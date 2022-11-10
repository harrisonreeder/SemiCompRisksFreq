#' Fit Parametric Univariate Survival Model
#'
#' @inheritParams nll_func
#' @param Formula a Formula object, with the outcome on the left of a
#'   \code{~}, and covariates on the right. It is of the form,
#'   \code{left truncation time | time to event +
#'   corresponding censoring indicator ~ covariates}.
#'   For example, \code{y_L | y + delta ~ x_1 + x_2 + x_3}.
#'   If there is no left truncation, then the lefthand side should only contain
#'   the two sets of outcome variables.
#' @param data a \code{data.frame} in which to interpret
#'   the variables named in Formula.
#' @param na.action how NAs are treated. See \code{\link[stats]{model.frame}}.
#' @param subset a specification of the rows to be used: defaults to all rows.
#'   See \code{\link[stats]{model.frame}}.
#' @param knots_vec Used for hazard specifications besides Weibull, a
#'   vector increasing sequence of integers, each corresponding to
#'   the knots for the flexible model on the baseline hazard.
#'   If \code{NULL}, will be created by \code{\link{get_default_knots}}
#'   according to the number of parameters specified by \code{nP0}.
#' @param p0 Integer indicating how many baseline hazard parameters
#'   should be specified for each of the three transition hazards.
#'   This input is only relevant when hazard is something other than
#'   \code{"weibull"} and is superceded by knots_vec.
#' @param startVals A numeric vector of parameter starting values,
#'   arranged as follows:
#'   the first \eqn{k} elements are the baseline hazard parameters,
#'   then the last\eqn{q} elements are the regression parameters.
#'   If set to \code{NULL}, generated internally using \code{\link{get_start}}.
#' @param control a list of control attributes passed directly
#'   into the \code{optim} function.
#' @param n_quad Scalar for number of Gaussian quadrature points used to
#'   evaluate numerical integral of B-spline.
#' @param quad_method String indicating which quadrature method to use to
#'   evaluate numerical integral of B-spline. Options are
#'   \code{"kronrod"} for Gauss-Kronrod quadrature or
#'   \code{"legendre"} for Gauss-Legendre quadrature.
#' @param optim_method a string naming which \code{optim}
#'   method should be used for optimization.
#' @param extra_starts Integer giving the number of extra
#'   starts to try when optimizing.
#' @param output_options List of named boolean elements specifying whether
#'   certain additional components should be included in
#'   the model object output. Options include
#'   \itemize{
#'     \item{Finv}{Variance-covariance matrix. Defaults to \code{TRUE}.}
#'     \item{grad_mat_return}{Matrix with rowwise
#'     score vectors for each individual evaluated at the MLE.
#'     Used to compute "cheese" or "meat" in robust standard error computation.
#'     Defaults to \code{FALSE}.}
#'     \item{cheese}{Sum of outer products of individual score vectors,
#'     used as the "cheese" or "meat" in robust standard error computation.
#'     Defaults to \code{TRUE}.}
#'     \item{data_return}{Original data frame used to fit model.
#'     Defaults to \code{FALSE}.}
#'   }
#'
#' @return An object of class \code{Freq_HReg}.
#' @import Formula
#'
#' @examples
#' #loading a data set
#' data(scrData)
#'
#' #fitting Weibull survival model on terminal event
#' form <- Formula::Formula(time2 + event2 ~ x1 + x2 + x3 | x1 + x2 | x1 + x2)
#' fit_WB	<- FreqSurv_HReg2(Formula = form, data=scrData,
#' extra_starts = 0,hazard = "weibull",optim_method = c("BFGS"))
#'
#' #exploring results
#' fit_WB
#' summ.fit_WB <- summary(fit_WB); names(summ.fit_WB)
#' summ.fit_WB
#' pred_WB <- predict(fit_WB, tseq=seq(from=0.1, to=30, length.out=100))
#' plot(pred_WB, plot.est="Haz")
#' plot(pred_WB, plot.est="Surv")
#'
#' @export
FreqSurv_HReg2 <- function(Formula, data, na.action="na.fail", subset=NULL,
                           weights=NULL, hazard=c("weibull"), knots_vec = NULL, p0=4,
                           startVals=NULL, control=NULL, quad_method="kronrod", n_quad=15,
                           optim_method = "BFGS", extra_starts=0, output_options=NULL){
  # browser()
  ##Check that chosen hazard is among available options
  if(!(tolower(hazard) %in% c("weibull","royston-parmar","bspline","piecewise","wb","pw","bs","rp"))){
    stop("valid choices of hazard are 'weibull', 'royston-parmar', or 'piecewise'")
  } else{
    hazard <- tolower(hazard)
  }

  ##INITIALIZE OPTIONS##
  ##******************##
  na.action <- match.arg(na.action) #technically na.fail I think is the only one currently implemented
  if (na.action != "na.fail" & na.action != "na.omit") {
    stop("na.action should be either na.fail or na.omit")
  }
  #if any control options have been given, override the defaults
  con=list(maxit=1000)
  nmsC <- names(con)
  namc <- names(control)
  con[namc] <- control

  #if any output options have been given, override the defaults
  out_options=list(Finv=TRUE,grad_mat_return=FALSE,
                   cheese=TRUE,data_return=FALSE)
  nmsO <- names(out_options)
  namO <- names(output_options)
  out_options[namO] <- output_options

  ##DATA PREPROCESSING##
  ##******************##

  #This line ensures that the formula is of type Formula, and not just formula
  #the below manipulations require special methods of Formula.
  form2 <- Formula::as.Formula(paste0(Formula[2], Formula[1], Formula[3]))
  data <- stats::model.frame(form2, data = data, na.action = na.action,
                      subset = subset)

  ##ACCOUNT FOR LEFT TRUNCATION WITH SUBSET OPTION (AND OTHER INPUTS TOO!)

  #to account for possible left truncation, look on the left side of the formula
  #if there are two pieces on the left side, the first is the left truncation variable
  if(length(form2)[1] == 2){
    yL <- model.part(form2, data = data, lhs = 1)[[1]]
    time1 <- Formula::model.part(form2, data = data, lhs = 2)
  } else if(length(form2)[1] == 1){
    yL <- 0
    time1 <- Formula::model.part(form2, data = data, lhs = 1)
  }
  anyLT <- as.numeric(any(yL>0))

  y <- time1[[1]]
  delta <- time1[[2]]
  Xmat <- as.matrix(stats::model.frame(stats::formula(form2, lhs = 0, rhs = 1),
                                       data = data))

  if(is.null(weights)){
    weights <- rep(1,length(y))
  } else{
    stopifnot(length(weights) == length(y))
  }

  ##PREPARE KNOTS AND BASIS FUNCTIONS FOR FLEXIBLE MODELS##
  ##*****************************************************##

  if(hazard %in% c("bspline","royston-parmar","piecewise","pw","rp","bs")){
    if(is.null(knots_vec)){
      knots_vec <- get_default_knots(y = y,delta = delta,p0 = p0,hazard = hazard)
    }

    # now, each spline specification has it's own details regarding generating
    # basis functions (esp. to account for left-truncation),
    # so we go through them one by one.
    if(hazard %in% c("piecewise","pw")){
      basis <- get_basis(y=y,knots_vec=knots_vec,hazard=hazard)
      dbasis <- get_basis(y=y,knots_vec=knots_vec,hazard=hazard,deriv=TRUE)
      #to account for left truncation, directly subtract off the basis of yL
      if(anyLT){
        basis <- basis - get_basis(y=yL,knots_vec=knots_vec,hazard=hazard)
      }
    } else if(hazard %in% c("bspline","bs")){
      if(anyLT){ #presence of left-truncation
        #note yL is lower bound of integral to account for left truncation
        y_quad <- c(y, transform_quad_points(n_quad=n_quad, a=yL,b=y,
                                             quad_method=quad_method))
      } else{
        y_quad <- c(y, transform_quad_points(n_quad=n_quad, a=0,b=y,
                                             quad_method=quad_method))
      }
      basis <- get_basis(y=y_quad,knots_vec=knots_vec,hazard=hazard)
      basis_yL <- dbasis <- NULL
      attr(x=basis,which="quad_method") <- tolower(quad_method)
    } else if(hazard %in% c("royston-parmar","rp")){
      basis <- get_basis(y=y,knots_vec=knots_vec,hazard=hazard)
      dbasis <- get_basis(y=y,knots_vec=knots_vec,hazard=hazard,deriv=TRUE)
      #royston-parmar model requires separate matrix for left truncation
      basis_yL <- if(anyLT) get_basis(y=yL,knots_vec=knots_vec,hazard=hazard) else NULL
    }
  } else{
    basis_yL <- basis <- dbasis <- NULL
    p0 <- 2 #must be weibull
  }

  #in bspline loglikelihood, passed in y is only used to define "weights"
  #for quadrature, which now must account for change of bounds above.
  #So, we just redefine y to be this difference.
  if(hazard %in% c("bspline","bs") & anyLT){ y <- y-yL }

  #finalize the labels of the parameters
  if(tolower(hazard) %in% c("weibull","wb")){
    myLabels <- c("log(kappa)", "log(alpha)")
  } else{
    myLabels <- c(paste0("phi",1:p0))
  }
  myLabels <- c(myLabels, colnames(Xmat))
  nP <- ncol(Xmat)
  nP0 <- p0


  ##GET MLE##
  ##*******##

  #if the user has not provided start values, we generate them here
  if(is.null(startVals)){
    startVals <- get_start_uni(y=y,delta=delta,yL=yL,anyLT=anyLT,Xmat=Xmat,
                               knots_vec=knots_vec,hazard=hazard,basis=basis,weights=weights)
  }

  grad1 <- ngrad_uni_func(para=startVals, y=y,delta=delta,yL=yL,anyLT=anyLT,Xmat=Xmat,
                          weights=weights, hazard=hazard,basis=basis, dbasis=dbasis,basis_yL=basis_yL)
  grad2 <- pracma::grad(f = nll_uni_func,x0 = startVals,y=y,delta=delta,yL=yL,
                        anyLT=anyLT,Xmat=Xmat,hazard=hazard, weights=weights,
                        basis=basis, dbasis=dbasis,basis_yL=basis_yL)
  grad3 <- colSums(ngrad_uni_mat_func(para=startVals, y=y,delta=delta,yL=yL,anyLT=anyLT,Xmat=Xmat,
                   weights=weights,hazard=hazard,basis=basis, dbasis=dbasis,basis_yL=basis_yL))
  if(max(abs(grad1-grad2)) >= 1e-4){stop("check gradient of non-frailty model")}
  if(max(abs(grad1-grad3)) >= 1e-4){stop("check gradient of non-frailty model")}
  cbind(grad1,grad2,grad3)

  #now, run the fitting function, which calls the correct optimization engine
  value <- get_fit_uni(startVals=startVals, y=y, delta=delta,
                      Xmat=Xmat,hazard=hazard, weights=weights,
                      basis=basis,dbasis=dbasis,
                      basis_yL=basis_yL,yL=yL,anyLT=anyLT, control=con,
                      hessian=out_options$Finv, optim_method=optim_method,
                      extra_starts=extra_starts)

  #if the fit fails, then return generic data
  if(!is.null(value$fail)){
    return(list(fail=TRUE,formula=form2,hazard=hazard,
                startVals=startVals,knots_vec=knots_vec,
                basis=basis,dbasis=dbasis,
                control=control,optim_method=optim_method,
                extra_starts=extra_starts))
  }

  #if requested, compute the inverse hessian (aka sandwich "bread")
  if(out_options$Finv){
    Finv <- tryCatch(MASS::ginv(value$nhess),
                    error=function(cnd){message(cnd);cat("\n");return(NA)})
  } else{ Finv <- NA }

  #if requested, compute gradient contributions for every subject
  if(out_options$grad_mat_return){
    #note division by n following `sandwich::meat' function
    grad_mat <- -ngrad_uni_mat_func(para = value$estimate,
                                         y=y, delta=delta, yL=yL, anyLT=anyLT,
                                         Xmat=Xmat, hazard=hazard,
                                         basis=basis, basis_yL=basis_yL,
                                         dbasis=dbasis, weights=weights)
  } else{ grad_mat <- NULL }

  #if requested, compute the sandwich variance "cheese" outer product of scores
  if(out_options$cheese){
    if(out_options$grad_mat_return){
      cheese <- crossprod(grad_mat)
    } else{
      cheese <- crossprod(ngrad_uni_mat_func(para = value$estimate,
                                           y=y, delta=delta, yL=yL, anyLT=anyLT,
                                           Xmat=Xmat, hazard=hazard,
                                           basis=basis, basis_yL=basis_yL,
                                           dbasis=dbasis, weights=weights))
      }
  } else{ cheese <- NA }

  #add the other quantities to the output
  value <- list(
    estimate=as.vector(value$estimate),
    logLike=value$logLike,
    grad=as.vector(value$grad),
    optim_details=value$optim_details,
    Finv = Finv,
    cheese = cheese,
    startVals=value$startVals,
    knots_vec=knots_vec,
    myLabels=myLabels,
    formula=form2,nP=nP,nP0=nP0,nobs=length(y),ymax=max(y),n_quad=n_quad,
    quad_method=quad_method,optim_method=optim_method,extra_starts=extra_starts,
    control=con,
    grad_mat=grad_mat,
    data=if(out_options$data_return) data else NULL)

  names(value$estimate) <- names(startVals)
  names(value$grad) <- names(startVals)

  value$class <- c("Freq_HReg2","Surv","Ind",
                   switch(tolower(hazard),
                      weibull="Weibull",wb="Weibull",
                      bspline="B-Spline",bs="B-Spline",
                      "royston-parmar"="Royston-Parmar",rp="Royston-Parmar",
                      piecewise="Piecewise Constant",pw="Piecewise Constant"))
  class(value) <- "Freq_HReg2"
  return(value)

}
