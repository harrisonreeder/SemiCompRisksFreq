#' Fit Parametric Univariate Survival Model
#'
#' @inheritParams nll_func
#' @param Formula a Formula object, with the outcome on the left of a
#'   \code{~}, and covariates on the right. It is of the form, \code{time to non-terminal
#'   event + corresponding censoring indicator | time to terminal event
#'   + corresponding censoring indicator ~ covariates for h_1 |
#'   covariates for h_2 | covariates for h_3}. For example, \code{y_1 + delta_1 | y_2 + delta_2 ~ x_1 | x_2 | x_3}.
#' @param data a \code{data.frame} in which to interpret the variables named in Formula.
#' @param na.action how NAs are treated. See \code{\link[stats]{model.frame}}.
#' @param subset a specification of the rows to be used: defaults to all rows. See \code{\link[stats]{model.frame}}.
#' @param knots_vec Used for hazard specifications besides Weibull, a
#'   vector increasing sequence of integers, each corresponding to
#'   the knots for the flexible model on the baseline hazard. If
#'   \code{NULL}, will be created by \code{\link{get_default_knots}}.
#' @param p0 Integer indicating how many baseline hazard parameters
#'   should be specified for each of the three transition hazards. This input is only relevant when
#'   hazard is something other than \code{"weibull"} and is superceded by knots_vec.
#' @param startVals A numeric vector of parameter starting values, arranged as follows:
#'   the first \eqn{k_1+k_2+k_3} elements correspond to the baseline hazard parameters,
#'   then the last\eqn{q_1+q_2+q_3} elements correspond with the regression parameters.
#'   If set to \code{NULL}, will be generated automatically using \code{\link{get_start}}.
#' @param hessian Boolean indicating whether the hessian (aka, the inverse of the covariance matrix)
#'   should be computed and returned.
#' @param control a list of control attributes passed directly into the \code{optim} function.
#' @param n_quad Scalar for number of Gaussian quadrature points used to evaluate numerical integral of B-spline.
#' @param quad_method String indicating which quadrature method to use to evaluate numerical integral of B-spline.
#'   Options are \code{"kronrod"} for Gauss-Kronrod quadrature or \code{"legendre"} for Gauss-Legendre quadrature.
#' @param optim_method a string naming which \code{optim} method should be used for optimization.
#'
#' @return \code{FreqID_HReg2} returns an object of class \code{Freq_HReg}.
#' @import Formula
#' @export
FreqSurv_HReg2 <- function(Formula, data, na.action="na.fail", subset=NULL,
                         hazard=c("weibull"), knots_vec = NULL, p0=4,
                         startVals=NULL, hessian=TRUE, control=NULL, quad_method="kronrod", n_quad=15,
                         optim_method = "BFGS"){
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

  con=list(maxit=1000)
  nmsC <- names(con)
  namc <- names(control)
  con[namc] <- control

  ##DATA PREPROCESSING##
  ##******************##
  ##MAKE THIS MORE EFFICIENT BY GOING DIRECTLY INTO NUMERICS

  #This line ensures that the formula is of type Formula, and not just formula
  #the below manipulations require special methods of Formula.
  form2 <- Formula::as.Formula(paste0(Formula[2], Formula[1], Formula[3]))
  data <- stats::model.frame(form2, data = data, na.action = na.action,
                      subset = subset)

  #to account for possible left truncation, look on the left side of the formula
  #if there are two pieces on the left side, the first is the left truncation variable
  if(length(form2)[1] == 2){
    yL <- model.part(Formula, data = data, lhs = 1)[[1]]
    time1 <- Formula::model.part(Formula, data = data, lhs = 2)
  } else if(length(form2)[1] == 1){
    yL <- 0
    time1 <- Formula::model.part(Formula, data = data, lhs = 1)
  }
  y <- time1[[1]]
  delta <- time1[[2]]
  Xmat <- as.matrix(stats::model.frame(stats::formula(Formula, lhs = 0, rhs = 1),
                                       data = data))
  anyLT <- as.numeric(any(yL>0))

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
      basis <- get_basis(x=y,knots=knots_vec,hazard=hazard)
      dbasis <- get_basis(x=y,knots=knots_vec,hazard=hazard,deriv=TRUE)
      #to account for left truncation, directly subtract off the basis of yL
      if(anyLT){
        basis <- basis - get_basis(x=yL,knots=knots_vec,hazard=hazard)
      }
    } else if(hazard %in% c("bspline","bs")){
      if(anyLT){ #presence of left-truncation
        #note yL is lower bound of integral to account for left truncation
        y_quad <- c(y, transform_quad_points(n_quad=n_quad,
                                             quad_method=quad_method,
                                             a=yL,b=y))
      } else{
        y_quad <- c(y, transform_quad_points(n_quad=n_quad,
                                             quad_method=quad_method,
                                             a=0,b=y))
      }
      basis <- get_basis(x=y_quad,knots=knots_vec,hazard=hazard)
      basis_yL <- dbasis <- NULL
      attr(x=basis,which="quad_method") <- tolower(quad_method)
    } else if(hazard %in% c("royston-parmar","rp")){
      basis <- get_basis(x=y,knots=knots_vec,hazard=hazard)
      dbasis <- get_basis(x=y,knots=knots_vec,hazard=hazard,deriv=TRUE)
      #royston-parmar model requires separate matrix for left truncation
      basis_yL <- if(anyLT) get_basis(x=yL,knots=knots_vec,hazard=hazard) else NULL
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
    startVals <- get_start_uni(y=y,delta=delta,yL=yL,anyLT=anyLT,Xmat=Xmat,knots=knots_vec,
                               hazard=hazard,basis=basis)
  }

  #now, run the fitting function, which calls the correct optimization engine
  value <- get_fit_uni(startVals=startVals, y=y, delta=delta,
                      Xmat=Xmat,hazard=hazard,
                      basis=basis,dbasis=dbasis,
                      basis_yL=basis_yL,yL=yL,anyLT=anyLT,
                      control=con, hessian=hessian, optim_method=optim_method)
  #if the fit fails, then return
  if(!is.null(value$fail)){
    return(list(fail=TRUE,formula=form2,hazard=hazard,
                startVals=startVals,knots_vec=knots_vec,
                basis=basis,dbasis=dbasis,
                control=control,optim_method=optim_method))
  }

  #add the other quantities to the output
  value <- list(
    estimate=value$estimate,
    logLike=value$logLike,
    grad=value$grad,
    optim_details=value$optim_details,
    Finv= if(hessian) MASS::ginv(value$nhess) else NA,
    startVals=startVals,
    knots_vec=knots_vec,
    myLabels=myLabels,
    formula=form2,nP=nP,nP0=nP0,nobs=length(y),ymax=max(y),n_quad=n_quad,
    quad_method=quad_method,optim_method=optim_method)

  value$class <- c("Freq_HReg2","Surv","Ind",
                   switch(tolower(hazard),
                          weibull="Weibull",wb="Weibull",
                          bspline="B-Spline",bs="B-Spline",
                          "royston-parmar"="Royston-Parmar",rp="Royston-Parmar",
                          piecewise="Piecewise Constant",pw="Piecewise Constant"))
  class(value) <- "Freq_HReg2"
  return(value)

}
