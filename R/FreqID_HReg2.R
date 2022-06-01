#' Fit Parametric Frailty Illness-Death Model for Semi-Competing Risks Data
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
#' @param knots_list Used for hazard specifications besides Weibull, a
#'   list of three increasing sequences of integers, each corresponding to
#'   the knots for the flexible model on the corresponding transition baseline hazard. If
#'   \code{NULL}, will be created by \code{\link{get_default_knots_list}}.
#' @param nP0 vector of length three of integers indicating how many baseline hazard parameters
#'   should be specified for each of the three transition hazards. This input is only relevant when
#'   hazard is something other than "weibull" and is superceded by knots_list.
#' @param startVals A numeric vector of parameter starting values, arranged as follows:
#'   the first \eqn{k_1+k_2+k_3} elements correspond to the baseline hazard parameters,
#'   then the \eqn{k_1+k_2+k_3+1} element corresponds to the gamma frailty log-variance parameter,
#'   then the last\eqn{q_1+q_2+q_3} elements correspond with the regression parameters.
#'   If set to \code{NULL}, will be generated automatically using \code{\link{get_start}}.
#' @param hessian Boolean indicating whether the hessian (aka, the inverse of the covariance matrix)
#'   should be computed and returned.
#' @param control a list of control attributes passed directly into the \code{optim} function.
#' @param n_quad Scalar for number of Gaussian quadrature points used to evaluate numerical integral of B-spline.
#' @param quad_method String indicating which quadrature method to use to evaluate numerical integral of B-spline.
#'   Options are 'kronrod' for Gauss-Kronrod quadrature or 'legendre' for Gauss-Legendre quadrature.
#' @param optim_method a string naming which \code{optim} method should be used.
#'
#' @return \code{FreqID_HReg2} returns an object of class \code{Freq_HReg}.
#' @import Formula
#' @export
FreqID_HReg2 <- function(Formula, data, na.action="na.fail", subset=NULL,
                        hazard=c("weibull"), frailty=TRUE, model, knots_list=NULL,
                        nP0=rep(4,3), startVals=NULL, hessian=TRUE,
                        quad_method="kronrod", n_quad=15,
                        optim_method="BFGS", control=NULL){
  browser()

  ##INITIALIZE OPTIONS##
  ##******************##
  #technically na.fail I think is the only one currently implemented
  na.action <- match.arg(na.action)
  if (na.action != "na.fail" & na.action != "na.omit") {
    stop("na.action should be either na.fail or na.omit")
  }
  #Check that chosen hazard is among available options
  if(!(tolower(hazard) %in% c("weibull","royston-parmar","bspline",
                              "piecewise","wb","pw","bs","rp"))){
    stop("valid choices of hazard are 'weibull', 'royston-parmar', or 'piecewise'")
  } else{
    hazard <- tolower(hazard)
  }
  #prepare "control" list that is passed to optim function
  #by default sets maxit to 1000, and otherwise overwrites defaults with user inputs
  con=list(maxit=1000)
  nmsC <- names(con)
  namc <- names(control)
  con[namc] <- control

  ##DATA PREPROCESSING##
  ##******************##

  #rearrange input Formula object (which stores the different pieces of the input formula)
  #This line ensures that the formula is of type Formula, and not just formula
  #the below manipulations require special methods of Formula.
  form2 <- Formula::as.Formula(paste0(Formula[2], Formula[1], Formula[3]))
  #arrange data into correct form, extracting subset of variables
  #and observations needed in model
  data <- stats::model.frame(form2,data=data,na.action=na.action,subset=subset)
  #create matrices storing two outcomes, and then component vectors
  time1 <- Formula::model.part(form2, data=data, lhs=1)
  time2 <- Formula::model.part(form2, data=data, lhs=2)
  y1 <- time1[[1]]
  delta1 <- time1[[2]]
  y2 <- time2[[1]]
  delta2 <- time2[[2]]
  #Create covariate matrices for each of three transition hazards
  Xmat1 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=1),
                                        data=data))
  Xmat2 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=2),
                                        data=data))
  Xmat3 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=3),
                                        data=data))

  ##PREPARE KNOTS AND BASIS FUNCTIONS FOR FLEXIBLE MODELS##
  ##*****************************************************##
  if(hazard %in% c("bspline","royston-parmar","piecewise","pw","rp","bs")){
    #if the knots (or changepoints) have not been provided, get the "defaults"
    if(is.null(knots_list)){
      if(length(nP0) != 3){stop("If hazard not equal to 'weibull' and knots_list set to NULL, then nP0 must be vector of 3 integers.")}
      p01 <- nP0[1]; p02 <- nP0[2]; p03 <- nP0[3]
      knots_list <- get_default_knots_list(y1,y2,delta1,delta2,
                                           p01,p02,p03,hazard,model)
    }

    #now, each spline specification has it's own details regarding generating
    #basis functions, so we go through them one by one.
    if(hazard %in% c("royston-parmar","rp")){
      basis1 <- get_basis(x=y1,knots=knots_list[[1]],hazard=hazard)
      basis2 <- get_basis(x=y1,knots=knots_list[[2]],hazard=hazard)
      #royston-parmar model also uses a set of derivative basis functions
      dbasis1 <- get_basis(x=y1,knots=knots_list[[1]],hazard=hazard,deriv=TRUE)
      dbasis2 <- get_basis(x=y1,knots=knots_list[[2]],hazard=hazard,deriv=TRUE)
      if(tolower(model)=="semi-markov"){
        basis3 <- get_basis(x=y2-y1,knots=knots_list[[3]],hazard=hazard)
        dbasis3 <- get_basis(x=y2-y1,knots=knots_list[[3]],hazard=hazard,deriv=TRUE)
        basis3_y1 <-  NULL
      } else{
        basis3 <- get_basis(x=y2,knots=knots_list[[3]],hazard=hazard)
        dbasis3 <- get_basis(x=y2,knots=knots_list[[3]],hazard=hazard,deriv=TRUE)
        #royston-parmar model requires separate basis3_y1 matrix
        basis3_y1 <-  get_basis(x=y1,knots=knots_list[[3]],hazard=hazard)
      }

    } else if(hazard %in% c("piecewise","pw")){
      basis1 <- get_basis(x=y1,knots=knots_list[[1]],hazard=hazard)
      basis2 <- get_basis(x=y1,knots=knots_list[[2]],hazard=hazard)
      #piecewise constant model also uses a set of derivative basis functions
      dbasis1 <- get_basis(x=y1,knots=knots_list[[1]],hazard=hazard,deriv=TRUE)
      dbasis2 <- get_basis(x=y1,knots=knots_list[[2]],hazard=hazard,deriv=TRUE)
      if(tolower(model)=="semi-markov"){
        basis3 <- get_basis(x=y2-y1,knots=knots_list[[3]],hazard=hazard)
        dbasis3 <- get_basis(x=y2-y1,knots=knots_list[[3]],hazard=hazard,deriv=TRUE)
      } else{
        #for markov piecewise constant model, we simplify by setting
        #basis3 as the difference of bases for y2 and y1
        basis3 <- get_basis(x=y2,knots=knots_list[[3]],hazard=hazard) -
          get_basis(x=y1,knots=knots_list[[3]],hazard=hazard)
        dbasis3 <- get_basis(x=y2,knots=knots_list[[3]],hazard=hazard,deriv=TRUE)
        basis3_y1 <- NULL
      }

    } else if(hazard %in% c("bspline","bs")){
      #for bspline, evaluate bases at observed times
      #as well as all quadrature points for numerical integration
      y1_quad <- c(y1,transform_quad_points(n_quad=n_quad,
                                            quad_method=quad_method,
                                            a=0,b=y1))
      basis1 <- get_basis(x=y1_quad, knots=knots_list[[1]],hazard=hazard)
      basis2 <- get_basis(x=y1_quad, knots=knots_list[[2]],hazard=hazard)
      if(tolower(model)=="semi-markov"){
        y2y1_quad <- c(y2-y1, transform_quad_points(n_quad=n_quad,
                                                    quad_method=quad_method,
                                                    a=0,b=y2-y1))
        basis3 <- get_basis(x=y2y1_quad,knots=knots_list[[3]],hazard=hazard)
        dbasis3 <- NULL
      } else{ #markov
        y2_quad <- c(y2, transform_quad_points(n_quad=n_quad,
                                                 quad_method=quad_method,
                                                 a=y1,b=y2))
        basis3 <- get_basis(x=y2_quad,knots=knots_list[[3]],hazard=hazard)
        basis3_y1 <- dbasis3 <- NULL
      }
      #lastly, we add attribute storing quadrature method
      attr(x=basis1,which="quad_method") <-
        attr(x=basis2,which="quad_method") <-
        attr(x=basis3,which="quad_method") <- tolower(quad_method)
    }
    p01 <- ncol(basis1); p02 <- ncol(basis2); p03 <- ncol(basis3)

  } else{ #if not a flexible method, must be weibull
    basis1 <- basis2 <- basis3 <- basis3_y1 <-
      dbasis1 <- dbasis2 <- dbasis3 <- NULL
    p01 <- p02 <- p03 <- 2
  }

  ##Store a few values for the output##
  ##*********************************##
  nP <- c(ncol(Xmat1), ncol(Xmat2), ncol(Xmat3))
  nP0 <- c(p01,p02,p03)
  if(tolower(hazard) %in% c("weibull","wb")){
    myLabels <- c("log(kappa1)", "log(alpha1)", "log(kappa2)",
                  "log(alpha2)", "log(kappa3)", "log(alpha3)")
  } else{
    myLabels <- c(paste0("phi1",1:p01),paste0("phi2",1:p02),paste0("phi3",1:p03))
  }

  #if the user has not provided start values, we generate them here
  #we start by fitting non-frailty model, so startVals is generated without frailty variance
  if(is.null(startVals)){
    startVals_nf <- get_start(y1=y1,y2=y2,delta1=delta1,delta2=delta2,
                   Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,
                   hazard=hazard,frailty=FALSE,model=model,knots_list=knots_list,
                   basis1=basis1,basis2=basis2,basis3=basis3)

    # nll_ltheta_func <- function(x){
    #   nll_func(para=c(startVals_nf[1:sum(nP0)],x,startVals_nf[-(1:sum(nP0))]),
    #          y1=y1, y2=y2,
    #          delta1=delta1, delta2=delta2,
    #          Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
    #          hazard=hazard,frailty=frailty,model=model,
    #          basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
    #          dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)}
    # start_ltheta <- stats::optimize(f = nll_ltheta_func,interval = c(-10,5))$minimum
    start_ltheta <- log(1)

    startVals <- c(startVals_nf[1:sum(nP0)],start_ltheta,startVals_nf[-(1:sum(nP0))])
  } else{
    startVals_nf <- if(frailty) startVals[-(nP0+1)] else startVals
  }

  ##FIT MODEL##
  ##*********##
  if(!frailty){
    #instead of fitting a non-frailty model directly, we fit it as three
    #univariate models for the three transitions.
    value <- get_nf_fit(startVals_nf=startVals_nf,
                y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                hazard=hazard,model=model,
                basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                control=con, hessian=hessian,optim_method=optim_method)
    #compute and store final gradient of log likelihood
    value$grad <- -ngrad_func(para=value$estimate,y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                      hazard=hazard,frailty=FALSE,model=model,
                      basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                      dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)

  } else{

    #create list that we will fill with results
    value <- list()


    #two options for optimization approaches: optim(), and nleqslv()
    if(tolower(optim_method)=="nleqslv"){

      fit0 <- tryCatch(nleqslv::nleqslv(x=startVals, fn=ngrad_func,
                                    method="Broyden",global = "dbldog",
                                    y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                    Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                    hazard=hazard,frailty=frailty,model=model,
                                    basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                    dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                    control=con, jacobian=hessian),
                       error=function(cnd){message(cnd); cat("\n")
                         return(list(fail=TRUE, formula=form2,
                                     hazard=hazard,frailty=frailty,model=model,
                                     startVals=startVals,knots_list=knots_list,
                                     basis1=basis1,basis2=basis2,basis3=basis3,
                                     dbasis1=dbasis1,dbasis2=dbasis2,dbasis3=dbasis3,
                                     control=control))})
      #If optimization fails, return list of inputs and flag "fail=TRUE"
      if(!is.null(fit0$fail)){ return(fit0)}
      if (!(fit0$termcd %in% c(1,2))){warning("check convergence.")}

      #prepare the results from the nleqslv function into the output format
      value$estimate <- fit0$x
      value$Finv <- if(hessian) MASS::ginv(fit0$jac) else NA
      value$logLike <- -nll_func(para=value$estimate,y1=y1, y2=y2,
                                 delta1=delta1, delta2=delta2,
                                 Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                 hazard=hazard,frailty=frailty,model=model,
                                 basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                 dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
      #WHAT TO DO IF THIS IS NaN, which has occurred for royston-parmar model

      value$optim_details <- list(counts=fit0$iter,
                                  convergence=fit0$termcd,
                                  message=fit0$message)
      value$grad <- -fit0$fvec

    } else{ #use method from optim()

      fit0 <- tryCatch(stats::optim(par=startVals, fn=nll_func, gr=ngrad_func,
                                    y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                    Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                    hazard=hazard,frailty=frailty,model=model,
                                    basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                    dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                    control=con, hessian=hessian,method=optim_method),
                       error=function(cnd){message(cnd); cat("\n")
                         return(list(fail=TRUE, formula=form2,
                                     hazard=hazard,frailty=frailty,model=model,
                                     startVals=startVals,knots_list=knots_list,
                                     basis1=basis1,basis2=basis2,basis3=basis3,
                                     dbasis1=dbasis1,dbasis2=dbasis2,dbasis3=dbasis3,
                                     control=control))})
      #If optimization fails, return list of inputs and flag "fail=TRUE"
      if(!is.null(fit0$fail)){ return(fit0)}
      if (!(fit0$convergence %in% c(0,1))){warning("check convergence.")}

      #prepare the results from the optim function into the output format
      value$estimate <- fit0$par
      value$Finv <- if(hessian) MASS::ginv(fit0$hessian) else NA
      value$logLike <- -fit0$value
      value$optim_details <- list(counts=fit0$counts,
                                  convergence=fit0$convergence,
                                  message=fit0$message)
      #compute and store final gradient of log likelihood
      value$grad <- as.vector(-ngrad_func(para=value$estimate,y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                hazard=hazard,frailty=frailty,model=model,
                                basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3))
    }
    value$startVals <- startVals
    myLabels <- c(myLabels, "log(theta)")



    ##FINALLY, OPTIONALLY COMPUTE NON-FRAILTY MODEL FOR COMPARISON##
    #First, maximize the likelihood without a frailty using univariate functions
    value$fit_nf <- get_nf_fit(startVals_nf=startVals_nf,
                               y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                               Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                               hazard=hazard,model=model,
                               basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                               dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                               control=con, hessian=hessian,optim_method=optim_method)
    #compute and store final gradient of log likelihood (non-frailty model)
    value$fit_nf$grad <- as.vector(-ngrad_func(para=value$fit_nf$estimate,y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                     Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                     hazard=hazard,frailty=FALSE,model=model,
                                     basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                     dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3))
    #compute a likelihood-ratio test for the frailty
    #if the frailty likelihood is below the non-frailty likelihood, set to 0
    frail_test <- max(0,2*(value$logLike - value$fit_nf$logLike))
    if(frail_test >= 0){
      value$frailty_test <- c(
        stat=frail_test,
        #corrected null distribution
        pval= 0.5 * (stats::pchisq(q = frail_test,df = 0,lower.tail = FALSE) +
                       stats::pchisq(q = frail_test,df = 1,lower.tail = FALSE))
      )
    }
  }

  ##PREPARE OUTPUTS##
  ##***************##
  myLabels <- c(myLabels,colnames(Xmat1),colnames(Xmat2),colnames(Xmat3))
  value$knots_list <- knots_list
  value$myLabels <- myLabels
  value$frailty <- frailty
  value$formula <- form2
  value$nP <- nP
  value$nP0 <- nP0
  value$nobs <- length(y1)
  value$ymax <- max(y2)
  value$n_quad <- n_quad
  value$quad_method <- quad_method


  value$class <- c("Freq_HReg2","ID","Ind",
                   switch(tolower(hazard),
                          weibull="Weibull",wb="Weibull",
                          bspline="B-Spline",bs="B-Spline",
                          "royston-parmar"="Royston-Parmar",rp="Royston-Parmar",
                          piecewise="Piecewise Constant",pw="Piecewise Constant"),
                   switch(tolower(model),
                          "semi-markov"="semi-Markov",
                          "markov"="Markov"))
  class(value) <- "Freq_HReg2"
  return(value)
}

