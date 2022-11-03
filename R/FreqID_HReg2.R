#' Fit Parametric Frailty Illness-Death Model for Semi-Competing Risks Data
#'
#' @inheritParams nll_func
#' @inheritParams FreqSurv_HReg2
#' @param Formula a Formula object, with the outcome on the left of a
#'   \code{~}, and covariates on the right. It is of the form,
#'   \code{left truncation time | time to non-terminal
#'   event + corresponding censoring indicator | time to terminal event
#'   + corresponding censoring indicator ~ covariates for h_1 |
#'   covariates for h_2 | covariates for h_3}.
#'   For example, \code{y_L | y_1 + delta_1 | y_2 + delta_2 ~ x_1 | x_2 | x_3}.
#'   If there is no left truncation, then the lefthand side should only contain
#'   the two sets of outcome variables.
#' @param knots_list Used for hazard specifications besides Weibull, a
#'   list of three increasing sequences of integers, each corresponding to
#'   the knots for the flexible model on the corresponding transition baseline hazard.
#'   If \code{NULL}, will be created by \code{\link{get_default_knots_list}}
#'   according to the number of parameters specified by \code{nP0}.
#' @param nP0 vector of length three of integers indicating how many
#'   baseline hazard parameters
#'   should be specified for each of the three transition hazards.
#'   This input is only relevant when
#'   hazard is something other than "weibull" and is superceded by knots_list.
#' @param startVals A numeric vector of parameter starting values,
#'   arranged as follows:
#'   the first \eqn{k_1+k_2+k_3} elements are the baseline hazard parameters,
#'   then the gamma frailty log-variance parameter,
#'   then the last\eqn{q_1+q_2+q_3} elements are the regression parameters.
#'   If set to \code{NULL}, generated internally using \code{\link{get_start}}.
#' @param output_options List of named boolean elements specifying whether
#'   certain additional components should be included in
#'   the model object output. Options include
#'   \itemize{
#'     \item{Finv}{Variance-covariance matrix. Defaults to \code{TRUE}.}
#'     \item{nf_fit}{If \code{frailty} is \code{TRUE}, include an
#'     additional fit object for the corresponding non-frailty model.
#'     Defaults to \code{TRUE}.}
#'     \item{eb_frailties}{Include vector of empirical Bayes
#'     predicted gamma frailty values. Defaults to \code{TRUE}.}
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
#' @return \code{FreqID_HReg2} returns an object of class \code{Freq_HReg}.
#' @import Formula
#' @export
FreqID_HReg2 <- function(Formula, data, na.action="na.fail", subset=NULL, weights=NULL,
                        hazard=c("weibull"), frailty=TRUE, model, knots_list=NULL,
                        nP0=rep(4,3), startVals=NULL, control=NULL,
                        quad_method="kronrod", n_quad=15,
                        optim_method="BFGS", extra_starts=0, output_options=NULL){
  # browser()
  ##INITIALIZE OPTIONS##
  ##******************##
  #technically na.fail I think is the only one currently implemented
  na.action <- match.arg(na.action)
  if (na.action != "na.fail" & na.action != "na.omit") {
    stop("na.action should be either na.fail or na.omit")
  }
  #Check that chosen hazard is among available options
  hazard_options <- c("weibull","wb","piecewise","pw","royston-parmar","rp","bspline","bs")
  if(!(tolower(hazard) %in% hazard_options)){
    stop(paste0("valid choices of hazard are '", paste(hazard_options,collapse = "', '"),"'"))
  } else{ hazard <- tolower(hazard) }
  #prepare "control" list that is passed to optim function
  #by default sets maxit to 1000, and otherwise overwrites defaults with user inputs
  con=list(maxit=1000)
  nmsC <- names(con)
  namc <- names(control)
  con[namc] <- control

  #if any output options have been given, override the defaults
  out_options=list(nf_fit=TRUE,Finv=TRUE,grad_mat_return=FALSE,
                   cheese=TRUE,eb_frailties=TRUE,data_return=FALSE)
  nmsO <- names(out_options)
  namO <- names(output_options)
  out_options[namO] <- output_options

  ##DATA PREPROCESSING##
  ##******************##
  #rearrange input Formula object (which stores the different pieces of the input formula)
  #This line ensures that the formula is of type Formula, and not just formula
  #the below manipulations require special methods of Formula.
  form2 <- Formula::as.Formula(paste0(Formula[2], Formula[1], Formula[3]))
  #arrange data into correct form, extracting subset of variables
  #and observations needed in model
  data <- stats::model.frame(form2,data=data,na.action=na.action,subset=subset)

  #to account for possible left truncation, look on the left side of the formula
  #if there are two pieces on the left side, the first is the left truncation variable
  if(length(form2)[1] == 3){
    yL <- model.part(form2, data = data, lhs = 1)[[1]]
    #create matrices storing two outcomes, and then component vectors
    time1 <- Formula::model.part(form2, data=data, lhs=2)
    time2 <- Formula::model.part(form2, data=data, lhs=3)
  } else if(length(form2)[1] == 2){
    yL <- 0
    time1 <- Formula::model.part(form2, data=data, lhs=1)
    time2 <- Formula::model.part(form2, data=data, lhs=2)
  }
  anyLT <- as.numeric(any(yL>0))
  y1 <- time1[[1]]
  delta1 <- time1[[2]]
  y2 <- time2[[1]]
  delta2 <- time2[[2]]
  #Create covariate matrices for each of three transition hazards
  Xmat1 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=1),data=data))
  Xmat2 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=2),data=data))
  Xmat3 <- as.matrix(stats::model.frame(stats::formula(form2, lhs=0, rhs=3),data=data))

  if(is.null(weights)){
    weights <- rep(1,length(y1))
  } else{
    stopifnot(length(weights) == length(y1))
  }

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
      basis1 <- get_basis(y=y1,knots_vec=knots_list[[1]],hazard=hazard)
      basis2 <- get_basis(y=y1,knots_vec=knots_list[[2]],hazard=hazard)
      if(anyLT){
        basis1_yL <- get_basis(y=yL,knots_vec=knots_list[[1]],hazard=hazard)
        basis2_yL <- get_basis(y=yL,knots_vec=knots_list[[2]],hazard=hazard)
      } else{
        basis1_yL <- NULL
        basis2_yL <- NULL
      }
      #royston-parmar model also uses a set of derivative basis functions
      dbasis1 <- get_basis(y=y1,knots_vec=knots_list[[1]],hazard=hazard,deriv=TRUE)
      dbasis2 <- get_basis(y=y1,knots_vec=knots_list[[2]],hazard=hazard,deriv=TRUE)
      if(tolower(model)=="semi-markov"){
        basis3 <- get_basis(y=y2-y1,knots_vec=knots_list[[3]],hazard=hazard)
        dbasis3 <- get_basis(y=y2-y1,knots_vec=knots_list[[3]],hazard=hazard,deriv=TRUE)
        basis3_y1 <-  NULL
      } else{
        basis3 <- get_basis(y=y2,knots_vec=knots_list[[3]],hazard=hazard)
        dbasis3 <- get_basis(y=y2,knots_vec=knots_list[[3]],hazard=hazard,deriv=TRUE)
        #royston-parmar model requires separate basis3_y1 matrix
        basis3_y1 <-  get_basis(y=y1,knots_vec=knots_list[[3]],hazard=hazard)
      }
    } else if(hazard %in% c("piecewise","pw")){
      basis1 <- get_basis(y=y1,knots_vec=knots_list[[1]],hazard=hazard)
      basis2 <- get_basis(y=y1,knots_vec=knots_list[[2]],hazard=hazard)
      if(anyLT){
        basis1_yL <- get_basis(y=yL,knots_vec=knots_list[[1]],hazard=hazard)
        basis2_yL <- get_basis(y=yL,knots_vec=knots_list[[2]],hazard=hazard)
      } else{
        basis1_yL <- NULL
        basis2_yL <- NULL
      }
      #piecewise constant model also uses a set of derivative basis functions
      dbasis1 <- get_basis(y=y1,knots_vec=knots_list[[1]],hazard=hazard,deriv=TRUE)
      dbasis2 <- get_basis(y=y1,knots_vec=knots_list[[2]],hazard=hazard,deriv=TRUE)
      if(tolower(model)=="semi-markov"){
        basis3 <- get_basis(y=y2-y1,knots_vec=knots_list[[3]],hazard=hazard)
        dbasis3 <- get_basis(y=y2-y1,knots_vec=knots_list[[3]],hazard=hazard,deriv=TRUE)
      } else{
        #for markov piecewise constant model, we simplify by setting
        #basis3 as the difference of bases for y2 and y1
        basis3 <- get_basis(y=y2,knots_vec=knots_list[[3]],hazard=hazard) -
          get_basis(y=y1,knots_vec=knots_list[[3]],hazard=hazard)
        dbasis3 <- get_basis(y=y2,knots_vec=knots_list[[3]],hazard=hazard,deriv=TRUE)
        basis3_y1 <- NULL
      }
    } else if(hazard %in% c("bspline","bs")){
      #for bspline, evaluate bases at observed times
      #as well as all quadrature points for numerical integration
      #plan is the following:
        #first, for any non-frailty model fits we will directly create basis1 and basis2
        #incorporating yL as a left truncation time
        #then, before fitting the frailty model, redefine basis1 and basis2 without left truncation
        #and also generate separate basis1_yL and basis2_yL accordingly (this is because of how frailty is computed).
      y1_quad <- c(y1,transform_quad_points(n_quad=n_quad,
                                            quad_method=quad_method,
                                            a=yL,b=y1)) #note that for a start, incorporate left truncation in here
      basis1 <- get_basis(y=y1_quad, knots_vec=knots_list[[1]],hazard=hazard)
      basis2 <- get_basis(y=y1_quad, knots_vec=knots_list[[2]],hazard=hazard)
      if(tolower(model)=="semi-markov"){
        y2y1_quad <- c(y2-y1, transform_quad_points(n_quad=n_quad,a=0,b=y2-y1,
                                                    quad_method=quad_method))
        basis3 <- get_basis(y=y2y1_quad,knots_vec=knots_list[[3]],hazard=hazard)
        dbasis3 <- NULL
      } else{ #markov (i.e., time to y2 treating y1 as left-truncation time)
        y2_quad <- c(y2, transform_quad_points(n_quad=n_quad,a=y1,b=y2,
                                                 quad_method=quad_method))
        basis3 <- get_basis(y=y2_quad,knots_vec=knots_list[[3]],hazard=hazard)
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
  if(tolower(hazard) %in% c("weibull","wb")){
    myLabels <- c("log(kappa1)", "log(alpha1)", "log(kappa2)",
                  "log(alpha2)", "log(kappa3)", "log(alpha3)")
  } else{
    myLabels <- c(paste0("phi1",1:p01),paste0("phi2",1:p02),paste0("phi3",1:p03))
  }
  nP <- c(ncol(Xmat1), ncol(Xmat2), ncol(Xmat3))
  nP0 <- c(p01,p02,p03)

  ##FIT MODEL##
  ##*********##
  #if the user has not provided start values, we generate them here
  #we start by fitting non-frailty model, so startVals is generated without frailty variance
  if(is.null(startVals)){
    startVals_nf <- get_start(y1=y1,y2=y2,delta1=delta1,delta2=delta2,yL=yL,anyLT=anyLT,
                              Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,weights=weights,
                              hazard=hazard,frailty=FALSE,model=model,knots_list=knots_list,
                              basis1=basis1,basis2=basis2,basis3=basis3)
    #In old code I also had a step where we initialized theta at coordinatewise max
    start_ltheta <- log(1)
    startVals <- c(startVals_nf[1:sum(nP0)],ltheta=start_ltheta,startVals_nf[-(1:sum(nP0))])
  } else{
    startVals_nf <- if(frailty) startVals[-(nP0+1)] else startVals
  }

  #check non-frailty gradient function
  grad1_nf <- ngrad_func(para = startVals_nf,y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                         Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                         hazard=hazard,model=model,frailty=FALSE,weights=weights,
                         basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                         basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                         dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
  grad2_nf <- pracma::grad(f = nll_func,x0 = startVals_nf,y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                           Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                           hazard=hazard,model=model,frailty=FALSE,weights=weights,
                           basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                           basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                           dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
  grad3_nf <- colSums(ngrad_mat_func(para = startVals_nf,y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                         Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                         hazard=hazard,model=model,frailty=FALSE,weights=weights,
                         basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                         basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                         dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3))
  # if(max(abs(grad1_nf-grad2_nf)) >= 1e-4){stop("check gradient of non-frailty model")}
  # if(max(abs(grad1_nf-grad3_nf)) >= 1e-4){stop("check gradient of non-frailty model")}
  # if(max(abs(grad1_nf-grad2_nf)) >= 1e-4){warning("check gradient of non-frailty model")}
  # cbind(grad1_nf,grad2_nf,grad3_nf)

  if(!frailty){
    #non-frailty model fit as three univariate models for the three transitions.
    value <- get_fit_nf(startVals_nf=startVals_nf,
                y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                hazard=hazard,model=model,
                basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                basis1_yL=basis1_yL, basis2_yL=basis2_yL,weights=weights,
                dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                control=con, hessian=out_options$Finv,optim_method=optim_method,
                extra_starts=extra_starts)
  } else{
    #Even for frailty model fit, optionally fit non-frailty model for comparison
    if(out_options$nf_fit){
      fit_nf <- get_fit_nf(startVals_nf=startVals_nf,
                           y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                           Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                           hazard=hazard,model=model,
                           basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                           basis1_yL=basis1_yL, basis2_yL=basis2_yL,weights=weights,
                           dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                           control=con, hessian=out_options$Finv,optim_method=optim_method,
                           extra_starts=extra_starts)
      #check gradient of individual non-frailty fits against gradient function
      grad3_nf <- ngrad_func(para = fit_nf$estimate,y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                             Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                             hazard=hazard,model=model,frailty=FALSE,
                             basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                             basis1_yL=basis1_yL, basis2_yL=basis2_yL,weights=weights,
                             dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
      grad4_nf <- pracma::grad(f = nll_func, x0 = fit_nf$estimate,y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                               Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                               hazard=hazard,model=model,frailty=FALSE,
                               basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                               basis1_yL=basis1_yL, basis2_yL=basis2_yL,weights=weights,
                               dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
      grad5_nf <- colSums(ngrad_mat_func(para = fit_nf$estimate,y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                                         Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                         hazard=hazard,model=model,frailty=FALSE,
                                         basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                         basis1_yL=basis1_yL, basis2_yL=basis2_yL,weights=weights,
                                         dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3))
      # if(max(abs(grad3_nf+fit_nf$grad)) >= 1e-4){stop("check gradient of non-frailty model")}
      # if(max(abs(grad3_nf-grad5_nf)) >= 1e-4){stop("check gradient of non-frailty model")}
      # if(max(abs(grad3_nf+fit_nf$grad)) >= 1e-4){warning("check gradient of non-frailty model")}
      # cbind(fit_nf$grad,grad3_nf,grad4_nf,grad5_nf)
    } else {
      nf_fit <- NULL
      frailty_lrtest <- c(stat=NA,pval=NA)
    }

    #now, update the basis matrices if we're
    #1. using a b-spline specification
    #2. have left truncation (and a frailty)
    if(hazard %in% c("bspline","bs") && anyLT && frailty){
      y1_quad <- c(y1,transform_quad_points(n_quad=n_quad,a=0,b=y1,
                                            quad_method=quad_method))
      basis1 <- get_basis(y=y1_quad, knots_vec=knots_list[[1]],hazard=hazard)
      basis2 <- get_basis(y=y1_quad, knots_vec=knots_list[[2]],hazard=hazard)
      yL_quad <- c(y1,transform_quad_points(n_quad=n_quad,a=0,b=yL,
                                            quad_method=quad_method))
      basis1_yL <- get_basis(y=yL_quad, knots_vec=knots_list[[1]],hazard=hazard)
      basis2_yL <- get_basis(y=yL_quad, knots_vec=knots_list[[2]],hazard=hazard)
      #lastly, we add attribute storing quadrature method
      attr(x=basis1,which="quad_method") <- attr(x=basis2,which="quad_method") <-
        attr(x=basis1_yL,which="quad_method") <- attr(x=basis2_yL,which="quad_method") <-
        tolower(quad_method)
    }

    #check frailty-based gradient functions
    grad1 <- ngrad_func(para = startVals,y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
               Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
               hazard=hazard,model=model,frailty=frailty,
               basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
               basis1_yL=basis1_yL, basis2_yL=basis2_yL,weights=weights,
               dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
    grad2 <- pracma::grad(f = nll_func,x0 = startVals,y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                 Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                 hazard=hazard,model=model,frailty=frailty,
                 basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                 basis1_yL=basis1_yL, basis2_yL=basis2_yL, weights=weights,
                 dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
    grad3 <- colSums(ngrad_mat_func(para = startVals,y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                        Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                        hazard=hazard,model=model,frailty=frailty,
                        basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                        basis1_yL=basis1_yL, basis2_yL=basis2_yL, weights=weights,
                        dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3))
    # if(max(abs(grad1-grad2)) >= 1e-4){stop("check gradient of frailty model")}
    # if(max(abs(grad1-grad3)) >= 1e-4){stop("check gradient of frailty model")}
    # if(max(abs(grad1-grad2)) >= 1e-4){warning("check gradient of non-frailty model")}
    # cbind(grad1,grad2,grad3)

    #Now, fit frailty model
    value <- get_fit_frail(startVals=startVals,
                y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                hazard=hazard,model=model,frailty=frailty,
                basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                basis1_yL=basis1_yL, basis2_yL=basis2_yL, weights=weights,
                dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                control=con, hessian=out_options$Finv,optim_method=optim_method,
                extra_starts=extra_starts)

    if(!is.null(value$fail)){
      return(list(fail=TRUE, formula=form2,
                  hazard=hazard,frailty=frailty,model=model,
                  startVals=startVals,knots_list=knots_list,
                  basis1=basis1,basis2=basis2,basis3=basis3,
                  basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                  dbasis1=dbasis1,dbasis2=dbasis2,dbasis3=dbasis3,
                  optim_method=optim_method,control=con,fit_nf=fit_nf))
    }

    if(out_options$nf_fit){
      #compute a likelihood-ratio test for the frailty
      #if the frailty likelihood is below the non-frailty likelihood, set to 0
      frail_lrtest_stat <- max(0,2*(value$logLike - fit_nf$logLike))
      if(!is.na(frail_lrtest_stat) & frail_lrtest_stat >= 0){
        frailty_lrtest <- c(stat=frail_lrtest_stat,
                            #corrected null distribution (mixture of chi-squareds)
                            pval= 0.5 * (stats::pchisq(q = frail_lrtest_stat,df = 0,lower.tail = FALSE) +
                                           stats::pchisq(q = frail_lrtest_stat,df = 1,lower.tail = FALSE)))
      } else{frailty_lrtest <- c(stat=NA,pval=NA)}
    }

    if(out_options$eb_frailties){
      eb_frailties <- pred_frailties(para = value$estimate,
                        hazard = hazard, model = model, knots_list = knots_list,
                        y1=y1, y2 = y2, delta1=delta1, delta2=delta2,
                        Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                        yL=yL, quad_method = quad_method, n_quad = n_quad)
    }

    myLabels <- c(myLabels, "log(theta)")
  }

  ##PREPARE OUTPUTS##
  ##***************##
  myLabels=c(myLabels,colnames(Xmat1),colnames(Xmat2),colnames(Xmat3))

  #if requested, compute the inverse hessian (aka sandwich "bread")
  if(out_options$Finv){
    Finv <- tryCatch(MASS::ginv(value$nhess),
            error=function(cnd){message(cnd);cat("\n");return(NA)})
  } else{ Finv <- NA }

  #if requested, compute gradient contributions for every subject
  if(out_options$grad_mat_return){
    grad_mat <- -ngrad_mat_func(para = value$estimate,
                 y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                 Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                 hazard=hazard,model=model,frailty=frailty,
                 basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                 basis1_yL=basis1_yL, basis2_yL=basis2_yL, weights=weights,
                 dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
  } else{ grad_mat <- NULL }

  #if requested, compute the sandwich variance "cheese"
  if(out_options$cheese){
    if(out_options$grad_mat_return){
      cheese <- crossprod(grad_mat)
    } else{
      cheese <- crossprod(ngrad_mat_func(para = value$estimate,
               y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
               Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
               hazard=hazard,model=model,frailty=frailty,
               basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
               basis1_yL=basis1_yL, basis2_yL=basis2_yL, weights=weights,
               dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3))
    }
  } else{ cheese <- NA }

  value <- list(
             estimate=as.vector(value$estimate),
             logLike=value$logLike,
             grad=as.vector(value$grad),
             optim_details=value$optim_details,
             Finv=Finv,
             cheese=cheese,
             eb_frailties=eb_frailties,
             startVals=value$startVals,
             knots_list=knots_list,
             myLabels=myLabels,
             formula=form2,nP=nP,nP0=nP0,nobs=length(y1),ymax=max(y2),n_quad=n_quad,
             quad_method=quad_method,optim_method=optim_method,
             extra_starts=extra_starts,control=con,
             grad_mat=grad_mat,
             data=if(out_options$data_return) data else NULL,
             frailty=frailty,
             frailty_lrtest = if(frailty) frailty_lrtest else NULL,
             fit_nf = if(frailty) fit_nf else NULL)

  names(value$estimate) <- if(frailty) names(startVals) else names(startVals_nf)
  names(value$grad) <- if(frailty) names(startVals) else names(startVals_nf)
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

