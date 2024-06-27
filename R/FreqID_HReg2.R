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
#'     \item{\code{Finv} Variance-covariance matrix. Defaults to \code{TRUE}.}
#'     \item{\code{fit_nf} If \code{frailty} is \code{TRUE}, include an
#'     additional fit object for the corresponding non-frailty model.
#'     Defaults to \code{TRUE}.}
#'     \item{\code{eb_frailties} Include vector of empirical Bayes
#'     predicted gamma frailty values. Defaults to \code{TRUE}.}
#'     \item{\code{grad_mat} Matrix with rowwise
#'     score vectors for each individual evaluated at the MLE.
#'     Used to compute "cheese" or "meat" in robust standard error computation.
#'     Defaults to \code{FALSE}.}
#'     \item{\code{cheese} Sum of outer products of individual score vectors,
#'     used as the "cheese" or "meat" in robust standard error computation.
#'     Defaults to \code{TRUE}.}
#'     \item{\code{data} Original data frame used to fit model.
#'     Defaults to \code{FALSE}.}
#'   }
#'
#' @return \code{FreqID_HReg2} returns an object of class \code{Freq_HReg}.
#' @import Formula
#' @importFrom pracma grad jacobian
#'
#' @examples
#' #loading a data set
#' data(scrData)
#'
#' #fitting Weibull semi-Markov illness-death model with gamma frailties
#' form <- Formula::Formula(time1 + event1 | time2 + event2 ~ x1 + x2 + x3 | x1 + x2 | x1 + x2)
#' fit_WB	<- FreqID_HReg2(Formula = form, data=scrData, model="semi-Markov",
#' extra_starts = 0,hazard = "weibull",frailty = TRUE,optim_method = c("BFGS"))
#'
#' #exploring results
#' fit_WB
#' summ.fit_WB <- summary(fit_WB); names(summ.fit_WB)
#' summ.fit_WB
#' pred_WB <- predict(fit_WB, tseq=seq(from=0, to=30, length.out=100))
#' plot(pred_WB, plot.est="Haz")
#' plot(pred_WB, plot.est="Surv")
#'
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

  #define list of components to include in model fit output list
  #if any output options have been given, override the defaults
  out_options=list(fit_nf=TRUE, Finv=TRUE, grad_mat=FALSE,
                   cheese=TRUE, eb_frailties=TRUE,
                   # beta_test=TRUE, #indicator for whether score and lr tests of each beta should also be performed, not currently used
                   data=FALSE)
  nmsO <- names(out_options)
  namO <- names(output_options)
  out_options[namO] <- output_options

  #initialize placeholder versions of possible objects to export
  eb_frailties <-
    # beta_test <- #for now, set aside the idea of adding score and lr tests of each beta
    grad_mat <-
    cheese <-
    fit_nf <- NULL
  frailty_lrtest <- c(frail_ll=NA,nonfrail_ll=NA,lrtest=NA,lrpvalue=NA)
  # frailty_scoretest <- c(sctest_onesided=NA,scpvalue_onesided=NA,
  #                        sctest_twosided=NA,scpvalue_twosided=NA)

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
  #if there are three pieces on the left side, the first is the left truncation variable
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

  #if weights have not been supplied, set all weights to 1.
  if(is.null(weights)){
    weights <- rep(1,length(y1))
  } else{
    stopifnot(length(weights) == length(y1))
  }


  #a few quick data checks
  if(any( (y2-y1)[delta1==1] == 0 )){
    warning("Some rows with nonterminal event (delta1==1) have sojourn time of 0 (y2-y1==0).\nThis will cause an error.")
  }
  if(any( (y2-y1)[delta1==0] != 0 )){
    warning("Some rows without the nonterminal event (delta1==0) have nonzero sojourn time (y2-y1 =/= 0).\nThis will cause incorrect estimation.")
  }



  ##PREPARE KNOTS AND BASIS FUNCTIONS FOR FLEXIBLE MODELS##
  ##*****************************************************##
  if(hazard %in% c("bspline","royston-parmar","piecewise","pw","rp","bs")){
    #if the knots (or changepoints) have not been provided, get the "defaults"
    if(is.null(knots_list)){
      if(length(nP0) != 3){
        stop("If hazard not equal to 'weibull' and knots_list set to NULL,
           then nP0 must be vector of 3 integers giving number of baseline parameters
           for each hazard submodel.")
      }
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
      basis1_yL <- basis2_yL <-
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
  #we start by fitting non-frailty model, so startVals_nf is generated without frailty variance
  if(is.null(startVals)){
    startVals_nf <- get_start(y1=y1,y2=y2,delta1=delta1,delta2=delta2,
                        yL=yL,anyLT=anyLT,
                        Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,weights=weights,
                        hazard=hazard,frailty=FALSE,model=model,
                        knots_list=knots_list,
                        basis1=basis1,basis2=basis2,basis3=basis3)

    #if the fitting fails, then return a bunch of information
    if(any(is.infinite(startVals_nf))){
      warning("Starting values generated by univariate transition model approximations were infinite. Try setting start values manually.")
      return(list(fail=TRUE, formula=form2,
                  hazard=hazard,frailty=frailty,model=model,
                  startVals_nf=startVals_nf,knots_list=knots_list,
                  basis1=basis1,basis2=basis2,basis3=basis3,
                  basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                  dbasis1=dbasis1,dbasis2=dbasis2,dbasis3=dbasis3,
                  optim_method=optim_method,control=con))
    }

    #In old code I also had a step where we initialized theta at coordinatewise max
    #but now, let's just start at 0 = log(1).
    start_ltheta <- log(1)
    startVals <- c(startVals_nf[1:sum(nP0)],ltheta=start_ltheta,startVals_nf[-(1:sum(nP0))])

  } else{
    startVals_nf <- if(frailty) startVals[-(nP0+1)] else startVals
  }


  #if the user requests a non-frailty model
  if(!frailty){
    #non-frailty model fit as three univariate models for the three transitions.
    value <- get_fit_nf(startVals=startVals_nf,
                y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                hazard=hazard,model=model,
                basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                basis1_yL=basis1_yL, basis2_yL=basis2_yL,weights=weights,
                dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                control=con, hessian=out_options$Finv,optim_method=optim_method,
                extra_starts=extra_starts)

    #if the fitting fails, then return a bunch of information
    if(!is.null(value$fail)){
      return(list(fail=TRUE, formula=form2,
                  hazard=hazard,frailty=frailty,model=model,
                  startVals_nf=startVals_nf,knots_list=knots_list,
                  basis1=basis1,basis2=basis2,basis3=basis3,
                  basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                  dbasis1=dbasis1,dbasis2=dbasis2,dbasis3=dbasis3,
                  optim_method=optim_method,control=con))
    }

  } else{
    #Even for frailty model fit, optionally fit non-frailty model for comparison
    if(out_options$fit_nf){
      fit_nf <- get_fit_nf(startVals=startVals_nf,
                 y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                 Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                 hazard=hazard,model=model,
                 basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                 basis1_yL=basis1_yL, basis2_yL=basis2_yL,weights=weights,
                 dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                 control=con, hessian=out_options$Finv,optim_method=optim_method,
                 extra_starts=extra_starts)
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

    if(out_options$fit_nf){
      #compute a likelihood-ratio test for the frailty
      #if the frailty likelihood is below the non-frailty likelihood, set to 0
      frail_lrtest_stat <- max(0,2*(value$logLike - fit_nf$logLike))
      if(!is.na(frail_lrtest_stat) & frail_lrtest_stat >= 0){
        frailty_lrtest <- c(frail_ll=value$logLike,
                            nonfrail_ll=fit_nf$logLike,
                            lrtest=frail_lrtest_stat,
                            #corrected null distribution (mixture of chi-squareds)
                            lrpvalue= 0.5 * (stats::pchisq(q = frail_lrtest_stat,df = 0,lower.tail = FALSE) +
                                           stats::pchisq(q = frail_lrtest_stat,df = 1,lower.tail = FALSE)))
      }


      #FOR NOW, I WILL COMMENT OUT THE POSSIBLE SCORE TESTING OF THE FRAILTY!


      # #here's what I think it would take to have a score test of the frailty:
      # #compute the score equation (and its jacobian)
      # #from the theta-parameterized likelihood using non-frailty ests
      # #and setting theta=0. Compute score test and compare to mixture null dist.
      #
      # #Technically, this looks like a two-sided test, but I think there's a deeper connection:
      # #if you imagined a statistic normally-distributed under the null (e.g., theta=0)
      # #then the two-sided statistic looks at theta =/= 0. Put another way, the
      # #statistic squared would follow a chi-squared with 1 df.
      # #But, if you first truncated that statistic at 0, and then squared it,
      # #the distribution would be different! It would be exactly a 50:50 mixture
      # #between a point mass at 0, and a chi-squared with 1 df.
      # #truncating the statistic is akin to "constraining" the alternative space?
      # if(hazard=="weibull"){
      #   #brute forcing numerical differentiation to get gradient and hessian under null
      #   ngrad_ID_theta <- function(x0){ pracma::grad(x0 = x0, f = nlogLikWB_ID_theta,
      #                y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
      #                X1=Xmat1, X2=Xmat2, X3=Xmat3, model=model, frailty_ind=1, weights=weights)}
      #   null_mle_para <- append(fit_nf$estimate,0,6)
      #   score_ngrad_temp <- ngrad_ID_theta(null_mle_para)
      #   score_njac_temp <- pracma::jacobian(x0 = null_mle_para, f = ngrad_ID_theta)
      #
      #   # Test defined by Silvapulle and Silvapulle, and clarified by Claeskens 2008
      #   # my notion is that this would be normally distributed under the null
      #   # (note the negative sign)
      #   frailty_scoretest_norm <- -score_ngrad_temp[7] * sqrt(MASS::ginv(score_njac_temp)[7,7])
      #   frailty_scoretest_twosided <- frailty_scoretest_norm^2
      #   # now, after truncating and squaring this will have "mixture" distribution of LRT above.
      #   frailty_scoretest_onesided <- max(0, frailty_scoretest_norm)^2
      #
      #   #separately, following discussion of Robertson (1988) p. 321,
      #   #we might actually be inclined to use a "score test" approach to
      #   #determine if the frailty is strictly positive, or not.
      #   temp_para <- value$estimate
      #   temp_para[7] <- exp(value$estimate[7])
      #   score_ngrad_alt <- ngrad_ID_theta(temp_para)
      #   score_njac_alt <- pracma::jacobian(x0 = temp_para, f = ngrad_ID_theta)
      #   sctest_negative <- (score_ngrad_alt[7])^2 * MASS::ginv(score_njac_alt)[7,7]
      #
      #   # Finally, another test defined by Robertson (1988) p. 321 is a score test
      #   # for theta=0 vs theta>=0
      #   # sctest_robertson <- t(score_ngrad_alt - score_ngrad_temp) %*%
      #   #   MASS::ginv(score_njac_alt) %*% (score_ngrad_alt - score_ngrad_temp)
      #   # sctest_robertson <- (score_ngrad_alt[7] - score_ngrad_temp[7])^2 * MASS::ginv(score_njac_alt)[7,7]
      #
      #   frailty_scoretest <- c(sctest_onesided=frailty_scoretest_onesided,
      #                          scpvalue_onesided= 0.5 *
      #                            (stats::pchisq(q = frailty_scoretest_onesided,df = 0,lower.tail = FALSE) +
      #                            stats::pchisq(q = frailty_scoretest_onesided,df = 1,lower.tail = FALSE)),
      #                          sctest_twosided=frailty_scoretest_twosided,
      #                          scpvalue_twosided=stats::pchisq(q = frailty_scoretest_twosided,df = 1,lower.tail = FALSE))
      # }

    }

    #if requested, empirical bayes estimation of the gamma frailties
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

  #if requested, compute gradient contributions for every subject
  if(out_options$grad_mat){
    grad_mat <- -ngrad_mat_func(para = value$estimate,
                 y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                 Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                 hazard=hazard,model=model,frailty=frailty,
                 basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                 basis1_yL=basis1_yL, basis2_yL=basis2_yL, weights=weights,
                 dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
  }

  #if requested, compute the sandwich variance "cheese" outer product of scores
  if(out_options$cheese){
    if(out_options$grad_mat){
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
  }


  #for now, comment out the score and likelihood ratio tests for every beta

  # #if requested, compute score and likelihood ratio tests for every beta
  # if(out_options$beta_test && sum(nP) > 0){
  #   para_mle <- as.vector(value$estimate)
  #   beta_test <- cbind(beta=para_mle[-(1:(sum(nP0) + as.numeric(frailty)))],
  #                        lrtest=NA, sctest=NA,
  #                      lrpvalue=NA,scpvalue=NA)
  #   for(i in 1:sum(nP)){
  #     #compute "null" negative log-likelihood fixing beta=0
  #     temp_prof <- nll_profile_func(fixed_param=0,
  #                      fixed_param_ind=sum(nP0)+i+as.numeric(frailty)-1,
  #                      para_mle=para_mle,
  #                      y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
  #                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
  #                      hazard=hazard, frailty=frailty, model=model, weights=weights,
  #                      basis1=basis1, basis2=basis2, basis3=basis3,
  #                      basis3_y1=basis3_y1, basis1_yL=basis1_yL, basis2_yL=basis2_yL,
  #                      dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
  #                      verbose=FALSE, control=control,
  #                      optim_method=optim_method, extra_starts=extra_starts)
  #     #compute likelihood ratio test statistic
  #     beta_test[i,"lrtest"] <- 2 * (value$logLike + temp_prof$nll)
  #     beta_test[i,"lrpvalue"] <- pchisq(q=beta_test[i,"lrtest"],df=1,lower.tail = FALSE)
  #     #implement score test of beta
  #     #compute gradient at MLE "under null" of beta=0
  #     score_ngrad <- ngrad_func(para = temp_prof$paramat[1,],
  #                    y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
  #                    Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
  #                    hazard=hazard, frailty=frailty, model=model, weights=weights,
  #                    basis1=basis1, basis2=basis2, basis3=basis3,
  #                    basis3_y1=basis3_y1, basis1_yL=basis1_yL, basis2_yL=basis2_yL,
  #                    dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
  #     #compute numerical jacobian of above gradient
  #     score_nhess <- pracma::jacobian(f = ngrad_func, x0 = temp_prof$paramat[1,],
  #                    y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
  #                    Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
  #                    hazard=hazard, frailty=frailty, model=model, weights=weights,
  #                    basis1=basis1, basis2=basis2, basis3=basis3,
  #                    basis3_y1=basis3_y1, basis1_yL=basis1_yL, basis2_yL=basis2_yL,
  #                    dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
  #     #compute score test statistic
  #     beta_test[i,"sctest"] <- score_ngrad[sum(nP0)+i+as.numeric(frailty)]^2 *
  #       MASS::ginv(score_nhess)[sum(nP0)+i+as.numeric(frailty),sum(nP0)+i+as.numeric(frailty)]
  #     # beta_test[i,"sctest"] <- t(score_ngrad) %*% MASS::ginv(score_nhess) %*% score_ngrad
  #     beta_test[i,"scpvalue"] <- pchisq(q=beta_test[i,"sctest"],df=1,lower.tail = FALSE)
  #   }
  # }

  class_temp <- c("Freq_HReg2","ID","Ind",
                   switch(tolower(hazard),
                          weibull="Weibull",wb="Weibull",
                          bspline="B-Spline",bs="B-Spline",
                          "royston-parmar"="Royston-Parmar",rp="Royston-Parmar",
                          piecewise="Piecewise Constant",pw="Piecewise Constant"),
                   switch(tolower(model),
                          "semi-markov"="semi-Markov",
                          "markov"="Markov"))

  value <- list(
             estimate=as.vector(value$estimate),
             logLike=value$logLike,
             grad=as.vector(value$grad),
             Finv=value$Finv, #automatically set to NULL above if out_list$hessian is FALSE
             cheese=cheese, grad_mat=grad_mat,
             eb_frailties=eb_frailties,
             startVals=value$startVals,
             knots_list=knots_list,
             myLabels=myLabels,
             formula=form2, nP=nP, nP0=nP0, nobs=length(y1), ymax=max(y2),
             optim_details=value$optim_details,
             optim_method=optim_method, extra_starts=extra_starts, control=con,
             n_quad=n_quad, quad_method=quad_method,
             data=if(out_options$data) data else NULL,
             hazard=hazard, model=model, frailty=frailty,
             frailty_test=frailty_lrtest,
             # frailty_test=c(frailty_lrtest,frailty_scoretest),
             # beta_test=beta_test,
             class=class_temp,
             fit_nf=fit_nf)

  names(value$estimate) <- if(frailty) names(startVals) else names(startVals_nf)
  names(value$grad) <- if(frailty) names(startVals) else names(startVals_nf)
  class(value) <- "Freq_HReg2"
  return(value)
}

