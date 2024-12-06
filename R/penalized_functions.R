##********************************************************##
####Main penalized illness-death model fitting functions####
##********************************************************##


#' Fit Penalized Parametric Frailty Illness-Death Model Solution Path
#'
#' @inheritParams solution_path_function
#' @inheritParams FreqID_HReg2
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
#' @param startVals A numeric vector of parameter starting values, arranged as follows:
#'   the first \eqn{k_1+k_2+k_3} elements correspond to the baseline hazard parameters,
#'   then the \eqn{k_1+k_2+k_3+1} element corresponds to the gamma frailty log-variance parameter,
#'   then the last\eqn{q_1+q_2+q_3} elements correspond with the regression parameters.
#'   If set to \code{NULL}, will be generated automatically using \code{\link{get_start}}.
#' @param lambda_target Final lambda value for parameterwise penalty. Ignored if lambda_path is specified.
#' @param N_path_steps number of steps generated for regularization path. Ignored if lambda_path is specified.
#' @param standardize Boolean for whether columns of covariate matrix should be standardized before estimation.
#' Estimates are then rescaled back to original scale when output.
#' @param fusion_tol Positive numeric value for thresholding estimates that are close
#'   to being considered fused, for the purposes of estimating degrees of freedom.
#'
#' @return A list.
#'
#' @export
FreqID_HReg_Rpath <- function(Formula, data, na.action="na.fail", subset=NULL, weights=NULL,
                              hazard=c("weibull"),frailty=TRUE,model, knots_list = NULL,
                              penalty=c("scad","mcp","lasso"),
                              lambda_path=NULL, lambda_target=0, N_path_steps = 40, #passing in lambda_path overrides automatic calculation of path
                              a=NULL, select_tol=1e-4, fusion_tol=1e-3,
                              penalty_fusedcoef=c("none","fusedlasso"), lambda_fusedcoef_path=0,
                              penweights_list=list(), mu_smooth_path=0,
                              fit_method="prox_grad",
                              nP0=rep(4,3), startVals=NULL,
                              quad_method="kronrod", n_quad=15,
                              ball_L2=Inf,
                              warm_start=TRUE, step_size_min=1e-6, step_size_max=1e6, step_size_init=1,
                              step_size_scale=0.5, #no checks implemented on these values!!
                              maxit=300, extra_starts=0,
                              conv_crit = "nll_pen_change", conv_tol=1e-6, standardize=TRUE,
                              verbose=0){

  # To start, I'm going to implement the PISTA algorithm of Wang et al (2014).  I know it is somewhat deficient compared to
  # the APISTA and PICASSO algorithms subsequently proposed by the same authors, but it's a good start with what I have previously implemented
  # I also think that relatively speaking, the cost of the nll and ngrad functions doesn't nicely decompose for coordinate methods in the same way? idk.
  # Gotta start somewhere, that's all.

  ##INITIALIZE OPTIONS##
  ##******************##
  #checks on the penalties and on 'a' happen in the underlying functions to minimize duplication

  na.action <- match.arg(na.action) #technically na.fail I think is the only one currently implemented
  if (na.action != "na.fail" & na.action != "na.omit") {
    stop("na.action should be either na.fail or na.omit")
  }

  #Check that chosen hazard is among available options
  hazard_options <- c("weibull","wb","piecewise","pw","royston-parmar","rp","bspline","bs")
  if(!(tolower(hazard) %in% hazard_options)){
    stop(paste0("valid choices of hazard are '", paste(hazard_options,collapse = "', '"),"'"))
  } else{ hazard <- tolower(hazard) }

  if(!(penalty %in% c("scad","mcp","lasso"))){
    stop("penalty must be either 'scad', 'mcp', 'lasso'")
  }
  if(!(penalty_fusedcoef %in% c("none","fusedlasso"))){
    stop("penalty_fusedcoef must be 'none' or 'fusedlasso'")
  }

  if (a <= 0)
    stop("a must be greater than 0")
  if (a <= 1 & penalty == "mcp")
    stop("a must be greater than 1 for the MC penalty")
  if (a <= 2 & penalty == "scad")
    stop("a must be greater than 2 for the SCAD penalty")

  # browser()

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
  n <- length(y1)

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

  if(standardize){
    Xmat1 <- scale(Xmat1)
    scaling_center1 <- attr(Xmat1,"scaled:center")
    scaling_factor1 <- attr(Xmat1,"scaled:scale")
    Xmat2 <- scale(Xmat2)
    scaling_center2 <- attr(Xmat2,"scaled:center")
    scaling_factor2 <- attr(Xmat2,"scaled:scale")
    Xmat3 <- scale(Xmat3)
    scaling_center3 <- attr(Xmat3,"scaled:center")
    scaling_factor3 <- attr(Xmat3,"scaled:scale")
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
  nP0tot <- sum(nP0, frailty)
  nPtot <- sum(nP0, nP, frailty)

  ##Set up algorithmic parameters##
  ##*****************************##

  #if the user has not provided start values, we generate them here
  #we start by fitting non-frailty model, so startVals_nf is generated without frailty variance
  if(is.null(startVals)){
    startVals_nf <- get_start(y1=y1,y2=y2,delta1=delta1,delta2=delta2,
                              yL=yL,anyLT=anyLT,
                              Xmat1=Xmat1,Xmat2=Xmat2,Xmat3=Xmat3,weights=weights,
                              hazard=hazard,frailty=FALSE,model=model,
                              knots_list=knots_list,
                              basis1=basis1,basis2=basis2,basis3=basis3, sparse_start = TRUE)

    #if the fitting fails, then return a bunch of information
    if(any(is.infinite(startVals_nf))){
      warning("Starting values generated by univariate transition model approximations were infinite. Try setting start values manually.")
      return(list(fail=TRUE, formula=form2,
                  hazard=hazard,frailty=frailty,model=model,
                  startVals_nf=startVals_nf,knots_list=knots_list,
                  basis1=basis1,basis2=basis2,basis3=basis3,
                  basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                  dbasis1=dbasis1,dbasis2=dbasis2,dbasis3=dbasis3))
    }

    #In old code I also had a step where we initialized theta at coordinatewise max
    #but now, let's just start at 0 = log(1).
    start_ltheta <- log(1)
    if(frailty){
      startVals <- c(startVals_nf[1:sum(nP0)],ltheta=start_ltheta,startVals_nf[-(1:sum(nP0))])
    } else{
      startVals <- startVals_nf
    }

  }

  if(length(startVals) != nPtot){
    stop("length of startVals vector provided does not match dimensionality.")
  }

  # nPtot <- length(startVals)
  # nP1 <- ncol(Xmat1)
  # nP2 <- ncol(Xmat2)
  # nP3 <- ncol(Xmat3)
  # nP0 <- nPtot - nP1 - nP2 - nP3
  # nP_vec <- c(nP0,nP1,nP2,nP3)

  #if lambda_path is null then compute it according to Wang (2014), otherwise go with what was provided
  if(is.null(lambda_path)){
    startVals_grad <- ngrad_func(para=startVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                 yL=yL, anyLT=anyLT,
                                 Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                 hazard=hazard, frailty=frailty, model=model, weights=weights,
                                 basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                 basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                                 dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)/n #notice the normalization by n

    # N_path_steps <- 40
    lambda0 <- 1.1 * max(abs(startVals_grad[-(1:nP0tot)])) #a tad bit larger than the largest gradient of the betas
    #using Wang (2014) 3.3, solving for 'eta' given the other three pieces
    #add the -1 to account for starting at 0, so that final results are of length N_path_steps
    grid_eta <- exp((log(lambda_target) - log(lambda0)) / (N_path_steps-1))
    if(grid_eta < 0.9 | grid_eta > 1){
      print(c(N_path_steps=N_path_steps,grid_eta=grid_eta,lambda_target=lambda_target,lambda0=lambda0))
      warning("provided lambda_target is too big, or provided N_path_steps is too small. recommend grid_eta > 0.9")
    }
    if(lambda_target > lambda0){
      print(c(N_path_steps=N_path_steps,grid_eta=grid_eta,lambda_target=lambda_target,lambda0=lambda0))
      warning("provided lambda_target is bigger than largest gradient, results may be overly regularized.")
    }
    #for now, we're just gonna let the same lambda govern all the betas, we'll work on the rest later.
    #note that we already 'have the solution' for the max value, because it is the startVal, FIX
    lambda_path_vec <- lambda0 * grid_eta^(0:(N_path_steps-1)) #
    lambda_path <- cbind(lambda_path_vec,lambda_path_vec,lambda_path_vec)
  } else{
    grid_eta <- NULL
    if(NCOL(lambda_path)==1){
      lambda_path <- cbind(lambda_path,lambda_path,lambda_path)
    }
  }

  colnames(lambda_path) <- paste0("lambda",1:ncol(lambda_path))

  ##ACTUALLY BEGIN THE PATH##
  ##***********************##

  # browser()

  solution_path_out <- solution_path_function(para=startVals, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                              yL=yL, anyLT=anyLT,
                                              Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                              hazard=hazard, frailty=frailty, model=model, weights=weights,
                                              basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                              basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                                              dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                              penalty=penalty, lambda_path=lambda_path, a=a,
                                              penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef_path=lambda_fusedcoef_path,
                                              #penalty_fusedbaseline="none", lambda_fusedbaseline=0, #idk what to do about the baseline fusing, for now I skip
                                              penweights_list=penweights_list, mu_smooth_path=mu_smooth_path, ball_L2=ball_L2,
                                              fit_method=fit_method, warm_start=warm_start,
                                              step_size_min=step_size_min, step_size_max=step_size_max, step_size_init=step_size_init,
                                              step_size_scale=step_size_scale, #no checks implemented on these values!!
                                              maxit=maxit,
                                              conv_crit = conv_crit,
                                              conv_tol=conv_tol,
                                              extra_starts=extra_starts,
                                              verbose=verbose)

  solution_path_out$lambda_path <- lambda_path
  solution_path_out$lambda_fusedcoef_path <- lambda_fusedcoef_path
  solution_path_out$penweights_list <- penweights_list

  solution_path_out$grid_eta <- grid_eta
  solution_path_out$knots_list <- knots_list
  solution_path_out$ball_L2 <- ball_L2
  solution_path_out$nP <- nP
  solution_path_out$nP0 <- nP0
  solution_path_out$detail <- list(
    quad_method=quad_method, n_quad=n_quad,
    select_tol=select_tol, fusion_tol=fusion_tol, fit_method=fit_method,
    mu_smooth_path=mu_smooth_path,
    warm_start=warm_start, step_size_scale=step_size_scale,
    maxit=maxit, extra_starts=extra_starts,
    conv_crit = conv_crit, conv_tol=conv_tol, standardize=standardize,
    verbose=verbose)

  if(standardize){
    solution_path_out$ests[,-(1:nP0tot)] <-
      sweep(x = solution_path_out$ests[,-(1:nP0tot)],
            STATS = c(scaling_factor1,scaling_factor2,scaling_factor3),
            MARGIN = 2, FUN = "/")
  }

  return(solution_path_out)
}






#' Estimate Penalized Illness-Death Model Solution Path
#'
#' This function estimates penalized illness-death model results along a range of
#'   penalty, fused penalty, and smoothing parameters.
#'
#' This is a function to loop through a path of lambda, lambda_fusedcoef, and mu_smooth values in a somewhat
#'   thoughtful way to maximize the pathwise connections between starting values and step sizes,
#'   with some adjustments tailored to each approach to optimization.
#' @inheritParams proximal_gradient_descent
#' @param lambda_path Numeric sequence of decreasing regularization parameters
#'   for the parameterwise penalties, along which the solution path runs.
#'   Assumes a single shared penalty across transitions.
#' @param lambda_fusedcoef_path Numeric sequence of increasing regularization parameters
#'   for the fusion penalties.
#' @param mu_smooth_path Numeric sequence of decreasing Nesterov smoothing parameters for the fusion penalties.
#' @param fit_method String indicating which optimization method should be used at each step.
#' @param warm_start Boolean indicating whether each step of the solution should start from
#'   ending of previous (\code{TRUE}) or from the original starting point (\code{FALSE}).
#' @param extra_starts numeric indicating how many additional optimization runs from random start values
#'   should be performed at each grid point.
#' @param fusion_tol Numeric value indicating when to consider fused parameters
#'   that are close to be considered the same value, for estimating degrees of freedom.
#'
#' @return A list.
#' @export
solution_path_function <- function(para, y1, y2, delta1, delta2, yL, anyLT,
                                   Xmat1 = matrix(nrow(length(y1)),ncol=0), #the default is a 'empty' matrix
                                   Xmat2 = matrix(nrow(length(y1)),ncol=0),
                                   Xmat3 = matrix(nrow(length(y1)),ncol=0),
                                   hazard, frailty, model, weights,
                                   basis1=NULL, basis2=NULL, basis3=NULL, basis3_y1=NULL,
                                   basis1_yL=NULL, basis2_yL=NULL,
                                   dbasis1=NULL, dbasis2=NULL, dbasis3=NULL,
                                   penalty, lambda_path, a,
                                   penalty_fusedcoef, lambda_fusedcoef_path,
                                   #penalty_fusedbaseline="none", lambda_fusedbaseline=0, #idk what to do about the baseline fusing, for now I skip
                                   penweights_list, mu_smooth_path, ball_L2=Inf,
                                   fit_method, warm_start=TRUE, extra_starts=0,
                                   select_tol=1e-4, fusion_tol = 1e-3,
                                   step_size_min=1e-6, step_size_max=1e6,
                                   step_size_init=1,
                                   step_size_scale=0.5, #no checks implemented on these values!!
                                   maxit=300,
                                   conv_crit = "nll_pen_change",
                                   conv_tol=1e-6,
                                   verbose){

  # browser()

  #If lambda_path is NULL, set it equal to 0 (aka no parameterwise regularization)
  #Otherwise, put it into matrix format
  # (note, this admits a 3-column matrix to let different transitions have different values)
  lambda_path <- if(is.null(lambda_path)) as.matrix(0) else as.matrix(lambda_path)
  lambda_length <- NROW(lambda_path)
  if(NCOL(lambda_path) == 1){
    lambda_path <- cbind(lambda_path,lambda_path,lambda_path)
  }
  colnames(lambda_path) <- paste0("lambda",1:NCOL(lambda_path))

  #If lambda_fusedcoef_path is NULL, set it equal to 0 (aka no fused regularization)
  #Otherwise, put it into matrix format
  # (note, this admits a 3-column matrix to let different transitions have different values)
  lambda_fusedcoef_path <- if(is.null(lambda_fusedcoef_path)) as.matrix(0) else as.matrix(lambda_fusedcoef_path)
  lambda_fusedcoef_length <- NROW(lambda_fusedcoef_path)
  if(NCOL(lambda_path) == 1){
    lambda_fusedcoef_path <- cbind(lambda_fusedcoef_path,lambda_fusedcoef_path,lambda_fusedcoef_path)
  }
  colnames(lambda_fusedcoef_path) <- paste0("lambda_fusedcoef",1:ncol(lambda_fusedcoef_path))

  #If mu_smooth_path is NULL, set it equal to 0 (aka no fused smoothing)
  #Otherwise, put it into matrix format
  mu_smooth_path <- if(is.null(mu_smooth_path)) as.matrix(0) else as.matrix(mu_smooth_path)
  mu_smooth_length <- nrow(mu_smooth_path)
  colnames(mu_smooth_path) <- "mu_smooth" #for now, implicit that there is only one column


  n <- length(y1)
  Xmat1 <- if(!is.null(Xmat1)) as.matrix(Xmat1) else matrix(nrow=n,ncol=0)
  Xmat2 <- if(!is.null(Xmat2)) as.matrix(Xmat2) else matrix(nrow=n,ncol=0)
  Xmat3 <- if(!is.null(Xmat3)) as.matrix(Xmat3) else matrix(nrow=n,ncol=0)

  nPtot <- length(para)
  nP1 <- ncol(Xmat1)
  nP2 <- ncol(Xmat2)
  nP3 <- ncol(Xmat3)
  nP0 <- nPtot - nP1 - nP2 - nP3
  nP_vec <- c(nP0,nP1,nP2,nP3)


  #to correctly count number of iterations, we sum up number of iterations with no
  # total_length <- lambda_length*sum(lambda_fusedcoef_path %in% 0) + #total number of steps with no fusion
  #   lambda_length*sum(!(lambda_fusedcoef_path %in% 0)) * mu_smooth_length #total number of fusion steps

  #in this new version, we're not storing every single step of the mu_smoothing loop, just the final one
  #therefore, we don't need to play around with how many iterations are smoothed vs not smoothed.
  #we just get a single output, regardless of whether there is any fusion or not.
  total_length <- lambda_length*lambda_fusedcoef_length

  out_starts <- out_pen_ngrads <- out_ests <- matrix(nrow=total_length,
                                                     ncol=nPtot,dimnames=list(NULL,names(para)))
  out_ics <- matrix(nrow=total_length, ncol=11,
                    dimnames=list(NULL,c("nll", "nll_pen",
                                         "nPtot", "nPtot_selected",
                                         "df_unique",
                                         "AIC", "AIC_unique",
                                         "BIC", "BIC_unique",
                                         "GCV", "GCV_unique")))
  out_info <- matrix(nrow=total_length,
                     ncol=ncol(lambda_path) + ncol(lambda_fusedcoef_path) + ncol(mu_smooth_path) + 1,
                     dimnames=list(NULL,c(colnames(lambda_path),
                                          colnames(lambda_fusedcoef_path),
                                          colnames(mu_smooth_path),"ball_L2")))
  out_conv_stats <- matrix(nrow=total_length,
                        ncol=9,dimnames=list(NULL,c("niter","fit_code","final_step_size","bad_step_count",
                                                    "best_start_iter",
                                                    "est_change_max","est_change_2norm",
                                                    "nll_pen_change","subopt")))

  nll_pen_trace_mat <-  matrix(nrow=total_length,
                               ncol=maxit)

  iter <- 1
  ##ACTUALLY BEGIN THE PATH##
  ##***********************##

  # browser()

  #keep running pathwise starting values and step sizes for the lambda, lambda_fused, and mu_smooth loops
  startVals_outer <- startVals_middle <- startVals_inner <- para
  step_size_init_outer <- step_size_init_middle <- step_size_init_inner <- step_size_init

  #OUTER LOOP: LOOPING THROUGH THE LAMBDA VALUES, GOVERNING PARAMETERWISE SPARSITY
  for(lambda_iter in 1:lambda_length){
    lambda <- lambda_path[lambda_iter,]

    #Wang (2014) paper advises specific weaker tolerance earlier in the path, according to lambda value:
    # conv_tol <- if(lambda_iter==lambda_length) lambda/8 else lambda/4 #gotta be less than lambda_target/4, so we just halve it again.

    startVals_middle <- startVals_outer
    step_size_init_middle <- step_size_init_outer

    #MIDDLE LOOP: LOOPING THROUGH THE LAMBDA_FUSEDCOEF VALUES, GOVERNING FUSION SPARSITY
    for(lambda_fusedcoef_iter in 1:lambda_fusedcoef_length){

      lambda_fusedcoef <- lambda_fusedcoef_path[lambda_fusedcoef_iter,]
      startVals_inner <- startVals_middle
      step_size_init_inner <- step_size_init_middle

      if(verbose >= 1){
        print(paste0("step: ",lambda_iter,
                     " lambda: ",paste(round(lambda,6),collapse = " "),
                     " lambda_fusedcoef: ",paste(round(lambda_fusedcoef,6),collapse = " ")))#, " mu_smooth_fused: ",mu_smooth_fused, " extra_start: ", extra_start_iter))
      }

      #monitor to track which random start achieves the lowest objective value
      best_nll_pen <- Inf

      #INNER LOOP: LOOPING THROUGH EXTRA START VALUES
      for(extra_start_iter in 1:(extra_starts+1)){

        #if this is the first time through, admit the warm start approach, otherwise randomize the start value
        #and reset the step size to the global default
        if(extra_start_iter == 1){
          startVals_innermost <- startVals_inner
          step_size_init_innermost <- step_size_init_inner
        } else{
          #Here I'm basing my perturbed start values based off of the 'middle' starting values
          startVals_innermost <- (startVals_inner + stats::rnorm(n = length(startVals_inner),mean = 0,sd=0.7)) *
            stats::runif(n = length(startVals_inner),min = 0.9, max = 1.1) #adding random multiplicative and additive noise
          step_size_init_innermost <- step_size_init
        }

        #INNERMOST LOOP: LOOPING THROUGH MU_SMOOTH VALUES, GOVERNING SMOOTH APPROXIMATION OF FUSION SPARSITY
        for(mu_smooth_iter in 1:mu_smooth_length){

          #if there is no fusion, there should be no smoothing
          mu_smooth_fused <- if(all(lambda_fusedcoef==0)) 0 else mu_smooth_path[mu_smooth_iter,]

          if(fit_method=="prox_grad"){
            tempfit <- proximal_gradient_descent(para=startVals_innermost,
                                                 y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                                 yL=yL, anyLT=anyLT,
                                                 Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                                 hazard=hazard, frailty=frailty, model=model, weights=weights,
                                                 basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                                 basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                                                 dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                                 penalty=penalty, lambda=lambda, a=a,
                                                 penalty_fusedcoef=penalty_fusedcoef, lambda_fusedcoef=lambda_fusedcoef,
                                                 #penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
                                                 penweights_list=penweights_list, mu_smooth_fused=mu_smooth_fused,
                                                 step_size_init=step_size_init_innermost,
                                                 step_size_min=step_size_min,step_size_max=step_size_max,
                                                 step_size_scale=step_size_scale, select_tol=select_tol,
                                                 maxit=maxit, conv_crit=conv_crit, conv_tol=conv_tol,
                                                 ball_L2=ball_L2, verbose=verbose)
          }

          #now, in this innermost loop we really just step through the iterations,
          #without really tracking the intermediate smoothing results
          startVals_innermost <- tempfit$estimate

          # Wang (2014) suggest in their pathwise approach to carry over the step size from each iteration
          # So we could do that here, but it's currently commented out, and even then
          # only written for the proximal algorithm.
          # The accelerated method re-estimates its step size
          # at each step anyways so the initial step being too big would be less of an issue.
          # if(fit_method %in%  c("prox_grad")){
          #   step_size_init_innermost <- tempfit$final_step_size
          # }

          #if there is no smoothing, then break out of the smoothing (innermost) loop
          #recall that above we set this to 0 if theres no fusion at all.
          if(mu_smooth_fused==0){
            break
          }
        } #END OF INNERMOST (SMOOTHING) LOOP

        #Now, having run the innermost loop though all of the mu_smooth steps (if there's only one, just the one)
        #We assess whether it has reached a better spot than the previous start

        #If this random start has reached a better penalized value  (this is a completely unsmoothed value, fyi)
        #than the previous start, then update all of the resulting outputs
        if(tempfit$final_nll_pen < best_nll_pen){
          if(verbose >= 1 && extra_start_iter > 1){
            print(paste0("start ",extra_start_iter,
                         "yielded better obj value ", tempfit$final_nll_pen,
                         " over ",best_nll_pen))
          }
          best_nll_pen <- tempfit$final_nll_pen

          ##Now, let's COMPUTE SOME USEFUL FIT STATISTICS based on this fit##
          ##***************************************************************##
          final_nll <- tempfit$final_nll

          beta1_selected <- if(nP1 != 0) tempfit$estimate[(1+nP0):(nP0+nP1)] else numeric(0)
          nP1_selected <- sum(beta1_selected != 0)
          beta2_selected <- if(nP2 != 0) tempfit$estimate[(1+nP0+nP1):(nP0+nP1+nP2)] else numeric(0)
          nP2_selected <- sum(beta2_selected != 0)
          beta3_selected <- if(nP3 != 0) tempfit$estimate[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)] else numeric(0)
          nP3_selected <- sum(beta3_selected != 0)

          nPtot_selected <- nP0 + nP1_selected + nP2_selected + nP3_selected

          tempaic <- 2*final_nll + 2* nPtot_selected
          tempbic <- 2*final_nll + log(n) * nPtot_selected
          tempgcv <- final_nll/(n^2*(1-nPtot_selected/n)^2)

          # The provisional idea for degrees of freedom is to take the number of unique nonzero estimates, e.g.,
          # nPtot_selected_unique <- sum(unique(finalVals_selected) != 0)
          #but here we use a tolerance because we can't get exact equality
          #NOTE THIS DOES NOT INCORPORATE POTENTIAL EQUALITY OF BASELINE PARAMETERS
          #ALSO IT ASSUMES THE SAME COVARIATES IN THE SAME ORDER ACROSS ALL THREE HAZARDSSSS!!!
          if(nP1==nP2 & nP2==nP3){
            df_unique <- nP0 + sum(beta1_selected != 0 &
                                     abs(beta1_selected - beta2_selected) > fusion_tol &
                                     abs(beta1_selected - beta3_selected) > fusion_tol) +
              sum(beta2_selected != 0 &
                    abs(beta2_selected - beta3_selected) > fusion_tol) +
              sum(beta3_selected != 0)

            #currently, these account for fusion in the 'degrees of freedom' by counting unique parameters
            tempaic_unique <- 2*final_nll + 2* df_unique
            tempbic_unique <- 2*final_nll + log(n) * df_unique
            tempgcv_unique <- final_nll/(n^2*(1-df_unique/n)^2)
          } else{
            df_unique <- tempaic_unique <- tempbic_unique <- tempgcv_unique <- NA
          }

          out_info[iter,] <- c(lambda,lambda_fusedcoef, mu_smooth_fused, ball_L2)
          out_ests[iter,] <- tempfit$estimate
          out_pen_ngrads[iter,] <- as.vector(tempfit$final_ngrad_pen)
          out_ics[iter,] <- c(tempfit$final_nll, tempfit$final_nll_pen,
                              nPtot, nPtot_selected, df_unique,
                              tempaic, tempaic_unique,
                              tempbic, tempbic_unique,
                              tempgcv, tempgcv_unique)
          out_conv_stats[iter,] <- c(tempfit$niter,tempfit$fit_code,
                                     tempfit$final_step_size,
                                     tempfit$bad_step_count,
                                     extra_start_iter,
                                     tempfit$conv_stats)
          out_starts[iter,] <- startVals_innermost
          nll_pen_trace_mat[iter, ] <- tempfit$nll_pen_trace

          #If this new best value is arising during the first time through the fused coef (middle) loop,
          #then, e.g., this is most likely the result of an unfused grid point, so
          #the estimates are relevant for the next iteration of the outer (parameterwise) loop,
          #which e.g., will also start unfused
          if(warm_start & lambda_fusedcoef_iter==1){
            startVals_outer <- tempfit$estimate
            #again, Wang (2014) suggests tempering the step size along the path but I'm setting that aside for now.
            # if(fit_method %in% c("prox_grad")){
            #   step_size_init_outer <- tempfit$final_step_size
            # }
          }

          #we also want to tell the middle (fused lasso) loop where to start at the next iteration,
          #so if we've found a new best value for this lambda/fusedlambda pair, then we want to
          #store it so that the next fused increment also starts from a warm spot for this fused value
          if(warm_start){
            startVals_middle <- tempfit$estimate
          }

        } #END OF IF STATEMENT HOUSING UPDATES IF NEW BEST VALUE IS REACHED

      } #END OF INNER (EXTRA STARTS) LOOP

      iter <- iter + 1

    } #END OF MIDDLE (FUSED COEFFICIENT LAMBDA) LOOP

  } #END OF OUTER (PARAMETERWISE LAMBDA) LOOP

  # browser()

  ##MAKE SOME SUMMARY PLOTS AND FINALIZE THE OUTPUT TABLES##
  ##******************************************************##

  rownames(out_info) <- NULL
  rownames(out_ics) <- NULL
  rownames(out_ests) <- NULL
  rownames(out_pen_ngrads) <- NULL
  rownames(out_starts) <- NULL
  rownames(out_conv_stats) <- NULL
  rownames(nll_pen_trace_mat) <- NULL

  return(list(info=as.data.frame(out_info),
              ests=as.data.frame(out_ests),
              ics=as.data.frame(out_ics),
              conv_stats=as.data.frame(out_conv_stats),
              pen_trace_mat=nll_pen_trace_mat,
              pen_ngrads=out_pen_ngrads,
              starts=out_starts,
              hazard=hazard, frailty=frailty, model=model,
              penalty=penalty, a=a, ball_L2=ball_L2,
              select_tol=select_tol, fusion_tol=fusion_tol,
              penalty_fusedcoef=penalty_fusedcoef,
              #penalty_fusedbaseline=penalty_fusedbaseline,
              #lambda_fusedbaseline=lambda_fusedbaseline,
              fit_method=fit_method))

}












#' Proximal Gradient Descent Algorithm
#'
#' This function runs a proximal gradient descent algorithm with backtracking similar to that presented by
#'   Wang et al. (2014).
#'
#' @inheritParams nll_func
#' @param penalty A string value indicating the form of parameterwise penalty
#'   to apply. "lasso", "scad", and "mcp" are the options.
#' @param lambda The strength of the parameterwise penalty. Either a single non-negative numeric value
#'   for all three transitions, or a length 3 vector with elements corresponding to the three transitions.
#' @param a For two-parameter penalty functions (e.g., scad and mcp), the second parameter.
#' @param penalty_fusedcoef A string value indicating the form of the fusion penalty to apply
#'   to the regression parameters. "none" and "fusedlasso" are the options.
#' @param lambda_fusedcoef The strength of the fusion penalty on the regression parameters.
#'   Either a single non-negative numeric value
#'   for all three transitions, or a length 3 vector with elements corresponding to the three transitions.
#' @param penweights_list A list of numeric vectors representing weights for each
#'   penalty term (e.g., for adaptive lasso.) Elements of the list should be indexed by the
#'   names "coef1", "coef2", "coef3", "fusedcoef12", "fusedcoef13", "fusedcoef23", "fusedbaseline12", "fusedbaseline13", and "fusedbaseline23"
#' @param mu_smooth_fused A non-negative numeric value for the Nesterov smoothing parameter applied to the fusion penalty.
#' @param step_size_init Positive numeric value for the initial step size.
#' @param step_size_min Positive numeric value for the minimum allowable step size to allow during backtracking.
#' @param step_size_max Positive numeric value for the maximum allowable step size to allow by size increase at each iteration.
#' @param step_size_scale Positive numeric value for the multiplicative change in step size at each step of backtracking.
#' @param ball_L2 Positive numeric value for \eqn{l_2} ball constraint around the origin for the regression parameters.
#'   Typically set to \code{Inf} indicating no constraint, otherwise equivalent to an extra \eqn{l_2} penalty.
#' @param maxit Positive integer maximum number of iterations.
#' @param conv_crit String (possibly vector) giving the convergence criterion.
#' @param conv_tol Positive numeric value giving the convergence tolerance for the chosen criterion.
#' @param verbose Numeric indicating the amount of iteration information should be printed to the user.
#'   Higher numbers provide more detailed information to user, but will slow down the algorithm.
#' @param select_tol Positive numeric value for thresholding estimates to be equal to zero.
#'
#' @return A list.
#' @export
proximal_gradient_descent <- function(para, y1, y2, delta1, delta2, yL, anyLT,
                                      Xmat1 = matrix(nrow(length(y1)),ncol=0), #the default is a 'empty' matrix
                                      Xmat2 = matrix(nrow(length(y1)),ncol=0),
                                      Xmat3 = matrix(nrow(length(y1)),ncol=0),
                                      hazard, frailty, model, weights,
                                      basis1=NULL, basis2=NULL, basis3=NULL, basis3_y1=NULL,
                                      basis1_yL=NULL, basis2_yL=NULL,
                                      dbasis1=NULL, dbasis2=NULL, dbasis3=NULL,
                                      penalty, lambda, a,
                                      penalty_fusedcoef, lambda_fusedcoef,
                                      #penalty_fusedbaseline, lambda_fusedbaseline,
                                      penweights_list, mu_smooth_fused,
                                      step_size_init=1, step_size_min = 1e-6, step_size_max = 1e6,
                                      step_size_scale=1/2, ball_L2=Inf, maxit=300, select_tol = 1e-4,
                                      conv_crit = "nll_pen_change", conv_tol=if(lambda>0) lambda/4 else 1e-6,
                                      verbose){


  #THIS IS STANDARD PROXIMAL GRADIENT DESCENT WITH A LINE SEARCH
  #THE BUILDING BLOCK OF THE PISTA ALGORITHM OF WANG ET AL. (2014) THAT I'M USING FOR A START

  # browser()
  ##Set up control pieces##
  ##*********************##

  #conv_crit determines criterion used to judge convergence
  #est_change_2norm' looks at l1 norm of change in parameter values (scaled by l1 norm of prior values): sum(abs(finalVals-prevVals))/sum(abs(prevVals))
  #est_change_max' looks at largest absolute change in a parameter value
  #nll_pen_change' looks at change in regularized log-likelihood
  #maj_grad_norm' looks at l1 norm of majorized gradient
  if(!all(conv_crit %in% c("est_change_2norm","est_change_max","nll_pen_change","suboptimality"))){
    stop("unknown convergence criterion.")
  }
  if(maxit <0){
    stop("maxit must be nonnegative.")
  }

  n <- length(y1)
  Xmat1 <- if(!is.null(Xmat1)) as.matrix(Xmat1) else matrix(nrow=n,ncol=0)
  Xmat2 <- if(!is.null(Xmat2)) as.matrix(Xmat2) else matrix(nrow=n,ncol=0)
  Xmat3 <- if(!is.null(Xmat3)) as.matrix(Xmat3) else matrix(nrow=n,ncol=0)

  nPtot <- length(para)
  nP1 <- ncol(Xmat1)
  nP2 <- ncol(Xmat2)
  nP3 <- ncol(Xmat3)
  nP0 <- nPtot - nP1 - nP2 - nP3

  ##CREATE CONTRAST MATRICES AND DECOMPOSITIONS FOR LATER USE IN PROXIMAL ALGORITHMS##
  ##I HAD PREVIOUSLY DEVISED A TWO-STAGE CONSTRUCTION TO MAKE UNWEIGHTED VERSION, COMPUTE WEIGHTS, AND THEN MAKE WEIGHTED VERSION
  ##BUT THEN I DECIDED THAT FOR NOW, COMPUTING WEIGHTS MIGHT HAPPEN OUTSIDE OF THIS FUNCTION
  #create list of unweighted contrast matrices
  # pen_mat_list <- contrast_mat_list(nP0=nP0,nP1 = nP1,nP2 = nP2,nP3 = nP3,
  #                                      penalty_fusedcoef = penalty_fusedcoef, lambda_fusedcoef = lambda_fusedcoef,
  #                                      penalty_fusedbaseline = penalty_fusedbaseline,lambda_fusedbaseline = lambda_fusedbaseline,
  #                                      hazard = hazard, penweights_list = NULL)
  # for(pen_name in names(D_list_noweight)){
  #create list of weights based on adaptive lasso inverse contrasts
  #   penweights_list[[pen_name]] <- penweights_internal(parahat=mle_optim,
  #                                                                D=pen_mat_list[[var]],
  #                                                                penweight_type="adaptive",addl_penweight=NULL)
  # }

  #create list of (adaptive) weighted contrast matrices
  #if there is no fusion, this should just be an empty list
  pen_mat_list_temp <- contrast_mat_list(nP0=nP0,nP1 = nP1,nP2 = nP2,nP3 = nP3,
                                         penalty_fusedcoef = penalty_fusedcoef,
                                         lambda_fusedcoef = lambda_fusedcoef,
                                         #penalty_fusedbaseline = penalty_fusedbaseline,
                                         #lambda_fusedbaseline = lambda_fusedbaseline,
                                         hazard = hazard, penweights_list = penweights_list)
  #append them into single weighted matrix
  #if there is no fusion, this should return NULL
  pen_mat_w <- do.call(what = rbind,args = pen_mat_list_temp[["pen_mat_list"]])

  #vector with length equal to the total number of rows of pen_mat_w, with each entry lambda corresponding to that contrast
  #if there is no fusion, this should return NULL
  lambda_f_vec <- pen_mat_list_temp[["lambda_f_vec"]]

  #now, the matrix with the penalty baked in for use in smooth approximation of the fused penalty
  if(is.null(pen_mat_w)){
    pen_mat_w_lambda <- NULL
  } else{
    pen_mat_w_lambda <- diag(lambda_f_vec) %*% pen_mat_w #consider shifting to Matrix package version, because it never goes to cpp
  }

  #compute eigenvector decomposition of this weighted contrast matrix
  #if there is no fusion, this should return a list with Q = as.matrix(0) and eigval = 0
  # pen_mat_w_eig <- pen_mat_decomp(pen_mat_w)
  pen_mat_w_eig <- NULL

  ##Run algorithm##
  ##*************##
  nll_pen_trace <- rep(NA,maxit)
  trace_mat <- xcurr <- para
  fit_code <- 4 #follows convention from 'nleqslv' L-BFGS that is 1 if converged, 4 if maxit reached, 3 if stalled
  bad_step_count <- restart_count <- 0


  #note the division by n to put it on the mean scale--gradient descent works better then!
  nll_pen_xcurr <- nll_pen_func(para=xcurr,
                                y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                yL=yL, anyLT=anyLT,
                                Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                hazard=hazard, frailty=frailty, model=model, weights=weights,
                                basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                                dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                penalty=penalty,lambda=lambda, a=a,
                                penweights_list=penweights_list,
                                pen_mat_w_lambda = pen_mat_w_lambda,
                                mu_smooth_fused = mu_smooth_fused)/n
  if(is.na(nll_pen_xcurr)){stop("Initial values chosen yield infinite likelihood values. Consider other initial values.")}

  #note the division by n to put it on the mean scale--gradient descent works better then!
  ngrad_xcurr <- smooth_obj_grad_func(para=xcurr,
                                      y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                      yL=yL, anyLT=anyLT,
                                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                      hazard=hazard, frailty=frailty, model=model, weights=weights,
                                      basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                      basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                                      dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                      penalty=penalty,lambda=lambda, a=a,
                                      penweights_list=penweights_list,
                                      pen_mat_w_lambda = pen_mat_w_lambda,
                                      mu_smooth_fused = mu_smooth_fused)/n


  #constants
  #define things using notation from Wang (2014), which also covers Zhao (2016) and Zhao (2018)
  #outputs
  #monitors

  #I'M TESTING SOMETHING HERE
  # step_L_ti <- step_L
  step_size_ti <- step_size_init #this is 1/L in the notation of Wang (2014)


  i <- 1 #this is instead of 'k' in Wang (2014)
  while(i <= maxit){
    if(verbose >= 2)print(i)
    # if(i==15){browser()}
    ##RUN ACTUAL STEP##
    ##***************##

    #slightly increase step size at each step (though, with line search this might get knocked back down.)
    step_size_ti <- min(step_size_max,step_size_ti/step_size_scale)

    #run prox function:
    #ADMM prox subroutine for fused lasso only if mu_smooth_fused=0
    #soft thresholding according to lambda*step_size*weight
    xnext <- prox_func(para=xcurr-ngrad_xcurr * step_size_ti, prev_para = xcurr,
                       nP1=nP1,nP2=nP2,nP3=nP3,step_size=step_size_ti,
                       penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                       pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig,
                       lambda_f_vec=lambda_f_vec,
                       mu_smooth_fused=mu_smooth_fused, ball_L2=ball_L2)

    #note the division by n to put it on the mean scale--gradient descent works better then!
    nll_pen_xnext <- nll_pen_func(para=xnext,
                                  y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                  yL=yL, anyLT=anyLT,
                                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                  hazard=hazard, frailty=frailty, model=model, weights=weights,
                                  basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                  basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                                  dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                  penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                  pen_mat_w_lambda = pen_mat_w_lambda,
                                  mu_smooth_fused = mu_smooth_fused)/n

    smooth_obj_lqa_pen_xnext <- smooth_obj_lqa_pen_func(para=xnext, prev_para=xcurr,
                                                        y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                                        yL=yL, anyLT=anyLT,
                                                        Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                                        hazard=hazard, frailty=frailty, model=model, weights=weights,
                                                        basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                                        basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                                                        dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                                        penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                                        pen_mat_w_lambda = pen_mat_w_lambda,
                                                        mu_smooth_fused = mu_smooth_fused,
                                                        step_size=step_size_ti)/n

    if(is.nan(nll_pen_xnext)){
      if(verbose >= 4)print("whoops, proposed step yielded infinite likelihood value, just fyi")
      nll_pen_xnext <- smooth_obj_lqa_pen_xnext <- Inf
    }

    while(is.infinite(nll_pen_xnext) || nll_pen_xnext > smooth_obj_lqa_pen_xnext){
      if(step_size_ti < step_size_min){
        if(verbose >= 4)print("step size got too small, accepting result anyways")
        bad_step_count <- bad_step_count + 1
        break
      }
      step_size_ti <- step_size_ti * step_size_scale
      if(verbose >= 4)print(paste0("nll_pen_xnext:",nll_pen_xnext," bigger than smooth_obj_lqa_pen_xnext:",smooth_obj_lqa_pen_xnext))
      if(verbose >= 4)print(paste0("effective step size reduced to: ",step_size_ti))

      xnext <- prox_func(para=xcurr-ngrad_xcurr*step_size_ti, prev_para = xcurr,
                         nP1=nP1,nP2=nP2,nP3=nP3,step_size=step_size_ti,
                         penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                         pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig,
                         lambda_f_vec=lambda_f_vec,
                         mu_smooth_fused=mu_smooth_fused, ball_L2=ball_L2)

      #note the division by n to put it on the mean scale--gradient descent works better then!
      nll_pen_xnext <- nll_pen_func(para=xnext,
                                    y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                    yL=yL, anyLT=anyLT,
                                    Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                    hazard=hazard, frailty=frailty, model=model, weights=weights,
                                    basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                    basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                                    dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                    penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                    pen_mat_w_lambda = pen_mat_w_lambda,
                                    mu_smooth_fused = mu_smooth_fused)/n

      smooth_obj_lqa_pen_xnext <- smooth_obj_lqa_pen_func(para=xnext, prev_para=xcurr,
                                                          y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                                          yL=yL, anyLT=anyLT,
                                                          Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                                          hazard=hazard, frailty=frailty, model=model, weights=weights,
                                                          basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                                          basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                                                          dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                                          penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                                          pen_mat_w_lambda = pen_mat_w_lambda,
                                                          mu_smooth_fused = mu_smooth_fused,step_size=step_size_ti)/n

      if(is.nan(nll_pen_xnext)){
        if(verbose >= 4)print("whoops, proposed step yielded infinite likelihood value, just fyi")
        nll_pen_xnext <- smooth_obj_lqa_pen_xnext <- Inf
      }
    }

    ##UPDATE MONITORS##
    ##***************##

    xprev <- xcurr
    xcurr <- xnext

    # trace_mat <- cbind(trace_mat,xnext)

    #RECORD TRACE RESULT WITHOUT SMOOTHING, TO PUT US ON A COMMON SCALE!
    nll_pen_trace[i] <- nll_pen_func(para=xnext,
                                     y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                     yL=yL, anyLT=anyLT,
                                     Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                     hazard=hazard, frailty=frailty, model=model, weights=weights,
                                     basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                     basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                                     dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                     penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                     pen_mat_w_lambda = pen_mat_w_lambda,mu_smooth_fused = 0)/n

    ##Check for convergence##
    ##*********************##

    #note the division by n to put it on the mean scale--gradient descent works better then!
    ngrad_xcurr <- smooth_obj_grad_func(para=xcurr,
                                        y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                        yL=yL, anyLT=anyLT,
                                        Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                        hazard=hazard, frailty=frailty, model=model, weights=weights,
                                        basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                        basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                                        dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                        penalty=penalty,lambda=lambda, a=a,
                                        penweights_list=penweights_list,
                                        pen_mat_w_lambda = pen_mat_w_lambda,mu_smooth_fused = mu_smooth_fused)/n

    #suboptimality convergence criterion given as omega in Wang (2014)
    #I actually think this might be incorrect IF we use inexact prox fused lasso, or if we have a ball constraint....
    subopt_t <- max(abs(prox_func(para=ngrad_xcurr, prev_para = xprev,
                                  nP1=nP1,nP2=nP2,nP3=nP3,step_size=1,
                                  penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                                  pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig,
                                  lambda_f_vec=lambda_f_vec,
                                  mu_smooth_fused = mu_smooth_fused, ball_L2=ball_L2)))

    est_change_max <- max(abs(xcurr-xprev))
    est_change_2norm <- sqrt(sum((xcurr-xprev)^2))
    nll_pen_change <- nll_pen_xnext-nll_pen_xcurr #these are both SMOOTHED, so they are direct comparison

    if(verbose >= 3){
      print(paste("max change in ests",est_change_max))
      print(paste("suboptimality (max norm of prox grad)", subopt_t))
      print(paste("estimate with max change",names(para)[abs(xcurr-xprev) == est_change_max]))
      print(paste("max norm of change", est_change_2norm)) #essentially a change in estimates norm
      print(paste("change in nll_pen", nll_pen_change))
      print(paste("new nll_pen", nll_pen_xnext))
    }

    if("suboptimality" %in% conv_crit){
      if(subopt_t <= conv_tol){
        fit_code <- 2
        break
      } #conv_tol is called 'epsilon' in the Wang (2014) paper
    }
    if("est_change_2norm" %in% conv_crit){
      if(est_change_2norm < conv_tol){
        fit_code <- 2
        break
      }
    }
    if("est_change_max" %in% conv_crit){
      if(est_change_max < conv_tol){
        fit_code <- 2
        break
      }
    }
    if("nll_pen_change" %in% conv_crit){
      if(abs(nll_pen_change) < conv_tol){
        fit_code <- 2
        break
      }
    }

    #after the updates, now update the nll monitor values
    nll_pen_xcurr <- nll_pen_xnext
    #iterate algorithm
    i <- i + 1
  }

  ##END ALGORITHM##
  ##*************##

  #if algorithm stopped because learning rate dropped too low, that is worth noting
  # if(lr < con[["min_lr"]]){
  #   fit_code <- 3
  # }

  ##Compute traits of final estimates##
  ##*********************************##

  finalVals <- as.numeric(xnext)

  #for some reason, the fusion prox step can sometimes induce very small nonzero
  #beta values in estimates that are supposed to be 0, so for now we threshold them back to 0.
  #replace any of the betas that are under the selection tolerance with 0, ignoring the baseline variables
  if(nP1+nP2+nP3 > 0){
    if(any(abs(finalVals[(1+nP0):nPtot]) < select_tol)) {
      finalVals[(1+nP0):nPtot][ abs(finalVals[(1+nP0):nPtot]) < select_tol] <- 0
    }
  }

  names(finalVals) <- names(para)

  #Here, report final nll on the SUM scale! No division by n
  final_nll <- nll_func(para=finalVals,
                        y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                        yL=yL, anyLT=anyLT,
                        Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                        hazard=hazard, frailty=frailty, model=model, weights=weights,
                        basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                        basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                        dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
  final_pen <- pen_func(para=finalVals,nP1=nP1,nP2=nP2,nP3=nP3,
                        penalty=penalty,lambda=lambda, a=a,penweights_list=penweights_list,
                        pen_mat_w_lambda = pen_mat_w_lambda)
  #Note we're also reporting the final results under the un-smoothed fused lasso

  final_nll_pen <- final_nll + (n*final_pen)

  #note the division by n to put it on the mean scale--gradient descent works better then!
  final_ngrad <- smooth_obj_grad_func(para=finalVals,
                                      y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                      yL=yL, anyLT=anyLT,
                                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                      hazard=hazard, frailty=frailty, model=model, weights=weights,
                                      basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                      basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                                      dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                      penalty=penalty,lambda=lambda, a=a, penweights_list=penweights_list,
                                      pen_mat_w_lambda = pen_mat_w_lambda,
                                      mu_smooth_fused = mu_smooth_fused) #here, we do still report the nesterov-smoothed results

  final_ngrad_pen <- prox_func(para=final_ngrad, prev_para = xprev,
                               nP1=nP1,nP2=nP2,nP3=nP3,step_size=1,
                               penalty=penalty,lambda=lambda, penweights_list=penweights_list,
                               pen_mat_w=pen_mat_w,pen_mat_w_eig=pen_mat_w_eig,
                               lambda_f_vec=lambda_f_vec,
                               mu_smooth_fused = mu_smooth_fused,
                               ball_L2=ball_L2)
  ##Return list of final objects##
  ##****************************##

  return(list(estimate=finalVals, niter=i, fit_code=fit_code,
              final_nll=final_nll, final_nll_pen=final_nll_pen,
              final_ngrad=final_ngrad, final_ngrad_pen=final_ngrad_pen,
              #penalty_fusedbaseline=penalty_fusedbaseline, lambda_fusedbaseline=lambda_fusedbaseline,
              nll_pen_trace = n*nll_pen_trace,
              conv_stats=c(est_change_max=est_change_max,est_change_2norm=est_change_2norm,
                           nll_pen_change=nll_pen_change,subopt=subopt_t),
              final_step_size=step_size_ti,
              bad_step_count=bad_step_count))
}


##***********************************##
####Define Regularization Functions####
##***********************************##

penweights_internal <- function(parahat,D,penweight_type="adaptive",addl_penweight){
  ####function to compute adaptive weights based on input
  ####written to be general enough to apply to either univariate or fused penalty
  if(penweight_type == "adaptive"){
    if(!is.null(D)){
      if(!is.null(parahat)){
        ada_weight <- abs(D %*% parahat)#^(-a) #originally, the adaptive lasso allows the weight to be scaled by an exponent, but that causes confusion I feel.
        if(any(!is.finite(ada_weight))){
          warning("at least one adaptive fusion weight is not finite, setting to 0") #wait, it may not make sense to set to 0
          ada_weight[!is.finite(ada_weight)] <- 0
        }
      } else {
        ada_weight <- rep(1,nrow(D))
      }
    } else if(!is.null(parahat)){
      ada_weight <- abs(parahat)#^(-a) #originally, the adaptive lasso allows the weight to be scaled by an exponent, but that causes confusion I feel.
    } else {
      ada_weight <- rep(1,length(parahat))
    }

    if(!is.null(addl_penweight) ){
      if(length(ada_weight)==length(addl_penweight)){
        ada_weight <- ada_weight * addl_penweight
      } else{
        stop("supplied additional weight vector is incorrect dimension")
      }
    }
    return(ada_weight)

  } else{ #if not adaptive weights, then just return whatever additional weights might exist
    return(addl_penweight)
  }
}

pen_internal <- function(para, penalty, lambda, a, D, penweights){
  #unified function that returns penalties--  ASSUMES CHECKS HAVE BEEN DONE
  #D is a matrix of contrasts, with number of columns equal to the total number of parameters, rows equal to the total number of differences being `fused'
  #e.g., to fuse first two transitions' betas, set D <- cbind(matrix(data=0,nrow=nP1,ncol=7), diag(nP1), -diag(nP2), matrix(data=0,nrow=nP1,ncol=nP1))

  if(is.null(D)){
    beta <- abs(para)
  } else{
    beta <- as.vector(abs(D %*% para))
  }

  if(penalty == "lasso"){
    out <- abs(beta) * lambda
  } else if(penalty == "scad"){
    beta <- abs(beta)
    ind1 <- ifelse(beta <= lambda,1,0)
    ind2 <- ifelse(beta > lambda & beta <= a*lambda ,1,0)
    ind3 <- ifelse(beta > a*lambda ,1,0)
    out <- ind1 * lambda*beta +
      ind2 * (2*a*lambda*beta - beta^2 - lambda^2)/(2*(a-1)) +
      ind3 * lambda^2*(a+1)/2
  } else if(penalty == "mcp"){
    beta <- abs(beta)
    ind1 <- rep(0, length(beta))
    ind1[beta <= a*lambda] <- 1
    out <- ind1 * (lambda*beta - beta^2/(2*a)) +
      (1-ind1)* (a*lambda^2)/2
  } else{
    out <- rep(0, length(beta)) #vector of 0's
  }

  #in principle, any penalty can have weights added, e.g., to set some terms to 0
  if(!is.null(penweights) && length(penweights)==length(beta)){
    return(penweights * out)
  } else{
    # warning("weights supplied to lasso function are NULL, or of incorrect length. ignoring weights.")
    return(out)
  }
}



pen_func <- function(para,nP1,nP2,nP3,
                     penalty, lambda, a,
                     penweights_list,
                     pen_mat_w_lambda){
  #function to compute complete penalty term for penalized (negative) log likelihood
  #this computes it on the 'mean' scale, i.e., the penalty is not scaled by the number of observations, which can happen in the calling function

  #checks on 'a' to ensure no errors
  if(is.null(a)){
    a <- switch(tolower(penalty),"scad"=3.7, "mcp"=3, "lasso"=1)
  }

  # check_pen_params(penalty,penalty_fusedcoef,penalty_fusedbaseline,
  #                  lambda,lambda_fusedcoef,lambda_fusedbaseline,a)

  #redefine lambda and lambda_fusedcoef depending on a single value or three separate values are given
  if(length(lambda)==1){
    lambda1 <- lambda2 <- lambda3 <- lambda
  } else if(length(lambda)==3){
    lambda1 <- lambda[1]; lambda2 <- lambda[2]; lambda3 <- lambda[3]
  } else{ stop("lambda is neither a single value or a 3-vector!!") }

  nPtot <- length(para)
  nP0 <- length(para) - nP1 - nP2 - nP3


  ##COMPUTE PENALTY##
  ##***************##

  pen <- 0

  #break out the beta vectors from the larger parameter vector, using nP0 to correctly pad out the baseline hazard and theta parameters
  if(nP1 != 0){
    # beta1 <- para[(1+nP0):(nP0+nP1)]
    pen <- pen + sum(pen_internal(para=para[(1+nP0):(nP0+nP1)],
                                  penalty=penalty,lambda=lambda1,a=a,D=NULL,penweights=penweights_list[["coef1"]]))
  }
  if(nP2 != 0){
    # beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
    pen <- pen + sum(pen_internal(para=para[(1+nP0+nP1):(nP0+nP1+nP2)],
                                  penalty=penalty,lambda=lambda2,a=a,D=NULL,penweights=penweights_list[["coef2"]]))
  }
  if(nP3 != 0){
    # beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    pen <- pen + sum(pen_internal(para=para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)],
                                  penalty=penalty,lambda=lambda3,a=a,D=NULL,penweights=penweights_list[["coef3"]]))
  }
  ##If fused lasso is being used, then the following comes into play##
  ##THIS MATCHES THE \|C\bbeta\| formulation of the fused penalty seen in Chen (2012)##
  if(!is.null(pen_mat_w_lambda)){
    pen <- pen + norm(pen_mat_w_lambda %*% para, type = "1")
  }

  return(pen)
}

nll_pen_func <- function(para, y1, y2, delta1, delta2, yL, anyLT,
                         Xmat1, Xmat2, Xmat3,
                         hazard, frailty, model, weights,
                         basis1, basis2, basis3, basis3_y1, basis1_yL, basis2_yL,
                         dbasis1, dbasis2, dbasis3,
                         penalty, lambda, a, penweights_list,
                         pen_mat_w_lambda, mu_smooth_fused){

  #general regularized ll function
  #this function should be able to return two possible things:
  #1. the standard penalized negative log likelihood, including parameter-wise penalties and fused penalties
  #2. the negative log likelihood, plus the original parameter-wise penalty, plus the nesterov-smoothed fused penalty
  #if the convex parameterwise lasso penalty is also being smoothed, then this function does not need to exist, because the 'smoothed_obj_func' above
  #is all you need.


  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)

  nll <- nll_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                  frailty=frailty, hazard=hazard, model=model, weights=weights,
                  basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                  basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                  dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)

  #here, we're not passing in the fused penalty, because whether mu_smooth_fused is 0 or nonzero, it is correctly estimated below
  pen <- pen_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,
                  penalty=penalty,lambda=lambda, a=a,
                  penweights_list=penweights_list,
                  pen_mat_w_lambda = NULL)

  out <- nll + (n * pen)

  #include nesterov smoothed result here instead of in above pen_func! This also returns the original fused penalty if mu_smooth_fused=0
  out <- out + n * smoothed_fused_lasso_func(para = para,
                                             pen_mat_w_lambda = pen_mat_w_lambda,
                                             mu_smooth_fused = mu_smooth_fused)

  return(out)
}





pen_prime_internal <- function(para, penalty, lambda, a, penweights){

  ##Function implementing the 'derivative' of various penalty terms (defined everywhere except where para=0)
  #para and lambda are in all of them, a is in SCAD or MCP

  if(penalty == "lasso"){
    out <- lambda
  } else if(penalty == "scad"){
    para <- abs(para)
    ind2 <- ind1 <- rep(0, length(para))
    ind1[para > lambda] <- 1
    ind2[para <= (lambda * a)] <- 1
    out <- lambda * (1 - ind1) +
      ((lambda * a) - para) * ind2/(a - 1) * ind1
  } else if(penalty == "mcp"){ #careful of the sign here! This is from breheny 2-29 slide 15
    ind1 <- rep(0, length(para))
    ind1[abs(para) <= a*lambda] <- 1
    # out <- ind1*(lambda - abs(para)/a)*sign(para)
    out <- ind1*(lambda - abs(para)/a)
  } else{
    out <- rep(0, length(para)) #vector of 0's
  }

  #in principle, any penalty can have weights added, e.g., to set some terms to 0
  if(!is.null(penweights) && length(penweights)==length(para)){
    return(penweights * out)
  } else{
    # warning("weights supplied to lasso function are NULL, or of incorrect length. ignoring weights.")
    return(out)
  }
}


pen_concave_part_prime_func <- function(para,nP1,nP2,nP3,
                                        penalty, lambda, a,
                                        penweights_list){

  #function to compute the gradient of the smooth part of the penalty, following Yao (2018)
  #this is also on the 'mean' rather than the 'sum' scale, so must be scaled appropriately if grad/hess depend on sample size
  #THIS DOES NOT TAKE ANYTHING ABOUT THE FUSED PENALTY, BECAUSE WE DO NOT CONSIDER NONCONVEX FUSED PENALTY AT PRESENT


  #redefine lambda and lambda_fusedcoef depending on a single value or three separate values are given
  if(length(lambda)==1){
    lambda1 <- lambda2 <- lambda3 <- lambda
  } else if(length(lambda)==3){
    lambda1 <- lambda[1]; lambda2 <- lambda[2]; lambda3 <- lambda[3]
  } else{ stop("lambda is neither a single value or a 3-vector!!") }

  nPtot <- length(para)
  nP0 <- length(para) - nP1 - nP2 - nP3

  grad_part <- numeric(nPtot)

  #break out the beta vectors from the larger parameter vector, using nP0 to correctly pad out the baseline hazard and theta parameters
  if(nP1 != 0){
    beta1 <- para[(1+nP0):(nP0+nP1)]
    grad_part[(1+nP0):(nP0+nP1)] <- sign(beta1)*(pen_prime_internal(para=beta1,lambda=lambda1,penalty=penalty,a=a,penweights=penweights_list[["coef1"]]) -
                                                   pen_prime_internal(para=beta1,lambda=lambda1,penalty="lasso",a=a,penweights=penweights_list[["coef1"]]))
  }
  if(nP2 != 0){
    beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
    grad_part[(1+nP0+nP1):(nP0+nP1+nP2)] <- sign(beta2)*(pen_prime_internal(para=beta2,lambda=lambda2,penalty=penalty,a=a,penweights=penweights_list[["coef2"]]) -
                                                           pen_prime_internal(para=beta2,lambda=lambda2,penalty="lasso",a=a,penweights=penweights_list[["coef2"]]))
  }
  if(nP3 != 0){
    beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    grad_part[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)] <- sign(beta3)*(pen_prime_internal(para=beta3,lambda=lambda3,penalty=penalty,a=a,penweights=penweights_list[["coef3"]]) -
                                                                   pen_prime_internal(para=beta3,lambda=lambda3,penalty="lasso",a=a,penweights=penweights_list[["coef3"]]))
  }

  # names(grad_part) <- names(para)
  return(grad_part)
}



##*****************************************##
####Nesterov-Smoothed Objective Functions####
##*****************************************##


#function for the nesterov smoothed approximation of the fused lasso
#Defined in Chen (2012) equation 3.5, with the solution alphastar defined
#in proposition 2
smoothed_fused_lasso_func <- function(para,pen_mat_w_lambda, mu_smooth_fused){
  if(is.null(pen_mat_w_lambda) | is.null(mu_smooth_fused)){
    return(0)
  }
  if(mu_smooth_fused==0){
    return(norm(pen_mat_w_lambda %*% para, type = "1"))
  }
  temp_vec <- pen_mat_w_lambda %*% para
  alphastar <-  sign(temp_vec/mu_smooth_fused)*pmin(1, abs(temp_vec/mu_smooth_fused))

  return( as.vector( t(alphastar) %*% temp_vec - mu_smooth_fused*sum(alphastar^2)/2 ) )
}


#function for the gradient of the nesterov smoothed approximation of the fused lasso
smoothed_fused_lasso_prime_func <- function(para,pen_mat_w_lambda, mu_smooth_fused){
  if(is.null(pen_mat_w_lambda) | is.null(mu_smooth_fused)){
    return(0)
  }
  if(mu_smooth_fused==0){
    return(0)
  }
  temp_vec <- pen_mat_w_lambda %*% para /mu_smooth_fused
  alphastar <-  sign(temp_vec)*pmin(1, abs(temp_vec))

  return( as.vector( t(pen_mat_w_lambda) %*% alphastar ) )
}


smooth_obj_func <- function(para, y1, y2, delta1, delta2, yL, anyLT,
                            Xmat1, Xmat2, Xmat3,
                            hazard, frailty, model, weights,
                            basis1, basis2, basis3, basis3_y1, basis1_yL, basis2_yL,
                            dbasis1, dbasis2, dbasis3,
                            penalty, lambda, a, penweights_list,
                            pen_mat_w_lambda, mu_smooth_fused){

  #function that combines the nll, the smooth concave part of any parameter-wise penalties, and optionally the nesterov-smoothed penalty
  #this function should be able to return two possible things:
  #1. The negative log-likelihood, plus the smooth nonconvex part of any penalties
  #this possibility is for the proximal gradient method that solves the fused lasso proximal operator,
  #so we just want the smooth part of the objective function
  #2. The above smooth function, plus the nesterov-smoothed lasso fused lasso components


  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)

  out <- nll_func(para=para, y1=y1, y2=y2, yL=yL, anyLT=anyLT,
                  delta1=delta1, delta2=delta2,
                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                  frailty=frailty, hazard=hazard, model=model, weights=weights,
                  basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                  basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                  dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)

  #note that this only accounts for the elementwise penalties, because the fused penalty is either convex, or smoothed in the next step
  if(penalty %in% c("scad","mcp")){
    componentwise_pen <- pen_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,
                                  penalty=penalty,lambda=lambda, a=a,
                                  penweights_list=penweights_list,
                                  pen_mat_w_lambda = NULL)

    neg_componentwise_pen <- pen_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,
                                      penalty="lasso",lambda=lambda, a=a,
                                      penweights_list=penweights_list,
                                      pen_mat_w_lambda = NULL)
    out <- out  + (n * componentwise_pen) - (n * neg_componentwise_pen)
  }

  #if we are also nesterov smoothing the fused penalty, include that!
  if(!is.null(mu_smooth_fused) && mu_smooth_fused > 0){
    out <- out + n * smoothed_fused_lasso_func(para = para, pen_mat_w_lambda = pen_mat_w_lambda, mu_smooth_fused = mu_smooth_fused)
  }

  return(out)
}




smooth_obj_grad_func <- function(para, y1, y2, delta1, delta2, yL, anyLT,
                                 Xmat1, Xmat2, Xmat3,
                                 hazard, frailty, model, weights,
                                 basis1, basis2, basis3, basis3_y1, basis1_yL, basis2_yL,
                                 dbasis1, dbasis2, dbasis3,
                                 penalty, lambda, a,
                                 penweights_list,
                                 pen_mat_w_lambda, mu_smooth_fused){
  #general gradient function of smooth part of penalized negative loglikelihood
  #corresponds to the smooth_obj_func from above

  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)

  out <- ngrad_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                    Xmat1=Xmat1, Xmat2=Xmat2,
                    Xmat3=Xmat3, frailty=frailty, hazard=hazard, model=model, weights=weights,
                    basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                    basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                    dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)

  if(tolower(penalty) %in% c("scad","mcp")){
    out <- out + n * pen_concave_part_prime_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,
                                                 penalty=penalty,lambda=lambda, a=a,
                                                 penweights_list=penweights_list)
  }

  #if we are also nesterov smoothing the fused lasso component, include that!
  if(!is.null(mu_smooth_fused) && mu_smooth_fused > 0){
    out <- out + n * smoothed_fused_lasso_prime_func(para=para, pen_mat_w_lambda=pen_mat_w_lambda, mu_smooth_fused=mu_smooth_fused)
  }

  return( out )

}





smooth_obj_lqa_pen_func <- function(para, prev_para,
                                    y1, y2, delta1, delta2, yL, anyLT,
                                    Xmat1, Xmat2, Xmat3,
                                    hazard, frailty, model, weights,
                                    basis1, basis2, basis3, basis3_y1,
                                    basis1_yL, basis2_yL,
                                    dbasis1, dbasis2, dbasis3,
                                    penalty, lambda, a, penweights_list,
                                    pen_mat_w_lambda, mu_smooth_fused,step_size){

  #function for the lasso-penalized local quadratic approximation of the smooth objective function
  #a quadratic expansion of the nll around prev_para, evaluated at para
  #following Wang (2014), equation 3.7


  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)

  #the original smooth objective at prev_para
  obj_val <- smooth_obj_func(para=prev_para, y1=y1, y2=y2, yL=yL, anyLT=anyLT,
                             delta1=delta1, delta2=delta2,
                             Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                             frailty=frailty, hazard=hazard, model=model, weights=weights,
                             basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                             basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                             dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                             penalty=penalty, lambda=lambda, a=a,
                             penweights_list=penweights_list,
                             pen_mat_w_lambda=pen_mat_w_lambda, mu_smooth_fused=mu_smooth_fused)

  #the corresponding grad of the smoothed objective function, evaluated at prev_para
  #oops this function is technically not defined until below but go look for it it's there!
  smooth_ngrad_temp <- smooth_obj_grad_func(para=prev_para,
                                            y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                                            Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                                            frailty=frailty, hazard=hazard, model=model, weights=weights,
                                            basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                                            basis1_yL=basis1_yL, basis2_yL=basis2_yL,
                                            dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                                            penalty=penalty, lambda=lambda, a=a,
                                            penweights_list=penweights_list,
                                            pen_mat_w_lambda=pen_mat_w_lambda, mu_smooth_fused=mu_smooth_fused)

  #the linear term of the taylor expansion
  out <- obj_val + t(smooth_ngrad_temp) %*% (para - prev_para)

  #the quadratic term of the expansion (notice I had to multiply by n, originally forgot that)
  out <- out + n/(2*step_size) * sum((para - prev_para)^2)

  #As a last step, add in any non-smooth penalties that remain after we've accounted for the smoothness above
  #Here, it is always lasso because SCAD and MCP have had smooth concave component absorbed above.
  pen <- pen_func(para=para,nP1=nP1,nP2=nP2,nP3=nP3,
                  penalty="lasso",lambda=lambda, a=a,
                  penweights_list=penweights_list,
                  pen_mat_w_lambda = NULL) #notice, we are separating out fused lasso piece, which is separately added below.
  out <- out + (n * pen)

  #Separately, we add any non-smooth fused lasso that remains.
  if(is.null(mu_smooth_fused) || mu_smooth_fused==0){
    out <- out + n * smoothed_fused_lasso_func(para = para,
                                               pen_mat_w_lambda = pen_mat_w_lambda,
                                               mu_smooth_fused = mu_smooth_fused)
  }

  return(out)
}



##***********************##
####Helper Proximal operator functions####
##***********************##

lasso_prox_internal <- function(beta,lambda,step_size,penweights){

  #function to soft-threshold a vector at lambda*penweights*step_size
  #step_size is step size

  if(!is.null(penweights) && length(penweights)==length(beta)){
    return( as.matrix(  sign(beta)*pmax(0, abs(beta)-(step_size*lambda*penweights)) ) )
  } else{
    # warning("weights supplied to lasso function are NULL, or of incorrect length. ignoring weights.")
    return( as.matrix(  sign(beta)*pmax(0, abs(beta)-(step_size*lambda)) ) )
  }
}

prox_func <- function(para, prev_para, nP1, nP2, nP3, step_size,
                      penalty, lambda, penweights_list,
                      pen_mat_w,pen_mat_w_eig=NULL,lambda_f_vec,
                      mu_smooth_fused, ball_L2=Inf){

  #perform proximal operator based on convex part of penalty term (following Yao (2018))

  # a actually never gets used here...
  # check_pen_params(penalty,penalty_fusedcoef,penalty_fusedbaseline,
  #                  lambda,lambda_fusedcoef,lambda_fusedbaseline,a)

  #redefine lambda and lambda_fusedcoef depending on a single value or three separate values are given
  if(length(lambda)==1){
    lambda1 <- lambda2 <- lambda3 <- lambda
  } else if(length(lambda)==3){
    lambda1 <- lambda[1]; lambda2 <- lambda[2]; lambda3 <- lambda[3]
  } else{ stop("lambda is neither a single value or a 3-vector!!") }

  nPtot <- length(para)
  nP0 <- length(para) - nP1 - nP2 - nP3

  ###COMPUTE THE FUSED LASSO PROXIMAL STEP USING ADMM
  eps_num <- min(sqrt(.Machine$double.eps), 100 * .Machine$double.eps) #definition from Smurf code
  if(!is.null(pen_mat_w) && mu_smooth_fused == 0){
    stop("proximal function for fused penalty requires nesterov smoothing, i.e., mu_smooth_fused>0.")

    # para_fl <- admm_po_cpp(beta_tilde = para,
    #                        slambda = lambda_f_vec * step_size,
    #                        penmat = pen_mat_w,
    #                        Q = if(!is.null(pen_mat_w_eig)) pen_mat_w_eig$Q else as.matrix(0),
    #                        eigval =  if(!is.null(pen_mat_w_eig)) pen_mat_w_eig$eigval else 0,
    #                        fast = if(!is.null(pen_mat_w_eig)) all(abs(pen_mat_w_eig$eigval) >= eps_num) else FALSE,
    #                        maxiter = 1e4, rho = 1,
    #                        beta_old = prev_para)

  } else{
    para_fl <- para
  }


  prox_out <- para_fl #i think this could equivalently be para but whatever

  #break out the beta vectors from the larger parameter vector, using nP0 to correctly pad out the baseline hazard and theta parameters
  if(nP1 != 0){
    beta1 <- para_fl[(1+nP0):(nP0+nP1)]
    prox_out[(1+nP0):(nP0+nP1)] <- lasso_prox_internal(beta=beta1,lambda=lambda1,step_size=step_size,penweights=penweights_list[["coef1"]])
  }
  if(nP2 != 0){
    beta2 <- para_fl[(1+nP0+nP1):(nP0+nP1+nP2)]
    prox_out[(1+nP0+nP1):(nP0+nP1+nP2)] <- lasso_prox_internal(beta=beta2,lambda=lambda2,step_size=step_size,penweights=penweights_list[["coef2"]])
  }
  if(nP3 != 0){
    beta3 <- para_fl[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    prox_out[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)] <- lasso_prox_internal(beta=beta3,lambda=lambda3,step_size=step_size,penweights=penweights_list[["coef3"]])
  }

  #now, add projection of the covariates onto the ball of radius R to potentially accommodate the constraints of Wang (2014)
  #this is equivalent to ridge regression
  if(nP1 + nP2 + nP3 > 0 && !is.null(ball_L2) && !is.infinite(ball_L2)){
    temp_norm2 <- sum(prox_out[(1+nP0):(nP0+nP1+nP2+nP3)]^2)
    if(temp_norm2 > ball_L2){
      prox_out[(1+nP0):(nP0+nP1+nP2+nP3)] <- prox_out[(1+nP0):(nP0+nP1+nP2+nP3)] * ball_L2 / temp_norm2
    }
  }

  # names(prox_out) <- names(para)
  return(prox_out)
}



#### Helper Fused Lasso Functions ####


#' Create List of Contrast Matrices for Parameter Fusion
#'
#' Creates list of numeric matrices corresponding to the fusion of blocks of parameters.
#'   The list is indexed by pairs of parameter blocks being fused, i.e.,
#'   \itemize{
#'   \item \code{fusedcoef12},\code{fusedcoef13},\code{fusedcoef23} list elements correspond to
#'     matrices that, when multiplied to the parameter vector, yield a vector of fused differences
#'     between corresponding regression parameters between the stated transition hazards.
#'   \item \code{fusedbaseline12},\code{fusedbaseline13},\code{fusedbaseline23} list elements
#'     are as above, but fusing the baseline hazard parameters.
#'   }
#'
#' @inheritParams proximal_gradient_descent
#' @param nP0,nP1,nP2,nP3 Number of parameters.
#'
#' @return Returns a list with elements named as above, each containing a numeric contrast matrix.
#' @export
contrast_mat_list <- function(nP0,nP1,nP2,nP3, #NEEDS UPDATE TO WORK WITH WEIBULL
                              penalty_fusedcoef,lambda_fusedcoef,
                              # penalty_fusedbaseline,lambda_fusedbaseline,
                              hazard, penweights_list){

  ##FUNCTIONS TO CREATE CONTRAST MATRICES FOR FUSED PENALIZATION##

  #create single (weighted) contrast matrix that applies all of the differences at once
  #some circularity is induced because we need matrix to get adaptive weights to get weighted matrix, but that's ok

  pen_mat_list <-list()
  lambda_f_vec <- NULL
  nPtot <- nP0 + nP1 + nP2 + nP3

  if(tolower(penalty_fusedcoef) %in% c("fusedlasso","adafusedlasso")){

    if(length(lambda_fusedcoef)==1){
      lambda_fusedcoef12 <- lambda_fusedcoef13 <- lambda_fusedcoef23 <- lambda_fusedcoef
    } else if(length(lambda_fusedcoef)==3){
      lambda_fusedcoef12 <- lambda_fusedcoef[1]; lambda_fusedcoef13 <- lambda_fusedcoef[2];lambda_fusedcoef23 <- lambda_fusedcoef[3]
    } else{ stop("lambda_fusedcoef is neither a single value or a 3-vector!!") }

    #create a difference matrix connecting first and second sets of betas, left padded by 7 columns of zeros for the hazards
    if(lambda_fusedcoef12 != 0){
      stopifnot(nP1==nP2)
      #in principle, any penalty can have weights added, e.g., to set some terms to 0
      if(!is.null(penweights_list[["fusedcoef12"]]) && length(penweights_list[["fusedcoef12"]])==nP1){
        pen_mat_list[["fusedcoef12"]] <- cbind(matrix(data=0,nrow=nP1,ncol=nP0),
                                               diag(as.vector(penweights_list[["fusedcoef12"]])),
                                               -diag(as.vector(penweights_list[["fusedcoef12"]])),
                                               matrix(data=0,nrow=nP1,ncol=nP1))
      } else{
        pen_mat_list[["fusedcoef12"]] <- cbind(matrix(data=0,nrow=nP1,ncol=nP0), diag(nP1), -diag(nP1), matrix(data=0,nrow=nP1,ncol=nP1))
      }
      lambda_f_vec <- c(lambda_f_vec, rep(lambda_fusedcoef12, nP1))
    }

    #create a difference matrix connecting first and third sets of betas, left padded by 7 columns of zeros for the hazards
    if(lambda_fusedcoef13 != 0){
      stopifnot(nP1==nP3)
      #in principle, any penalty can have weights added, e.g., to set some terms to 0
      if(!is.null(penweights_list[["fusedcoef13"]]) && length(penweights_list[["fusedcoef13"]])==nP1){
        pen_mat_list[["fusedcoef13"]] <- cbind(matrix(data=0,nrow=nP1,ncol=nP0),
                                               diag(as.vector(penweights_list[["fusedcoef13"]])),
                                               matrix(data=0,nrow=nP1,ncol=nP1),
                                               -diag(as.vector(penweights_list[["fusedcoef13"]])))
      } else{
        pen_mat_list[["fusedcoef13"]] <- cbind(matrix(data=0,nrow=nP1,ncol=nP0), diag(nP1),matrix(data=0,nrow=nP1,ncol=nP1),-diag(nP1))
      }
      lambda_f_vec <- c(lambda_f_vec, rep(lambda_fusedcoef13, nP1))
    }

    #create a difference matrix connecting second and third sets of betas, left padded by 7 columns of zeros for the hazards
    if(lambda_fusedcoef23 != 0){
      stopifnot(nP2==nP3)
      #in principle, any penalty can have weights added, e.g., to set some terms to 0
      if(!is.null(penweights_list[["fusedcoef23"]]) && length(penweights_list[["fusedcoef23"]])==nP2){
        pen_mat_list[["fusedcoef23"]] <- cbind(matrix(data=0,nrow=nP2,ncol=nP0),
                                               matrix(data=0,nrow=nP2,ncol=nP2),
                                               diag(as.vector(penweights_list[["fusedcoef23"]])),
                                               -diag(as.vector(penweights_list[["fusedcoef23"]])))
      } else{
        pen_mat_list[["fusedcoef23"]] <- cbind(matrix(data=0,nrow=nP2,ncol=nP0), matrix(data=0,nrow=nP2,ncol=nP2),diag(nP2),-diag(nP2))
      }
      lambda_f_vec <- c(lambda_f_vec, rep(lambda_fusedcoef23, nP2))
    }

  }

  # #Finally, add penalty of fusion for baseline parameters
  # if(tolower(penalty_fusedbaseline) != "none" & !(is.null(hazard))){
  #   if(!(tolower(hazard) %in% c("weibull","wb"))){
  #     stop("non-Weibull not yet implemented")
  #   }
  #
  #   if(length(lambda_fusedbaseline)==1){
  #     lambda_fusedbaseline12 <- lambda_fusedbaseline13 <- lambda_fusedbaseline23 <- lambda_fusedbaseline
  #   } else if(length(lambda_fusedbaseline)==3){
  #     lambda_fusedbaseline12 <- lambda_fusedbaseline[1]
  #     lambda_fusedbaseline13 <- lambda_fusedbaseline[2]
  #     lambda_fusedbaseline23 <- lambda_fusedbaseline[3]
  #   } else{ stop("lambda_fusedbaseline is neither a single value or a 3-vector!!") }
  #
  #   #create a difference matrix connecting first and second sets of baseline parameters
  #   if(lambda_fusedbaseline12 != 0){
  #     D_temp <- matrix(data=0,nrow=2,ncol=nPtot)
  #
  #     #FIGURE OUT A BETTER WAY TO INCORPORATE THESE WEIGHTS
  #     if(!is.null(penweights_list[["fusedbaseline12"]]) && length(penweights_list[["fusedbaseline12"]])==2){
  #       D_temp[1,1] <- D_temp[2,2] <- penweights_list[["fusedbaseline12"]][1]
  #       D_temp[1,3] <- D_temp[2,4] <- -penweights_list[["fusedbaseline12"]][2]
  #     } else{
  #       D_temp[1,1] <- D_temp[2,2] <- 1
  #       D_temp[1,3] <- D_temp[2,4] <- -1
  #     }
  #     pen_mat_list[["fusedbaseline12"]] <- D_temp
  #     lambda_f_vec <- c(lambda_f_vec, rep(lambda_fusedbaseline12, 2))
  #   }
  #   #create a difference matrix connecting first and third sets of baseline parameters
  #   if(lambda_fusedbaseline13 != 0){
  #     D_temp <- matrix(data=0,nrow=2,ncol=nPtot)
  #
  #     if(!is.null(penweights_list[["fusedbaseline13"]]) && length(penweights_list[["fusedbaseline13"]])==2){
  #       D_temp[1,1] <- D_temp[2,2] <- penweights_list[["fusedbaseline13"]][1]
  #       D_temp[1,5] <- D_temp[2,6] <- -penweights_list[["fusedbaseline13"]][2]
  #
  #     } else{
  #       D_temp[1,1] <- D_temp[2,2] <- 1
  #       D_temp[1,5] <- D_temp[2,6] <- -1
  #     }
  #     pen_mat_list[["fusedbaseline13"]] <- D_temp
  #     lambda_f_vec <- c(lambda_f_vec, rep(lambda_fusedbaseline13, 2))
  #   }
  #   #create a difference matrix connecting second and third sets of baseline parameters
  #   if(lambda_fusedbaseline23 != 0){
  #     D_temp <- matrix(data=0,nrow=2,ncol=nPtot)
  #
  #     if(!is.null(penweights_list[["fusedbaseline23"]]) && length(penweights_list[["fusedbaseline23"]])==2){
  #       D_temp[1,3] <- D_temp[2,4] <- penweights_list[["fusedbaseline23"]]
  #       D_temp[1,5] <- D_temp[2,6] <- -penweights_list[["fusedbaseline23"]]
  #     } else{
  #       D_temp[1,3] <- D_temp[2,4] <- 1
  #       D_temp[1,5] <- D_temp[2,6] <- -1
  #     }
  #     pen_mat_list[["fusedbaseline23"]] <- D_temp
  #     lambda_f_vec <- c(lambda_f_vec, rep(lambda_fusedbaseline23, 2))
  #   }
  # }

  return(list(pen_mat_list=pen_mat_list,lambda_f_vec=lambda_f_vec))
}



pen_mat_decomp <- function(pen_mat) {

  #FROM SMURF CODE "PENALTY_MATRICES"
  # Compute eigenvalue decomposition of t(pen_mat[[j]]) %*% pen_mat[[j]] for fused penalties
  # except "none", "lasso" and "grouplasso"
  #
  # pen_mat: (weighted) penalty matrix


  pen_mat_aux <- list()
  # Return NULL if error
  tmp <- tryCatch(eigen(t(pen_mat) %*% pen_mat), error = function(e) NULL)

  if (!is.null(tmp)) {
    # Get eigenvectors and -values if eigen did not give an error
    pen_mat_aux$Q <- tmp$vectors
    pen_mat_aux$eigval <- tmp$values

  } else {
    # eigen gave an error, use slower ADMM version (check happens when calling C++ code)
    pen_mat_aux$Q <- as.matrix(0)
    pen_mat_aux$eigval <- 0
  }
  return(pen_mat_aux)
}

