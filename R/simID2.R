#' The function that simulates independent/cluster-correlated semi-competing
#'   risks data under Markov/semi-Markov Weibull models.
#'
#' @inheritParams FreqID_HReg2
#' @param id A vector of cluster information for \code{n} subjects.
#' The cluster membership must be set to consecutive positive integers, \code{1:J}.
#' Required only when generating clustered data.
#' @param x1,x2,x3 Covariate matrices with \code{n} rows.
#' @param beta1.true,beta2.true,beta3.true Vectors of true regression
#' parameter values. The length of each vector should equal
#' the number of columns in the corresponding covariate matrix.
#' @param alpha1.true,alpha2.true,alpha3.true,kappa1.true,kappa2.true,kappa3.true scalar true baseline parameter values.
#' @param phi1.true,phi2.true,phi3.true Vectors of true baseline parameter values.
#' @param theta.true True value for \eqn{\theta}.
#' @param SigmaV.true True value for covariance matrix of MVN
#' cluster-level random effects. Required only when generating clustered data.
#' Should be a numeric \eqn{J\times J} matrix.
#' @param frailty_type string denoting "gamma" for gamma-distributed frailty
#' with variance \code{theta}, or "lognormal" for lognormal distributed
#' frailty with log-frailty variance \code{theta}
#' @param beta2frail.true,beta3frail.true scalar for the coefficient
#' of the shared frailty in the \eqn{h_2} and \eqn{h_3} submodels.
#' @param beta3tv.true Vectors of true regression parameter values
#' for direct effects of \eqn{T_1} on the \eqn{h_3} submodel.
#' @param h3tv_degree either the string "cs" indicating restricted
#' cubic spline, or an integer for degree of time-varying hazard/odds ratio
#' B-spline basis. (0 is piecewise constant)
#' @param LT_interval A numeric vector of two elements. The left truncation
#'  times are generated from Uniform(\eqn{LT_interval[1]}, \eqn{LT_interval[2]}).
#'  Setting to \code{c(0,0)} corresponds with no left-truncation.
#' @param cens A numeric vector of two elements. The right censoring
#'  times are generated from Uniform(\eqn{cens[1]}, \eqn{cens[2]}).
#'  Setting to \code{c(0,0)} corresponds with no independent right censoring.
#'
#' @return returns a data.frame containing semi-competing risks
#'   outcomes from \code{n} subjects, with columns as follows: \cr
#'   \describe{
#'   \item{y1}{a vector of \code{n} times for the non-terminal event}
#'   \item{delta1}{a vector of \code{n} censoring indicators for the non-terminal event time (1=event occurred, 0=censored)}
#'   \item{y2}{a vector of \code{n} times for the terminal event}
#'   \item{delta2}{a vector of \code{n} censoring indicators for the terminal event time (1=event occurred, 0=censored)}
#'   \item{yL}{a vector of \code{n} times of left censoring}
#'   \item{gamma.true}{a vector of true frailty values}
#'   \item{x3tv}{a matrix of basis functions corresponding to the direct
#'   time-varying effect of \eqn{T_1} on \eqn{h_3}, if applicable.}
#'   \item{id}{a vector of cluster membership labels}
#'   \item{Vmat}{an \eqn{n \times 3} matrix of true cluster-specific random effects}
#'   }
#'   \cr
#'   Additionally, if a direct time-varying effect of \eqn{T_1} on \eqn{h_3} is
#'   specified, the returned data frame is given the following attributes:
#'   \describe{
#'   \item{p3tv}{total number of parameters corresponding to effect}
#'   \item{h3tv_degree}{degree of spline effect, as specified by user}
#'   \item{h3tv_knots}{vector of knots corresponding with effect}
#'   }
#'
#' @export
simID2 <- function(id = NULL, x1, x2, x3,
                  beta1.true, beta2.true, beta3.true,
                  alpha1.true, alpha2.true, alpha3.true,
                  kappa1.true, kappa2.true, kappa3.true,
                  phi1.true, phi2.true, phi3.true,
                  theta.true, SigmaV.true = NULL, hazard="weibull", knots_list,
                  model="semi-markov", frailty_type="gamma",
                  beta2frail.true=1, beta3frail.true=1,
                  beta3tv.true=NULL, h3tv_degree=3,
                  LT_interval = c(0,0), cens = c(0,0)) {
  # browser()

  if (!is.null(id) & is.null(SigmaV.true)) {
    stop("SigmaV.true must be given to simulate correlated data")
  }

  n <- dim(x1)[1]
  p1 <- dim(x1)[2]
  p2 <- dim(x2)[2]
  p3 <- dim(x3)[2]
  anyLT <- max(LT_interval) > 0

  if(theta.true >0) {
    if(tolower(frailty_type)=="gamma"){
      gamma.true <- stats::rgamma(n, 1/theta.true, 1/theta.true)
    } else {
      gamma.true <- stats::rlnorm(n, meanlog = 0, sdlog = sqrt(theta.true))
    }
  } else if(theta.true == 0){
    gamma.true <- rep(1, n)
  }

  LP1	<- if(p1>0) as.vector(beta1.true %*% t(x1)) else 0
  LP2	<- if(p2>0) as.vector(beta2.true %*% t(x2)) else 0
  LP3	<- if(p3>0) as.vector(beta3.true %*% t(x3)) else numeric(n) #made a vector bc it is subset automatically by yesR below

  #incorporate clustering random effects
  if (!is.null(id)) {
    J <- length(unique(id))
    nj <- as.vector(table(id))
    Vmat <- MASS::mvrnorm(J, rep(0, 3), SigmaV.true)
    Vmat <- cbind(V1=rep(Vmat[,1], nj),V2=rep(Vmat[,2], nj),V3=rep(Vmat[,3], nj))
    LP1 <- LP1 + as.vector(Vmat[,1])
    LP2 <- LP2 + as.vector(Vmat[,2])
    LP3 <- LP3 + as.vector(Vmat[,3])
  } else{ Vmat <- NULL }

  #if not weibull, must be piecewise, so set that up
  if(!(tolower(hazard) %in% c("weibull","wb"))){
    #make it so that index can be made rowwise whether there is 1 col or 3
    #function already accounts for whether or not knots vectors begin with 0
    if(is.list(knots_list)){
      if(length(knots_list) == 3){
        knots01 <- knots_list[[1]]
        knots02 <- knots_list[[2]]
        knots03 <- knots_list[[3]]
      } else if(length(knots_list)==1){
        knots01 <- knots02 <- knots03 <- knots_list[[1]]
      } else{
        stop("knots must be either vector of knots lambdas,
           or list with three elements corresponding to hazard-specific knot vectors")
      }
    } else{
      if(is.vector(knots_list)){
        knots03 <- knots02 <- knots01 <- knots_list
      } else{
        stop("knots must be either vector of knots lambdas,
           or list with three elements corresponding to hazard-specific knot vectors")
      }
    }
  }

  if(anyLT){
    yL <- stats::runif(n,min=LT_interval[1],max=LT_interval[2])
    if(tolower(hazard) %in% c("weibull","wb")){
      R_bound <- stats::pweibull(yL, lower.tail = FALSE, shape = alpha1.true,
                  scale = exp(-(log(kappa1.true) + LP1 + log(gamma.true))/alpha1.true))
      D_bound <- stats::pweibull(yL, lower.tail = FALSE, shape = alpha2.true,
                  scale = exp(-(log(kappa2.true) + LP2 + beta2frail.true * log(gamma.true))/alpha2.true))
      R_prob <- stats::runif(n,min=0, max=R_bound)
      D_prob <- stats::runif(n,min=0, max=D_bound)
      R <- stats::qweibull(R_prob, lower.tail = FALSE, shape = alpha1.true,
                  scale = exp(-(log(kappa1.true) + LP1 + log(gamma.true))/alpha1.true))
      D <- stats::qweibull(D_prob, lower.tail = FALSE, shape = alpha2.true,
                  scale = exp(-(log(kappa2.true) + LP2 + beta2frail.true * log(gamma.true))/alpha2.true))
    } else{
      R_bound <- ppwexp(yL, lower.tail = FALSE, log.p=FALSE,
                        phi=phi1.true, knots_vec=knots01, eta = LP1 + log(gamma.true))
      D_bound <- ppwexp(yL, lower.tail = FALSE, log.p=FALSE,
                        phi=phi2.true, knots_vec=knots02, eta = LP2 + beta2frail.true * log(gamma.true))
      R_prob <- stats::runif(n,min=0, max=R_bound)
      D_prob <- stats::runif(n,min=0, max=D_bound)
      R <- qpwexp(R_prob, lower.tail = FALSE, phi=phi1.true, knots_vec=knots01,
                  eta = LP1 + log(gamma.true))
      D <- qpwexp(D_prob, lower.tail = FALSE, phi=phi2.true, knots_vec=knots02,
                  eta = LP2 + beta2frail.true * log(gamma.true))
    }
  } else{
    yL <- numeric(n)
    if(tolower(hazard) %in% c("weibull","wb")){
      R <- stats::rweibull(n, shape = alpha1.true,
            scale = exp(-(log(kappa1.true) + LP1 + log(gamma.true))/alpha1.true))
      D <- stats::rweibull(n, shape = alpha2.true,
            scale = exp(-(log(kappa2.true) + LP2 + beta2frail.true * log(gamma.true))/alpha2.true))
    } else{
      #still generate data by inverse probability
      R <- qpwexp(stats::runif(n), lower.tail = FALSE, phi=phi1.true,
                         knots_vec=knots01, eta = LP1 + log(gamma.true))
      D <- qpwexp(stats::runif(n), lower.tail = FALSE, phi=phi2.true,
                         knots_vec=knots02, eta = LP2 + beta2frail.true * log(gamma.true))
    }
  }

  yesR <- R < D

  #now, incorporate a possibly time-varying component into LP3 (semi-markov only)
  #if we set total number of parameters to 0, then we have no time-varying component.
  if( !(tolower(model) %in% c("markov","m")) && !is.null(beta3tv.true)){
    p3tv <- length(beta3tv.true)
    if(h3tv_degree == "linear"){ #linear
      stopifnot(p3tv==1)
      x3tv <- as.matrix(pmin(R,D,Cen))
      h3_knots <- c(0,Inf)
    } else if(h3tv_degree == "log1p") {
      stopifnot(p3tv==1)
      x3tv <- as.matrix(log1p(pmin(R,D,Cen)))
      h3_knots <- c(0,Inf)
    } else if(h3tv_degree == "cs"){ #cubic spline model
      #in cubic spline model, boundary knots are set directly at min/max endpoints,
      #so no need to fix at 0
      h3_quantile_seq <- seq(from = 0,to = 1, length.out = p3tv+1)
      h3_knots <- stats::quantile(pmin(R,Cen)[yesR==1 & R<Cen], h3_quantile_seq)
      x3tv <- splines::ns(x = pmin(R,D,Cen), knots = h3_knots[-c(1,length(h3_knots))],
                          Boundary.knots = h3_knots[c(1,length(h3_knots))],
                          intercept = FALSE)
    } else { #if we don't use restricted cubic, then we are using a regular b-spline with specified degree
      #this also includes piecewise constant if degree is 0
      h3tv_degree <- as.numeric(h3tv_degree)
      stopifnot(p3tv>=h3tv_degree)
      h3_quantile_seq <- seq(from = 0,to = 1, length.out = p3tv+2-h3tv_degree)[-c(1,p3tv+2-h3tv_degree)]
      #fixing piecewise endpoint at maximum is ok, because splines2 prediction will extrapolate beyond it
      h3_knots <- c(0,stats::quantile(pmin(R,Cen)[yesR==1 & R<Cen],
                                      h3_quantile_seq),max(pmin(R,D,Cen)))
      x3tv <- splines2::bSpline(x = pmin(R,D,Cen), intercept = FALSE, degree = h3tv_degree,
                                knots = h3_knots[-c(1,length(h3_knots))],
                                Boundary.knots = h3_knots[c(1,length(h3_knots))])
    }
    colnames(x3tv) <- paste0("h3tv",1:p3tv)
    LP3 <- LP3 + x3tv %*% beta3tv.true
  } else{
    p3tv <- 0
  }

  #check the lower.bound stuff to make sure that I have the right probs and not 1-probs
  if(tolower(model) %in% c("markov","m")){
    if(tolower(hazard) %in% c("weibull","wb")){
      M_bound <- stats::pweibull(R[yesR],lower.tail = FALSE, shape = alpha3.true,
                        scale = exp(-(log(kappa3.true) + LP3[yesR] +
                                      beta3frail.true * log(gamma.true[yesR]))/alpha3.true),)
      M_prob <- stats::runif(sum(yesR),min=0, max=M_bound)
      D[yesR] <- stats::qweibull(p = M_prob, lower.tail = FALSE, shape = alpha3.true,
                  scale = exp(-(log(kappa3.true) + LP3[yesR] +
                                beta3frail.true * log(gamma.true[yesR]))/alpha3.true))
    } else{
      M_bound <- ppwexp(R[yesR], lower.tail = FALSE, log.p=FALSE,
                        phi=phi3.true, knots_vec=knots03,
                        eta = LP3[yesR] + beta3frail.true * log(gamma.true[yesR]))
      M_prob <- stats::runif(sum(yesR),min=0, max=M_bound)
      D[yesR] <- qpwexp(M_prob, lower.tail = FALSE, phi=phi3.true,
                        knots_vec=knots03, eta = LP3[yesR] + beta3frail.true * log(gamma.true[yesR]))
    }
  } else{
    if(tolower(hazard) %in% c("weibull","wb")){
      D[yesR] <- R[yesR] + stats::rweibull(sum(yesR), shape = alpha3.true,
                            scale = exp(-(log(kappa3.true) + LP3[yesR] +
                                          beta3frail.true * log(gamma.true[yesR]))/alpha3.true))
    } else {
      D[yesR] <- R[yesR] + qpwexp(stats::runif(sum(yesR)), lower.tail = FALSE,
                                  phi=phi3.true, knots_vec=knots03,
                                  eta =  LP3[yesR] + beta3frail.true * log(gamma.true[yesR]))
    }
  }
  delta1 <- rep(NA, n)
  delta2 <- rep(NA, n)
  y1 <- R
  y2 <- D

  if(cens[2] == 0){
    Cen <- rep(Inf,n)
  } else{
    #fix what to do about intersection of left truncation and censoring interval!

    # if(cens[1] < LT_interval[2]){
    #   warning(paste0("Censoring distribution cannot overlap truncation distribution.",
    #                  "Setting minimum censoring time to ", LT_interval[2],"."))
    #   cens[1] <- LT_interval[2]
    # }
    Cen <- stats::runif(n, cens[1], cens[2])
  }

  #cases where terminal occurs before non-terminal and censoring
  ind01 <- which(D < R & D < Cen)
  y1[ind01] <- D[ind01]
  delta1[ind01] <- 0
  delta2[ind01] <- 1

  #cases where nonterminal occurs, then censoring before terminal
  ind10 <- which(R < D & R < Cen & D >= Cen)
  y2[ind10] <- Cen[ind10]
  delta1[ind10] <- 1
  delta2[ind10] <- 0

  #cases where censoring occurs first
  ind00 <- which(R >= Cen & D >= Cen)
  y1[ind00] <- Cen[ind00]
  y2[ind00] <- Cen[ind00]
  delta1[ind00] <- 0
  delta2[ind00] <- 0

  #cases where nonterminal occurs, then terminal, then censoring
  ind11 <- which(R < Cen & D < Cen & R < D)
  delta1[ind11] <- 1
  delta2[ind11] <- 1

  #this is a gut-check that the values I build the basis for h3tv on above
  #match the values of y1 for all observations that matter (e.g., those with delta1==1)
  stopifnot(all((y1==pmin(R,Cen))[delta1==1]))

  ret <- data.frame(cbind(y1, delta1, y2, delta2, yL, gamma.true,
                          if(p3tv > 0) x3tv,
                          id, Vmat))
  if(!is.null(beta3tv.true)){
    attr(ret,which = "p3tv") <- p3tv
    attr(ret,which = "h3tv_degree") <- h3tv_degree
    attr(ret,which = "h3tv_knots") <- h3_knots
  }

  # if(!(tolower(hazard) %in% c("weibull","wb"))){
  #   attr(ret,which = "knots01") <- knots01
  #   attr(ret,which = "knots02") <- knots02
  #   attr(ret,which = "knots03") <- knots03
  # }

  return(ret)
}
