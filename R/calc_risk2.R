#' Calculate absolute risk profiles
#'
#' This function calculates absolute risk profiles under weibull baseline hazard specifications.
#'
#' @inheritParams nll_func
#' @inheritParams FreqID_HReg2
#' @param t_cutoff Numeric vector indicating the time(s) to compute the risk profile.
#' @param t_start Numeric scalar indicating the dynamic start time to compute the risk profile. Set to 0 by default.
#' @param tol Numeric value for the tolerance of the numerical integration procedure.
#' @param type String either indicating 'marginal' for population-averaged probabilities,
#'   or 'conditional' for probabilities computed at the specified gamma
#' @param gamma Numeric value indicating the fixed level of the frailty assumed for predicted probabilities,
#'   if 'type' is set to 'conditional'
#' @param h3_tv String indicating whether there is an effect of t1 on hazard 3.
#' @param tv_knots for piecewise effect of t1 in h3, these are the knots at which the effect jumps
#'
#' @return if Xmat has only one row, and t_cutoff is a scalar, then returns a 4 element row matrix
#'   of probabilities. If Xmat has \code{n} rows, then returns an \code{n} by 4 matrix of probabilities.
#'   If Xmat has \code{n} rows and t_cutoff is a vector of length \code{s}, then returns an \code{s} by 4 by \code{n} array.
#' @export
calc_risk2 <- function(para, Xmat1, Xmat2, Xmat3,hazard,knots_list=NULL,
                      t_cutoff, t_start=0, tol=1e-3, frailty=TRUE,
                      type="marginal", gamma=1, model="semi-markov",
                      h3_tv="none",tv_knots=NULL,
                      quad_method="legendre",n_quad=15){
  # browser()

  n <- max(1,nrow(Xmat1),nrow(Xmat2),nrow(Xmat3))
  t_length <- length(t_cutoff)
  #standardize namings for the use of "switch" below
  stopifnot(tolower(hazard) %in% c("wb","weibull","pw","piecewise",
                                   "bs","bspline","rp","royston-parmar"))
  hazard <- switch(tolower(hazard),
                   wb="weibull",weibull="weibull",
                   pw="piecewise",piecewise="piecewise",
                   bs="bspline",bspline="bspline",
                   rp="rp","royston-parmar"="rp")
  stopifnot(tolower(type) %in% c("c","conditional","m","marginal"))
  type <- switch(tolower(type),
                 c="conditional",conditional="conditional",
                 m="marginal",marginal="marginal")
  stopifnot(tolower(model) %in% c("sm","semi-markov","m","markov"))
  model <- switch(tolower(model),
                  sm="semi-markov","semi-markov"="semi-markov",
                  m="markov",markov="markov")

  if(hazard == "bspline"){
    quad_weights <- get_quad_pointsweights(n_quad=n_quad,
                                           quad_method=quad_method)$weights
  }


  #set up parameters
  if(hazard == "weibull"){
    nP01 <- nP02 <- nP03 <- 2
  } else{
    stopifnot(!is.null(knots_list))

    if(hazard!="rp"){
      #left pad with a zero if it is not already present
      if(knots_list[[1]][1] != 0){knots_list[[1]] <- c(0,knots_list[[1]])}
      if(knots_list[[2]][1] != 0){knots_list[[2]] <- c(0,knots_list[[2]])}
      if(knots_list[[3]][1] != 0){knots_list[[3]] <- c(0,knots_list[[3]])}
    }
    knots1 <- knots_list[[1]]
    knots2 <- knots_list[[2]]
    knots3 <- knots_list[[3]]

    if(hazard!="bspline"){ #pw and rp have knots vec of length nP
      nP01 <- length(knots1)
      nP02 <- length(knots2)
      nP03 <- length(knots3)
    } else{ #bs has knots vec of length nP - 2
      nP01 <- length(knots1) + 2
      nP02 <- length(knots2) + 2
      nP03 <- length(knots3) + 2
    }
  }

  #define baseline parameter subvectors
  phi1 <- para[1:nP01]
  phi2 <- para[(1+nP01):(nP01+nP02)]
  phi3 <- para[(1+nP01+nP02):(nP01+nP02+nP03)]

  if(frailty){
    nP0 <- nP01 + nP02 + nP03 + 1
    theta <- exp(para[nP0])
    if(type=="conditional" & length(gamma)==1){
      gamma <- rep(gamma,n)
    }
  } else{
    nP0 <- nP01 + nP02 +nP03
    type <- "conditional"
    gamma <- rep(1,n)
  }

  if(!is.null(Xmat1) && !(ncol(Xmat1)==0)){
    nP1 <- ncol(Xmat1)
    beta1 <- para[(1+nP0):(nP0+nP1)]
    eta1 <- as.vector(Xmat1 %*% beta1)
  } else{
    nP1 <- 0; eta1 <- 0
  }
  if(!is.null(Xmat2) && !(ncol(Xmat2)==0)){
    nP2 <- ncol(Xmat2)
    beta2 <- para[(1+nP0+nP1):(nP0+nP1+nP2)]
    eta2 <- as.vector(Xmat2 %*% beta2)
  } else{
    nP2 <- 0; eta2 <- 0
  }
  if(!is.null(Xmat3) && !(ncol(Xmat3)==0)){
    nP3 <- ncol(Xmat3)
    beta3 <- para[(1+nP0+nP1+nP2):(nP0+nP1+nP2+nP3)]
    eta3 <- as.vector(Xmat3 %*% beta3)
  } else{
    nP3 <- 0; eta3 <- 0
  }

  #DO MORE WITH THIS!!
  #specify different forms by which t1 can be incorporated into h3
  if(tolower(h3_tv) == "linear"){
    if(tolower(model) != "semi-markov"){stop("must be semi-markov to have t1 in h3.")}
    n_tv <- 1
    beta3_tv_linear <- utils::tail(para,n = n_tv)
  } else if(tolower(h3_tv) %in% c("pw","piecewise")){
    if(tolower(model) != "semi-markov"){stop("must be semi-markov to have t1 in h3.")}
    stopifnot(!is.null(tv_knots))
    if(tv_knots[1] != 0){tv_knots <- c(0,tv_knots)}
    if(utils::tail(tv_knots, n=1) != Inf){tv_knots <- c(tv_knots,Inf)}
    n_tv <- length(tv_knots) - 2
    beta3_tv <- c(0,utils::tail(para,n=n_tv))
    beta3_tv_linear <- 0
  } else{
    n_tv <- 0
    beta3_tv_linear <- 0
  }

  #if the size of the parameter vector doesn't match the expected size, throw a fuss
  stopifnot(length(para) == nP0 + nP1 + nP2 + nP3 + n_tv)

  #NEW PLAN
  #1. first, define each "integrand" function for each baseline hazard specification for marginal + conditional
  #   do this rather than specifying hazard and Hazard functions and building up, because
  #   it might be ever so slightly more efficient / straightforward
  #2. either adapt vectorized integration functions from below, or implement numerical integration a la bspline function
  #gamma (vector) and theta (scalar) are implicitly defined

  p_neither_func <- switch(hazard,
    weibull=function(t){
      if(type == "marginal"){
        (1 + theta*(t^exp(phi1[2]) * exp(phi1[1] + eta1) +
                      t^exp(phi2[2]) * exp(phi2[1] + eta2)))^(-theta^(-1))
      } else{
        exp(-gamma*(t^exp(phi1[2]) * exp(phi1[1] + eta1) +
                             t^exp(phi2[2]) * exp(phi2[1] + eta2)))
      }
    },
    piecewise=function(t){
      basis1 <- get_basis(y = t,knots_vec = knots1,hazard = "piecewise",deriv = FALSE)
      basis2 <- get_basis(y = t,knots_vec = knots2,hazard = "piecewise",deriv = FALSE)
      Lambda01 <- as.vector(basis1 %*% exp(phi1))
      Lambda02 <- as.vector(basis2 %*% exp(phi2))
      if(type == "marginal"){
        (1 + theta*(Lambda01 * exp(eta1) + Lambda02 * exp(eta2)))^(-theta^(-1))
      } else{
        exp(-gamma*(Lambda01 * exp(eta1) + Lambda02 * exp(eta2)))
      }
    },
    rp=function(t){
      basis1 <- get_basis(y = t,knots_vec = knots1,hazard = "royston-parmar",deriv = FALSE)
      basis2 <- get_basis(y = t,knots_vec = knots2,hazard = "royston-parmar",deriv = FALSE)
      s1 <- as.vector(basis1 %*% phi1)
      s2 <- as.vector(basis2 %*% phi2)
      if(type == "marginal"){
        (1 + theta*(exp(s1 + eta1) + exp(s2 + eta2)))^(-theta^(-1))
      } else{
        exp(-gamma*(exp(s1 + eta1) + exp(s2 + eta2)))
      }
    },
    bspline=function(t){
      quad_points <- transform_quad_points(n_quad = n_quad,
                                           quad_method=quad_method, a = 0,b = t)
      basis1_quad <- get_basis(y=quad_points, knots_vec=knots1,hazard="bspline")
      basis2_quad <- get_basis(y=quad_points, knots_vec=knots2,hazard="bspline")
      lambda01 <- as.vector(exp(basis1_quad %*% phi1))
      lambda02 <- as.vector(exp(basis2_quad %*% phi2))
      #reshape lambda0 from a n*n_quad length vector
      #to an n by n_quad matrix, then multiply with n_quad length weights
      #to get final Lambda0
      Lambda01 <- t/2 * as.vector(matrix(lambda01,ncol=n_quad,byrow = TRUE) %*% quad_weights)
      Lambda02 <- t/2 * as.vector(matrix(lambda02,ncol=n_quad,byrow = TRUE) %*% quad_weights)
      if(type == "marginal"){
        (1 + theta*(Lambda01 * exp(eta1) + Lambda02 * exp(eta2)))^(-theta^(-1))
      } else{
        exp(-gamma*(Lambda01 * exp(eta1) + Lambda02 * exp(eta2)))
      }
    })

  t2_only_integrand <- switch(hazard,
     weibull=function(t,index){
       if(type == "marginal"){
         exp(phi2[1] + phi2[2] + eta2[index]) * t^expm1(phi2[2]) * #h2
           (1 + theta*(t^exp(phi1[2]) * exp(phi1[1] + eta1[index]) +
                         t^exp(phi2[2]) * exp(phi2[1] + eta2[index])))^(-theta^(-1)) #S12
       } else{
         gamma[index] * exp(phi2[1] + phi2[2] + eta2[index]) * t^expm1(phi2[2]) * #h2
           exp(-gamma[index]*(t^exp(phi1[2]) * exp(phi1[1] + eta1[index]) +
                                t^exp(phi2[2]) * exp(phi2[1] + eta2[index]))) #S12
       }
     },
     piecewise=function(t,index){
       basis1 <- get_basis(y = t,knots_vec = knots1,hazard = "piecewise",deriv = FALSE)
       basis2 <- get_basis(y = t,knots_vec = knots2,hazard = "piecewise",deriv = FALSE)
       Lambda01 <- as.vector(basis1 %*% exp(phi1))
       Lambda02 <- as.vector(basis2 %*% exp(phi2))
       dbasis2 <- get_basis(y = t,knots_vec = knots2,hazard = "piecewise",deriv = TRUE)
       lambda02 <- as.vector(dbasis2 %*% exp(phi2))
       if(type == "marginal"){
         lambda02 * exp(eta2[index]) * #h2
           (1 + theta*(Lambda01 * exp(eta1[index]) + Lambda02 * exp(eta2[index])))^(-theta^(-1)) #S12
       } else{
         gamma[index] * lambda02 * exp(eta2[index]) * #h2
           exp(-gamma[index]*(Lambda01 * exp(eta1[index]) + Lambda02 * exp(eta2[index]))) #S12
       }
     },
     rp=function(t,index){
       basis1 <- get_basis(y = t,knots_vec = knots1,hazard = "royston-parmar",deriv = FALSE)
       basis2 <- get_basis(y = t,knots_vec = knots2,hazard = "royston-parmar",deriv = FALSE)
       dbasis2 <- get_basis(y = t,knots_vec = knots2,hazard = "royston-parmar",deriv = TRUE)
       s1 <- as.vector(basis1 %*% phi1)
       s2 <- as.vector(basis2 %*% phi2)
       s2prime <- as.vector(dbasis2 %*% phi2)
       if(type == "marginal"){
         s2prime * exp(s2 + eta2[index]) / t * #h2
           (1 + theta*(exp(s1 + eta1[index]) + exp(s2 + eta2[index])))^(-theta^(-1)) #S12
       } else{
         gamma[index] * s2prime * exp(s2 + eta2[index]) / t * #h2
           exp(-gamma[index]*(exp(s1 + eta1[index]) + exp(s2 + eta2[index]))) #S12
       }
     },
     bspline=function(t,index){
       quad_points <- transform_quad_points(n_quad = n_quad,
                                            quad_method=quad_method, a = 0,b = t)
       basis2 <- get_basis(y=t, knots_vec=knots2,hazard="bspline")
       basis1_quad <- get_basis(y=quad_points, knots_vec=knots1,hazard="bspline")
       basis2_quad <- get_basis(y=quad_points, knots_vec=knots2,hazard="bspline")
       lambda01 <- as.vector(exp(basis1_quad %*% phi1))
       lambda02 <- as.vector(exp(basis2_quad %*% phi2))
       #reshape lambda0 from a n_t*n_quad length vector
       #to an n_t by n_quad matrix, then multiply with n_quad length weights
       #to get final Lambda0
       Lambda01 <- t/2 * as.vector(matrix(lambda01,ncol=n_quad,byrow = TRUE) %*% quad_weights)
       Lambda02 <- t/2 * as.vector(matrix(lambda02,ncol=n_quad,byrow = TRUE) %*% quad_weights)
       if(type == "marginal"){
         as.vector(exp(basis2 %*% phi2 + eta2[index])) *
          (1 + theta*(Lambda01 * exp(eta1[index]) + Lambda02 * exp(eta2[index])))^(-theta^(-1))
       } else{
         gamma[index] * as.vector(exp(basis2 %*% phi2 + eta2[index])) *
           exp(-gamma[index]*(Lambda01 * exp(eta1[index]) + Lambda02 * exp(eta2[index])))
       }
     })

  #this one is trickier because it depends on semi-markov vs markov, and effect of t1 on h3
  t1_only_integrand <- switch(hazard,
    weibull=function(t,t_bound,index){
      if(model=="markov"){
        Lambda03 <- (t_bound^exp(phi3[2]) - t^exp(phi3[2])) * exp(phi3[1])
      } else{
        Lambda03 <- (t_bound-t)^exp(phi3[2]) * exp(phi3[1])
      }
      if(type == "marginal"){
        exp(phi1[1] + phi1[2] + eta1[index]) * t^expm1(phi1[2]) * #h2
          (1 + theta*(t^exp(phi1[2]) * exp(phi1[1] + eta1[index]) +
                        t^exp(phi2[2]) * exp(phi2[1] + eta2[index]) +
                        Lambda03 * exp(eta3[index])))^(-theta^(-1)) #S12
      } else{
        gamma[index] * exp(phi1[1] + phi1[2] + eta1[index]) * t^expm1(phi1[2]) * #h2
          exp(-gamma[index]*(t^exp(phi1[2]) * exp(phi1[1] + eta1[index]) +
                               t^exp(phi2[2]) * exp(phi2[1] + eta2[index]) +
                               Lambda03 * exp(eta3[index]))) #S12
      }
    },
    piecewise=function(t,t_bound,index){
      basis1 <- get_basis(y = t,knots_vec = knots1,hazard = "piecewise",deriv = FALSE)
      basis2 <- get_basis(y = t,knots_vec = knots2,hazard = "piecewise",deriv = FALSE)
      Lambda01 <- as.vector(basis1 %*% exp(phi1))
      Lambda02 <- as.vector(basis2 %*% exp(phi2))
      dbasis1 <- get_basis(y = t,knots_vec = knots1,hazard = "piecewise",deriv = TRUE)
      lambda01 <- as.vector(dbasis1 %*% exp(phi1))

      if(model=="markov"){
        basis3 <- get_basis(y = t_bound,knots_vec = knots3,hazard = "piecewise",deriv = FALSE)
        Lambda03 <- as.vector(basis3 %*% exp(phi3))
        basis3 <- get_basis(y = t,knots_vec = knots3,hazard = "piecewise",deriv = FALSE)
        Lambda03 <- Lambda03 - as.vector(basis3 %*% exp(phi3))
      } else{
        basis3 <- get_basis(y = t_bound-t,knots_vec = knots3,hazard = "piecewise",deriv = FALSE)
        Lambda03 <- as.vector(basis3 %*% exp(phi3))
      }

      if(type == "marginal"){
        lambda01 * exp(eta1[index]) * #h2
          (1 + theta*(Lambda01 * exp(eta1[index]) +
                        Lambda02 * exp(eta2[index]) +
                        Lambda03 * exp(eta3[index])))^(-theta^(-1)) #S12
      } else{
        gamma[index] * lambda01 * exp(eta1[index]) * #h2
          exp(-gamma[index]*(Lambda01 * exp(eta1[index]) +
                               Lambda02 * exp(eta2[index]) +
                               Lambda03 * exp(eta3[index]))) #S12
      }
    },
    rp=function(t,t_bound,index){
      basis1 <- get_basis(y = t,knots_vec = knots1,hazard = "royston-parmar",deriv = FALSE)
      basis2 <- get_basis(y = t,knots_vec = knots2,hazard = "royston-parmar",deriv = FALSE)
      dbasis1 <- get_basis(y = t,knots_vec = knots1,hazard = "royston-parmar",deriv = TRUE)
      s1 <- as.vector(basis1 %*% phi1)
      s2 <- as.vector(basis2 %*% phi2)
      s1prime <- as.vector(dbasis1 %*% phi1)

      if(model=="markov"){
        basis3 <- get_basis(y = t_bound,knots_vec = knots3,hazard = "royston-parmar",deriv = FALSE)
        Lambda03 <- exp(as.vector(basis3 %*% phi3))
        basis3 <- get_basis(y = t,knots_vec = knots3,hazard = "royston-parmar",deriv = FALSE)
        Lambda03 <- Lambda03 - exp(as.vector(basis3 %*% phi3))
      } else{
        basis3 <- get_basis(y = t_bound-t,knots_vec = knots3,hazard = "royston-parmar",deriv = FALSE)
        Lambda03 <- exp(as.vector(basis3 %*% phi3))
      }
      if(type == "marginal"){
        s1prime * exp(s1 + eta1[index]) / t * #h1
          (1 + theta*(exp(s1 + eta1[index]) +
                        exp(s2 + eta2[index]) +
                        Lambda03*exp(eta3[index])))^(-theta^(-1))
      } else{
        gamma[index] * s1prime * exp(s1 + eta1[index]) / t * #h1
          exp(-gamma[index]*(exp(s1 + eta1[index]) +
                               exp(s2 + eta2[index]) +
                               Lambda03*exp(eta3[index])))
      }
    },
    bspline=function(t,t_bound,index){
      quad_points <- transform_quad_points(n_quad = n_quad,
                                           quad_method=quad_method, a = 0,b = t)
      basis1 <- get_basis(y=t, knots_vec=knots1,hazard="bspline")
      basis1_quad <- get_basis(y=quad_points, knots_vec=knots1,hazard="bspline")
      basis2_quad <- get_basis(y=quad_points, knots_vec=knots2,hazard="bspline")
      lambda01 <- as.vector(exp(basis1_quad %*% phi1))
      lambda02 <- as.vector(exp(basis2_quad %*% phi2))
      #reshape lambda0 from a n_t*n_quad length vector
      #to an n_t by n_quad matrix, then multiply with n_quad length weights
      #to get final Lambda0
      Lambda01 <- t/2 * as.vector(matrix(lambda01,ncol=n_quad,byrow = TRUE) %*% quad_weights)
      Lambda02 <- t/2 * as.vector(matrix(lambda02,ncol=n_quad,byrow = TRUE) %*% quad_weights)

      if(model=="markov"){
        quad_points <- transform_quad_points(n_quad = n_quad,
                                             quad_method=quad_method, a = t,b = t_bound)
        basis3_quad <- get_basis(y=quad_points, knots_vec=knots3,hazard="bspline")
      } else{
        quad_points <- transform_quad_points(n_quad = n_quad,
                                             quad_method=quad_method, a = 0,b = t_bound-t)
        basis3_quad <- get_basis(y=quad_points, knots_vec=knots3,hazard="bspline")
      }
      lambda03 <- as.vector(exp(basis3_quad %*% phi3))
      Lambda03 <- (t_bound-t)/2 * as.vector(matrix(lambda03,ncol=n_quad,byrow = TRUE) %*% quad_weights)

      if(type == "marginal"){
        as.vector(exp(basis1 %*% phi1 + eta1[index])) *
          (1 + theta*(Lambda01 * exp(eta1[index]) +
                        Lambda02 * exp(eta2[index]) +
                        Lambda03 * exp(eta3[index])))^(-theta^(-1))
      } else{
        gamma[index] * as.vector(exp(basis1 %*% phi1 + eta1[index])) *
          exp(-gamma[index]*(Lambda01 * exp(eta1[index]) +
                               Lambda02 * exp(eta2[index]) +
                               Lambda03 * exp(eta3[index])))
      }
    })

  ##ACTUALLY COMPUTING THE PROBABILITIES##
  ##************************************##

  #If we are computing 'dynamic probabilities' updated to some later timepoint t_start,
  #then we need to compute the probability of having experienced the event by t_start.
  #This is easier for 'neither' and 'tonly' outcomes because the inner integral does not depend on t_start.
  if(t_start > 0){
    p_neither_start <- p_neither_func(t=t_start,type=type)
    p_tonly_start <- sapply(1:n,function(x){
      tryCatch(stats::integrate(t2_only_integrand, index=x,
                                lower=0, upper=t_start)$value,
                error=function(cnd){return(NA)}) })
  } else{
    p_tonly_start <- 0
    p_neither_start <- 1
  }

  #this function allows inputs with multiple subjects, multiple time points, or both
  #therefore, we need to create the right data structure to contain the output.
  if(n > 1){
    if(t_length > 1){
      out_mat <- array(dim=c(t_length,4,
                             n),dimnames = list(paste0("t",t_cutoff),
                                                c("p_ntonly","p_both","p_tonly","p_neither"),paste0("i",1:n)))
    } else{
      out_mat <- matrix(nrow=n,ncol=4,
                        dimnames = list(paste0("i",1:n),
                                        c("p_ntonly","p_both","p_tonly","p_neither")))
    }
  } else{
    out_mat <- matrix(nrow=t_length,ncol=4,
                      dimnames = list(paste0("t",t_cutoff),
                                      c("p_ntonly","p_both","p_tonly","p_neither")))
  }

  #loop through each time point, and compute the predicted probability at that time point for all subjects
  #each probability is an n-length vector
  for(t_ind in 1:t_length){
    t_temp <- t_cutoff[t_ind]

    p_ntonly <- sapply(1:n,function(x){
      tryCatch(stats::integrate(t1_only_integrand, index=x, t_bound = t_temp,
                                lower=0, upper=t_temp)$value,
               error=function(cnd){return(NA)}) })
    if(t_start > 0){
      p_ntonly_start <- sapply(1:n,function(x){
        tryCatch(stats::integrate(t1_only_integrand, index=x, t_bound = t_temp,
                                  lower=0, upper=t_start)$value,
                 error=function(cnd){return(NA)}) })
      p_both_start <- 1 - p_tonly_start - p_ntonly_start - p_neither_start
    } else{
      p_both_start <- p_ntonly_start <- 0
    }

    p_tonly <- sapply(1:n,function(x){
      tryCatch(stats::integrate(t2_only_integrand, index=x,
                                lower=0, upper=t_temp)$value,
               error=function(cnd){return(NA)}) })
    p_neither <- p_neither_func(t=t_temp)
    p_both <- 1 - p_tonly - p_neither - p_ntonly

    out_temp <- cbind(p_ntonly=(p_ntonly-p_ntonly_start)/p_neither_start,
                      p_both=(p_both-p_both_start)/p_neither_start,
                      p_tonly=(p_tonly-p_tonly_start)/p_neither_start,
                      p_neither=(p_neither)/p_neither_start)

    if(n > 1){
      if(t_length > 1){
        out_mat[t_ind,,] <- t(out_temp)
      } else{
        out_mat <- out_temp
      }
    } else{
      out_mat[t_ind,] <- out_temp
    }
  }

  return(out_mat)
}







#' Get matrix of observed outcome categories
#'
#' This function returns a matrix giving the observed outcome categories of each observation at various
#'   time cutoffs.
#'
#' @inheritParams nll_func
#' @inheritParams calc_risk2
#'
#' @return a matrix or array.
#' @export
get_outcome_mat <- function(y1, y2, delta1, delta2, t_cutoff){

  n <- length(y1)
  t_length <- length(t_cutoff)

  if(n > 1){
    if(t_length > 1){
      out_mat <- array(dim=c(t_length,4,n),dimnames = list(paste0("t",t_cutoff),c("ntonly","both","tonly","neither"),paste0("i",1:n)))
    } else{
      out_mat <- matrix(nrow=n,ncol=4,dimnames = list(paste0("i",1:n),c("ntonly","both","tonly","neither")))
    }
  } else{
    out_mat <- matrix(nrow=t_length,ncol=4,dimnames = list(paste0("t",t_cutoff),c("ntonly","both","tonly","neither")))
  }

  for(t_ind in 1:t_length){

    #For cases where y=t_cutoff, I consider events that happened exactly at t_cutoff in categorization.
    neither <- t_cutoff[t_ind] < y1 | #neither
      y2 <= t_cutoff[t_ind] & delta1 == 0 & delta2 == 0 #neither
    ntonly <- y1 <= t_cutoff[t_ind] & t_cutoff[t_ind] < y2 | #ntonly
      y2 <= t_cutoff[t_ind] & delta1 == 1 & delta2 == 0 #ntonly
    tonly <- y2 <= t_cutoff[t_ind] & delta1 == 0 & delta2 == 1 #tonly
    both <- y2 <= t_cutoff[t_ind] & delta1 == 1 & delta2 == 1 #both



    out_temp <- cbind(ntonly=as.numeric(ntonly),
                      both=as.numeric(both),
                      tonly=as.numeric(tonly),
                      neither=as.numeric(neither))

    if(n > 1){
      if(t_length > 1){
        out_mat[t_ind,,] <- t(out_temp)
      } else{
        out_mat <- out_temp
      }
    } else{
      out_mat[t_ind,] <- out_temp
    }
  }
  out_mat
}

#' Get inverse probability of censoring weights
#'
#' This function returns a vector of inverse probability of censoring weights from an unadjusted Cox model
#'   for censoring times.
#'
#' @inheritParams nll_func
#' @inheritParams calc_risk2
#'
#' @return a vector.
#' @export
get_ipcw_mat <- function(y2,delta2,t_cutoff){

  # browser()
  n <- length(y2)
  t_length <- length(t_cutoff)

  #this is Ghat, a non-parametric model of the 'survival distribution' of censoring var C.
  sfcens <- survival::survfit(survival::Surv(y2, delta2==0) ~ 1)

  #* want to get Ghat(z) where (following Graf 1999 'three categories')
  #* Category 1:
  #* z=y2- if y2<=s and delta2=1,
  #* Category 2:
  #* z=s if s<y2
  #* Category 3:
  #* z=Inf if y2<s and delta2=0, (aka, 1/Ghat(z)=0, I know they're subtly different)
  #* z=s if s=y2 and delta2=0, #this one situation is what I'm unsure about
  #* because graf and spitoni would say that if s=y2 and delta2=0, then z=Inf and result should be tossed.
  #* below I defer to them, but I'm just noting that to me there is logic in the other direction.
  #*
  #* so, we define
  #* z = min(y2,s)
  #* then change z=y2- if y2<=s and delta2=1
  #* then change z=Inf if y2<=s and delta2=0 (again, technically change 1/Ghat(z)=0 for these obs but still)
  #*
  #* then change z=Inf if y2< s and delta2=0 (again, technically change 1/Ghat(z)=0 for these obs but still)
  #* and that should do it! (I hope)

  ipcw_mat <- matrix(nrow=n,ncol=t_length,dimnames = list(paste0("i",1:n),paste0("t",t_cutoff)))
  for(t_ind in 1:t_length){
    #vector of min(ttilde,t)
    y_last <- pmin(y2,t_cutoff[t_ind])
    if(sum(y2 <= t_cutoff[t_ind] & delta2==1)>0){
      y_last[y2 <= t_cutoff[t_ind] & delta2==1] <- y_last[y2 <= t_cutoff[t_ind] & delta2==1] - 1e-8
    }
    y_last_cens <- rep(NA,n)
    y_last_cens[order(y_last)] <- summary(sfcens, times = y_last, extend=TRUE)$surv

    #now, to eliminate the the 'censored' observations
    if(sum(y_last <= t_cutoff[t_ind] & delta2==0) > 0){
      y_last_cens[y2 <= t_cutoff[t_ind] & delta2==0] <- Inf
    }
    #below is my former definition, but I will defer to graf and spitoni above
    # if(sum(y_last < t_cutoff[t_ind] & delta2==0) > 0){
    #   y_last_cens[y2 < t_cutoff[t_ind] & delta2==0] <- Inf
    # }
    ipcw_mat[,t_ind] <- 1/y_last_cens
  }
  ipcw_mat
}


#' Compute prediction performance score
#'
#' This function takes in all of the ingredients needed for prediction validation,
#'   and returns the corresponding scores.
#'
#' @param outcome_mat Output from get_outcome_mat function
#' @param pred_mat Output from calc_risks function
#' @param ipcw_mat Output from get_ipcw_mat function
#' @param score String indicating whether 'brier' score, or 'entropy' should be computed.
#'
#' @return a vector.
#' @export
compute_score <- function(outcome_mat, pred_mat, ipcw_mat, score="brier"){
  #this function is for brier and kl scores (aka cross entropy)
  #one weird thing is that in spitoni, the authors divide by n instead of by the sum of weights, is that right?
  # browser()
  if(length(dim(outcome_mat))==3){
    if(tolower(score) %in% c("brier")){
      out <- apply( t(apply((outcome_mat - pred_mat)^2,
                            MARGIN = c(1,3),
                            FUN = sum)) *
                      ipcw_mat, MARGIN = 2, FUN = mean)
    } else{
      out <- apply( t(apply(-outcome_mat*log(pred_mat),
                            MARGIN = c(1,3),
                            FUN = sum)) *
                      ipcw_mat, MARGIN = 2, FUN = mean)
    }
  } else{ #this must mean that there is only a single time cutoff, so the input mats are n by 4 and weights is just a vector
    if(tolower(score) %in% c("brier")){
      out <- colMeans( apply((outcome_mat - pred_mat)^2,
                             MARGIN = c(1),
                             FUN = sum) * ipcw_mat)
    } else{
      out <- colMeans( apply(-outcome_mat*log(pred_mat),
                             MARGIN = c(1),
                             FUN = sum) * ipcw_mat)
    }
  }
  out
}



