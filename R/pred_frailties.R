#' Predict Gamma Frailties
#'
#' This function predicts estimated gamma frailties via empirical bayes
#'
#' @inheritParams nll_func
#' @inheritParams FreqID_HReg2
#'
#' @return if Xmat has only one row, and t_cutoff is a scalar, then returns a 1 by 4 element row matrix
#'   of the following values:
#'   \itemize{
#'     \item{\code{mean} of subject's gamma empirical bayes posterior, i.e., shape/rate.)}
#'     \item{\code{median} of subject's gamma empirical bayes posterior.)}
#'     \item{\code{shape} parameter (as used in, e.g., \code{dgamma}, of subject's gamma empirical bayes posterior.)}
#'     \item{\code{rate} parameter (as used in, e.g., \code{dgamma}, of subject's gamma empirical bayes posterior.)}
#'   }
#'   If Xmat has \code{n} rows, then returns an \code{n} by 4 matrix.
#'   If Xmat has \code{n} rows and t_cutoff is a vector of length \code{s}, then returns an \code{s} by 4 by \code{n} array.
#' @export
pred_frailties <- function(para, y1, delta1, y2, delta2, yL=NULL,
                           Xmat1, Xmat2, Xmat3,
                           hazard, knots_list=NULL, model,
                           quad_method="legendre",n_quad=15){
  # browser()

  n <- max(1,nrow(Xmat1),nrow(Xmat2),nrow(Xmat3))
  anyLT <- as.numeric(any(yL>0))

  #standardize namings for the use of "switch" below
  stopifnot(tolower(hazard) %in% c("wb","weibull","pw","piecewise",
                                   "bs","bspline","rp","royston-parmar"))
  hazard <- switch(tolower(hazard),
                   wb="weibull",weibull="weibull",
                   pw="piecewise",piecewise="piecewise",
                   bs="bspline",bspline="bspline",
                   rp="rp","royston-parmar"="rp")
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

  nP0 <- nP01 + nP02 + nP03 + 1
  invtheta <- exp(-para[nP0])

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

  #Compute Lambda01, Lambda02 and Lambda03
  if(hazard %in% c("weibull","wb")){
    Lambda01 <- y1^exp(phi1[2]) * exp(phi1[1])
    Lambda02 <- y1^exp(phi2[2]) * exp(phi2[1])
    if(model=="markov"){
      Lambda03 <- (y2^exp(phi3[2]) - y1^exp(phi3[2])) * exp(phi3[1])
    } else{
      Lambda03 <- (y2-y1)^exp(phi3[2]) * exp(phi3[1])
    }
    if(anyLT){
      Lambda01_yL <- yL^exp(phi1[2]) * exp(phi1[1])
      Lambda02_yL <- yL^exp(phi2[2]) * exp(phi2[1])
    }
  } else if(hazard %in% c("piecewise","pw")) {
    basis1 <- get_basis(y = y1,knots_vec = knots1,hazard = "piecewise",deriv = FALSE)
    basis2 <- get_basis(y = y1,knots_vec = knots2,hazard = "piecewise",deriv = FALSE)
    Lambda01 <- as.vector(basis1 %*% exp(phi1))
    Lambda02 <- as.vector(basis2 %*% exp(phi2))
    if(model=="markov"){
      basis3 <- get_basis(y = y2,knots_vec = knots3,hazard = "piecewise",deriv = FALSE)
      Lambda03 <- as.vector(basis3 %*% exp(phi3))
      basis3 <- get_basis(y = y1,knots_vec = knots3,hazard = "piecewise",deriv = FALSE)
      Lambda03 <- Lambda03 - as.vector(basis3 %*% exp(phi3))
    } else{
      basis3 <- get_basis(y = y2-y1,knots_vec = knots3,hazard = "piecewise",deriv = FALSE)
      Lambda03 <- as.vector(basis3 %*% exp(phi3))
    }
    if(anyLT){
      #reuse the basis objects
      basis1 <- get_basis(y = yL,knots_vec = knots1,hazard = "piecewise",deriv = FALSE)
      basis2 <- get_basis(y = yL,knots_vec = knots2,hazard = "piecewise",deriv = FALSE)
      Lambda01_yL <- as.vector(basis1 %*% exp(phi1))
      Lambda02_yL <- as.vector(basis2 %*% exp(phi2))
    }
  } else if(hazard %in% c("royston-parmar","rp")){
    basis1 <- get_basis(y = y1,knots_vec = knots1,hazard = "royston-parmar",deriv = FALSE)
    basis2 <- get_basis(y = y1,knots_vec = knots2,hazard = "royston-parmar",deriv = FALSE)
    Lambda01 <- exp(as.vector(basis1 %*% phi1))
    Lambda02 <- exp(as.vector(basis2 %*% phi2))
    if(model=="markov"){
      basis3 <- get_basis(y = y2,knots_vec = knots3,hazard = "royston-parmar",deriv = FALSE)
      Lambda03 <- exp(as.vector(basis3 %*% phi3))
      basis3 <- get_basis(y = y1,knots_vec = knots3,hazard = "royston-parmar",deriv = FALSE)
      Lambda03 <- Lambda03 - exp(as.vector(basis3 %*% phi3))
    } else{
      basis3 <- get_basis(y = y2-y1,knots_vec = knots3,hazard = "royston-parmar",deriv = FALSE)
      Lambda03 <- exp(as.vector(basis3 %*% phi3))
    }
    if(anyLT){
      #reuse basis objects
      basis1 <- get_basis(y = yL,knots_vec = knots1,hazard = "royston-parmar",deriv = FALSE)
      basis2 <- get_basis(y = yL,knots_vec = knots2,hazard = "royston-parmar",deriv = FALSE)
      Lambda01_yL <- exp(as.vector(basis1 %*% phi1))
      Lambda02_yL <- exp(as.vector(basis2 %*% phi2))
    }
  } else if(hazard %in% c("bspline","bs")){
    quad_points <- transform_quad_points(n_quad = n_quad,
                                         quad_method=quad_method, a = 0,b = y1)
    basis1_quad <- get_basis(y=quad_points, knots_vec=knots1,hazard="bspline")
    basis2_quad <- get_basis(y=quad_points, knots_vec=knots2,hazard="bspline")
    lambda01 <- as.vector(exp(basis1_quad %*% phi1))
    lambda02 <- as.vector(exp(basis2_quad %*% phi2))
    #reshape lambda0 from a n_t*n_quad length vector
    #to an n_t by n_quad matrix, then multiply with n_quad length weights
    #to get final Lambda0
    Lambda01 <- y1/2 * as.vector(matrix(lambda01,ncol=n_quad,byrow = TRUE) %*% quad_weights)
    Lambda02 <- y1/2 * as.vector(matrix(lambda02,ncol=n_quad,byrow = TRUE) %*% quad_weights)
    if(model=="markov"){
      quad_points <- transform_quad_points(n_quad = n_quad,
                                           quad_method=quad_method, a = y1,b = y2)
      basis3_quad <- get_basis(y=quad_points, knots_vec=knots3,hazard="bspline")
    } else{
      quad_points <- transform_quad_points(n_quad = n_quad,
                                           quad_method=quad_method, a = 0,b = y2-y1)
      basis3_quad <- get_basis(y=quad_points, knots_vec=knots3,hazard="bspline")
    }
    lambda03 <- as.vector(exp(basis3_quad %*% phi3))
    Lambda03 <- (y2-y1)/2 * as.vector(matrix(lambda03,ncol=n_quad,byrow = TRUE) %*% quad_weights)
    if(anyLT){
      quad_points <- transform_quad_points(n_quad = n_quad,
                                           quad_method=quad_method, a = 0,b = yL)
      basis1_quad <- get_basis(y=quad_points, knots_vec=knots1,hazard="bspline")
      basis2_quad <- get_basis(y=quad_points, knots_vec=knots2,hazard="bspline")
      lambda01 <- as.vector(exp(basis1_quad %*% phi1))
      lambda02 <- as.vector(exp(basis2_quad %*% phi2))
      #reshape lambda0 from a n_t*n_quad length vector
      #to an n_t by n_quad matrix, then multiply with n_quad length weights
      #to get final Lambda0
      Lambda01_yL <- yL/2 * as.vector(matrix(lambda01,ncol=n_quad,byrow = TRUE) %*% quad_weights)
      Lambda02_yL <- yL/2 * as.vector(matrix(lambda02,ncol=n_quad,byrow = TRUE) %*% quad_weights)
    }
  }


  #compute frailties as means of "posterior" gamma distribution with parameters:
  #shape: delta1 + delta2 + invtheta
  #rate: H1 + H2 + H3 - H1_yL - H2_yL + invtheta
  #mean is then shape / rate

  numer <- delta1 + delta2 + invtheta
  denom <- Lambda01 * exp(eta1) + Lambda02 * exp(eta2) +
    delta1 * Lambda03 * exp(eta3) + invtheta #note extra delta1 to zero out unnecessary contributions
  if(anyLT){
    denom <- denom - Lambda01_yL * exp(eta1) - Lambda02_yL * exp(eta2)
  }
  cbind(mean = numer/denom,
        median = stats::qgamma(p=0.5,shape=numer,rate=denom),
        shape=numer, rate=denom)
}

