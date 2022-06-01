#' Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#' @param para A numeric vector of parameters, arranged as follows:
#'   the first \eqn{k} elements correspond to the baseline hazard parameters,
#'   then the last\eqn{q} elements correspond with the regression parameters.
#' @param y Numeric vector of length \eqn{n} with (possibly censored) non-terminal and terminal event times
#' @param delta Numeric vector of length \eqn{n}  with indicator of 1 if the event was observed and 0 otherwise
#' @param yL Numeric vector of length \eqn{n} with left-truncation times.
#' @param anyLT Boolean (or 1/0) value that is \code{TRUE} (or 1) if there are any non-zero left truncation times, and \code{FALSE} (or 0) otherwise
#' @param Xmat Numeric matrices with \eqn{n} rows and \eqn{q_1,q_2,q_3} columns containing covariates.
#' @param hazard String (not case sensitive) specifying the form of the baseline hazard.
#'   Options are \code{"weibull"} for Weibull,
#'   \code{"piecewise"} for piecewise constant on the hazard scale,
#'   \code{"bspline"} for cubic B-spline on the log-hazard scale,
#'   and \code{"royston-parmar"} for restricted cubic spline on the log-cumulative hazard scale.
#'   Aliases for these are \code{"wb"}, \code{"pw"}, \code{"bs"}, and \code{"rp"} respectively.
#' @param basis Numeric matrix with \eqn{n} rows and \eqn{k} columns
#'   with piecewise/spline basis function values at the corresponding \code{y} values.
#'   Not used under Weibull model.
#' @param dbasis Numeric matrix with \eqn{n} rows and \eqn{k} columns
#'   with piecewise/spline basis function derivative values at the corresponding \code{y} values.
#'   Used only under Royston-Parmar model.
#' @param basis_yL Numeric matrix with \eqn{n} rows and \eqn{k} columns
#'   with piecewise/spline basis function values at the corresponding \code{yL} values.
#'   Not used under Weibull model.
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
nll_uni_func <- function(para, y, delta, yL, anyLT, Xmat, hazard, basis, dbasis, basis_yL){
  # browser()
  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP <- if(!is.null(Xmat)) ncol(Xmat) else 0
  n <- length(y)

  if(tolower(hazard) %in% c("weibull","wb")){
    nll <- nlogLikWB_uni(para, y=y,delta=delta,
                         yL=if(anyLT) yL else numeric(0),
                         anyLT = as.numeric(anyLT),
                         X=if(nP>0) as.matrix(Xmat) else matrix(nrow = n, ncol=0))
  } else if(tolower(hazard) %in%  c("bspline","bs")){
    # assume that basis has entries 1:n of original data, then
    # n blocks of n_quad rows corresponding the quadrature points for each subject
    n_quad <- NROW(basis)/NROW(y) - 1
    quad_method <- attr(basis,which="quad_method")
    nll <- nlogLikBS_uni(para, y=y,delta=delta,
                         X=if(nP>0) as.matrix(Xmat) else matrix(nrow = n, ncol=0),
                         basis=basis,quad_weights = get_quad_pointsweights(n_quad = n_quad,quad_method = quad_method)$weights)
  } else if(tolower(hazard) %in% c("piecewise","pw")){
    nll <- nlogLikPW_uni(para, y=y,delta=delta,
                         X=if(nP>0) as.matrix(Xmat) else matrix(nrow = n, ncol=0),
                         basis=basis,dbasis=dbasis)
  } else if(tolower(hazard) %in% c("royston-parmar","rp")){
    nll <- nlogLikRP_uni(para, y=y,delta=delta,
                         X=if(nP>0) as.matrix(Xmat) else matrix(nrow = n, ncol=0),
                         basis=basis,dbasis=dbasis,
                         basis_yL=if(anyLT) basis_yL else matrix(nrow = n, ncol=0),
                         anyLT = as.numeric(anyLT))
  } else{ stop("please choose hazard of 'weibull', 'bspline', 'royston-parmar', or 'piecewise'")}
  return(nll)
}

#' Gradient of Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#' @inheritParams nll_uni_func
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
ngrad_uni_func <- function(para, y, delta, yL, anyLT, Xmat, hazard, basis, dbasis, basis_yL){
  # browser()
  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP <- if(!is.null(Xmat)) ncol(Xmat) else 0
  n <- length(y)

  if(tolower(hazard) %in% c("weibull","wb")){
    ngrad <- ngradWB_uni(para, y=y,delta=delta,
                         yL=if(anyLT) yL else numeric(0),
                         anyLT = as.numeric(anyLT),
                         X=if(nP>0) as.matrix(Xmat) else matrix(nrow = n, ncol=0))
  } else if(tolower(hazard) %in%  c("bspline","bs")){
    # assume that basis has entries 1:n of original data, then
    # n blocks of n_quad rows corresponding the quadrature points for each subject
    n_quad <- NROW(basis)/NROW(y) - 1
    quad_method <- attr(basis,which="quad_method")
    ngrad <- ngradBS_uni(para, y=y,delta=delta, X=if(nP>0) as.matrix(Xmat) else matrix(nrow = n, ncol=0),
                         basis=basis,quad_weights = get_quad_pointsweights(n_quad = n_quad,quad_method = quad_method)$weights)
  } else if(tolower(hazard) %in% c("piecewise","pw")){
    ngrad <- ngradPW_uni(para, y=y,delta=delta,
                         X=if(nP>0) as.matrix(Xmat) else matrix(nrow = n, ncol=0),
                         basis=basis,dbasis=dbasis)
  } else if(tolower(hazard) %in% c("royston-parmar","rp")){
    ngrad <- ngradRP_uni(para, y=y,delta=delta, X=if(nP>0) as.matrix(Xmat) else matrix(nrow = n, ncol=0),
                         basis=basis,dbasis=dbasis,
                         basis_yL=if(anyLT) basis_yL else matrix(nrow = n, ncol=0),
                         anyLT = as.numeric(anyLT))
  } else{ stop("please choose hazard of 'weibull', 'bspline', 'royston-parmar', or 'piecewise'")}
  return(ngrad)
}







#' Combined function with Negative Log-Likelihood Function and Gradient for Univariate Survival Model
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#' @inheritParams nll_uni_func
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
nll_ngrad_uni_func <- function(para, y, delta, yL, anyLT, Xmat, hazard,
                               basis, dbasis, basis_yL){
  value <- nll_uni_func(para=para, y=y, delta=delta, yL=yL,
              anyLT=anyLT, Xmat=Xmat, hazard=hazard,
              basis=basis, dbasis=dbasis, basis_yL=basis_yL)
  attr(value,"gradient") <- ngrad_uni_func(para=para, y=y, delta=delta, yL=yL,
                               anyLT=anyLT, Xmat=Xmat, hazard=hazard,
                               basis=basis, dbasis=dbasis, basis_yL=basis_yL)
  value
}



#' Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#' @param para A numeric vector of parameters, arranged as follows:
#'   the first \eqn{k_1+k_2+k_3} elements correspond to the baseline hazard parameters,
#'   then the \eqn{k_1+k_2+k_3+1} element corresponds to the gamma frailty log-variance parameter,
#'   then the last\eqn{q_1+q_2+q_3} elements correspond with the regression parameters.
#' @param y1,y2 Numeric vectors of length \eqn{n} with (possibly censored) non-terminal and terminal event times
#' @param delta1,delta2 Numeric vectors of length \eqn{n}  with indicators of 1 if the event was observed and 0 otherwise
#' @param Xmat1,Xmat2,Xmat3 Numeric matrices with \eqn{n} rows and \eqn{q_1,q_2,q_3} columns containing covariates.
#' @param hazard String (not case sensitive) specifying the form of the baseline hazard.
#'   Options are \code{"weibull"} for Weibull,
#'   \code{"piecewise"} for piecewise constant on the hazard scale,
#'   \code{"bspline"} for cubic B-spline on the log-hazard scale,
#'   and \code{"royston-parmar"} for restricted cubic spline on the log-cumulative hazard scale.
#'   Aliases for these are \code{"wb"}, \code{"pw"}, \code{"bs"}, and \code{"rp"} respectively.
#' @param frailty Boolean indicating whether a gamma distributed subject-specific frailty should
#'   be included. Currently this must be set to \code{TRUE}.
#' @param model String (not case sensitive) specifying the transition assumption: either \code{"semi-Markov"}
#' @param basis1,basis2,basis3,basis3_y1 Numeric matrices with \eqn{n} rows and \eqn{k_1,k_2,k_3} columns
#'   with piecewise/spline basis function values at the corresponding \code{y1} and \code{y2} values.
#'   Under semi-Markov model, basis3 represents basis derived from \eqn{y_2-y_1} and \code{basis3_y1} is unused,
#'   while under Markov model, basis3 represents basis derived from \eqn{y_2} and \code{basis3_y1} is from \eqn{y_1}
#'   Not used under Weibull model.
#' @param dbasis1,dbasis2,dbasis3 Numeric matrices with \eqn{n} rows and \eqn{k_1,k_2,k_3} columns
#'   with piecewise/spline basis function derivative values at the corresponding \code{y1} and \code{y2} values.
#'   Used only under Royston-Parmar model.
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
nll_func <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                     hazard, frailty, model,
                     basis1, basis2, basis3, basis3_y1,
                     dbasis1, dbasis2, dbasis3){

  # browser()

  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  nP01 <- if(!is.null(basis1)) ncol(basis1) else 0
  nP02 <- if(!is.null(basis2)) ncol(basis2) else 0
  nP03 <- if(!is.null(basis3)) ncol(basis3) else 0
  n <- length(y1)

  if(tolower(hazard) %in% c("weibull","wb")){
    # if(frailty){
      nP0 <- 6 + as.numeric(frailty)
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
      nll <- nlogLikWB_ID(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                          X1=if(nP1>0) Xmat1 else matrix(nrow = n, ncol=0),
                          X2=if(nP2>0) Xmat2 else matrix(nrow = n, ncol=0),
                          X3=if(nP3>0) Xmat3 else matrix(nrow = n, ncol=0),
                          model = tolower(model), frailty_ind = as.numeric(frailty))
    # } else{stop("non-frailty weibull model not yet implemented")}
  }
  else if(tolower(hazard) %in%  c("bspline","bs")){
    # assume that basis has entries 1:n of original data, then
    # n blocks of n_quad rows corresponding the quadrature points for each subject
    n_quad <- NROW(basis1)/NROW(y1) - 1
    quad_method <- attr(basis1,which="quad_method")
    # if(frailty){
      nP0 <- nP01 + nP02 + nP03 + as.numeric(frailty)
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
      nll <- nlogLikBS_ID(para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                X1=if(nP1>0) Xmat1 else matrix(nrow = n, ncol=0),
                X2=if(nP2>0) Xmat2 else matrix(nrow = n, ncol=0),
                X3=if(nP3>0) Xmat3 else matrix(nrow = n, ncol=0),
                basis1 = basis1, basis2 = basis2, basis3=basis3,
                quad_weights = get_quad_pointsweights(n_quad = n_quad,quad_method = quad_method)$weights,
                frailty_ind = as.numeric(frailty))
    # } else{
    #   stop("non-frailty bspline model not yet implemented")
    # }
  }
  else if(tolower(hazard) %in% c("piecewise","pw")){
    # if(frailty){
      nP0 <- nP01 + nP02 + nP03 + as.numeric(frailty)
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
      nll <- nlogLikPW_ID(para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                          X1=if(nP1>0) Xmat1 else matrix(nrow = n, ncol=0),
                          X2=if(nP2>0) Xmat2 else matrix(nrow = n, ncol=0),
                          X3=if(nP3>0) Xmat3 else matrix(nrow = n, ncol=0),
                          basis1=basis1, basis2=basis2, basis3=basis3,
                          dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                          frailty_ind=as.numeric(frailty))
    # } else{
    #   stop("non-frailty piecewise constant model not yet implemented")
    # }
  } else if(tolower(hazard) %in% c("royston-parmar","rp")){
    # if(frailty){
      nP0 <- nP01 + nP02 + nP03 + as.numeric(frailty)
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss

      nll <- nlogLikRP_ID(para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                 X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                 X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                 X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0),
                 basis1=basis1, basis2=basis2, basis3=basis3,
                 basis3_y1=if(!is.null(basis3_y1)) basis3_y1 else matrix(nrow = n, ncol=0),
                 dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                 model = tolower(model), frailty_ind=as.numeric(frailty))

      # if(tolower(model)=="semi-markov"){
      #   nll <- nlogLikRP_ID_frail_SM(para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
      #                                X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
      #                                X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
      #                                X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0),
      #                                basis1=basis1, basis2=basis2, basis3=basis3,
      #                                dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
      #                                frailty_ind=as.numeric(frailty))
      # } else{
      #   nll <- nlogLikRP_ID_frail_M(para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
      #                               X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
      #                               X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
      #                               X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0),
      #                               basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1 = basis3_y1,
      #                               dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
      #                               frailty_ind=as.numeric(frailty))
      # }
    # } else {"non-frailty royston-parmar model not yet implemented"}
  } else{ stop("please choose hazard of 'weibull', 'bspline', 'royston-parmar', or 'piecewise'")}

  return(nll)
}



#' Gradient of Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the gradient of the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#'
#' @inheritParams nll_func
#'
#' @return Returns numeric vector of same length and arrangement as \code{para}
#'   with sum of gradient contributions for the negative log likelihood.
#' @export
ngrad_func <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                       hazard, frailty, model,
                       basis1, basis2, basis3, basis3_y1,
                       dbasis1, dbasis2, dbasis3){

  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  nP01 <- if(!is.null(basis1)) ncol(basis1) else 0
  nP02 <- if(!is.null(basis2)) ncol(basis2) else 0
  nP03 <- if(!is.null(basis3)) ncol(basis3) else 0
  n <- length(y1)

  if(tolower(hazard) %in% c("weibull","wb")){
    # if(frailty){
      nP0 <- 6 + as.numeric(frailty)
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss

      ngrad <- ngradWB_ID(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                   X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                   X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                   X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0),
                   model = tolower(model),frailty_ind=as.numeric(frailty))

      # if(tolower(model)=="semi-markov"){
      #   ngrad <- ngradWB_ID_frail_SM(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
      #                                X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
      #                                X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
      #                                X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0))
      # } else{ #markov model
      #   ngrad <- ngradWB_ID_frail_M(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
      #                               X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
      #                               X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
      #                               X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0))
      # }
    # } else{stop("non-frailty not yet implemented")}
  }
  else if(tolower(hazard) %in% c("bspline","bs")){
    n_quad <- NROW(basis1)/NROW(y1) - 1
    quad_method <- attr(basis1,which="quad_method")
    # if(frailty){
    nP0 <- nP01 + nP02 + nP03 + as.numeric(frailty)
    stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
    ngrad <- ngradBS_ID(para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
              X1=if(nP1>0) Xmat1 else matrix(nrow = n, ncol=0),
              X2=if(nP2>0) Xmat2 else matrix(nrow = n, ncol=0),
              X3=if(nP3>0) Xmat3 else matrix(nrow = n, ncol=0),
              basis1 = basis1, basis2 = basis2, basis3=basis3,
              quad_weights = get_quad_pointsweights(n_quad = n_quad,quad_method = quad_method)$weights,
              frailty_ind=as.numeric(frailty))
    # } else{
    #   stop("non-frailty bspline model not yet implemented.")
    # }
  }
  else if(tolower(hazard) %in% c("piecewise","pw")){
    # if(frailty){
      nP0 <- nP01 + nP02 + nP03 + as.numeric(frailty)
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
        ngrad <- ngradPW_ID(para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                            X1=if(nP1>0) Xmat1 else matrix(nrow = n, ncol=0),
                            X2=if(nP2>0) Xmat2 else matrix(nrow = n, ncol=0),
                            X3=if(nP3>0) Xmat3 else matrix(nrow = n, ncol=0),
                            basis1=basis1, basis2=basis2, basis3=basis3,
                            dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                            frailty_ind=as.numeric(frailty))
    # } else{
    #   stop("non-frailty piecewise model not yet implemented.")
    # }
  } else if(tolower(hazard) %in% c("royston-parmar","rp")){
    # if(frailty){
      nP0 <- nP01 + nP02 + nP03 + as.numeric(frailty)
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss

      ngrad <- ngradRP_ID(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                  X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                  X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                  X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0),
                  basis1=basis1, basis2=basis2, basis3=basis3,
                  basis3_y1=if(!is.null(basis3_y1)) basis3_y1 else matrix(nrow = n, ncol=0),
                  dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                  model = tolower(model),frailty_ind = as.numeric(frailty))

    #   if(tolower(model)=="semi-markov"){
    #     ngrad <- ngradRP_ID_frail_SM(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
    #                                  X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
    #                                  X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
    #                                  X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0),
    #                                  basis1=basis1, basis2=basis2, basis3=basis3,
    #                                  dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
    #   } else{
    #     ngrad <- ngradRP_ID_frail_M(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
    #                                 X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
    #                                 X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
    #                                 X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0),
    #                                 basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
    #                                 dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
    #   }
    # } else {"non-frailty royston-parmar model not yet implemented"}
  } else {stop("please choose hazard of 'weibull', 'bspline', 'royston-parmar', or 'piecewise'")}
  return( ngrad )
}


#' Combined function with Negative Log-Likelihood Function and Gradient for Univariate Survival Model
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#' @inheritParams nll_func
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
nll_ngrad_func <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                           hazard, frailty, model,
                           basis1, basis2, basis3, basis3_y1,
                           dbasis1, dbasis2, dbasis3){
  value <- nll_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                    Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                    hazard=hazard, frailty=frailty, model=model,
                    basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                    dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
  attr(value,"gradient") <- ngrad_func(para=para, y1=y1, y2=y2,
                             delta1=delta1, delta2=delta2,
                             Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                             hazard=hazard, frailty=frailty, model=model,
                             basis1=basis1, basis2=basis2, basis3=basis3, basis3_y1=basis3_y1,
                             dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
  value
}







#' Gradient Matrix of Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning a matrix of contributions to the gradient of the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#'
#' @inheritParams nll_func
#'
#' @return Returns numeric matrix with \eqn{n} rows and width same as length of \code{para} with gradient contributions
#'   for the negative log likelihood.
#' @export
ngrad_mat_func <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                           hazard, frailty, model){

  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0

  if(tolower(hazard) %in% c("weibull","wb")){
    if(frailty){
      nP0 <- 7
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
      if(tolower(model)=="semi-markov"){
        ngrad_mat <- ngradWB_ID_frail_mat_SM(para = para,y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                             X1=as.matrix(Xmat1), X2=as.matrix(Xmat2), X3=as.matrix(Xmat3))
      } else{
        ngrad_mat <- ngradWB_ID_frail_mat_M(para = para,
                                            y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                            X1=as.matrix(Xmat1), X2=as.matrix(Xmat2), X3=as.matrix(Xmat3))
      }
    } else{stop("non-frailty not yet implemented")}
  } else{stop("non-Weibull not yet implemented")}
  return(ngrad_mat)

}




#' Hessian of Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the Hessian of the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption.
#'   Typically, this function will not be used directly by the user, but as part of a
#'   larger estimation procedure.
#'
#' @inheritParams nll_func
#'
#' @return Returns numeric square matrix with dimensions the same length as \code{para}
#'   with sum of gradient contributions for the negative log likelihood.
#' @export
nhess_func <- function(para, y1, y2, delta1, delta2, Xmat1, Xmat2, Xmat3,
                       frailty, hazard, model){

  #number of parameters in each arm dictated by number of covariate columns in each matrix
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)

  if(tolower(hazard) %in% c("weibull","wb")){
    if(frailty){
      nP0 <- 7
      stopifnot(length(para) == nP0 + nP1 + nP2 + nP3) #if the size of the parameter vector doesn't match the expected size, throw a fuss
      if(tolower(model)=="semi-markov"){
        nhess <- nhessWB_ID_frail_SM(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                     X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                                     X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                                     X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0))
      } else{ #markov model
        nhess <- nhessWB_ID_frail_M(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                                    X1=if(nP1>0) as.matrix(Xmat1) else matrix(nrow = n, ncol=0),
                                    X2=if(nP2>0) as.matrix(Xmat2) else matrix(nrow = n, ncol=0),
                                    X3=if(nP3>0) as.matrix(Xmat3) else matrix(nrow = n, ncol=0))
      }
    } else{stop("non-frailty not yet implemented")}
  } else{stop("non-Weibull not yet implemented")}

  return(nhess)

}
