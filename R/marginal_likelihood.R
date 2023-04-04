#' Negative Log-Likelihood Function for Illness-Death Model (Numerically Marginalized)
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption.
#'   Typically, this function will not be used directly by the user,
#'   but as part of a larger estimation procedure.
#' @inheritParams nll_func
#' @inheritParams FreqSurv_HReg2
#' @param frailtyparam_num count of parameters associated with frailty
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @importFrom statmod gauss.quad.prob
#' @export
nll_marg_func <- function(para, y1, y2, delta1, delta2, yL, anyLT,
                          Xmat1, Xmat2, Xmat3,
                          hazard, frailtyparam_num, model, weights, n_quad){
  # browser()
  nP1 <- if(!is.null(Xmat1)) ncol(Xmat1) else 0
  nP2 <- if(!is.null(Xmat2)) ncol(Xmat2) else 0
  nP3 <- if(!is.null(Xmat3)) ncol(Xmat3) else 0
  n <- length(y1)

  if(tolower(hazard) %in% c("weibull","wb")){
    #assume logtheta is in the 7th position for now!!
    gauss_list <- statmod::gauss.quad.prob(n = n_quad, dist = "gamma",
                                           alpha=exp(-para[7]),beta=exp(para[7]))

    nP0 <- 6 + frailtyparam_num
    stopifnot(length(para) == nP0 + nP1 + nP2 + nP3)
    nll <- nlogLikWB_ID_marg( para=para,
                              y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                              X1=if(nP1>0) Xmat1 else matrix(nrow = n, ncol=0),
                              X2=if(nP2>0) Xmat2 else matrix(nrow = n, ncol=0),
                              X3=if(nP3>0) Xmat3 else matrix(nrow = n, ncol=0),
                              yL=yL, anyLT=as.numeric(anyLT), weights=weights,
                              model = tolower(model), frailty_ind = as.numeric(frailtyparam_num),
                              gauss_nodes=gauss_list$nodes, gauss_weights=gauss_list$weights)
  } else{
    stop("please choose hazard of 'weibull'")}
  return(nll)
}

#' Conditional Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption, fixing some parameters
#'   Typically, this function will not be used directly by the user,
#'   but as part of a larger estimation procedure.
#' @inheritParams nll_func
#' @param other_para vector of values for remaining parameters.
#' @param fixed_param value of parameter being fixed (can be vector if fixed parameters are contiguous)
#' @param fixed_param_ind index of parameter after which fixed parameters would be inserted
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
nll_marg_profile_helper_func <- function(other_para, fixed_param, fixed_param_ind,
                                    y1, y2, delta1, delta2, yL, anyLT,
                                    Xmat1, Xmat2, Xmat3,
                                    hazard, model, weights,
                                    frailtyparam_num, n_quad){
  #fixed_param can be multidimensional as long as the "fixed" parameters are contiguous
  #i.e., fixed_param_ind is a scalar (this is fine bc use case is like theta+betafrail)
  para <- append(x = other_para,values = fixed_param,after = fixed_param_ind)
  nll <- nll_marg_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                  yL=yL, anyLT=anyLT, weights=weights,
                  frailtyparam_num = frailtyparam_num,
                  hazard=hazard, model=model, n_quad = n_quad)
  return(nll)
}



#' Negative Profile Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the negative profile log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption, fixing some parameters
#'   Typically, this function will not be used directly by the user,
#'   but as part of a larger estimation procedure.
#' @inheritParams nll_profile_helper_func
#' @inheritParams get_fit_uni
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
nll_marg_profile_func <- function(fixed_param, fixed_param_ind, para_mle,
                             y1, y2, delta1, delta2, yL, anyLT,
                             Xmat1, Xmat2, Xmat3, weights,
                             hazard, model, frailtyparam_num, n_quad,
                             verbose=FALSE, control, hessian,
                             optim_method, extra_starts){
  # browser()
  fixed_param <- as.matrix(fixed_param)
  n_fixed <- NCOL(fixed_param)
  n_points <- NROW(fixed_param)
  #if multiple params are being fixed, then fixed_param should be matrix with different values
  start_vec <- other_para <-
    para_mle[-((1+fixed_param_ind):(fixed_param_ind+n_fixed))]

  out_list <- list()
  out_list$paramat <- matrix(nrow=n_points,
                             ncol=length(para_mle),
                             dimnames = list(NULL,names(para_mle)))
  if(hessian) out_list$semat <- out_list$paramat
  out_list$nll <- numeric(n_points)
  out_list$counts <- matrix(nrow=n_points,ncol=2)
  out_list$conv <- numeric(n_points)
  for(i in 1:n_points){
    if(verbose) print(paste0(i," of ",n_points))
    temp_fit <- stats::optim(par = start_vec, fn = nll_marg_profile_helper_func,
                      fixed_param=fixed_param[i,],fixed_param_ind=fixed_param_ind,
                      y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                      hazard=hazard, model=model, weights=weights,
                      frailtyparam_num=frailtyparam_num, n_quad=n_quad,
                      method=optim_method, control=control, hessian=hessian)
    #warm start next iteration at previous
    start_vec <- temp_fit$par
    out_list$paramat[i,] <- append(temp_fit$par, values = fixed_param[i,],
                                   after = fixed_param_ind)
    if(hessian) out_list$semat[i,] <- append(sqrt(diag(MASS::ginv(temp_fit$hessian))),
                          values = numeric(length(fixed_param[i,])),after=fixed_param_ind)
    out_list$nll[i] <- temp_fit$value
    out_list$counts[i,] <- temp_fit$counts
    out_list$conv[i] <- temp_fit$convergence
  }
  return(out_list)
}






#' Conditional Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption, fixing some parameters
#'   Typically, this function will not be used directly by the user,
#'   but as part of a larger estimation procedure.
#' @inheritParams nll_func
#' @param other_para vector of values for remaining parameters.
#' @param fixed_param value of parameter being fixed (can be vector if fixed parameters are contiguous)
#' @param fixed_param_ind index of parameter after which fixed parameters would be inserted
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
nll_profile_helper_func <- function(other_para, fixed_param, fixed_param_ind,
                                    y1, y2, delta1, delta2, yL, anyLT,
                                    Xmat1, Xmat2, Xmat3,
                                    hazard, frailty, model, weights,
                                    basis1, basis2, basis3, basis3_y1,
                                    basis1_yL,basis2_yL,
                                    dbasis1, dbasis2, dbasis3){
  #fixed_param can be multidimensional as long as the "fixed" parameters are contiguous
  #i.e., fixed_param_ind is a scalar (this is fine bc use case is like theta+betafrail)
  para <- append(x = other_para,values = fixed_param,after = fixed_param_ind)
  nll <- nll_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                  hazard=hazard, frailty=frailty, model=model, weights=weights,
                  basis1=basis1, basis2=basis2, basis3=basis3,
                  basis3_y1=basis3_y1, basis1_yL=basis1_yL,basis2_yL=basis2_yL,
                  dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)
  return(nll)
}


#' Conditional Negative Log-Likelihood Gradient Function for Illness-Death Model
#'
#' Function returning the gradient of the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption, fixing some parameters.
#'   Typically, this function will not be used directly by the user,
#'   but as part of a larger estimation procedure.
#' @inheritParams nll_func
#' @param other_para vector of values for remaining parameters.
#' @param fixed_param value of parameter being fixed (can be vector if fixed parameters are contiguous)
#' @param fixed_param_ind index of parameter after which fixed parameters would be inserted
#'
#' @return Returns numeric sum of negative log likelihood gradient contributions.
#' @export
nll_profile_helper_grad_func <- function(other_para, fixed_param, fixed_param_ind,
                                    y1, y2, delta1, delta2, yL, anyLT,
                                    Xmat1, Xmat2, Xmat3,
                                    hazard, frailty, model, weights,
                                    basis1, basis2, basis3, basis3_y1,
                                    basis1_yL,basis2_yL,
                                    dbasis1, dbasis2, dbasis3){
  #fixed_param can be multidimensional as long as the "fixed" parameters are contiguous
  #i.e., fixed_param_ind is a scalar (this is fine bc use case is like theta+betafrail)
  n_fixed <- length(fixed_param)
  para <- append(x = other_para,values = fixed_param,after = fixed_param_ind)
  #note that then gradient for "added" component must be again removed.
  ngrad <- ngrad_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                  Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                  hazard=hazard, frailty=frailty, model=model, weights=weights,
                  basis1=basis1, basis2=basis2, basis3=basis3,
                  basis3_y1=basis3_y1, basis1_yL=basis1_yL,basis2_yL=basis2_yL,
                  dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3)[-((1+fixed_param_ind):(fixed_param_ind+n_fixed))]
  return(ngrad)
}





#' Conditional Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption, fixing some parameters
#'   Typically, this function will not be used directly by the user,
#'   but as part of a larger estimation procedure.
#' @inheritParams nll_profile_helper_func
#' @inheritParams get_fit_uni
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
nll_profile_func <- function(fixed_param, fixed_param_ind, para_mle,
                             y1, y2, delta1, delta2, yL, anyLT,
                             Xmat1, Xmat2, Xmat3,
                             hazard, frailty, model, weights,
                             basis1, basis2, basis3, basis3_y1, basis1_yL,basis2_yL,
                             dbasis1, dbasis2, dbasis3,verbose=FALSE,
                             control,optim_method,extra_starts){
  # browser()
  fixed_param <- as.matrix(fixed_param)
  n_fixed <- NCOL(fixed_param)
  n_points <- NROW(fixed_param)
  #if multiple params are being fixed, then fixed_param should be matrix with different values
  start_vec <- other_para <-
    para_mle[-((1+fixed_param_ind):(fixed_param_ind+n_fixed))]

  out_list <- list()
  out_list$paramat <- matrix(nrow=n_points,
                             ncol=length(para_mle),
                             dimnames = list(NULL,names(para_mle)))
  out_list$nll <- numeric(n_points)
  out_list$counts <- matrix(nrow=n_points,ncol=2)
  out_list$conv <- numeric(n_points)
  for(i in 1:n_points){
    if(verbose) print(paste0(i," of ",n_points))
    temp_fit <- stats::optim(par = start_vec, fn = nll_profile_helper_func,
                      gr = nll_profile_helper_grad_func,
                      fixed_param=fixed_param[i,],fixed_param_ind=fixed_param_ind,
                      y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
                      hazard=hazard, frailty=frailty, model=model, weights=weights,
                      basis1=basis1, basis2=basis2, basis3=basis3,
                      basis3_y1=basis3_y1, basis1_yL=basis1_yL,basis2_yL=basis2_yL,
                      dbasis1=dbasis1, dbasis2=dbasis2, dbasis3=dbasis3,
                      method=optim_method, control=control, hessian=FALSE)
    #warm start next iteration at previous
    start_vec <- temp_fit$par
    out_list$paramat[i,] <- append(temp_fit$par, values = fixed_param[i,],
                                   after = fixed_param_ind)
    out_list$nll[i] <- temp_fit$value
    out_list$counts[i,] <- temp_fit$counts
    out_list$conv[i] <- temp_fit$convergence
  }
  return(out_list)
}













#' Conditional Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption, fixing some parameters
#'   Typically, this function will not be used directly by the user,
#'   but as part of a larger estimation procedure.
#' @inheritParams nll_func
#' @param other_para vector of values for remaining parameters.
#' @param fixed_param value of parameter being fixed (can be vector if fixed parameters are contiguous)
#' @param fixed_param_ind index of parameter after which fixed parameters would be inserted
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
nll_profile_theta_helper_func <- function(other_para, fixed_param, fixed_param_ind,
                                          y1, y2, delta1, delta2, yL, anyLT,
                                          Xmat1, Xmat2, Xmat3, model, weights){
  #fixed_param can be multidimensional as long as the "fixed" parameters are contiguous
  #i.e., fixed_param_ind is a scalar (this is fine bc use case is like theta+betafrail)
  para <- append(x = other_para,values = fixed_param,after = fixed_param_ind)
  nll <- nlogLikWB_ID_theta(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                            X1=Xmat1, X2=Xmat2, X3=Xmat3, model=model, frailty_ind=1, weights=weights)
  return(nll)
}

#' Conditional Negative Log-Likelihood Function for Illness-Death Model
#'
#' Function returning the negative log-likelihood for the illness-death model,
#'   under specified baseline hazard, and specified frailty,
#'   and specified Markov/semi-Markov transition assumption, fixing some parameters
#'   Typically, this function will not be used directly by the user,
#'   but as part of a larger estimation procedure.
#' @inheritParams nll_profile_helper_func
#' @inheritParams get_fit_uni
#'
#' @return Returns numeric sum of negative log likelihood contributions.
#' @export
nll_profile_theta_func <- function(fixed_param, fixed_param_ind, para_mle,
                                   y1, y2, delta1, delta2, yL, anyLT,
                                   Xmat1, Xmat2, Xmat3,
                                   hazard, frailty, model, weights,
                                   basis1, basis2, basis3, basis3_y1, basis1_yL,basis2_yL,
                                   dbasis1, dbasis2, dbasis3,verbose=FALSE,
                                   control,hessian,optim_method,extra_starts){
  # browser()
  fixed_param <- as.matrix(fixed_param)
  n_fixed <- NCOL(fixed_param)
  n_points <- NROW(fixed_param)
  #if multiple params are being fixed, then fixed_param should be matrix with different values
  start_vec <- other_para <-
    para_mle[-((1+fixed_param_ind):(fixed_param_ind+n_fixed))]

  out_list <- list()
  out_list$paramat <- matrix(nrow=n_points,
                             ncol=length(para_mle),
                             dimnames = list(NULL,names(para_mle)))
  out_list$nll <- numeric(n_points)
  out_list$counts <- matrix(nrow=n_points,ncol=2)
  out_list$conv <- numeric(n_points)
  for(i in 1:n_points){
    if(verbose) print(paste0(i," of ",n_points))
    temp_fit <- stats::optim(par = start_vec, fn = nll_profile_theta_helper_func,
                      fixed_param=fixed_param[i,],fixed_param_ind=fixed_param_ind,
                      y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
                      Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3, model=model,
                      weights=weights,
                      method=optim_method, control=control, hessian=FALSE)
    #warm start next iteration at previous
    start_vec <- temp_fit$par
    out_list$paramat[i,] <- append(temp_fit$par, values = fixed_param[i,],
                                   after = fixed_param_ind)
    out_list$nll[i] <- temp_fit$value
    out_list$counts[i,] <- temp_fit$counts
    out_list$conv[i] <- temp_fit$convergence
  }
  return(out_list)
}












#eventually replace the above functions with this "overarching" function, but
#for now set it aside so I don't break the existing code I have



#' #' Conditional Negative Log-Likelihood Function for Illness-Death Model
#' #'
#' #' Function returning the negative log-likelihood for the illness-death model,
#' #'   under specified baseline hazard, and specified frailty,
#' #'   and specified Markov/semi-Markov transition assumption, fixing some parameters
#' #'   Typically, this function will not be used directly by the user,
#' #'   but as part of a larger estimation procedure.
#' #' @inheritParams nll_func
#' #' @param other_para vector of values for remaining parameters.
#' #' @param fixed_param value of parameter being fixed (can be vector if fixed parameters are contiguous)
#' #' @param fixed_param_ind index of parameter after which fixed parameters would be inserted
#' #'
#' #' @return Returns numeric sum of negative log likelihood contributions.
#' #' @export
#' nll_profile_helper_func <- function(other_para, fixed_param, fixed_param_ind,
#'                                     y1, y2, delta1, delta2, yL, anyLT,
#'                                     Xmat1, Xmat2, Xmat3,
#'                                     hazard, frailtyparam_num, model, weights,
#'                                     n_quad, nll_type){
#'   stopifnot(tolower(hazard) %in% c("weibull","wb")) #restrict to weibull for now
#'
#'   #fixed_param can be multidimensional as long as the "fixed" parameters are contiguous
#'   #i.e., fixed_param_ind is a scalar (this is fine bc use case is like theta+betafrail)
#'   para <- append(x = other_para,values = fixed_param,after = fixed_param_ind)
#'
#'   if(tolower(nll_type) %in% c("m","marg","marginal")){
#'     nll <- nll_marg_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2,
#'                          Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
#'                          yL=yL, anyLT=anyLT, weights=weights,
#'                          frailty_ind = frailtyparam_num,
#'                          hazard=hazard, model=model, n_quad = n_quad)
#'   } else{
#'     stopifnot(frailtyparam_num %in% c(0,1))
#'     nll <- nll_func(para=para, y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
#'                     Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
#'                     hazard=hazard, frailty=frailty, model=model, weights=weights,
#'                     basis1=NULL, basis2=NULL, basis3=NULL,
#'                     basis3_y1=NULL, basis1_yL=NULL,basis2_yL=NULL,
#'                     dbasis1=NULL, dbasis2=NULL, dbasis3=NULL)
#'   }
#'   return(nll)
#' }
#'
#'
#' #' Negative Profile Log-Likelihood Function for Illness-Death Model
#' #'
#' #' Function returning the negative profile log-likelihood for the illness-death model,
#' #'   under specified baseline hazard, and specified frailty,
#' #'   and specified Markov/semi-Markov transition assumption, fixing some parameters
#' #'   Typically, this function will not be used directly by the user,
#' #'   but as part of a larger estimation procedure.
#' #' @inheritParams nll_profile_helper_func
#' #' @inheritParams get_fit_uni
#' #'
#' #' @return Returns numeric sum of negative log likelihood contributions.
#' #' @export
#' nll_profile_func <- function(fixed_param, fixed_param_ind, para_mle,
#'                              y1, y2, delta1, delta2, yL, anyLT,
#'                              Xmat1, Xmat2, Xmat3,
#'                              hazard, model, frailtyparam_num, n_quad, nll_type,
#'                              verbose=FALSE, control, hessian,
#'                              optim_method, extra_starts){
#'   # browser()
#'   fixed_param <- as.matrix(fixed_param)
#'   n_fixed <- NCOL(fixed_param)
#'   n_points <- NROW(fixed_param)
#'   #if multiple params are being fixed, then fixed_param should be matrix with different values
#'   start_vec <- other_para <-
#'     para_mle[-((1+fixed_param_ind):(fixed_param_ind+n_fixed))]
#'
#'   out_list <- list()
#'   out_list$paramat <- matrix(nrow=n_points,
#'                              ncol=length(para_mle),
#'                              dimnames = list(NULL,names(para_mle)))
#'   if(hessian) out_list$semat <- out_list$paramat
#'   out_list$nll <- numeric(n_points)
#'   out_list$counts <- matrix(nrow=n_points,ncol=2)
#'   out_list$conv <- numeric(n_points)
#'   for(i in 1:n_points){
#'     if(verbose) print(paste0(i," of ",n_points))
#'     temp_fit <- optim(par = start_vec, fn = nll_profile_helper_func,
#'                       fixed_param=fixed_param[i,],fixed_param_ind=fixed_param_ind,
#'                       y1=y1, y2=y2, delta1=delta1, delta2=delta2, yL=yL, anyLT=anyLT,
#'                       Xmat1=Xmat1, Xmat2=Xmat2, Xmat3=Xmat3,
#'                       hazard=hazard, model=model, weights=weights,
#'                       frailtyparam_num=frailtyparam_num, n_quad=n_quad,
#'                       method=optim_method, control=control, hessian=hessian)
#'     #warm start next iteration at previous
#'     start_vec <- temp_fit$par
#'     out_list$paramat[i,] <- append(temp_fit$par, values = fixed_param[i,],
#'                                    after = fixed_param_ind)
#'     if(hessian) out_list$semat[i,] <- append(sqrt(diag(MASS::ginv(temp_fit$hessian))),
#'                                              values = numeric(length(fixed_param[i,])),after=fixed_param_ind)
#'     out_list$nll[i] <- temp_fit$value
#'     out_list$counts[i,] <- temp_fit$counts
#'     out_list$conv[i] <- temp_fit$convergence
#'   }
#'   return(out_list)
#' }

