
# #idea adapted from https://github.com/lgaborini/bayessource/blob/master/src/statistical_functions.cpp
# logcumsumexp <- function(x){
#   xmax <- max(x)
#   return(xmax + log(cumsum(exp(x-xmax))))
# }
# logsumexp <- function(x,y){
#   xmax <- pmax(x,y)
#   return(xmax + log(exp(x-xmax) + exp(y-xmax)))
# }
# logdiffexp <- function(x,y){
#   xmax <- pmax(x,y)
#   return(xmax + log(exp(x-xmax) - exp(y-xmax)))
# }


#need derivative with respect to theta!!!
#double check derivative wrt logtheta or theta
#https://www.wolframalpha.com/input?i=derivative+of+%281+%2B+exp%28x%29+*+a%29%5E%28-exp%28-x%29%29
#https://www.wolframalpha.com/input?i=derivative+of+%281+%2B+x+*+a%29%5E%28-1%2Fx%29
marg_logtheta_deriv <- function(a,ltheta){
  (a * exp(ltheta) + 1)^(-exp(-ltheta)) * (
    exp(-ltheta) * log1p(a * exp(ltheta)) -
      a / (a * exp(ltheta) + 1)
  )
}


#in this script I'm going to build up the following tools
#1. a function to predict (log) hazard or (log) cumulative hazard
#2. a function to compute jacobian matrix for (log) hazard or (log) cumulative hazard
#3. a function to output complete prediction based on (1) and (2) of
#h1 h2 h3, H1 H2 H3, S1 S2 S3
#4. a function for illness-death models to (optionally) compute allllll the stuff
#CIFs and risk profiles, cond and marg
#all with delta method standard errors
#by sharing as much redundant intermediate product as possible
#with a function argument saying which things you do and don't want output

#function strictly to predict (log or natural) hazard and cumulative hazard
#using outer() function, also admits vectors of tseq and eta
#in order to generate matrix output of dimension length(tseq) by length(eta)

#' @export
get_haz <- function(tseq, eta,
                    hazard, phi, basis, dbasis, basis_quad, quad_weights,
                    func_type, log_out){
  # browser()
  stopifnot(func_type %in% c("h","H"))
  # stopifnot(tseq[1]>0)
  eta <- as.vector(eta)

  if(hazard %in% c("weibull","wb")){
    if(func_type=="h"){
      #Weibull loghazard is logalpha + logkappa + (exp(logalpha)-1) * log(t) + eta
      out_temp <- outer(as.vector(phi[2] + phi[1] + expm1(phi[2]) * log(tseq)),
                        eta, "+")
    } else {
      #out_temp is log-Lambda, may be transformed to survival at the end
      out_temp <- outer(as.vector(phi[1] + exp(phi[2]) * log(tseq)),
                        eta, "+")
    }
    #for weibull, it is natural to compute on log scale, so exponentiate as needed
    if(!log_out) out_temp <- exp(out_temp)

  } else if(hazard %in% c("piecewise","pw")){
    if(func_type=="h"){
      out_temp <- outer(as.vector(dbasis %*% phi), eta, "+")
      if(!log_out) out_temp <- exp(out_temp)
    } else {
      out_temp <- outer(as.vector(basis %*% exp(phi)),
                        exp(eta), "*")
      if(log_out) out_temp <- log(out_temp)
    }

  } else if(hazard %in% c("royston-parmar","rp")){
    #this is log-cumulative hazard
    out_temp <- outer(as.vector(basis %*% phi), eta, "+")

    #if hazard is requested, update it
    if(func_type=="h"){
      #implicitly, this is columnwise addition, which is what we want
      out_temp <- as.vector(log(dbasis %*% phi) - log(tseq)) + out_temp
    }

    if(!log_out) out_temp <- exp(out_temp)

  } else if(hazard %in% c("bspline","bs")){
    if(func_type=="h"){
      out_temp <- outer(as.vector(basis %*% phi),
                        eta, "+")
      if(!log_out) out_temp <- exp(out_temp)

    } else {
      n_quad <- dim(basis_quad)[1] / dim(basis)[1]
      # lambda0 <- as.vector(exp(basis_quad %*% phi))
      #reshape lambda0 from a n*n_quad length vector
      #to an n by n_quad matrix, then multiply with n_quad length weights
      #to get final Lambda0
      Lambda0 <- tseq/2 * as.vector(
        matrix(as.vector(exp(basis_quad %*% phi)),
               ncol=n_quad,byrow = TRUE) %*%
          quad_weights)

      if(!log_out){
        out_temp <- outer(as.vector(Lambda0), exp(eta), "*")
      } else{
        out_temp <- outer(as.vector(log(Lambda0)), eta, "+")
      }

    }
  } else stop("'hazard' must be 'Weibull', 'Piecewise Constant', 'Royston-Parmar' or 'BSpline'")

  return(out_temp)
}


#function that "gets" jacobian matrix from inputs for univariate logh and logH
#to reduce overall size of output, this function only takes vector xnew rather than matrix

#' @export
get_jac <- function(tseq, xnew, beta, eta,
                    hazard, phi, basis, dbasis, basis_quad, quad_weights,
                    func_type, h, H, log_out){
  # browser()
  stopifnot(func_type %in% c("logh","h","logH","H","S"))
  # stopifnot(tseq[1]>0)

  #compute eta if it is not already provided
  if(is.null(eta) & length(beta) > 0){
    eta <- as.vector( xnew %*% beta)
  }

  #general formulae:

  #for a generic parameter theta,

  #dh/dtheta is some time by parameter matrix
  # dlog(h)/dtheta = 1/h * dh/dtheta columnwise

  #dH/dtheta is some time by parameter matrix
  # dlog(H)/dtheta = 1/H * dH/dtheta columnwise
  # dexp(-H)/dtheta = -H * dH/dtheta columnwise


  #now, generate the appropriate Jacobian matrix
  if(hazard %in% c("weibull","wb")){
    if(func_type %in% c("h","logh")){
      #conditional hazard jacobian for computing SE
      #log(haz) = logalpha + logkappa + (exp(logalpha)-1)*log(t) + xtbeta
      #Use delta method on log(haz) to compute baseline hazard confidence intervals
      #each row of J is gradient wrt logkappa, logalpha, beta1, ..., betap
      J <- cbind(1, #partial derivative of logkappa
                 (1 + exp(phi[2]) * log(tseq)), #partial derivative of logalpha
                 if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow=length(tseq),
                                                              ncol=length(xnew), byrow=T))
    } else {
      #under weibull, log(-log(S(t))) = log(H(t)) = logkappa + xtbeta + exp(logalpha)*log(t)
      #matrix with as many columns as parameters in S1, and as many rows as times in tseq
      #so, each row of J is the gradient at a particular tseq wrt logkappa, logalpha, beta1, ..., betap
      J <- cbind(1, #partial derivative of logkappa
                 exp(phi[2]) * log(tseq), #partial derivative of logalpha
                 if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow=length(tseq),
                                                              ncol=length(xnew), byrow=T))
    }
  } else if(hazard %in% c("piecewise","pw")){
    if(func_type %in% c("h","logh")){
      J <- cbind(t(t(dbasis) * exp(phi)) / as.vector(dbasis %*% exp(phi)),
                 if (!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                               ncol = length(xnew), byrow = T))
    } else {
      #under piecewise, log(-log(S(t))) = log( basis %*% exp(phi)) + xtbeta
      J <- cbind(t(t(basis) * exp(phi)) / as.vector(basis %*% exp(phi)), #partial deriv of phi
                 if (!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                               ncol = length(xnew), byrow = T))
    }
  } else if(hazard %in% c("royston-parmar","rp")){
    if(func_type %in% c("h","logh")){
      #log(h(t)) = log(s'(z)) + s(z) + xtbeta where z=log(t)
      J <- cbind(dbasis/as.vector(dbasis %*% phi) + basis,
                 if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                              ncol = length(xnew), byrow = T))
    } else {
      J <- cbind(basis,
                 if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                              ncol = length(xnew), byrow = T))
    }
  } else if(hazard %in% c("bspline","bs")){
    if(func_type %in% c("h","logh")){
      #log(haz) is Btphi + xtbeta
      J <- cbind(basis, #partial derivative of phi
                 if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                              ncol = length(xnew), byrow = T))
    } else {
      n_quad <- dim(basis_quad)[1] / dim(basis)[1]
      #matrix with as many columns as parameters in S1, and as many rows as times in tseq
      #under bspline, log(-log(S(t))) = log( [numerical integral of exp(Btphi)] ) + xtbeta
      #so, each row of J is the gradient at a particular tseq wrt phi1, ..., phiK, beta1, ..., betap
      #specifically, for gradient of phi we take each basis column and multiply it by lambda,
      #evaluate the numerical integral of the product, and then divide off Lambda0
      lambda0 <- as.vector(exp(basis_quad %*% phi))
      Lambda0 <- tseq/2 * as.vector(
        matrix(lambda0,
               ncol=n_quad,byrow = TRUE) %*%
          quad_weights)

      J_phi <- apply(X = basis_quad * lambda0, MARGIN = 2,
                     FUN = function(x) tseq/2 * as.vector(matrix(x,ncol=n_quad,byrow = TRUE) %*% quad_weights)) / Lambda0
      if(is.null(nrow(J_phi))){ J_phi <- t(J_phi)} #if J_phi is a vector, then make it a one row matrix instead...
      J <- cbind( J_phi,
                  if (!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                                ncol = length(xnew), byrow = T))
    }
  } else stop("'hazard' must be 'Weibull', 'Piecewise Constant', 'Royston-Parmar' or 'B-Spline'")

  #if needed, compute h or H
  if(func_type == "h"){
    if(is.null(h)){
      h <- get_haz(tseq=tseq, eta=eta, hazard=hazard, phi=phi,
                   basis=basis, dbasis=dbasis, basis_quad=basis_quad,
                   quad_weights=quad_weights, func_type="h", log_out=FALSE)
    }
    #update Jacobian to be on scale of h
    J <- h * J
  }
  if(func_type %in% c("H","S")){
    if(is.null(H)){
      H <- get_haz(tseq=tseq, eta=eta,hazard=hazard, phi=phi,
                   basis=basis, dbasis=dbasis, basis_quad=basis_quad,
                   quad_weights=quad_weights, func_type="H", log_out=FALSE)
    }
    #update Jacobian to be on scale of H
    if(func_type == "H") J <- H * J
    #update Jacobian to be on scale of S
    if(func_type == "S") J <- -H^2 * J
  }

  if(tseq[1] == 0){
    J[1,] <- 0
  }

  return(J)
}

#' @export
pred_helper_uni <- function(tseq, xnew=NULL, beta, phi, hazard, knots_vec,
                            n_quad, quad_method, vcov_sub, alpha, offset = 0){
  # browser()
  eta <- if(!is.null(xnew) & length(beta) > 0) as.vector(xnew %*% beta) + offset else offset

  basis <- get_basis(y = tseq,knots_vec = knots_vec,hazard = hazard, deriv = FALSE)
  dbasis <- get_basis(y = tseq,knots_vec = knots_vec,hazard = hazard, deriv = TRUE)
  quad_weights <- get_quad_pointsweights(n_quad=n_quad,
                     quad_method=quad_method)$weights
  quad_points <- transform_quad_points(n_quad = n_quad,
                    quad_method=quad_method, a = 0,b = tseq)
  basis_quad <- get_basis(y=quad_points, knots_vec=knots_vec,hazard=hazard,deriv = FALSE)

  #preallocate list (some of which will remain null, if excluded from out_options)
  value <- list(tseq=tseq, h=NULL,se_h=NULL,se_logh=NULL,ll_h=NULL,ul_h=NULL,
                H=NULL,se_H=NULL,se_logH=NULL,ll_H=NULL,ul_H=NULL,
                S=NULL,se_S=NULL,ll_S=NULL,ul_S=NULL, class=NULL)

  for(type in c("h","H")){
    value[[type]] <- get_haz(tseq=tseq, eta=eta,
                             hazard=hazard,
                             phi=phi, basis=basis, dbasis=dbasis,
                             basis_quad=basis_quad, quad_weights=quad_weights,
                             func_type=type, log_out=FALSE)

    for(log_out_temp in c("","log")){
      if(!is.null(vcov_sub)){
        value[[paste0("se_",log_out_temp,type)]] <- sapply(X = 1:length(eta),
             function(i){
               J <- get_jac(tseq=tseq, xnew=xnew[i,,drop=FALSE], beta=beta, eta=eta[i],
                            hazard=hazard,
                            phi=phi, basis=basis, dbasis=dbasis,
                            basis_quad=basis_quad, quad_weights=quad_weights,
                            func_type=paste0(log_out_temp,type),
                            log_out = (log_out_temp=="log"),
                            h=value[["h"]][,i], H=value[["H"]][,i])
               #https://stackoverflow.com/questions/27157127/efficient-way-of-calculating-quadratic-forms-avoid-for-loops
               sqrt(rowSums(J * (J %*% vcov_sub)))
             })
      }
    }
  }

  #avoid any weirdness that might happen if you start at 0
  if(tseq[1]==0){
    value[["H"]][1,] <- 0
  }


  #fill in survivor function using H, or if needed, computing it freshly
  value[["S"]] <- exp(-value$H)

  if(!is.null(vcov_sub)){
    #second, take the "alpha" given and generate confidence bounds
    #use log-scale CI by default
    value$ll_h <- value$h * exp(stats::qnorm(alpha/2) * value$se_logh)
    value$ul_h <- value$h * exp(-stats::qnorm(alpha/2) * value$se_logh)
    # value$ll_h <- value$h - stats::qnorm(alpha/2) * value$se_h
    # value$ul_h <- value$h + stats::qnorm(alpha/2) * value$se_h

    #use log-log-scale CI by default
    value$ll_H <- value$H * exp(stats::qnorm(alpha/2) * value$se_logH)
    value$ul_H <- value$H * exp(-stats::qnorm(alpha/2) * value$se_logH)
    value$ll_S <- exp(-value$ll_H)
    value$ul_S <- exp(-value$ul_H)
    # #Former log-scale CIs
    # value$ll_H <- value$H - stats::qnorm(alpha/2) * value$se_H
    # value$ul_H <- value$H + stats::qnorm(alpha/2) * value$se_H
    # value$ll_S <- exp(-value$ll_H)
    # value$ul_S <- exp(-value$ul_H)
    # #Former direct CI for survival
    # value$ll_S <- value$S - stats::qnorm(alpha/2) * value$se_S
    # value$ul_S <- value$S + stats::qnorm(alpha/2) * value$se_S
  }

  value
}


#First, a function that just predicts those things that do not require any amount of integration and minimial extra fuss.
#so h, H, S, and p_neither (without an SE).
#That way, there's minimal confusion or unexpected behavior based on
#the nature of tseq, and also it should be very fast

#Another reason for this is that if you want to show, say, h3 for particular choices of t1,
#you can directly include that in x3new, and this function doesn't need to `know` anything
#about how h3 depends on t1.
#this is in contrast to the risk profile function below, which will explicitly depend
#on the structure of dependence of h3 on t1 under the semi-markov model
#however, x3new needs to have the full set of columns corresponding to the effect of t1

#' @export
pred_helper_ID_simple <- function(tseq,
                                  x1new=NULL, x2new=NULL, x3new=NULL, frailnew=NULL,
                                  para, nP0, nP, frailty, model,
                                  hazard, knots_list, n_quad, quad_method,
                                  Finv, alpha){
  # browser()
  marg_flag <- is.null(frailnew)
  if(is.null(frailnew)){
    frailnew <- 1
  }
  se_fit_flag <- !is.null(Finv)

  #setting up some preliminary quantities
  nP0_tot <- if(frailty) sum(nP0) + 1 else sum(nP0)
  nP0_start <- 1 + c(0,nP0[1],nP0[1]+nP0[2])
  nP0_end <- c(nP0[1],nP0[1]+nP0[2],nP0[1]+nP0[2]+nP0[3])
  nP_start <- 1 + c(nP0_tot,nP0_tot+nP[1],nP0_tot+nP[1]+nP[2])
  nP_end <- c(nP0_tot+nP[1],nP0_tot+nP[1]+nP[2],nP0_tot+nP[1]+nP[2]+nP[3])

  quad_weights <- get_quad_pointsweights(n_quad=n_quad,
                                                             quad_method=quad_method)$weights
  quad_points <- transform_quad_points(n_quad = n_quad,
                                                           quad_method=quad_method, a = 0,b = tseq)

  #prefill the list with all potential inputs
  value_temp <- list()
  value <- list(tseq=tseq, tseq_h3_marg=NULL,
                h1=NULL,se_h1=NULL, se_logh1=NULL, ll_h1=NULL, ul_h1=NULL,
                h2=NULL,se_h2=NULL, se_logh2=NULL, ll_h2=NULL, ul_h2=NULL,
                h3=NULL,se_h3=NULL, se_logh3=NULL, ll_h3=NULL, ul_h3=NULL,
                H1=NULL,se_H1=NULL, se_logH1=NULL, ll_H1=NULL, ul_H1=NULL,
                H2=NULL,se_H2=NULL, se_logH2=NULL, ll_H2=NULL, ul_H2=NULL,
                H3=NULL,se_H3=NULL, se_logH3=NULL, ll_H3=NULL, ul_H3=NULL,
                S1=NULL,se_S1=NULL, ll_S1=NULL, ul_S1=NULL,
                S2=NULL,se_S2=NULL, ll_S2=NULL, ul_S2=NULL,
                S3=NULL,se_S3=NULL, ll_S3=NULL, ul_S3=NULL,
                x1new=x1new, x2new=x2new, x3new=x3new, frailnew=frailnew, class=NULL)
  #fill in the predicted quantities from the univariate prediction function
  for(i in 1:3){
    xtemp <- switch(i, x1new, x2new, x3new)

    if(se_fit_flag){
      Finv_temp <- Finv[c(nP0_start[i]:nP0_end[i],
                          if(!is.null(xtemp)) nP_start[i]:nP_end[i]),
                        c(nP0_start[i]:nP0_end[i],
                          if(!is.null(xtemp)) nP_start[i]:nP_end[i])]
    } else Finv_temp <- NULL

    value_temp <- pred_helper_uni(tseq = tseq,
                                  xnew = xtemp,
                                  phi = para[nP0_start[i]:nP0_end[i]],
                                  beta = para[nP_start[i]:nP_end[i]],
                                  hazard = hazard, knots_vec = knots_list[[i]],
                                  n_quad = n_quad, quad_method = quad_method,
                                  vcov_sub = Finv_temp, offset = log(frailnew),
                                  alpha=alpha) #as in, give me everything!
    #loop through and copy over the results from the univariate fit into the output list

    for(temp_name in names(value_temp)){
      if(temp_name=="tseq"){next}
      value[[paste0(temp_name,"",i)]] <- value_temp[[temp_name]]
    }
  }

  #if tseq=0, just fix these.
  #Then, it can be the default for tseq to start with 0
  if(tseq[1]==0){
    value$H1[1,] <- 0; value$H2[1,] <- 0; value$H3[1,] <- 0
    value$S1[1,] <- 1; value$S2[1,] <- 1; value$S3[1,] <- 1
    if(se_fit_flag){
      value$ll_H1[1,] <- 0; value$ll_H2[1,] <- 0; value$ll_H3[1,] <- 0
      value$ul_H1[1,] <- 0; value$ul_H2[1,] <- 0; value$ul_H3[1,] <- 0
      value$ll_S1[1,] <- 1; value$ll_S2[1,] <- 1; value$ll_S3[1,] <- 1
      value$ul_S1[1,] <- 1; value$ul_S2[1,] <- 1; value$ul_S3[1,] <- 1
    }
  }

  value
}


#' @export
predict.Freq_HReg2 <- function (object, xnew = NULL,
                                x1new = NULL, x2new = NULL, x3new = NULL,
                                tseq = seq(0,object$ymax,length.out = 100),
                                se.fit = TRUE, frailnew = NULL,
                                alpha = 0.05, n_quad=15, quad_method="kronrod", ...) {
  # browser()
  yLim <- NULL
  nP = object$nP
  nP0 = object$nP0
  if (object$class[2] == "Surv") {
    if (!(is.null(x1new) & is.null(x2new) & is.null(x3new))) {
      stop("'x1new','x2new', and 'x3new' are for semi-competing risks models and must be specified as NULL for univariate models")
    }

    value <- pred_helper_uni(tseq=tseq, xnew=xnew,
                             beta=object$estimate[-(1:nP0)],
                             phi=object$estimate[1:nP0],
                             hazard=object$hazard,
                             knots_vec=object$knots_vec,
                             n_quad=n_quad, quad_method=quad_method,
                             vcov_sub= if(se.fit) object$Finv[c(1:nP0, if(!is.null(xnew)) (nP0+1):(nP0+nP)),
                                                              c(1:nP0, if(!is.null(xnew)) (nP0+1):(nP0+nP))] else NULL,
                             alpha = alpha)
    value$xnew <- xnew #not in the univariate helper function so needs to be added

  } else if(object$class[2]=="ID"){
    if (!is.null(xnew)) {
      stop("'xnew' is for univariate models so it must be specified as NULL for semi-competing risks models")
    }

    if(!is.null(x1new) & !is.null(x2new)) stopifnot(NROW(x1new) == NROW(x2new))
    if(!is.null(x1new) & !is.null(x3new)) stopifnot(NROW(x1new) == NROW(x3new))
    if(!is.null(x2new) & !is.null(x3new)) stopifnot(NROW(x2new) == NROW(x3new))

    value <- pred_helper_ID_simple(tseq=tseq, x1new=x1new, x2new=x2new, x3new=x3new, frailnew=frailnew,
                                   para=object$estimate, nP0=object$nP0, nP=object$nP,
                                   frailty=object$frailty, model=object$model,
                                   hazard=object$hazard, knots_list=object$knots_list,
                                   n_quad=n_quad, quad_method=quad_method,
                                   Finv= if(se.fit) object$Finv else NULL, alpha = alpha)

  }

  value$class <- object$class
  class(value) <- "pred.Freq_HReg2"
  return(value)
}

#' @export
plot.pred.Freq_HReg2 <- function (x, plot.est = "Haz",
                                  xlab = NULL, ylab = NULL, ci=TRUE, ...){
  #color blind friendly colors from here: https://davidmathlogic.com/colorblind/#%23648FFF-%23785EF0-%23DC267F-%23FE6100-%23FFB000
  cb_blue <- "#648FFF"; cb_red <- "#DC267F"; cb_purple <- "#785EF0";
  cb_orange <- "#FE6100"; cb_grey <- "#CACACA"

  # browser()

  if (x$class[2] == "Surv") {
    if (is.null(ylab)){
      ylab <- switch(tolower(plot.est),
                     "surv"="Survival", "s"="Survival",
                     "haz"="Hazard", "h"="Hazard",
                     "cumhaz"="Cumulative Hazard","ch"="Cumulative Hazard")
    }
    if(is.null(xlab)){xlab <- "Time"}

    if (tolower(plot.est) %in% c("surv","s")) {
      yLim <- c(0,1)
      main_temp <- expression(paste("Estimated ", S(t), ""))
    }
    if (tolower(plot.est) %in% c("haz","h")) {
      yLim <- if(ci) c(0,max(x$ul_h[is.finite(x$ul_h)])) else c(0,max(x$h[is.finite(x$h)]))
      main_temp <- expression(paste("Estimated ", h(t), ""))
    }
    if (tolower(plot.est) %in% c("cumhaz","ch")){
      yLim <- if(ci) c(0,max(x$ul_H[is.finite(x$ul_H)])) else c(0,max(x$H[is.finite(x$H)]))
      main_temp <- expression(paste("Estimated ", H(t), ""))
    }

    plot_vec <- switch(tolower(plot.est), "surv"=x$S, "s"=x$S,
                       "haz"=x$h, "h"=x$h, "cumhaz"=x$H, "ch"=x$H)
    graphics::matplot(x$tseq, plot_vec, type="l", lwd=3, lty=1,
            main = main_temp,
            xlim=c(0,max(x$tseq)), ylim = yLim, xlab=xlab, ylab=ylab)
    if(ci){
      plot_vec <- switch(tolower(plot.est), "surv"=x$ll_S, "s"=x$ll_S,
                         "haz"=x$ll_h, "h"=x$ll_h, "cumhaz"=x$ll_H, "ch"=x$ll_H)
      graphics::matplot(x$tseq, plot_vec, type="l", lwd=3, lty=3, add=T)
      plot_vec <- switch(tolower(plot.est), "surv"=x$ul_S, "s"=x$ul_S,
                         "haz"=x$ul_h, "h"=x$ul_h, "cumhaz"=x$ul_H, "ch"=x$ul_H)
      graphics::matplot(x$tseq, plot_vec, type="l", lwd=3, lty=3, add=T)
    }
  }

  if (x$class[2] == "ID") {

    if (is.null(xlab)) {
      xlab <- c("Time", "Time", "Time")
      if (x$class[5] == "semi-Markov") {
        xlab[3] <- "Time since non-terminal event"
      }
    }
    if (is.null(ylab)){
      ylab <- switch(tolower(plot.est),
                     "surv"=c("Cause-Specific Survival","Cause-Specific Survival","Conditional Survival"),
                     "s"=c("Cause-Specific Survival","Cause-Specific Survival","Conditional Survival"),
                     "haz"=c("Cause-Specific Hazard","Cause-Specific Hazard","Conditional Hazard"),
                     "h"=c("Cause-Specific Hazard","Cause-Specific Hazard","Conditional Hazard"),
                     "cumhaz"=c("Cause-Specific Cumulative Hazard","Cause-Specific Cumulative Hazard","Conditional Cumulative Hazard"),
                     "ch"=c("Cause-Specific Cumulative Hazard","Cause-Specific Cumulative Hazard","Conditional Cumulative Hazard"))
    }

    main_temp <- vector(mode = "list", length = 3)
    if (tolower(plot.est) %in% c("surv","s")) {
      yLim <- c(0,1)
      main_temp[[1]] <- expression(paste("Estimated ", S[1](t), ""))
      main_temp[[2]] <- expression(paste("Estimated ", S[2](t), ""))
      main_temp[[3]] <- expression(paste("Estimated ", S[3](t), ""))
    }
    if (tolower(plot.est) %in% c("haz","h")) {
      if(ci){
        maxy <- max(x$ul_h1[is.finite(x$ul_h1)],
                    x$ul_h2[is.finite(x$ul_h2)],
                    x$ul_h3[is.finite(x$ul_h3)])
      }  else{
        maxy <- max(x$h1[is.finite(x$h1)],
                    x$h2[is.finite(x$h2)],
                    x$h3[is.finite(x$h3)])
      }
      yLim <- c(0,maxy)
      main_temp[[1]] <- expression(paste("Estimated ", h[1](t), ""))
      main_temp[[2]] <- expression(paste("Estimated ", h[2](t), ""))
      main_temp[[3]] <- expression(paste("Estimated ", h[3](t), ""))
    }
    if (tolower(plot.est) %in% c("cumhaz","ch")){
      if(ci){
        maxy <- max(x$ul_H1[is.finite(x$ul_H1)],
                    x$ul_H2[is.finite(x$ul_H2)],
                    x$ul_H3[is.finite(x$ul_H3)])
      }  else{
        maxy <- max(x$H1[is.finite(x$H1)],
                    x$H2[is.finite(x$H2)],
                    x$H3[is.finite(x$H3)])
      }
      yLim <- c(0,maxy)
      main_temp[[1]] <- expression(paste("Estimated ", H[1](t), ""))
      main_temp[[2]] <- expression(paste("Estimated ", H[2](t), ""))
      main_temp[[3]] <- expression(paste("Estimated ", H[3](t), ""))
    }

    graphics::par(mfrow = c(1, 3))

    plot_vec <- switch(tolower(plot.est), "surv"=x$S1, "s"=x$S1,
                       "haz"=x$h1, "h"=x$h1, "cumhaz"=x$H1, "ch"=x$H1)
    graphics::matplot(x$tseq, plot_vec, type="l", lwd=3, lty=1,
            main = main_temp[[1]],
            xlim=c(0,max(x$tseq)), ylim = yLim, xlab=xlab[1], ylab=ylab[1])
    if(ci){
      plot_vec <- switch(tolower(plot.est), "surv"=x$ll_S1, "s"=x$ll_S1,
                         "haz"=x$ll_h1, "h"=x$ll_h1, "cumhaz"=x$ll_H1, "ch"=x$ll_H1)
      graphics::matplot(x$tseq, plot_vec, type="l", lwd=3, lty=3, add=T)
      plot_vec <- switch(tolower(plot.est), "surv"=x$ul_S1, "s"=x$ul_S1,
                         "haz"=x$ul_h1, "h"=x$ul_h1, "cumhaz"=x$ul_H1, "ch"=x$ul_H1)
      graphics::matplot(x$tseq, plot_vec, type="l", lwd=3, lty=3, add=T)
    }

    plot_vec <- switch(tolower(plot.est), "surv"=x$S2, "s"=x$S2,
                       "haz"=x$h2, "h"=x$h2, "cumhaz"=x$H2, "ch"=x$H2)
    graphics::matplot(x$tseq, plot_vec, type="l", lwd=3, lty=1,
            main = main_temp[[2]],
            xlim=c(0,max(x$tseq)), ylim = yLim, xlab=xlab[2], ylab=ylab[2])
    if(ci){
      plot_vec <- switch(tolower(plot.est), "surv"=x$ll_S2, "s"=x$ll_S2,
                         "haz"=x$ll_h2, "h"=x$ll_h2, "cumhaz"=x$ll_H2, "ch"=x$ll_H2)
      graphics::matplot(x$tseq, plot_vec, type="l", lwd=3, lty=3, add=T)
      plot_vec <- switch(tolower(plot.est), "surv"=x$ul_S2, "s"=x$ul_S2,
                         "haz"=x$ul_h2, "h"=x$ul_h2, "cumhaz"=x$ul_H2, "ch"=x$ul_H2)
      graphics::matplot(x$tseq, plot_vec, type="l", lwd=3, lty=3, add=T)
    }

    plot_vec <- switch(tolower(plot.est), "surv"=x$S3, "s"=x$S3,
                       "haz"=x$h3, "h"=x$h3, "cumhaz"=x$H3, "ch"=x$H3)
    graphics::matplot(x$tseq, plot_vec, type="l", lwd=3, lty=1,
            main = main_temp[[3]],
            xlim=c(0,max(x$tseq)), ylim = yLim, xlab=xlab[3], ylab=ylab[3])
    if(ci){
      plot_vec <- switch(tolower(plot.est), "surv"=x$ll_S3, "s"=x$ll_S3,
                         "haz"=x$ll_h3, "h"=x$ll_h3, "cumhaz"=x$ll_H3, "ch"=x$ll_H3)
      graphics::matplot(x$tseq, plot_vec, type="l", lwd=3, lty=3, add=T)
      plot_vec <- switch(tolower(plot.est), "surv"=x$ul_S3, "s"=x$ul_S3,
                         "haz"=x$ul_h3, "h"=x$ul_h3, "cumhaz"=x$ul_H3, "ch"=x$ul_H3)
      graphics::matplot(x$tseq, plot_vec, type="l", lwd=3, lty=3, add=T)
    }

    graphics::par(mfrow = c(1, 1))
  }
  invisible()
}


#ok, so here is the last piece: the full on risk prediction function, which
#depends on the full sequence of tseq values, allows h3 to vary with t1, and
#admits left truncation.


#ok, so I have the univariate prediction function done!
#now, the motherlode: illness-death model prediction
#the only thing missing from this now is deciding how to incorporate
#dependence of h3 on t1.

#' @export
pred_risk_ID <- function(tseq,
                            x1new=NULL, x2new=NULL, x3new=NULL, frailnew=NULL,
                            para, nP0, nP, frailty, model,
                            p3tv, h3tv_basis_func,
                            hazard, knots_list, n_quad, quad_method,
                            Finv, alpha){
  # browser()
  marg_flag <- is.null(frailnew)
  if(is.null(frailnew)){
    frailnew <- 1
  }

  #it's crazy but I think I'll specify things related to the effect of t1 on h3 as follows:
  #you tell me how many total parameters are used in estimating the h3tv effect
  #and you give me the function that defines the basis for that effect.
  #Could be linear, piecewise, spline, I don't care. All I know is that it's an effect on the log-hazard scale.
  #also, assume that it's the last thing in the vector of parameters

  #setting up some preliminary quantities
  nP0_tot <- if(frailty) sum(nP0) + 1 else sum(nP0)
  nP0_start <- 1 + c(0,nP0[1],nP0[1]+nP0[2])
  nP0_end <- c(nP0[1],nP0[1]+nP0[2],nP0[1]+nP0[2]+nP0[3])
  nP_start <- 1 + c(nP0_tot,nP0_tot+nP[1],nP0_tot+nP[1]+nP[2])
  nP_end <- c(nP0_tot+nP[1],nP0_tot+nP[1]+nP[2],nP0_tot+nP[1]+nP[2]+nP[3]-p3tv) #leave off the h3tv effects!!!

  eta_list <- beta_list <- vector(mode="list",3)
  beta_list[[1]] <- if(!is.null(x1new)) para[nP_start[1]:nP_end[1]] else numeric(0)
  beta_list[[2]] <- if(!is.null(x2new)) para[nP_start[2]:nP_end[2]] else numeric(0)
  beta_list[[3]] <- if(!is.null(x3new)) para[nP_start[3]:nP_end[3]] else numeric(0)
  eta_list[[1]] <- if(length(beta_list[[1]]) > 0) as.vector(x1new %*% beta_list[[1]] + log(frailnew)) else 0
  eta_list[[2]] <- if(length(beta_list[[2]]) > 0) as.vector(x2new %*% beta_list[[2]] + log(frailnew)) else 0
  eta_list[[3]] <- if(length(beta_list[[3]]) > 0) as.vector(x3new %*% beta_list[[3]] + log(frailnew)) else 0

  h3tv_phi <- if(p3tv>0) utils::tail(para,n=p3tv) else numeric(0)

  se_fit_flag <- !is.null(Finv)
  if(frailty){
    theta <- exp(para[sum(nP0) + 1])
  }

  #prefill the list with all potential inputs
  value_temp <- list()
  value <- list(tseq=tseq,
                p_term_only=NULL, se_p_term_only=NULL,
                p_nonterm=NULL, se_p_nonterm=NULL,
                p_nonterm_only=NULL, se_p_nonterm_only=NULL,
                p_both=NULL, se_p_both=NULL,
                p_neither=NULL, se_p_neither=NULL, se_logp_neither=NULL,
                p_term_only_marg=NULL, se_p_term_only_marg=NULL,
                p_nonterm_marg=NULL, se_p_nonterm_marg=NULL,
                p_nonterm_only_marg=NULL, se_p_nonterm_only_marg=NULL,
                p_both_marg=NULL, se_p_both_marg=NULL,
                p_neither_marg=NULL, se_p_neither_marg=NULL,
                x1new=x1new, x2new=x2new, x3new=x3new, frailnew=frailnew, class=NULL)

  quad_weights <- SemiCompRisksFreq:::get_quad_pointsweights(n_quad=n_quad,
                                                             quad_method=quad_method)$weights
  quad_points <- SemiCompRisksFreq:::transform_quad_points(n_quad = n_quad,
                                                           quad_method=quad_method, a = 0,b = tseq)
  H_list <-  vector(mode="list",3)
  for(i in 1:3){
    xtemp <- switch(i, x1new, x2new, x3new)

    basis <- get_basis(y = tseq,knots_vec = knots_list[[i]],hazard = hazard, deriv = FALSE)
    dbasis <- get_basis(y = tseq,knots_vec = knots_list[[i]],hazard = hazard, deriv = TRUE)
    basis_quad <- get_basis(y=quad_points, knots_vec=knots_list[[i]],hazard=hazard,deriv = FALSE)

    #don't bother with h3_tv here because H_list[[3]] is only used by semi-markov setting
    H_list[[i]] <- get_haz(tseq=tseq, eta=eta_list[[i]],
                           hazard=hazard, phi=para[nP0_start[i]:nP0_end[i]],
                           basis=basis, dbasis=dbasis,
                           basis_quad=basis_quad, quad_weights=quad_weights,
                           func_type="H", log_out=FALSE)
  }

  # browser()

  #if tseq=0, just fix these.
  #Then, it can be the default for tseq to start with 0
  if(tseq[1]==0){
    H_list[[1]][1,] <- 0
    H_list[[2]][1,] <- 0
    H_list[[3]][1,] <- 0
  }

  #next, compute conditional probabilities
  #note, for the moment these are length(tseq) by length(eta) matrices!
  value$p_neither <- exp(-H_list[[1]] - H_list[[2]])
  #make it so that you've "started integrating" from time tseq[1]
  #in practice, this is no different if tseq = 0
  p_neither0 <- exp(-H_list[[1]][1,] - H_list[[2]][1,])
  value$p_neither <- t(t(value$p_neither) / p_neither0)

  #instead of assuming that there's an implicit 0, just directly use the first thing!
  #this allows for left truncation by just beginning the integral at an arbitrary time
  S1_b <- exp(-H_list[[1]][-1,,drop=FALSE])
  S1_a <- exp(-H_list[[1]][-length(tseq),,drop=FALSE])
  S2_b <- exp(-H_list[[2]][-1,,drop=FALSE])
  S2_a <- exp(-H_list[[2]][-length(tseq),,drop=FALSE])

  #and then in turn, we know that wherever the integral is starting from, it's
  #gotta start at 0, so just add that back
  value$p_term_only <- rbind(0, 0.5 *
    apply( (S1_a + S1_b) * (S2_a - S2_b), MARGIN = 2, cumsum))
  value$p_nonterm <- rbind(0, 0.5 *
    apply( (S1_a - S1_b) * (S2_a + S2_b), MARGIN = 2, cumsum))

  #if there is left truncation, then correct for it by rescaling by the probability
  #of no event by the starting timepoint
  value$p_nonterm <- t(t(value$p_nonterm) / p_neither0)
  value$p_term_only <- t(t(value$p_term_only) / p_neither0)

  if(frailty & marg_flag){
    #Following Xu (2010), marginal hazards reflect interplay of H1 and H2
    #this could be computed on log-scale, fyi! could be useful for some reason
    value$p_neither_marg <- exp( (-1/theta) * log1p(theta * (H_list[[1]] + H_list[[2]])) )
    #make it so that you've "started integrating" from time tseq[1]
    #in practice, this is no different if tseq = 0
    p_neither0_marg <- exp( (-1/theta) * log1p(theta * (H_list[[1]][1,] + H_list[[2]][1,])) )
    value$p_neither_marg <- t(t(value$p_neither_marg) / p_neither0_marg)

    H1_b <- H_list[[1]][-1,,drop=FALSE]
    H1_a <- H_list[[1]][-length(tseq),,drop=FALSE]
    H2_b <- H_list[[2]][-1,,drop=FALSE]
    H2_a <- H_list[[2]][-length(tseq),,drop=FALSE]

    #Here, we interchange the integral of gamma with the summation of trapezoid rule
    #to get sum of "integrated integrands"
    value$p_term_only_marg <- rbind(0, 0.5 * apply(
      (1 + theta * (H1_a + H2_a))^(-1/theta) +
        (1 + theta * (H1_b + H2_a))^(-1/theta) -
        (1 + theta * (H1_a + H2_b))^(-1/theta) -
        (1 + theta * (H1_b + H2_b))^(-1/theta), MARGIN=2, cumsum))
    value$p_nonterm_marg <- rbind(0, 0.5 * apply(
      (1 + theta * (H1_a + H2_a))^(-1/theta) +
        (1 + theta * (H1_a + H2_b))^(-1/theta) -
        (1 + theta * (H1_b + H2_a))^(-1/theta) -
        (1 + theta * (H1_b + H2_b))^(-1/theta), MARGIN=2, cumsum))

    value$p_term_only_marg <- t(t(value$p_term_only_marg) / p_neither0_marg)
    value$p_nonterm_marg <- t(t(value$p_nonterm_marg) / p_neither0_marg)

    # #Below is a log-scale implementation of these same computations,
    # #I'm sure we could do this throughout but I'm not sure how essential it is...
    # value$p_term_only_marg <- 0.5 * exp(apply(
    #   logdiffexp(
    #     logsumexp((-1/theta)*log1p(theta * (H1_a + H2_a)),
    #                  (-1/theta)*log1p(theta * (H1_b + H2_a))),
    #     logsumexp((-1/theta)*log1p(theta * (H1_a + H2_b)),
    #                  (-1/theta)*log1p(theta * (H1_b + H2_b)))
    #   ), MARGIN=2, logcumsumexp))
    # value$p_term_only_marg <- 0.5 * exp(apply(
    #   logdiffexp(
    #     logsumexp((-1/theta)*log1p(theta * (H1_a + H2_a)),
    #                  (-1/theta)*log1p(theta * (H1_a + H2_b))),
    #     logsumexp((-1/theta)*log1p(theta * (H1_b + H2_a)),
    #                  (-1/theta)*log1p(theta * (H1_b + H2_b)))
    #   ), MARGIN=2, logcumsumexp))

  }

  #now compute SEs for these quantities
  #loop through "individuals" (I know this is woefully inefficient but idk what to say...)
  if(se_fit_flag){
    Finv_sub12 <- Finv[c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                         nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2]),
                       c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                         nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2])]

    #while I'm here, set up the matrix for the more complicated  probabilities in the loop far below
    Finv_sub123 <- Finv[c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                          nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2],
                          nP0_start[3]:nP0_end[3], if(!is.null(x3new)) nP_start[3]:(nP_end[3] + p3tv)),
                        c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                          nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2],
                          nP0_start[3]:nP0_end[3], if(!is.null(x3new)) nP_start[3]:(nP_end[3] + p3tv))]

    if(frailty & marg_flag){
      Finv_sub12_marg <- Finv[c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                                nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2],
                                nP0_tot), #add frailty
                              c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                                nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2],
                                nP0_tot)] #add frailty

      #while I'm here, set up the matrix for the more complicated probabilities in the loop far below
      Finv_sub123_marg <- Finv[c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                                 nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2],
                                 nP0_start[3]:nP0_end[3], if(!is.null(x3new)) nP_start[3]:(nP_end[3] + p3tv),
                                 nP0_tot),
                               c(nP0_start[1]:nP0_end[1], if(!is.null(x1new)) nP_start[1]:nP_end[1],
                                 nP0_start[2]:nP0_end[2], if(!is.null(x2new)) nP_start[2]:nP_end[2],
                                 nP0_start[3]:nP0_end[3], if(!is.null(x3new)) nP_start[3]:(nP_end[3] + p3tv),
                                 nP0_tot)]
    }
    value$se_p_neither_marg <-
      value$se_p_neither <- value$se_logp_neither <-
      value$se_p_nonterm <- value$se_p_term_only <-
      value$se_p_nonterm_marg <- value$se_p_term_only_marg <-
      matrix(data=NA, nrow=NROW(H_list[[1]]),ncol=NCOL(H_list[[1]]))

    #again set up containers for the p values to be computed in the loop far below
    value$se_p_both_marg <- value$se_p_nonterm_only_marg <-
      value$se_p_both <- value$se_p_nonterm_only <-
      matrix(data = NA, nrow = length(tseq), ncol=NCOL(H_list[[1]]))


    basis1 <- get_basis(y = tseq,knots_vec = knots_list[[1]],hazard = hazard,deriv = FALSE)
    dbasis1 <- get_basis(y = tseq,knots_vec = knots_list[[1]],hazard = hazard,deriv = TRUE)
    basis1_quad <- get_basis(y=quad_points, knots_vec=knots_list[[1]],hazard=hazard,deriv = FALSE)
    basis2 <- get_basis(y = tseq,knots_vec = knots_list[[2]],hazard = hazard,deriv = FALSE)
    dbasis2 <- get_basis(y = tseq,knots_vec = knots_list[[2]],hazard = hazard,deriv = TRUE)
    basis2_quad <- get_basis(y=quad_points, knots_vec=knots_list[[2]],hazard=hazard,deriv = FALSE)

    for(j in 1:NCOL(H_list[[1]])){
      dH1b <- get_jac(tseq=tseq, xnew=x1new[j,,drop=FALSE], beta=beta_list[[1]],
                      eta=eta_list[[1]][j], #WILL THIS CAUSE A PROBLEM IF ONE X HAS NO ELEMENTS BUT ANOTHER HAS SEVERAL?
                      hazard=hazard, phi=para[nP0_start[1]:nP0_end[1]],
                      basis=basis1, dbasis=dbasis1,
                      basis_quad=basis1_quad, quad_weights=quad_weights,
                      func_type="H", H = H_list[[1]][,j],
                      log_out = FALSE )
      dH2b <- get_jac(tseq=tseq, xnew=x2new[j,,drop=FALSE], beta=beta_list[[2]],
                      eta=eta_list[[2]][j], #WILL THIS CAUSE A PROBLEM IF ONE X HAS NO ELEMENTS BUT ANOTHER HAS SEVERAL?
                      hazard=hazard, phi=para[nP0_start[2]:nP0_end[2]],
                      basis=basis2, dbasis=dbasis2,
                      basis_quad=basis2_quad, quad_weights=quad_weights,
                      func_type="H", H = H_list[[2]][,j],
                      log_out = FALSE )

      # does this cause a problem depending on left truncation??
      if(tseq[1]==0){
        dH1b[1,] <- 0
        dH2b[1,] <- 0
      }
      dH1a <- dH1b[-NROW(dH1b),,drop=FALSE]
      dH2a <- dH2b[-NROW(dH2b),,drop=FALSE]
      dH1b <- dH1b[-1,,drop=FALSE]
      dH2b <- dH2b[-1,,drop=FALSE]

      #see scratchwork for this derivation, basically it's leveraging the summation of trapezoid rule
      #for p_term_only (aka CIF_t)
      #apply(S1_a*S2_a + S1_b*S2_a - S1_a*S2_b - S1_b*S2_b, MARGIN = 2, cumsum)
      J <- 0.5 * cbind(
        dH1a*S1_a[,j]*(S2_a[,j] - S2_b[,j]) + dH1b*S1_b[,j]*(S2_a[,j] - S2_b[,j]),
        dH2a*S2_a[,j]*(S1_a[,j] + S1_b[,j]) - dH2b*S2_b[,j]*(S1_a[,j] + S1_b[,j])
      )
      value$se_p_term_only[,j] <- c(0,sqrt(diag(apply(apply(J %*% Finv_sub12 %*% t(J),1,cumsum),1,cumsum))))

      #for p_nonterm (aka CIF_nt)
      #apply(S2_a*S1_a + S2_b*S1_a - S2_a*S1_b - S2_b*S1_b, MARGIN = 2, cumsum)
      J <- 0.5 * cbind(
        dH1a*S1_a[,j]*(S2_a[,j] + S2_b[,j]) - dH1b*S1_b[,j]*(S2_a[,j] + S2_b[,j]),
        dH2a*S2_a[,j]*(S1_a[,j] - S1_b[,j]) + dH2b*S2_b[,j]*(S1_a[,j] - S1_b[,j])
      )
      value$se_p_nonterm[,j] <- c(0, sqrt(diag(apply(apply(J %*% Finv_sub12 %*% t(J),1,cumsum),1,cumsum))))

      #SE of "neither" probability
      J <- -cbind(dH1b, dH2b)
      value$se_logp_neither[,j] <- c(0, sqrt(rowSums(J * (J %*% Finv_sub12))))
      J <- value$p_neither[-1,j] * J
      value$se_p_neither[,j] <- c(0, sqrt(rowSums(J * (J %*% Finv_sub12))))

      if(frailty & marg_flag){
        #for p_term_only (aka CIF_t)
        #apply(S1_a*S2_a + S1_b*S2_a - S1_a*S2_b - S1_b*S2_b, MARGIN = 2, cumsum)
        J <- 0.5 * cbind(
          dH1a * ((1 + theta * (H1_a[,j] + H2_a[,j]))^(-1/theta - 1) -
                  (1 + theta * (H1_a[,j] + H2_b[,j]))^(-1/theta - 1)) +
          dH1b * ((1 + theta * (H1_b[,j] + H2_a[,j]))^(-1/theta - 1) -
                  (1 + theta * (H1_b[,j] + H2_b[,j]))^(-1/theta - 1)),
          dH2a * ((1 + theta * (H1_a[,j] + H2_a[,j]))^(-1/theta - 1) +
                  (1 + theta * (H1_b[,j] + H2_a[,j]))^(-1/theta - 1)) -
          dH2b * ((1 + theta * (H1_a[,j] + H2_b[,j]))^(-1/theta - 1) +
                  (1 + theta * (H1_b[,j] + H2_b[,j]))^(-1/theta - 1)),
          -marg_logtheta_deriv(a = H1_a[,j] + H2_a[,j], ltheta = log(theta)) +
            -marg_logtheta_deriv(a = H1_b[,j] + H2_a[,j], ltheta = log(theta)) -
            -marg_logtheta_deriv(a = H1_a[,j] + H2_b[,j], ltheta = log(theta)) -
            -marg_logtheta_deriv(a = H1_b[,j] + H2_b[,j], ltheta = log(theta)))
        value$se_p_term_only_marg[,j] <- c(0, sqrt(diag(apply(apply(J %*% Finv_sub12_marg %*% t(J),1,cumsum),1,cumsum))))

        #for p_nonterm (aka CIF_nt)
        #apply(S2_a*S1_a + S2_b*S1_a - S2_a*S1_b - S2_b*S1_b, MARGIN = 2, cumsum)
        J <- 0.5 * cbind(
          dH1a * ((1 + theta * (H1_a[,j] + H2_a[,j]))^(-1/theta - 1) +
                  (1 + theta * (H1_a[,j] + H2_b[,j]))^(-1/theta - 1)) -
          dH1b * ((1 + theta * (H1_b[,j] + H2_a[,j]))^(-1/theta - 1) +
                  (1 + theta * (H1_b[,j] + H2_b[,j]))^(-1/theta - 1)),
          dH2a * ((1 + theta * (H1_a[,j] + H2_a[,j]))^(-1/theta - 1) -
                  (1 + theta * (H1_b[,j] + H2_a[,j]))^(-1/theta - 1)) +
          dH2b * ((1 + theta * (H1_a[,j] + H2_b[,j]))^(-1/theta - 1) -
                  (1 + theta * (H1_b[,j] + H2_b[,j]))^(-1/theta - 1)),
          -marg_logtheta_deriv(a = H1_a[,j] + H2_a[,j], ltheta = log(theta)) +
            -marg_logtheta_deriv(a = H1_a[,j] + H2_b[,j], ltheta = log(theta)) -
            -marg_logtheta_deriv(a = H1_b[,j] + H2_a[,j], ltheta = log(theta)) -
            -marg_logtheta_deriv(a = H1_b[,j] + H2_b[,j], ltheta = log(theta)))
        value$se_p_nonterm_marg[,j] <- c(0, sqrt(diag(apply(apply(J %*% Finv_sub12_marg %*% t(J),1,cumsum),1,cumsum))))

        #SE of marginal "neither" probability
        J <- cbind(dH1b * (1 + theta * (H1_b[,j] + H2_b[,j]))^(-1/theta - 1),
                   dH2b * (1 + theta * (H1_b[,j] + H2_b[,j]))^(-1/theta - 1),
                   -marg_logtheta_deriv(a = H1_b[,j] + H2_b[,j], ltheta = log(theta)))
        value$se_p_neither_marg[,j] <- c(0, sqrt(rowSums(J * (J %*% Finv_sub12_marg))))
      }
    }
  }

  # browser()

  #now, to handle the quantities that depend on h3...

  #I think now, we'll have to use loops to incorporate h3 correctly
  value$p_both_marg <- value$p_nonterm_only_marg <-
    value$p_both <- value$p_nonterm_only <-
    matrix(data = NA, nrow = length(tseq),
           ncol=NCOL(H_list[[1]]))

  # if(abs(max(diff(tseq)) - min(diff(tseq))) > 1e-6){
  #   warning("tseq points are not equally spaced. Some predicted functions may be incorrect.")
  # }

  value$p_nonterm_only[1,] <- 0
  value$p_both[1,] <- 0
  if(frailty & marg_flag){
    value$p_nonterm_only_marg[1,] <- 0
    value$p_both_marg[1,] <- 0
  }

  # #and further assume that t1 is not included as a covariate...
  # if(tolower(model) %in% c("m","markov")){
  #   H3_b <- value$H3[-1,,drop=FALSE]
  #   H3_a <- value$H3[-length(tseq),,drop=FALSE]
  #   value$p_nonterm_only <- exp(value$H3) * rbind(0, 0.5 *
  #             apply(S2_a*exp(-H3_a)*S1_a + S2_b*exp(-H3_b)*S1_a -
  #                     S2_a*exp(-H3_a)*S1_b - S2_b*exp(-H3_b)*S1_b, MARGIN = 2, cumsum))
  # }



  if(p3tv > 0){
    h3_tv_basis <- h3tv_basis_func(tseq)
    h3_tv_mult <- exp( as.vector(h3_tv_basis %*% h3tv_phi) )
  } else{
    h3_tv_mult <- rep(1,length(tseq))
    h3_tv_basis <- NULL
  }




  for(i in 1:(length(tseq)-1)){

    # browser()

    #and further assume that t1 is not included as a covariate...
    if(tolower(model) %in% c("m","markov")){
      tseq_temp <- tseq[1:(i+1)]
      #last row of H3[1:i] is H3(t2i), so
      #compute matrix of H3(t2) - H3(t1) for every t1 between 0 and t1i
      #by "sweeping" out H3(t2i) from every row of -H3(t1) matrix

      # H3_a_sub <- t( -t(value$H3[1:i,,drop=FALSE]) + value$H3[i,])

      #YAY THIS UPDATE WORKS!! EVEN FOR LEFT TRUNCATED and IRREGULARLY SPACE DATA!
      H3_temp <- t( -t(H_list[[3]][1:(i+1),,drop=FALSE]) + H_list[[3]][(i+1),])
      H3_b_sub <- H3_temp[-1,,drop=FALSE]
      H3_a_sub <- H3_temp[-(i+1),,drop=FALSE]

    } else{
      tseq_temp <- tseq[(i+1)] - tseq[1:(i+1)]

      #semi-markov is trickier and I think may in fact need to be re-predicted in general
      #both to account for irregular spacing and for dependence on t_1

      #take the first i entries, and then reverse their order
      #to get a vector that increments from H3(t2i-0) to H3(t2i-t2i)

      #OK THIS WORKS FOR NON-TRUNCATED CASE! But not truncated case
      #because I think we'd need to compute new values of H3 regardless! i.e., H_3(5.5-5)
      #is not already computed if we are only computing H_3 on the range of 5 to 6
      # temp <- value$H3[1:(i+1),,drop=FALSE][(i+1):1,,drop=FALSE]
      # H3_b_sub <- temp[-1,,drop=FALSE]
      # H3_a_sub <- temp[-(i+1),,drop=FALSE]

      #reprediction solves the above issue
      #only lingering issue is how to bring in dependence on t1
      H3_temp <- pred_helper_uni(tseq = tseq_temp,
                                 xnew = x3new,
                                 phi = para[nP0_start[3]:nP0_end[3]],
                                 beta = para[nP_start[3]:nP_end[3]],
                                 hazard = hazard, knots_vec = knots_list[[3]],
                                 n_quad = n_quad, quad_method = quad_method,
                                 vcov_sub = NULL, offset = log(frailnew),
                                 alpha=alpha)$H * h3_tv_mult[1:(i+1)] #NOTE COLUMNWISE MULTIPLICATION BY T1 EFFECT FACTOR

      #the final one should always be zero
      H3_temp[i+1,] <- 0
      H3_b_sub <- H3_temp[-1,,drop=FALSE]
      H3_a_sub <- H3_temp[-(i+1),,drop=FALSE]
    }

    S1_a_sub <- S1_a[1:i,,drop=FALSE]
    S1_b_sub <- S1_b[1:i,,drop=FALSE]
    integrand_a_sub <- S2_a[1:i,,drop=FALSE] * exp(-H3_a_sub)
    integrand_b_sub <- S2_b[1:i,,drop=FALSE] * exp(-H3_b_sub)
    value$p_nonterm_only[i+1,] <- 0.5 *
      colSums( (integrand_a_sub + integrand_b_sub) * (S1_a_sub - S1_b_sub) )

    integrand_a_sub <- S2_a[1:i,,drop=FALSE] * (1-exp(-H3_a_sub))
    integrand_b_sub <- S2_b[1:i,,drop=FALSE] * (1-exp(-H3_b_sub))
    value$p_both[i+1,] <- 0.5 *
      colSums( (integrand_a_sub + integrand_b_sub) * (S1_a_sub - S1_b_sub) )

    if(frailty & marg_flag){
      H1_a_sub <- H1_a[1:i,,drop=FALSE]
      H1_b_sub <- H1_b[1:i,,drop=FALSE]
      H2_a_sub <- H2_a[1:i,,drop=FALSE]
      H2_b_sub <- H2_b[1:i,,drop=FALSE]
      value$p_nonterm_only_marg[i+1,] <- 0.5 * colSums(
        (1 + theta * (H1_a_sub + H2_a_sub + H3_a_sub))^(-1/theta) +
          (1 + theta * (H1_a_sub + H2_b_sub + H3_b_sub))^(-1/theta) -
          (1 + theta * (H1_b_sub + H2_a_sub + H3_a_sub))^(-1/theta) -
          (1 + theta * (H1_b_sub + H2_b_sub + H3_b_sub))^(-1/theta))
    }

    #RETURN TO THESE ONCE I FINISH SOLVING THE MATTER OF THE SEMI-MARKOV!!
    if(se_fit_flag){
      #now, to generate standard errors for these predictions

      basis <- get_basis(y = tseq_temp,knots_vec = knots_list[[3]],hazard = hazard,deriv = FALSE)
      dbasis <- get_basis(y = tseq_temp,knots_vec = knots_list[[3]],hazard = hazard,deriv = TRUE)
      quad_points_sub <- transform_quad_points(n_quad = n_quad,
                                               quad_method=quad_method, a = 0,b = tseq_temp)
      basis_quad <- get_basis(y=quad_points_sub, knots_vec=knots_list[[3]], hazard = hazard,deriv = FALSE)

      dH1a_sub <- dH1a[1:i,,drop=FALSE]
      dH1b_sub <- dH1b[1:i,,drop=FALSE]
      dH2a_sub <- dH2a[1:i,,drop=FALSE]
      dH2b_sub <- dH2b[1:i,,drop=FALSE]

      # browser()

      for(j in 1:NCOL(H_list[[1]])){
        #jacobian for H3 has to follow the same structure as H3 itself does above
        #this dH3a will be "updated" below absed on the model type
        dH3a <- get_jac(tseq=tseq_temp, xnew=x1new[j,,drop=FALSE], beta=beta_list[[3]],
                        eta=eta_list[[1]][j], #WILL THIS CAUSE A PROBLEM IF ONE X HAS NO ELEMENTS BUT ANOTHER HAS SEVERAL?
                        hazard=hazard, phi=para[nP0_start[3]:nP0_end[3]],
                        basis=basis, dbasis=dbasis,
                        basis_quad=basis_quad, quad_weights=quad_weights,
                        func_type="H", H = H3_temp[,j],
                        log_out = FALSE )
        if(tolower(model) %in% c("m","markov")){
          #last row of H3[1:i] is H3(t2i), so
          #compute matrix of H3(t2) - H3(t1) for every t1 between 0 and t1i
          #by "sweeping" out H3(t2i) from every row of H3(t1) matrix
          dH3a <- t( -t(dH3a) + dH3a[(i+1),])
        } else if(p3tv>0){
          #if semi-markov and some kind of t1-dependence, just include the corresponding columns here!
          dH3a <- cbind(dH3a, H3_temp[,j] * h3_tv_basis[1:(i+1),,drop=FALSE])
        }
        #we know that by definition, the last element of tseq_temp is supposed to be zero
        dH3a[(i+1),] <- 0

        dH3b <- dH3a[-1,,drop=FALSE]
        dH3a <- dH3a[-(i+1),,drop=FALSE]

        #for p_term_only (aka CIF_t)
        #apply(S1_a*S2_a + S1_b*S2_a - S1_a*S2_b - S1_b*S2_b, MARGIN = 2, cumsum)
        S2_a_sub <- S2_a[1:i,,drop=FALSE]
        S2_b_sub <- S2_b[1:i,,drop=FALSE]
        J <- 0.5 * cbind(
          dH1a_sub * (S1_a_sub[,j] * (S2_a_sub[,j] * exp(-H3_a_sub[,j]) +
                                      S2_b_sub[,j] * exp(-H3_b_sub[,j]))) -
          dH1b_sub * (S1_b_sub[,j] * (S2_a_sub[,j] * exp(-H3_a_sub[,j]) +
                                      S2_b_sub[,j] * exp(-H3_b_sub[,j]))),
          dH2a_sub * (S2_a_sub[,j] * exp(-H3_a_sub[,j]) * (S1_a_sub[,j] - S1_b_sub[,j])) +
          dH2b_sub * (S2_b_sub[,j] * exp(-H3_b_sub[,j]) * (S1_a_sub[,j] -S1_b_sub[,j])),
          dH3a * (S2_a_sub[,j] * exp(-H3_a_sub[,j]) * (S1_a_sub[,j] - S1_b_sub[,j])) +
          dH3b * (S2_b_sub[,j] * exp(-H3_b_sub[,j]) * (S1_a_sub[,j] - S1_b_sub[,j])))
        value$se_p_nonterm_only[i+1,j] <- sqrt(sum(J %*% Finv_sub123 %*% t(J)))

        J <- 0.5 * cbind(
          dH1a_sub * (S1_a_sub[,j] * (S2_a_sub[,j] * (1-exp(-H3_a_sub[,j])) +
                                      S2_b_sub[,j] * (1-exp(-H3_b_sub[,j])))) -
          dH1b_sub * (S1_b_sub[,j] * (S2_a_sub[,j] * (1-exp(-H3_a_sub[,j])) +
                                      S2_b_sub[,j] * (1-exp(-H3_b_sub[,j])))),
          dH2a_sub * (S2_a_sub[,j] * ((1-exp(-H3_a_sub[,j])) * S1_a_sub[,j] -
                                      (1-exp(-H3_a_sub[,j])) * S1_b_sub[,j])) +
          dH2b_sub * (S2_b_sub[,j] * ((1-exp(-H3_b_sub[,j])) * S1_a_sub[,j] -
                                      (1-exp(-H3_b_sub[,j])) * S1_b_sub[,j])),
          dH3a * (S2_a_sub[,j] * exp(-H3_a_sub[,j]) * (S1_b_sub[,j] - S1_a_sub[,j])) +
          dH3b * (S2_b_sub[,j] * exp(-H3_b_sub[,j]) * (S1_b_sub[,j] - S1_a_sub[,j])))
        value$se_p_both[i+1,j] <- sqrt(sum(J %*% Finv_sub123 %*% t(J)))

        #last piece, se for marginal probabilities of "both" and "nonterm only"
        if(frailty & marg_flag){
          #for p_nonterm (aka CIF_nt)
          #apply(S2_a*S1_a + S2_b*S1_a - S2_a*S1_b - S2_b*S1_b, MARGIN = 2, cumsum)
          J <- 0.5 * cbind(
            dH1a_sub * ((1 + theta * (H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1) +
                        (1 + theta * (H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1)) -
            dH1b_sub * ((1 + theta * (H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1) +
                        (1 + theta * (H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1)),
            dH2a_sub * ((1 + theta * (H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1) -
                        (1 + theta * (H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1)) +
            dH2b_sub * ((1 + theta * (H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1) -
                        (1 + theta * (H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1)),
            dH3a * ((1 + theta * (H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1) -
                    (1 + theta * (H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1)) +
            dH3b * ((1 + theta * (H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1) -
                    (1 + theta * (H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1)),
            -marg_logtheta_deriv(a = H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j], ltheta = log(theta)) +
              -marg_logtheta_deriv(a = H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j], ltheta = log(theta)) -
              -marg_logtheta_deriv(a = H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j], ltheta = log(theta)) -
              -marg_logtheta_deriv(a = H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j], ltheta = log(theta)))
          value$se_p_nonterm_only_marg[i+1,j] <- sqrt(sum(J %*% Finv_sub123_marg %*% t(J)))

          #for p_both
          J <- 0.5 * cbind(
            dH1a_sub * (((1 + theta * (H1_a_sub[,j] + H2_a_sub[,j]))^(-1/theta - 1) -
                         (1 + theta * (H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1)) +
                        ((1 + theta * (H1_a_sub[,j] + H2_b_sub[,j]))^(-1/theta - 1) -
                         (1 + theta * (H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1))) -
            dH1b_sub * (((1 + theta * (H1_b_sub[,j] + H2_a_sub[,j]))^(-1/theta - 1) -
                         (1 + theta * (H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1)) +
                        ((1 + theta * (H1_b_sub[,j] + H2_b_sub[,j]))^(-1/theta - 1) -
                         (1 + theta * (H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1))),
            dH2a_sub * (((1 + theta * (H1_a_sub[,j] + H2_a_sub[,j]))^(-1/theta - 1) -
                         (1 + theta * (H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1)) -
                        ((1 + theta * (H1_b_sub[,j] + H2_a_sub[,j]))^(-1/theta - 1) -
                         (1 + theta * (H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1))) +
            dH2b_sub * (((1 + theta * (H1_a_sub[,j] + H2_b_sub[,j]))^(-1/theta - 1) -
                         (1 + theta * (H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1)) -
                        ((1 + theta * (H1_b_sub[,j] + H2_b_sub[,j]))^(-1/theta - 1) -
                         (1 + theta * (H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1))),
            dH3a * ((1 + theta * (H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1) -
                    (1 + theta * (H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j]))^(-1/theta - 1)) +
            dH3b * ((1 + theta * (H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1) -
                    (1 + theta * (H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j]))^(-1/theta - 1)),
            (marg_logtheta_deriv(a = H1_a_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j], ltheta = log(theta)) -
               marg_logtheta_deriv(a = H1_a_sub[,j] + H2_a_sub[,j], ltheta = log(theta))) +
              (marg_logtheta_deriv(a = H1_a_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j], ltheta = log(theta)) -
                 marg_logtheta_deriv(a = H1_a_sub[,j] + H2_b_sub[,j], ltheta = log(theta))) -
              (marg_logtheta_deriv(a = H1_b_sub[,j] + H2_a_sub[,j] + H3_a_sub[,j], ltheta = log(theta)) -
                 marg_logtheta_deriv(a = H1_b_sub[,j] + H2_a_sub[,j], ltheta = log(theta))) -
              (marg_logtheta_deriv(a = H1_b_sub[,j] + H2_b_sub[,j] + H3_b_sub[,j], ltheta = log(theta)) -
                 marg_logtheta_deriv(a = H1_b_sub[,j] + H2_b_sub[,j], ltheta = log(theta))))
          value$se_p_both_marg[i+1,j] <- sqrt(sum(J %*% Finv_sub123_marg %*% t(J)))
        }
      }

    }

  }

  value$p_nonterm_only <- t( t(value$p_nonterm_only) / p_neither0)
  value$p_both <- t( t(value$p_both) / p_neither0)
  if(frailty & marg_flag){
    value$p_nonterm_only_marg <- t( t(value$p_nonterm_only_marg) / p_neither0_marg)

    #The marginal expression for "both" probability is literally
    #the difference of these two quantities...
    value$p_both_marg <- value$p_nonterm_marg - value$p_nonterm_only_marg

  }

  value
}


