#' Extract robust "sandwich" estimate of variance-covariance matrix
#'
#' This function computes the robust "sandwich" variance estimate from frequentist
#' model object, provided that the model object includes the "cheese" or "meat" matrix comprising
#' the summed outer products of each individual subject's scores evaluated at the MLE.
#'
#' @param object Model object
#' @param adjust Boolean for whether to include "finite sample"
#'  degrees-of-freedom adjustment by multiplication with \eqn{n/(n-k)}
#'  where \eqn{n} is the number of observations and
#'  \eqn{k} is the number of estimated parameters.
#' @param ... Other arguments (currently unused).
#'
#' @export
vcov_sandwich <- function (object, adjust = FALSE, ...){
  if(is.null(object$cheese) || is.null(object$Finv)){stop("fitted model does not have 'cheese' for sandwich variance.")}
  val <- object$Finv %*% object$cheese %*% object$Finv
  if(adjust) val <- val * object$nobs / (object$nobs - NROW(object$Finv))
  rownames(val) <- colnames(val) <- names(object$estimate)
  val
}

#' @export
vcov.Freq_HReg2 <- function (object, ...){
  val <- object$Finv
  rownames(val) <- colnames(val) <- names(object$estimate)
  val
}

#' @export
logLik.Freq_HReg2 <- function (object, ...){
  val <- object$logLike; class(val) <- "logLik"
  attr(x = val,which = "df") <- length(object$estimate)
  attr(x = val,which = "nobs") <- object$nobs
  val
}

#' @export
extractAIC.Freq_HReg2 <- function (fit, scale, k=2, ...){
  c(fit$nobs,stats::AIC(fit,k=k))
}

#' @export
coef.Freq_HReg2 <- function (object, ...){
  object$estimate
}

#' @export
print.Freq_HReg2 <- function (x, digits = 3, alpha = 0.05, ...)
{
  conf.level = alpha
  logEst <- x$estimate
  logSE <- if(all(is.na(x$Finv))) NA else sqrt(diag(x$Finv))
  value <- cbind(logEst, logSE,
                 logEst - abs(stats::qnorm(conf.level/2, 0, 1)) * logSE,
                 logEst + abs(stats::qnorm(conf.level/2, 0, 1)) * logSE)
  dimnames(value) <- list(x$myLabels, c("Estimate", "SE", "LL", "UL"))
  if (x$class[2] == "Surv") {
    cat("\nAnalysis of independent univariate time-to-event data \n")
    cat(x$class[4], "baseline hazard specification\n")
    cat("Confidence level: ", conf.level, "\n", sep = "")
    if (sum(x$nP) != 0) {
      cat("\nRegression coefficients:\n")
      print(round(value[-c(1:(sum(x$nP0))), ], digits = digits))
    }
  }
  if (x$class[2] == "ID") {
    cat("\nAnalysis of independent semi-competing risks data \n")
    cat(x$class[4], "baseline hazard specification\n")
    cat(x$class[5], "specification for h3\n")
    cat("Confidence level: ", conf.level, "\n", sep = "")
    cat("\nVariance of frailties, theta:\n")
    if (x$frailty == TRUE)
      value_theta <- matrix(exp(value[sum(x$nP0)+1, ]), ncol = 4)
      dimnames(value_theta) <- list("", c("Estimate", "SE", "LL", "UL"))
      #update theta SE with delta method
      value_theta[1, 2] <- value[sum(x$nP0)+1, 2] * exp(value[sum(x$nP0)+1, 1])

      print(round(value_theta, digits = digits))
    if (x$frailty == FALSE)
      cat("NA")
    if (sum(x$nP) != 0) {
      cat("\nRegression coefficients:\n")
      if (x$frailty == TRUE)
        print(round(value[-c(1:(sum(x$nP0)+1)), ], digits = digits))
      if (x$frailty == FALSE)
        print(round(value[-c(1:(sum(x$nP0))), ], digits = digits))
    }
    cat("\nNote: Covariates are arranged in order of transition number, 1->3.\n")
  }
  invisible()
}

#' @export
summary.Freq_HReg2 <- function (object, digits = 3, alpha = 0.05, ...) {
  # browser()
  conf.level = alpha
  obj <- object
  logEst <- obj$estimate
  logSE <- if(all(is.na(obj$Finv))) NA else sqrt(diag(obj$Finv))
  results <- cbind(logEst,
                   logEst - abs(stats::qnorm(conf.level/2, 0, 1)) * logSE,
                   logEst + abs(stats::qnorm(conf.level/2, 0, 1)) * logSE)
  if (obj$class[2] == "Surv") {
    output.coef <- results[-c(1:obj$nP0),,drop=FALSE]
    dimnames(output.coef) <- list(unique(obj$myLabels[-c(1:obj$nP0)]),
                                  c("beta", "LL", "UL"))
    output.h0 <- results[1:obj$nP0,,drop=FALSE]
    if(obj$class[4]=="Weibull"){
      dimnames(output.h0) <- list(c("Weibull: log-kappa", "Weibull: log-alpha"),
                                  c("h-PM", "LL", "UL"))
      knots_mat <- NULL
    } else{
      dimnames(output.h0) <- list(paste0(obj$class[4],": phi",1:obj$nP0),
                               c("h-PM", "LL", "UL"))
      knots_mat <- as.matrix(obj$knots_vec)
      rownames(knots_mat) <- paste0("knot",1:length(obj$knots_vec))
    }
    value <- list(coef = output.coef, h0 = output.h0,
                  logLike = obj$logLike,
                  nP = nrow(results), class = obj$class,
                  conf.level = conf.level, knots_mat=knots_mat)
  }
  if (obj$class[2] == "ID") {
    nP.0 <- ifelse(obj$frailty, sum(obj$nP0)+1, sum(obj$nP0))
    nP.1 <- obj$nP[1]; nP.2 <- obj$nP[2]; nP.3 <- obj$nP[3]
    beta.names <- unique(obj$myLabels[-c(1:nP.0)])
    nP <- length(beta.names)
    output <- matrix(NA, nrow = nP, ncol = 9)
    dimnames(output) <- list(beta.names, c("beta1", "LL", "UL",
                                           "beta2", "LL", "UL",
                                           "beta3", "LL", "UL"))
    for (i in 1:nP) {
      if (nP.1 != 0) {
        for (j in 1:nP.1) if (obj$myLabels[nP.0+j] == beta.names[i])
          output[i,1:3] <- results[nP.0+j,]
      }
      if (nP.2 != 0) {
        for (j in 1:nP.2) if (obj$myLabels[nP.0+nP.1+j] == beta.names[i])
          output[i,4:6] <- results[nP.0+nP.1+j,]
      }
      if (nP.3 != 0) {
        for (j in 1:nP.3) if (obj$myLabels[nP.0+nP.1+nP.2+j] == beta.names[i])
          output[i,7:9] <- results[nP.0+nP.1+nP.2+j,]
      }
    }
    output.coef <- output
    output <- matrix(NA, nrow = 1, ncol = 3)
    dimnames(output) <- list(c("theta"), c("Estimate", "LL", "UL"))
    if (obj$frailty == TRUE)
      output[1, ] <- exp(results[nP.0, ])
    if (obj$frailty == FALSE)
      output[1, ] <- rep(NA, 3)
    output.theta <- output

    knots_mat <- NULL
    if(obj$class[4]=="Weibull"){
      output <- matrix(NA, nrow = 2, ncol = 9)
      dimnames(output) <- list(c("Weibull: log-kappa", "Weibull: log-alpha"),
                               c("h1-PM", "LL", "UL",
                                 "h2-PM", "LL", "UL",
                                 "h3-PM", "LL", "UL"))
      output[1, 1:3] <- results[1, ]
      output[1, 4:6] <- results[3, ]
      output[1, 7:9] <- results[5, ]
      output[2, 1:3] <- results[2, ]
      output[2, 4:6] <- results[4, ]
      output[2, 7:9] <- results[6, ]
    } else{ #this covers piecewise and spline
      p01 <- obj$nP0[1]; p02 <- obj$nP0[2]; p03 <- obj$nP0[3]
      p0max <- max(obj$nP0)

      #generate "wide" matrix of baseline parameters by padding with 0s so all are same height
      output <- cbind(
        rbind(results[1:p01,],matrix(data=0,ncol=3,nrow=(p0max-p01))),
        rbind(results[(1+p01):(p01+p02),],matrix(data=0,ncol=3,nrow=(p0max-p02))),
        rbind(results[(1+p01+p02):(p01+p02+p03),],matrix(data=0,ncol=3,nrow=(p0max-p03)))
      )
      dimnames(output) <- list(paste0(obj$class[4],": phi",1:p0max),
                               c("h1-PM", "LL", "UL",
                                 "h2-PM", "LL", "UL",
                                 "h3-PM","LL", "UL"))

      #lastly, make a matrix with the knot locations, padded with NAs
      knotmax <- max(sapply(obj$knots_list,length))
      knots_mat <- sapply(obj$knots_list,FUN = function(x) c(x,rep(NA,knotmax-length(x))))
      dimnames(knots_mat) <- list(paste0("knot",1:knotmax), c("h1","h2","h3"))
    }

    output.h0 <- output
    value <- list(coef = output.coef, theta = output.theta,
                  h0 = output.h0, logLike = obj$logLike,
                  nP = nrow(results), class = obj$class,
                  conf.level = conf.level, knots_mat=knots_mat)
  }
  class(value) <- "summ.Freq_HReg2"
  return(value)
}

#' @export
print.summ.Freq_HReg2 <- function (x, digits = 3, ...)
{
  obj <- x
  if (obj$class[2] == "Surv") {
    cat("\nAnalysis of independent univariate time-to-event data \n")
  }
  if (obj$class[2] == "ID") {
    cat("\nAnalysis of independent semi-competing risks data \n")
    cat(obj$class[4], "baseline hazard specification\n")
    cat(obj$class[5], "specification for h3\n")
  }
  cat("Confidence level: ", x$conf.level, "\n", sep = "")
  if (!is.null(obj$coef)) {
    cat("\nHazard ratios:\n")
    print(round(exp(obj$coef), digits = digits))
  }
  if (obj$class[2] == "ID") {
    cat("\nVariance of frailties:\n")
    print(round(obj$theta, digits = digits))
  }
  cat("\nBaseline hazard function components:\n")
  print(round(obj$h0, digits = digits))
  if (obj$class[4] != "Weibull") {
    cat("\nKnots:\n")
    print(round(obj$knots_mat, digits = digits))
  }
  invisible()
}








pred_helper <- function(tseq, xnew, func_type, conf.level,
                        phi, beta, logtheta, pred_type,
                        hess_sub, hazard_lab, knots_vec,
                        n_quad, quad_method){
  # browser()
  eta <- if(!is.null(xnew) & length(beta) > 0) as.vector(xnew %*% beta) else 0

  #derivations for marginal curves come mostly from Balan & Putter (2020)
  if(pred_type=="marginal"){ theta <- exp(logtheta); invtheta <- exp(-logtheta)}
  #marginal survival is S(t) = (1 + theta * Haz)^(-invtheta)
  #marginal hazard takes the form: h(t) = h(t) * E[gamma | T>t]
  #for gamma dist, marginal hazard is: haz / (1 + theta * Haz)

  #Originally I had considered delta method computation of CIs for marginal curves,
  #but the results I was getting seemed anticonservative so I've set those aside.

  if(hazard_lab=="Weibull"){
    kappa <- exp(phi[1])
    logalpha <- phi[2]
    alpha <- exp(logalpha)
    if(func_type=="S"){
      if(pred_type=="marginal"){
        H <- exp(eta)*kappa*(tseq)^alpha
        denom <- (1+theta*H)
        out_temp <- denom^(-invtheta)
        # #under weibull and gamma frailty, marginal
        # #log(S(t)) = -exp(-logtheta) * log(1 + exp(logtheta + logkappa + xtbeta) * t^(exp(logalpha)))
        # #Use log delta method to compute baseline survival confidence intervals
        # #matrix with as many columns as parameters in S1, and as many rows as times in tseq
        # J <- cbind(-H/denom, #partial derivative of logkappa
        #            -H/denom*alpha*log(tseq), #partial derivative of logalpha
        #            invtheta * log1p(theta*exp(eta)*kappa*tseq^alpha) - H/denom, #partial derivative of logtheta
        #            if(!is.null(xnew) & length(beta) > 0) -H/denom * matrix(xnew, nrow=length(tseq),
        #                                                                    ncol=length(xnew), byrow=T))
      } else{
        out_temp <- exp(-kappa * (tseq)^alpha * exp(eta))
        #under weibull, log(-log(S(t))) = logkappa + xtbeta + exp(logalpha)*log(t)
        #Use log-log delta method to compute baseline survival confidence intervals
        #matrix with as many columns as parameters in S1, and as many rows as times in tseq
        #so, each row of J is the gradient at a particular tseq wrt logkappa, logalpha, beta1, ..., betap
        J <- cbind(1, #partial derivative of logkappa
                   exp(logalpha) * log(tseq), #partial derivative of logalpha
                   if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow=length(tseq),
                                                                ncol=length(xnew), byrow=T))
      }
    } else if(func_type=="h"){
      #Weibull hazard is exp(logalpha) * exp(logkappa) * t^(exp(logalpha)-1) * exp(xtbeta)
      out_temp <- alpha * kappa * (tseq)^(alpha - 1) * exp(eta)
      if(pred_type=="marginal"){
        H <- exp(eta)*kappa*(tseq)^alpha
        denom <- (1+theta*H)
        out_temp <- out_temp / denom

        #here I tested out a few equivalent ways of computing marginal hazard CIs,
        #and they gave the same (seemingly anticonservative) results..

        # #under weibull and gamma frailty, marginal
        # #log(haz) = logalpha + logkappa + (exp(logalpha)-1)*log(t) + xtbeta
        # #           - log(1 + exp(logtheta + logkappa + xtbeta) * t^(exp(logalpha))
        # #Use delta method on log(haz) to compute baseline hazard confidence intervals
        # J <- cbind(1 - theta*H/denom, #partial derivative of logkappa
        #            1 + alpha*log(tseq) - theta*H*alpha*log(tseq)/denom, #partial derivative of logalpha
        #            -theta*H/denom, #partial derivative of logtheta
        #            if(!is.null(xnew) & length(beta) > 0) (1 - theta*H/denom) * matrix(xnew, nrow=length(tseq),
        #                                                                               ncol=length(xnew), byrow=T))

        # #equivalently, compute jacobian of log marginal hazard numerically
        # marg_loghaz <- function(para,t){
        #   #1 is alpha, 2 is kappa, 3 is logtheta assume no covariates for now
        #   para[1] + para[2]  + (exp(para[2])-1) * log(t) - log1p(exp(para[1]+para[3])*t^exp(para[2]))
        # }
        # out_temp <- exp(marg_loghaz(para = c(phi,logtheta), t=tseq))
        # J <- NULL
        # for(i in 1:length(tseq)){
        #   J <- rbind(J,pracma::grad(f = marg_loghaz,x0 = c(phi,logtheta),t=tseq[i]))
        # }

        # #this one also used a different delta method, in terms of hazard directly
        # marg_haz <- function(para,t){
        #   #1 is alpha, 2 is kappa, 3 is logtheta assume no covariates for now
        #   exp(para[1] + para[2]) * t^(exp(para[2])-1) / (1 + exp(para[1]+para[3])*t^exp(para[2]))
        # }
        # out_temp <- marg_haz(para = c(phi,logtheta), t=tseq)
        # J <- NULL
        # for(i in 1:length(tseq)){
        #   J <- rbind(J,pracma::grad(f = marg_haz,x0 = c(phi,logtheta),t=tseq[i]))
        # }
      } else{ #conditional hazard jacobian for computing SE
        #log(haz) = logalpha + logkappa + (exp(logalpha)-1)*log(t) + xtbeta
        #Use delta method on log(haz) to compute baseline hazard confidence intervals
        #each row of J is gradient wrt logkappa, logalpha, beta1, ..., betap
        J <- cbind(1, #partial derivative of logkappa
                   (1 + alpha * log(tseq)), #partial derivative of logalpha
                   if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow=length(tseq),
                                                                ncol=length(xnew), byrow=T))
      }
    } else stop("'func_type' must be 'Surv' or 'Haz'")
  } else if(hazard_lab=="Piecewise Constant"){
    if(func_type=="S"){
      basis <- get_basis(y = tseq,knots_vec = knots_vec,hazard = "piecewise")
      Lambda0 <- as.vector(basis %*% exp(phi))
      if(pred_type=="marginal"){
        denom <- (1+theta*Lambda0*exp(eta))
        out_temp <- denom^(-invtheta)
        # #under piecewise, marginal
        # #log(S(t)) = -exp(-logtheta)*log(1 + exp(logtheta+xtbeta)* (basis %*% exp(phi)))
        # #For phi partial derivs, first, multiply every row of basis by exp(phi),
        # #then, multiply every column of the resulting matrix by exp(eta)/denom
        # J <- cbind(-t(t(basis) * exp(phi)) * exp(eta) / denom, #partial deriv of phi
        #            invtheta*log1p(theta*Lambda0*exp(eta)) - Lambda0*exp(eta)/denom, #partial deriv of logtheta
        #            if (!is.null(xnew) & length(beta) > 0) -Lambda0*exp(eta)/denom * matrix(xnew, nrow = length(tseq),
        #                                                          ncol = length(xnew), byrow = T))
      } else{
        out_temp <- exp(-Lambda0*exp(eta))
        #under piecewise, log(-log(S(t))) = log( basis %*% exp(phi)) + xtbeta
        J <- cbind(t(t(basis) * exp(phi)) / Lambda0, #partial deriv of phi
                   if (!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                                 ncol = length(xnew), byrow = T))
      }
    } else if(func_type=="h"){
      dbasis <- get_basis(y = tseq,knots_vec = knots_vec,hazard = "piecewise",deriv = TRUE)
      out_temp <- as.vector(dbasis %*% exp(phi) * exp(eta))
      if(pred_type=="marginal"){

        #here I tested out a few equivalent ways of computing marginal hazard
        #they gave same haz but slightly different (all anticonservative) CIs

        basis <- get_basis(y = tseq,knots_vec = knots_vec,hazard = "piecewise")
        Lambda0 <- as.vector(basis %*% exp(phi))
        denom <- (1+theta*Lambda0*exp(eta))
        # #For phi partial derivs, first, multiply every row of basis by exp(phi),
        # #then, multiply every column of the resulting matrix by exp(eta)/denom
        # J <- cbind(t(t(dbasis) * exp(phi)) * exp(eta) / out_temp -
        #              t(t(basis) * exp(phi)) * exp(eta+logtheta) / denom, #partial derivative of phi
        #            -theta*Lambda0*exp(eta)/denom, #partial derivative of logtheta
        #            if (!is.null(xnew) & length(beta) > 0) (1-theta*Lambda0*exp(eta)/denom) * matrix(xnew, nrow = length(tseq),
        #                                                                                             ncol = length(xnew), byrow = T))
        out_temp <- out_temp / denom

        # #compute jacobian of log marginal hazard numerically (slightly diff CI)
        # marg_loghaz <- function(para,t){
        #   basis <- get_basis(y = t,knots_vec = knots_vec,hazard = "piecewise")
        #   dbasis <- get_basis(y = t,knots_vec = knots_vec,hazard = "piecewise",deriv=TRUE)
        #   log(as.vector(dbasis %*% exp(para[-1]))) - log1p(exp(para[1])*as.vector(basis %*% exp(para[-1])))
        # }
        # out_temp <- exp(marg_loghaz(para = c(logtheta,phi), t=tseq))
        # J <- NULL
        # for(i in 1:length(tseq)){
        #   J <- rbind(J,pracma::grad(f = marg_loghaz,x0 = c(logtheta,phi),t=tseq[i]))
        # }

        # #this one also used a different delta method, in terms of hazard directly
        # marg_haz <- function(para,t){
        #   basis <- get_basis(y = t,knots_vec = knots_vec,hazard = "piecewise")
        #   dbasis <- get_basis(y = t,knots_vec = knots_vec,hazard = "piecewise",deriv=TRUE)
        #   as.vector(dbasis %*% exp(para[-1])) / (1 + exp(para[1])*as.vector(basis %*% exp(para[-1])))
        # }
        # out_temp <- marg_haz(para = c(logtheta,phi), t=tseq)
        # J <- NULL
        # for(i in 1:length(tseq)){
        #   J <- rbind(J,pracma::grad(f = marg_haz,x0 = c(logtheta,phi),t=tseq[i]))
        # }
      } else{
        #Use delta method on log(haz) to compute baseline hazard confidence intervals
        J <- cbind(t(t(dbasis) * exp(phi)) * exp(eta) / out_temp,
                   if (!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                                 ncol = length(xnew), byrow = T))
        # #Use delta method on hazard to compute baseline hazard confidence intervals
        # J <- cbind(t(t(dbasis) * exp(phi)) * exp(eta),
        #            if (!is.null(xnew) & length(beta) > 0) out_temp * matrix(xnew, nrow = length(tseq),
        #                                                                     ncol = length(xnew), byrow = T))
      }
    } else stop("'func_type' must be 'Surv' or 'Haz'")
  } else if(hazard_lab=="Royston-Parmar"){
    basis <- get_basis(y = tseq,knots_vec = knots_vec,hazard = "royston-parmar")
    Lambda0 <- as.vector(exp(basis %*% phi))
    if(func_type=="S"){
      if(pred_type=="marginal"){
        denom <- (1+theta*Lambda0*exp(eta))
        out_temp <- denom^(-invtheta)
      } else{
        out_temp <- exp(-Lambda0*exp(eta))
        J <- cbind(basis,
                   if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                                ncol = length(xnew), byrow = T))
      }
    } else if(func_type=="h"){
      dbasis <- get_basis(y = tseq,knots_vec = knots_vec,hazard = "royston-parmar",deriv = TRUE)
      out_temp <- as.vector(dbasis %*% phi / tseq * exp(basis %*% phi + eta))
      if(pred_type=="marginal"){
        basis <- get_basis(y = tseq,knots_vec = knots_vec,hazard = "royston-parmar")
        Lambda0 <- as.vector(exp(basis %*% phi))
        denom <- (1+theta*Lambda0*exp(eta))
        out_temp <- out_temp / denom
      } else {
        #Use delta method on log(haz) to compute baseline hazard confidence intervals
        #log(h(t)) = log(s'(z)) + s(z) + xtbeta where z=log(t)
        J <- cbind(dbasis/as.vector(dbasis %*% phi) + basis,
                   if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                                ncol = length(xnew), byrow = T))
        # #Use delta method on hazard to compute baseline hazard confidence intervals
        # #h(t) = s'(z)/t * exp(s(z) + xtbeta) where z=log(t)
        # #so, each row of J is gradient wrt phi1, ..., phik, beta1, ..., betap
        # J <- cbind(dbasis*Lambda0*exp(eta)/tseq + out_temp*basis,
        #            if(!is.null(xnew) & length(beta) > 0) out_temp * matrix(xnew, nrow = length(tseq),
        #                                                          ncol = length(xnew), byrow = T))
      }
    } else stop("'func_type' must be 'Surv' or 'Haz'")
  } else if(hazard_lab=="B-Spline"){
    quad_weights <- get_quad_pointsweights(n_quad=n_quad,
                                           quad_method=quad_method)$weights
    quad_points <- transform_quad_points(n_quad = n_quad,
                                         quad_method=quad_method, a = 0,b = tseq)
    basis <- get_basis(y=tseq,knots_vec=knots_vec,hazard="bspline")
    basis_quad <- get_basis(y=quad_points, knots_vec=knots_vec,hazard="bspline")
    lambda0 <- as.vector(exp(basis_quad %*% phi))
    if(func_type=="S"){
      #reshape lambda0 from a n*n_quad length vector
      #to an n by n_quad matrix, then multiply with n_quad length weights
      #to get final Lambda0
      Lambda0 <- tseq/2 * as.vector(matrix(lambda0,ncol=n_quad,byrow = TRUE) %*% quad_weights)
      if(pred_type=="marginal"){
        denom <- (1+theta*Lambda0*exp(eta))
        out_temp <- denom^(-invtheta)
      } else{
        out_temp <- exp(-Lambda0*exp(eta))
        #matrix with as many columns as parameters in S1, and as many rows as times in tseq
        #under bspline, log(-log(S(t))) = log( [numerical integral of exp(Btphi)] ) + xtbeta
        #so, each row of J is the gradient at a particular tseq wrt phi1, ..., phiK, beta1, ..., betap
        #specifically, for gradient of phi we take each basis column and multiply it by lambda,
        #evaluate the numerical integral of the product, and then divide off Lambda0
        J_phi <- apply(X = basis_quad * lambda0, MARGIN = 2,
                       FUN = function(x) tseq/2 * as.vector(matrix(x,ncol=n_quad,byrow = TRUE) %*% quad_weights)) / Lambda0
        J <- cbind( J_phi,
                    if (!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                                  ncol = length(xnew), byrow = T))
      }
    } else if(func_type=="h"){
      out_temp <- as.vector(exp(basis %*% phi) * exp(eta))
      if(pred_type=="marginal"){
        Lambda0 <- tseq/2 * as.vector(matrix(lambda0,ncol=n_quad,byrow = TRUE) %*% quad_weights)
        denom <- (1+theta*Lambda0*exp(eta))
        out_temp <- out_temp / denom
      } else{
        #Use delta method on log(haz) to compute baseline hazard confidence intervals
        #log(haz) is Btphi + xtbeta
        J <- cbind(basis, #partial derivative of phi
                   if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                                ncol = length(xnew), byrow = T))
        # #Use delta method on hazard to compute baseline hazard confidence intervals
        # #hazard is exp(Btphi) * exp(xtbeta)
        # #so, each row of J is gradient wrt phi1, ..., phik, beta1, ..., betap
        # J <- cbind(out_temp * basis, #partial derivative of phi
        #            if(!is.null(xnew) & length(beta) > 0) out_temp * matrix(xnew, nrow = length(tseq),
        #                                                          ncol = length(xnew), byrow = T))
      }
    } else stop("'func_type' must be 'Surv' or 'Haz'")
  } else stop("'hazard_lab' must be 'Weibull', 'Piecewise Constant', 'Royston-Parmar' or 'B-Spline'")

  #don't compute CI for marginal hazard/survivor curves, only conditional ones
  if(pred_type=="marginal"){
    LL <- UL <- rep(NA,length(out_temp))
  } else{
    out_var <- if(all(is.na(hess_sub))) NA else J %*% hess_sub %*% t(J)
    out_se <- if(all(is.na(out_var))) NA else sqrt(diag(out_var))
    out_se[is.nan(out_se)] <- 0
    if(func_type=="S"){
      #formula based on se of log(-log(S(t)))
      LL <- out_temp^exp(-stats::qnorm(conf.level/2) * out_se)
      UL <- out_temp^exp(stats::qnorm(conf.level/2) * out_se)

      # #when I was experimenting with marginal CI, I used log(S(t)) for those
      # if(pred_type=="marginal"){
      #   LL <- out_temp * exp(stats::qnorm(conf.level/2) * out_se)
      #   UL <- out_temp * exp(-stats::qnorm(conf.level/2) * out_se)
      # }
    } else if(func_type=="h"){
      #computing CI from log-scale hazard se
      LL <- out_temp * exp(stats::qnorm(conf.level/2) * out_se)
      UL <- out_temp * exp(-stats::qnorm(conf.level/2) * out_se)
      # #old formula from direct calculation of se on hazard scale
      # LL <- out_temp + stats::qnorm(conf.level/2) * out_se #sign reversed because 0.025 quantile is negative
      # UL <- out_temp - stats::qnorm(conf.level/2) * out_se
      # LL[LL < 0] <- 0
    } else stop("'func_type' must be 'Surv' or 'Haz'")
  }

  if (tseq[1] == 0) {
    tseq <- tseq[-1]
    out_temp <- out_temp[-1]
    LL <- LL[-1]
    UL <- UL[-1]
  }

  return(data.frame(time=tseq,out=out_temp,LL=LL,UL=UL))
}

#' @export
predict.Freq_HReg2 <- function (object, xnew = NULL,
                                x1new = NULL, x2new = NULL, x3new = NULL,
                                tseq = seq(0,object$ymax,length.out = 100),
                                alpha = 0.05, pred_type="conditional",
                                n_quad=15, quad_method="kronrod", ...) {
  # browser()
  conf.level = alpha
  obj <- object
  yLim <- NULL
  nP = obj$nP
  nP0 = obj$nP0
  if (obj$class[2] == "Surv") {
    if (!(is.null(x1new) & is.null(x2new) & is.null(x3new))) {
      stop("'x1new','x2new', and 'x3new' are for semi-competing risks models and must be specified as NULL for univariate models")
    }
    value <- list()
    for(type in c("h","S")){
      value[[type]] <- pred_helper(tseq = tseq,xnew = xnew,
                            func_type = type,conf.level = conf.level,
                            phi = obj$estimate[1:nP0],
                            beta = obj$estimate[-(1:nP0)],
                            logtheta=0, pred_type = "conditional", #by default
                            hess_sub = obj$Finv[c(1:nP0, if(!is.null(xnew)) (nP0+1):(nP0+nP)),
                                                c(1:nP0, if(!is.null(xnew)) (nP0+1):(nP0+nP))],
                            knots_vec = obj$knots_vec,
                            hazard_lab = obj$class[4],
                            n_quad=n_quad,quad_method=quad_method)
      colnames(value[[type]]) <- c("time",type,"LL","UL")
    }
  } else if(obj$class[2]=="ID"){
    #if there is no frailty, silently switch prediction type to "conditional"
    if(pred_type=="marginal" & !obj$frailty){ pred_type <- "conditional" }
    if (!is.null(xnew)) {
      stop("'xnew' is for univariate models so it must be specified as NULL for semi-competing risks models")
    }
    nP0_tot <- if (obj$frailty == TRUE) sum(nP0) + 1 else sum(nP0)
    nP0_start <- 1 + c(0,nP0[1],nP0[1]+nP0[2])
    nP0_end <- c(nP0[1],nP0[1]+nP0[2],nP0[1]+nP0[2]+nP0[3])
    nP_start <- 1 + c(nP0_tot,nP0_tot+nP[1],nP0_tot+nP[1]+nP[2])
    nP_end <- c(nP0_tot+nP[1],nP0_tot+nP[1]+nP[2],nP0_tot+nP[1]+nP[2]+nP[3])
    x_list <- list(x1new,x2new,x3new)
    value <- list()
    for(type in c("h","S")){
      for(i in 1:3){
        value[[paste0(type,".",i)]] <-
          pred_helper(tseq = tseq,xnew = x_list[[i]],
            func_type = type,conf.level = conf.level,
            phi = obj$estimate[nP0_start[i]:nP0_end[i]],
            beta = obj$estimate[nP_start[i]:nP_end[i]],
            logtheta=if(obj$frailty) obj$estimate[nP0_tot] else NA,
            pred_type = pred_type,
            hess_sub = obj$Finv[c(nP0_start[i]:nP0_end[i],
                                  if(pred_type=="marginal") nP0_tot,
                                  if(!is.null(x_list[[i]])) nP_start[i]:nP_end[i]),
                                c(nP0_start[i]:nP0_end[i],
                                  if(pred_type=="marginal") nP0_tot,
                                  if(!is.null(x_list[[i]])) nP_start[i]:nP_end[i])],
            hazard_lab = obj$class[4], knots_vec=obj$knots_list[[i]],
            n_quad=n_quad,quad_method=quad_method)
        colnames(value[[paste0(type,".",i)]]) <- c("time",paste0(type,".",i),
                                                   paste0("LL.",i),
                                                   paste0("UL.",i))
      }
    }
  }
  value$xnew <- xnew
  value$x1new <- x1new; value$x2new <- x2new; value$x3new <- x3new
  value$tseq <- tseq
  value$class <- obj$class
  class(value) <- "pred.Freq_HReg2"
  return(value)
}


#' @export
plot.pred.Freq_HReg2 <- function (x, plot.est = "Haz", xlab = NULL, ylab = NULL, ...)
{
  # browser()
  obj <- x
  T2seq <- x$tseq
  yLim <- NULL

  if (is.null(ylab)){
    ylab <- switch(plot.est, "Surv"="Survival","Haz"="Hazard")
  }

  if (obj$class[2] == "Surv") {
    if (is.null(yLim)) {
      if (plot.est == "Surv") {
        yLim <- seq(from = 0, to = 1, by = 0.2)
      }
      if (plot.est == "Haz") {
        grid <- (max(obj$h$UL) - min(obj$h$LL))/5
        yLim <- seq(from = min(obj$h$LL), to = max(obj$h$UL),
                    by = grid)
      }
    }
    if (is.null(xlab)){ xlab <- "Time" }
    if (plot.est == "Surv") {
      plot(range(T2seq), range(yLim), xlab = xlab, ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               S(t), "")), axes = FALSE)
      graphics::axis(1, at = T2seq)
      graphics::axis(2, at = yLim)
      graphics::lines(obj$S$time, obj$S$S, col = "red", lwd = 3)
      graphics::lines(obj$S$time, obj$S$LL, col = "red", lwd = 3, lty = 3)
      graphics::lines(obj$S$time, obj$S$UL, col = "red", lwd = 3, lty = 3)
    }
    if (plot.est == "Haz") {
      plot(range(T2seq), range(yLim), xlab = xlab, ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               h(t), "")), axes = FALSE)
      graphics::axis(1, at = T2seq)
      graphics::axis(2, at = round(yLim, 4))
      graphics::lines(obj$h$time, obj$h$h, col = "red", lwd = 3)
      graphics::lines(obj$h$time, obj$h$LL, col = "red", lwd = 3, lty = 3)
      graphics::lines(obj$h$time, obj$h$UL, col = "red", lwd = 3, lty = 3)
    }
  } else if (obj$class[2] == "ID") {
    if (is.null(xlab)) {
      xlab <- c("Time", "Time", "Time")
      if (obj$class[5] == "semi-Markov") {
        xlab[3] <- "Time since non-terminal event"
      }
    }
    if (is.null(yLim)) {
      if (plot.est == "Surv") {
        yLim <- seq(from = 0, to = 1, by = 0.2)
      }
      if (plot.est == "Haz") {
        ygrid <- (max(x$h.1$h.1,x$h.2$h.2,x$h.3$h.3,
                      x$h.1$UL.1, x$h.2$UL.2, x$h.3$UL.3,na.rm = TRUE) -
                    0)/5
        yLim <- seq(from = 0, to = max(x$h.1$h.1,x$h.2$h.2,x$h.3$h.3,
                                       x$h.1$UL.1, x$h.2$UL.2, x$h.3$UL.3,na.rm = TRUE), by = ygrid)
      }
    }
    if (plot.est == "Surv") {
      graphics::par(mfrow = c(1, 3))
      plot(range(T2seq), range(yLim), xlab = xlab[1], ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               S[1](t), "")), axes = FALSE)
      graphics::axis(1, at = T2seq)
      graphics::axis(2, at = yLim)
      graphics::lines(obj$S.1$time, obj$S.1$S.1, col = "blue", lwd = 3)
      graphics::lines(obj$S.1$time, obj$S.1$LL.1, col = "blue", lwd = 3,
                      lty = 3)
      graphics::lines(obj$S.1$time, obj$S.1$UL.1, col = "blue", lwd = 3,
                      lty = 3)
      plot(range(T2seq), range(yLim), xlab = xlab[2], ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               S[2](t), "")), axes = FALSE)
      graphics::axis(1, at = T2seq)
      graphics::axis(2, at = yLim)
      graphics::lines(obj$S.2$time, obj$S.2$S.2, col = "red", lwd = 3)
      graphics::lines(obj$S.2$time, obj$S.2$LL.2, col = "red", lwd = 3,
                      lty = 3)
      graphics::lines(obj$S.2$time, obj$S.2$UL.2, col = "red", lwd = 3,
                      lty = 3)
      plot(range(T2seq), range(yLim), xlab = xlab[3], ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               S[3](t), "")), axes = FALSE)
      graphics::axis(1, at = T2seq)
      graphics::axis(2, at = yLim)
      graphics::lines(obj$S.3$time, obj$S.3$S.3, col = "red", lwd = 3)
      graphics::lines(obj$S.3$time, obj$S.3$LL.3, col = "red", lwd = 3,
                      lty = 3)
      graphics::lines(obj$S.3$time, obj$S.3$UL.3, col = "red", lwd = 3,
                      lty = 3)
      graphics::par(mfrow = c(1, 1))
    }
    if (plot.est == "Haz") {
      graphics::par(mfrow = c(1, 3))
      plot(range(T2seq), range(yLim), xlab = xlab[1], ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               h[1](t), "")), axes = FALSE)
      graphics::axis(1, at = T2seq)
      graphics::axis(2, at = round(yLim, 4))
      graphics::lines(obj$h.1$time, obj$h.1$h.1, col = "blue", lwd = 3)
      graphics::lines(obj$h.1$time, obj$h.1$LL.1, col = "blue", lwd = 3,
                      lty = 3)
      graphics::lines(obj$h.1$time, obj$h.1$UL.1, col = "blue", lwd = 3,
                      lty = 3)
      plot(range(T2seq), range(yLim), xlab = xlab[2], ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               h[2](t), "")), axes = FALSE)
      graphics::axis(1, at = T2seq)
      graphics::axis(2, at = round(yLim, 4))
      graphics::lines(obj$h.2$time, obj$h.2$h.2, col = "red", lwd = 3)
      graphics::lines(obj$h.2$time, obj$h.2$LL.2, col = "red", lwd = 3,
                      lty = 3)
      graphics::lines(obj$h.2$time, obj$h.2$UL.2, col = "red", lwd = 3,
                      lty = 3)
      plot(range(T2seq), range(yLim), xlab = xlab[3], ylab = ylab,
           type = "n", main = expression(paste("Estimated ",
                                               h[3](t), "")), axes = FALSE)
      graphics::axis(1, at = T2seq)
      graphics::axis(2, at = round(yLim, 4))
      graphics::lines(obj$h.3$time, obj$h.3$h.3, col = "red", lwd = 3)
      graphics::lines(obj$h.3$time, obj$h.3$LL.3, col = "red", lwd = 3,
                      lty = 3)
      graphics::lines(obj$h.3$time, obj$h.3$UL.3, col = "red", lwd = 3,
                      lty = 3)
      graphics::par(mfrow = c(1, 1))
    }
  }
  invisible()
}

