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

pred_helper <- function(tseq, xnew, func_type,
                        phi, beta, hazard_lab, knots_vec,
                        n_quad, quad_method, vcov_sub, conf.level){
  # browser()
  stopifnot(func_type %in% c("S","h","H"))
  get_se <- !is.null(vcov_sub)
  if (tseq[1] == 0) {
    tseq <- tseq[-1]
  }

  #helper function for univariate hazard/survivor function estimation
  eta <- if(!is.null(xnew) & length(beta) > 0) as.vector(xnew %*% beta) else 0
  if(hazard_lab=="Weibull"){
    kappa <- exp(phi[1])
    logalpha <- phi[2]
    alpha <- exp(logalpha)
    if(func_type=="h"){
      #Weibull hazard is exp(logalpha) * exp(logkappa) * t^(exp(logalpha)-1) * exp(xtbeta)
      out_temp <- alpha * kappa * (tseq)^(alpha - 1) * exp(eta)
      if(get_se){
        #conditional hazard jacobian for computing SE
        #log(haz) = logalpha + logkappa + (exp(logalpha)-1)*log(t) + xtbeta
        #Use delta method on log(haz) to compute baseline hazard confidence intervals
        #each row of J is gradient wrt logkappa, logalpha, beta1, ..., betap
        J <- cbind(1, #partial derivative of logkappa
                   (1 + alpha * log(tseq)), #partial derivative of logalpha
                   if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow=length(tseq),
                                                                ncol=length(xnew), byrow=T))
      }
    } else {
      #out_temp is Lambda, may be transformed to survival at the end
      out_temp <- kappa * (tseq)^alpha * exp(eta)
      if(get_se){
        #under weibull, log(-log(S(t))) = logkappa + xtbeta + exp(logalpha)*log(t)
        #Use log-log delta method to compute baseline survival confidence intervals
        #matrix with as many columns as parameters in S1, and as many rows as times in tseq
        #so, each row of J is the gradient at a particular tseq wrt logkappa, logalpha, beta1, ..., betap
        J <- cbind(1, #partial derivative of logkappa
                   exp(logalpha) * log(tseq), #partial derivative of logalpha
                   if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow=length(tseq),
                                                                ncol=length(xnew), byrow=T))
      }
    }
  } else if(hazard_lab=="Piecewise Constant"){
    if(func_type=="h"){
      dbasis <- get_basis(y = tseq,knots_vec = knots_vec,hazard = "piecewise",deriv = TRUE)
      out_temp <- as.vector(dbasis %*% exp(phi) * exp(eta))
      if(get_se){
        #Use delta method on log(haz) to compute baseline hazard confidence intervals
        J <- cbind(t(t(dbasis) * exp(phi)) * exp(eta) / out_temp,
                   if (!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                                 ncol = length(xnew), byrow = T))
      }
    } else {
      basis <- get_basis(y = tseq,knots_vec = knots_vec,hazard = "piecewise")
      Lambda0 <- as.vector(basis %*% exp(phi))
      out_temp <- Lambda0*exp(eta)
      if(get_se){
        #under piecewise, log(-log(S(t))) = log( basis %*% exp(phi)) + xtbeta
        J <- cbind(t(t(basis) * exp(phi)) / Lambda0, #partial deriv of phi
                   if (!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                                 ncol = length(xnew), byrow = T))
      }
    }
  } else if(hazard_lab=="Royston-Parmar"){
    basis <- get_basis(y = tseq,knots_vec = knots_vec,hazard = "royston-parmar")
    Lambda0 <- as.vector(exp(basis %*% phi))
    if(func_type=="h"){
      dbasis <- get_basis(y = tseq,knots_vec = knots_vec,hazard = "royston-parmar",deriv = TRUE)
      out_temp <- as.vector(dbasis %*% phi / tseq * exp(basis %*% phi + eta))
      if(get_se){
        #Use delta method on log(haz) to compute baseline hazard confidence intervals
        #log(h(t)) = log(s'(z)) + s(z) + xtbeta where z=log(t)
        J <- cbind(dbasis/as.vector(dbasis %*% phi) + basis,
                   if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                                ncol = length(xnew), byrow = T))
      }
    } else {
      out_temp <- Lambda0*exp(eta)
      if(get_se){
        J <- cbind(basis,
                   if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                                ncol = length(xnew), byrow = T))
      }
    }
  } else if(hazard_lab=="B-Spline"){
    quad_weights <- get_quad_pointsweights(n_quad=n_quad,
                                           quad_method=quad_method)$weights
    quad_points <- transform_quad_points(n_quad = n_quad,
                                         quad_method=quad_method, a = 0,b = tseq)
    basis <- get_basis(y=tseq,knots_vec=knots_vec,hazard="bspline")
    basis_quad <- get_basis(y=quad_points, knots_vec=knots_vec,hazard="bspline")
    lambda0 <- as.vector(exp(basis_quad %*% phi))
    if(func_type=="h"){
      out_temp <- as.vector(exp(basis %*% phi) * exp(eta))
      if(get_se){
        #Use delta method on log(haz) to compute baseline hazard confidence intervals
        #log(haz) is Btphi + xtbeta
        J <- cbind(basis, #partial derivative of phi
                   if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                                ncol = length(xnew), byrow = T))
      }
    } else {
      #reshape lambda0 from a n*n_quad length vector
      #to an n by n_quad matrix, then multiply with n_quad length weights
      #to get final Lambda0
      Lambda0 <- tseq/2 * as.vector(matrix(lambda0,ncol=n_quad,byrow = TRUE) %*% quad_weights)
      out_temp <- Lambda0*exp(eta)
      if(get_se){
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
    }
  } else stop("'hazard_lab' must be 'Weibull', 'Piecewise Constant', 'Royston-Parmar' or 'B-Spline'")

  if(get_se){
    out_var <- if(all(is.na(vcov_sub))) NA else J %*% vcov_sub %*% t(J)
    out_se <- if(all(is.na(out_var))) NA else sqrt(diag(out_var))
    out_se[is.nan(out_se)] <- 0
    #computing CI from log-scale se
    LL <- out_temp * exp(stats::qnorm(conf.level/2) * out_se)
    UL <- out_temp * exp(-stats::qnorm(conf.level/2) * out_se)
  } else{
    LL <- UL <- rep(NA,length(tseq))
  }

  #transform survivor function
  if(func_type == "S"){
    out_temp <- exp(-out_temp)
    LL <- exp(-LL)
    UL <- exp(-UL)
  }

  return(as.data.frame(cbind(time=tseq,out=out_temp,
                             LL=if(get_se) LL,
                             UL=if(get_se) UL)))
}

#' @export
predict.Freq_HReg2 <- function (object, xnew = NULL,
                                x1new = NULL, x2new = NULL, x3new = NULL,
                                tseq = seq(0,object$ymax,length.out = 100),
                                alpha = 0.05, pred_type="conditional",
                                n_quad=15, quad_method="kronrod", ...) {
  browser()
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
                            vcov_sub = obj$Finv[c(1:nP0, if(!is.null(xnew)) (nP0+1):(nP0+nP)),
                                                c(1:nP0, if(!is.null(xnew)) (nP0+1):(nP0+nP))],
                            knots_vec = obj$knots_vec,
                            hazard_lab = obj$class[4],
                            n_quad=n_quad,quad_method=quad_method)
      colnames(value[[type]]) <- c("time",type,"LL","UL")
    }
  } else if(obj$class[2]=="ID"){
    if (!is.null(xnew)) {
      stop("'xnew' is for univariate models so it must be specified as NULL for semi-competing risks models")
    }
    nP0_tot <- if(obj$frailty) sum(nP0) + 1 else sum(nP0)
    nP0_start <- 1 + c(0,nP0[1],nP0[1]+nP0[2])
    nP0_end <- c(nP0[1],nP0[1]+nP0[2],nP0[1]+nP0[2]+nP0[3])
    nP_start <- 1 + c(nP0_tot,nP0_tot+nP[1],nP0_tot+nP[1]+nP[2])
    nP_end <- c(nP0_tot+nP[1],nP0_tot+nP[1]+nP[2],nP0_tot+nP[1]+nP[2]+nP[3])
    x_list <- list(x1new,x2new,x3new)
    value <- list()
    for(i in 1:3){
      for(type in c("h","H")){
        value[[paste0(type,".",i)]] <-
          pred_helper(tseq = tseq,xnew = x_list[[i]],
            func_type = type,conf.level = conf.level,
            phi = obj$estimate[nP0_start[i]:nP0_end[i]],
            beta = obj$estimate[nP_start[i]:nP_end[i]],
            vcov_sub = obj$Finv[c(nP0_start[i]:nP0_end[i],
                                  if(!is.null(x_list[[i]])) nP_start[i]:nP_end[i]),
                                c(nP0_start[i]:nP0_end[i],
                                  if(!is.null(x_list[[i]])) nP_start[i]:nP_end[i])],
            hazard_lab = obj$class[4], knots_vec=obj$knots_list[[i]],
            n_quad=n_quad,quad_method=quad_method)
        colnames(value[[paste0(type,".",i)]]) <- c("time",paste0(type,".",i),
                                                   paste0("LL.",i),
                                                   paste0("UL.",i))
      }
      #finally, get the survivor curves from the cumulative hazard curves, and relabel
      value[[paste0("S.",i)]] <- as.data.frame(cbind(time=value[[paste0("H.",i)]][,"time"],
                                       out=exp(-value[[paste0("H.",i)]][,paste0("H.",i)]),
                                       LL=exp(-value[[paste0("H.",i)]][,paste0("UL.",i)]),
                                       UL=exp(-value[[paste0("H.",i)]][,paste0("LL.",i)])))
      colnames(value[[paste0("S.",i)]]) <- c("time",paste0("S.",i),
                                             paste0("LL.",i),
                                             paste0("UL.",i))
    }

    #Following Xu (2010), marginal hazards reflect interplay of H1 and H2
    if(obj$frailty){
      theta <- exp(obj$estimate[nP0_tot])
      value$h.1$h.1.marg <- value$h.1$h.1 / (1 + theta*(value$H.1$H.1 + value$H.2$H.2))
      value$h.2$h.2.marg <- value$h.2$h.2 / (1 + theta*(value$H.1$H.1 + value$H.2$H.2))

      #Marginal hazard of h3 also depends on t1, so we calculate it for every
      #combination.
      if(obj$class[5] == "Markov"){
        t_len_new <- length(value$h.1$h.1) #length without 0
        index_mat <- expand.grid(1:length(value$h.3$h.3),
                                 1:length(value$h.1$h.1))
        index_mat <- index_mat[index_mat[,1]>=index_mat[,2],]
        value$h.3.marg <- as.data.frame(cbind(
          time1=value$h.1$time[index_mat[,2]],
          time2=value$h.3$time[index_mat[,1]],
          h.3.marg=value$h.3$h.3[index_mat[,1]] * (1+theta) / #note extra (1+theta) term
            (1 + theta * (value$H.1$H.1[index_mat[,2]] +
                          value$H.2$H.2[index_mat[,2]] +
                          value$H.3$H.3[index_mat[,1]] -
                          value$H.3$H.3[index_mat[,2]]))
          ))

        #add on a "time zero" set showing the marginal h3
        #assuming non-terminal event occurred "immediately"
        value$h.3.marg <- rbind(
          cbind(time1=0,
                time2=value$h.3$time,
                h.3.marg=value$h.3$h.3 * (1+theta) / (1+theta*value$H.3$H.3)),
          value$h.3.marg)

        #for convenience, add functions where the non-terminal event
        #perpetually "just" happened to the h.3 data frame, and where the non-terminal
        #event happened at time 0
        value$h.3$h.3.marg0 <- value$h.3.marg$h.3.marg[value$h.3.marg$time1==value$h.3.marg$time2]
        value$h.3$h.3.marg.curr <- value$h.3.marg$h.3.marg[value$h.3.marg$time1==0]

      } else{
        index_mat <- expand.grid(1:length(value$h.3$h.3),
                                 1:length(value$h.1$h.1))
        value$h.3.marg <- as.data.frame(cbind(
          time1=value$h.1$time[index_mat[,2]],
          time2=value$h.1$time[index_mat[,2]] + value$h.3$time[index_mat[,1]],
          h.3.marg=value$h.3$h.3[index_mat[,1]] * (1+theta) / #note extra (1+theta) term
            (1 + theta * (value$H.1$H.1[index_mat[,2]] +
                            value$H.2$H.2[index_mat[,2]] +
                            value$H.3$H.3[index_mat[,1]]))

        ))

        #add on a "time zero" set showing the marginal h3
        #assuming non-terminal event occurred "immediately"
        value$h.3.marg <- rbind(
          cbind(time1=0,
                time2=value$h.3$time,
                h.3.marg=value$h.3$h.3 * (1+theta) / (1+theta*value$H.3$H.3)),
          value$h.3.marg)

        #for convenience, add functions where the non-terminal event
        #event happened at time 0
        value$h.3$h.3.marg0 <- value$h.3.marg$h.3.marg[value$h.3.marg$time1==0]
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

