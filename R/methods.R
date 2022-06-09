
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
vcov.Freq_HReg2 <- function (object, ...){
  val <- object$Finv
  rownames(val) <- colnames(val) <- names(object$estimate)
  val
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








pred_helper <- function(tseq, xnew, pred_type, conf.level,
                        phi, beta, hess_sub, hazard_lab, knots_vec,
                        n_quad, quad_method){
  eta <- if(!is.null(xnew) & length(beta) > 0) as.vector(xnew %*% beta) else 0
  if(hazard_lab=="Weibull"){
    kappa <- exp(phi[1])
    log_alpha <- phi[2]
    alpha <- exp(log_alpha)
    if(pred_type=="S"){
      out_temp <- exp(-kappa * (tseq)^alpha * exp(eta))
      #Use delta method to compute baseline survival confidence intervals
      #matrix with as many columns as parameters in S1, and as many rows as times in tseq
      #under weibull, log(-log(S(t))) = log_alpha + log_kappa + xtbeta + alpha*log(t)
      #so, each row of J is the gradient at a particular tseq wrt logkappa, logalpha, beta1, ..., betap
      J <- cbind(1, exp(log_alpha) * log(tseq),
                 if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow=length(tseq),
                                                   ncol=length(xnew), byrow=T))
    } else if(pred_type=="h"){
      #hazard is exp(logalpha) * exp(logkappa) * t^(exp(logalpha)-1) * exp(xtbeta)
      out_temp <- alpha * kappa * (tseq)^(alpha - 1) * exp(eta)
      #Use delta method to compute baseline hazard confidence intervals
      #each row of J is gradient wrt logkappa, logalpha, beta1, ..., betap
      J <- cbind(out_temp, out_temp * (1 + alpha * log(tseq)),
                 if(!is.null(xnew) & length(beta) > 0) out_temp * matrix(xnew, nrow=length(tseq),
                                                              ncol=length(xnew), byrow=T))
    } else stop("'pred_type' must be 'Surv' or 'Haz'")
  } else if(hazard_lab=="Piecewise Constant"){
    if(pred_type=="S"){
      basis <- get_basis(x = tseq,knots = knots_vec,hazard = "piecewise")
      Lambda0 <- as.vector(basis %*% exp(phi))
      out_temp <- exp(-Lambda0*exp(eta))
      #matrix with as many columns as parameters in S1, and as many rows as times in tseq
      #under piecewise, log(-log(S(t))) = log( basis1 %*% exp(phi1)) + xtbeta
      #so, each row of J is the gradient at a particular tseq wrt phi1, ..., phiK, beta1, ..., betap
      #For phi partial derivs, first, multiply every row of basis1 by exp(phi1),
      #then, divide every column of the resulting matrix by the vector Lambda01
      J <- cbind(t(t(basis) * exp(phi)) / Lambda0,
                 if (!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                     ncol = length(xnew), byrow = T))
    } else if(pred_type=="h"){
      dbasis <- get_basis(x = tseq,knots = knots_vec,hazard = "piecewise",deriv = TRUE)
      out_temp <- as.vector(dbasis %*% exp(phi) * exp(eta))
      J <- cbind(t(t(dbasis) * exp(phi)) * exp(eta),
                 if (!is.null(xnew) & length(beta) > 0) out_temp * matrix(xnew, nrow = length(tseq),
                                                              ncol = length(xnew), byrow = T))
    } else stop("'pred_type' must be 'Surv' or 'Haz'")
  } else if(hazard_lab=="Royston-Parmar"){
    basis <- get_basis(x = tseq,knots = knots_vec,hazard = "royston-parmar")
    Lambda0 <- as.vector(exp(basis %*% phi))
    if(pred_type=="S"){
      out_temp <- exp(-Lambda0*exp(eta))
      J <- cbind(basis,
                 if(!is.null(xnew) & length(beta) > 0) matrix(xnew, nrow = length(tseq),
                                                    ncol = length(xnew), byrow = T))
    } else if(pred_type=="h"){
      dbasis <- get_basis(x = tseq,knots = knots_vec,hazard = "royston-parmar",deriv = TRUE)
      #Use delta method compute baseline hazard confidence intervals
      # h(t) = s'(z)/t * exp(s(z) + xtbeta) where z=log(t)
      out_temp <- as.vector(dbasis %*% phi / tseq * exp(basis %*% phi + eta))
      #so, each row of J is gradient wrt phi1, ..., phik, beta1, ..., betap
      J <- cbind(dbasis*Lambda0*exp(eta)/tseq + out_temp*basis,
                 if(!is.null(xnew) & length(beta) > 0) out_temp * matrix(xnew, nrow = length(tseq),
                                                               ncol = length(xnew), byrow = T))
    } else stop("'pred_type' must be 'Surv' or 'Haz'")
  } else if(hazard_lab=="B-Spline"){
    quad_weights <- get_quad_pointsweights(n_quad=n_quad,
                                           quad_method=quad_method)$weights
    quad_points <- transform_quad_points(n_quad = n_quad,
                                         quad_method=quad_method, a = 0,b = tseq)
    basis <- get_basis(x=tseq,knots=knots_vec,hazard="bspline")
    basis_quad <- get_basis(x=quad_points, knots=knots_vec,hazard="bspline")
    lambda0 <- as.vector(exp(basis_quad %*% phi))
    if(pred_type=="S"){
      #reshape lambda0 from a n*n_quad length vector
      #to an n by n_quad matrix, then multiply with n_quad length weights
      #to get final Lambda0
      Lambda0 <- tseq/2 * as.vector(matrix(lambda0,ncol=n_quad,byrow = TRUE) %*% quad_weights)
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
    } else if(pred_type=="h"){
      #Use delta method compute baseline hazard confidence intervals
      out_temp <- as.vector(exp(basis %*% phi) * exp(eta))
      #hazard is exp(Btphi) * exp(xtbeta)
      #so, each row of J is gradient wrt phi1, ..., phik, beta1, ..., betap
      J <- cbind(out_temp * basis,
                 if(!is.null(xnew) & length(beta) > 0) out_temp * matrix(xnew, nrow = length(tseq),
                                                               ncol = length(xnew), byrow = T))
    } else stop("'pred_type' must be 'Surv' or 'Haz'")
  } else stop("'hazard_lab' must be 'Weibull', 'Piecewise Constant', 'Royston-Parmar' or 'B-Spline'")

  out_var <- if(all(is.na(hess_sub))) NA else J %*% hess_sub %*% t(J)
  out_se <- if(all(is.na(out_var))) NA else sqrt(diag(out_var))
  out_se[is.nan(out_se)] <- 0
  if(pred_type=="S"){
    LL <- out_temp^exp(-stats::qnorm(conf.level/2) * out_se)
    UL <- out_temp^exp(stats::qnorm(conf.level/2) * out_se)
  } else if(pred_type=="h"){
    LL <- out_temp + stats::qnorm(conf.level/2) * out_se #sign reversed because 0.025 quantile is negative
    UL <- out_temp - stats::qnorm(conf.level/2) * out_se
    LL[LL < 0] <- 0
    if (tseq[1] == 0) {
      tseq <- tseq[-1]
      out_temp <- out_temp[-1]
      LL <- LL[-1]
      UL <- UL[-1]
    }
  } else stop("'pred_type' must be 'Surv' or 'Haz'")
  return(data.frame(time=tseq,out=out_temp,LL=LL,UL=UL))
}

#' @export
predict.Freq_HReg2 <- function (object, xnew = NULL,
                                x1new = NULL, x2new = NULL, x3new = NULL,
                                tseq = seq(0,object$ymax,length.out = 100),
                                alpha = 0.05, n_quad=15, quad_method="kronrod", ...) {
  # browser()
  conf.level = alpha
  obj <- object
  yLim <- NULL
  nP = obj$nP
  nP0 = obj$nP0
  if (obj$class[2] == "Surv") {
    if (!(is.null(x1new) & is.null(x2new) & is.null(x3new))) {
      stop("'x1new','x2new', and 'x1new' are for semi-competing risks models and must be specified as NULL for univariate models")
    }
    value <- list()
    for(type in c("h","S")){
      value[[type]] <- pred_helper(tseq = tseq,xnew = xnew,
                            pred_type = type,conf.level = conf.level,
                            phi = obj$estimate[1:nP0],
                            beta = obj$estimate[-(1:nP0)],
                            hess_sub = obj$Finv[c(1:nP0, if(!is.null(xnew)) (nP0+1):(nP0+nP)),
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
            pred_type = type,conf.level = conf.level,
            phi = obj$estimate[nP0_start[i]:nP0_end[i]],
            beta = obj$estimate[nP_start[i]:nP_end[i]],
            hess_sub = obj$Finv[c(nP0_start[i]:nP0_end[i],
                                  if(!is.null(x_list[[i]])) nP_start[i]:nP_end[i]),
                                c(nP0_start[i]:nP0_end[i],
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





# #
# #' @export
# predict2.Freq_HReg2 <- function (object, xnew = NULL,
#                                  x1new = NULL, x2new = NULL, x3new = NULL,
#                                  tseq = c(0, 5, 10), alpha = 0.05, ...) {
#   # browser()
#   conf.level = alpha
#   obj <- object
#   yLim <- NULL
#   if (obj$class[2] == "ID") {
#     if (!is.null(xnew)) {
#       stop("'xnew' is for univariate models so it must be specified as NULL for semi-competing risks models")
#     }
#     nP = obj$nP
#     nP0 = obj$nP0
#     nP0_tot <- if (obj$frailty == TRUE) sum(nP0) + 1 else sum(nP0)
#     T2 <- tseq
#
#     eta1 <- if (!is.null(x1new) & nP[1] > 0) {
#       x1new %*% obj$estimate[(1 + nP0_tot):(nP0_tot + nP[1])]
#     } else 0
#     eta2 <- if (!is.null(x2new) & nP[2] > 0) {
#       x2new %*% obj$estimate[(1 + nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])]
#     } else 0
#     eta3 <- if (!is.null(x3new) & nP[3] > 0) {
#       x3new %*% obj$estimate[(1 + nP0_tot + nP[1] + nP[2]):(nP0_tot + nP[1] + nP[2] + nP[3])]
#     } else 0
#
#     #FIRST TRANSITION
#     if(obj$class[4] == "Weibull"){
#       kappa <- exp(obj$estimate[1])
#       alpha <- exp(obj$estimate[2])
#       log_alpha <- obj$estimate[2]
#       S.1 <- exp(-kappa * (T2)^alpha) * exp(eta1)
#       #Use delta method to compute baseline survival confidence intervals
#       if(all(is.na(obj$Finv))){
#         Var.loglogS.1 <- NA
#       } else if (!is.null(x1new) & nP[1] > 0) {
#         #matrix with as many columns as parameters in S1, and as many rows as times in T2
#         #under weibull, log(-log(S(t))) = log_alpha + log_kappa + xtbeta + alpha*log(t)
#         #so, each row of J is the gradient at a particular T2 wrt logkappa, logalpha, beta1, ..., betap
#         J <- cbind(1, exp(log_alpha) * log(T2),
#                    matrix(x1new, nrow = length(T2), ncol = length(x1new), byrow = T))
#         #vector with as many elements as times in T2
#         Var.loglogS.1 <- J %*% obj$Finv[c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1])),
#                                         c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1]))] %*% t(J)
#       } else if (is.null(x1new) | nP[1] == 0) {
#         J <- cbind(1, exp(log_alpha) * log(T2))
#         Var.loglogS.1 <- J %*% obj$Finv[1:nP0[1], 1:nP0[1]] %*% t(J)
#       }
#
#       h.1 <- alpha * kappa * (T2)^(alpha - 1) * exp(eta1)
#       #hazard is exp(logalpha) * exp(logkappa) * t^(exp(logalpha)-1) * exp(xtbeta)
#       #so, each row of J is gradient wrt logkappa, logalpha, beta1, ..., betap
#       if(all(is.na(obj$Finv))){
#         Var.h.1 <- NA
#       } else if (!is.null(x1new) & nP[1] > 0) {
#         J <- cbind(h.1, h.1 * (1 + alpha * log(T2)), h.1 *
#                      matrix(x1new, nrow = length(T2), ncol = length(x1new), byrow = T))
#         Var.h.1 <- J %*% obj$Finv[c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1])),
#                                   c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1]))] %*% t(J)
#       }
#       else if (is.null(x1new) | nP[1] == 0) {
#         J <- cbind(h.1, h.1 * (1 + alpha * log(T2)))
#         Var.h.1 <- J %*% obj$Finv[1:nP0[1], 1:nP0[1]] %*% t(J)
#       }
#     } else if(obj$class[4] == "Piecewise Constant"){
#       basis1 <- get_basis(x = T2,knots = obj$knots_list[[1]],hazard = "piecewise")
#       phi1 <- obj$estimate[1:nP0[1]]
#       Lambda01 <- as.vector(basis1 %*% exp(phi1))
#       S.1 <- exp(-Lambda01*exp(eta1))
#       #Use delta method to compute baseline survival confidence intervals
#       if(all(is.na(obj$Finv))){
#         Var.loglogS.1 <- NA
#       } else if (!is.null(x1new) & nP[1] > 0) {
#         #matrix with as many columns as parameters in S1, and as many rows as times in T2
#         #under piecewise, log(-log(S(t))) = log( basis1 %*% exp(phi1)) + xtbeta
#         #so, each row of J is the gradient at a particular T2 wrt phi1, ..., phiK, beta1, ..., betap
#         J <- cbind(
#           #first, this multiplies every row of basis1 by exp(phi1),
#           #then, this divides every column of the resulting matrix by the vector Lambda01
#           t(t(basis1) * exp(phi1)) / Lambda01,
#           matrix(x1new, nrow = length(T2), ncol = length(x1new), byrow = T))
#         #vector with as many elements as times in T2
#         Var.loglogS.1 <- J %*% obj$Finv[c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1])),
#                                         c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1]))] %*%
#           t(J)
#       } else if (is.null(x1new) | nP[1] == 0) {
#         J <- t(t(basis1) * exp(phi1)) / Lambda01
#         Var.loglogS.1 <- J %*% obj$Finv[1:nP0[1], 1:nP0[1]] %*% t(J)
#       }
#
#       #Use delta method compute baseline hazard confidence intervals
#       #vector saying which interval each time falls into
#       cut_cats1 <- rowSums(basis1!=0)
#       if(T2[1]==0){cut_cats1[1] <- 1}
#       stopifnot(length(cut_cats1)==length(T2))
#       h.1 <- exp(phi1)[cut_cats1] * exp(eta1)
#       #build a matrix with 0's everywhere, except on ith row, set column that T2i falls in to 1
#       temp_mat <- matrix(data=0,nrow=length(T2),ncol=nP0[1])
#       temp_mat[cbind(1:length(T2),cut_cats1)] <- 1
#       temp_mat <- t(t(temp_mat) * exp(phi1))
#
#       if(all(is.na(obj$Finv))){
#         Var.h.1 <- NA
#       } else if (!is.null(x1new) & nP[1] > 0) {
#         J <- cbind(temp_mat * exp(eta1), h.1 *
#                      matrix(x1new, nrow = length(T2), ncol = length(x1new), byrow = T))
#         Var.h.1 <- J %*% obj$Finv[c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1])),
#                                   c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1]))] %*% t(J)
#       }
#       else if (is.null(x1new) | nP[1] == 0) {
#         J <- temp_mat * exp(eta1)
#         Var.h.1 <- J %*% obj$Finv[1:nP0[1], 1:nP0[1]] %*% t(J)
#       }
#     } else if(obj$class[4] == "Royston-Parmar"){
#       basis1 <- get_basis(x = T2,knots = obj$knots_list[[1]],hazard = "royston-parmar")
#       dbasis1 <- get_basis(x = T2,knots = obj$knots_list[[1]],hazard = "royston-parmar",deriv = TRUE)
#       phi1 <- obj$estimate[1:nP0[1]]
#       Lambda01 <- as.vector(exp(basis1 %*% phi1))
#       S.1 <- exp(-Lambda01*exp(eta1))
#
#       #Use delta method to compute baseline survival confidence intervals
#       if(all(is.na(obj$Finv))){
#         Var.loglogS.1 <- NA
#       } else if (!is.null(x1new) & nP[1] > 0) {
#         #matrix with as many columns as parameters in S1, and as many rows as times in T2
#         #under bspline, log(-log(S(t))) = Btphi + xtbeta
#         #so, each row of J is the gradient at a particular T2 wrt phi1, ..., phiK, beta1, ..., betap
#         #specifically, for gradient of phi we take each basis column and multiply it by lambda,
#         #evaluate the numerical integral of the product, and then divide off Lambda01
#         J <- cbind(
#           #first, this multiplies every row of basis1 by exp(phi1),
#           #then, this divides every column of the resulting matrix by the vector Lambda01
#           basis1,
#           matrix(x1new, nrow = length(T2), ncol = length(x1new), byrow = T))
#         #vector with as many elements as times in T2
#         Var.loglogS.1 <- J %*% obj$Finv[c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1])),
#                                         c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1]))] %*%
#           t(J)
#       } else if (is.null(x1new) | nP[1] == 0) {
#         J <- basis1
#         Var.loglogS.1 <- J %*% obj$Finv[1:nP0[1], 1:nP0[1]] %*% t(J)
#       }
#
#       #Use delta method compute baseline hazard confidence intervals
#       # h(t) = s'(z)/t * exp(s(z) + xtbeta) where z=log(t)
#       h.1 <- as.vector(dbasis1 %*% phi1 / T2 * exp(basis1 %*% phi1 + eta1))
#       #so, each row of J is gradient wrt phi1, ..., phik, beta1, ..., betap
#       if(all(is.na(obj$Finv))){
#         Var.h.1 <- NA
#       } else if (!is.null(x1new) & nP[1] > 0) {
#         J <- cbind(dbasis1*Lambda01*exp(eta1)/T2 + h.1*basis1, h.1 *
#                      matrix(x1new, nrow = length(T2), ncol = length(x1new), byrow = T))
#         Var.h.1 <- J %*% obj$Finv[c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1])),
#                                   c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1]))] %*% t(J)
#       }
#       else if (is.null(x1new) | nP[1] == 0) {
#         J <- dbasis1*Lambda01*exp(eta1)/T2 + h.1*basis1
#         Var.h.1 <- J %*% obj$Finv[1:nP0[1], 1:nP0[1]] %*% t(J)
#       }
#     } else if(obj$class[4] == "B-Spline"){
#       basis1 <- get_basis(x = T2,knots = obj$knots_list[[1]],hazard = "bspline")
#       basis1_quad <- get_basis(x = transform_quad_points(n_quad = obj$n_quad,
#                                                          quad_method=obj$quad_method,
#                                                          a = 0,b = T2),
#                                knots = obj$knots_list[[1]],hazard = "bspline")
#       phi1 <- obj$estimate[1:nP0[1]]
#       lambda01 <- as.vector(exp(basis1_quad %*% phi1))
#
#       #reshape lambda01 from a n*n_quad length vector
#       #to an n by n_quad matrix, then multiply with n_quad length weights
#       #to get final Lambda01
#       Lambda01 <- T2/2 * as.vector(matrix(lambda01,ncol=obj$n_quad,byrow = TRUE) %*% get_quad_pointsweights(n_quad=obj$n_quad)$weights)
#       S.1 <- exp(-Lambda01*exp(eta1))
#
#       #Use delta method to compute baseline survival confidence intervals
#       if(all(is.na(obj$Finv))){
#         Var.loglogS.1 <- NA
#       } else if (!is.null(x1new) & nP[1] > 0) {
#         #matrix with as many columns as parameters in S1, and as many rows as times in T2
#         #under bspline, log(-log(S(t))) = log( [numerical integral of exp(Btphi)] ) + xtbeta
#         #so, each row of J is the gradient at a particular T2 wrt phi1, ..., phiK, beta1, ..., betap
#         #specifically, for gradient of phi we take each basis column and multiply it by lambda,
#         #evaluate the numerical integral of the product, and then divide off Lambda01
#         J_phi <- apply(X = basis1_quad * lambda01, MARGIN = 2,
#                        FUN = function(x) T2/2 * as.vector(matrix(x,ncol=obj$n_quad,byrow = TRUE) %*% get_quad_pointsweights(n_quad=obj$n_quad)$weights))
#         J <- cbind(
#           #first, this multiplies every row of basis1 by exp(phi1),
#           #then, this divides every column of the resulting matrix by the vector Lambda01
#           J_phi / Lambda01,
#           matrix(x1new, nrow = length(T2), ncol = length(x1new), byrow = T))
#         #vector with as many elements as times in T2
#         Var.loglogS.1 <- J %*% obj$Finv[c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1])),
#                                         c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1]))] %*%
#           t(J)
#       } else if (is.null(x1new) | nP[1] == 0) {
#         J_phi <- apply(X = basis1_quad * lambda01, MARGIN = 2,
#                        FUN = function(x) T2/2 * as.vector(matrix(x,ncol=obj$n_quad,byrow = TRUE) %*% get_quad_pointsweights(n_quad=obj$n_quad)$weights))
#         J <- J_phi / Lambda01
#         Var.loglogS.1 <- J %*% obj$Finv[1:nP0[1], 1:nP0[1]] %*% t(J)
#       }
#
#       #Use delta method compute baseline hazard confidence intervals
#       h.1 <- as.vector(exp(basis1 %*% phi1) * exp(eta1))
#       #hazard is exp(Btphi) * exp(xtbeta)
#       #so, each row of J is gradient wrt phi1, ..., phik, beta1, ..., betap
#       if(all(is.na(obj$Finv))){
#         Var.h.1 <- NA
#       } else if (!is.null(x1new) & nP[1] > 0) {
#         J <- cbind(h.1 * basis1, h.1 *
#                      matrix(x1new, nrow = length(T2), ncol = length(x1new), byrow = T))
#         Var.h.1 <- J %*% obj$Finv[c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1])),
#                                   c(1:nP0[1], (1+nP0_tot):(nP0_tot + nP[1]))] %*% t(J)
#       }
#       else if (is.null(x1new) | nP[1] == 0) {
#         J <- h.1 * basis1
#         Var.h.1 <- J %*% obj$Finv[1:nP0[1], 1:nP0[1]] %*% t(J)
#       }
#     } else{stop("class must be Weibull, Piecewise Constant, or B-Spline")}
#
#     if(all(is.na(Var.loglogS.1))){
#       se.loglogS.1 <- NA
#     } else{
#       se.loglogS.1 <- sqrt(diag(Var.loglogS.1))
#       se.loglogS.1[is.nan(se.loglogS.1)] <- 0
#     }
#     LL.1 <- S.1^exp(-stats::qnorm(conf.level/2) * se.loglogS.1)
#     UL.1 <- S.1^exp(stats::qnorm(conf.level/2) * se.loglogS.1)
#     if(all(is.na(Var.h.1))){
#       se.h.1 <- NA
#     } else{
#       se.h.1 <- sqrt(diag(Var.h.1))
#       se.h.1[is.nan(se.h.1)] <- 0
#     }
#     LLh.1 <- h.1 + stats::qnorm(conf.level/2) * se.h.1 #sign reversed because 0.025 quantile is negative
#     ULh.1 <- h.1 - stats::qnorm(conf.level/2) * se.h.1
#     LLh.1[LLh.1 < 0] <- 0
#
#     #SECOND TRANSITION
#     if(obj$class[4] == "Weibull"){
#       kappa <- exp(obj$estimate[3])
#       alpha <- exp(obj$estimate[4])
#       log_alpha <- obj$estimate[4]
#       S.2 <- exp(-kappa * (T2)^alpha) * exp(eta2)
#       #Use delta method to compute baseline survival confidence intervals
#       if(all(is.na(obj$Finv))){
#         Var.loglogS.2 <- NA
#       } else if (!is.null(x2new) & nP[2] > 0) {
#         #matrix with as many columns as parameters in S2, and as many rows as times in T2
#         #under weibull, log(-log(S(t))) = log_alpha + log_kappa + xtbeta + alpha*log(t)
#         #so, each row of J is the gradient at a particular T2 wrt logalpha, logkappa, beta1, ..., betap
#         J <- cbind(1, exp(log_alpha) * log(T2),
#                    matrix(x2new, nrow = length(T2), ncol = length(x2new), byrow = T))
#         #vector with as many elements as times in T2
#         Var.loglogS.2 <- J %*% obj$Finv[c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])),
#                                         c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2]))] %*% t(J)
#       } else if (is.null(x2new) | nP[2] == 0) {
#         J <- cbind(1, exp(log_alpha) * log(T2))
#         Var.loglogS.2 <- J %*% obj$Finv[(1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0[1]):(nP0[1]+nP0[2])] %*% t(J)
#       }
#
#       h.2 <- alpha * kappa * (T2)^(alpha - 1) * exp(eta2)
#       #hazard is exp(logalpha) * exp(logkappa) * t^(exp(logalpha)-1) * exp(xtbeta)
#       #so, each row of J is gradient wrt logkappa, logalpha, beta1, ..., betap
#       #WHY ISNT EVERYTHING MULTIPLIED BY exp(xtbeta) ? seems like it should be. I'm gonna add it
#       if(all(is.na(obj$Finv))){
#         Var.h.2 <- NA
#       } else if (!is.null(x2new) & nP[2] > 0) {
#         J <- cbind(h.2, h.2 * (1 + alpha * log(T2)), h.2 *
#                      matrix(x2new, nrow = length(T2), ncol = length(x2new), byrow = T))
#         Var.h.2 <- J %*% obj$Finv[c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])),
#                                   c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2]))] %*% t(J)
#       }
#       else if (is.null(x2new) | nP[2] == 0) {
#         J <- cbind(h.2, h.2 * (1 + alpha * log(T2)))
#         Var.h.2 <- J %*% obj$Finv[(1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0[1]):(nP0[1]+nP0[2])] %*% t(J)
#       }
#     } else if(obj$class[4] == "Piecewise Constant"){ #for now the only alternative we entertain is piecewise constant
#       basis2 <- get_basis(x = T2,knots = obj$knots_list[[2]],hazard = "piecewise")
#       phi2 <- obj$estimate[(1+nP0[1]):(nP0[1]+nP0[2])]
#       Lambda02 <- as.vector(basis2 %*% exp(phi2))
#       S.2 <- exp(-Lambda02*exp(eta2))
#       #Use delta method to compute baseline survival confidence intervals
#       if(all(is.na(obj$Finv))){
#         Var.loglogS.2 <- NA
#       } else if (!is.null(x2new) & nP[2] > 0) {
#         #matrix with as many columns as parameters in S1, and as many rows as times in T2
#         #under piecewise, log(-log(S(t))) = log( basis1 %*% exp(phi1)) + log(xtbeta)
#         #so, each row of J is the gradient at a particular T2 wrt phi1, ..., phiK, beta1, ..., betap
#         J <- cbind(
#           #first, this multiplies every row of basis1 by exp(phi1),
#           #then, this divides every column of the resulting matrix by the vector Lambda01
#           t(t(basis2) * exp(phi2)) / Lambda02,
#           matrix(x2new, nrow = length(T2), ncol = length(x2new), byrow = T))
#         #vector with as many elements as times in T2
#         Var.loglogS.2 <- J %*% obj$Finv[c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])),
#                                         c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2]))] %*%
#           t(J)
#       } else if (is.null(x2new) | nP[2] == 0) {
#         J <- t(t(basis2) * exp(phi2)) / Lambda02
#         Var.loglogS.2 <- J %*% obj$Finv[(1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0[1]):(nP0[1]+nP0[2])] %*% t(J)
#       }
#
#       #Use delta method compute baseline hazard confidence intervals
#       #vector saying which interval each time falls into
#       cut_cats2 <- rowSums(basis2!=0)
#       if(T2[1]==0){cut_cats2[1] <- 1}
#       stopifnot(length(cut_cats2)==length(T2))
#       h.2 <- exp(phi2)[cut_cats2] * exp(eta2)
#       #build a matrix with 0's everywhere, except on ith row, set column that T2i falls in to 1
#       temp_mat <- matrix(data=0,nrow=length(T2),ncol=nP0[2])
#       temp_mat[cbind(1:length(T2),cut_cats2)] <- 1
#       temp_mat <- t(t(temp_mat) * exp(phi2))
#
#       if(all(is.na(obj$Finv))){
#         Var.h.2 <- NA
#       } else if (!is.null(x2new) & nP[2] > 0) {
#         J <- cbind(temp_mat * exp(eta2), h.2 *
#                      matrix(x2new, nrow = length(T2), ncol = length(x2new), byrow = T))
#         Var.h.2 <- J %*% obj$Finv[c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])),
#                                   c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2]))] %*% t(J)
#       }
#       else if (is.null(x2new) | nP[2] == 0) {
#         J <- temp_mat * exp(eta2)
#         Var.h.2 <- J %*% obj$Finv[(1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0[1]):(nP0[1]+nP0[2])] %*% t(J)
#       }
#     } else if(obj$class[4] == "Royston-Parmar"){
#       basis2 <- get_basis(x = T2,knots = obj$knots_list[[2]],hazard = "royston-parmar")
#       dbasis2 <- get_basis(x = T2,knots = obj$knots_list[[2]],hazard = "royston-parmar",deriv = TRUE)
#       phi2 <- obj$estimate[(1+nP0[1]):(nP0[1]+nP0[2])]
#       Lambda02 <- as.vector(exp(basis2 %*% phi2))
#       S.2 <- exp(-Lambda02*exp(eta2))
#
#       #Use delta method to compute baseline survival confidence intervals
#       if(all(is.na(obj$Finv))){
#         Var.loglogS.2 <- NA
#       } else if (!is.null(x2new) & nP[2] > 0) {
#         #matrix with as many columns as parameters in S1, and as many rows as times in T2
#         #under bspline, log(-log(S(t))) = Btphi + xtbeta
#         #so, each row of J is the gradient at a particular T2 wrt phi1, ..., phiK, beta1, ..., betap
#         #specifically, for gradient of phi we take each basis column and multiply it by lambda,
#         #evaluate the numerical integral of the product, and then divide off Lambda01
#         J <- cbind(
#           #first, this multiplies every row of basis1 by exp(phi1),
#           #then, this divides every column of the resulting matrix by the vector Lambda01
#           basis2,
#           matrix(x2new, nrow = length(T2), ncol = length(x2new), byrow = T))
#         #vector with as many elements as times in T2
#         Var.loglogS.2 <- J %*% obj$Finv[c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])),
#                                         c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2]))] %*%
#           t(J)
#       } else if (is.null(x2new) | nP[2] == 0) {
#         J <- basis2
#         Var.loglogS.2 <- J %*% obj$Finv[(1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0[1]):(nP0[1]+nP0[2])] %*% t(J)
#       }
#
#       #Use delta method compute baseline hazard confidence intervals
#       # h(t) = s'(z)/t * exp(s(z) + xtbeta) where z=log(t)
#       h.2 <- as.vector(dbasis2 %*% phi2 / T2 * exp(basis2 %*% phi2 + eta2))
#       #so, each row of J is gradient wrt phi1, ..., phik, beta1, ..., betap
#       if(all(is.na(obj$Finv))){
#         Var.h.2 <- NA
#       } else if (!is.null(x2new) & nP[2] > 0) {
#         J <- cbind(dbasis2*Lambda02*exp(eta2)/T2 + h.2*basis2, h.2 *
#                      matrix(x2new, nrow = length(T2), ncol = length(x2new), byrow = T))
#         Var.h.2 <- J %*% obj$Finv[c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])),
#                                   c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2]))] %*% t(J)
#       }
#       else if (is.null(x2new) | nP[2] == 0) {
#         J <- dbasis2*Lambda02*exp(eta2)/T2 + h.2*basis2
#         Var.h.2 <- J %*% obj$Finv[(1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0[1]):(nP0[1]+nP0[2])] %*% t(J)
#       }
#     } else if(obj$class[4] == "B-Spline"){
#       basis2 <- get_basis(x = T2,knots = obj$knots_list[[2]],hazard = "bspline")
#       basis2_quad <- get_basis(x = transform_quad_points(n_quad = obj$n_quad,
#                                                          quad_method = obj$quad_method,
#                                                          a = 0,b = T2),
#                                knots = obj$knots_list[[2]],hazard = "bspline")
#       phi2 <- obj$estimate[(1+nP0[1]):(nP0[1]+nP0[2])]
#       lambda02 <- as.vector(exp(basis2_quad %*% phi2))
#
#       #reshape lambda01 from a n*n_quad length vector
#       #to an n by n_quad matrix, then multiply with n_quad length weights
#       #to get final Lambda01
#       Lambda02 <- T2/2 * as.vector(matrix(lambda02,ncol=obj$n_quad,byrow = TRUE) %*% get_quad_pointsweights(n_quad=obj$n_quad)$weights)
#       S.2 <- exp(-Lambda02*exp(eta2))
#
#       #Use delta method to compute baseline survival confidence intervals
#       if(all(is.na(obj$Finv))){
#         Var.loglogS.2 <- NA
#       } else if (!is.null(x2new) & nP[2] > 0) {
#         #matrix with as many columns as parameters in S1, and as many rows as times in T2
#         #under piecewise, log(-log(S(t))) = log( basis1 %*% exp(phi1)) + xtbeta
#         #so, each row of J is the gradient at a particular T2 wrt phi1, ..., phiK, beta1, ..., betap
#         J_phi <- apply(X = basis2_quad * lambda02, MARGIN = 2,
#                        FUN = function(x) T2/2 * as.vector(matrix(x,ncol=obj$n_quad,byrow = TRUE) %*% get_quad_pointsweights(n_quad=obj$n_quad)$weights))
#         J <- cbind(
#           #first, this multiplies every row of basis1 by exp(phi1),
#           #then, this divides every column of the resulting matrix by the vector Lambda01
#           J_phi / Lambda02,
#           matrix(x2new, nrow = length(T2), ncol = length(x2new), byrow = T))
#         #vector with as many elements as times in T2
#         Var.loglogS.2 <- J %*% obj$Finv[c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])),
#                                         c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2]))] %*%
#           t(J)
#       } else if (is.null(x2new) | nP[2] == 0) {
#         J_phi <- apply(X = basis2_quad * lambda02, MARGIN = 2,
#                        FUN = function(x) T2/2 * as.vector(matrix(x,ncol=obj$n_quad,byrow = TRUE) %*% get_quad_pointsweights(n_quad=obj$n_quad)$weights))
#         J <- J_phi / Lambda02
#         Var.loglogS.2 <- J %*% obj$Finv[(1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0[1]):(nP0[1]+nP0[2])] %*% t(J)
#       }
#
#       #Use delta method compute baseline hazard confidence intervals
#       h.2 <- as.vector(exp(basis2 %*% phi2) * exp(eta2))
#       #hazard is exp(Btphi) * exp(xtbeta)
#       #so, each row of J is gradient wrt phi1, ..., phik, beta1, ..., betap
#       if(all(is.na(obj$Finv))){
#         Var.h.2 <- NA
#       } else if (!is.null(x2new) & nP[2] > 0) {
#         J <- cbind(h.2 * basis2, h.2 *
#                      matrix(x2new, nrow = length(T2), ncol = length(x2new), byrow = T))
#         Var.h.2 <- J %*% obj$Finv[c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2])),
#                                   c((1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0_tot + nP[1]):(nP0_tot + nP[1] + nP[2]))]  %*% t(J)
#       }
#       else if (is.null(x2new) | nP[2] == 0) {
#         J <- h.2 * basis2
#         Var.h.2 <- J %*% obj$Finv[(1+nP0[1]):(nP0[1]+nP0[2]), (1+nP0[1]):(nP0[1]+nP0[2])] %*% t(J)
#       }
#     } else{stop("class must be Weibull, Piecewise Constant, or B-Spline")}
#
#
#     if(all(is.na(Var.loglogS.2))){
#       se.loglogS.2 <- NA
#     } else{
#       se.loglogS.2 <- sqrt(diag(Var.loglogS.2))
#       se.loglogS.2[is.nan(se.loglogS.2)] <- 0
#     }
#     LL.2 <- S.2^exp(-stats::qnorm(conf.level/2) * se.loglogS.2)
#     UL.2 <- S.2^exp(stats::qnorm(conf.level/2) * se.loglogS.2)
#     if(all(is.na(Var.h.2))){
#       se.h.2 <- NA
#     } else{
#       se.h.2 <- sqrt(diag(Var.h.2))
#       se.h.2[is.nan(se.h.2)] <- 0
#     }
#     LLh.2 <- h.2 + stats::qnorm(conf.level/2) * se.h.2 #sign reversed because 0.025 quantile is negative
#     ULh.2 <- h.2 - stats::qnorm(conf.level/2) * se.h.2
#     LLh.2[LLh.2 < 0] <- 0
#
#     #THIRD TRANSITION
#     if(obj$class[4] == "Weibull"){
#       kappa <- exp(obj$estimate[5])
#       alpha <- exp(obj$estimate[6])
#       log_alpha <- obj$estimate[6]
#       S.3 <- exp(-kappa * (T2)^alpha) * exp(eta3)
#       #Use delta method to compute baseline survival confidence intervals
#       if(all(is.na(obj$Finv))){
#         Var.loglogS.3 <- NA
#       } else if(!is.null(x3new) & nP[3] > 0) {
#         #matrix with as many columns as parameters in S3, and as many rows as times in T2
#         #under weibull, log(-log(S(t))) = log_alpha + log_kappa + xtbeta + alpha*log(t)
#         #so, each row of J is the gradient at a particular T2 wrt logalpha, logkappa, beta1, ..., betap
#         J <- cbind(1, exp(log_alpha) * log(T2),
#                    matrix(x3new, nrow = length(T2), ncol = length(x3new), byrow = T))
#         #vector with as many elements as times in T2
#         Var.loglogS.3 <- J %*% obj$Finv[c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3])),
#                                         c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3]))] %*% t(J)
#       } else if (is.null(x3new) | nP[2] == 0) {
#         J <- cbind(1, exp(log_alpha) * log(T2))
#         Var.loglogS.3 <- J %*% obj$Finv[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])] %*% t(J)
#       }
#
#       h.3 <- alpha * kappa * (T2)^(alpha - 1) * exp(eta3)
#       #hazard is exp(logalpha) * exp(logkappa) * t^(exp(logalpha)-1) * exp(xtbeta)
#       #so, each row of J is gradient wrt logkappa, logalpha, beta1, ..., betap
#       #kyu ha's code does not multiply the variances by exp(xtbeta), which it seems like it should be. I'm gonna add it
#       if(all(is.na(obj$Finv))){
#         Var.h.3 <- NA
#       } else if (!is.null(x3new) & nP[3] > 0) {
#         J <- cbind(h.3, h.3 * (1 + alpha * log(T2)), h.3 *
#                      matrix(x3new, nrow = length(T2), ncol = length(x3new), byrow = T))
#         Var.h.3 <- J %*% obj$Finv[c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3])),
#                                   c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3]))] %*% t(J)
#       }
#       else if (is.null(x3new) | nP[2] == 0) {
#         J <- cbind(h.3, h.3 * (1 + alpha * log(T2)))
#         Var.h.3 <- J %*% obj$Finv[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])] %*% t(J)
#       }
#     } else if(obj$class[4] == "Piecewise Constant"){
#       basis3 <- get_basis(x = T2,knots = obj$knots_list[[3]],hazard = "piecewise")
#       phi3 <- obj$estimate[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])]
#       Lambda03 <- as.vector(basis3 %*% exp(phi3))
#       S.3 <- exp(-Lambda03*exp(eta3))
#       #Use delta method to compute baseline survival confidence intervals
#       if(all(is.na(obj$Finv))){
#         Var.loglogS.3 <- NA
#       } else if (!is.null(x3new) & nP[3] > 0) {
#         #matrix with as many columns as parameters in S1, and as many rows as times in T2
#         #under piecewise, log(-log(S(t))) = log( basis1 %*% exp(phi1)) + log(xtbeta)
#         #so, each row of J is the gradient at a particular T2 wrt phi1, ..., phiK, beta1, ..., betap
#         J <- cbind(
#           #first, this multiplies every row of basis1 by exp(phi1),
#           #then, this divides every column of the resulting matrix by the vector Lambda01
#           t(t(basis3) * exp(phi3)) / Lambda03,
#           matrix(x3new, nrow = length(T2), ncol = length(x3new), byrow = T))
#         #vector with as many elements as times in T2
#         Var.loglogS.3 <- J %*% obj$Finv[c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3])),
#                                         c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3]))] %*%
#           t(J)
#       } else if (is.null(x3new) | nP[3] == 0) {
#         J <- t(t(basis3) * exp(phi3)) / Lambda03
#         Var.loglogS.3 <- J %*% obj$Finv[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])] %*% t(J)
#       }
#
#       #Use delta method compute baseline hazard confidence intervals
#       #vector saying which interval each time falls into
#       cut_cats3 <- rowSums(basis3!=0)
#       if(T2[1]==0){cut_cats3[1] <- 1}
#       stopifnot(length(cut_cats3)==length(T2))
#       h.3 <- exp(phi3)[cut_cats3] * exp(eta3)
#       #build a matrix with 0's everywhere, except on ith row, set column that T2i falls in to 1
#       temp_mat <- matrix(data=0,nrow=length(T2),ncol=nP0[3])
#       temp_mat[cbind(1:length(T2),cut_cats3)] <- 1
#       temp_mat <- t(t(temp_mat) * exp(phi3))
#
#       if(all(is.na(obj$Finv))){
#         Var.h.3 <- NA
#       } else if (!is.null(x3new) & nP[3] > 0) {
#         J <- cbind(temp_mat * exp(eta3), h.3 *
#                      matrix(x3new, nrow = length(T2), ncol = length(x3new), byrow = T))
#         Var.h.3 <- J %*% obj$Finv[c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3])),
#                                   c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3]))] %*% t(J)
#       }
#       else if (is.null(x3new) | nP[3] == 0) {
#         J <- temp_mat * exp(eta3)
#         Var.h.3 <- J %*% obj$Finv[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])] %*% t(J)
#       }
#     } else if(obj$class[4] == "Royston-Parmar"){
#       basis3 <- get_basis(x = T2,knots = obj$knots_list[[3]],hazard = "royston-parmar")
#       dbasis3 <- get_basis(x = T2,knots = obj$knots_list[[3]],hazard = "royston-parmar",deriv = TRUE)
#       phi3 <- obj$estimate[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])]
#       Lambda03 <- as.vector(exp(basis3 %*% phi3))
#       S.3 <- exp(-Lambda03*exp(eta3))
#
#       #Use delta method to compute baseline survival confidence intervals
#       if(all(is.na(obj$Finv))){
#         Var.loglogS.3 <- NA
#       } else if (!is.null(x3new) & nP[3] > 0) {
#         #matrix with as many columns as parameters in S1, and as many rows as times in T2
#         #under bspline, log(-log(S(t))) = Btphi + xtbeta
#         #so, each row of J is the gradient at a particular T2 wrt phi1, ..., phiK, beta1, ..., betap
#         #specifically, for gradient of phi we take each basis column and multiply it by lambda,
#         #evaluate the numerical integral of the product, and then divide off Lambda01
#         J <- cbind(
#           #first, this multiplies every row of basis1 by exp(phi1),
#           #then, this divides every column of the resulting matrix by the vector Lambda01
#           basis3,
#           matrix(x3new, nrow = length(T2), ncol = length(x3new), byrow = T))
#         #vector with as many elements as times in T2
#         Var.loglogS.3 <- J %*% J %*% obj$Finv[c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3])),
#                                               c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3]))] %*%
#           t(J)
#       } else if (is.null(x3new) | nP[3] == 0) {
#         J <- basis3
#         Var.loglogS.3 <- J %*% obj$Finv[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])] %*% t(J)
#       }
#
#       #Use delta method compute baseline hazard confidence intervals
#       # h(t) = s'(z)/t * exp(s(z) + xtbeta) where z=log(t)
#       h.3 <- as.vector(dbasis3 %*% phi3 / T2 * exp(basis3 %*% phi3 + eta3))
#       #so, each row of J is gradient wrt phi1, ..., phik, beta1, ..., betap
#       if(all(is.na(obj$Finv))){
#         Var.h.3 <- NA
#       } else if (!is.null(x3new) & nP[2] > 0) {
#         J <- cbind(dbasis3*Lambda03*exp(eta3)/T2 + h.3*basis3, h.3 *
#                      matrix(x3new, nrow = length(T2), ncol = length(x3new), byrow = T))
#         Var.h.3 <- J %*% J %*% obj$Finv[c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3])),
#                                         c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3]))] %*% t(J)
#       }
#       else if (is.null(x3new) | nP[3] == 0) {
#         J <- dbasis3*Lambda03*exp(eta3)/T2 + h.3*basis3
#         Var.h.3 <- J %*% obj$Finv[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])] %*% t(J)
#       }
#     } else if(obj$class[4] == "B-Spline"){
#       basis3 <- get_basis(x = T2,knots = obj$knots_list[[3]],hazard = "bspline")
#       basis3_quad <- get_basis(x = transform_quad_points(n_quad = obj$n_quad,
#                                                          quad_method = obj$quad_method,
#                                                          a = 0,b = T2),
#                                knots = obj$knots_list[[3]],hazard = "bspline")
#       phi3 <- obj$estimate[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])]
#       lambda03 <- as.vector(exp(basis3_quad %*% phi3))
#
#       #reshape lambda01 from a n*n_quad length vector
#       #to an n by n_quad matrix, then multiply with n_quad length weights
#       #to get final Lambda01
#       Lambda03 <- T2/2 * as.vector(matrix(lambda03,ncol=obj$n_quad,byrow = TRUE) %*% get_quad_pointsweights(n_quad=obj$n_quad)$weights)
#       S.3 <- exp(-Lambda03*exp(eta3))
#
#       #Use delta method to compute baseline survival confidence intervals
#       if(all(is.na(obj$Finv))){
#         Var.loglogS.3 <- NA
#       } else if (!is.null(x3new) & nP[3] > 0) {
#         #matrix with as many columns as parameters in S1, and as many rows as times in T2
#         #under piecewise, log(-log(S(t))) = log( basis1 %*% exp(phi1)) + xtbeta
#         #so, each row of J is the gradient at a particular T2 wrt phi1, ..., phiK, beta1, ..., betap
#         J_phi <- apply(X = basis3_quad * lambda03, MARGIN = 2,
#                        FUN = function(x) T2/2 * as.vector(matrix(x,ncol=obj$n_quad,byrow = TRUE) %*% get_quad_pointsweights(n_quad=obj$n_quad)$weights))
#         J <- cbind(
#           #first, this multiplies every row of basis1 by exp(phi1),
#           #then, this divides every column of the resulting matrix by the vector Lambda01
#           J_phi / Lambda03,
#           matrix(x3new, nrow = length(T2), ncol = length(x3new), byrow = T))
#         #vector with as many elements as times in T2
#         Var.loglogS.3 <- J %*% obj$Finv[c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3])),
#                                         c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3]))] %*%
#           t(J)
#       } else if (is.null(x3new) | nP[3] == 0) {
#         J_phi <- apply(X = basis3_quad * lambda03, MARGIN = 2,
#                        FUN = function(x) T2/2 * as.vector(matrix(x,ncol=obj$n_quad,byrow = TRUE) %*% get_quad_pointsweights(n_quad=obj$n_quad)$weights))
#         J <- J_phi / Lambda03
#         Var.loglogS.3 <- J %*% obj$Finv[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])] %*% t(J)
#       }
#
#       #Use delta method compute baseline hazard confidence intervals
#       h.3 <- as.vector(exp(basis3 %*% phi3) * exp(eta3))
#       #hazard is exp(Btphi) * exp(xtbeta)
#       #so, each row of J is gradient wrt phi1, ..., phik, beta1, ..., betap
#       if(all(is.na(obj$Finv))){
#         Var.h.3 <- NA
#       } else if (!is.null(x3new) & nP[3] > 0) {
#         J <- cbind(h.3 * basis3, h.3 *
#                      matrix(x3new, nrow = length(T2), ncol = length(x3new), byrow = T))
#         Var.h.3 <- J %*% obj$Finv[c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3])),
#                                   c((1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0_tot+nP[1]+nP[2]):(nP0_tot+nP[1]+nP[2]+nP[3]))]  %*% t(J)
#       }
#       else if (is.null(x3new) | nP[3] == 0) {
#         J <- h.3 * basis3
#         Var.h.3 <- J %*% obj$Finv[(1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3]), (1+nP0[1]+nP0[2]):(nP0[1]+nP0[2]+nP0[3])] %*% t(J)
#       }
#     } else{stop("class must be Weibull, Piecewise Constant, or B-Spline")}
#
#     if(all(is.na(Var.loglogS.3))){
#       se.loglogS.3 <- NA
#     } else{
#       se.loglogS.3 <- sqrt(diag(Var.loglogS.3))
#       se.loglogS.3[is.nan(se.loglogS.3)] <- 0
#     }
#     LL.3 <- S.3^exp(-stats::qnorm(conf.level/2) * se.loglogS.3)
#     UL.3 <- S.3^exp(stats::qnorm(conf.level/2) * se.loglogS.3)
#     if(all(is.na(Var.h.3))){
#       se.h.3 <- NA
#     } else{
#       se.h.3 <- sqrt(diag(Var.h.3))
#       se.h.3[is.nan(se.h.3)] <- 0
#     }
#     LLh.3 <- h.3 + stats::qnorm(conf.level/2) * se.h.3 #sign reversed because 0.025 quantile is negative
#     ULh.3 <- h.3 - stats::qnorm(conf.level/2) * se.h.3
#     LLh.3[LLh.3 < 0] <- 0
#
#     T2h <- T2
#     if (T2[1] == 0) {
#       T2h <- T2h[-1]
#       h.1 <- h.1[-1]
#       LLh.1 <- LLh.1[-1]
#       ULh.1 <- ULh.1[-1]
#       h.2 <- h.2[-1]
#       LLh.2 <- LLh.2[-1]
#       ULh.2 <- ULh.2[-1]
#       h.3 <- h.3[-1]
#       LLh.3 <- LLh.3[-1]
#       ULh.3 <- ULh.3[-1]
#     }
#     BH1_tbl <- data.frame(time = T2h, h.1 = h.1, LL.1 = LLh.1,
#                           UL.1 = ULh.1)
#     BH2_tbl <- data.frame(time = T2h, h.2 = h.2, LL.2 = LLh.2,
#                           UL.2 = ULh.2)
#     BH3_tbl <- data.frame(time = T2h, h.3 = h.3, LL.3 = LLh.3,
#                           UL.3 = ULh.3)
#     BS1_tbl <- data.frame(time = T2, S.1 = S.1, LL.1 = LL.1,
#                           UL.1 = UL.1)
#     BS2_tbl <- data.frame(time = T2, S.2 = S.2, LL.2 = LL.2,
#                           UL.2 = UL.2)
#     BS3_tbl <- data.frame(time = T2, S.3 = S.3, LL.3 = LL.3,
#                           UL.3 = UL.3)
#     value <- list(h.1 = BH1_tbl, h.2 = BH2_tbl, h.3 = BH3_tbl,
#                   S.1 = BS1_tbl, S.2 = BS2_tbl, S.3 = BS3_tbl)
#   }
#   value$xnew <- xnew
#   value$x1new <- x1new
#   value$x2new <- x2new
#   value$x3new <- x3new
#   value$tseq <- tseq
#   value$setup$model <- obj$setup$model
#   value$class <- obj$class
#   class(value) <- "pred.Freq_HReg2"
#   return(value)
# }

