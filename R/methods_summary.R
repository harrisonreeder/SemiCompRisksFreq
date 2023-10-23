#' Compute estimates of variance-covariance matrix
#'
#' This function computes different variance estimate from frequentist
#' model object, using the observed information and/or the summed outer products
#' of each individual subject's scores evaluated at the MLE.
#'
#' @param Finv inverse of observed information matrix
#' @param cheese summed score function outer products
#' @param n number of observations
#' @param var_type String giving type of variance covariance matrix estimate to return.
#'  Options are "modelbased" which returns the inverse of the observed information,
#'  "outer" which returns the inverse of the summed score function outer products,
#'  or "sandwich" which returns the sandwich estimator
#' @param name_vec String vector giving the desired names of each row/column of the matrix
#' @param df_adjust Boolean for whether to include "finite sample"
#'  degrees-of-freedom adjustment by multiplication with \eqn{n/(n-k)}
#'  where \eqn{n} is the number of observations and
#'  \eqn{k} is the number of estimated parameters.
#'
#' @export
vcov_helper <- function (Finv=NULL, cheese=NULL,
                         var_type=c("modelbased","sandwich","outer"),
                         name_vec=NULL,
                         df_adjust = FALSE, n){
  var_type <- match.arg(var_type)

  if(var_type %in% c("modelbased","sandwich") & is.null(Finv)){
    warning("'Finv' variance-covariance matrix not provided. Returning NULL.")
    return(NULL)
  }
  if(var_type %in% c("outer","sandwich") & is.null(cheese)){
    warning("fitted model does not have gradient outer-product 'cheese' matrix stored. Returning NULL")
    return(NULL)
  }

  val <- NULL
  if(var_type == "modelbased"){ val <- Finv }
  if(var_type == "sandwich"){ val <- Finv %*% cheese %*% Finv }
  if(var_type == "outer"){
    val <- tryCatch(MASS::ginv(cheese),
             error=function(cnd){message(cnd);cat("\n");return(NULL)})
  }

  if(!is.null(val) & !all(is.na(val))){
    if(df_adjust) val <- val * n / (n - NROW(val))
    rownames(val) <- colnames(val) <- name_vec
  }

  val
}


#Now, a bunch of standard methods defined for model fit object classes

#returns model-based estimate using inverse of observed information.
#' @export
vcov.Freq_HReg2 <- function (object, ...){
  if(is.null(object$Finv)){
    warning("object does not have 'Finv' variance-covariance matrix stored.
            Returning NULL.")
    return(NULL)
  }
  val <- object$Finv
  rownames(val) <- colnames(val) <- names(object$estimate)
  val
}

#logLik is also used by such methods as stats:::AIC.default
#note that for now, I'm defining degrees of freedom excluding exact 0's
#because in practice, this will be something I've "set" to zero for a reason

#' @export
logLik.Freq_HReg2 <- function (object, ...){
  val <- object$logLike; class(val) <- "logLik"
  attr(x = val,which = "df") <- sum(object$estimate != 0) #length(object$estimate)
  attr(x = val,which = "nobs") <- object$nobs
  val
}


#' @export
extractAIC.Freq_HReg2 <- function (fit, scale, k=2, ...){
  c(sum(fit$estimate != 0), #length(fit$estimate),
    stats::AIC(fit,k=k))
}

#' @export
coef.Freq_HReg2 <- function (object, ...){
  object$estimate
}

#' @export
print.Freq_HReg2 <- function (x, digits = 3, alpha = 0.05, ...)
{
  logEst <- x$estimate
  logSE <- if(all(is.na(x$Finv))) NA else sqrt(diag(x$Finv))
  value <- cbind(beta=logEst,SE=logSE,
                 LL=logEst - abs(stats::qnorm(alpha/2, 0, 1)) * logSE,
                 UL=logEst + abs(stats::qnorm(alpha/2, 0, 1)) * logSE,
                 z=logEst/logSE,
                 pvalue=stats::pchisq((logEst/logSE)^2,df = 1,lower.tail = FALSE))
  dimnames(value) <- list(x$myLabels, c("beta", "SE", "LL", "UL","z","pvalue"))
  if (x$class[2] == "Surv") {
    cat("\nAnalysis of independent univariate time-to-event data \n")
    cat(x$class[4], "baseline hazard specification\n")
    cat("Confidence level: ", round(1-alpha,3)*100, "%\n", sep = "")
    if (sum(x$nP) != 0) {
      cat("\nRegression coefficients:\n")
      print(round(value[-c(1:(sum(x$nP0))), ], digits = digits))
    }
  }
  if (x$class[2] == "ID") {
    cat("\nAnalysis of independent semi-competing risks data \n")
    cat(x$class[4], "baseline hazard specification\n")
    cat(x$class[5], "specification for h3\n")
    cat("Confidence level: ", round(1-alpha,3)*100, "%\n", sep = "")

    #print the frailty variance results
    if (x$frailty == TRUE){
      cat("\nVariance of frailties, theta:\n")
      theta_ests <- c(exp(value[sum(x$nP0)+1,1:4]))
      names(theta_ests) <- c("theta", "SE", "LL", "UL")
      #update theta SE using delta method
      #(accounting for the fact that we exponentiated everything just now)
      theta_ests["SE"] <- log(theta_ests["SE"]) * theta_ests["theta"]

      print(round(c(theta_ests,x$frailty_test[c("lrtest","lrpvalue")]),
                  digits = digits))
      cat("SE computed from SE(log(theta)) via delta method.")
      cat("\nBounds formed for log(theta) and exponentiated.")
      cat("\nLikelihood ratio test of theta=0 vs. theta>0 using mixture of chi-squareds null.\n")
    } else{
      #no frailty variance estimate, so just say so.
      cat("\nNon-frailty model fit\n")
    }
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
summary.Freq_HReg2 <- function (object, alpha = 0.05,
                                var_type=c("modelbased","sandwich","outer"), ...) {
  # browser()
  var_type <- match.arg(var_type)
  logEst <- object$estimate
  temp_vcov <- vcov_helper(Finv = object$Finv, cheese = object$cheese,
                           n = object$nobs, name_vec = names(object$estimate),
                           var_type = var_type, df_adjust = FALSE)
  logSE <- if(!is.null(temp_vcov) & !all(is.na(temp_vcov))) sqrt(diag(temp_vcov)) else NA
  results <- cbind(beta=logEst,SE=logSE,
                LL=logEst - abs(stats::qnorm(alpha/2, 0, 1)) * logSE,
                UL=logEst + abs(stats::qnorm(alpha/2, 0, 1)) * logSE,
                z=logEst/logSE,
                pvalue=stats::pchisq((logEst/logSE)^2,df = 1,lower.tail = FALSE))
  if (object$class[2] == "Surv") {
    output.coef <- results[-c(1:object$nP0),,drop=FALSE]
    dimnames(output.coef) <- list(unique(object$myLabels[-c(1:object$nP0)]),
                                  c("beta", "SE", "LL", "UL","z","pvalue"))
    output.HR <- exp(output.coef[,c("beta","LL","UL"),drop=FALSE])
    colnames(output.HR) <- c("exp(beta)","LL","UL")

    output.h0 <- results[1:object$nP0,1:4,drop=FALSE]
    if(object$class[4]=="Weibull"){
      dimnames(output.h0) <- list(c("Weibull: log-kappa", "Weibull: log-alpha"),
                                  c("h-PM", "SE", "LL", "UL"))
      knots_mat <- NULL
    } else{
      dimnames(output.h0) <- list(paste0(object$class[4],": phi",1:object$nP0),
                               c("h-PM", "SE", "LL", "UL"))
      knots_mat <- as.matrix(object$knots_vec)
      rownames(knots_mat) <- paste0("knot",1:length(object$knots_vec))
    }

    value <- list(HR = output.HR,
                  coef = output.coef, h0 = output.h0,
                  logLike = object$logLike,
                  nP = nrow(results), class = object$class,
                  alpha=alpha, #conf.level = conf.level,
                  knots_mat=knots_mat, var_type=var_type)
  }
  if (object$class[2] == "ID") {
    nP.0 <- ifelse(object$frailty, sum(object$nP0)+1, sum(object$nP0))
    nP.1 <- object$nP[1]; nP.2 <- object$nP[2]; nP.3 <- object$nP[3]
    beta.names <- unique(object$myLabels[-c(1:nP.0)])
    nP <- length(beta.names)
    output.coef <- matrix(NA, nrow = nP, ncol = 18)
    dimnames(output.coef) <- list(beta.names,
        c("beta1", "beta1_SE", "beta1_LL", "beta1_UL","beta1_z","beta1_pvalue",
          "beta2", "beta2_SE", "beta2_LL", "beta2_UL","beta2_z","beta2_pvalue",
          "beta3", "beta3_SE", "beta3_LL", "beta3_UL","beta3_z","beta3_pvalue"))
    for (i in 1:nP) {
      if (nP.1 != 0) {
        for (j in 1:nP.1) if (object$myLabels[nP.0+j] == beta.names[i])
          output.coef[i,1:6] <- results[nP.0+j,]
      }
      if (nP.2 != 0) {
        for (j in 1:nP.2) if (object$myLabels[nP.0+nP.1+j] == beta.names[i])
          output.coef[i,7:12] <- results[nP.0+nP.1+j,]
      }
      if (nP.3 != 0) {
        for (j in 1:nP.3) if (object$myLabels[nP.0+nP.1+nP.2+j] == beta.names[i])
          output.coef[i,13:18] <- results[nP.0+nP.1+nP.2+j,]
      }
    }
    output.HR <- exp(output.coef[,c("beta1","beta1_LL","beta1_UL",
                                    "beta2","beta2_LL","beta2_UL",
                                    "beta3","beta3_LL","beta3_UL"),drop=FALSE])
    colnames(output.HR) <- c("exp(beta1)", "LL", "UL",
                             "exp(beta2)", "LL", "UL",
                             "exp(beta3)", "LL", "UL")

    output.theta <- output.ltheta <- rep(NA,6)
    if(object$frailty){
      output.ltheta <- c(results[nP.0,1:4],
                         object$frailty_test[c("lrtest","lrpvalue")])
      output.theta <- c(exp(results[nP.0,1:4]),
                        object$frailty_test[c("lrtest","lrpvalue")])
      #replace with standard error of theta from delta method
      output.theta[2] <- results[nP.0,2] * exp(results[nP.0,1])
    }
    names(output.ltheta) <- c("log(theta)", "SE", "LL", "UL", "lrtest", "lrpvalue")
    names(output.theta) <- c("theta", "SE", "LL", "UL", "lrtest", "lrpvalue")

    knots_mat <- NULL
    if(object$class[4]=="Weibull"){
      output.h0 <- matrix(NA, nrow = 2, ncol = 12,
                    dimnames=list(c("Weibull: log-kappa", "Weibull: log-alpha"),
                               c("h1-PM", "SE", "LL", "UL",
                                 "h2-PM", "SE", "LL", "UL",
                                 "h3-PM", "SE", "LL", "UL")))
      output.h0[1, 1:4] <- results[1,1:4]
      output.h0[1, 5:8] <- results[3,1:4]
      output.h0[1, 9:12] <- results[5,1:4]
      output.h0[2, 1:4] <- results[2,1:4]
      output.h0[2, 5:8] <- results[4,1:4]
      output.h0[2, 9:12] <- results[6,1:4]
    } else{ #this covers piecewise and spline models
      p01 <- object$nP0[1]; p02 <- object$nP0[2]; p03 <- object$nP0[3]
      p0max <- max(object$nP0)

      #generate "wide" matrix of baseline parameters by padding with 0s so all are same height
      output.h0 <- cbind(
        rbind(results[1:p01,1:4],matrix(data=0,ncol=4,nrow=(p0max-p01))),
        rbind(results[(1+p01):(p01+p02),1:4],matrix(data=0,ncol=4,nrow=(p0max-p02))),
        rbind(results[(1+p01+p02):(p01+p02+p03),1:4],matrix(data=0,ncol=4,nrow=(p0max-p03)))
      )
      dimnames(output.h0) <- list(paste0(object$class[4],": phi",1:p0max),
                               c("h1-PM", "SE", "LL", "UL",
                                 "h2-PM", "SE", "LL", "UL",
                                 "h3-PM", "SE", "LL", "UL"))

      #lastly, make a matrix with the knot locations, padded with NAs
      knotmax <- max(sapply(object$knots_list,length))
      knots_mat <- sapply(object$knots_list,FUN = function(x) c(x,rep(NA,knotmax-length(x))))
      dimnames(knots_mat) <- list(paste0("knot",1:knotmax), c("h1","h2","h3"))
    }

    value <- list(HR = output.HR,
                  coef = output.coef,
                  coef_long = results[-(1:nP.0),,drop=FALSE],
                  theta = output.theta,
                  ltheta = output.ltheta,
                  h0 = output.h0,
                  h0_long = results[1:sum(object$nP0),1:4],
                  logLike = object$logLike,
                  nP = nrow(results), class = object$class,
                  frailty = object$frailty,
                  alpha = alpha, #conf.level = conf.level,
                  knots_mat = knots_mat, var_type=var_type)
  }
  class(value) <- "summ.Freq_HReg2"
  return(value)
}

#' @export
print.summ.Freq_HReg2 <- function (x, digits = 3, ...)
{

  if (x$class[2] == "Surv") {
    cat("\nAnalysis of independent univariate time-to-event data \n")
  }
  if (x$class[2] == "ID") {
    cat("\nAnalysis of independent semi-competing risks data \n")
    cat(x$class[4], "baseline hazard specification\n")
    cat(x$class[5], "specification for h3\n")
  }
  cat("Confidence level: ", round(1-x$alpha,3)*100, "%\n", sep = "")

  if(x$var_type != "modelbased"){
    if(x$var_type=="sandwich"){
      cat("Variances estimated by robust sandwich estimator.")
    } else if(x$var_type=="outer"){
      cat("Variances estimated by outer product of gradients.")
    }
  }

  if (x$class[2] == "ID") {
    if(x$frailty){
      cat("\nVariance of frailties:\n")
      print(round(x$theta, digits = digits))
      cat("SE computed from log(theta) via delta method. Bounds exponentiated from log(theta).")
      cat("\nLikelihood ratio test of theta=0 vs. theta>0 using mixture of chi-squareds null.\n")
    } else{
      #no frailty variance estimate, so just say so.
      cat("\nNon-frailty model fit\n")
    }
  }

  if (!is.null(x$coef)) {
    cat("\nHazard ratios:\n")
    print(round(x$HR, digits = digits))
  }

  cat("\nBaseline hazard function components:\n")
  print(round(x$h0, digits = digits))
  if (x$class[4] != "Weibull") {
    cat("\nKnots:\n")
    print(round(x$knots_mat, digits = digits))
  }
  invisible()
}
