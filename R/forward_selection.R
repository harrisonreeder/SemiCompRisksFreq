#' Forward Selection Procedure for Parametric Illness-Death Model
#'
#' This method is meant as a comparator for the penalized estimation procedure.
#'
#' @inheritParams FreqID_HReg2
#' @param vars character string with all of the possible variable names to be searched through.
#' @param steps the maximum number of steps to be considered. The default is 1000 (essentially as many as required). It is typically used to stop the process early.
#' @param k the multiple of the number of degrees of freedom used for the penalty. Only k = 2 gives the genuine AIC: k = log(n) is sometimes referred to as BIC or SBC.
#' @param verbose if positive, information is printed during the running of step. Larger values may give more detailed information.
#' @param fixed1 a string vector for parameters you want guaranteed included in arm 1
#' @param fixed2 a string vector for parameters you want guaranteed included in arm 2
#' @param fixed3 a string vector for parameters you want guaranteed included in arm 3
#'
#' @return A list.
#' @export
FreqID_HReg2_forward <- function(time1_name, time2_name, event1_name, event2_name, trunc_name=NULL,
                                 vars,
                                 data, na.action="na.fail", subset=NULL, weights=NULL,
                               hazard=c("weibull"), frailty=TRUE, model, knots_list=NULL,
                               fixed1=character(0),fixed2=character(0),fixed3=character(0),
                               k = 2, steps=1000, control=NULL, verbose = 1,
                               quad_method="kronrod", n_quad=15,
                               optim_method="BFGS", extra_starts=0, output_options=NULL){


  out_options_loop=list(fit_nf=FALSE, Finv=FALSE, grad_mat=FALSE,
                   cheese=FALSE, eb_frailties=FALSE,
                   # beta_test=TRUE, #indicator for whether score and lr tests of each beta should also be performed, not currently used
                   data=FALSE)

  n <- nrow(data)

  #empty vectors/lists used to track the strings that have been 'used' aka are active in the model.
  used1 <- used2 <- used3 <- character(0)
  used1_list <- used2_list <- used3_list <- list()
  used1_vec <- used2_vec <- used3_vec <- character(0)

  #monitors of the search process.
  modelcount <- 1
  best_crit_round <- Inf
  continue_flag <- TRUE
  model_crit_vec <- numeric(0)
  counter <- 1

  if(is.null(trunc_name)) {
    lefthand_string <- paste0(time1_name," + ",event1_name," | ",time2_name," + ",event2_name)
  } else{
    lefthand_string <- paste0(trunc_name," | ",time1_name," + ",event1_name," | ",time2_name," + ",event2_name)
  }

  #loop that iterates the forward selection process
  while(continue_flag & counter <= steps){
    best_crit_round <- Inf
    if(verbose>=1){
      print(counter)
    }

    ##loop that iterates through the covariates
    for(currvar in unique(vars)) {
      if(verbose >=3){
        print(currvar)
      }

      run1 <- run2 <- run3 <- FALSE

      #During earlier project, had more complicated logic here to accomodate
      #interactions and squared terms, but for now things are simplified
      if(!(currvar %in% c(fixed1,used1))){
        run1 <- TRUE
      }
      if(!(currvar %in% c(fixed2,used2))){
        run2 <- TRUE
      }
      if(!(currvar %in% c(fixed3,used3))){
        run3 <- TRUE
      }

      # print(paste0("run1: ",run1,", run2: ",run2,", run3: ",run3))
      if(run1){
        form <- Formula::as.Formula(paste0(lefthand_string," ~ ",
                                           paste(c(1,fixed1,used1,currvar),collapse = "+")," |",
                                           paste(c(1,fixed2,used2),collapse = "+")," |",
                                           paste(c(1,fixed3,used3),collapse = "+")," "))
        fit_temp <- FreqID_HReg2(Formula = form,
                     data=data, na.action=na.action, subset=subset, weights=weights,
                     hazard=hazard, frailty=frailty, model=model, knots_list=knots_list,
                     control=control,
                     quad_method=quad_method, n_quad=n_quad,
                     optim_method=optim_method, extra_starts=extra_starts, output_options=out_options_loop)

        #rewrite to allow different criteria. Maybe make a specific function that computes different ones?
        model_crit_vec[modelcount] <- extractAIC.Freq_HReg2(fit=fit_temp, k=k)[2]
        if(model_crit_vec[modelcount] < best_crit_round) {best_crit_round <- model_crit_vec[modelcount]}
        used1_list[[modelcount]] <- c(used1,currvar)
        used2_list[[modelcount]] <- c(used2)
        used3_list[[modelcount]] <- c(used3)
        used1_vec[modelcount] <- paste0(c(used1,currvar),collapse=", ")
        used2_vec[modelcount] <- paste0(c(used2),collapse=", ")
        used3_vec[modelcount] <- paste0(c(used3),collapse=", ")

        if(verbose>= 4){
          print(paste0("used1: ",used1_vec[modelcount]))
          print(paste0("used2: ",used2_vec[modelcount]))
          print(paste0("used3: ",used3_vec[modelcount]))
        }
        modelcount = modelcount + 1
      }

      if(run2){
        form <- Formula::as.Formula(paste0(lefthand_string," ~ ",
                                  paste(c(1,fixed1,used1),collapse = "+")," |",
                                  paste(c(1,fixed2,used2,currvar),collapse = "+")," |",
                                  paste(c(1,fixed3,used3),collapse = "+")," "))
        fit_temp <- FreqID_HReg2(Formula = form,
           data=data, na.action=na.action, subset=subset, weights=weights,
           hazard=hazard, frailty=frailty, model=model, knots_list=knots_list,
           control=control,
           quad_method=quad_method, n_quad=n_quad,
           optim_method=optim_method, extra_starts=extra_starts, output_options=out_options_loop)
        #rewrite to allow different criteria. Maybe make a specific function that computes different ones?
        model_crit_vec[modelcount] <- extractAIC.Freq_HReg2(fit=fit_temp, k=k)[2]
        if(model_crit_vec[modelcount] < best_crit_round) {best_crit_round <- model_crit_vec[modelcount]}
        used1_list[[modelcount]] <- c(used1)
        used2_list[[modelcount]] <- c(used2,currvar)
        used3_list[[modelcount]] <- c(used3)
        used1_vec[modelcount] <- paste0(c(used1),collapse=", ")
        used2_vec[modelcount] <- paste0(c(used2,currvar),collapse=", ")
        used3_vec[modelcount] <- paste0(c(used3),collapse=", ")

        if(verbose>= 4){
          print(paste0("used1: ",used1_vec[modelcount]))
          print(paste0("used2: ",used2_vec[modelcount]))
          print(paste0("used3: ",used3_vec[modelcount]))
        }
        modelcount = modelcount + 1
      }

      if(run3){
        form <- Formula::as.Formula(paste0(lefthand_string," ~ ",
                                  paste(c(1,fixed1,used1),collapse = "+")," |",
                                  paste(c(1,fixed2,used2),collapse = "+")," |",
                                  paste(c(1,fixed3,used3,currvar),collapse = "+")," "))
        fit_temp <- FreqID_HReg2(Formula = form,
           data=data, na.action=na.action, subset=subset, weights=weights,
           hazard=hazard, frailty=frailty, model=model, knots_list=knots_list,
           control=control,
           quad_method=quad_method, n_quad=n_quad,
           optim_method=optim_method, extra_starts=extra_starts, output_options=out_options_loop)
        model_crit_vec[modelcount] <- extractAIC.Freq_HReg2(fit=fit_temp, k=k)[2]
        if(model_crit_vec[modelcount] < best_crit_round) {best_crit_round <- model_crit_vec[modelcount]}
        used1_list[[modelcount]] <- c(used1)
        used2_list[[modelcount]] <- c(used2)
        used3_list[[modelcount]] <- c(used3,currvar)
        used1_vec[modelcount] <- paste0(c(used1),collapse=", ")
        used2_vec[modelcount] <- paste0(c(used2),collapse=", ")
        used3_vec[modelcount] <- paste0(c(used3,currvar),collapse=", ")

        if(verbose>= 4){
          print(paste0("used1: ",used1_vec[modelcount]))
          print(paste0("used2: ",used2_vec[modelcount]))
          print(paste0("used3: ",used3_vec[modelcount]))
        }
        modelcount = modelcount + 1
      }
    }


    #tail used because if there's a tie for some reason, just pick the last
    best_model <- utils::tail(which(min(model_crit_vec) == model_crit_vec),n=1)
    best_crit <- model_crit_vec[best_model]
    used1 <- used1_list[[best_model]]
    used2 <- used2_list[[best_model]]
    used3 <- used3_list[[best_model]]

    if(verbose>= 2){
      print(paste("best model at step", counter,"is:"))
      print(paste0("h1: "))
      print(c(fixed1,used1))
      print(paste0("h2: "))
      print(c(fixed2,used2))
      print(paste0("h3: "))
      print(c(fixed3,used3))
      print(paste0("with criterion value: ",best_crit))
    }

    if(best_crit < best_crit_round){
      continue_flag <- FALSE
    }
    counter <- counter + 1

  }

  if(verbose >= 3){
    print(paste0("done! Best model is (by hazard): "))
    print(paste(c(used1,fixed1),collapse = ", "))
    print(paste(c(used2,fixed2),collapse = ", "))
    print(paste(c(used3,fixed3),collapse = ", "))
  }

  crit_tab <- data.frame(model_crit_vec,used1_vec,used2_vec,used3_vec)


  out_options=list(fit_nf=TRUE, Finv=TRUE, grad_mat=FALSE,
                   cheese=TRUE, eb_frailties=TRUE,
                   # beta_test=TRUE, #indicator for whether score and lr tests of each beta should also be performed, not currently used
                   data=FALSE)
  nmsO <- names(out_options)
  namO <- names(output_options)
  out_options[namO] <- output_options

  #FIT BEST MODEL
  form <- Formula::as.Formula(paste0(lefthand_string," ~ ",
                                     paste(c(1,used1,fixed1),collapse = "+")," |",
                                     paste(c(1,used2,fixed2),collapse = "+")," |",
                                     paste(c(1,used3,fixed3),collapse = "+")," "))
  fit_temp <- FreqID_HReg2(Formula = form,
                            data, na.action=na.action, subset=subset,
                            hazard=hazard,frailty=frailty,
                            model=model, knots_list=knots_list,
                            optim_method = optim_method)

  return(list(final_fit=fit_temp,
              crit_tab=crit_tab,
              fixed1=fixed1,fixed2=fixed2,fixed3=fixed3,
              used1=used1,used2=used2,used3=used3,
              best_crit=best_crit,
              used1_list=used1_list,used2_list=used2_list,used3_list=used3_list,
              model_crit_vec=model_crit_vec,
              k=k))
}

