#' @export
ppwexp <- function(t, phi, knots_vec, eta=0, lower.tail = TRUE, log.p = FALSE){
  #first correct any simple flubs in the knots vector
  if(knots_vec[1] != 0){knots_vec <- c(0,knots_vec)} #starts with 0
  if(knots_vec[length(knots_vec)]==Inf){knots_vec <- knots_vec[-length(knots_vec)]}
  stopifnot(all(diff(knots_vec)> 0)) #increasing
  stopifnot(length(knots_vec)==length(phi))

  basis <- get_basis(y=t, knots_vec = knots_vec,hazard = "pw",deriv = FALSE)
  Lambda <- as.vector(basis %*% exp(phi)) * exp(eta)
  if(lower.tail){
    if(log.p){
      log1p(-exp(Lambda))
    } else{
      1 - exp(-Lambda)
    }
  } else{
    if(log.p){
      -Lambda
    } else{
      exp(-Lambda)
    }
  }
}

#' @export
invHaz_pwexp <- function(t, phi, knots_vec, eta=numeric(length(t))){
  # browser()
  #first correct any simple flubs in the knots vector
  if(knots_vec[1] != 0){knots_vec <- c(0,knots_vec)} #starts with 0
  if(knots_vec[length(knots_vec)]==Inf){knots_vec <- knots_vec[-length(knots_vec)]}
  stopifnot(all(diff(knots_vec)> 0)) #increasing
  stopifnot(length(knots_vec)==length(phi))
  if(length(t) > 1 & length(eta)==1){eta <- rep(eta,length(t))}

  #slope of first interval from 0=knots_vec[1] to knots_vec[2]) is exp(phi[1])
  #so inverse is line on interval from
  #0=knots_vec[1] to knots_vec[2]*exp(phi[1]) with slope exp(-phi[1])

  #slope of second interval (from knots_vec[2] to knots_vec[3]) is exp(phi[2])
  #so inverse is line on interval from
  #knots_vec[2]*exp(phi[1]) to
  #knots_vec[2]*exp(phi[1]) + (knots_vec[3]-knots_vec[2]) * exp(phi[2])
  #with slope exp(-phi[2])

  #and so on, with last interval going on until the end. So, let's make set of
  #"inverted" knots based on the baseline. (These will be further shifted by covariates)
  temp_knots <- c(0, cumsum(diff(knots_vec) * exp(phi[-length(phi)])))
  temp_coef <- exp(-phi)

  out <- sapply(1:length(t),
                function(i) {
                  # temp_x <- -log(p[i])
                  basis_i <- get_basis(y = t[i], knots_vec = temp_knots * exp(eta[i]),
                                       hazard="pw",deriv=FALSE)
                  Lambda0i <- basis_i %*% temp_coef
                  Lambda0i
                })
  as.vector(out) * exp(-eta)
}

#' @export
qpwexp <- function(p, phi, knots_vec, eta, lower.tail = TRUE, log.p = FALSE){
  temp_t <- if(lower.tail) -log1p(-p) else -log(p)
  out <- invHaz_pwexp(t = temp_t,phi = phi,knots_vec=knots_vec,eta = eta)
  if(log.p) log(out) else out
}
