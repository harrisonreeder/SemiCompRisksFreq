% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal_likelihood.R
\name{nll_marg_profile_func}
\alias{nll_marg_profile_func}
\title{Negative Profile Log-Likelihood Function for Illness-Death Model}
\usage{
nll_marg_profile_func(
  fixed_param,
  fixed_param_ind,
  para_mle,
  y1,
  y2,
  delta1,
  delta2,
  yL,
  anyLT,
  Xmat1,
  Xmat2,
  Xmat3,
  weights,
  hazard,
  model,
  fraildist = "gamma",
  frailtyparam_num,
  n_quad,
  verbose = FALSE,
  control,
  hessian,
  optim_method,
  extra_starts
)
}
\arguments{
\item{fixed_param}{value of parameter being fixed (can be vector if fixed parameters are contiguous)}

\item{fixed_param_ind}{index of parameter after which fixed parameters would be inserted}

\item{y1, y2}{Numeric vectors of length \eqn{n} with (possibly censored) non-terminal and terminal event times}

\item{delta1, delta2}{Numeric vectors of length \eqn{n}  with indicators of 1 if the event was observed and 0 otherwise}

\item{yL}{Numeric vector of length \eqn{n} with left-truncation times.}

\item{anyLT}{Boolean (or 1/0) value that is \code{TRUE} (or 1) if there are any non-zero left truncation times, and \code{FALSE} (or 0) otherwise}

\item{Xmat1, Xmat2, Xmat3}{Numeric matrices with \eqn{n} rows and \eqn{q_1,q_2,q_3} columns containing covariates.}

\item{weights}{Vector of numeric weights of length \eqn{n} corresponding
to each observed individual.}

\item{hazard}{String (not case sensitive) specifying the form of the baseline hazard.
Options are \code{"weibull"} for Weibull,
\code{"piecewise"} for piecewise constant on the hazard scale,
\code{"bspline"} for cubic B-spline on the log-hazard scale,
and \code{"royston-parmar"} for restricted cubic spline on the log-cumulative hazard scale.
Aliases for these are \code{"wb"}, \code{"pw"}, \code{"bs"}, and \code{"rp"} respectively.}

\item{model}{String (not case sensitive) specifying the transition assumption: either \code{"semi-Markov"}}

\item{control}{a list of control attributes passed directly
into the \code{optim} function.}

\item{hessian}{Boolean indicating whether the hessian
(inverse of the covariance matrix) should be computed and returned.}

\item{optim_method}{a string naming which \code{optim}
method should be used for optimization.}

\item{extra_starts}{Integer giving the number of extra
starts to try when optimizing.}
}
\value{
Returns numeric sum of negative log likelihood contributions.
}
\description{
Function returning the negative profile log-likelihood for the illness-death model,
under specified baseline hazard, and specified frailty,
and specified Markov/semi-Markov transition assumption, fixing some parameters
Typically, this function will not be used directly by the user,
but as part of a larger estimation procedure.
}
