% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal_likelihood.R
\name{nll_profile_helper_grad_func}
\alias{nll_profile_helper_grad_func}
\title{Conditional Negative Log-Likelihood Gradient Function for Illness-Death Model}
\usage{
nll_profile_helper_grad_func(
  other_para,
  fixed_param,
  fixed_param_ind,
  y1,
  y2,
  delta1,
  delta2,
  yL,
  anyLT,
  Xmat1,
  Xmat2,
  Xmat3,
  hazard,
  frailty,
  model,
  weights,
  basis1,
  basis2,
  basis3,
  basis3_y1,
  basis1_yL,
  basis2_yL,
  dbasis1,
  dbasis2,
  dbasis3
)
}
\arguments{
\item{other_para}{vector of values for remaining parameters.}

\item{fixed_param}{value of parameter being fixed (can be vector if fixed parameters are contiguous)}

\item{fixed_param_ind}{index of parameter after which fixed parameters would be inserted}

\item{y1, y2}{Numeric vectors of length \eqn{n} with (possibly censored) non-terminal and terminal event times}

\item{delta1, delta2}{Numeric vectors of length \eqn{n}  with indicators of 1 if the event was observed and 0 otherwise}

\item{yL}{Numeric vector of length \eqn{n} with left-truncation times.}

\item{anyLT}{Boolean (or 1/0) value that is \code{TRUE} (or 1) if there are any non-zero left truncation times, and \code{FALSE} (or 0) otherwise}

\item{Xmat1, Xmat2, Xmat3}{Numeric matrices with \eqn{n} rows and \eqn{q_1,q_2,q_3} columns containing covariates.}

\item{hazard}{String (not case sensitive) specifying the form of the baseline hazard.
Options are \code{"weibull"} for Weibull,
\code{"piecewise"} for piecewise constant on the hazard scale,
\code{"bspline"} for cubic B-spline on the log-hazard scale,
and \code{"royston-parmar"} for restricted cubic spline on the log-cumulative hazard scale.
Aliases for these are \code{"wb"}, \code{"pw"}, \code{"bs"}, and \code{"rp"} respectively.}

\item{frailty}{Boolean indicating whether a gamma distributed subject-specific frailty should
be included. Currently this must be set to \code{TRUE}.}

\item{model}{String (not case sensitive) specifying the transition assumption: either \code{"semi-Markov"}}

\item{weights}{Vector of numeric weights of length \eqn{n} corresponding
to each observed individual.}

\item{basis1, basis2, basis3, basis3_y1, basis1_yL, basis2_yL}{Numeric matrices with \eqn{n} rows and \eqn{k_1,k_2,k_3} columns
with piecewise/spline basis function values at the corresponding \code{y1} and \code{y2} values.
Under semi-Markov model, basis3 represents basis derived from \eqn{y_2-y_1} and \code{basis3_y1} is unused,
while under Markov model, basis3 represents basis derived from \eqn{y_2} and \code{basis3_y1} is from \eqn{y_1}
Not used under Weibull model.}

\item{dbasis1, dbasis2, dbasis3}{Numeric matrices with \eqn{n} rows and \eqn{k_1,k_2,k_3} columns
with piecewise/spline basis function derivative values at the corresponding \code{y1} and \code{y2} values.
Used only under Royston-Parmar model.}
}
\value{
Returns numeric sum of negative log likelihood gradient contributions.
}
\description{
Function returning the gradient of the negative log-likelihood for the illness-death model,
under specified baseline hazard, and specified frailty,
and specified Markov/semi-Markov transition assumption, fixing some parameters.
Typically, this function will not be used directly by the user,
but as part of a larger estimation procedure.
}
