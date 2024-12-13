% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/penalized_functions.R
\name{FreqID_HReg_Rpath}
\alias{FreqID_HReg_Rpath}
\title{Fit Penalized Parametric Frailty Illness-Death Model Solution Path}
\usage{
FreqID_HReg_Rpath(
  Formula,
  data,
  na.action = "na.fail",
  subset = NULL,
  weights = NULL,
  hazard = c("weibull"),
  frailty = TRUE,
  model,
  knots_list = NULL,
  penalty = c("scad", "mcp", "lasso"),
  lambda_path = NULL,
  lambda_target = 0,
  N_path_steps = 40,
  a = NULL,
  select_tol = 1e-04,
  fusion_tol = 0.001,
  penalty_fusedcoef = c("none", "fusedlasso"),
  lambda_fusedcoef_path = 0,
  penweights_list = list(),
  mu_smooth_path = 0,
  fit_method = "prox_grad",
  nP0 = rep(4, 3),
  startVals = NULL,
  quad_method = "kronrod",
  n_quad = 15,
  ball_L2 = Inf,
  warm_start = TRUE,
  step_size_min = 1e-06,
  step_size_max = 1e+06,
  step_size_init = 1,
  step_size_scale = 0.5,
  maxit = 300,
  extra_starts = 0,
  conv_crit = "nll_pen_change",
  conv_tol = 1e-06,
  standardize = TRUE,
  verbose = 0
)
}
\arguments{
\item{Formula}{a Formula object, with the outcome on the left of a
\code{~}, and covariates on the right. It is of the form, \code{time to non-terminal
  event + corresponding censoring indicator | time to terminal event
  + corresponding censoring indicator ~ covariates for h_1 |
  covariates for h_2 | covariates for h_3}. For example, \code{y_1 + delta_1 | y_2 + delta_2 ~ x_1 | x_2 | x_3}.}

\item{data}{a \code{data.frame} in which to interpret the variables named in Formula.}

\item{na.action}{how NAs are treated. See \code{\link[stats]{model.frame}}.}

\item{subset}{a specification of the rows to be used: defaults to all rows. See \code{\link[stats]{model.frame}}.}

\item{weights}{Vector of numeric weights of length \eqn{n} corresponding
to each observed individual.}

\item{hazard}{String (not case sensitive) specifying the form of the baseline hazard.
Options are \code{"weibull"} for Weibull,
\code{"piecewise"} for piecewise constant on the hazard scale,
\code{"bspline"} for cubic B-spline on the log-hazard scale,
and \code{"royston-parmar"} for restricted cubic spline on the log-cumulative hazard scale.
Aliases for these are \code{"wb"}, \code{"pw"}, \code{"bs"}, and \code{"rp"} respectively.}

\item{frailty}{Boolean indicating whether a gamma distributed subject-specific frailty should
be included. Currently this must be set to \code{TRUE}.}

\item{model}{String (not case sensitive) specifying the transition assumption: either \code{"semi-Markov"}}

\item{knots_list}{Used for hazard specifications besides Weibull, a
list of three increasing sequences of integers, each corresponding to
the knots for the flexible model on the corresponding transition baseline hazard. If
\code{NULL}, will be created by \code{\link{get_default_knots_list}}.}

\item{penalty}{A string value indicating the form of parameterwise penalty
to apply. "lasso", "scad", and "mcp" are the options.}

\item{lambda_path}{Numeric sequence of decreasing regularization parameters
for the parameterwise penalties, along which the solution path runs.
Assumes a single shared penalty across transitions.}

\item{lambda_target}{Final lambda value for parameterwise penalty. Ignored if lambda_path is specified.}

\item{N_path_steps}{number of steps generated for regularization path. Ignored if lambda_path is specified.}

\item{a}{For two-parameter penalty functions (e.g., scad and mcp), the second parameter.}

\item{select_tol}{Positive numeric value for thresholding estimates to be equal to zero.}

\item{fusion_tol}{Positive numeric value for thresholding estimates that are close
to being considered fused, for the purposes of estimating degrees of freedom.}

\item{penalty_fusedcoef}{A string value indicating the form of the fusion penalty to apply
to the regression parameters. "none" and "fusedlasso" are the options.}

\item{lambda_fusedcoef_path}{Numeric sequence of increasing regularization parameters
for the fusion penalties.}

\item{penweights_list}{A list of numeric vectors representing weights for each
penalty term (e.g., for adaptive lasso.) Elements of the list should be indexed by the
names "coef1", "coef2", "coef3", "fusedcoef12", "fusedcoef13", "fusedcoef23", "fusedbaseline12", "fusedbaseline13", and "fusedbaseline23"}

\item{mu_smooth_path}{Numeric sequence of decreasing Nesterov smoothing parameters for the fusion penalties.}

\item{fit_method}{String indicating which optimization method should be used at each step.}

\item{nP0}{vector of length three of integers indicating how many
baseline hazard parameters
should be specified for each of the three transition hazards.
This input is only relevant when
hazard is something other than "weibull" and is superceded by knots_list.}

\item{startVals}{A numeric vector of parameter starting values, arranged as follows:
the first \eqn{k_1+k_2+k_3} elements correspond to the baseline hazard parameters,
then the \eqn{k_1+k_2+k_3+1} element corresponds to the gamma frailty log-variance parameter,
then the last\eqn{q_1+q_2+q_3} elements correspond with the regression parameters.
If set to \code{NULL}, will be generated automatically using \code{\link{get_start}}.}

\item{quad_method}{String indicating which quadrature method to use to
evaluate numerical integral of B-spline. Options are
\code{"kronrod"} for Gauss-Kronrod quadrature or
\code{"legendre"} for Gauss-Legendre quadrature.}

\item{n_quad}{Scalar for number of Gaussian quadrature points used to
evaluate numerical integral of B-spline.}

\item{ball_L2}{Positive numeric value for \eqn{l_2} ball constraint around the origin for the regression parameters.
Typically set to \code{Inf} indicating no constraint, otherwise equivalent to an extra \eqn{l_2} penalty.}

\item{warm_start}{Boolean indicating whether each step of the solution should start from
ending of previous (\code{TRUE}) or from the original starting point (\code{FALSE}).}

\item{step_size_min}{Positive numeric value for the minimum allowable step size to allow during backtracking.}

\item{step_size_max}{Positive numeric value for the maximum allowable step size to allow by size increase at each iteration.}

\item{step_size_init}{Positive numeric value for the initial step size.}

\item{step_size_scale}{Positive numeric value for the multiplicative change in step size at each step of backtracking.}

\item{maxit}{Positive integer maximum number of iterations.}

\item{extra_starts}{numeric indicating how many additional optimization runs from random start values
should be performed at each grid point.}

\item{conv_crit}{String (possibly vector) giving the convergence criterion.}

\item{conv_tol}{Positive numeric value giving the convergence tolerance for the chosen criterion.}

\item{standardize}{Boolean for whether columns of covariate matrix should be standardized before estimation.
Estimates are then rescaled back to original scale when output.}

\item{verbose}{Numeric indicating the amount of iteration information should be printed to the user.
Higher numbers provide more detailed information to user, but will slow down the algorithm.}
}
\value{
A list.
}
\description{
Fit Penalized Parametric Frailty Illness-Death Model Solution Path
}
