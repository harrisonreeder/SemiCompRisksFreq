% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_helper_functions.R
\name{transform_quad_points}
\alias{transform_quad_points}
\title{Return transformed Gaussian quadrature points based on the specified integral limits}
\usage{
transform_quad_points(n_quad, quad_method, a, b)
}
\arguments{
\item{n_quad}{Scalar for number of Gaussian quadrature points used to
evaluate numerical integral of B-spline.}

\item{quad_method}{String indicating which quadrature method to use to
evaluate numerical integral of B-spline. Options are
\code{"kronrod"} for Gauss-Kronrod quadrature or
\code{"legendre"} for Gauss-Legendre quadrature.}

\item{a}{vector lower limit of the integral}

\item{b}{vector upper limit of the integral}
}
\description{
Return transformed Gaussian quadrature points based on the specified integral limits
}
