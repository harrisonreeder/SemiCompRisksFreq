#' A simulated clustered semi-competing risks data set
#'
#' Simulated semi-competing risks data
#'
#' @format ## `scrData`
#' A data frame with 2000 rows and 8 columns:
#' \describe{
#'   \item{time1}{the time to non-terminal event}
#'   \item{event1}{the censoring indicators for the
#'   non-terminal event time; 1=event observed, 0=censored/truncated}
#'   \item{time2}{the time to terminal event}
#'   \item{event2}{the censoring indicators for the
#'   terminal event time; 1=event observed, 0=censored}
#'   \item{cluster}{cluster numbers}
#'   \item{x1,x2,x3}{vectors of continuous covarates}
#' }
#'
"scrData"
