% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linear_regression_mcr.r
\name{fit_lm_regularized}
\alias{fit_lm_regularized}
\title{Can be used to fit either kernel models or regularized linear models}
\usage{
fit_lm_regularized(suff_stats, alpha = NA, reg_threshold = Inf,
  tol = NA)
}
\arguments{
\item{suff_stats}{output from \code{\link{get_suff_stats_kernel}}, or from the internal function \code{get_suff_stats_lm}. See also \code{\link{precompute_mcr_objects_and_functions}}}

\item{alpha}{regularization parameter}

\item{reg_threshold}{regularization threshold}

\item{tol}{passed to \code{\link{solve_QP1QC}}}
}
\value{
The argmin over values of \eqn{\beta}.
}
\description{
Minimizes \eqn{||Y_i - X \beta ||^2 +}\code{alpha}\eqn{\beta'(}\code{reg_matrix}\eqn{)\beta} across \eqn{\beta}. Alternatively, if \code{reg_threshold} is finite,  this function minimizes \eqn{||Y_i - X \beta ||^2} across \eqn{\beta} satisfying \eqn{\beta'(}\code{reg_matrix}\eqn{)\beta\le}\code{reg_threshold}.
}
