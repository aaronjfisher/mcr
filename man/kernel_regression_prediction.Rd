% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linear_regression_mcr.r
\name{kernel_regression_prediction}
\alias{kernel_regression_prediction}
\title{Get prediction from a (Nadaraya–Watson) kernel regression model}
\usage{
kernel_regression_prediction(X, X_ref, y_ref, kern_fun)
}
\arguments{
\item{X}{covariate matrix for which predictions are desired}

\item{X_ref}{covariate matrix for existing data}

\item{y_ref}{corresponding output vector for existing data}

\item{kern_fun}{kernel function to use}
}
\description{
Get prediction from a (Nadaraya–Watson) kernel regression model
}
