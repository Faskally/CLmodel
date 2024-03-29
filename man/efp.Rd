% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Likelihoods_and_gradients.R
\name{efp}
\alias{efp}
\title{Estimate capture probabilites from electrofishing data}
\usage{
efp(
  formula,
  data = NULL,
  passes = NULL,
  verbose = TRUE,
  init = "0",
  hessian = FALSE
)
}
\arguments{
\item{formula}{a formula object}

\item{data}{a data.frame containing all relavent info}

\item{passes}{todo}

\item{verbose}{show verbose output}

\item{init}{todo}

\item{hessian}{return hessian or not}
}
\value{
glm type object
}
\description{
This function uses the marginal likelihood of capture probabilities
to estimate model parameters
}
