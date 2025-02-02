% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblogreg.R
\name{blblogreg}
\alias{blblogreg}
\title{blblogreg}
\usage{
blblogreg(formula, data, m = 10, B = 5000, ncpu = 1)
}
\arguments{
\item{formula}{the formula of choice - logistic linear regression with weight in this case.}

\item{data}{data of choice}

\item{m}{a positive number, the number of items to choose from.}

\item{B}{B a posive integer. It represents a number of bootstrap samples.}

\item{ncpu}{Parallel is set as TRUE in default which means the algorithm is implementing multiple CPUs, 4. The user can select FALSE to implement only one CPU.}
}
\description{
blblogreg
}
