% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gallery.R
\name{lauchli}
\alias{lauchli}
\title{Create Lauchli Matrix}
\usage{
lauchli(n, mu = NULL, sparse = F)
}
\arguments{
\item{n}{number of columns}

\item{mu}{constant applied to identity}

\item{sparse}{whether matrix should be sparse}
}
\value{
Lauchli matrix.
}
\description{
the (N + 1) x (N) matrix [ones(1,n); mu*eye(n)]. Well-known example in least squares of the
  danger of forming t(A) %*% A (due to inexact arithmetic, gives singular matrix)
}
