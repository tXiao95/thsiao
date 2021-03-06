% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gallery.R
\name{dramadah}
\alias{dramadah}
\title{Create anti-Hadamard matrix \code{A}}
\usage{
dramadah(n, k = 1)
}
\arguments{
\item{n}{order of matrix}

\item{k}{decides type of matrix returned.}
}
\description{
Returns a \code{n}-by-\code{n} nonsingular matrix of 0's and 1's.
  With large determinant or inverse. If \code{k=1}, \code{A} is Toeplitz and \code{abs(det(A))=1}.
  If \code{k=2}, \code{A} is upper triangular and Toeplitz. If \code{k=3}, \code{A} has maximal
  determinant among (0,1) lower Hessenberg matrices. Also is Toeplitz.

  Also known as an anti-Hadamard matrix.
}
