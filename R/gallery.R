#' @import pracma
#' @import Matrix

#' @name binomial_matrix
#' @title Create binomial matrix
#'
#' @description Binomial matrix: an N-by-N matrix with
#' integer entries such that A^2 = 2^(N-1)*EYE(N)
#' Thus B = A*2^((1-N)/2) is involutory, that is B^2 = EYE(N)
#'
#' @param n - row dimension
#'
#' @export
binomial_matrix <- function(n){
  L <- abs(pracma::pascal(n,1))
  U <- L[n:1,n:1]
  D <- diag((-2)^(0:(n-1)))
  L %*% D %*% U
}

#' @name cauchy_matrix
#' @title Create cauchy matrix
#'
#' @description
#'
#' @param x vector of length n
#' @param y vector of length n
#'
#' @export
cauchy_matrix <- function(x,y=NULL){
  n <- length(x)
  if(n==1){
    n <- x
    x <- 1:n
  }

  if(is.null(y)){
    y <- x
  }

  if(length(x) != length(y)){
    stop("cauchy:ParamLengthMismatch")
  }

  1 / (matrix(x,nrow=n,ncol=n) + matrix(y,nrow=n,ncol=n,byrow=T))
}

#' @name spdiags
#' @title Create sparse diagonal matrix
#'
#' @description Creates a sparse representation of multipel diagonal matrix
#'
#' @param A - matrix where columns correspond to the desired diagonals
#' @param d - indices of the diagonals to be filled in. 0 is main diagonal. -1
#' is first subdiagonal and +1 is first superdiagonal.
#' @param m - row dim
#' @param n - col dim
#'
#' @return dgcMatrix sparse diagonal
#'
#' @export
spdiags <- function(A, d, m, n){
  rows <- unlist(sapply(d, function(d){
    if(d == 0){
      i <- 1:m
    } else if(d > 0){
      i <- 1:(m-d)
    } else if(d < 0){
      i <- (1-d):m
    }
    i
  }))

  cols <- unlist(sapply(d, function(d){
    if(d == 0){
      j <- 1:n
    } else if(d > 0){
      j <- (1+d):n
    } else if(d < 0){
      j <- 1:(n+d)
    }
    j
  }))

  B <- Matrix::sparseMatrix(i=rows, j=cols, x = as.vector(A),dims = c(m,n))
  B
}

#' @description Create a sparse tridiagonal matrix
#'
#' @param n - dimension of the square matrix
#' @param x - subdiagonal
#' @param y - diagonal
#' @param z - superdiagonal
#'
#' @return Sparse tridiagonal matrix
tridiag <- function(n, x=NULL, y=NULL, z=NULL){


  if(is.null(x)){
    x <- -1; y <- 2; z <- 1
  }
  if(is.null(z)){
    z <- y; y <- x; x <- n
  }

  if(max(c(length(x), length(y), length(z))) == 1){
    x <- x * rep(1,n-1)
    z <- z * rep(1,n-1)
    y <- y * rep(1,n)
  } else{
    nx <- length(x)
    ny <- length(y)
    nz <- length(z)
    if((ny - nx - 1) | (ny - nz - 1)){
      stop("tridiag:InvalidVectorArgDim")
    }
  }
  n <- length(y)
  # Create sparse representation of tridiagonal. Below is MATLAB notation
  # T = spdiags([ [x;0] y [0;z] ], -1:1, n, n)
  spdiags(matrix(c(x,0, y, 0,z),nrow=n), -1:1, n, n)
}
