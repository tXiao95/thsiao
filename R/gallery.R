#' @importFrom pracma pascal
#' @importFrom Matrix sparseMatrix norm

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
  L <- abs(pascal(n, 1))
  U <- L[n:1,n:1]
  D <- diag((-2)^(0:(n-1)))
  L %*% D %*% U
}

#' @name cauchy_matrix
#' @title Create cauchy matrix
#'
#' @description Arguments \code{x} and \code{y} are vectors of length \code{n}. \code{C[i,j] = 1 / (x[i] + y[j])}
#'
#' @param x vector of length n
#' @param y vector of length n
#'
#' @export
cauchy_matrix <- function(x,y=NULL){
  n <- length(x)
  if(n == 1){
    n <- x
    x <- 1:n
  }

  if(is.null(y)){
    y <- x
  }

  if(length(x) != length(y)){
    stop("cauchy:ParamLengthMismatch")
  }

  1 / (matrix(x, nrow = n, ncol = n) + matrix(y, nrow = n,ncol = n, byrow = T))
}

#' @name spdiags
#' @title Create sparse diagonal matrix
#'
#' @description Creates a sparse representation of multiple diagonal matrix
#'
#' @param A matrix where columns correspond to the desired diagonals
#' @param d indices of the diagonals to be filled in. 0 is main diagonal. -1
#' is first subdiagonal and +1 is first superdiagonal.
#' @param m row dim
#' @param n col dim
#'
#' @return dgcMatrix sparse diagonal
#'
#' @export
spdiags <- function(A, d, m, n){
  num_diags <- length(d)
  A_vec <- rows <- cols <- vector(mode = "list", length = num_diags)

  for(k in 1:num_diags){
    d_k <- d[k]
    if(d_k == 0){
      i <- 1:m
      j <- 1:n
    } else if(d_k > 0){
      i <- 1:(m-d_k)
      j <- (1 + d_k):n
    } else if(d_k < 0){
      i <- (1 - d_k):m
      j <- 1:(n + d_k)
    }
    rows[[k]] <- i
    cols[[k]] <- j
    A_vec[[k]] <- A[j,k]
  }

  rows  <- unlist(rows)
  cols  <- unlist(cols)
  A_vec <- unlist(A_vec)
  B <- sparseMatrix(i=rows, j=cols, x = A_vec,dims = c(m,n))
  return(B)
}

#' @name tridiag
#' @title Create a sparse tridiagonal matrix
#'
#' @description Create a sparse tridiagonal matrix of dgcMatrix class.
#'
#' @param n dimension of the square matrix
#' @param x subdiagonal (-1)
#' @param y diagonal (0)
#' @param z superdiagonal (+1)
#'
#' @return Sparse tridiagonal matrix
tridiag <- function(n, x=NULL, y=NULL, z=NULL){
  if(is.null(x)){
    x <- -1; y <- 2; z <- -1
  }
  if(is.null(z)){
    z <- y; y <- x; x <- n
  }

  if(max(c(length(x), length(y), length(z))) == 1){
    x <- x * rep(1, n-1)
    z <- z * rep(1, n-1)
    y <- y * rep(1, n)
  } else{
    nx <- length(x)
    ny <- length(y)
    nz <- length(z)
    if((ny - nx - 1) | (ny - nz - 1)){
      stop("tridiag:InvalidVectorArgDim")
    }
  }
  n <- length(y)
  spdiags(matrix(c(x, 0, y, 0, z), nrow = n), -1:1, n, n)
}
