#' @importFrom pracma pascal Toeplitz Diag
#' @importFrom Matrix sparseMatrix norm triu tril

# binomial_matrix: binomial matrix ---------------------------------------------------------

#' @name binomial_matrix
#' @title Create binomial matrix
#'
#' @description Binomial matrix: an N-by-N multiple of an involuntory matrix with
#' integer entries such that $A^2 = 2^(N-1)*I_N$
#' Thus B = A*2^((1-N)/2) is involutory, that is B^2 = EYE(N)
#'
#' @param n - row dimension
#'
#' @export
binomial_matrix <- function(n){
  L <- abs(pascal(n, 1))
  U <- L[n:1,n:1]
  D <- diag((-2)^(0:(n-1)))

  return(L %*% D %*% U)
}

# cauchy_matrix: Cauchy matrix ------------------------------------------------------------------

#' @name cauchy_matrix
#' @title Create Cauchy matrix
#'
#' @description Arguments \code{x} and \code{y} are vectors of length \code{n}.
#'   \code{C[i,j] = 1 / (x[i] + y[j])}
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

  return(1 / (matrix(x, nrow = n, ncol = n) + matrix(y, nrow = n, ncol = n, byrow = T)))
}

# chebspec: Chebyshev spectral differentiation matrix ---------------------------------------------------------

#' @name chebspec
#' @title Create Chebyshev spectral differentiation matrix
#'
#' @description Chebyshev spectral differentiation matrix of order \code{n}. \code{k} determines
#'   the character of the output matrix. For either form, the eigenvector matrix is ill-conditioned.
#'
#' @param n order of the matrix.
#' @param k \code{k=0} is the default, no boundary conditions. The matrix is similar to a Jordan block
#'   of size \code{n} with eigenvalue 0. If \code{k=1}, the matrix is nonsingular
#'   and well-conditioned, and its eigenvalues have negative real parts.
#'
#' @export
chebspec <- function(n, k = NULL){
  if(is.null(k)){
    k <- 0
  }
  if(k == 1){
    n <- n + 1
  }
  if(!(k %in% 0:1)){
    stop("k must be 0 or 1")
  }

  n <- n - 1
  C <- matrix(0, nrow = n + 1, ncol = n + 1)

  x <- cos(matrix(0:n, ncol = 1) * (pi/n))
  d <- matrix(1, nrow = n + 1, ncol = 1)
  d[1] <- 2; d[n + 1] <- 2

  X <- matrix(x, nrow = n + 1, ncol = n + 1) - matrix(x, nrow = n + 1, ncol = n + 1, byrow = T)
  C <- (d %*% t(matrix(1, nrow = n + 1, ncol = 1) / d)) / (X + diag(nrow(C)))

  # Now fix diagonal and signs
  C[1,1] <- (2*n^2 + 1) / 6
  for(i in 2:(n + 1)){
    if((i %% 2) == 0){
      C[,i] <- -C[,i]
      C[i,] <- -C[i,]
    }
    if(i < (n+1)){
      C[i,i] <- -x[i] / (2 * (1 - x[i]^2))
    } else{
      C[n + 1,n + 1] <- -C[1,1]
    }
  }

  if(k == 1){
    C <- C[2:(n + 1),2:(n + 1)]
  }

  return(C)
}

# chebvand: Vandermonde-like matrix for Chebyshev polynomials ---------------------------------------------------------

#' @name chebvand
#' @title Creating Vandermonde-like matrix for the Chebyshev polynomials
#'
#' @description Produces the (primal) Chebyshev Vandermonde matrix based on
#'   the points \code{p}. \code{C[i,j] = T_{i-1}p[j]}, where \code{T_{i-1}} is the Chebyshev
#'   polynomial of degree \code{i-1}
#'
#' @param p points to evaluate. If a scalar, then \code{p} equally spaced points
#'   on \code{[0,1]} are used.
#' @param m number of rows of the matrix. \code{chebvand(p, m)} is the rectangular version of
#'   \code{chebvand(p)} with \code{m} rows.
#'
#' @export
chebvand <- function(p, m = NULL){
  # If no rectangular row argument, make square matrix. Else, rectangular and m != n
  if(is.null(m)){
    #m <- p
    square <- 1
  } else{
    square <- 0
  }
  n <- length(p)

  # If p is scalar, make ncols = p and change p to vector within [0,1] of length p
  if(n == 1){
    n <- p
    p <- seq(0, 1, length.out = n)
  }

  # if will be square matrix, make nrows = ncols
  if(square == 1){
    m <- n
  }

  # make p a row vector
  p <- matrix(p, ncol = length(p))
  # Create m x n matrix, and since first row will be all 1's.
  C <- matrix(1, nrow = m, ncol = n)
  if(m == 1){
    return(C)
  }
  C[2,] <- p
  if(m == 2){
    return(C)
  }
  for(i in 3:m){
    C[i,] <- 2*p*C[i-1,] - C[i-2,]
  }

  return(C)
}

# chow: Singular Toeplitz lower Hessenberg matrix -------------------------------

#' @name chow
#' @title Creating singular Toeplitz lower Hessenberg matrix
#'
#' @description returns matrix \code{A = H(alpha) + delta * EYE}, such that
#'   \code{H[i,j] = alpha^(i-j+1)}.
#'
#' @param n order of the matrix
#' @param alpha defaults to 1
#' @param delta defaults to 0
#'
#' @export
chow <- function(n, alpha = 1, delta = 0){
  row <- alpha^(1:n)
  col <- c(alpha, 1, rep(0,n-2))
  A <- Toeplitz(row, col) + delta * diag(n)

  return(A)
}

# circul: Circulant matrix --------------------------------------------------------

#' @name circul
#' @title Create circulant matrix
#'
#' @description Each row is obtained from the previous by cyclically permuting the
#'   entries one step forward. A special Toeplitz matrix in which diagonals "wrap around"
#'
#' @param v first row of the matrix. If \code{v} is a scalar, then \code{C = circul(1:v)}
#'
#' @return a circulant matrix whose first row is the vector \code{v}
#'
#' @export
circul <- function(v){
  if(length(v) == 1){
    v <- 1:v
  }
  n <- length(v)
  A <- pracma::Toeplitz(c(v[1], v[n:2]), v)

  return(A)
}


# clement: Tridiagonal matrix with zero diagonal entries ------------------

clement <- function(n, k = 0){
  n <- n - 1
  z <- 1:n
  x <- n:1

  if(k == 0){
    A <- Diag(x, -1) + Diag(z, 1)
  } else{
    y <- sqrt(x * z)
    A <- Diag(y, -1) + Diag(y, 1)
  }

  return(A)
}

# compar: Comparison matrices ---------------------------------------------

compar <- function(A, k = 0){
  m <- nrow(A)
  n <- ncol(A)

  if(k == 0){
    C <- -abs(A)
  }
}

# Sparse Diagonal Matrix --------------------------------------------------

#' @name spdiags
#' @title Create sparse diagonal matrix
#'
#' @description Creates a sparse representation of multiple diagonal matrix
#'
#' @param A matrix where columns correspond to the desired diagonals
#' @param d indices of the diagonals to be filled in. 0 is main diagonal. -1
#'   is first subdiagonal and +1 is first superdiagonal.
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
  B <- sparseMatrix(i = rows, j = cols, x = A_vec,dims = c(m, n))

  return(B)
}

# Sparse tridiagonal matrix -----------------------------------------------

#' @name tridiag
#' @title Create sparse tridiagonal matrix
#'
#' @description Create a sparse tridiagonal matrix of dgcMatrix class.
#'
#' @param n dimension of the square matrix
#' @param x subdiagonal (-1)
#' @param y diagonal (0)
#' @param z superdiagonal (+1)
#'
#' @return Sparse tridiagonal matrix
#'
#' @export
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

  return(spdiags(matrix(c(x, 0, y, 0, z), nrow = n), -1:1, n, n))
}

# Fiedler Symmetric matrix ------------------------------------------------

#' @name fiedler
#' @title Create Fiedler matrix
#'
#' @description Fiedler matrix that has a dominant positive eigenvalue and all others are negative
#'
#' @param c N-vector. If \code{c} is a scalar, then returns fiedler(1:c)
#'
#' @return a symmetric dense matrix A with a dominant positive eigenvalue and all others are negative.
#'
#' @export
fiedler <- function(c){
  if (!is.vector(c)){
    stop("'c' is not a vector")
  }
  if(length(c) == 1){
    c <- 1:c
  }
  n <- length(c)
  A <- matrix(raw(), n, n)
  A <- matrix(data = abs(c[col(A)] - c[row(A)]), nrow = n, ncol = n)

  return(A)
}

# Toeplitz matrix with sensitive eigenvalues ------------------------------

#' @name grcar
#' @title Create Toeplitz matrix with sensitive eigenvalues
#'
#' @description Eigenvalues are sensitive.
#'
#' @param n dimension of the square matrix
#' @param k number of superdiagonals of ones
#'
#' @return n by n Toeplitz matrix with -1 on subdiagonal, 1 on diagonal, and k superdiagionals of 1s.
#'
#' @export
grcar <- function(n, k=NULL){
  if(is.null(k)){
    k <- 3
  }
  A <- tril(triu(matrix(1, nrow = n, ncol = n)), k)
  i <- row(A)
  j <- col(A)
  A[i == j + 1] <- -1

  return(A)
}

# Leslie matrix -----------------------------------------------------------

#' @name leslie
#' @title Create Leslie population model matrix
#'
#' @description N by N matrix from Leslie population model with average birth and survival rates.
#'
#' @param a average birth numbers (first row)
#' @param b survival rates (subdiagonal)
#' @param sparse whether to return a sparse matrix
#'
#' @return N by N Leslie population model matrix
#'
#' @export
leslie <- function(a, b=NULL, sparse = F){
  if(is.null(b)){
    n <- a
    a <- rep(1, n)
    b <- rep(1, n-1)
  }
  if(length(a) != length(b) + 1){
    stop("a must have length = length(b) + 1")
  }
  if(sparse){
    n <- length(a)
    i <- c(rep(1, n), 2:n)
    j <- c(1:n, 1:(n-1))
    L <- sparseMatrix(i = i, j = j, x = c(a, b))
  } else{
    L <- Diag(b, -1)
    L[1,] <- a
  }

  return(L)
}

# Lauchli matrix ----------------------------------------------------------

#' @name lauchli
#' @title Create Lauchli Matrix
#'
#' @description the (N + 1) x (N) matrix [ones(1,n); mu*eye(n)]. Well-known example in least squares of the
#'   danger of forming t(A) %*% A (due to inexact arithmetic, gives singular matrix)
#'
#' @param n number of columns
#' @param mu constant applied to identity
#' @param sparse whether matrix should be sparse
#'
#' @return Lauchli matrix.
#'
#' @export
lauchli <- function(n, mu = NULL, sparse = F){
  if(is.null(mu)){
    mu <- sqrt(.Machine$double.eps)
  }
  if(sparse){
    rows <- c(rep(1, n), 2:(n + 1))
    cols <- c(1:n, 1:n)
    A <- sparseMatrix(i = rows, j = cols, x = c(rep(1, n), rep(mu, length(rows) - n)))
  } else{
    A <- mu * diag(n)
    A <- rbind(rep(1, n), A)
  }

  return(A)
}

# Lehmer matrix -----------------------------------------------------------

#' @name lehmer
#' @title Create Lehmer matrix
#'
#' @description the symmetric positive-definite matrix such that A[i,j] = i/j, for j >= i
#'
#' @param n order of matrix
#'
#' @export
lehmer <- function(n){
  A <- matrix(0, n, n)
  i <- row(A)
  j <- col(A)
  I <- i / j
  J <- j / i
  A[j >= i] <- I[j >= i]
  A[j < i] <- J[j < i]

  return(A)
}
