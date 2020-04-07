context("Gallery matrices")

library(pracma)

# binomial_matrix: Binomial matrix ---------------------------------------------------------

test_that("Binomial matrix", {
  for(n in 10:13){
    A <- binomial_matrix(n)
    expect_equal(norm(A %*% A - 2^(n-1) * diag(n)), 0)

    # B is involutory
    B <- A*2^((1-n)/2)
    expect_equal(norm(B %*% B - diag(n)), 0)
  }
  # Test for n = 4
  A <- binomial_matrix(4)
  B <- matrix(c(1,1,1,1,3,1,-1,-3,3,-1,-1,3,1,-1,1,-1), nrow=4)
  expect_equal(norm(A-B),0)
})

# cauchy_matrix: Cauchy matrix --------------------------------------------

test_that("Cauchy matrix is cauchy", {
  set.seed(1)
  x <- rnorm(10)
  y <- rnorm(10)

  A <- cauchy_matrix(x, y)

  # Check identity of each matrix cell is correct for random x and y
  for(i in 1:nrow(A)){
    for(j in 1:ncol(A)){
      expect_equal(A[i,j], 1 / (x[i] + y[j]))
    }
  }

  # Two equal length vectors
  x <- 1:3
  y <- 2:4
  A <- cauchy_matrix(x, y)
  B <- matrix(c(1/3, 1/4,1/5,1/4,1/5,1/6,1/5,1/6,1/7),nrow=3)
  expect_equal(norm(A-B),0)

  # Scalar x
  x <- 3
  A <- cauchy_matrix(x)
  B <- matrix(c(1/2, 1/3,1/4,1/3,1/4,1/5,1/4,1/5,1/6),nrow=3)
  expect_equal(norm(A-B),0)

  # unequal lengths
  x <- 2:3
  y <- 1:100
  expect_error(cauchy_matrix(x, y))
})

# chebspec: Chebyshev spectral differentiation matrix ---------------------------------------------------------

test_that("Chebspec matrix is nilpotent", {
  C <- chebspec(4)
  # C is nilpotent, k=0
  expect_equal(norm(C %*% C %*% C %*% C), 0)

  # C has nullspace dim 1
  sigmas <- svd(C)$d
  smallest <- sigmas[length(sigmas)]
  next_smallest <- sigmas[length(sigmas)-1]
  expect_equal(smallest, 0)
  expect_true(next_smallest > sqrt(.Machine$double.eps))

  # and the null vector is all ones
  expect_equal(norm(C %*% c(1,1,1,1)), 0)

  # Eigenvalues have all negative real parts when k=1
  C <- chebspec(4, 1)
  eigs <- eigen(C)$values
  expect_true(all(Re(eigs) < 0))
})

# chebvand: Vandermonde-like matrix for Chebyshev polynomials ----------------------------------------------------------------

# test_that()

# chow: Singular Toeplitz lower Hessenberg matrix -------------------------
# circul: Circulant matrix --------------------------------------------------------

test_that("Circulant matrix is circular", {
  # Check every diagonal in the matrix is of the same value (Toeplitz)
  n <- 100
  A <- circul(n)
  B <- circul(rnorm(n))

  i <- row(A)
  j <- col(A)

  bound <- (2*n - 2) / 2
  for (k in -bound:bound){
    # Check that matrix is Toeplitz (diagonals all equal)
    expect_equal(length(unique(A[i == (j + k)])), 1)
    expect_equal(length(unique(B[i == (j + k)])), 1)

    # Check that matrix is circulant: ith diagonal equal to (i+n)th diagonal
    if(k < 0){
      expect_equal(unique(Diag(A, k = k)), unique(Diag(A, k = k + n)))
      expect_equal(unique(Diag(B, k = k)), unique(Diag(B, k = k + n)))
    }
  }
})

# clement: Tridagonal matrix with 0 diagonal ------------------------------

test_that("clement matrix", {
  # Diagonal is 0
  A <- clement(10)
  B <- clement(10,1)
  expect_equal(unique(diag(A)), 0)
  expect_equal(unique(diag(B)), 0)

  # Symmetric if k=1
  expect_true(isSymmetric(B))
})

# dorr: Dorr matrix -------------------------------------------------------

test_that("Dorr matrix is correct", {
  # theta 0.01
  Atrue <- matrix(c(1.32, -1.16, 0,
                    -.16, .32, -.16,
                    0, -1.16, 1.32), nrow=3, byrow=T)
  Atest <- as.matrix(dorr(3))
  expect_equal(norm(Atest - Atrue), 0)

  # negative theta
  Btrue <- matrix(c(-31, 15, 0,
                    16, -32, 16,
                    0, 15, -31), nrow=3, byrow=T)
  Btest <- as.matrix(dorr(3, -1))
  expect_equal(norm(Btest - Btrue), 0)
})

# dramadah: Anti-hadamard matrix ------------------------------------------

test_that("Dramadah matrix is correct", {
  A1 <- dramadah(4, 1)
  A2 <- dramadah(4, 2)
  A3 <- dramadah(4, 3)

  refA1 <- matrix(c(1,0,0,1,1,1,0,0,0,1,1,0,1,0,1,1), nrow=4)
  refA2 <- matrix(c(1,0,0,0,1,1,0,0,0,1,1,0,1,0,1,1), nrow=4)
  refA3 <- matrix(c(1,0,1,0,1,1,0,1,0,1,1,0,0,0,1,1), nrow=4)
  expect_equal(norm(A1 - refA1), 0)
  expect_equal(norm(A2 - refA2), 0)
  expect_equal(norm(A3 - refA3), 0)
})

# jordbloc ----------------------------------------------------------------

test_that("Jordan block is correct", {
  A1 <- jordbloc(4)
  A2 <- jordbloc(4, 2)

  refA1 <- matrix(c(1,0,0,0,1,1,0,0,0,1,1,0,0,0,1,1), nrow=4)
  refA2 <- refA1
  diag(refA2) <- 2

  expect_equal(norm(A1 - refA1), 0)
  expect_equal(norm(A2 - refA2), 0)
})

# tridiag: Sparse tridiagonal matrix ------------------------------------------------------

test_that("Tridiagonal matrix eigenvalues", {
  # Are tridiagonal matrix eigenvalues satisfied by identity
  n <- 50
  x <- 10
  y <- 2
  z <- 3

  eigs_analytic <- sort(y + 2*sqrt(x*z)*cos((1:n)*pi/(n+1)))
  A <- tridiag(n, x, y, z)
  eigs <- sort(eigen(A)$values)

  expect_equal(eigs, eigs_analytic, tolerance = 0.00001)
})

test_that("Tridiagonal matrix is tridiagonal", {
  n <- 1000
  x <- 10
  y <- 2
  z <- 3
  A <- tridiag(n, x, y, z)

  for(i in 1:(n-1)){
    A[i,i] <- A[i+1,i] <- A[i,i+1] <- 0
  }
  A[n,n] <- 0
  expect_equal(norm(A), 0, tolerance = 0.00001)
})

# fiedler: Fiedler Symmetric matrix ----------------------------------------------------------

test_that("Fiedler matrix has one dominant positive eigenvalue (rest negative) and symmetric", {
  A <- fiedler(10)
  eigs <- eigen(A)$values
  expect_more_than(eigs[1], 0)
  expect_equal(sum(eigs[2:length(eigs)] > 0), 0)
  expect_identical(isSymmetric(A), TRUE)

  v <- rnorm(10)
  A <- fiedler(v)
  eigs <- eigen(A)$values
  expect_more_than(eigs[1], 0)
  expect_equal(sum(eigs[2:length(eigs)] > 0), 0)
  expect_identical(isSymmetric(A), TRUE)
})

# lehmer: Lehmer matrix -----------------------------------------------------------

test_that("Lehmer matrix is correct", {
  A2 <- matrix(c(1, 1/2, 1/2, 1), nrow = 2)
  B2 <- lehmer(2)
  expect_equal(norm(A2 - B2), 0)

  A3 <- matrix(c(1, 1/2, 1/3, 1/2, 1, 2/3, 1/3, 2/3, 1), nrow = 3)
  B3 <- lehmer(3)
  expect_equal(norm(A3 - B3), 0)

  A4 <- matrix(c(1, 1/2, 1/3, 1/4, 1/2, 1, 2/3, 1/2, 1/3, 2/3, 1, 3/4, 1/4, 1/2, 3/4, 1), nrow = 4)
  B4 <- lehmer(4)
  expect_equal(norm(A4 - B4), 0)
})
