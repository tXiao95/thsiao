context("Gallery matrices")

# Binomial matrix ---------------------------------------------------------

test_that("Binomial matrix", {
  for(n in 2:30){
    A <- binomial_matrix(n)
    expect_equal(norm(A %*% A - 2^(n-1) * diag(n)), 0)
    B <- A*2^((1-n)/2)
    expect_equal(norm(B %*% B - diag(n)), 0)
  }
})

# Tridiagonal matrix ------------------------------------------------------

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
  expect_equal(Matrix::norm(A), 0, tolerance = 0.00001)
})

# Fiedler matrix ----------------------------------------------------------

test_that("Fiedler matrix has one dominant positive eigenvalue and symmetric", {
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

# Circulant matrix --------------------------------------------------------

test_that("Circulant matrix is circular", {
  # Check every diagonal in the matrix is of the same value
  n <- 100
  A <- circul(n)
  B <- circul(rnorm(n))

  i <- row(A)
  j <- col(A)

  bound <- (2*n - 2) / 2
  for (k in -bound:bound){
    expect_equal(length(unique(A[i == (j + k)])), 1)
    expect_equal(length(unique(B[i == (j + k)])), 1)
  }
})
