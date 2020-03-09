context("Numeric functions")

test_that("Changing directories correctly", {
  x <- 1:10

  # Scaling 1:10 to 0:1
  xScale <- linScale(x, 0, 1)
  correct <- (x - 1) / 9
  expect_equal(sum((correct - xScale)^2), 0)

  # Scaling 1:10 to 1:10
  expect_equal(sum((linScale(x, 1, 10) - x)^2), 0)

  # Changing order of 1:10 reverse
  expect_equal(sum((linScale(10:1, 1, 10) - 10:1)^2), 0)
})

test_that("derivative Test test", {
  # See if test can detect correct derivative of x^2
  f <- function(x){x^2}
  f_grad <- function(x){2*x}
  f_grad_wrong <- function(x){2.01*x}

  h <- c(2e-1, 2e-2, 2e-3, 2e-4, 2e-5)

  W <- 3
  D <- 1

  derivs <- derivativeTest(f, f_grad, h, W, D)
  derivsWrong <- derivativeTest(f, f_grad_wrong, h, W, D)
  log_h <- log(h)
  log_d2 <- log(derivs$d2)
  log_d2_wrong <- log(derivsWrong$d2)

  # Calculate slopes of each d1 and d2
  slope_d2 <- diff(log_d2) / diff(log_h)
  slope_d2_wrong <- diff(log_d2_wrong) / diff(log_h)

  # Correct derivative - slope == 2
  expect_equal(sum((slope_d2 - 2)^2), 0)
  # Wrong derivative - slope != 2
  expect_true(sum((slope_d2_wrong - 2)^2) != 0)
})
