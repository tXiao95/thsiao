context("simNew")

test_that("simNew works for linear models", {
  model <- lm(Sepal.Length ~ Sepal.Width, iris)
  new <- rbind(iris, iris)

  draws <- simNew(model, newdata = new, nsim = 10)
  pred <- draws$Y
  Beta <- draws$B

  X <- cbind(rep(1, nrow(new)), new[,2])
  Y <- X %*% t(Beta)

  expect_equal(norm(pred - Y), 0)
})
