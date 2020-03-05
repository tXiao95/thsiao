# Model fitting
#' @importFrom MASS mvrnorm

#' @title Generate Monte Carlo draws of a linear model applied to new observations. Works for
#' lm and glm objects. Does not work for lmer random and mixed effects models.
create_draws <- function(model, newdata, nsim){
  betas <- coef(model)
  V <- vcov(model)

  nsim = 100
  coef_draws <- MASS::mvrnorm(n = nsim, mu = betas, Sigma = V)

  coef_draws %*% newdata
}
