# Model fitting
#' @importFrom MASS mvrnorm
#' @importFrom stats .checkMFClasses delete.response model.frame model.matrix na.pass terms vcov

# simulateDraws generic ---------------------------------------------------

#' @name simulateDraws
#' @title Generate simulations of a fitted model object on new data
#'
#' @description \code{simulateDraws} is a generic function for generating simulation predictions
#' from a model object and newdata.
#'
#' @param object model object for which draws of prediction is desired
#' @param ... additional arguments affecting the prediction
#'
#' @export
simulateDraws <- function(object, ...){
  UseMethod("simulateDraws", object)
}

# simulateDraws.lm  -------------------------------------------------------

#' @name simulateDraws.lm
#' @title Generate simulations of a fitted \code{lm} object on new data
#'
#' @description Returns a matrix of simulated responses from a fitted lm object, applied to a new dataset
#'   separate from the original data used to fit the model.
#'
#' @param object object of class inheriting from "lm"
#' @param newdata a data frame where to look for variables to predict
#' @param nsim the number of simulations to return for each observation in the new dataset. Must be greater than 1.
#' @param ... further arguments passed to or from other methods
#'
#' @return a \code{nrow(newdata) x nsim} sized matrix. Each row corresponds to an observation in the new dataset,
#'   and each colulmn corresponds to a draw/simulation.
#'
#' @export
simulateDraws.lm <- function(object, newdata, nsim = 100, ...){
  tt <- terms(object)
  # Check whether model object inherits lm ... otherwise warning
  if (!inherits(object, "lm"))
    stop("Not a lm object")

  # Create design matrix for newdata
  Terms <- delete.response(tt)
  m <- model.frame(Terms, newdata, na.action = na.pass,
                   xlev = object$xlevels)
  if (!is.null(cl <- attr(Terms, "dataClasses")))
    .checkMFClasses(cl, m)
  X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)

  # Generate offset
  offset <- rep(0, nrow(X))
  if (!is.null(off.num <- attr(tt, "offset")))
    for (i in off.num) offset <- offset + eval(attr(tt,
                                                      "variables")[[i + 1]], newdata)
  if (!is.null(object$call$offset))
    offset <- offset + eval(object$call$offset, newdata)

  # Check if original design matrix needed pivoting, and apply to newdata
  p <- object$rank
  p1 <- seq_len(p)
  piv <- if (p)
    qr(object)$pivot[p1]
  if (p < ncol(X) && !(missing(newdata) || is.null(newdata)))
    warning("prediction from a rank-deficient fit may be misleading")

  # Generate draws of the coefficients from MVN, and add on offset
  beta <- object$coefficients
  V <- vcov(object)
  B <- mvrnorm(n = nsim, mu = beta, Sigma = V)
  predictor <- drop(X[, piv, drop = FALSE] %*% t(B[,piv]))
  if (!is.null(offset))
    # Broadcasting offset, which will sum offset to each column of predictor matrix
    predictor <- predictor + offset

  return(predictor)
}