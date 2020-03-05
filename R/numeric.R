#' @importFrom assertthat is.scalar
#'
#' @name linScale
#' @title Linear scaling of numeric vectors
#'
#' @description Linear transformation of a numeric vector to a
#' newly set minimum and maximum.
#'
#' @param x numeric vector
#' @param newmin new minimum
#' @param newmax new maximum
#'
#' @return A numeric vector linearly scaled of the input to be between \code{newmin} and \code{newmax}.
#'
#' @export
#'
#' @examples
#' # Scale samples from standard normal between 0 and 1
#' x <- rnorm(100)
#' scaled_x <- linScale(x, 0, 1)
#' plot(x, scaled_x)
linScale <- function(x, newmin = 0, newmax = 1){
  oldmin <- min(x)
  oldmax <- max(x)
  return((x - oldmin)/(oldmax - oldmin) * (newmax - newmin) + newmin)
}

#' @name derivativeTest
#' @title Testing derivatives based on the Taylor Theorem
#'
#' @description Derivative check for gradients, Jacobians, and Hessians based on Taylor's Theorem.
#'
#' @param f scalar-valued function to be evaluated
#' @param f_grad function for derivative of \code{f}. Should take exact same arguments as \code{f}, and be a row vector.
#' @param h step sizes: must be (0,1]
#' @param W fixed constant to evaluate the function at. If multivariate, should be a column vector.
#' @param D direction to move in from W. If multivariate, should be a column vector
#' @param ... Any other arguments to evaluate \code{f(x)} and \code{f_grad(x)}
#'
#' @details In many optimizations, automatic differentiation and other numerical methods are preferred when
#' pursuing derivative-based optimization. However, derivatives are assumed to be wrong unless proven otherwise. This
#' derivative test makes use of properties of the Taylor expansion for any differentiable function to identify whether
#' a computed gradient is correct, given that the function is the correct one.
#'
#' @return list of two vectors, \code{d1} holds the first order differences which have slope close to 1, \code{d2} should
#' be the second order differences which have slope close to 2.
#'
#' @export

derivativeTest <- function(f, f_grad, h, W, D, ...){
  num_steps <- length(h)
  d1 <- numeric(num_steps)
  d2 <- numeric(num_steps)

  for(i in 1:num_steps){
    E_hD <- f(W + h[i]*D, ...)
    E_W <- f(W, ...)
    d1[i] <- E_hD - E_W
    grad_E <- f_grad(W, ...)
    if(is.scalar(W)){
      d2[i] <- E_hD - E_W - (h[i]*D)*grad_E
    } else{
      d2[i] <- E_hD - E_W - h[i]*sum(diag(D %*% grad_E))
    }
  }
  return(list(d1=abs(d1), d2=abs(d2)))
}




