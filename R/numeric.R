#' @importFrom assertthat is.scalar
#' @importFrom data.table data.table CJ
#' @importFrom directlabels geom_dl
#' @import ggplot2

# linScale: linear scaling of a vector ------------------------------------

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


# derivativeTest: Taylor's Theorem based derivative test ------------------

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

# plotPseudospectra: plot pseudospectrum for a matrix ---------------------

#' @name plotPseudospectra
#' @title Visualize pseudospectra of a matrix
#'
#' @description Plots contour plot on the complex plane for the
#' pseudospectra for the matrix A, under the given epsilon.
#'
#' @param A a square matrix
#' @param title title of the pseudospectra plot
#' @param eps epsilon values
#'
#' @return a ggplot object
#'
#' @export
plotPseudospectra <- function(A, title = "e-pseudospectrum", eps = c(1,.5,.1, .01, .001)){
  if(nrow(A) != ncol(A)){stop("A must be a square matrix")}

  n <- nrow(A)

  true_eigs <- eigen(A)$values
  eig_table <- data.table(a = Re(true_eigs), b = Im(true_eigs))

  min_Re <- min(eig_table$a)
  min_Im <- min(eig_table$b)
  max_Re <- max(eig_table$a)
  max_Im <- max(eig_table$b)

  # Try to narrow down where the contours are...ideally you would fix 'a' and the
  # minimum singular value at epsilon, and solve for the complex part 'b' but
  # I didn't know how to do that. So just finding all the computed minimum singular values
  # that are below some tolerance distance from the contour I want. Then let the contour
  #software handle the rest

  # Creating grid of scalars on the complex plane. Then calculate minimum singular value for
  # all matrices zI - A, where z \in C. data.table doesn't support complex numbers so need to make vector
  # outside of the data frame
  grid <- CJ(a = seq(min_Re - 5, max_Re + 5,.1),
             b = seq(min_Im - 5, max_Im + 5,.1))
  c    <-  complex(real = grid$a, imaginary =grid$b)
  # Method 1: calculate all singular values and take minimum
  grid[,sigma := sapply(c, function(c) min(svd(diag(c, n) - A)$d))]

  # Want to know if there's a way to use svds() or similar to calculate smallest singular values?
  # I don't want to calculate the inverse first...should experiment

  # Method 2: run svds to get 1st singular value and take reciprocal
  # grid[, sigma := sapply(c, function(c){
  #   B <- diag(c,n) - A
  #   X <- t(B) %*% B
  #   eigs_sym(X,k=1)$values
  #   #1 / eigs(Conj(t(B)) %*% B,k=1,sigma=0)$values
  # })]

  p <- ggplot(grid, aes(a, b)) +
    geom_point(data = eig_table, aes(a,b),col = "red",size=4,alpha = 0.5) +
    geom_contour(aes(z = sigma, colour = stat(level)),
                 size = 1.5,
                 breaks = c(eps)) +
    theme_dark() +
    scale_colour_distiller(palette = "Spectral", direction = 1) +
    xlab("Real") +
    ylab("Imaginary") +
    ggtitle(title) +
    labs(colour = "epsilon") +
    geom_dl(aes(label=..level.., z = sigma), method = "top.pieces",
            stat="contour",breaks = eps) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    coord_fixed()
  p
}


