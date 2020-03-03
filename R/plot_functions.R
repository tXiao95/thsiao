#' @importFrom grDevices dev.off pdf
#'
#' @name ggplot_to_pdf
#' @title Printing ggplot objects to PDF's
#'
#' @description Printing plots to a pdf
#'
#' @param plot a plot object. ggplot most common
#' @param ... Arguments to pass to \code{pdf}
#'
#' @export
ggplot_to_pdf <- function(plot, ...){
  pdf(...)
  print(plot)
  dev.off()
}
