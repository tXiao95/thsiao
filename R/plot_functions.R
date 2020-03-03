# Plotting
#' @name ggplot_to_pdf
#' @title Printing ggplot objects to PDF's
#'
#' @description Printing plots to a pdf
#'
#' @export
ggplot_to_pdf <- function(plot, ...){
  pdf(...)
  print(plot)
  dev.off()
}
