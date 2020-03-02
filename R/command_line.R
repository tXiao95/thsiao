#' @name winPath
#' @title Change Windows filepaths to R format
#'
#' @description Convert Windows backslash filepaths to forward slashes to
#' be compatible with R
#'
#' @export
winPath <- function(path){
  return(chartr("\\", "/", path))
}

#' @name setwd_jump
#' @title setwd_jump
#'
#' @description Changes working directory to folder up the chain, without
#' needing to type full filepath or series of "../../.."
#'
#' @param dir -- name of directory up the chain of the current working directory
#' to change to
#'
#' @return None
#'
#' @export
setwd_jump <- function(dir){
  path <- getwd()
  path_arr <- strsplit(path, "/")[[1]]
  len <- length(path_arr)
  if (path_arr[len] == dir){
    message("Current working directory is already this folder")
  } else{
    for(i in 1:len){
      d <- path_arr[i]
      if(d == dir){
        string <- paste0(rep("..", len-i), collapse = "/")
        setwd(string)
        message("Changed working dir to ", getwd())
        return(invisible(NULL))
      }
    }
    stop("Directory does not exist in this chain")
  }
}
