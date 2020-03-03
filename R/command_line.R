#' @name winPath
#' @title Change Windows filepaths to R format
#'
#' @description Convert Windows backslash filepaths to forward slashes to
#' be compatible with R
#'
#' @param path Filepath with \\ slashes that needs to be formatted. If left blank, will pull from the clipboard.
#'
#' @return Filepath with forward slashes
#'
#' @examples
#' # If Windows path on my clipboard
#' winPath()
#' # Input path
#' winPath("C:\\Users\\user")
#'
#' @export
winPath <- function(path="clipboard"){
  if(path=="clipboard"){
    path <- readClipboard()
  }
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
#' @examples
#' \dontrun{
#' setwd_jump("Users")
#'}
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
