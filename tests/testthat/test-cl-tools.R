context("Command line tools")

test_that("Changing directories correctly", {
  sample_path <- strsplit(getwd(), "/")[[1]]
  num_dir <- length(sample_path)
  for(i in num_dir:1){
    setwd_jump(sample_path[i])
    if(i > 1){
      testthat::expect_identical(getwd(), paste0(sample_path[1:i], collapse="/"))
    } else{
      testthat::expect_identical(getwd(), paste0(sample_path[i], "/"))
    }
  }
  # Make sure to change working directory back to original - leave User's space untouched
  setwd(paste0(sample_path, collapse="/"))
})
