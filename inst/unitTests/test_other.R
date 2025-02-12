test_other <- function () {
  
  # covers helper and other internal functions
  
  data(ramr)
  
  .data <- ramr:::.preprocessData(
    data.ranges=ramr.data, data.samples=ramr.samples,
    data.coverage=data.frame(), transform="identity",
    exclude.range=c(2,0), ncores=1, verbose=TRUE
  )

  ramr:::rcpp_extract_out(.data)
  t(ramr:::rcpp_extract_coef(.data))
}