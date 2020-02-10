#' @importClassesFrom S4Vectors DataFrame

#' @rdname DE_Results
#' @export
setClass("DE_Results", contains = "DataFrame")

#' DE_Results object and constructor
#'
#'DE_Results extends the DataFrame class of S4Vectors. It is used by \code{\link{lm_results}}
#' to wrap up the results table. This constructor function would not typically be used by "end users".
#'
#'
#' @param DataFrame a DataFrame of lm_results.
#'
#' @return a DE_Results object
#' @docType class
#' @aliases DE_Results-class
#' @rdname DE_Results
#' @export
DE_Results <- function(DataFrame) {
  methods::new("DE_Results", DataFrame)
}



