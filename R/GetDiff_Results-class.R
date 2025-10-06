#' GetDiff_Results Class
#'
#' This class is a specialized subclass of \code{\link{DE_Results}},
#' intended for results that have undergone a specific "get_diff" transformation
#' or analysis step. It inherits all slots and methods from \code{DE_Results}.
#'
#' @import methods
#' @exportClass GetDiff_Results
setClass("GetDiff_Results", contains = "DE_Results")


#' GetDiff_Results Constructor
#'
#' Creates an instance of the \code{GetDiff_Results} class. It first coerces the
#' input to a \code{DE_Results} object and then re-casts it as the specialized
#' subclass.
#'
#' @param x A data frame or any object coercible to a \code{DE_Results} object.
#' @return A valid \code{GetDiff_Results} object.
#'
#' @export
#' @rdname GetDiff_Results-class
GetDiff_Results <- function(x) {
  # 1. Ensure the input is first converted to the immediate parent class
  #    (This handles input validation and structure from DE_Results)
  parent_obj <- DE_Results(x)

  # 2. Use the S4 constructor 'new()' to create the instance of the subclass
  obj <- methods::new("GetDiff_Results", parent_obj)

  # 3. Final validation check
  methods::validObject(obj)

  return(obj)
}