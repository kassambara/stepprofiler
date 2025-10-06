#' DE_Results Class
#'
#' This class extends the \code{\link[S4Vectors]{DFrame}} class from the
#' S4Vectors package, designed specifically to hold results from
#' differential expression analysis, ensuring specialized methods can be
#' applied.
#'
#' @import methods
#' @importClassesFrom S4Vectors DFrame
#' @importFrom S4Vectors DataFrame mcols 'mcols<-' rownames 'rownames<-'
#' @exportClass DE_Results
setClass("DE_Results", contains = "DFrame")


#' DE_Results Constructor
#'
#' Creates an instance of the \code{DE_Results} class. This function
#' validates the input and ensures proper construction of the S4 object.
#'
#' @param x A data frame or any object coercible to a \code{DFrame}.
#' @return A \code{DE_Results} object.
#'
#' @export
#' @rdname DE_Results-class
DE_Results <- function(x) {
  if (!is(x, "DFrame")) x <- S4Vectors::DataFrame(x)

  # Use the standard S4 constructor 'new()' to create the instance
  obj <- methods::new("DE_Results", x)

  methods::validObject(obj)
  return(obj)
}


#' @rdname DE_Results-class
#' @importFrom methods callNextMethod
#' @export
setMethod("$<-", signature(x = "DE_Results", value = "ANY"),
function(x, name, value) {
  # Call the method for the parent class (DFrame) to perform the actual operation
  res <- callNextMethod()

  # Use the constructor function to re-cast the modified DFrame back to DE_Results
  res <- DE_Results(res)

  return(res)
})