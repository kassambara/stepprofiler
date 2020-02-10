#' @include summarize.R
NULL
#' Expressed genes in each cell population
#'
#' @description Identification of expressed genes in each cell population.
#' Expressed genes in a given cell population is defined as genes with average expression above a given cutoff.
#' @param x matrix or numeric data frame containing gene expression data. Rows are genes and columns are samples.
#' @param groups a factor containing the group assignement of samples
#' @param cutoff a numeric value specifying the threshold above which a gene is declared as expressed.
#' @return
#' \itemize{
#' \item expressed() returns an object of class "expressed". Columns are samples or groups and rows are genes. Values can be 1 (if expressed)
#' or 0 (if not expressed).
#' \item summarize() returns a list containing the following components:
#'  \itemize{
#'    \item exprs_in_all: The number of expressed genes in all groups
#'    \item exprs_in_one: The number of expressed genes in at least one group
#'    \item exprs_in_each: a vector containing the number of genes expressed in each group
#'    \item exprs_in_specific: the number of genes expressed only in specific group
#'  }
#' }
#' @name expressed
#' @rdname expressed
#' @export
expressed <- function(x, groups, cutoff = 1){

if(!inherits(x, c("matrix", "data.frame")))
  stop("x must be a matrix or data frame containing gene expression data")
  if(ncol(x)!=length(groups))
    stop( "The number of columns in x is different to the length of groups")
 if(!inherits(groups, "factor")) groups <- as.factor(groups)
 if(length(levels(groups)>=2)) x <- summarizeby(x, groups, mean)
 x <- apply(x, 2, function(x){ifelse(x >= cutoff, 1, 0)})
 class(x) <- c("expressed", class(x))

 x
}

#' @param object an object of class expressed
#' @param ... additional arguments.
#' @rdname expressed
#' @export
summary.expressed <- function(object, ...){
  if(!inherits(object, "expressed"))
    stop("An object of class expressed is required.")

  res = list()
  row_min <- apply(object, 1, min)
  row_max <- apply(object, 1, max)
  res$exprs_in_all <- length(which(row_min == 1))
  res$exprs_in_one <- length(which(row_max == 1))
  res$exprs_in_each <- apply(object, 2, sum)
  exprs_in_specific = NULL
  for(i in 1:ncol(object)){
    index <- which(object[, i] == 1)
    if(length(index) == 0) exprs_in_specific <- c(exprs_in_specific, 0)
    else{
      res2 <- apply(object[index,  -i, drop = FALSE], 1, sum)
      exprs_in_specific <- c(exprs_in_specific, length(which(res2 == 0)))
    }
  }
  names(exprs_in_specific) <- colnames(object)
  res$exprs_in_specific <- exprs_in_specific

  res
}
