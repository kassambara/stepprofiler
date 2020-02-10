#' Summarize gene expression data
#' @description
#' Summarize gene expression data by computing some descriptive statistics.
#' \itemize{
#' \item summarize(): computes descriptive statistics for expression data.
#' \item summarizeby(): summarizes expression data by groups and by \bold{one descriptive statistic}.
#' For example computes, the
#' average expression of each gene by sample groups.
#' \item summarizebys(): summarizes expression data by groups
#' and by \bold{multiple descriptive statistics} at the same time. For example computes, the
#' mean expression and the SD of each gene by sample groups.
#' }
#' @return
#' \itemize{
#' \item summarize() returns a data frame containing some descriptive statistics for each genes including:
#' the minimum (Min), the first quartile (Q1), the median (Med), the mean (Mean), the third quartile (Q3),
#' the maximum (Max), the standard deviation (SD), the interquartile range (IQR), the median absolute deviation (MAD).
#' \item summarizeby() retuns a data frame containing the summarized expression data by groups.
#' }
#' @param x matrix or numeric data frame containing gene expression data. Rows are genes and columns are samples.
#' @name summarize
#' @rdname summarize
#' @export
summarize <- function(x){
  distributions <- data.frame(
    Min = apply(x, 2, min), # minimum
    Q1 = apply(x,2, quantile, 1/4), # First quartile
    Med = apply(x,2, median), # median
    Mean = apply(x, 2, mean), # mean
    Q3 = apply(x,2, quantile, 3/4), # Third quartile
    Max = apply(x, 2, max), # Maximum
    SD = apply(x, 2, sd),
    IQR = apply(x, 2, IQR),
    MAD = apply(x, 2, mad)
  )
  distributions <- round(distributions, 2)
  distributions
}

#' @param fun a statistical function such as mean, median etc
#' @param groupby a factor variable specifying sample groups
#' @param melt logical value. If TRUE, then the data is summarized by groups
#'  and melted using reshape2 package
#' @rdname summarize
#' @export
summarizeby <- function(x, groupby, fun, melt = FALSE){
  val <- apply(x, 1,
               function(x, groupby){ tapply(unlist(x), groupby, fun)},
               groupby)
 val <- t(val)

 if(melt){
   val <- as.data.frame(t(val))
   val$ids <- rownames(val)
   val <- reshape2::melt(val, id.vars = "ids")
 }
 val
}

#' @rdname summarize
#' @export
#' @param stat a vector containing the name of the statistic summaries you want to compute.
#' Possible values are the combination of the following statistics:
#' c("mean", "sd", "min", "max", "median", "mad", "iqr")
summarizebys <- function(x, groupby,
                         stat = c("mean", "sd"))
{
  if(!inherits(x, "data.frame"))
    x <- as.data.frame(x)
  if(!inherits(groupby, "factor")) groupby <- as.factor(groupby)

  res <- summarizeby(x, groupby, mean, melt = TRUE)
  if(!"mean" %in% stat) res <- res[, -3]
  else colnames(res)[3] <- "mean"
  if("sd" %in% stat) res$sd <- summarizeby(x, groupby, stats::sd, melt = TRUE)$value
  if("min" %in% stat) res$min <- summarizeby(x, groupby, base::min, melt = TRUE)$value
  if("median" %in% stat) res$median <- summarizeby(x, groupby, median, melt = TRUE)$value
  if("max" %in% stat) res$max <- summarizeby(x, groupby, base::max, melt = TRUE)$value
  if("iqr" %in% stat) res$iqr <- summarizeby(x, groupby, stats::IQR, melt = TRUE)$value
  if("mad" %in% stat) res$mad <- summarizeby(x, groupby, stats::mad, melt = TRUE)$value
  res$ids <- factor(res$ids, levels = levels(groupby))
  res
}
