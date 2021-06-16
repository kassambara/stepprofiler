#' @include utilities.R expressed.R utilities_diffexprs.R
#' @importFrom dplyr %>%
NULL
#' Differential Expression Between Groups
#'
#' @description Identify differentially expressed genes between groups:
#' \itemize{
#' \item \bold{get_diff_1()}: Common up/down regulated genes in one group compared to each
#' of the other groups in pairwise analysis. It's an alias of get_diff_pwc() when length(grp1) = 1.
#' \item \bold{get_diff_2()}: Differentially expressed genes between two groups in a simple two class analysis.
#' Alias of get_diff_pwc() when length(grp1) = length(grp2) = 1.
#' \item \bold{get_diff_n()}: Differentially expressed genes between multiple groups
#' using pairewise analysis or not (see details section).
#' \item \bold{get_diff_pwc()}: Pairwise comparisons (pwc) between sample classes in grp1 (e.g.: grp1 = c("A", "B", ...))
#' and sample classes in grp2 (e.g.: grp2 = c("C", "D", ...) ). Can do all the above mentioned type of analyses.
#' It's an alias of get_diff_n(), get_diff_1( ) and get_diff_2().
#' \item \bold{significant()}: Take an object of class "get_diff", "lm_results" or "DESeqResult" and returns
#' significantly differentially expressed genes with a given cutoff for false discovery rate
#' and fold change.
#' }
#' For more information about these functions, read details section.
#' @param dds an object of class:
#' \itemize{
#'  \item DESeqDataSet containing the results of the function DESeq(), or
#'  \item lm_fit
#' }
#' dds must contain the information about samples.
#' @param factor the name of a factor in the design formula.
#' @param grp1,grp2 a list of 2 character vectors containing the groups to be compared.
#' The argument grp1 is used as numerator and grp2 as the denominator for computing the fold change (See details  section).
#' If grp2 is NULL, then grp1 is compared to all the other groups in the data.
#' @param alpha The significance cutoff used for optimizing the independent filtering
#' (by default 0.05); see ?DESeq2::results.
#' @param active_exprs a numeric value specifying detected expression cutoff.
#' If not NULL, any gene with an average expression above this cutoff, in at least one group,
#' is declared as actively expressed. Used for defining presence and absence detection call.
#' If not NULL, a column named detection_call is added to the result. Possible values for detection_call
#' is 0 (the gene is not actively expressed) or 1 (the gene is actively expressed).
#' @param count.norm a data frame containing the normalized count as provided by DESeq2.
#' Only used for computing the detection call when active_expression is not NULL.
#' @param exprs.norm a data frame containing the normalized gene expression data. An alias of count.norm.
#' @param merge a logical value. Default is FALSE. Used only in multiclass comparison cases (seee details section).
#' @param verbose logical value; if TRUE the progression is shown. Default is TRUE.
#'
#'
#' @details
#' Definition of general terms contained in the results returned by get_diff_*() functions:
#' \itemize{
#' \item The \bold{log2FoldChange} is calculated using grp1 as numerator and grp2 as denominator.
#' \item The \bold{baseMean} is calculated as mean of normalized counts for all samples.
#' \item The \bold{detection_call} defines gene presence/absence detection calls which possible values are:
#' 0 (the gene is not actively expressed) or 1 (the gene is actively expressed,
#' i.e. normalize expression >= active_exprs cutoff).
#' \item \bold{significance}: Deregulation significance: Up = 1; Down = -1; Not significant = 0. A gene is
#' significantly deregulated if a) the detection_call = 1, b) the fold change >= fc cutoff
#' and c) the padj < = alpha cutoff.\cr\cr
#' }
#' Possible values for the arguments grp1 and grp2 can be of the forms:
#' \itemize{
#' \item \bold{grp1 = "A" and grp2 = "B"}: Simple two class comparisons.
#' Recommended functions: \bold{get_diff_2()} or \bold{get_diff_pwc()}.
#' The \bold{adjusted p-values} (padj) and the \bold{log2 fold changes} (log2FoldChange) are those returned by
#' DESeq2::results.
#' \item \bold{grp1 = "A" and grp2 = c("B", "C", ...)}. One class compared to more than one other groups.
#' Recommended functions: \bold{get_diff_1()} or \bold{get_diff_pwc()} for identifying common up/down regulated genes
#' in group A compared to each of groups B and C using pairwise analysis (e.g.: A vs B and A vs C).
#' \itemize{
#' \item \bold{Adjusted p-value}: The maximum of all the adjusted p-values (padj) obtained from the different pairwise comparisons
#' is used as the final padj in the results. For a given gene, a max(padj) < 0.05 indicates that each of the padj is < 0.05 in
#' all pairwise comparisons. If a gene is upregulated in one comparison but downregulated in another comparison,
#' then the max(padj) is set to NA. Only genes commonly up- or down- regulated in all the comparisons are considered to be
#' specific gene signatures of grp1.
#' \item \bold{log2 fold change}: The minimum of all log2FoldChange (lfc) obtained from the different pairwise comparisons
#' is used as the final log2FoldChange. For a given gene, min(lfc) > 1.5 in absolute value (abs) indicates that
#' each of the abs(lfcs) is > 1.5 in all pairwise comparisons. If a gene is upregulated in one comparison but downregulated in another comparison,
#' then the min(lfc) is set to 0. Only genes commonly up- or down- regulated in all the comparisons are considered to be
#' specific gene signatures of grp1.
#' \item \bold{lfcSE, stat and pvalue}: These columns can be ignored in the result as they are
#' obtained by performing a global comparison (not a pairwise comparison) between grp1 and grp2.
#' We kept them only for conveniance with DESeq2 result formats.
#' }
#' \item \bold{grp1 = c("A", "B", ...) and grp2 = c("D", "E", ...)}; Multiple class comparisons.
#' Recommended functions: \bold{get_diff_n()} or \bold{get_diff_pwc()}.
#'  \itemize{
#'    \item If merge = TRUE, then A+B are merged and compared to merged D+E in a single two class analysis.
#'    The \bold{adjusted p-values} (padj) and the \bold{log2 fold changes} (log2FoldChange)
#'    are those returned by DESeq2::results.
#'    \item If merge = FALSE, then separated pairewise comparisons are performed using \bold{get_diff_pwc()},
#'    for example as follow: "A" vs c("D", "E", ...) and "B" vs c("D", "E", ...). In this particular case, the
#'    different comparisons are  A vs D, A vs E, B vs D and B vs E, etc.
#'    Kept genes in the final results, are only those commonly up/down regulated
#'    in both A and B compared to D and E.
#'    The returned results have the same format as those returned by get_diff_1().
#'    Read the paragraphs above concerning \bold{adjusted p-value} and {log2 fold change}.
#' }
#' }
#'
#' @return
#' \itemize{
#' \item The functions get_diff_*() returns an object of class DESeqResults (see ?DESeq2::results) and "get_diff" which contains
#' the following columns: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, detection_call and significance.
#' To see the description of the columns use the function mcols(res, use.names = TRUE), where res is the
#' result of get_diff_*().\cr\cr
#' The function plot() can be used to draw the results. (See ?plot_MA).
#' \item significant() returns an object of class "significant" containing
#' significantly differentially expressed genes
#' with a given cutoff for false discovery rate (e.g.: 0.05) and fold change (e.g.: 1.5) and,
#' a detection_call =1 (if available in the argument object).\cr\cr
#' The object "significant" is a list with the following components:
#' \itemize{
#'  \item up: data frame containing upregulated genes in grp1 compared to grp2
#'  \item down: data frame containing downregulated genes in grp1 compared to grp2
#'  \item nup: number of upregulated genes
#'  \item ndown: number of downregulated genes
#' }
#' }
#' @name get_diff
#' @rdname get_diff
#' @export
get_diff_1 <- function(dds, factor, grp1, grp2 = NULL,
                     alpha = 0.05, fc = 1.5, active_exprs = 1, count.norm = NULL, exprs.norm = NULL,
                     verbose = TRUE, ...)
{

  if (!missing(count.norm)) {
    warning("argument count.norm is deprecated; please use exprs.norm instead.",
            call. = FALSE)
    exprs.norm <- count.norm
  }

  if(length(grp1)!=1)
    stop("grp1 must be of length 1")
  grps <- .check_groups(dds, factor, grp1, grp2)
  grp1 <- grps$grp1
  grp2 <- grps$grp2

  if(is.null(exprs.norm)) exprs.norm <- .data_norm(dds)
  if(length(grp2) == 1){
    res <- get_diff_2(dds, factor = factor, grp1 = grp1, grp2 = grp2,
               alpha = alpha, fc = fc, active_exprs = active_exprs, exprs.norm = exprs.norm,
               verbose = verbose, ...)
    res$significance <- .get_significance(res, alpha = alpha, fc = fc)
  }
  else res <- get_diff_1 (dds, factor = factor, grp1 = grp1, grp2 = grp2,
                            alpha = alpha, fc = fc, active_exprs = active_exprs,
                            exprs.norm = exprs.norm, verbose = verbose, ...)

  res <- .set_mcols(res) # description of column names
  res
}



#' @param ... other arguments to be passed to the functions DESeq2:results() or lm_results()
#' @rdname get_diff
#' @export
get_diff_2 <- function(dds, factor, grp1, grp2,
                       alpha = 0.05, fc = 1.5, active_exprs = 1, count.norm = NULL, exprs.norm = NULL,
                       verbose = TRUE, ...)
{

  # Deprecated arguments
  if (!missing(count.norm)) {
    warning(
      "argument count.norm is deprecated; ",
      "please use exprs.norm instead.", call. = FALSE
      )
    exprs.norm <- count.norm
  }
  . <- NULL

  # Check parameters
  if(length(grp1)!=1) stop("grp1 must be of length 1")
  if(length(grp2)!=1) stop("grp2 must be of length 1")
  if(fc < 1) stop("fc must be >= 1.")
  grps <- dds %>% .check_groups(factor, grp1, grp2)
  grp1 <- grps$grp1
  grp2 <- grps$grp2

  # Create a name for the current comparison
  step_name <- paste0(grp2, "_vs_", grp1) %>% tolower()
  if(verbose) message(
    "- Comparison ", step_name, ": Progress = ",
    1, "/", length(grp1)
    )

  # Differential expression between grp1 and grp2
  res <- dds %>% .results(c(factor, grp1, grp2), alpha = 0.05, ...)
  gene.id <- rownames(res)
  # Mark genes that are actively expressed
  res$detection_call <- dds %>% .get_detection_call(
    factor, grp1, grp2,active_exprs, exprs.norm ) %>%
    .[gene.id]

  # Mark significant genes
  res <- res %>% mark_significant(log2.foldchange = log2(fc), pvalue = alpha)
  res <- .set_mcols(res) # description of column names
  structure(res, class = c(class(res), "get_diff"))
}


#' @rdname get_diff
#' @export
get_diff_n <- function(dds, factor, grp1, grp2 = NULL, alpha = 0.05, fc = 1.5,
                        active_exprs = 1, count.norm = NULL, exprs.norm = NULL, merge = FALSE,
                       verbose = TRUE, ...)
{

  if (!missing(count.norm)) {
    warning("argument count.norm is deprecated; please use exprs.norm instead.",
            call. = FALSE)
    exprs.norm <- count.norm
  }

  grps <- .check_groups(dds, factor, grp1, grp2)
  grp1 <- grps$grp1
  grp2 <- grps$grp2

  if(length(grp1) == 1 & length(grp2) == 1)
    res <- get_diff_2(dds, factor = factor, grp1 = grp1, grp2 = grp2,
                      alpha = alpha,  fc = fc, active_exprs = active_exprs,
                      exprs.norm = exprs.norm, verbose = verbose, ...)
  else if(merge){
    res <- .results(dds, alpha = alpha,
                    contrast = list(paste0(factor, grp1), paste0(factor, grp2) ), ...)

    res$detection_call <- .get_detection_call(dds, factor, grp1, grp2,
                                              active_exprs, exprs.norm )[rownames(res)]
    res$significance <- .get_significance(res, alpha = alpha, fc = fc)
    res <- structure(res, class = c(class(res), "get_diff"))
  }
  else if(!merge) res <- get_diff_pwc (dds, factor = factor, grp1 = grp1, grp2 = grp2,
                        alpha = alpha, fc = fc, active_exprs = active_exprs,
                        exprs.norm = exprs.norm, verbose = verbose, ...)

  res <- .set_mcols(res) # description of column names
  res
}


#' @rdname get_diff
#' @export
get_diff_pwc <- function(dds, factor, grp1, grp2 = NULL, alpha = 0.05, fc = 1.5,
                         active_exprs = 1, count.norm = NULL, exprs.norm = NULL, merge = FALSE,
                         verbose = TRUE, ...)
{

  if (!missing(count.norm)) {
    warning("argument count.norm is deprecated; please use exprs.norm instead.",
            call. = FALSE)
    exprs.norm <- count.norm
  }

  grps <- .check_groups(dds, factor, grp1, grp2)
  grp1 <- grps$grp1
  grp2 <- grps$grp2

  if(merge)
    res <- get_diff_n(dds, factor, grp1, grp2, alpha = alpha, fc = fc,
                      active_exprs = active_exprs,
                      exprs.norm = exprs.norm, merge = TRUE, verbose = verbose, ...)

  else if(length(grp1) == 1){
    if(length(grp2) == 1)
      res <- get_diff_2(dds, factor = factor, grp1 = grp1, grp2 = grp2,
                        alpha = alpha, fc = fc, active_exprs = active_exprs,
                        exprs.norm = exprs.norm, verbose = verbose, ...)
    else if(is.null(grp2) | length(grp2) >= 2)
      res <- get_diff_pwc(dds, factor = factor, grp1 = grp1, grp2 = grp2,
                        alpha = alpha, fc = fc, active_exprs = active_exprs,
                        exprs.norm = exprs.norm, verbose = verbose, ...)
  }
  else{
    ids <- rownames(dds)
    # Global comparison to have global base mean, log2fc, detection call and other statistics
    # the adjusted pvalue is replaced
    res <- get_diff_n(dds, factor, grp1, grp2, alpha = alpha,
                      active_exprs = active_exprs, fc = fc,
                      exprs.norm = exprs.norm, merge = TRUE, verbose = verbose, ...)
    res <- res[ids, , drop = FALSE]

    # Pairwise comparison: Compare each element in grp1 to each of elements in grp2
    pwc <- .pairewise_compare2(dds = dds, factor = factor, grp1 = grp1, grp2 = grp2,
                               alpha = alpha, verbose = verbose, ...)

    # changing padj and log2foldchange
    padj <- apply(pwc$padj, 1, max )
    status <- .deregulation_status(pwc$fcs, fc = 1) # up = 1, down = -1 and ns = 0
    log2fcs <- apply(abs(pwc$fcs), 1, min) * status
    padj[which(status == 0)] = NA

    res$padj <- padj[rownames(res)]
    res$log2FoldChange <- log2fcs
    res$significance <- .get_significance(res, alpha = alpha, fc = fc)
    res <- structure(res, class = c(class(res), "get_diff"))
  }
  res <- .set_mcols(res) # description of column names
  res
}




#' @rdname get_diff
#' @param object an object of class get_diff, DE_Results, DESeqResults or a data frame containing the
#' columns padj, log2FoldChange and optionally the column detection_call.
#' @param fdr Accepted false discovery rate for considering genes as differentially expressed.
#' @param fc the fold change threshold. Only genes with an absolute fold change >= fc and adjusted p-value <= fdr
#' and a detection_call = 1, are considered as significantly differentially expressed.
#' @param orderby possible values are "padj" and "lfc" for ordering
#' significant genes by the adjusted p-value or the log2FoldChange
#' @export
significant <- function(object, fdr = 0.05, fc = 1.5,
                        orderby = c("padj", "lfc"))
{

  if(!inherits(object, c("get_diff", "DESeqResults", "DE_Results")))
    stop("The argument object must be an object of class get_diff, DE_Results or DESeqResults.")
  object <- as.data.frame(object)

  object$significance <- .get_significance(object, alpha = fdr, fc = fc)
  signifs <- subset(object, object$significance!=0)

  if(nrow(signifs)>0) {
    orderby <- orderby[1]
    if(orderby == "lfc") signifs <- signifs[order(signifs$log2FoldChange, decreasing = TRUE), , drop = FALSE]
    else signifs <- signifs[order(signifs$padj, decreasing = FALSE), , drop = FALSE]
  }

  up <- subset(signifs, signifs$log2FoldChange > 0 )
  down <- subset(signifs, signifs$log2FoldChange < 0 )

  res <- list(up = up, down = down, nup = nrow(up), ndown = nrow(down))

  class(res) <- c("significant", class(res))
  res
}


# Mark significant genes in res
#::::::::::::::::::::::::::::::::::::::::::::
# res: an object of class DESeqDataSet
# log2.foldchange and pvalue cutoff are specified
mark_significant <- function(res, log2.foldchange = 0, pvalue = 0.05){
  # Check if lfc is above cutoff
  # Consider lfc sign to decide up and down genes
  lfc.ok <- (abs(res$log2FoldChange)) >= log2.foldchange
  lfc.sign <- sign(res$log2FoldChange)
  # Check if pvalue is below 0.05
  pval.ok <- ifelse(
    is.na(res$padj), FALSE,
    res$padj <= pvalue
  )
  # Detect significant genes
  if(!is.null(res$detection_call))
    gn.significant <- (lfc.ok & pval.ok) %>%
    as.numeric()
  else
  gn.significant <- (lfc.ok & pval.ok & res$detection_call == 1) %>%
    as.numeric()

  # Significance
  # 1 = up, -1 = down & 0 = NS
  significance <- gn.significant*lfc.sign
  res$significance <- significance
  res
}



