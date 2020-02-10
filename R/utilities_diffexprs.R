# Utilities for differential expression

# Check contrasts
# resNames : as returned by results_names()
.check_contrast <- function (contrast, resNames)
{
  if (!(is.numeric(contrast) | !is.character(contrast) | !is.list(contrast))) {
    stop("'contrast' vector should be either a character vector of length 3,\na list of length 2 containing character vectors,\nor a numeric vector, see the argument description in ?results")
  }
  if (is.character(contrast)) {
    if (length(contrast) != 3) {
      stop("'contrast', as a character vector of length 3, should have the form:\ncontrast = c('factorName','numeratorLevel','denominatorLevel'),\nsee the manual page of ?results for more information")
    }
    if (contrast[2] == contrast[3]) {
      stop(paste(contrast[2], "and", contrast[3], "should be different level names"))
    }

    contrast <- c(paste0(contrast[1], contrast[2]), paste0(contrast[1], contrast[3]))
  }

  if (is.list(contrast)) {
    if (length(contrast) == 1) contrast <- list(contrast[[1]], character())
    if (length(contrast) != 2) {
      stop("'contrast', as a list, should have length 2,\nor, if length 1, an empty vector will be added for the second element.\nsee the manual page of ?results for more information")
    }
    if (!(is.character(contrast[[1]]) & is.character(contrast[[2]]))) {
      stop("'contrast', as a list of length 2, should have character vectors as elements,\nsee the manual page of ?lm_results for more information")
    }
    if (!all(c(contrast[[1]], contrast[[2]]) %in% resNames)) {
      stop("all elements of the contrast as a list of length 2 should be elements of 'results_names(object)'")
    }
    if (length(intersect(contrast[[1]], contrast[[2]])) >
        0) {
      stop("elements in the contrast list should only appear in the numerator (first element of contrast list)\nor the denominator (second element of contrast list), but not both")
    }
    if (length(c(contrast[[1]], contrast[[2]])) == 0) {
      stop("one of the two elements in the list should be a character vector of non-zero length")
    }
    contrast<- c(paste(contrast[[1]], collapse="+"), paste(contrast[[2]], collapse="+"))
  }

  if (is.numeric(contrast)) {
    stop("Can't handle numeric contrast")
  }
  return(contrast)
}

# function for independent filtering
#++++++++++++++++++++++++++++
# res: an object of class lm_results or DSeqResults
# independentFiltering: logical value. If true independent filtering is applied
# for the others arguments: see ?lm_results
.pvalue_adjustment<- function (res, independentFiltering, filter, theta, alpha, pAdjustMethod){
  if (independentFiltering) {
    if (missing(filter)) {
      filter <- res$baseMean
    }
    if (missing(theta)) {
      lowerQuantile <- mean(filter == 0)
      if (lowerQuantile < 0.95)
        upperQuantile <- 0.95
      else upperQuantile <- 1
      theta <- seq(lowerQuantile, upperQuantile, length = 20)
    }
    stopifnot(length(theta) > 1)
    stopifnot(length(filter) == nrow(res))
    filtPadj <- genefilter::filtered_p(filter = filter, test = res$pvalue,
                                       theta = theta, method = pAdjustMethod)
    numRej <- colSums(filtPadj < alpha, na.rm = TRUE)
    j <- which.max(numRej)
    padj <- filtPadj[, j, drop = TRUE]
    cutoffs <- quantile(filter, theta)
    filterThreshold <- cutoffs[j]
    filterNumRej <- data.frame(theta = theta, numRej = numRej)
    return(list(padj = padj, filterThreshold = filterThreshold,
                filterNumRej = filterNumRej))
  }
  else {
    padj <- p.adjust(res$pvalue, method = pAdjustMethod)
    return(list(padj = padj))
  }
}


# Get the results from lm_fit or DEseq
# dds : lm_fit or DESeqDataSet
# contrast: c(factor, grp1, grp2)
.results <- function(dds, contrast, alpha = 0.05, ...){
  if(inherits(dds, "DESeqDataSet"))
    res <- DESeq2::results(dds, contrast = contrast, alpha = alpha, ...)
  else if("lm_fit" %in% class(dds))
    res <- lm_results(dds, contrast = contrast, alpha = alpha, ...)
  else stop("Can't handle an object of class ", class(dds))
  res
}

# Get sample annotation
.samples <- function(dds){
  if(inherits(dds, "DESeqDataSet"))
    res <- as.data.frame(SummarizedExperiment::colData(dds))
  else if("lm_fit" %in% class(dds))
    res <- as.data.frame(dds$samples)
  else stop("Can't handle an object of class ", class(dds))
  res
}

# Get the normalized data
.data_norm <- function(dds){
  if(inherits(dds, "DESeqDataSet"))
    res<- as.data.frame(round(DESeq2::counts(dds, normalized=TRUE),2))
  else if("lm_fit" %in% class(dds))
    res <- dds$data
  else stop("Can't handle an object of class ", class(dds))
  res
}

# Data into log scale.
# if DESeq --> rlog
# else log2()
.rlog_data <- function(dds){
  if(inherits(dds, "DESeqDataSet"))
    res <- as.data.frame(SummarizedExperiment::assay(DESeq2::rlog(dds)))
  else if("lm_fit" %in% class(dds))
    res <- log2(dds$data)
  else stop("Can't handle an object of class ", class(dds))
  res
}


# Get  available groups
.available_groups <- function(dds, factor){

  if(inherits(dds, "DESeqDataSet"))
    available.grps <- levels(SummarizedExperiment::colData(dds)[, factor])
  else if("lm_fit" %in% class(dds))
    available.grps <- levels(dds$samples[, factor])
  else stop("dds must be an object of class DESeqDataSet or lm_fit.")

  available.grps
}

# Check dds and if groups exist in dds
# If grp2 is null, the remaining groups in dds (exluding grp1) are used as grp2
# returns a list containing the elements grp1 and grp2
.check_groups <- function(dds, factor, grp1, grp2){

  available.grps <- .available_groups(dds, factor)

  if(is.null(grp2)){
    grp2 <- available.grps
    grp2 <-setdiff(grp2, grp1)
  }
  dif <- setdiff(c(grp1, grp2), available.grps )
  if(length(dif) > 0) stop("Can't found the following group in the data: ",
                           paste(dif, collapse = ", "), ". ",
                           "Available groups are: ", paste(available.grps, collapse = ", "), ".")
  if(length(grp2)==0)
    stop("The data contains only one group: ", paste(grp1, collapse = ", "))

  list(grp1 = grp1, grp2 = grp2)
}




# Pairwise comparisons between grp1 and each of the elements in grp2
# Compare grp1 to each of the elements in grp2
# grp1 and grp2 must be as follow: grp1 = c("A"), grp2 = c("A", "B", "C", ...)
# returns a list(fcs, padj, significance)
.pairewise_compare <- function(dds, factor, grp1, grp2, alpha = 0.05, fc = 1,  verbose = TRUE, ...)
{

  ids <- rownames(dds)
  fcs <- padj <- significance <-data.frame(id = rownames(dds))
  rownames(significance) <- rownames(fcs) <- rownames(padj) <-  ids

  # Pairwise comparison
  for(i in 1:length(grp2)){
    step_name <- tolower(paste0(grp2[i], "_vs_", grp1))
    if(verbose) message("- Comparison ", step_name, ": Progress = ", i, "/", length(grp2))
    res <- as.data.frame(.results(dds, contrast = c(factor, grp1, grp2[i]), alpha = alpha, ...))[ids, ]
    fcs <- cbind(fcs, res[ids, "log2FoldChange", drop = FALSE])
    padj <- cbind(padj, res[ids, "padj", drop = FALSE])
    # significance
    if (fc == 1) fcsOk <- abs(res$log2FoldChange) > log2(fc) else fcsOk <- abs(res$log2FoldChange) >= log2(fc)
    sig <- as.numeric(ifelse(is.na(res$padj), FALSE, res$padj <= alpha & fcsOk))*sign(res$log2FoldChange)
    names(sig) <- ids
    significance  <- cbind(significance, sig = sig)
    colnames(fcs)[ncol(fcs)] <- step_name
    colnames(padj)[ncol(padj)] <- step_name
    colnames(significance)[ncol(significance)] <- step_name
  }
  padj <- as.data.frame(padj)[ids, -1, drop = FALSE]
  fcs <- as.data.frame(fcs)[ids, -1, drop = FALSE]
  significance <- as.data.frame(significance)[ids, -1, drop = FALSE]
  list(padj = padj, fcs = fcs, significance = significance)
}

# Pairwise comparisons between each element of grp1 and each of the elements in grp2
# grp1 and grp2 must be as follow: grp1 = c("A", "B", ...), grp2 = c("C", "D",  ...)
.pairewise_compare2 <- function(dds, factor, grp1, grp2, alpha = 0.05, fc = 1,  verbose = TRUE, ...){
  ids <- rownames(dds)
  pwc <- list(fcs = data.frame(id = ids), padj = data.frame(id = ids), significance =data.frame(id = ids) )
  for(i in 1:length(grp1)){
    # Pairwise comparison: Compare grp1 to each of elements in grp2
    pwci <- .pairewise_compare(dds = dds, factor = factor, grp1 = grp1[i], grp2 = grp2,
                               alpha = alpha, fc = fc, verbose = verbose, ...)
    pwc$fcs <- cbind(pwc$fcs, pwci$fcs[ids, , drop = FALSE])
    pwc$padj <- cbind(pwc$padj, pwci$padj[ids, , drop = FALSE])
    pwc$significance <- cbind(pwc$significance, pwci$significance[ids, , drop = FALSE])
  }
  pwc$fcs <- pwc$fcs[, -1, drop = FALSE]
  pwc$padj <- pwc$padj[, -1, drop = FALSE]
  pwc$significance <- pwc$significance[, -1, drop = FALSE]
  pwc
}




# fcs a data frame containing fold changes of different comparisons
# fc fold change cutoff
# the result is a vector with the following values : -1 (if downregulated in all columns);
# 1 (if upregulated in all columns); 0 if not changed
.deregulation_status <- function(fcs, fc = 1.5){
  if(fc < 1) stop("fc must be >= 1.")
  # Sign to identify if up or down
  fcs.sign <- apply(fcs, 2, sign)
  fcs.sign <- apply(fcs.sign, 1,
                    function(x, ncol){
                      s <- sum(x)
                      if(s == ncol) 1
                      else if(s == -ncol) -1
                      else 0
                    },
                    ncol(fcs))
  # Levels of deregulation
  fcs.lev <- apply(abs(fcs), 1, min)
  if(fc ==1) fcs.lev <- fcs.lev > 0 else fcs.lev <- fcs.lev >= log2(fc)
  fcs.sign*fcs.lev
}

# Tell us if a gene is significant or not
# res is a result of DESeq2
.get_significance <- function(res, alpha = 0.05, fc = 1.5){

  if(fc == 1) fcsOk <- abs(res$log2FoldChange) > log2(fc) else fcsOk <- abs(res$log2FoldChange) >= log2(fc)
  if("detection_call" %in% colnames(res)) detection_call <- res$detection_call else detection_call <- rep(1, nrow(res))
  significance <- as.numeric(ifelse(is.na(res$padj), FALSE, res$padj <= alpha & fcsOk &
                                      detection_call == 1))*sign(res$log2FoldChange)
  significance
}

# detection call
# dds an object of class DESeqResults as provided by DESeq2
# factor, grp1, grp2, active_exprs, count.norm: see ?get_diff
# exprs.norm an alias of count.norm
.get_detection_call <- function(dds, factor, grp1, grp2 = NULL,
                                active_exprs = 1, count.norm = NULL, exprs.norm = NULL){

  if(is.null(count.norm))  count.norm <- exprs.norm

  samples <- .samples(dds)
  detection_call <- rep(1, nrow(dds))
  names(detection_call) <- rownames(dds)
  if(!is.null(active_exprs)){
    # Detection call: Expressed genes in at least one groups
    if(is.null(count.norm))
      count.norm <- .data_norm(dds)

    sples <- subset(samples, samples[, factor] %in% c(grp1, grp2) )
    grps <- droplevels(sples[, factor])
    exprs <- count.norm[, rownames(sples), drop = FALSE]

    detection_call <- apply(expressed(exprs,  grps, cutoff = active_exprs), 1, max)[rownames(dds)]
    detection_call <- detection_call
  }
  return(detection_call)
}

# Description of column names
.set_mcols <- function(res){

  S4Vectors::mcols(res)$type[names(res) == "detection_call"] <- "results"
  S4Vectors::mcols(res)$description[names(res) == "detection_call"] <-
    "Detection call: Active expression = 1; Low expression = 0"
  S4Vectors::mcols(res)$type[names(res) == "significance"] <- "results"
  S4Vectors::mcols(res)$description[names(res) == "significance"] <-
    "Deregulation significance: Up = 1; Down = -1; Not significant = 0"
  res
}

# used in lm_results()
.lm_set_mcols <- function(res){
  tpes <- rep("results", 6 )
  description <- c("Estimate of the log2-fold-change",
                   "Average log2-expression for the probe over all samples",
                   "moderated t-statistic", "raw p-value", "adjusted p-value or q-value",
                   "log-odds that the gene is differentially expressed")
  S4Vectors::mcols(res) <- S4Vectors::DataFrame(type = tpes, description = description,
                                                row.names = c("log2FoldChange", "baseMean", "t", "pvalue", "padj", "B"))
  res
}


