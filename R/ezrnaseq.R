#' @include utilities.R
#' @importFrom  DESeq2 DESeqDataSetFromMatrix counts DESeq
#' @importFrom tibble rownames_to_column as_tibble
NULL
#' Import RNAseq Data
#'
#' @description Imports RNAseq data contained in a directory.
#' @param directory path to the directory containing the RNAseq data. This directory
#'   should contain the following files:
#'   \itemize{
#'   \item{count.txt}: containsraw count data
#'   \item{samples.txt}: contains sample groups
#'   \item {annotation.txt}: contains the gene annotation file.
#'   }
#' @describeIn ezrnaseq Import RNAseq data
#' @export
import_rnaseq <- function(directory = "data"){
  options(readr.num_columns = 0)
  raw_count <- readr::read_tsv(file.path(directory, "count.txt"))
  samples <- readr::read_tsv(file.path(directory, "samples.txt"))
  gene_annotation <- readr::read_tsv(file.path(directory, "annotation.txt"))
  res <- list(
    count = raw_count,
    annotation = gene_annotation,
    samples = samples,
    directory = directory
  )
  res$samples$group <- factor(res$samples$group, levels = unique(res$samples$group))
  res
}


#' @describeIn ezrnaseq Compare two groups
#' @param data a list containing gene expression data plus samples and gene annotation.
#' @param grp1 group 1 used as numerator for calculating the fold change.
#' @param grp2 group 2 used as denominator for calculating the fold change.
#' @param fold_change fold change to be considered
#' @export
compare_two_groups <- function(data, grp1, grp2, fold_change = 2){

  # Create DESeq dataset
  dds <- DESeqDataSetFromMatrix(
    countData = .as_dseq_input(data$count),
    colData = .as_dseq_input(data$samples),
    design = ~group
    )
  gene_annotation <- .as_dseq_input(data$annotation)
  # Minimal pre-filtering to remove rows that have 0 read.
 dds <- dds[ rowSums(counts(dds)) > 0, ]
 # Differential expression
 BiocParallel::register( BiocParallel::MulticoreParam() )
 dds <- DESeq(dds, parallel = TRUE)
 # Normalize count
 count.norm <- as.data.frame(round(counts(dds, normalized=TRUE),2))
 # Significant results
 comparison <- get_diff_2(dds, factor = "group", grp1 = grp1,  grp2 = grp2,
                   active_exprs = 64, fc = fold_change,
                   exprs.norm = count.norm)

 res.signif <- significant(comparison, fdr = 0.05, fc = fold_change)
 res <- NULL
 if(res.signif$nup!=0) res <- rbind(res, res.signif$up)
 if(res.signif$ndown!=0) res <- rbind(res, res.signif$down)

 # Annotation of significant results; and order by p-values
 if(!is.null(res)){
   res <- cbind(gene_annotation[rownames(res),   "name", drop = FALSE],
                res[, c("baseMean", "log2FoldChange", "padj")],
                count.norm[rownames(res), ])
   res <- res[order(res$padj), ] %>%
     rownames_to_column("id") %>%
     as_tibble()
   colnames(res)[1] <- colnames(data$count)[1]

 }

 # Annotation of all the comparisons
 comparison <- as.data.frame(comparison, stringsAsFactors = FALSE)
 comparison <- cbind(gene_annotation[rownames(comparison),   "name", drop = FALSE],
               comparison[, c("baseMean", "log2FoldChange", "padj")],
              count.norm[rownames(comparison), ])
 comparison <- comparison[order(comparison$padj), ] %>%
   rownames_to_column("id") %>%
   as_tibble()
 colnames(comparison)[1] <- colnames(data$count)[1]

 results <- list(
   comparison = paste0(grp1, "_vs_", grp2),
   significant = res,
   all_results = comparison
 ) %>%
   c(data)
 results
}

.as_dseq_input <- function(data){
  data <- data %>%
    as.data.frame(stringsAsFactors = FALSE)
  rownames(data) <- data[[1]]
  data <- data[, -1, drop = FALSE]
  data
}


#' @describeIn ezrnaseq Export significant results
#' @param results results returned by compare_groups
#' @export
export_significant <- function(results){
  result_file = paste0("results_", results$comparison, ".txt")
  readr::write_tsv(results$significant, file.path(results$directory, result_file))
}
