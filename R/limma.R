#' @include utilities.R expressed.R utilities_diffexprs.R AllClasses.R
NULL
#' Differential expression using limma
#'
#' @description Differential expression using limma:
#' \itemize{
#' \item \bold{lm_fit()}: Fit easily linear model for gene expression data using limma::lmFit().
#' \item \bold{lm_results()}: Extract results from lm_fit analysis. Compared to the function topTable(),
#' it performs automatically an independent filtering as described by Bourgon et al., PNAS 2010.
#' \item \bold{results_names()}: Returns the names of the estimated effects (coefficents) of the model
#' }
#'
#' @return
#' \itemize{
#' \item \bold{lm_fit()}: Returns an object of class "MArrayLM" (see limma:?lmFit) and "lm_fit".
#' Comprared to lmFit(), the object returned by lm_fit() contains supplementary variables named:
#' \itemize{
#' \item design_formula: which holds the design formula.
#' \item samples: annotation of each samples as provided in eset.
#' \item data: unlogged gene expression data
#' }
#' \item \bold{lm_results()}: Returns an object of class "DE_Results" which is a simple subclass of DataFrame (in S4Vectors)
#' containing log2 fold changes, base means across samples, t-test statistics, p-values
#' adjusted p-values. The function plot_MA() can be used to draw the results. (See ?plot_MA).
#' \item \bold{results_names()}: Returns the names of the estimated effects (coefficents) of the model.
#' }
#' @examples
#' \donttest{
#' # Load and prepare data
#' data(ALL_bcrneg)
#' ALL_bcrneg$mol.biol <- factor(make.names(ALL_bcrneg$mol.biol),
#' levels = c("NEG", "BCR.ABL"))
#'
#' # Fit linear model
#' fit <- lm_fit(ALL_bcrneg, design = ~0+mol.biol, logged2 = TRUE)
#'
#' # result names
#' results_names(fit)
#'
#' # Results of differential expression
#' res <- lm_results(fit)
#'
#' # Plot MA
#' plot_MA(res)
#'
#' # Independent filtering threshold and rejection plot
#'  attr(res,"filterThreshold")
#'
#'  plot(attr(res,"filterNumRej"),type="b",
#'  xlab="quantiles of 'baseMean'",
#'  ylab="number of rejections")
#'
#' }
#' @name limma
#' @rdname limma
#' @param eset an object of class ExpressionSet or a list with two components including data and samples.
#' eset must contain the information about samples.
#' @param design a formula which expresses how the expression for each gene depend on the variables in the sample data.
#' e.g.: ~ 0 + group (design with multiple groups ),  ~ group + condition (design with multiple variables)
#' and ~ genotype + treatment + genotype:treatment (design with interactions).
#' @param logged2 logical value. Specify if the data are in log2 scale or not. Allowed values are TRUE or FALSE.
#' @export
lm_fit <- function(eset, design, logged2){

  if(!inherits(eset, c("ExpressionSet", "list")))
    stop("eset must be an object of class ExpressionSet or a list with two components including data and samples.")
  # Data and samples
  if(inherits(eset, "ExpressionSet")){
    data <- Biobase::exprs(eset)
    samples <- Biobase::pData(eset)
  }
  else{
    data <- eset[[1]]
    samples <- eset[[2]]
  }
  if (any(is.na(data)))
    stop("NA values are not allowed in the expression matrix")
  if (any(duplicated(rownames(data))))
    stop("Duplicated row names in the data are not allowed.")

  if(!logged2) {
    unlogged_data <- data
    data <- log2(data+0.5)
  }
  else unlogged_data <- 2^data

  if(inherits(samples, "matrix")) samples <- as.data.frame(samples)

  designVars <- all.vars(design)
  if (!all(designVars %in% names(samples))) {
    stop("all variables in design formula must be columns in samples data")
  }

  designVarsClass <- sapply(designVars, function(v) class(samples[[v]]))
  if (any(designVarsClass != "factor")) {
    warning("some variables in design formula are not a 'factor', converting to factors")
    for (v in designVars[designVarsClass != "factor"]) {
      samples[[v]] <- factor(samples[[v]])
    }
  }
  designVarsClass <- sapply(designVars, function(v) class(samples[[v]]))

  if (length(designVars) == 1) {
    var <- samples[[designVars]]
    if (all(var == var[1])) {
      stop("design has a single variable, with all samples having the same value.\n")
    }
  }

  designFactors <- designVars[designVarsClass == "factor"]
  missingLevels <- sapply(designFactors, function(v) any(table(samples[[v]]) == 0))
  if (any(missingLevels)) {
    message("factor levels were dropped which had no samples")
    for (v in designFactors[missingLevels]) {
      samples[[v]] <- droplevels(samples[[v]])
    }
  }

  modelMatrix <- stats::model.matrix(design, data = samples)
  fit <- limma::lmFit(data, modelMatrix)
  fit$design_formula <- design
  fit$samples <- samples
  fit$data <- unlogged_data
  structure(fit, class = c(class(fit), "lm_fit"))
}


#' @param object an object of class "lm_fit" generated from the function lm_fit()
#' @param contrast specifies what comparison to extract from the object to build a results table.
#' Possible values include:
#' \itemize{
#' \item a character vector with three elements:
#' the name of a factor in the design formula, the name of group 1 and group 2 to be compared;
#' e.g.: c("group", "grp1", "grp2"). "grp1" is used as numerator and "grp2" as denominator for
#' computing the fold change.
#' \item a list of 2 character vectors to be compared. These names should be elements of lm_results_names(object).
#' e.g.: list(c("grp1","grp2"), c("grp3", "group4"))
#' }
#' @param independentFiltering logical, whether independent filtering should be applied automatically
#' @param alpha the significance cutoff used for optimizing the independent filtering (by default 0.1).
#' @param filter the vector of filter statistics over which the independent filtering will be optimized.
#' @param theta the quantiles at which to assess the number of rejections from independent filtering
#' @param pAdjustMethod method used to adjust the p-values for multiple testing.
#' Allowed values include "none", "BH", "BY" and "holm" (See ?p.adjust).
#' @param ... other arguments to be passed to the function limma::topTable
#' @rdname limma
#' @export
lm_results <- function (object,  contrast,
                         independentFiltering = TRUE, alpha = 0.05,
                         filter, theta, pAdjustMethod="BH", ...)
{
  stopifnot(length(alpha) == 1)
  stopifnot(length(pAdjustMethod) == 1)

  if(missing(contrast)){
    designVars <- all.vars(object$design_formula)
    lastVarName <- designVars[length(designVars)]
    lastVar <- object$samples[[lastVarName]]
    if (is.factor(lastVar)) {
      nlvls <- nlevels(lastVar)
      contrast <- c(lastVarName, levels(lastVar)[nlvls],
                    levels(lastVar)[1])
    }
  }

  resNames <- results_names(object)
  contrast <- .check_contrast(contrast, resNames)

  # comparison
  grps <- paste0("(", contrast[1], ")", "-", "(", contrast[2], ")")
  contrast.matrix <- limma::makeContrasts(contrasts = grps, levels = object$design)
  fit <- limma::contrasts.fit(object, contrast.matrix)
  fit <- limma::eBayes(fit)

  # top table
  res <- limma::topTable(fit, adjust.method=pAdjustMethod, sort.by="p", number = length(fit$p.value), ...)
  res$AveExpr <- 2^res$AveExpr
  # Renaming some column names
  nn <- colnames(res)
  nn[which(nn == "logFC")] <- "log2FoldChange"
  nn[which(nn == "AveExpr")] <- "baseMean"
  nn[which(nn == "P.Value")] <- "pvalue"
  nn[which(nn == "adj.P.Val")] <- "padj"
  colnames(res) <- nn

  paRes <- .pvalue_adjustment(res, independentFiltering, filter,
                              theta, alpha, pAdjustMethod)
  res$padj <- paRes$padj
  res <- S4Vectors::DataFrame(res)
  res <- .lm_set_mcols(res)
  res <- DE_Results(res)

  if (independentFiltering) {
    attr(res, "filterThreshold") <- paRes$filterThreshold
    attr(res, "filterNumRej") <- paRes$filterNumRej
  }
  res
}


#' @rdname limma
#' @export
results_names <- function (object){
  if(inherits(object, "MArrayLM")) colnames(object$design)
  else if(inherits(object, "DESeqDataSet")) DESeq2::resultsNames(object)
  else stop("Can't handle an object of class ", class(object))
}
