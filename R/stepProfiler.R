#' @include utilities.R get_diff.R utilities_stepProfiler.R summarize.R
NULL
#' StepProfiler: Extract temporally Up- and Down- regulated genes using DESeq2.
#'
#' @description
#' stepProfiler can be used to identify four temporal gene expression patterns that we refer to as
#' one-step-up, one-step-down, impulsed-up, impulsed-down.
#' \itemize{
#' \item pairewise_compare(): paire-wise comparison between all groups.
#' \item one_step_up(): Low before a step i and increase gradually, at step i+1, during all the next transition steps.
#' \item one_step_down(): High before a step i and  decrease gradually, from step i+1,  during all the next transition steps.
#' \item one_step_up_sustained(): Low before a step i, upregulated, from step i+1, and remain stable during all the next transition steps:
#' the maximum fold change between the step i+1 and the following steps must be in [-1.4, 1.4] (stability fold change).
#' \item one_step_down_sustained(): Highly expressed before step i, downregulated at step i+1 and remain stable
#' during all the next transition steps: the maximum fold change between the current step and the next step must be in [-1.4, 1.4].
#' \item impulses_up(): upregulated gradually and significantly during all the previous transitions (i--) until the current step i,
#' then decrease significantly and gradually during all the next transition steps. There is an impulse signal at
#' step i with a minimum fc = 1.5. This expression profile is also called as two-step-down-up.
#' \item impulse_up_sustained(): Low and stable during the previous and the next transitions,
#' significant impulse at step i with a minimum fc = 1.5.
#' \item impulse_down(): downregulated significantly and gradually during all the previous transitions (i--)
#' until the current step i, then increase significantly and gradually during all the next transition steps.
#' There is an impulse signal at step i with a minimum fc = 1.5. This expression profile is also called as two-step-down-up.
#' \item impulse_down_sustained(): High and stable during the previous and the next transitions,
#' significant impulse at step i with a minimum fc = 1.5.
#' \item plot(): Draw step genes using ggplot2
#' }
#' @param dds an object of class DESeqDataSet containing the diff_exprs of the function DESeq().
#' @param factor the name of a factor in the design formula (e.g factor = "group")
#' @param grps a vector containing at list 2 groups to be compared.
#' @param alpha false discovery rate
#' @param fc positive value specifying fold change threshold. Must be at least 1.
#' Genes with fc > 1 are kept.
#' @param active_exprs a numeric value specifying detected expression cutoff.
#' If not NULL, any gene with an average expression above this cutoff, in at least one group,
#' is declared as actively expressed. Used for defining presence and absence detection call.
#' If not NULL, a column named detection_call is added to the result. Possible values for detection_call
#' is 0 (the gene is not actively expressed) or 1 (the gene is actively expressed).
#' @param count.norm a data frame containing the normalized count as provided by DESeq2.
#' Only used for computing the detection call when active_expression is not NULL.
#' @param exprs.norm a data frame containing the normalized gene expression data. An alias of count.norm.
#' @param verbose logical value; if TRUE the progression is shown.
#' @param pwc an object of class "pairewise_compare" as returned by the function pairewise_compare().
#' @return
#' \bold{pairewise_compare()} returns a list containing the following elements:
#' \itemize{
#'
#' \item significance: a data frame containing the pairewise comparison results.
#' Each column corresponds to the significance of a given comparison. Rows are genes wich can have -1, 0 or 1 values for
#' each column. -1 = significantly downregulated; 0 = not significant; 1 = significantly upregulated.
#' A significant deregulated gene is a gene with active expression AND padj < = alpha AND absolute fold change >= fc.
#' \item foldchanges: a data frame containing log2 fold changes from the different comparisons.
#' \item padj: a data frame containing the adjusted p-values from the different comparisons.
#' \item groups:  the differents groups used for the pairewise comparisons.
#' This corresponds to levels(samples$group).
#' \item steps: transition steps from group_i to group_i+1
#' \item nstep: the number of steps
#' \item steps_mean_fc: a numeric vector containing the mean of the absolute fold change during transitions.
#' Can be useful for sorting genes.
#' \item dds: the dds used during pairewise comparison.
#' \item count.norm: the normalized count as provided by DESeq2
#' \item rlog.data: 'regularized log' transformed data as provided by DESeq2
#' \item samples: samples data as provided by colData(dds)
#' \item refs: a numeric vector containing the mean normalized expression of each gene in the first group.
#' This is used when plotting the results with argument transformby = "firstgroup" in order to display the
#' log2 fold change expression in the other groups compared to the first group as reference samples.
#' \item detection_call a vector containing 1 or 0 specifying whether the gene is actively expressed or not in at
#' least one group.
#' This computed according to the value of the argument active_exprs.
#' \item steps_mean_fc the mean of the absolute fold change during transitions. Can be useful for sorting genes.
#' \item baseMean the mean normalized count for all samples.
#' }
#' \bold{one_step_*(), impulsed_*()} returns an object of class stepgenes which is a list of sub-lists.
#' Each sub-list corresponds to a given step and contains the following elements:
#' \itemize{
#' \item exprs: a data frame containing the normalized gene expression values.
#' Columns are samples and rows are genes.
#' \item foldchanges, padj, refs, baseMean: see the description above.
#' \item group: a factor containing the group assignment of samples
#' }
#' Genes are ordered according to the mean fold change expression during all transitions,
#' so that the most relevant genes are shown at the top the data table.
#'
#' \bold{plot()}: Plots step genes using ggplot2.\cr\cr
#'
#' \bold{summary()}: Returns a data frame which columns are:
#' \itemize{
#' \item id: gene IDs
#' \item status: the expression status of a gene at a given step
#' \item step: the step at which the expression status changes
#' \item group: Numeric coding of the step column
#' \item max_fcs and min_fcs: the maximum and the minimum fold change expression at each transition, respectively.
#' }
#'
#' @name stepProfiler
#' @rdname stepProfiler
#' @export
pairewise_compare <- function(dds, factor, grps, alpha = 0.05, fc = 1,
                            active_exprs = 1, count.norm = NULL, exprs.norm = NULL, verbose = TRUE, ...)
{
  if (!missing(count.norm)) {
    warning("argument count.norm is deprecated; please use exprs.norm instead.",
            call. = FALSE)
    exprs.norm <- count.norm
  }

  if(length(grps)==1) stop("The argument groups must contain at least 2 groups.")
  available.grps <- .available_groups(dds, factor)

  dif <- setdiff(grps, available.grps )
  if(length(dif) > 0) stop("Can't found the following group in the data: ",
                           paste(dif, collapse = ", "), ". ",
                           "Available groups are: ", paste(available.grps, collapse = ", "), ".")
  if(fc < 1) stop("fc must be >= 1")

  # rlog data
  if(verbose) message("- Computing rlog data")
   rlog.data <- .rlog_data(dds)

  # Normalized count
  if(verbose) message("- Computing normalized data")
  if(is.null(exprs.norm)) exprs.norm <- .data_norm(dds)+1
  baseMean <- apply(exprs.norm, 1, mean)
  samples <- .samples(dds)
  grp_levels <- c(grps, setdiff(unique(samples[, factor]), grps))
  samples[, factor] <- factor(samples[, factor], levels = grp_levels)

  # The first group is used as reference
  refs <- which(samples[, factor] == grps[1])
  refs <- exprs.norm[, refs, drop = FALSE]
  refs <- apply(refs, 1, mean)

  # Identify actively expressed genes
  detection_call <- .get_detection_call(dds, factor = factor, grp1 = grps, grp2 = NULL,
                                            active_exprs = active_exprs, exprs.norm = exprs.norm )[rownames(dds)]

  n.comparison <- length(grps)*2
  # Differentially expressed genes between a given step and the next step
  ids <- rownames(dds)
  significance <- fcs <- padj <- data.frame(id = rownames(dds))
  rownames(significance) <- rownames(fcs) <- rownames(padj) <-  ids
  for(i in 1:(length(grps)-1)){
    grp_i <- grps[i]
    grp_js <- grps[(i+1):length(grps)]
    pwc <- .pairewise_compare2(dds = dds, factor = factor, grp1 = grp_js, grp2 = grp_i,
                       alpha = alpha,  fc = fc, verbose = verbose, ...)

    pwc$significance$id <- rownames(pwc$significance)
    pwc$padj$id <- rownames(pwc$padj)
    pwc$fcs$id <- rownames(pwc$fcs)

    significance <- merge(significance, pwc$significance, by.x = "id", by.y = "id", all.x = TRUE, sort = FALSE)
    fcs <- merge(fcs, pwc$fcs, by.x = "id", by.y = "id", all.x = TRUE, sort = FALSE)
    padj <- merge(padj, pwc$padj, by.x = "id", by.y = "id", all.x = TRUE, sort = FALSE)
  }
  rownames(significance) <- rownames(padj) <- rownames(fcs) <- as.vector(fcs$id)

  if(verbose) message("- Computing significance adjusted by the detection call")
  significance <- significance[ids, -1, drop = FALSE]
  fcs <- fcs[ids, -1, drop = FALSE]
  fcs <- apply(fcs, 2, function(x, detection_call) as.vector(x)*detection_call, detection_call)
  padj <- padj[ids, -1, drop = FALSE]

  # Transition step
  steps <- NULL
  for(i in 1:(length(grps)-1)){
    steps <- c(steps, paste0(grps[i], "_vs_", grps[i+1]))
  }
  steps <- tolower(steps)
  nstep <- length(steps)


  # Mean fold change expression in all comparison
  mean_fc <- apply(fcs[, steps, drop = FALSE], 1, mean)
  res <- list(significance = significance, foldchanges = fcs, padj = padj, factor = factor, groups = grps,
              steps = steps, nstep = nstep, steps_mean_fc = mean_fc,
              dds = dds, count.norm = exprs.norm, rlog.data = rlog.data,
              samples = samples, refs = refs, detection_call = detection_call,
              baseMean = baseMean)

  structure(res, class = c("list", "pairewise_compare"))
}



#' @rdname stepProfiler
#' @export
one_step_up <- function(pwc, fc = 1.5, verbose = TRUE){

  # Upregulated at step i
  # Stable at the next step: fold change between the current step and the next step in -1.4-1.4

  if(!inherits(pwc, "pairewise_compare"))
    stop("pwc must be an object of class pairewise_compare.")

  all <- pwc$significance
  foldchanges <- pwc$foldchanges
  groups <- pwc$groups
  max_stability_fc = 1.4

  results <- list()
  # For the first transition 1
  # +++++++++++++++++++++++
  if(verbose) message("- Identifying genes for step ", pwc$groups[1])
  grp_i <- groups[2]
  current_trans <- apply(all[, pwc$steps[1], drop = FALSE], 1, sum)
  current_fcs <- abs(apply(foldchanges[, pwc$steps[1], drop = FALSE], 1, sum))
  # Stability in the next steps
  isStable_nextsteps <- .is_upSustained_in_nextsteps(pwc, from = 2,
                                                     max_stability_fc = max_stability_fc)
  # Stability in the previous steps
  isStable_prevsteps <- .is_stable_in_prevsteps(pwc, to = 2,
                                                max_stability_fc = max_stability_fc )
  # Result
  index <- which(current_trans == 1 & isStable_nextsteps == 1 & current_fcs >= log2(fc) & pwc$detection_call == 1)
  if(length(index) == 0) warnings("No genes for step ", grp_i)
  res <- .subset_index(pwc, index, sort_ = "decreasing")
  res$title <- paste0( "Up-at-", grp_i)
  results[[1]] <- res

  # For the remaining steps
  # ++++++++++++++++++++++++++
  for(i in 2:(length(groups)-1)){
    if(verbose) message("- Identifying genes for step ", pwc$groups[i])
    current_trans <- apply(all[, pwc$steps[i], drop = FALSE], 1, sum)
    current_fcs <- abs(apply(foldchanges[, pwc$steps[i], drop = FALSE], 1, sum))
    # Stability in the next steps
    isStable_nextsteps <- .is_upSustained_in_nextsteps(pwc, from = i+1)
    # Stability in the previous steps
    isStable_prevsteps <- .is_stable_in_prevsteps(pwc, to = i, max_stability_fc = max_stability_fc )
    # Result
    index <- which(current_trans ==1 &
                     isStable_nextsteps ==1 & isStable_prevsteps ==1 &
                     current_fcs >= log2(fc) & pwc$detection_call == 1)
    if(length(index) == 0) warnings("No genes for step ", groups[i])
    res <- .subset_index(pwc, index, sort_ = "decreasing")
    res$title <- paste0( "Up-at-", groups[i+1])
    results[[i]] <- res
  }

  names(results) <- pwc$groups[2:(length(pwc$groups))]
  structure(results, class = c("list", "stepgenes", "onestepup"))
}




#' @rdname stepProfiler
#' @export
one_step_down <- function(pwc, fc = 1.5, verbose = TRUE){

  # Highly expressed at previous step
  # downregulated at step i and remain stable at the next step: fold change between the current step and the next step in -1.4-1.4

  if(!inherits(pwc, "pairewise_compare"))
    stop("pwc must be an object of class pairewise_compare.")

  all <- pwc$significance
  foldchanges <- pwc$foldchanges
  groups <- pwc$groups
  max_stability_fc = 1.4

  results <- list()
  # For the first transition 1
  # +++++++++++++++++++++++
  if(verbose) message("- Identifying genes for step ", pwc$groups[1])
  grp_i <- groups[2]
  current_trans <- apply(all[, pwc$steps[1], drop = FALSE], 1, sum)
  current_fcs <- abs(apply(foldchanges[, pwc$steps[1], drop = FALSE], 1, sum))
  # Stability in the next steps
  isStable_nextsteps <- .is_downSustained_in_nextsteps(pwc, from = 2,
                                                max_stability_fc = max_stability_fc )
  # Stability in the previous steps
  isStable_prevsteps <- .is_downSustained_in_prevsteps(pwc, to = 2,
                                                max_stability_fc = max_stability_fc )
  # Result
  index <- which(current_trans == -1 & isStable_nextsteps == 1 & current_fcs >= log2(fc) & pwc$detection_call == 1)
  if(length(index) == 0) warnings("No genes for step ", grp_i)
  res <- .subset_index(pwc, index, sort_ = "decreasing")
  res$title <- paste0( "Down-at-", grp_i)
  results[[1]] <- res



  # For the remaining steps
  # ++++++++++++++++++++++++++
  for(i in 2:(length(groups)-1)){
    if(verbose) message("- Identifying genes for step ", pwc$groups[i])
    current_trans <- apply(all[, pwc$steps[i], drop = FALSE], 1, sum)
    current_fcs <- abs(apply(foldchanges[, pwc$steps[i], drop = FALSE], 1, sum))
    # Stability in the next steps
    isStable_nextsteps <- .is_downSustained_in_nextsteps(pwc, from = i+1,
                                                  max_stability_fc = max_stability_fc )
    # Stability in the previous steps
    isStable_prevsteps <- .is_stable_in_prevsteps(pwc, to = i,
                                                  max_stability_fc = max_stability_fc )
    # Result
    index <- which(current_trans == -1 &
                     isStable_nextsteps ==1 & isStable_prevsteps ==1 &
                     current_fcs >= log2(fc) & pwc$detection_call == 1)
    if(length(index) == 0) warnings("No genes for step ", groups[i])
    res <- .subset_index(pwc, index, sort_ = "decreasing")
    res$title <- paste0( "Down-at-", groups[i+1])
    results[[i]] <- res
  }
  names(results) <- pwc$groups[2:(length(pwc$groups))]
  structure(results, class = c("list", "stepgenes", "onestepdown"))
}



#' @rdname stepProfiler
#' @export
impulsed_up<- function(pwc, fc = 1.5, verbose = TRUE){
  results <- list()
  for(i in 2:(length(pwc$groups)-1)){
    if(verbose) message("- Identifying genes for step ", pwc$groups[i])
    # Define all impulsed genes
    impulsed <- .is_impulsed_up(pwc, step_index = i, fc = fc)
    # Stability in the next steps
    isStable_nextsteps <- .is_downSustained_in_nextsteps(pwc, from = i+1,max_stability_fc = 1.4 )
    # Stability in the previous steps
    isStable_prevsteps <- .is_upSustained_in_prevsteps(pwc, to = i-1, max_stability_fc = 1.4)
    index <- which(impulsed == 1 & isStable_nextsteps == 1 & isStable_prevsteps == 1 &  pwc$detection_call == 1)
    if(length(index) == 0) warnings("No genes for step ", pwc$groups[i])
    res <- .subset_index(pwc, index, sort_ = "decreasing")
    res$title <- paste0( "Impulsed-up-at-", pwc$groups[i])
    results[[i-1]] <- res
  }
  names(results) <- pwc$groups[2:(length(pwc$groups)-1)]
  structure(results, class = c("list", "stepgenes", "impulsedup"))
}


#' @rdname stepProfiler
#' @export
impulsed_down <- function(pwc, fc = 1.5, verbose = TRUE){
  results <- list()
  for(i in 2:(length(pwc$groups)-1)){
    if(verbose) message("- Identifying genes for step ", pwc$groups[i])
    # Define all impulsed genes
    impulsed <- .is_impulsed_down(pwc, step_index = i, fc = fc)
    # Stability in the next steps
    isStable_nextsteps <- .is_upSustained_in_nextsteps(pwc, from = i+1,max_stability_fc = 1.4 )
    # Stability in the previous steps
    isStable_prevsteps <- .is_downSustained_in_prevsteps(pwc, to = i-1, max_stability_fc = 1.4)
    index <- which(impulsed == 1 & isStable_nextsteps == 1 & isStable_prevsteps == 1 &  pwc$detection_call == 1)
    if(length(index) == 0) warnings("No genes for step ", pwc$groups[i])
    res <- .subset_index(pwc, index, sort_ = "decreasing")
    res$title <- paste0( "Impulsed-down-at-", pwc$groups[i])
    results[[i-1]] <- res
  }
  names(results) <- pwc$groups[2:(length(pwc$groups)-1)]
  structure(results, class = c("list", "stepgenes", "impulseddown"))
}

#' @param x an object of class step_genes to be plotted
#' @param ... other arguments
#' @param transformby a character vector specifying how to transform the expression signal before plotting.
#' Possible values are:
#' \itemize{
#' \item standardize: the rlog data is standardized, before plotting, using the function scale().
#' \item firstep: the log2 fold change expression compared to the first group is displayed
#' \item none: the mean rlog expression is displayed
#' }
#' @param size a numeric vector of length 2 containing the size of lines and points, respectively.
#' @param getRegFunc a function for identifying regulators (transcription factor/epigenetic enzymes) in a list of genes.
#' A possible value is get_human_regulators.
#' @param print_plot a logical value used in the function plot.stepgenes().
#' If TRUE ggplots are printed. If FALSE, a list, which elements are ggplots, is returned instead of being printed.
#' @rdname stepProfiler
#' @method plot stepgenes
#' @export
plot.stepgenes <- function(x, transformby = c("firststep", "basemean", "standardize",  "none"),
                           size = c(0.2, 1), getRegFunc = NULL, print_plot = TRUE, ...){

  if(!inherits(x, "stepgenes"))
    stop("x must be an object of class stepgenes containing 'exprs' and 'group' elements")
  plots <- list()

  transformby <- transformby[1]

  ylab <- "log2 fold change vs the first group"
  for( i in 1:length(x)){
      xx <- x[[i]]
      exprs <- xx$exprs
    if(!is.null(xx$rlog.data)){
      exprs <- xx$exprs
      titl <- paste0(xx$title, " (", nrow(exprs), ")")
      # Transcription regulators
      if(!is.null(getRegFunc)){
        regulators <- getRegFunc(rownames(exprs), keep.all = FALSE)
        titl <- paste0(xx$title, " (", nrow(exprs), ", ", nrow(regulators), ")")
      }
      # Transform
      log_it = TRUE
      if(transformby %in% c("firststep", "firstgroup")){
        if(nrow(exprs) == 1) exprs[1, ] <- exprs[1, ]/xx$refs
        else exprs <- apply(xx$exprs, 2, function(x, ref_vals){x/ref_vals}, xx$refs)
      }
      if(transformby =="basemean") {
        if(nrow(exprs) == 1) exprs[1, ] <- exprs[1, ]/xx$baseMean
        else exprs <- apply(xx$exprs, 2, function(x, ref_vals){x/ref_vals}, xx$baseMean)
        ylab <- "log2 fold change vs base mean"
      }
      else if(transformby == "none") {
        exprs <- xx$rlog.data
        ylab = "Mean rlog expression"
        log_it = FALSE
      }
      else if(transformby == "standardize") {
        exprs <- xx$rlog.data
        exprs <- t(scale(t(exprs)))
        log_it = FALSE
        ylab = "Standardized rlog expression"
      }

      exprs <- summarizeby(exprs, xx$group, mean, melt = TRUE)
      colnames(exprs) <- c("group", "variable", "mean")
      exprs$group <- factor(exprs$group, levels = levels(xx$group))
      if(log_it) exprs$mean <- log2(exprs$mean)
      mapping <-  ggplot2::aes_string(x = "group", y= "mean", group = "variable")
      p <- ggplot2::ggplot(exprs, mapping) +
        ggplot2::geom_line(size = size[1], ...) +  ggplot2::geom_point(size = size[2], ...)+
        ggplot2::labs(y = ylab, x = "Groups", title = titl) +
        ggplot2::theme_classic()
      plots[[i]] <- p
      if(print_plot) print(p)
    }
    else warnings("No one step genes at step ", i)
  }
  invisible(plots)
}

#' @rdname stepProfiler
#' @param object an object of class step_genes to be plotted
#' @method summary stepgenes
#' @export
summary.stepgenes <- function(object, ...)
  {
  if(!inherits(object, "stepgenes"))
    stop("object must be an object of class stepgenes")
  id <- titl <- max_fcs <- min_fcs <- step <- group <- NULL
  nn <- names(object)
  for( i in 1:length(object)){
    xx <- object[[i]]
    if(!is.null(xx$exprs)){
      id <- c(id, rownames(xx$exprs))
      titl <- c(titl, rep(xx$title, nrow(xx$exprs)))
      max_fcs <- c(max_fcs, apply(xx$foldchanges, 1, max))
      min_fcs <- c(min_fcs, apply(xx$foldchanges, 1, min))
      step <- c(step, rep(nn[i], nrow(xx$exprs)))
      group <- c(group, rep(i, nrow(xx$exprs)))
    }
  }

  res <- NULL
  if(!is.null(id)) res <- data.frame(id = id, status = titl, step = step,group = group,
                                     max_fcs = round(max_fcs,2), min_fcs = round(min_fcs,2))
  res
}


#' @rdname stepProfiler
#' @param gn_id_type the type of gene ids. Allowed values are "ENSEMBL"
#' "ENTREZID", "SYMBOL", "UNIGENE", "GENENAME" or "UCSCKG".
#' @param pvalueCutoff Cutoff value of pvalue
#' @param OrgDb annotation database. Default is "org.Hs.eg.db" for human.
#' @export
pathways <- function(object, gn_id_type = "ENSEMBL", pvalueCutoff = 0.05,
                     OrgDb ="org.Hs.eg.db", verbose = TRUE){

  if(!inherits(object, "stepgenes"))
    stop("object must be an object of class stepgenes")
  if(!(gn_id_type %in% c("ENTREZID", "SYMBOL", "UNIGENE", "GENENAME","UCSCKG", "ENSEMBL")))
     stop('Allowed values for gn_id_type are "ENTREZID", "SYMBOL", "UNIGENE", "GENENAME" or "UCSCKG"')
  .load_package(c("clusterProfiler", "ReactomePA", OrgDb, "DOSE"))

  res <- NULL
  if(inherits(object, "onestepup"))
    res <- .stepup_paths(object = object, gn_id_type = gn_id_type,
                         pvalueCutoff = pvalueCutoff, OrgDb = OrgDb, verbose = verbose)
  else if(inherits(object, "onestepdown"))
    res <- .stepdown_paths(object = object, gn_id_type = gn_id_type,
                           pvalueCutoff = pvalueCutoff, OrgDb = OrgDb, verbose = verbose)
  else if(inherits(object, c("impulsedup")))
    res <- .impulsedup_paths(object = object, gn_id_type = gn_id_type,
                             pvalueCutoff = pvalueCutoff, OrgDb = OrgDb,  verbose = verbose)
  else if(inherits(object, c("impulseddown"))){
    res <- .impulsedup_paths(object = object, gn_id_type = gn_id_type,
                             pvalueCutoff = pvalueCutoff, OrgDb = OrgDb,  verbose = verbose)
    if(ncol(res$merged)>0){
      res$merged[,2:ncol(res$merged)] <- apply(res$merged[,2:ncol(res$merged)], 2, function(x){ifelse(x==1, 0, 1)})
    }
  }

  res
}




# step_specific_up(): Increased specifically at a given step but not actively expressed in the others
# (average expression in groups < active_exprs)
step_specific_up <- function(pwc, fc = 1.5, active_exprs = 20, verbose = TRUE){
  max_stability_fc = 1.4
  results <- list()
  for(i in 2:(length(pwc$groups)-1)){
    if(verbose) message("- Identifying genes for step ", pwc$groups[i])
    # Define all impulsed genes
    impulsed <- .is_impulsed_up(pwc, step_index = i, fc = fc)
    # Detection call in previous and next steps
    detection_call <- .get_detection_call(pwc$dds, pwc$factor,
                                          grp1 = pwc$groups[1:(i-1)], grp2 = pwc$groups[(i+1):length(pwc$groups)],
                                          active_exprs = active_exprs, count.norm = pwc$count.norm)
    index <- which(impulsed == 1 & detection_call == 0 & pwc$detection_call == 1)
    if(length(index) == 0) warnings("No genes for step ", pwc$groups[i])
    res <- .subset_index(pwc, index, sort_ = "decreasing")
    res$title <- paste0( "Impulsed-up-at-", pwc$groups[i])
    results[[i-1]] <- res
  }
  names(results) <- pwc$groups[2:(length(pwc$groups)-1)]
  structure(results, class = c("list", "stepgenes"))
}

