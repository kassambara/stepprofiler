
# Subset result for an index of genes
# +++++++++++++++++++++++++++++++++++++++++
## pwc an object of class pairewise compare
## index a vector of gene index
# sort_: order by mean_fc, possible values for order are "decreasing", or "increasing"
.subset_index <- function(pwc, index, sort_ = "decreasing"){

  fcs <- padj <- exprs <- rlog.data  <- refs <- baseMean <- NULL
  samples <- subset(pwc$samples, pwc$samples[, pwc$factor] %in% pwc$groups)
  if(length(index) > 0){
    fcs <- pwc$foldchanges[index, pwc$steps, drop = FALSE]
    padj <- pwc$padj[index, pwc$steps , drop = FALSE]
    exprs <- as.data.frame(pwc$count.norm[index, rownames(samples), drop = FALSE ])
    rlog.data <- as.data.frame(pwc$rlog.data[index, rownames(samples), drop = FALSE ])
    refs <- pwc$refs[rownames(exprs)]
    baseMean <- pwc$baseMean[rownames(exprs)]

    # Order
    mean_fc <- pwc$steps_mean_fc[rownames(exprs)]
    if(sort_ == "decreasing") ss <- order(mean_fc, decreasing = TRUE)
    else if(sort_ == "increasing") ss <- order(mean_fc, increasing = TRUE)
    fcs <- fcs[ss, , drop = FALSE]
    padj <- padj[ss, , drop = FALSE]
    exprs <- exprs[ss, , drop = FALSE]
    rlog.data <- rlog.data[ss, , drop = FALSE]
    refs <- refs[ss]
    baseMean <- baseMean[ss]
  }

  res <- list(exprs = exprs, rlog.data = rlog.data, foldchanges = fcs, padj = padj,
              group = factor(samples[, pwc$factor], levels = pwc$groups),
              refs = refs, baseMean = baseMean)
}

# UP at each transition step
#+++++++++++++++++++++++++
# pwc pairewise_compare object
# All step up
.all_step_up <- function(pwc){

  if(!inherits(pwc, "pairewise_compare"))
    stop("pwc must be an object of class pairewise_compare.")

  freq <- apply(pwc$significance[, pwc$steps, drop =FALSE], 1, sum)
  index <- which(freq == pwc$nstep)
  res <- .subset_index(pwc, index)
  res$title <- paste0("Up at ", pwc$groups[1])
  res
}


# All step down
# +++++++++++++++++++++++++++
.all_step_down <- function(pwc){

  if(!inherits(pwc, "pairewise_compare"))
    stop("pwc must be an object of class pairewise_compare.")

  freq <- apply(pwc$significance[, pwc$steps, drop =FALSE], 1, sum)
  index <- which(freq == -pwc$nstep)
  res <- .subset_index(pwc, index)
  res$title <- paste0("Down at ", pwc$groups[1])
  res
}

# Possible comparison between group 1 and the others
#++++++++++++++++++++++++++++++++++++++++
.possible_groups <- function(groups){tolower(paste0(groups[1], "_vs_", groups[2:length(groups)]))}


# Define genes that are stable or upregulated in the next step
# +++++++++++++++++++++++++++++++++++++++++++
.is_upSustained_in_nextsteps <- function(pwc, from, max_stability_fc = 1.4 ){
  if(max_stability_fc  < 1) stop("The max_stability_fcmust be >= 1")
  fc <- log2(max_stability_fc)
  res <- NULL
  # if i = ngroups, there is no further steps for evaluation of stable signal => all the genes are considered as stable
  if(from == length(pwc$groups)) res <- rep(1, nrow(pwc$count.norm))
  else{
    transition <- pwc$significance[, pwc$steps[from:length(pwc$steps)], drop = FALSE]
    transition <- apply(transition, 1,
                 function(x){
                   x <- as.numeric(.is_inRange(x, c(0, 1)))
                   sum(x) == length(x)})
    fcs <- pwc$foldchanges[, pwc$steps[from:length(pwc$steps)], drop = FALSE]
    fcs <- apply(fcs, 1,
                 function(x){
                   x <- as.numeric(.is_inRange(x, c(-fc, Inf)))
                   sum(x) == length(x)})
    res  <- (transition == 1 & fcs == 1)
  }
  res
}

# Define genes that are stable or upregulated in the prev step
# +++++++++++++++++++++++++++++++++++++++++++
.is_upSustained_in_prevsteps <- function(pwc, to, max_stability_fc = 1.4 ){
  if(max_stability_fc  < 1) stop("The max_stability_fcmust be >= 1")
  fc <- log2(max_stability_fc)
  res <- NULL
  # if to = 1, there is no previous steps for evaluation of stable signal => all the genes are considered as stable
  if(to == 1) res <- rep(1, nrow(pwc$count.norm))
  else{
    transition <- pwc$significance[, pwc$steps[1:(to-1)], drop = FALSE]
    transition <- apply(transition, 1,
                        function(x){
                          x <- as.numeric(.is_inRange(x, c(0, 1)))
                          sum(x) == length(x)})
    fcs <- pwc$foldchanges[, pwc$steps[1:(to-1)], drop = FALSE]
    fcs <- apply(fcs, 1,
                 function(x){
                   x <- as.numeric(.is_inRange(x, c(-fc, Inf)))
                   sum(x) == length(x)})
    res  <- (transition == 1 & fcs == 1)
  }
  res
}


# Define genes that are stable or downregulated in the next step
# +++++++++++++++++++++++++++++++++++++++++++
.is_downSustained_in_nextsteps <- function(pwc, from, max_stability_fc = 1.4 ){
  if(max_stability_fc  < 1) stop("The max_stability_fcmust be >= 1")
  fc <- log2(max_stability_fc)
  res <- NULL
  # if i = ngroups, there is no further steps for evaluation of stable signal => all the genes are considered as stable
  if(from == length(pwc$groups)) res <- rep(1, nrow(pwc$count.norm))
  else{
    transition <- pwc$significance[, pwc$steps[from:length(pwc$steps)], drop = FALSE]
    transition <- apply(transition, 1,
                        function(x){
                          x <- as.numeric(.is_inRange(x, c(-1, 0)))
                          sum(x) == length(x)})
    fcs <- pwc$foldchanges[, pwc$steps[from:length(pwc$steps)], drop = FALSE]
    fcs <- apply(fcs, 1,
                 function(x){
                   x <- as.numeric(.is_inRange(x, c(-Inf, fc)))
                   sum(x) == length(x)})
    res  <- (transition == 1 & fcs == 1)
  }
  res
}

# Define genes that are stable or downregulated in the prev step
# +++++++++++++++++++++++++++++++++++++++++++
.is_downSustained_in_prevsteps <- function(pwc, to, max_stability_fc = 1.4 ){
  if(max_stability_fc  < 1) stop("The max_stability_fcmust be >= 1")
  fc <- log2(max_stability_fc)
  res <- NULL
  # if i = ngroups, there is no further steps for evaluation of stable signal => all the genes are considered as stable
  if(to == length(pwc$groups)) res <- rep(1, nrow(pwc$count.norm))
  else{
    transition <- pwc$significance[, pwc$steps[1:(to-1)], drop = FALSE]
    transition <- apply(transition, 1,
                        function(x){
                          x <- as.numeric(.is_inRange(x, c(-1, 0)))
                          sum(x) == length(x)})
    fcs <- pwc$foldchanges[, pwc$steps[1:(to-1)], drop = FALSE]
    fcs <- apply(fcs, 1,
                 function(x){
                   x <- as.numeric(.is_inRange(x, c(-Inf, fc)))
                   sum(x) == length(x)})
    res  <- (transition == 1 & fcs == 1)
  }
  res
}




# Define if the expression of the genes is stable from step i to step i+1
#++++++++++++++++++++++++++++
## pwc: an object of class pairewise compare
## from: group index from which the next steps are defined
## max_stability_fc define the maximum fold change between step i and step i +1
#  a stable gene is gene which are not differentially expressed between
# step i to step i+1 and with a maximum absolute fold change of max_stability_fc
.is_stable_in_nextsteps<- function(pwc, from, max_stability_fc = 1.4 ){
  if(max_stability_fc  < 1) stop("The max_stability_fcmust be >= 1")
  fc <- log2(max_stability_fc)
  res <- NULL
  # if i = ngroups, there is no further steps for evaluation of stable signal => all the genes are considered as stable
  if(from == length(pwc$groups)) res <- rep(1, nrow(pwc$count.norm))
  else{
  possible_comparison <- .get_possible_parwise_comparison(pwc$groups[from:length(pwc$groups)])

  transition <- pwc$significance[, possible_comparison, drop = FALSE]
  transition <- apply(abs(transition), 1, sum)
  fcs <- pwc$foldchanges[, .possible_groups(pwc$groups[from:length(pwc$groups)]), drop = FALSE]
  fcs <- apply(fcs, 1,
               function(x){
                 x <- as.numeric(.is_inRange(x, c(-fc, fc)))
                 sum(x) == length(x)})
  res  <- transition == 0 & fcs == 1
  }
  res
}

# Check if a gene is stable in previous steps
#+++++++++++++++++++++++++++
## pwc pairewise compare
## to: group index to which previous steps are defined
#  max_stability_fc = 1.4  not used
.is_stable_in_prevsteps<- function(pwc, to, max_stability_fc = 1.4 ){
  if(max_stability_fc  < 1) stop("The max_stability_fc must be >= 1")
  fc <- log2(max_stability_fc)

  res <- NULL
  # if to = 1, there is no previous steps for evaluation of stable signal => all the genes are considered as stable
  if(to == 1) res <- rep(1, nrow(pwc$count.norm))
  else{
    possible_comparison <- .get_possible_parwise_comparison(pwc$groups[1:to])

    transition <- pwc$significance[, possible_comparison, drop = FALSE]
    transition <- apply(abs(transition), 1, sum)

    fcs <- pwc$foldchanges[, .possible_groups(pwc$groups[1:to]), drop = FALSE]
    fcs <- apply(fcs, 1,
                 function(x){
                   x <- as.numeric(.is_inRange(x, c(-fc, fc)))
                   sum(x) == length(x)})
   # res <- transition == 0 & fcs == 1
    res <- transition == 0
  }
  res
}


# Impulsed-up genes
# ++++++++++++++++++++++++++++++
# Impulsed-up genes are genes that are upregulated from the transition i-1 to i and
# then downregulated from i to i+1 transition
# return a vector of logical values with TRUE and FALSE values
# fc minimum fold change of up/downregulation
.is_impulsed_up <- function(pwc, step_index, fc = 1.5){
  fc <- log2(fc)
  # Define previous and next groups
  grp_i <- pwc$groups[step_index]
  prev_grp <- pwc$groups[step_index-1]
  next_grp <- pwc$groups[step_index+1]
  # Define transitions
  prev_transition <- tolower(paste0(prev_grp, "_vs_", grp_i))
  next_transition <- tolower(paste0(grp_i, "_vs_", next_grp))
  # Fold change expression during next and previous transitions
  prev_trans_fc <- pwc$foldchanges[, prev_transition, drop = TRUE]
  next_trans_fc <- pwc$foldchanges[, next_transition, drop = TRUE]
  # Check if foldchanges are in ranges
  prev_trans_fc_isOK <- .is_inRange(prev_trans_fc, c(fc, Inf)) # upregulation?
  next_trans_fc_isOK <- .is_inRange(next_trans_fc, c(-Inf, -fc)) # downregulation?
  # Final results
  return(next_trans_fc_isOK == 1 & prev_trans_fc_isOK == 1)
}

# Impulsed-down genes are genes that are downregulated from the transition i-1 to i and
# then upregulated from i to i+1 transition
# return a vector of logical values with TRUE and FALSE values
# fc minimum fold change of up/downregulation
.is_impulsed_down <- function(pwc, step_index, fc = 1.5){
  fc <- log2(fc)
  # Define previous and next groups
  grp_i <- pwc$groups[step_index]
  prev_grp <- pwc$groups[step_index-1]
  next_grp <- pwc$groups[step_index+1]
  # Define transitions
  prev_transition <- tolower(paste0(prev_grp, "_vs_", grp_i))
  next_transition <- tolower(paste0(grp_i, "_vs_", next_grp))
  # Fold change expression during next and previous transitions
  prev_trans_fc <- pwc$foldchanges[, prev_transition, drop = TRUE]
  next_trans_fc <- pwc$foldchanges[, next_transition, drop = TRUE]
  # Check if foldchanges are in ranges
  prev_trans_fc_isOK <- .is_inRange(prev_trans_fc, c(-Inf, -fc)) # downregulation?
  next_trans_fc_isOK <- .is_inRange(next_trans_fc, c(fc, Inf)) # upregulation?
  # Final results
  return(next_trans_fc_isOK == 1 & prev_trans_fc_isOK == 1)
}


# genes that are significantly upregulated during all previous transition steps
# and then significantly downregulated during all next transition steps
# minimum fold change expression during transitions are 1.5
.is_up_down <- function(pwc, step_index){
  fc <- log2(1)
  i <- step_index
  # Define transitions
  prev_transitions <- pwc$steps[1:(i-1)]
  nprev <- length(prev_transitions)
  next_transitions <- pwc$steps[i:length(pwc$steps)]
  nnext <- length(next_transitions)
  # Significance of differential expression in all transistions
  prev_trans_sig <- apply(pwc$significance[, prev_transitions, drop = FALSE],
                          1, sum)
  next_trans_sig <- apply(pwc$significance[, next_transitions, drop = FALSE],
                          1, sum)
  # Fold change expression during next and previous transitions
  prev_trans_fcs <- pwc$foldchanges[, prev_transitions, drop = TRUE]
  next_trans_fcs <- pwc$foldchanges[, next_transitions, drop = TRUE]
  # Check if foldchanges are in ranges
  prev_trans_fcs_isOK <- .is_inRange(prev_trans_fcs, c(fc, Inf)) # upregulation?
  next_trans_fcs_isOK <- .is_inRange(next_trans_fcs, c(-Inf, -fc)) # downregulation?
  # Final result
  return(prev_trans_sig == nprev & next_trans_sig == -nnext &  next_trans_fcs_isOK == 1 & prev_trans_fcs_isOK == 1)
}

# genes that are significantly downregulated during all previous transition steps
# and then significantly upregulated during all next transition steps
# minimum fold change expression during transitions are 1.5
.is_down_up <- function(pwc, step_index){
  fc <- log2(1)
  i <- step_index
  # Define transitions
  prev_transitions <- pwc$steps[1:(i-1)]
  nprev <- length(prev_transitions)
  next_transitions <- pwc$steps[i:length(pwc$steps)]
  nnext <- length(next_transitions)
  # Significance of differential expression in all transistions
  prev_trans_sig <- apply(pwc$significance[, prev_transitions, drop = FALSE],
                          1, sum)
  next_trans_sig <- apply(pwc$significance[, next_transitions, drop = FALSE],
                          1, sum)
  # Fold change expression during next and previous transitions
  prev_trans_fcs <- pwc$foldchanges[, prev_transitions, drop = TRUE]
  next_trans_fcs <- pwc$foldchanges[, next_transitions, drop = TRUE]

  # Check if foldchanges are in ranges
  prev_trans_fcs_isOK <- .is_inRange(prev_trans_fcs,c(-Inf, -fc))  # upregulation?
  next_trans_fcs_isOK <- .is_inRange(next_trans_fcs, c(fc, Inf)) # downregulation?

  # Final result
  return(prev_trans_sig == -nprev & next_trans_sig == nnext &  next_trans_fcs_isOK == 1 & prev_trans_fcs_isOK == 1)
}



# grps a vector of groups
.get_possible_parwise_comparison <- function(grps){
  res <- NULL
  for(i in 1:(length(grps)-1)){
    grp_i <- grps[i]
    for(j in (i+1):length(grps)){
      grp_j <- grps[j]
      res<-c(res, tolower(paste0(grp_i, "_vs_", grp_j)))
    }
  }
  res
}

# Check if a value is in a range
# val can be numeric vector, matrix or data frame
.is_inRange <- function(val, .range){
  if(inherits(val, c("matrix", "data.frame"))){

    apply(val, 1,
          function(x){
            x <- as.numeric(.is_inRange(x, .range ))
            sum(x) == length(x) })
  }
  else min(.range) <= val & val <= max(.range)
}

# Pathways
# object: an object of class onestepup
# gn_id_type: the type of gene ids. Allowed values are "ENSEMBL"
# "ENTREZID", "SYMBOL", "UNIGENE", "GENENAME" or "UCSCKG".
# for each step i, the gene set used is the genes in i and all i--
.stepup_paths <- function(object, gn_id_type = "ENSEMBL", pvalueCutoff = 0.05, OrgDb = "org.Hs.eg.db", verbose = TRUE){

  if(!inherits(object, "onestepup"))
    stop("object must be an object of class onestepup")
  if(!(gn_id_type %in% c("ENTREZID", "SYMBOL", "UNIGENE", "GENENAME","UCSCKG", "ENSEMBL")))
    stop('Allowed values for gn_id_type are "ENTREZID", "SYMBOL", "UNIGENE", "GENENAME" or "UCSCKG"')
  .install_pkgs_ifnot(c("clusterProfiler", "ReactomePA", OrgDb, "DOSE"))
  ss <- summary(object)
  grps <- unique(ss$group)

  # Step names
  step_names <- NULL
  g <- as.vector(ss$group)
  step <- as.vector(ss$step)
  for(grp in grps) step_names <- c(step_names, unique(step[g == grp]))
  paths = list()

  # Pathways
  for( i in 1: length(grps)) {
    if(verbose) message("Identifying pathways for step ", i)
    ids <- rownames(subset(ss, ss$group <= grps[i]))
    entrez <- clusterProfiler::bitr(ids, fromType = gn_id_type,
                                    toType = "ENTREZID",  OrgDb=OrgDb)$ENTREZID
    paths[[i]] <- ReactomePA::enrichPathway(gene= entrez, pvalueCutoff = pvalueCutoff, readable=T)
  }

  # Merge paths:
  path.annot <- NULL
  for(i in 1:length(paths)){
    path.annot <- rbind(path.annot, DOSE::summary(paths[[i]])[, c("ID", "Description", "pvalue")])
  }
  path.annot <- path.annot[order(path.annot$pvalue),]
  path.annot <- path.annot[!duplicated(path.annot$ID), ]

  all <- DOSE::summary(paths[[1]])[, "ID", drop = FALSE]
  all[, step_names[1]] <- rep(1, nrow(all))
  for(i in 2:length(paths)){
    p <- DOSE::summary(paths[[i]])[, "ID",  drop = FALSE]
    p[, step_names[i]] <- rep(1, nrow(p))
    colnames(p)[ncol(p)] <- step_names[i]
    all <- merge(all, p, by.x = "ID", by.y = "ID", all.y = TRUE)
  }
  all <- merge(path.annot, all, by.x = "ID", by.y = "ID", all.y = TRUE )
  all <- all[!duplicated(all$Description), ]
  rownames(all) <- all$Description
  all <- all[, c(-1, -2), drop = FALSE] # Remove ID and descriptioncolumns
  all[, 2:ncol(all)] <- apply(all[, 2:ncol(all)], 2, function(x) ifelse(is.na(x), 0, 1))
  rsum <- apply(all[, 2:ncol(all)], 1, sum)
  all <- all[order(-rsum, all$pvalue), ]

  names(paths) <- step_names
  list(paths = paths, merged = all)
}



.stepdown_paths <- function(object, gn_id_type = "ENSEMBL", pvalueCutoff = 0.05, OrgDb = "org.Hs.eg.db", verbose = TRUE){

  if(!inherits(object, "onestepdown"))
    stop("object must be an object of class onestepdown")
  if(!(gn_id_type %in% c("ENTREZID", "SYMBOL", "UNIGENE", "GENENAME","UCSCKG", "ENSEMBL")))
    stop('Allowed values for gn_id_type are "ENTREZID", "SYMBOL", "UNIGENE", "GENENAME" or "UCSCKG"')
  .install_pkgs_ifnot(c("clusterProfiler", "ReactomePA", OrgDb, "DOSE"))
  ss <- summary(object)
  grps <- unique(ss$group)

  # Step names
  step_names <- NULL
  g <- as.vector(ss$group)
  step <- as.vector(ss$step)
  for(grp in grps) step_names <- c(step_names, unique(step[g == grp]))
  paths = list()

  # Pathways
  for( i in 1: length(grps)) {
    if(verbose) message("Identifying pathways for step ", i)
    ids <- rownames(subset(ss, ss$group %in% c(grps[i]:length(grps))))
    entrez <- clusterProfiler::bitr(ids, fromType = gn_id_type,
                                    toType = "ENTREZID",  OrgDb=OrgDb)$ENTREZID
    paths[[i]] <- ReactomePA::enrichPathway(gene= entrez, pvalueCutoff = pvalueCutoff, readable=T)
  }

  # Merge paths:
  path.annot <- NULL
  for(i in 1:length(paths)){
    path.annot <- rbind(path.annot, DOSE::summary(paths[[i]])[, c("ID", "Description", "pvalue")])
  }
  path.annot <- path.annot[order(path.annot$pvalue),]
  path.annot <- path.annot[!duplicated(path.annot$ID), ]

  all <- DOSE::summary(paths[[1]])[, "ID", drop = FALSE]
  all[, step_names[1]] <- rep(1, nrow(all))
  for(i in 2:length(paths)){
    p <- DOSE::summary(paths[[i]])[, "ID",  drop = FALSE]
    p[, step_names[i]] <- rep(1, nrow(p))
    colnames(p)[ncol(p)] <- step_names[i]
    all <- merge(all, p, by.x = "ID", by.y = "ID", all.x = TRUE)
  }
  all <- merge(path.annot, all, by.x = "ID", by.y = "ID", all.y = TRUE )
  all <- all[!duplicated(all$Description), ]
  rownames(all) <- all$Description

  all <- all[, c(-1, -2), drop = FALSE] # Remove ID and description columns
  all[, 2:ncol(all)] <- apply(all[, 2:ncol(all)], 2, function(x) ifelse(is.na(x), 0, 1))
  rsum <- apply(all[, 2:ncol(all)], 1, sum)
  all <- all[order(rsum, all$pvalue), ]

  names(paths) <- step_names
  list(paths = paths, merged = all)
}


.impulsedup_paths <- function(object, gn_id_type = "ENSEMBL",
                              pvalueCutoff = 0.05, OrgDb = "org.Hs.eg.db", verbose = TRUE){

  if(!inherits(object, c("impulsedup", "impulseddown")))
    stop("object must be an object of class impulsedup or impulseddown")
  if(!(gn_id_type %in% c("ENTREZID", "SYMBOL", "UNIGENE", "GENENAME","UCSCKG", "ENSEMBL")))
    stop('Allowed values for gn_id_type are "ENTREZID", "SYMBOL", "UNIGENE", "GENENAME" or "UCSCKG"')
  .install_pkgs_ifnot(c("clusterProfiler", "ReactomePA", OrgDb, "DOSE"))
  ss <- summary(object)
  grps <- unique(ss$group)

  # Step names
  step_names <- NULL
  g <- as.vector(ss$group)
  step <- as.vector(ss$step)
  for(grp in grps) step_names <- c(step_names, unique(step[g == grp]))
  paths = list()

  # Pathways
  for( i in 1:length(grps)){
    if(verbose) message("Identifying pathways for step ", i)
    ids <- rownames(subset(ss, ss$group == grps[i]))
    entrez <- clusterProfiler::bitr(ids, fromType = gn_id_type,
                                    toType = "ENTREZID",  OrgDb=OrgDb)$ENTREZID
    paths[[i]] <- ReactomePA::enrichPathway(gene= entrez, pvalueCutoff = pvalueCutoff, readable=T)

  }

  # Merge paths:
  path.annot <- NULL
  for(i in 1:length(paths)){
    path.annot <- rbind(path.annot,
                        DOSE::summary(paths[[i]])[, c("ID", "Description", "pvalue")]
                        )
  }
  path.annot <- path.annot[order(path.annot$pvalue),]
  path.annot <- path.annot[!duplicated(path.annot$ID), ]

  all <- DOSE::summary(paths[[1]])[, "ID", drop = FALSE]
  all[, step_names[1]] <- rep(1, nrow(all))
  for(i in 2:length(paths)){
    p <- DOSE::summary(paths[[i]])[, "ID",  drop = FALSE]
    p[, step_names[i]] <- rep(1, nrow(p))
    colnames(p)[ncol(p)] <- step_names[i]
    all <- merge(all, p, by.x = "ID", by.y = "ID", all.y = TRUE, all.x = TRUE)
  }
  all <- merge(path.annot, all, by.x = "ID", by.y = "ID", all.y = TRUE )
  all <- all[!duplicated(all$Description), ]
  rownames(all) <- all$Description

  all <- all[, c(-1, -2), drop = FALSE] # Remove ID and descriptioncolumns
  all[, 2:ncol(all)] <- apply(all[, 2:ncol(all)], 2, function(x) ifelse(is.na(x), 0, 1))

  # Ordering
  weight <- 1:(ncol(all)-1)
  weight <- t(apply(all[, 2:ncol(all)], 1, "*" , weight))
  weight <- apply(weight, 1, sum)
  all <- all[order(weight, all$pvalue), ]

  names(paths) <- step_names
  list(paths = paths, merged = all)
}







