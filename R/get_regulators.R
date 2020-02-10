#' Identification of a transcription regulators in a list of genes
#' @description
#' Identification of DNA-binding transcription factors and/or epigenetic enzymes in a list of genes.
#' Only data for human are available (see ?human.tf and ?human.ee).
#' @param ids a vector containing gene Ensembl Ids ("NSG00000130816") or gene names ("DNMT1").
#' if ids is a matrix or a data.frame, then the row names are used as ids.
#' @param id_type the type of the ids. Possible values are "ensembl" or "name".
#' @param keep.all a logical value. If TRUE, all the genes are kept in the final result. If FALSE, only
#' transcription factors/epigenetic enzymes are returned
#' @return
#' Returns an object of class data.frame containing the following columns:
#' \itemize{
#' \item id: gene ids
#' \item TF: with a value of 1 (the genes is transcription factor) or 0 (the gene is not a transcription factor).
#' \item EE_class: specifies the epigenetic enzyme class of the genes.
#' }
#' @export
get_human_regulators <- function(ids, id_type = "ensembl", keep.all = TRUE){

  d <- ids
  if(!id_type %in% c("ensembl", "name"))
    stop("Allowed values for id_type are 'ensembl' or 'name'.")
  if(inherits(ids, c("matrix", "data.frame")))
    ids <- data.frame(id = rownames(ids))
  else ids <- data.frame(id = unique(ids))

  by.y <- "ensembl"
  if(id_type == "ensembl"){
    found <- intersect(ids$id, rownames(biomart.hs.genes))
    if(length(found) == 0) stop("Are you sure that your ids are Ensembl Ids? ",
                                "Can't find any of your ids in human Ensembl data base (>= GRCh38). ",
                                "Try id_type = 'name', if your ids are gene names.")
    by.y <- "ensembl"
  }
  else if(id_type == "name"){
    found <- intersect(ids$id, as.vector(biomart.hs.genes$name))
    if(length(found) == 0) stop("Are you sure that your ids are gene names? ",
                                "Can't find any of your ids in human gene names (Ensembl data base >= GRCh38). ",
                                "Try id_type = 'ensembl', if your ids are Ensembl ids.")
    by.y <- "name"
  }

  # Transcription factors
  tf <- human.tf
  tf$TF <-  tf$name
  res <- merge(ids, tf[, c(by.y, "TF")], by.x = "id", by.y = by.y, sort =  FALSE, all.x = TRUE)
  # res$TF <- ifelse(is.na(res$TF), 0, 1)

  # Epigenetic Enzyme
  ee <- human.ee
  colnames(ee)[which(colnames(ee) == "class")] <- "EE_class"
  res <- merge(res, ee[, c(by.y, "EE_class")], by.x = "id", by.y = "ensembl", sort = FALSE, all.x = TRUE)
  rownames(res) <- as.vector(res$id)

  if(inherits(d, c("matrix", "data.frame"))){
    diff_column <- setdiff(colnames(res), colnames(d))
    res <- cbind.data.frame(d, res[rownames(d), diff_column, drop = FALSE])
  }
  if(!keep.all)  res <- subset(res, (!is.na(res$TF) | !is.na(res$EE_class)))

  res
}
