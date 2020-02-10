#' @include utilities_diffexprs.R
NULL
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_path
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_blank
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_alpha
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 element_blank
#' @importFrom grid arrow
#' @importFrom grid unit
#' @importFrom utils data
#' @importFrom utils head
#' @importFrom utils read.delim
#' @importFrom utils write.table
#' @importFrom stats IQR
#' @importFrom stats mad
#' @importFrom stats median
#' @importFrom stats p.adjust
#' @importFrom stats quantile
#' @importFrom stats sd

# Create a directory
create_dir <- function(path){
  if(!dir.exists(path)) dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

# Get the list of files in a given directory
list_files<-function(dir=getwd(), ext="txt", recursive = TRUE, full.names = TRUE, ...){
  list.files(path=dir, pattern=(paste(ext, "$", sep="")),
             recursive=recursive, full.names = full.names, ...)
}

# Remove extension in files
remove_extension<-function(file_list){
  for(i in 1:length(file_list)) file_list[i]= unlist(strsplit(file_list[i],"\\."))[1]
  file_list
}

# Install missing packages
# pkgs: vector of packages
# check: if TRUE, check if packages exist
install_pkgs <- function(pkgs, check = TRUE){

  if(check) pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  else pkgs_miss <- pkgs

  if (length(pkgs_miss) > 0) {
    message('Installing missing packages:\n', pkgs_miss)
    source("http://bioconductor.org/biocLite.R")
    all_repos=c('http://cran.r-project.org', BiocInstaller::biocinstallRepos())
    install.packages(pkgs_miss, repos=all_repos)
  }
}

# Save file with overwritting option
# file a .RDATA file
save_file <- function( ..., file, overwrite = FALSE ) {
  if (file.exists(file) & !overwrite ){
   bname<- remove_extension(basename(file))
   bname <- paste0(bname, "_", get_random_string(14), ".RDATA" )
   file2 <- file.path(dirname(file), bname)
   warnings("The file ", file, " already exists.", "The name ", file2, "has been used.")
   file <- file2
  }
  save(..., file = file)
}


# generate random string
get_random_string <- function(length = 14){
  letter<-c("a","b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "x", "y", "z")
  length <- round(length/2)
  txt  <- paste(letter[sample(1:26,length)], collapse ="")
  numb <- paste(sample(1:9,length), collapse = "")
  paste0(txt, numb)
}


.load_package <- function(pkgs){
  for(i in 1:length(pkgs)){
    if (!requireNamespace(pkgs[i], quietly = TRUE)) {
      install_pkgs(pkgs[i], check = FALSE)
      requireNamespace(pkgs[i], quietly = TRUE)
    }
  }
}

.install_pkgs_ifnot <- function(pkgs){
  install_pkgs(pkgs)
}

