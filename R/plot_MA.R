#' MA-plot from means and log fold changes
#'
#' @description
#' Make MA-plot which is a scatter plot of log2 fold changes (on the y-axis)
#' versus the mean expression signal (on the x-axis).
#' @param x an object of class DESeqResults, get_diff, DE_Results, matrix or data frame containing the columns
#' baseMean, log2FoldChange, and padj. Rows are genes.
#' \itemize{
#' \item baseMean: the mean expression of genes in the two groups.
#' \item log2FoldChange: the log2 fold changes of group 2 compared to group 1
#' \item padj: the adjusted p-value of the used statiscal test.
#' }
#' @param fdr Accepted false discovery rate for considering genes as differentially expressed.
#' @param fc the fold change threshold. Only genes with a fold change >= fc and padj <= fdr are
#' considered as significantly differentially expressed.
#' @param detection_call a numeric vector with length = nrow(x),
#' specifying if the genes is expressed or not. Default is NULL. If detection_call column is available in x,
#' it will be used.
#' @param ylim The limits for the y-axis.
#' If missing, an attempt is made to choose a sensible value.
#' Dots exceeding the limits will be displayed as triangles at the limits, pointing outwards.
#' @param colNonSig,colUp,colDown colors to be used for non-significant,
#' up-regulated and down-regulated data points, respectively.
#' @param colLine colour to use for the horizontal (y=0) line.
#' @param log which axis/axes should be logarithmic; will be passed to plot.
#' @param cex point size
#' @param xlab,ylab x and y axis labels
#' @param addlegend a logical value. If TRUE a legend is added to the plot.
#' @param highlight highlight top genes. Default value is 10 for top 10 genes.
#' @param ... other parameters to be passed to the function plot
#' @return returns a base plot
#' @export
plot_MA <- function (x, fdr = 0.05, fc = 1.5, detection_call = NULL,
                     ylim = NULL, colNonSig = "darkgray",
                     colUp = "#BC002A",  colDown ="#16614E",
                     colLine = "black", log = "x", cex = 0.45,
                     xlab = "Mean expression", ylab = "Log2 fold change",
                     addlegend = TRUE, highlight = 10, ...)
{

 if(!base::inherits(x, c("matrix", "data.frame", "DataFrame", "DE_Results", "DESeqResults")))
   stop("x must be an object of class matrix, data.frame, DataFrame, DE_Results or DESeqResults")
  if(!is.null(detection_call)){
    if(nrow(x)!=length(detection_call))
      stop("detection_call must be a numeric vector of length = nrow(x)")
  }
  else if("detection_call" %in% colnames(x)){
    detection_call <- as.vector(x$detection_call)
  }
  else detection_call = rep(1, nrow(x))

  ss <- base::setdiff(c("baseMean", "log2FoldChange", "padj"), colnames(x))

  if(length(ss)>0) stop("The colnames of x must contain: ",
                        paste(ss, collapse = ", "))

  x <- data.frame(name = rownames(x), mean = x$baseMean, lfc = x$log2FoldChange,
                       sig = ifelse(is.na(x$padj), FALSE,
                                    x$padj <= fdr & abs(x$log2FoldChange) >= log2(fc) &
                                      detection_call ==1 ),
                  padj = x$padj)

  colnames(x) <- c("name", "mean", "lfc", "sig", "padj")

  py <- x$lfc
  if (is.null(ylim)) ylim <- 1.1 * range(py[is.finite(py)])
  pch <- ifelse(py < ylim[1], 6, ifelse(py > ylim[2],  2, 16))

  col <- rep(colNonSig, nrow(x))
  col[which(x$sig == TRUE & x$lfc < 0)] = colDown
  col[which(x$sig == TRUE & x$lfc > 0)] = colUp

  plot(x$mean, pmax(ylim[1], pmin(ylim[2], py)), log = log,
       pch = pch, cex = cex, col = col, xlab = xlab, ylab = ylab,
       ylim = ylim, frame = FALSE,
       font.main=2, font.lab=2, font.axis = 2,  ...)
  abline(h = 0, lwd = 1.5, col = colLine, lty = 2)

  if(addlegend){
    up <- which(x$sig == TRUE & x$lfc > 0)
    up_legend <- paste0("Up: ", length(up))
    down <- which(x$sig == TRUE & x$lfc < 0)
    down_legend <- paste0("Down: ", length(down))
    legend("topleft", legend=c(up_legend,  down_legend),
           col=c(colUp, colDown), pch= c(pch, pch), box.lty=0, cex=0.8)
  }

  # top genes
  if(highlight > 0){
    x <- x[order(abs(x$lfc), decreasing = TRUE), ]
    x <- head(subset(x, x$sig==1), highlight)
    text(x$mean, pmax(ylim[1], pmin(ylim[2], x$lfc)), labels = x$name, pos = 4, cex = cex)
  }
}

#' @rdname plot_MA
#' @method plot DESeqResults
#' @export
plot.DESeqResults <- function(x, ...){
  plot_MA(x, ...)
}

#' @rdname plot_MA
#' @method plot DE_Results
#' @export
plot.DE_Results <- function(x, ...){
  plot_MA(x, ...)
}



