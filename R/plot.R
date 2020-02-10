#' @include utilities.R summarize.R
NULL
#' Graphs
#' @description
#' Create graphs
#' @param x a numeric vector
#' @param title,xlab,ylab labels for the graph
#' @param fill a fill color
#' @param color for line or point
#' @param sort.value a string specifying whether x should be sorted or not.
#' Allowed values are "none" (no sorting), "asc" (for ascending) or "desc" (for descending)
#' @param top a numeric value specifing the top elements to be shown
#' @param ... additional arguments
#' @examples
#' \donttest{
#' data(mtcars)
#' ggbarplot(mtcars[, "mpg"], top = 10)
#' }
#' @name plot
#' @rdname plot
#' @export
ggbarplot <- function(x,  title ="", xlab ="", ylab="",
                       fill="steelblue", color = "steelblue",
                       sort.value = c("none", "desc", "asc"), top = Inf, ... ){
  # top elements
  if(top!=Inf & top < length(x))
    x <- sort(x, decreasing=TRUE)[1:top]
  # sorting
  if(sort.value[1]=="desc") x <- sort(x, decreasing = TRUE)
  else if(sort.value[1]=="asc") x <- sort(x, decreasing = FALSE)
  # bar names
  if(is.null(names(x))) names(x) <- 1:length(x)

  #data frame for ggplot2
  d <- cbind.data.frame(name = factor(names(x), levels = names(x)), val = x)
  # plot
  p <- ggplot(d, aes_string("name", "val")) +
    geom_bar(stat="identity", fill=fill, color = color)+
    labs(title = title,  x =xlab, y = ylab)+
    theme(axis.text.x = element_text(angle=45), axis.title.x = element_blank())

  return(p)
}

#' @param data a data frame
#' @param xName,yName variables to be used for x and y, respectively
#' @rdname plot
#' @export
ggboxplot <- function(data, xName, yName,
                      fill="steelblue", color = "black")
{
  data <- as.data.frame(data)
  ggplot(data, aes_string(xName, yName))+
    ggplot2::geom_boxplot(fill = fill, color = color)+
    theme(axis.text.x = element_text(angle=45), axis.title.x = element_blank())
}


#' @param exprs gene expression data (rows are genes and columns are samples)
#' @param groups a vector containing the group of each samples
#' @param labels logical value. If TRUE point labels are shown.
#' @param log_it a logical value. If TRUE the data are log2 transformed.
#' @rdname plot
#' @export
gglineplot_by <- function(exprs, groups, log_it = TRUE, title = "", labels = FALSE, xlab = "Groups",
                            ylab = "log2 fold change compared to the control", ...){

  if(!inherits(data, "list"))
    stop("data must be an object of class list containing 'exprs' and 'group' elements")
  if(!inherits(groups, "factor")) groups <- as.factor(groups)
  exprs <- summarizeby(exprs, groups, mean, melt = TRUE)
  colnames(exprs) <- c("group", "variable", "mean")
  exprs$group <- factor(exprs$group, levels = levels(groups))
  if(log_it) exprs$mean <- log2(exprs$mean)
  mapping <-  ggplot2::aes_string(x = "group", y= "mean", group = "variable")
  if(labels) mapping <-  ggplot2::aes_string(x = "group", y= "mean", group = "variable", color = "variable")
  p <- ggplot2::ggplot(exprs, mapping) +
    ggplot2::geom_line(...) +  ggplot2::geom_point(...)+
    ggplot2::labs(y = ylab, x = "Groups", title = title) +
    ggplot2::theme_classic()
  p
}
