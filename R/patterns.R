#' Available Temporal Patterns
#' @importFrom dplyr mutate
#' @importFrom dplyr %>%
#' @importFrom ggplot2 ggplot
#' @description Shows available gene expression patterns
#' @return a ggplot
#' @examples
#' patterns()
#' @export
patterns <- function(){

  # Data: one step
  x <- seq(0,5,0.5)
  data.one.step.down <- data.frame(x = x) %>%
    mutate(y = 2^(-x*2), type = "One-step-down")
  data.one.step.up <- data.one.step.down %>%
    mutate(y = -y, type = "One-step-up")

  # Data: two step
  x <- seq(-5,5,0.5)
  data.up.down <- data.frame(x = x) %>%
    mutate(y = -2*(x^2), type = "Up-down")
  data.down.up <- data.frame(x = x) %>%
    mutate(y = 2*(x^2), type = "Down-Up")

  df <- rbind(data.one.step.up, data.one.step.down,
              data.up.down,data.down.up )
  df$type <- factor(df$type, levels = c("One-step-up", "One-step-down",  "Up-down", "Down-Up" ))

  x <- y <- NULL
  blank <- ggplot2::element_blank
  ggplot(df, aes(x, y)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~type, scale = "free", ncol = 4)+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.ticks = blank(), axis.text = blank(),
                   axis.title = blank(), strip.background = blank())
}
