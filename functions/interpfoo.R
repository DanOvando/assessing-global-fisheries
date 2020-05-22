interpfoo <- function(data) {
  out <-
    approx(data$x, data$mean, xout = ceiling(min(data$x)):floor(max(data$x)))

  out <- data_frame(year = out$x, effort = out$y)
}