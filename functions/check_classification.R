check_classification <- function(bhat, b, breaks = c(0.8,1.2)) {

  breaks <- c(0,breaks,Inf)

  labels <- c("over","fully","under")

  comp <- data_frame(year = seq_along(b),b = b, bhat = bhat) %>%
    gather(metric, value, -year) %>%
    mutate(bin = cut(value, breaks = breaks, labels = labels)) %>%
    select(-value) %>%
    spread(metric,bin) %>%
    mutate(correct_class = b == bhat)

  accuracy <- mean(comp$correct_class)

}