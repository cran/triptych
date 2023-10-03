recalibrate_mean <- function(x, y) {
  out <- numeric(length(x))
  ord <- order(x, -y)
  out[ord] <- monotone::monotone(y[ord])
  out
}

#' @importFrom class knn
recalibrate_mean2 <- function(x, y, k = round(length(y)^(1/3))) {
  knnfit <- class::knn(matrix(x), matrix(x), y, k = k, prob = TRUE)
  abs((knnfit == 0) - attr(knnfit, "prob"))
}
