context("cpls")

suppressPackageStartupMessages(library(dplyr))
library(kernlab)
data(spam)

test_that("cpls", {
  x <- ifelse(spam %>% select(-starts_with("capital"), -type) > 0, 1, 0)
  y <- as.integer(spam$type) - 1 # 1: spam, 0: nonspam

  data <- data.frame(y, x)
  data <- data %>% sample_n(100)

  result <- cpls(y ~ ., data = data, ncomp = 3, verbose = FALSE)

  pred <- t(result$alpha) %*% t(result$w) %*% t(data %>% select(-y)) + result$center
  act <- data$y

  print(table(ifelse(pred >= 0.5, 1, 0), act) %>% prop.table %>% round(2))
})
