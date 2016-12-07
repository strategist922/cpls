context("cpls")

suppressPackageStartupMessages(library(dplyr))
library(kernlab)
data(spam)

test_that("cpls", {
  set.seed(71)
  x <- ifelse(spam %>% select(-starts_with("capital"), -type) > 0, 1, 0)
  y <- as.integer(spam$type) - 1 # 1: spam, 0: nonspam

  data <- data.frame(y, x)
  data <- data %>% sample_n(100)

  model <- cpls(y ~ ., data = data, ncomp = 3, verbose = FALSE)

  pred <- predict(model, data)
  act <- data$y

  print(table(ifelse(pred >= 0.5, 1, 0), act) %>% prop.table %>% round(2))
})
