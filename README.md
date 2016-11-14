# An R Package for Compression-based Partial Least Squares Regression (cPLS)
Koji MAKIYAMA (@hoxo_m)  



[![Travis-CI Build Status](https://travis-ci.org/hoxo-m/cpls.svg?branch=master)](https://travis-ci.org/hoxo-m/cpls)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/cpls)](https://cran.r-project.org/package=cpls)
[![Coverage Status](https://coveralls.io/repos/github/hoxo-m/cpls/badge.svg?branch=master)](https://coveralls.io/github/hoxo-m/cpls?branch=master)

## An Example

### Prepare Data


```r
library(kernlab)
library(dplyr)
set.seed(71)
data(spam)

x <- ifelse(spam %>% select(-starts_with("capital"), -type) > 0, 1, 0)
y <- as.integer(spam$type) - 1 # 1: spam, 0: nonspam

data <- data.frame(y, x)
head(data %>% select(y:receive))
#>   y make address all num3d our over remove internet order mail receive
#> 1 1    0       1   1     0   1    0      0        0     0    0       0
#> 2 1    1       1   1     0   1    1      1        1     0    1       1
#> 3 1    1       0   1     0   1    1      1        1     1    1       1
#> 4 1    0       0   0     0   1    0      1        1     1    1       1
#> 5 1    0       0   0     0   1    0      1        1     1    1       1
#> 6 1    0       0   0     0   1    0      0        1     0    0       0

data <- data %>% sample_n(200)

train_test_indicator <- sample(c("train", "test"), size = nrow(data),
                               replace = TRUE, prob = c(0.6, 0.4))
d <- split(data, train_test_indicator)
```

### Execute cPLS


```r
library(cpls)

result <- cpls(y ~ ., data = d$train, ncomp = 3, verbose = FALSE)

pred <- t(result$alpha) %*% t(result$w) %*% t(d$train %>% select(-y)) + result$center
act <- d$train$y
RMSE <- sqrt(sum((as.numeric(pred) - as.numeric(act))^2))
RMSE
#> [1] 3.750782

table(ifelse(pred >= 0.5, 1, 0), act) %>% prop.table %>% round(2)
#>    act
#>        0    1
#>   0 0.49 0.04
#>   1 0.05 0.42

pred <- t(result$alpha) %*% t(result$w) %*% t(d$test %>% select(-y)) + result$center
act <- d$test$y
RMSE <- sqrt(sum((as.numeric(pred) - as.numeric(act))^2))
RMSE
#> [1] 3.145719

table(ifelse(pred >= 0.5, 1, 0), act) %>% prop.table %>% round(2)
#>    act
#>        0    1
#>   0 0.67 0.01
#>   1 0.06 0.26
```

## References

[1] Yasuo Tabei, Hiroto Saigo, Yoshihiro Yamanishi, and Simon J. Puglisi, "Scalable partial least squares regression on grammar-compressed data matrices,"  in *Proc. 22nd KDD*, 2016, pp. 1875-1884.
