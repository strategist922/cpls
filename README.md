# An R Package for Compression-based Partial Least Squares Regression (cPLS)
Koji MAKIYAMA (@hoxo_m)  



### This package is under construction.

## Examples


```r
library(kernlab)
library(dplyr)

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

d <- split(data, sample(0:1, size = nrow(data), replace = TRUE, prob = c(9, 1)))
names(d) <- c("train", "test")
```


```r
library(cpls)

result <- cpls(y ~ ., data = data, ncomp = 4)

pred <- t(result$alpha) %*% t(result$w) %*% t(d$train %>% select(-y)) + result$center
act <- d$train$y
RMSE <- sqrt(sum((as.numeric(pred) - as.numeric(act))^2))
RMSE
#> [1] 18.60759

table(ifelse(pred >= 0.5, 1, 0), act) %>% prop.table %>% round(2)
#>    act
#>        0    1
#>   0 0.57 0.04
#>   1 0.03 0.36

pred <- t(result$alpha) %*% t(result$w) %*% t(d$test %>% select(-y)) + result$center
act <- d$test$y
RMSE <- sqrt(sum((as.numeric(pred) - as.numeric(act))^2))
RMSE
#> [1] 5.936161

table(ifelse(pred >= 0.5, 1, 0), act) %>% prop.table %>% round(2)
#>    act
#>        0    1
#>   0 0.57 0.04
#>   1 0.04 0.35
```

## References

[1] Yasuo Tabei, Hiroto Saigo, Yoshihiro Yamanishi, and Simon J. Puglisi, "Scalable partial least squares regression on grammar-compressed data matrices,"  in *Proc. 22nd KDD*, 2016, pp. 1875--1884.
