
# jackknifeR

<!-- badges: start -->
<!-- badges: end -->

The goal of jackknifeR is to perform Jackknife on regression or correlation estimation.

``` r
library(jackknifeR)
## basic example code
j.cor <- jackknife.cor(cars$speed, cars$dist, d = 2)
j.cor$jackknife.summary
j.cor$biased_cor
```

