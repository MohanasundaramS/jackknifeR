
# jackknifeR

<!-- badges: start -->
<!-- badges: end -->

The package jackknifeR includes the functions to perform delete-d Jackknife.

``` r
library(jackknifeR)
## basic example code
j.cor <- jackknife.cor(cars[, c("speed", "dist")], d = 2)
j.cor
```

