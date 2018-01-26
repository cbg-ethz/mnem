# M&NEM

Install:
--------

Open R and input:

```r
install.packages("devtools")

library(devtools)

install_github("cbg-ethz/mnem")

library(mnem)
```

Small toy example on 1000 simulated cells and 100 E-genes.

```r

data <- matrix(rnorm(100*1000), 100, 1000)

colnames(data) <- sample(paste0("S", 1:5), 1000, replace = TRUE)

result <- mnem(data, k = 3)

plot(result)

```

