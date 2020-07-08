sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
data <- (sim$data - 0.5)/0.5
data <- data + rnorm(length(data), 0, 1)
result <- mnem(data, k = 2, starts = 2)
resp <- getAffinity(result$probs, mw = result$mw, data = data)
resp1 <- apply(resp, 2, function(x) {
    xmax <- which(x == max(x))
    y <- x*0
    y[xmax] <- 1
    return(y)
    })
resp2 <- getAffinity(result$probs, affinity = 1, mw = result$mw, data = data)
checkEquals(resp1, resp2, checkNames = TRUE)
