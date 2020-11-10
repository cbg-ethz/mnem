#' Originally imported from the package 'nem'.
#' @noRd
#' @importFrom e1071 bincombinations
#' @author Florian Markowetz
enumerate.models <- function (x, name = NULL, trans.close = TRUE,
                              verbose = TRUE) {
    if (length(x) == 1) {
        n <- as.numeric(x)
        if (is.null(name)) 
            name <- letters[seq_len(n)]
    }
    else {
        n <- length(x)
        name <- x
    }
    if (n == 1) 
        stop("nem> choose n>1!")
    if (n > 5) 
        stop(paste0("nem> exhaustive enumeration not ",
                    "feasible with more than 5 perturbed genes"))
    if (n == 5) 
        cat("nem> this will take a while ... \n")
    bc <- e1071::bincombinations(n * (n - 1))
    fkt1 <- function(x, n, name) {
        M <- diag(n)
        M[which(M == 0)] <- x
        dimnames(M) <- list(name, name)
        if (trans.close) 
            M <- transitive.closure(M)
        return(list(M))
    }
    models <- apply(bc, 1, fkt1, n, name)
    models <- unique(matrix(unlist(models), ncol = n * n, byrow = TRUE))
    fkt2 <- function(x, n, name) {
        M <- matrix(x, n)
        dimnames(M) <- list(name, name)
        return(list(M))
    }
    models <- unlist(apply(models, 1, fkt2, n, name), recursive = FALSE)
    if (verbose) 
        cat("Generated", length(models), "unique models ( out of", 
            2^(n * (n - 1)), ")\n")
    return(models)
}
#' Originally imported from the package 'nem'.
#' @noRd
#' @author Holger Froehlich, Cordula Zeller
sampleRndNetwork <- function (Sgenes, scaleFree = TRUE, gamma = 2.5,
                              maxOutDegree = length(Sgenes), 
    maxInDegree = length(Sgenes), trans.close = TRUE, DAG = FALSE) {
    n = length(Sgenes)
    S = diag(n)
    maxOutDegree = min(n, maxOutDegree)
    degprob = (0:maxOutDegree)^(-gamma)
    degprob[1] = 1
    degprob = degprob/sum(degprob)
    for (i in seq_len(n)) {
        if (scaleFree) 
            outdeg = sample(0:maxOutDegree, 1, prob = degprob)
        else outdeg = sample(0:maxOutDegree, 1)
        if (outdeg > 0) {
            if (!DAG) 
                idx0 = which(S[i, ] == 0)
            else idx0 = which(S[i, ] == 0 & seq_len(n) < i)
            if (length(idx0) > 0) {
                idx = which(colSums(S[, idx0, drop = FALSE]) <= 
                  maxInDegree)
                if (length(idx) > 0) {
                  idx = sample(idx0[idx], min(outdeg, length(idx0[idx])), 
                    replace = TRUE)
                  S[i, idx] = 1
                }
            }
        }
    }
    if (trans.close) 
        S = transitive.closure(S)
    diag(S) = 0
    colnames(S) = Sgenes
    rownames(S) = Sgenes
    S
}
#' @noRd
#' @importFrom methods is
bigphi <- function(x) {
    if (is(x, "mnemsim")) {
        resfull <- NULL
        for (l in seq_len(length(x$Nem))) {
            tmp <- transitive.closure(x$Nem[[l]])
            resfull <- cbind(resfull, t(tmp))
        }
    } else {
        resfull <- NULL
        for (l in seq_len(length(x$comp))) {
            tmp <- transitive.closure(x$comp[[l]]$phi)
            colnames(tmp) <- rownames(tmp) <- seq_len(nrow(tmp))
            resfull <- cbind(resfull, t(tmp))
        }
    }
    return(resfull)
}
#' @noRd
Kratio <- function(x, y) {
    xn <- length(unlist(x$comp, recursive = FALSE))/2
    yn <- length(y$Nem)
    z <- xn/yn
    return(z)
}
#' @noRd
uniques <- function(x) {
    for (i in seq_len(length(x$comp))) {
        x$comp[[i]]$theta <- NULL
    }
    y <- length(unique(x$comp))
    return(y)
}
#' @noRd
mnemh.rec <- function(data, k = 2, logtype = 2, ...) {
    tmp <- mnemk(data, ks=seq_len(k), logtype = logtype, ...)
    cluster <- apply(getAffinity(tmp$best$probs,
                                 logtype = logtype,
                                 mw = tmp$best$mw), 2, which.max)
    if (length(tmp$best$comp) == 1) {
        return(tmp$best)
    } else {
        tmp1 <- NULL
        for (i in seq_len(length(tmp$best$comp))) {
            if (length(unique(colnames(data[, which(cluster == i)]))) > 1) {
                tmp2 <- mnemh.rec(data[, which(cluster == i)],
                                  k=k, logtype=logtype, ...)
            } else {
                tmp2 <- NULL
            }
            tmp1 <- c(tmp1, tmp2)
        }
        return(tmp1)
    }
}
#' @noRd
random_probs <- function(k, data, full = FALSE, logtype = 2) {
    samplefun <- function(n,p) {
        x <- sample(c(0.5,1.5), n,
               replace = TRUE,
               prob = p)

        return(x)
    }
    probs <- matrix(log(samplefun(k*ncol(data),
                                  c(0.9, 0.1)))/log(logtype), k,
                    ncol(data))
    if (full) {
        for (i in seq_len(k)) {
            if (i == 1) { next() }
            infcells <- which(apply(probs, 2, function(x) {
                bad <- FALSE
                if (all(is.infinite(x))) {
                    bad <- TRUE
                }
                return(bad)
            }))
            if (i == k) {
                probs[i, infcells] <- log(1)/log(logtype)
            } else {
                probs[i, infcells] <-
                    log(samplefun(length(infcells),
                                  c(1-1/k, 1/k)))/log(logtype)
            }
        }
    }
    while(any(apply(probs, 1, function(x) {
        bad <- FALSE
        if (all(is.infinite(x)) | all(x == 0)) {
            bad <- TRUE
        }
        return(bad)
    }))) {
        probs <- matrix(log(samplefun(k*ncol(data),
                                      c(0.9, 0.1)))/log(logtype), k,
                        ncol(data))
        if (full) {
            for (i in seq_len(k)) {
                if (i == 1) { next() }
                infcells <- which(apply(probs, 2, function(x) {
                    bad <- FALSE
                    if (all(is.infinite(x))) {
                        bad <- TRUE
                    }
                    return(bad)
                }))
                if (i == k) {
                    probs[i, infcells] <- log(1)/log(logtype)
                } else {
                    probs[i, infcells] <-
                        log(samplefun(length(infcells),
                                      c(1-1/k, 1/k)))/log(logtype)
                }
            }
        }
    }
    return(probs)
}
#' @noRd
random_probs2 <- function(k, data, full = FALSE, logtype = 2) {
    samplefun <- function(n,p) {
        x <- sample(c(0.5,1.5), n,
                    replace = TRUE,
                    prob = p)

        return(x)
    }
    probs <- matrix(0, k, ncol(data))
    for (i in seq_len(k)) {
        if (i == 1) {
            probs[i, ] <- samplefun(ncol(data), c(1-1/k,1/k))
        } else {
            probs[i, which(probs[i-1, ] == 0)] <-
                samplefun(sum(probs[i-1, ] == 0), c(1-1/(k+1-i),1/(k+1-i)))
        }
    }
    return(probs)
}
#' @noRd
sortAdj <- function(res, list = FALSE) {
    resmat <- NULL
    for (i in seq_len(length(res))) {
        if (list) {
            resmat <-
                rbind(resmat,
                      as.vector(mytc(res[[i]])))
        } else {
            resmat <-
                rbind(resmat,
                      as.vector(mytc(res[[i]]$adj)))
        }
    }
    d <- as.matrix(dist(resmat))
    dsum <- apply(d, 1, sum)
    resorder <- which.max(dsum)
    d[resorder, resorder] <- Inf
    for (i in 2:length(dsum)) {
        resorder <- c(resorder, which.min(d[, resorder[length(resorder)]]))
        d[resorder, resorder] <- Inf
    }
    res2 <- list()
    for (i in seq_len(length(res))) {
        res2[[i]] <- res[[resorder[i]]]
    }
    return(list(res = res2, order = resorder))
}
#' @noRd
calcEvopen <- function(res, list = TRUE) {
    evopen <- 0
    for (i in seq_len(length(res)-1)) {
        if (list) {
            evopen <- evopen + sum(abs(res[[i]] -
                                       res[[(i+1)]]))/length(res[[i]])
        } else {
            evopen <- evopen + sum(abs(res[[i]]$adj -
                                       res[[(i+1)]]$adj))/length(res[[i]]$adj)
        }
    }
    evopen <- -evopen#/(k-1)
    return(evopen)
}
#' @noRd
transProb <- function(x, y, lambda = 0.5) {
    n <- nrow(x)
    exp <- sum(abs(x-y))
    enum <- lambda*(1-lambda)^(exp-1)
    denom <- 0
    ne <- (n*(n-1))
    for (i in 0:ne) {
        denom <- denom + choose(ne, i)*lambda*(1-lambda)^(i-1)
    }
    p <- enum/denom
    return(p)
}
#' @noRd
fullTransProb <- function(x, ...) {
    p <- 1
    for (i in seq_len(length(x)-1)) {
        p <- p*transProb(x[[i]], x[[i+1]], ...)
    }
    return(p)
}
#' @noRd
modAdj <- function(adj, D) {
    Sgenes <- getSgenes(D)
    SgeneN <- getSgeneN(D)
    for (i in seq_len(SgeneN)) {
        colnames(adj) <- rownames(adj) <- gsub(i, Sgenes[i], colnames(adj))
    }
    return(adj)
}
#' @noRd
getRho <- function(data) {

    Sgenes <- getSgenes(data)

    Rho <- matrix(0, length(Sgenes), ncol(data))

    for (i in seq_len(length(Sgenes))) {
        Rho[i, grep(paste0("^", Sgenes[i], "_|_", Sgenes[i],
                           "$|_", Sgenes[i], "_|^", Sgenes[i], "$"),
                    colnames(data))] <- 1
    }

    rownames(Rho) <- Sgenes
    colnames(Rho) <- colnames(data)
    Rho <- Rho[naturalsort(rownames(Rho)), ]
    return(Rho)
}
#' @noRd
initComps <- function(data, k=2, starts=1, verbose = FALSE, meanet = NULL,
                      linets = TRUE) {
    n <- getSgeneN(data)
    init <- list()
    for (i in seq_len(starts)) {
        init[[i]] <- list()
        for (j in seq_len(k)) {
            if (linets) {
                tmp <- matrix(0, n, n)
                tmp[upper.tri(tmp)] <- 1
                diag(tmp) <- 1
                rownames(tmp) <- colnames(tmp) <- sample(seq_len(n), n)
                tmp <- tmp[naturalorder(rownames(tmp)),
                           naturalorder(colnames(tmp))]
            } else {
                tmp <- sampleRndNetwork(seq_len(n), DAG = TRUE)
            }
            init[[i]][[j]] <- tmp
        }
    }
    return(init)
}
#' @noRd
initps <- function(data, ks, k, starts = 3, ksel = "dist") {
    clusters <- list()
    multi <- grep("_", colnames(data))
    if (length(multi) > 0) {
        data2 <- data[, -multi]
    } else {
        data2 <- data
    }
    for (i in seq_len(length(unique(colnames(data2))))) {
        if ("cor" %in% ksel) {
            d <- (1 - cor(data[, which(colnames(data) %in% i)]))/2
            d <- as.dist(d)
        } else {
            d <- dist(t(data[, which(colnames(data) %in% i)]))
        }
        if (length(d) > 1) {
            hc <- hclust(d)
            clusters[[i]] <- cutree(hc, min(ks[i], length(hc$labels)))
        } else {
            clusters[[i]] <- 1
        }
    }
    probscl <- list()

    n <- getSgeneN(data)

    count <- 0
    llstr <- NULL
    resstr <- NULL
    counter <- 0
    takes <- NULL
    takesdone <- NULL
    while(count < starts & counter < prod(ks)) {
        counter <- counter + 1
        tmp <- matrix(0.5, k, ncol(data))
        if (ks[1] < k) {
            takes <- as.matrix(sample(seq_len(ks[1]), k, replace = TRUE), k, 1)
        } else {
            takes <- as.matrix(sample(seq_len(ks[1]), k), k, 1)
        }
        for (i in 2:length(unique(colnames(data)))) {
            if (ks[i] < k) {
                takes <- cbind(takes, sample(seq_len(ks[i]), k, replace = TRUE))
            } else {
                takes <- cbind(takes, sample(seq_len(ks[i]), k))
            }
        }
        for (i in seq_len(k)) {
            for (j in seq_len(n)) {
                tmp[i, which(colnames(data) == j)[
                           which(clusters[[j]] == takes[i, j])]] <- 1.5
            }
        }
        takestmp <- paste(takes, collapse = "")
        if (!(takestmp %in% takesdone)) {
            count <- count + 1
            takesdone <- c(takesdone, takestmp)
            probscl[[count]] <- log2(tmp)
        }
    }
    return(probscl)
}
#' @noRd
modData <- function(D) {
    SgeneN <- getSgeneN(D)
    Sgenes <- getSgenes(D)
    if (!all(is.numeric(Sgenes))) {
        colnamesD <- colnames(D)
        for (i in seq_len(SgeneN)) {
            colnamesD <- gsub(paste0("^", Sgenes[i], "$"), i, colnamesD)
            colnamesD <- gsub(paste0("_", Sgenes[i], "$"),
                              paste0("_", i), colnamesD)
            colnamesD <- gsub(paste0("^", Sgenes[i], "_"),
                              paste0(i, "_"), colnamesD)
            colnamesD <- gsub(paste0("_", Sgenes[i], "_"),
                              paste0("_", i, "_"), colnamesD)
        }
        colnames(D) <- colnamesD
    }
    rownames(D) <- as.numeric(seq_len(nrow(D)))
    return(D)
}
#' @noRd
#' @importFrom cluster silhouette
#' @importFrom stats kmeans
learnk <- function(data, kmax = 10, ksel = c("hc", "silhouette", "cor"),
                   starts = 10, mono = TRUE, verbose = FALSE) {
    ks <- numeric(length(unique(colnames(data))))
    lab <- list()
    index <- numeric(ncol(data))
    for (i in naturalsort(as.numeric(unique(colnames(data))))) {
        if (verbose) {
            print(i)
        }
        if (sum(colnames(data) %in% i) <= 1) { k <- 1; next() }
        if ("cor" %in% ksel) {
            d <- (1 - cor(data[, which(colnames(data) %in% i)]))/2
            d <- as.dist(d)
        } else {
            d <- dist(t(data[, which(colnames(data) %in% i)]))
        }
        if ("hc" %in% ksel) {
            hc <- hclust(d)
            ks[i] <- 2
            lab[[i]] <- rep(1, sum(colnames(data) %in% i))
            if (length(d) > 1) {
                silavg <- 0
                silavgs <- numeric(length(hc$order)-1)
                clusters <- list()
                for (j in 2:(length(hc$order)-1)) {
                    cluster <- cutree(hc, j)
                    clusters[[j]] <- cluster
                    if ("BIC" %in% ksel | "AIC" %in% ksel) {
                        m <- nrow(data)
                        n <- length(cluster)
                        k <- length(table(cluster))
                        D <- log(mean(silhouette(cluster, d)[, 3]))
                        if ("AIC" %in% ksel) {
                            silavgs[j] <- 2*m*k - 2*D
                        } else {
                            silavgs[j] <- log(n)*m*k - 2*D
                        }
                    }
                    if ("silhouette" %in% ksel) {
                        silavgs[j] <- mean(silhouette(cluster, d)[, 3])
                    }
                    if (verbose) {
                        print(silavgs[j])
                    }
                    if (silavgs[j] < silavgs[(j-1)]) {
                        break()
                    }
                    if (silavg < silavgs[j]) {
                        silavg <- silavgs[j]
                        ks[i] <- j
                        lab[[i]] <- cluster
                    }
                }
            }
            index[which(colnames(data) %in% i)] <- lab[[i]]
        }
        if ("kmeans" %in% ksel) {
            silavg <- 0
            silavgs <- numeric(kmax)
            clusters <- list()
            for (j in seq_len(kmax)) {
                if (j == 1) { next() }
                if (ncol(as.matrix(d)) <= j) { next() }
                kc <- kmeans(d, centers = j, nstart = starts)
                if ("BIC" %in% ksel | "AIC" %in% ksel) {
                    m <- ncol(kc$centers)
                    n <- length(kc$cluster)
                    k <- nrow(kc$centers)
                    D <- kc$tot.withinss
                    if ("AIC" %in% ksel) {
                        silavgs[j] <- 2*m*k + D
                    } else {
                        silavgs[j] <- log(n)*m*k + D
                    }
                }
                if ("silhouette" %in% ksel) {
                    silavgs[j] <- mean(silhouette(kc$cluster, d)[, 3])
                }
                if (verbose) {
                    print(silavgs[j])
                }
                if (j == 1) {
                    j2 <- j
                } else {
                    j2 <- j - 1
                }
                if (silavgs[j] < silavgs[j] & mono) {
                    break()
                }
                if (silavg < silavgs[j]) {
                    silavg <- silavgs[j]
                    ks[i] <- j
                    lab[[i]] <- kc$cluster
                }
            }
            index[which(colnames(data) %in% i)] <- lab[[i]]
        }
    }
    k <- min(kmax, max(ks))
    return(list(ks = ks, k = k, lab = lab, cluster = index))
}
#' @noRd
getLL <- function(x, logtype = 2, mw = NULL, data = NULL, complete = FALSE) {
    if (is.null(mw)) { mw = rep(1, nrow(x))/nrow(x) }
    if (complete) {
        Z <- getAffinity(x, logtype = logtype, mw = mw, data = data,
                         complete = complete)
        l <- sum(apply(Z*(x + log(mw)/log(logtype)), 2, sum))
    } else {
        x <- logtype^x
        x <- x*mw
        l <- sum(log(rep(1,nrow(x))%*%x)/log(logtype))
    }
    return(l)
}
#' @noRd
estimateSubtopo <- function(data, fun = which.max) {
    if (length(grep("_", colnames(data))) > 0) {
        data <- data[, -grep("_", colnames(data))]
    }
    effectsums <- effectsds <- matrix(0, nrow(data),
                                      length(unique(colnames(data))))
    n <- getSgeneN(data)
    for (i in seq_len(n)) {
        ilen <- length(grep(i, colnames(data)))
        if (ilen > 1) {
            effectsds[, i] <- apply(data[, grep(i, colnames(data))], 1, sd)
            effectsums[, i] <- apply(data[, grep(i, colnames(data))], 1, sum)
        } else if (ilen == 1) {
            effectsds[, i] <- 1
            effectsums[, i] <- data[, grep(i, colnames(data))]
        }
    }
    subtopoX <- as.numeric(apply(effectsums/effectsds, 1, fun))
    subtopoX[which(is.na(subtopoX) == TRUE)] <- 1
    return(subtopoX)
}
#' @noRd
getProbs <- function(probs, k, data, res, method = "llr", affinity = 0,
                     converged = 10^-2, subtopoX = NULL, ratio = TRUE,
                     logtype = 2, mw = NULL, fpfn = fpfn, Rho = NULL,
                     complete = FALSE) {
    n <- ncol(res[[1]]$adj)
    bestprobs <- probsold <- probs
    time0 <- TRUE
    count <- 0
    if (ncol(data) <= 100) {
        max_count <- 100
    } else {
        max_count <- 10
    }
    ll0 <- llold <- -Inf
    stop <- FALSE
    mw <- apply(getAffinity(probsold, affinity = affinity, norm = TRUE,
                            mw = mw, logtype = logtype, data = data,
                            complete = complete), 1, sum)
    mw <- mw/sum(mw)
    if (any(is.na(mw))) { mw <- rep(1, k)/k }
    while((!stop | time0) & count < max_count) {
        time0 <- FALSE
        probsold <- probs
        subtopo0 <- matrix(0, k, nrow(data))
        subweights0 <- matrix(0, nrow(data), n+1)
        postprobsold <- getAffinity(probsold, affinity = affinity, norm = TRUE,
                                    logtype = logtype, mw = mw, data = data,
                                    complete = complete)
        probs0 <- probsold*0
        for (i in seq_len(k)) {
            adj1 <- mytc(res[[i]]$adj)
            if (is.null(Rho)) {
                adj1 <- adj1[colnames(data), ]
            } else {
                adj1 <- t(Rho)%*%adj1
                adj1[which(adj1 > 1)] <- 1
            }
            adj1 <- cbind(adj1, "0" = 0)
            if (is.null(res[[i]]$subtopo)) {
                subtopo <- maxCol_row(t(t(data)*postprobsold[i, ])%*%adj1)
            } else {
                subtopo <- res[[i]]$subtopo
                subtopo[which(subtopo == 0)] <- ncol(adj1)
            }
            adj2 <- adj1[, subtopo]
            if (method %in% "llr") {
                tmp <- colSums(data*t(adj2))
            }
            probs0[i, ] <- tmp
        }
        ll0 <- getLL(probs0, logtype = logtype, mw = mw,
                     data = data, complete = complete)
        if (max(ll0) - llold > 0) {
            bestprobs <- probs0
        }
        probs <- probs0
        if (max(ll0) - llold <= converged) {
            stop <- TRUE
        }
        mw <- apply(getAffinity(probs, affinity = affinity, norm = TRUE,
                                logtype = logtype, mw = mw, data = data,
                                complete = complete), 1,
                    sum)
        mw <- mw/sum(mw)
        llold <- max(ll0)
        count <- count + 1
    }
    return(list(probs = bestprobs, subtopoX = subtopoX,
                mw = mw, ll = llold))
}
#' @noRd
annotAdj <- function(adj, data) {
    Sgenes <- getSgenes(data)
    colnames(adj) <- rownames(adj) <- naturalsort(Sgenes)
    return(adj)
}
#' @noRd
nemEst <- function(data, maxiter = 100, start = "null",
                   cut = 0, monoton = TRUE, logtype = 2,
                   method = "llr",
                   weights = NULL, fpfn = c(0.1, 0.1), Rho = NULL,
                   close = TRUE, domean = TRUE, modified = FALSE,
                   hierarchy = "totaleffect", ...) {
    if (method %in% "disc") {
        D <- data
        D[which(D == 1)] <- log((1-fpfn[2])/fpfn[1])/log(logtype)
        D[which(D == 0)] <- log(fpfn[2]/(1-fpfn[1]))/log(logtype)
        method <- "llr"
        data <- D
    }
    if (!modified) {
        data <- modData(data)
    }
    if (sum(duplicated(colnames(data)) == TRUE) > 0 &
        method %in% "llr" & domean) {
        if (!is.null(weights)) {
            D <- data
            data <- data*rep(weights, rep(nrow(data), ncol(data)))
        }
        data2 <- doMean(data, weights = weights, Rho = Rho, logtype = logtype)
        data <- D
    } else {
        data2 <- data
    }
    if (is.null(weights)) { weights <- rep(1, ncol(data2)) }
    R <- data2[, naturalsort(colnames(data2))]
    Rho <- diag(ncol(R))
    N <- rowSums(Rho)
    n <- getSgeneN(R)
    phibest <- phi <- matrix(0, n, n)
    Sgenes <- getSgenes(R)
    rownames(phi) <- colnames(phi) <- Sgenes
    if (hierarchy == "totaleffect") {
        E0 <- apply(t(t(R%*%t(Rho))/N), 2, sum)
        phi <- phi[order(E0, decreasing = TRUE), order(E0, decreasing = TRUE)]
        phi[upper.tri(phi)] <- 1
        phi <- phi[naturalsort(rownames(phi)), naturalsort(colnames(phi))]
        E <- phi
    } else if (hierarchy == "pairwise.greedy") {
        for (i in seq_len(n-1)) {
            for (j in (i+1):n) {
                pwg <- nem(R[, c(i,j)], logtype = logtype,
                             weights = weights[c(i,j)], trans.close = close,
                             domean = FALSE)
                phi[i, j] <- pwg$adj[1, 2]
                phi[j, i] <- pwg$adj[2, 1]
            }
        }
        E0 <- E <- phi
    } else {
        Rpos <- t(t(R%*%t(Rho))/N)
        Rpos[which(Rpos < 0)] <- 0
        E0 <- t(Rpos)%*%(R%*%t(Rho))
        E <- E0  - t(E0)
        parents <- which(E <= 0)
        children <- which(E > 0)
        E[parents] <- 1
        E[children] <- 0
        E <- E[naturalorder(rownames(E)), naturalorder(colnames(E))]
        phi <- phi[naturalsort(rownames(phi)), naturalsort(colnames(phi))]
    }
    if ("full" %in% start) {
        phi <- phi
        diag(phi) <- 1
    } else if ("rand" %in% start) {
        phi <- phi*0
        diag(phi) <- 1
        phi[seq_len(length(phi))] <- sample(c(0,1), length(phi), replace = TRUE)
        phi[lower.tri(phi)] <- 0
        phi <- phi[sample(seq_len(nrow(phi)), nrow(phi)),
                   sample(seq_len(nrow(phi)), nrow(phi))]
        phi <- phi[naturalsort(rownames(phi)), naturalsort(colnames(phi))]
    } else if ("null" %in% start) {
        phi <- phi*0
        diag(phi) <- 1
    } else {
        phi <- start
    }
    O <- phi*0
    iter <- Oold <- 0
    lls <- NULL
    llbest <- -Inf
    stop <- FALSE
    while(!stop & iter < maxiter) {
        iter <- iter + 1
        ll <- scoreAdj(R, phi,
                       weights = weights, dotopo = TRUE,
                       trans.close = close, Rho = Rho)
        P <- ll$subweights
        theta <- theta2theta(ll$subtopo, phi)
        ll <- ll$score
        Oold <- O
        if (ll %in% lls | all(phi == phibest)) {
            stop <- TRUE
        }
        if (monoton & iter > 1) {
            if (ll < lls[length(lls)]) { stop <- TRUE }
        }
        if (llbest < ll) {
            phibest <- phi
            thetabest <- theta
            llbest <- ll
            numbest <- iter
            Obest <- O
        }
        lls <- c(lls, ll)
        nogenes <- which(apply(theta, 1, sum) == 0)
        nozeros <- which(t(P) > 0, arr.ind = TRUE)
        nozeros <- nozeros[which(nozeros[, 1] %in% nogenes), ]
        theta[nozeros] <- 1
        O <- (t(R%*%t(Rho))*weights)%*%t(theta)
        cutoff <- cut*max(abs(O))
        phi[which(O > cutoff & E == 1)] <- 1
        phi[which(O <= cutoff | E == 0)] <- 0
        if (close) {
            phi <- mytc(phi)
        }
    }
    phintc <- phibest
    if (close) {
        phibest <- mytc(phibest)
    }
    ll <- scoreAdj(R, phibest,
                   weights = weights, dotopo = TRUE,
                       trans.close = close, Rho = Rho)
    P <- ll$subweights
    theta <- theta2theta(ll$subtopo, phibest)
    llbest <- ll$score
    nem <- list(phi = phibest, theta = thetabest, iter = iter,
                ll = llbest, lls = lls, num = numbest,
                O = Obest, E = E0, phintc = phintc)
    class(nem) <- "nemEst"
    return(nem)
}
#' @noRd
doMean <- function(D, weights = NULL, Rho = NULL, logtype = 2) {
    man <- FALSE
    if (man) {
        mD <- matrix(0, nrow(D), length(unique(colnames(D))))
        if (!is.null(weights)) {
            doodds <- FALSE
            if (doodds) {
                A <- exp(D)/(1+exp(D))
                D <- log(t(t(A)*weights +
                           (1-weights)*0.5)/t(t(1 - A)*weights +
                                              (1-weights)*0.5))/log(logtype)
            } else {
            D <- D*rep(weights, rep(nrow(D), ncol(D)))
            }
        }
        for (i in seq_len(length(unique(colnames(D))))) {
            mD[, i] <-
                apply(D[, which(colnames(D) %in% unique(colnames(D))[i]),
                        drop = FALSE], 1, mean)
        }
        colnames(mD) <- unique(colnames(D))
    } else {
        if (!is.null(weights)) {
            D <- D*rep(weights, rep(nrow(D), ncol(D)))
        }
        if (!is.null(Rho)) {
            Rho <- apply(Rho, 2, function(x) {
                if (sum(x) > 1) {
                    x <- x/sum(x)
                }
                return(x)
            })
            Rho[is.na(Rho)] <- 0
            mD <- D%*%t(Rho)
            mD <- mD[, naturalorder(colnames(mD))]
            colnames(mD) <- seq_len(ncol(mD))
        } else {
            global <- TRUE
            if (global) {
                mD <- t(rowsum(t(D), colnames(D)))/
                    (ncol(D)/length(unique(colnames(D))))
            } else {
                mD <- t(rowsum(t(D), colnames(D))/
                        as.numeric(table(colnames(D))))
            }
        }
    }
    mD <- mD[, naturalorder(colnames(mD))]
    return(mD)
}
#' @noRd
modules <- function(D, method = "llr", weights = NULL, reduce = FALSE,
                    start = NULL,
                    verbose = FALSE, trans.close = TRUE, redSpace = NULL,
                    subtopo = NULL, ratio = TRUE, parallel = NULL,
                    prior = NULL, fpfn = c(0.1, 0.1),
                    modulesize = 4, search = "exhaustive", domean = TRUE,
                    Rho = NULL, logtype = 2) {
    D <- data <- modData(D)
    n <- getSgeneN(D)
    Sgenes <- getSgenes(D)
    if (domean) {
        D <- doMean(D, weights = weights, Rho = Rho, logtype = logtype)
        weights <- rep(1, ncol(D))
        sumdata <- data <- D
    } else {
        sumdata <- matrix(0, nrow(data), n)
        if (!is.null(weights)) {
            D <- D*rep(weights, rep(nrow(D), ncol(D)))
            weights <- rep(1, ncol(sumdata))
        }
        for (i in seq_len(n)) {
            sumdata[, i] <-
                apply(D[, which(colnames(D) %in% i), drop = FALSE], 1, sum)
        }
        colnames(sumdata) <- seq_len(n)
        rownames(sumdata) <- rownames(D)
    }
    D <- NULL
    n <- getSgeneN(data)
    cordata <- cor(sumdata)
    cordata[is.na(cordata)] <- -1
    d <- as.dist((1 - cordata)/2)
    for (i in 2:n) {
        hc <- hclust(d)
        hcut <- cutree(hc, i)
        if (max(table(hcut)) <= modulesize) {
            break()
        }
    }
    groups <- table(hcut)
    fullnet <- NULL
    for (i in seq_len(length(groups))) {
        subset <- which(hcut == i)
        if (verbose) {
            print(paste(c("calculating module", subset), collapse = " "))
        }
        if (length(subset) > 1) {
            subdata <- data[, which(colnames(data) %in% subset)]
            if (is.null(start)) {
                start2 <- start
            } else {
                start2 <- start[which(rownames(start) %in% subset),
                                which(colnames(start) %in% subset)]
            }
            if (!is.null(Rho)) { Rho <- getRho(subdata) }
            tmp <- nem(subdata, search = search, method = method,
                         start = start2,
                         parallel = parallel, reduce = reduce,
                         weights = weights[which(colnames(data) %in% subset)],
                         verbose = verbose,
                         redSpace = redSpace, trans.close = trans.close,
                         subtopo = subtopo, prior = prior, ratio = ratio,
                         domean = FALSE, fpfn = fpfn, Rho = Rho,
                         logtype = logtype)
            if (is.null(fullnet)) {
                fullnet <- tmp$adj
            } else {
                tmpnames <- c(colnames(fullnet), colnames(tmp$adj))
                fullnet <-
                    rbind(cbind(fullnet,
                                matrix(0, nrow(fullnet), ncol(tmp$adj))),
                          cbind(matrix(0, nrow(tmp$adj), ncol(fullnet)),
                                tmp$adj))
                colnames(fullnet) <- rownames(fullnet) <- as.numeric(tmpnames)
            }
        } else {
            if (is.null(dim(fullnet))) {
                fullnet <- matrix(1, 1, 1)
                colnames(fullnet) <- rownames(fullnet) <- subset
            } else {
                fullnet <- rbind(cbind(fullnet, 0), 0)
                colnames(fullnet)[ncol(fullnet)] <-
                    rownames(fullnet)[nrow(fullnet)] <-
                    subset
            }
        }
    }
    fullnet <- transitive.reduction(fullnet)
    fullnet <- fullnet[order(as.numeric(rownames(fullnet))),
                       order(as.numeric(colnames(fullnet)))]
    return(fullnet)
}
#' @noRd
getSgeneN <- function(data) {
    Sgenes <- length(unique(unlist(strsplit(colnames(data), "_"))))
    return(Sgenes)
}
#' @noRd
getSgenes <- function(data) {
    Sgenes <- naturalsort(unique(unlist(strsplit(colnames(data), "_"))))
    return(Sgenes)
}
#' @noRd
get.rev.tc <-function (Phi) {
    Phitr <- transitive.reduction(Phi)
    idx = which(Phitr + t(Phitr) == 1, arr.ind = TRUE)
    models = list()
    nn <- dim(Phi)
    if (NROW(idx) > 0) {
        for (i in seq_len(NROW(idx))) {
            Phinew = Phi
            Phinew[idx[i, 1], idx[i, 2]] = 1 - Phinew[idx[i, 1], idx[i, 2]]
            Phinew[idx[i, 2], idx[i, 1]] = 1 - Phinew[idx[i, 2], idx[i, 1]]
            diag(Phinew) = 1
            if (Phinew[idx[i, 1], idx[i, 2]] == 1) {
                uv <- idx[i, ]
            } else {
                uv <- rev(idx[i, ])
            }
            Phinew <- mytc(Phinew, uv[1], uv[2])
            models[[i]] <- Phinew
        }
    }
    models
}
#' @noRd
get.ins.fast <- function (Phi, trans.close = TRUE) {
    idx = which(Phi == 0)
    models = list()
    nn <- dim(Phi)
    if (length(idx) > 0) {
        for (i in seq_len(length(idx))) {
            uv <- arrayInd(idx[i], nn)
            Phinew = Phi
            Phinew[idx[i]] = 1
            if (trans.close) {
                Phinew = mytc(Phinew, uv[1], uv[2])
            }
            models[[i]] <- Phinew
        }
    }
    models
}
#' @noRd
get.del.tc <- function (Phi) {
    Phi = Phi - diag(ncol(Phi))
    Phi2 <- transitive.reduction(Phi)
    idx = which(Phi2 == 1)
    models = list()
    if (length(idx) > 0) {
        for (i in seq_len(length(idx))) {
            Phinew = Phi
            Phinew[idx[i]] = 0
            diag(Phinew) = 1
            models[[i]] <- Phinew
        }
    }
    models
}
#' @noRd
theta2theta <- function(x, y) {
    if (is.matrix(x)) {
        z <- apply(x, 2, function(x) {
            if (max(x) == 0) {
                y <- 0
            } else {
                y <- which.max(x)
            }
            return(y)
        })
    } else {
        z <- matrix(0, nrow(y), length(x))
        z[cbind(x, seq_len(ncol(z)))] <- 1
        rownames(z) <- seq_len(nrow(z))
        colnames(z) <- seq_len(ncol(z))
    }
    return(z)
}
#' @noRd
adj2dnf <- function(A) {

    dnf <- NULL

    for (i in seq_len(ncol(A))) {
        for (j in seq_len(nrow(A))) {
            if (A[i, j] == 1) {
                dnf <- c(dnf, paste(colnames(A)[i], rownames(A)[j], sep = "="))
            }
            if (A[i, j] == -1) {
                dnf <- c(dnf, paste("!", colnames(A)[i], "=",
                                    rownames(A)[j], sep = ""))
            }
        }
    }

    dnf <- unique(dnf)

    return(dnf)

}
#' @noRd
#' @importFrom methods new
plot.adj <- function(x, ...) {
    adj2graph <- function(adj.matrix) {
        V   <- rownames(adj.matrix)
        edL <- vector("list", length=nrow(adj.matrix))
        names(edL) <- V
        for (i in seq_len(nrow(adj.matrix))) {
            edL[[i]] <- list(edges=which(!adj.matrix[i,]==0),
                             weights=adj.matrix[i,!adj.matrix[i,]==0])
        }
        gR <- new("graphNEL",nodes=V,edgeL=edL,edgemode="directed")
        return(gR)
    }
    g <- adj2graph(x)
    plot(g)
}
#' @noRd
graph2adj <- function(gR) {
    adj.matrix <- matrix(0,
                         length(nodes(gR)),
                         length(nodes(gR))
                         )
    rownames(adj.matrix) <- nodes(gR)
    colnames(adj.matrix) <- nodes(gR)
    for (i in seq_len(length(nodes(gR)))) {
        adj.matrix[nodes(gR)[i],adj(gR,nodes(gR)[i])[[1]]] <- 1
    }
    return(adj.matrix)
}
#' @noRd
#' @importFrom flexclust dist2
#' @import Rcpp
#' @import RcppEigen
llrScore <- function(data, adj, weights = NULL, ratio = TRUE) {
    if (is.null(weights)) {
        weights <- rep(1, ncol(data))
    }
    if (ratio) {
        score <- eigenMapMatMult(data, adj*weights)
    } else {
        if (max(data) == 1) {
            score <- -dist2(data, t(adj)*weights)
        } else {
            score <- -dist2(data, t((adj*mean(data))*weights))
        }
    }
    return(score)
}
#' @noRd
adj2dnf <- function(A) {
    dnf <- NULL
    for (i in seq_len(ncol(A))) {
        dnf <- c(dnf, rownames(A))
        for (j in seq_len(nrow(A))) {
            if (A[i, j] == 1) {
                dnf <- c(dnf, paste(colnames(A)[i], rownames(A)[j], sep = "="))
            }
            if (A[i, j] == -1) {
                dnf <- c(dnf, paste("!", colnames(A)[i], "=", rownames(A)[j],
                                    sep = ""))
            }
        }
    }

    dnf <- unique(dnf)

    return(dnf)

}
#' @noRd
#' @useDynLib mnem
#' @importFrom Rcpp sourceCpp
mytc <- function(x, u = NULL, v = NULL) {
    diag(x) <- 1
    if (is.null(u) | is.null(v)) {
        y <- matrix(transClose_W(x), nrow(x), ncol(x))
    } else {
        if (x[u,v] == 1) {
            y <- matrix(transClose_Ins(x, u, v), nrow(x), ncol(x))
        }
        if (x[u,v] == 0) {
            y <- matrix(transClose_Del(x, u, v), nrow(x), ncol(x))
        }
    }
    y[which(y != 0)] <- 1
    rownames(y) <- rownames(x)
    colnames(y) <- colnames(x)
    return(y)
}
#' @noRd
simulateDnf <- function(dnf, stimuli = NULL, inhibitors = NULL) {
    getStateDnf <- function(node, signalStates, graph, children = NULL) {
        graphCut <- graph[grep(paste("=", node, "$", sep = ""), graph)]
        if (length(graphCut) == 0) {
            signalStates[, node] <- 0
        } else {
            sop <- numeric(nrow(signalStates))
            children2 <- gsub("!", "", children)
            for (i in graphCut) {
                parents <- gsub("=.*$", "", unlist(strsplit(i, "\\+")))
                pob <- rep(1, nrow(signalStates))
                for (j in parents) {
                    j2 <- gsub("!", "", j)
                    if (sum(is.na(signalStates[, j2]) == TRUE) ==
                        length(signalStates[, j2])) {
                        if (j %in% j2) {
                            node2 <- node
                            add1 <- 0
                        } else {
                            node2 <- paste("!", node, sep = "")
                            add1 <- 1
                        }
                        if (j2 %in% children2) {
                            subGraph <- graph[
                                -grep(paste(".*=", node, "|.*",
                                            j2, ".*=.*", sep = ""), graph)]
                            signalStatesTmp <- getStateDnf(
                                node = j2,
                                signalStates = signalStates,
                                graph = subGraph,
                                children = NULL)
                            ifa <- children[
                                which(children2 %in% j2):length(children2)]
                            ifb <- (length(grep("!", ifa)) + add1)/2
                            ifc <- children[
                                which(children2 %in% j2):length(children2)]
                            if (ifb != ceiling((length(grep("!", ifc)) +
                                                add1)/2)) {
                            } else {
                            }
                            if (add1 == 0) {
                                pobMult <- signalStatesTmp[, j2]
                            } else {
                                pobMult <- add1 - signalStatesTmp[, j2]
                            }
                        } else {
                            signalStates <-
                                getStateDnf(node = j2,
                                            signalStates = signalStates,
                                            graph = graph,
                                            children = unique(c(children,
                                                                node2)))
                            if (add1 == 0) {
                                pobMult <- signalStates[, j2]
                            } else {
                                pobMult <- add1 - signalStates[, j2]
                            }
                        }
                        pob <- pob*pobMult
                    } else {
                        if (j %in% j2) {
                            add1 <- 0
                        } else {
                            add1 <- 1
                        }
                        if (add1 == 0) {
                            pobMult <- signalStates[, j2]
                        } else {
                            pobMult <- add1 - signalStates[, j2]
                        }
                        pob <- pob*pobMult
                    }
                    if (max(pob, na.rm = TRUE) == 0) { break() }
                }
                sop <- sop + pob
                if (min(sop, na.rm = TRUE) > 0) { break() }
            }
            sop[sop > 0] <- 1
            if (node %in% inhibitors) {
                sop <- sop*0
            }
            if (node %in% stimuli) {
                sop <- max(sop, 1)
            }
            signalStates[, node] <- sop
        }
        return(signalStates)
    }
    signals <-
        unique(gsub("!", "",
                    unlist(strsplit(
                        unlist(strsplit(dnf, "=")), "\\+"))))
    graph <- dnf
    signalStates <- matrix(NA, nrow = 1, ncol = length(signals))
    rownames(signalStates) <- paste(c("stimuli:", stimuli, "inhibitors:",
                                      inhibitors), collapse = " ")
    colnames(signalStates) <- signals
    signalStates[which(signals %in% stimuli)] <- 1
    for (k in signals) {
        if (is.na(signalStates[, k]) == TRUE) {
            signalStates <- getStateDnf(node = k, signalStates = signalStates,
                                        graph = graph, children = NULL)
        }
    }
    namestmp <- colnames(signalStates)
    signalStates <- as.vector(signalStates)
    names(signalStates) <- namestmp
    return(signalStates = signalStates)
}
#' @noRd
#' @importFrom graphics boxplot axis
myboxplot <- function(x, box = TRUE, dens = TRUE, scatter = "no",
                      polygon = TRUE, sd = 0.2, dcol = NULL,
                      scol = NULL, dlty = 1,
                      dlwd = 1, spch = 1, xaxt = "y", ...) {
    paras <- list(...)
    n <- ncol(x)
    if (box) {
        boxplot(x, xaxt = "n", ...)
    }
    if (dens) {
        if (is.null(dcol) & !is.null(paras$col)) {
            dcol <- paras$col
        } else {
            dcol <- rgb(0,0,0,0.5)
        }
        if (!box) {
            if (is.null(paras$ylim)) {
                yrange <- apply(x, 2, function(y) {
                    d <- density(y, na.rm = TRUE)
                    return(d$x)
                })
                plot(0, 0, xlim = c(0.5, n+0.5), ylim = c(min(yrange),
                                                          max(yrange)),
                     xaxt = "n", ...)
            } else {
                plot(0, 0, xlim = c(0.5, n+0.5), xaxt = "n", ...)
            }
        }
        if (length(dcol) == 1) { dcol <- rep(dcol, n) }
        for (i in seq_len(n)) {
            if (any(x[, i] != x[1, i])) {
                d <- density(x[, i], na.rm = TRUE)
                dy <- d$y
                dx <- d$x
                dy <- dy/max(dy)*0.5
                lines(c(dy+i,-dy+i), rep(dx, 2), lty = dlty, lwd = dlwd,
                      col = dcol[i])
                if (polygon) {
                    polygon(c(dy+i,-dy+i), rep(dx, 2), col = dcol[i])
                }
            }
        }
    }
    if ("random" %in% scatter) {
        if (is.null(scol) & !is.null(paras$col)) {
            scol <- paras$col
        } else {
            scol <- rgb(0,0,0,0.5)
        }
        lines(rep(seq_len(n), each = nrow(x))+rnorm(nrow(x)*n, 0, sd),
              as.vector(x), type = "p",
              pch = spch, col = rep(scol, each = nrow(x)))
    }
    if (xaxt != "n") {
        axis(1, seq_len(n), seq_len(n))
    }
}
#' @noRd
#' @importFrom graphics abline
addgrid <- function(x = c(0,1,0.1), y = c(0,1,0.1), lty = 2,
                    col = rgb(0.5,0.5,0.5,0.5)) {
    abline(h=seq(x[1], x[2], x[3]), col = col, lty = lty)
    abline(v=seq(y[1], y[2], y[3]), col = col, lty = lty)
}
