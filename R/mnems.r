#' @noRd
#' @export
annotAdj <- function(adj, data) {
    Sgenes <- sort(unique(colnames(data)))
    colnames(adj) <- rownames(adj) <- sort(Sgenes)
    return(adj)
}

#' plot mnem result
#' @param x mnem object
#' @param man logical. manual data penalty
#' @param degree different degree of complexity
#' @param logtype logarithm type of the data
#' @param pen penalty weight ofr the data
#' @param useF use F (see publication) as complexity
#' @param Fnorm normalize complexity of F
#' @author Martin Pirkl
#' @return penlized likelihood
#' @export
#' @examples
#' sim <- simData(Sgenes = 5, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' pen <- numeric(3)
#' result <- list()
#' for (k in 1:3) {
#' result[[k]] <- mnem(data, k = k, starts = 10)
#' pen[k] <- getIC(result[[k]])
#' }
#' print(pen)
getIC <- function(x, man = FALSE, degree = 4, logtype = 2, pen = 2, useF = FALSE, Fnorm = FALSE) {
    n <- ncol(x$data)
    if (useF) {
        for (i in 1:length(x$comp)) {
            tmp <- transitive.closure(x$comp[[i]]$phi, mat = TRUE)
            tmp2 <- matrix(0, nrow = nrow(tmp), ncol = length(x$comp[[i]]$theta))
            tmp3 <- x$comp[[i]]$theta
            tmp3[which(tmp3 > nrow(tmp2))] <- 0
            tmp2[cbind(tmp3, 1:ncol(tmp2))] <- 1
            tmp4 <- tmp%*%tmp2
            if (i == 1) {
                fpar <- tmp4
            } else {
                fpar <- fpar + tmp4
            }
        }
        if (Fnorm) {
            fpar[which(fpar > 0)] <- 1
        }
        fpar <- sum(fpar)
        fpar <- fpar + length(x$comp) - 1
    } else {
        if (degree != 5) {
            fpar <- 0   
            for (i in 1:length(x$comp)) {
                tmp <- transitive.reduction(x$comp[[i]]$phi)
                if (degree > 2) {
                    tmp <- transitive.closure(x$comp[[i]]$phi, mat = TRUE)
                }
                diag(tmp) <- 0
                fpar <- fpar + sum(tmp != 0)
            }
            if (degree > 1 & (degree != 3)) {
                fpar <- fpar + length(x$comp)*length(x$comp[[1]]$theta)
            }
        } else {
            fpar <- (length(x$comp[[1]]$phi)+length(x$comp[[1]]$theta))*length(x$comp)
        }
        if (degree > 0) {
            fpar <- fpar + length(x$comp) - 1
        }
    }
    LL <- max(x$ll)*log(logtype)
    if (man) {
        ic <- pen*fpar - 2*LL
    } else {
        ic <- log(n)*fpar - 2*LL
    }
    return(ic)
}
#' @noRd
#' @export
getOmega <- function(data) {
    
    Sgenes <- unique(unlist(strsplit(colnames(data), "_")))
    
    Omega <- matrix(0, length(Sgenes), ncol(data))
    
    for (i in seq_len(length(Sgenes))) {
        Omega[i, grep(Sgenes[i], colnames(data))] <- 1
    }
    
    rownames(Omega) <- Sgenes
    colnames(Omega) <- colnames(data)

    return(Omega)
}
#' @noRd
#' @export
initComps <- function(data, k=2, starts=1, verbose = FALSE, meanet = NULL) {
    n <- getSgeneN(data)
    nets <- list()
    for (i in seq_len(starts*k)) {
        tmp <- matrix(sample(c(0,1), replace = TRUE), n, n)
        tmp[lower.tri(tmp)] <- 0
        colnames(tmp) <- rownames(tmp) <- sample(1:n, n)
        tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
        nets[[i]] <- tmp
    }
    nets <- sortAdj(nets, list = TRUE)$res
    for (j in seq_len(starts)) {
        init[[j]] <- list()
        for (i in seq_len(k)) {
            do <- (i+starts*(i-1)+j-1-(i-1)*1)
            init[[j]][[i]] <- nets[[do]]
        }
    }
    return(init)
}
#' @noRd
#' @export
initps <- function(data, ks, k, starts = 3) {
    ## based on the clustering for each knock-down, we estimate non-random membership values for each start of the em:
    clusters <- list()
    for (i in seq_len(length(unique(colnames(data))))) {
        d <- dist(t(data[, which(colnames(data) %in% i)]))
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
        tmp <- matrix(0, k, ncol(data))
        if (ks[1] < k) {
            takes <- as.matrix(sample(1:ks[1], k, replace = TRUE), k, 1)
        } else {
            takes <- as.matrix(sample(1:ks[1], k), k, 1)
        }
        for (i in 2:length(unique(colnames(data)))) {
            if (ks[i] < k) {
                takes <- cbind(takes, sample(1:ks[i], k, replace = TRUE))
            } else {
                takes <- cbind(takes, sample(1:ks[i], k))
            }
        }
        for (i in 1:k) {
            for (j in 1:n) {
                tmp[i, which(colnames(data) == j)[which(clusters[[j]] == takes[i, j])]] <- 1
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
#' @export
modData <- function(D) {
    SgeneN <- getSgeneN(D)
    Sgenes <- naturalsort(unique(colnames(D)))
    if (!all(is.numeric(Sgenes))) {
        colnamesD <- numeric(ncol(D))
        for (i in 1:SgeneN) {
            colnamesD[which(colnames(D) %in% Sgenes[i])] <- i
        }
        colnames(D) <- as.numeric(colnamesD)
    }
    rownames(D) <- as.numeric(1:nrow(D))
    return(D)
}
#' @noRd
#' @export
learnk <- function(data, kmax = 10, verbose = FALSE) {
    ks <- numeric(length(unique(colnames(data))))
    lab <- list()
    for (i in naturalsort(as.numeric(unique(colnames(data))))) {
        if (verbose) {
            print(i)
        }
        if (sum(colnames(data) %in% i) <= 1) { k <- 1; next() }
        d <- dist(t(data[, which(colnames(data) %in% i)]))
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
                sil <- silhouette(cluster, d)
                silavgs[j] <- mean(sil[, 3])
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
    }
    k <- min(kmax, max(ks))
    return(list(ks = ks, k = k, lab = lab))
}
#' @noRd
#' @export
getLL <- function(x, logtype = 2, mw = NULL, data = NULL) {
    if (is.null(mw)) { mw = rep(1, nrow(x))/nrow(x) }
    if (any(is.infinite(logtype^apply(data, 2, function(x) return(sum(x[which(x>0)])))))) {
        Z <- getAffinity(x, logtype = logtype, mw = mw, data = data)
        l <- sum(apply(Z*(x + log(mw)/log(logtype)), 2, sum))
    } else {
        x <- logtype^x
        x <- x*mw
        l <- sum(log(apply(x, 2, sum))/log(logtype))
    }
    return(l)
}
#' @noRd
#' @export
getAffinity <- function(x, affinity = 0, norm = TRUE, logtype = 2, mw = NULL, data = matrix(0, 2, ncol(x))) {
    if (is.null(mw)) { mw <- rep(1, nrow(x))/nrow(x) }
    if (affinity == 1) {
        y <- logtype^x
        y <- apply(y, 2, function(x) {
            xnmax <- which(x < max(x))
            if (length(xnmax) > 0) {
                x[xnmax] <- 0
            }
            return(x)
        })
        if (is.null(dim(y))) {
            y <- matrix(y, nrow = 1)
        }
        if (norm) {
            y[which(y != 0)] <- 1
        }
    } else {
        if (nrow(x) > 1) {
            y <- x
            if (norm) {
                if (any(is.infinite(logtype^apply(data, 2, function(x) return(sum(x[which(x>0)])))))) {
                    y <- apply(y, 2, function(x) {
                        xmax <- max(x)
                        maxnum <- 2^1023
                        shrinkfac <- log(maxnum)/log(logtype)
                        x <- x - (xmax - shrinkfac)
                        return(x)
                    })
                }
                y <- logtype^y
                y <- y*mw
                y <- apply(y, 2, function(x) return(x/sum(x)))
            }
        } else {
            y <- matrix(1, nrow(x), ncol(x))
        }
    }
    y[which(is.na(y) == TRUE)] <- 0
    return(y)
}
#' @noRd
#' @export
estimateSubtopo <- function(data) {
    effectsums <- effectsds <- matrix(0, nrow(data), length(unique(colnames(data))))
    for (i in 1:length(unique(colnames(data)))) {
        if (length(grep(i, colnames(data))) > 1) {
            effectsds[, i] <- apply(data[, grep(i, colnames(data))], 1, sd)
            effectsums[, i] <- apply(data[, grep(i, colnames(data))], 1, sum)
        } else {
            effectsds[, i] <- 1
            effectsums[, i] <- data[, grep(i, colnames(data))]
        }
    }
    subtopoX <- as.numeric(apply(effectsums/effectsds, 1, which.max))
    subtopoX[which(is.na(subtopoX) == TRUE)] <- 1
    return(subtopoX)
}
#' @noRd
#' @export
getProbs <- function(probs, k, data, res, method = "llr", n, affinity = 0, converged = 10^-2, subtopoX = NULL, ratio = TRUE, logtype = 2) {
    if (is.null(subtopoX)) {
        subtopoX <- estimateSubtopo(data)
    }
    subtopoY <- bestsubtopoY <- subtopoX
    bestprobs <- probsold <- probs
    time0 <- TRUE
    count <- 0
    max_count <- 100
    ll0 <- 0
    stop <- FALSE
    mw <- apply(getAffinity(probsold, affinity = affinity, norm = TRUE, logtype = logtype, data = data), 1, sum)
    mw <- mw/sum(mw)
    if (any(is.na(mw))) { mw <- rep(1, k)/k }
    while((!stop | time0) & count < max_count) {
        llold <- max(ll0)
        time0 <- FALSE
        probsold <- probs
        subtopo0 <- matrix(0, k, nrow(data))
        subweights0 <- matrix(0, nrow(data), n+1) # account for null node
        postprobsold <- getAffinity(probsold, affinity = affinity, norm = TRUE, logtype = logtype, mw = mw, data = data)
        align <- list()
        for (i in 1:k) {
            n <- getSgeneN(data)
            dataF <- matrix(0, nrow(data), n)
            colnames(dataF) <- 1:n
            nozero <- which(postprobsold[i, ] != 0)
            if (length(nozero) != 0) {
                dataR <- cbind(data[, nozero, drop = FALSE], dataF)
                postprobsoldR <- c(postprobsold[i, nozero], rep(0, n))
            } else {
                dataR <- dataF
                postprobsoldR <- rep(0, n)
            }
            align[[i]] <- scoreAdj(dataR, res[[i]]$adj,
                                   method = method, ratio = ratio, weights = postprobsoldR)
            subtopo0[i, ] <- align[[i]]$subtopo
            subweights0 <- subweights0 + align[[i]]$subweights
        }
        subtopoMax <- apply(subweights0, 1, which.max)
        subtopoMax[which(subtopoMax > n)] <- subtopoMax[which(subtopoMax > n)] - n
        subtopo0 <- rbind(subtopoMax, subtopo0, subtopoX, subtopoY)
        ## go through all individual subtopolgies, the last subtopolgy and the max subtopology over all components
        probs0 <- list()
        ll0 <- numeric(nrow(subtopo0)+1)
        for (do in 1:(nrow(subtopo0)+1)) {
            probs0[[do]] <- probsold*0
            if (do > 1) {
                subtopo <- subtopo0[do-1, ]
            }
            for (i in 1:k) {
                if (do == 1) {
                    subtopo <- align[[i]]$subtopo
                }
                adj1 <- transitive.closure(res[[i]]$adj, mat = TRUE)
                adj1 <- cbind(adj1, "0" = 0)
                adj2 <- adj1[, subtopo]
                tmp <- llrScore(t(data), t(adj2), ratio = ratio)
                probs0[[do]][i, ] <- tmp[cbind(1:nrow(tmp), as.numeric(rownames(tmp)))]
            }
            ll0[do] <- getLL(probs0[[do]], logtype = logtype, mw = mw, data = data)
        }
        if (which.max(ll0) == 1) {
            sdo <- 1
        } else {
            sdo <- which.max(ll0) - 1
        }
        if (max(ll0) - llold > 0) {
            bestprobs <- probs0[[which.max(ll0)]]
            bestsubtopoY <- subtopo0[sdo, ]
        }
        probs <- probs0[[which.max(ll0)]]
        subtopoY <- subtopo0[sdo, ]
        if (max(ll0) - llold <= converged) {
            stop <- TRUE
        }
        mw <- apply(getAffinity(probs, affinity = affinity, norm = TRUE, logtype = logtype, mw = mw, data = data), 1, sum)
        mw <- mw/sum(mw)
        count <- count + 1
        ## if (verbose) { print(max(ll0)) }
    }
    ## if (count == max_count) { print(paste("no cell score convergence (", max_count, " iterations)", sep = "")) }
    return(list(probs = bestprobs, subtopoX = bestsubtopoY))
}

#' learn mixture of networks from a single-ceel knock-down experiment
#' @param D data with cells indexing the columns and features (E-genes) indexing the rows
#' @param inference inference method "em" for expectation maximization, "greedy" or "genetic"
#' @param search search method for single network inference "greedy", "exhaustive" or "modules"
#' @param start A list of n lists of k networks for n starts of the EM and k components
#' @param method "llr" for log ratios or foldchanges as input (see ratio)
#' @param parallel number of threads
#' @param reduce logical - reduce search space for exhaustive search
#' @param runs how many runs for greedy search
#' @param starts how many starts for em
#' @param type initialize with "random" probabilities or "cluster" the data
#' @param p initial probabilities
#' @param k number of components
#' @param kmax maximum number of components when k=NULL is inferred
#' @param verbose output
#' @param max_iter maximum iteration if likelihood does not converge
#' @param parallel2 if parallel=NULL, number of threads for component optimization
#' @param converged absolute distance for convergence
#' @param redSpace space for exhaustive search
#' @param affinity 0 is default, 1 is for hard clustering
#' @param evolution logical, penalty on difference of components
#' @param subtopoX hard prior on theta
#' @param ratio logical, if true data is log ratios, if false foldchanges
#' @param logtype logarithm of the data
#' @param initnets initial networks in a list
#' @param popSize population size
#' @param stallMax maximum number of stall generations till convergence
#' @param elitism number of networks to keep during evolution
#' @param maxGens maximum numbero f generations
#' @param domean average the data (speed imporvment)
#' @param modulesize max number of S-genes per module in module search
#' @param compress compress networks after search
#' @author Martin Pirkl
#' @return optimized network for data fit
#' @export
#' @import
#' epiNEM
#' cluster
#' nem
#' graph
#' Rgraphviz
#' naturalsort
#' flexclust
#' snowfall
#' grid
#' lattice
#' modeltools
#' stats4
#' @examples
#' sim <- simData(Sgenes = 5, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' result <- mnem(data, k = 2, starts = 10)
#' plot(result)
mnem <- function(D, inference = "em", search = "modules", start = NULL, method = "llr",
                 parallel = NULL, reduce = FALSE, runs = 1, starts = 3, type = "random",
                 p = NULL, k = NULL, kmax = 10, verbose = FALSE,
                 max_iter = 100, parallel2 = NULL, converged = 10^-1,
                 redSpace = NULL, affinity = 0, evolution = FALSE,
                 subtopoX = NULL, ratio = TRUE, logtype = 2, initnets = FALSE,
                 popSize = 10, stallMax = 2, elitism = NULL, maxGens = Inf,
                 domean = TRUE, modulesize = 5, compress = FALSE, increase = TRUE) {
    if (reduce & search %in% "exhaustive" & is.null(redSpace)) {
        redSpace <- mynem(data[, -which(duplicated(colnames(data)) == TRUE)], search = "exhaustive", reduce = TRUE, verbose = verbose, parallel = c(parallel, parallel2), subtopo = subtopoX, ratio = ratio, domean = FALSE, modulesize = modulesize)$redSpace
    }
    if (!is.null(parallel)) { if (parallel == 1) { parallel <- NULL } }
    D.backup <- D
    D <- modData(D)
    Sgenes <- getSgenes(D)
    data <- D
    D <- NULL
    ## learn k:
    if (is.null(k) & is.null(start) & is.null(p)) {
        tmp <- learnk(data, kmax = kmax)
        k <- tmp$k
        print(paste("components detected: ", k, sep = ""))
        ks <- tmp$ks
    } else {
        ks <- rep(k, length(Sgenes))
    }
    if (is.null(subtopoX)) {
        subtopoX <- estimateSubtopo(data)
    }
    if (!is.null(k)) {
        if (initnets & is.null(start)) {
            meanet <- mynem(D = data, search = search, start = start, method = method,
                               parallel = parallel2, reduce = reduce, runs = runs,
                               verbose = verbose, redSpace = redSpace, ratio = ratio, domean = domean, modulesize = modulesize)
            init <- initComps(data, k, starts, verbose, meanet)
            probscl <- 0
        } else {
            if (type %in% "cluster") {
                probscl <- initps(data, ks, k, starts = starts)
                ## starts <- length(probscl)
            } else {
                probscl <- 0
            }
        }
    } else {
        probscl <- 0
    }
    ## if we start with specific membership values or networks we do not do several runs:
    if (!is.null(parallel)) { parallel2 <- NULL }
    if (!is.null(start) | !is.null(p)) { starts <- length(start); k <- max(c(nrow(p), length(start[[1]]))) }
    init <- start
    res1 <- NULL
    mw <- rep(1, k)/k
    if (starts <= 1) { parallel2 <- parallel; parallel <- NULL }
    if (!is.null(parallel)) { parallel2 <- NULL }
    ## learn "mixture" of k = 1 or k > 1:
    if (k == 1) {
        if (!is.null(init)) {
            start <- start[[1]][[1]]
        }
        if (!is.null(parallel) & is.null(parallel2)) { parallel2 <- parallel }
        print("no evidence for sub populations or k set to 1")
        limits <- list()
        limits[[1]] <- list()
        limits[[1]]$res <- list()
        limits[[1]]$res[[1]] <- mynem(D = data, search = search, start = start, method = method,
                               parallel = parallel2, reduce = reduce, runs = runs,
                               verbose = verbose, redSpace = redSpace, ratio = ratio, domean = domean, modulesize = modulesize)
        limits[[1]]$res[[1]]$D <- NULL
        res <- list()
        res[[1]] <- list()
        res[[1]]$adj <- limits[[1]]$res[[1]]$adj
        n <- ncol(res[[1]]$adj)
        mw <- 1
        probs0 <- list()
        probs0$probs <- matrix(0, k, ncol(data))
        for (i in 1:2) { # i == 1 seems to win most of the time
            probs <- matrix(i - 1, k, ncol(data))
            probs <- getProbs(probs, k, data, res, method, n, affinity, converged, subtopoX, ratio)
            if (getLL(probs$probs, logtype = logtype, mw = mw, data = data) > getLL(probs0$probs, logtype = logtype, mw = mw, data = data)) {
                probs0 <- probs
            }
        }
        subtopoX <- probs0$subtopoX
        probs <- probs0$probs
        limits[[1]]$ll <- getLL(probs, logtype = logtype, mw = mw, data = data)
        limits[[1]]$probs <- probs
    } else {
        if (inference %in% "em") {
            if (!is.null(parallel)) {
                transitive.closure <- nem::transitive.closure
                transitive.reduction <- nem::transitive.reduction
                sfInit(parallel = TRUE, cpus = parallel)
                sfExport("modules", "mw", "ratio", "getSgeneN", "modData", "sortAdj", "calcEvopen", "evolution", "transitive.reduction", "getSgenes", "estimateSubtopo", "getLL", "getAffinity", "get.insertions", "get.reversions", "get.deletions", "D", "mynem", "scoreAdj", "transitive.closure", "max_iter", "verbose", "llrScore", "search", "redSpace", "affinity", "getProbs", "probscl", "method")#, "start", "better", "traClo", "method", "scoreAdj", "weights", "transitive.closure")
            }
            do_inits <- function(s) {
                ## initial with random membership values if not input is given:
                if (!is.null(init)) {
                    res1 <- list()
                    for (i in 1:length(init[[s]])) {
                        res1[[i]] <- list()
                        res1[[i]][["adj"]] <- init[[s]][[i]]
                    }
                    if (is.list(subtopoX)) {
                        subtopoX <- estimateSubtopo(data)
                    }
                    k <- length(res1)
                    n <- ncol(res1[[1]]$adj)
                    if (is.null(p)) {
                        probs0 <- list()
                        probs0$probs <- matrix(0, k, ncol(data))
                        for (i in 1:2) {
                            probs <- matrix(i - 1, k, ncol(data))
                            probs <- getProbs(probs, k, data, res1, method, n, affinity, converged, subtopoX, ratio)
                            if (getLL(probs$probs, logtype = logtype, mw = mw, data = data1) > getLL(probs0$probs, logtype = logtype, mw = mw, data = data1)) {
                                probs0 <- probs
                            }
                        }
                        subtopoX <- probs0$subtopoX
                        probs <- probs0$probs
                    } else {
                        probs <- p
                    }
                }
                if (is.null(p)) {
                    if (length(probscl) >= s & type %in% "cluster") {
                        probs <- probscl[[s]]
                    } else {
                        probs <- matrix(log2(sample(c(0,1), k*ncol(data), replace = TRUE, prob = c(0.9, 0.1))), k, ncol(data))
                    }
                    mw <- apply(getAffinity(probs, affinity = affinity, norm = TRUE, logtype = logtype, mw = mw, data = data), 1, sum)
                    mw <- mw/sum(mw)
                } else {
                    k <- nrow(p)
                    probs <- p
                }
                mw <- apply(getAffinity(probs, affinity = affinity, norm = TRUE, logtype = logtype, mw = mw, data = data), 1, sum)
                mw <- mw/sum(mw)
                if (any(is.na(mw))) { mw <- rep(1, k)/k }
                if (verbose) {
                    print(paste("start...", s))
                }
                limits <- list()
                ll <- getLL(probs, logtype = logtype, mw = mw, data = data)
                llold <- bestll <- -Inf
                lls <- NULL
                count <- 0
                time0 <- TRUE
                probsold <- probs
                while ((((ll - llold > converged & increase) | (abs(ll - llold) > converged & !increase) & count < max_iter)) | time0) {
                    if (!time0) {
                        if (ll - bestll > 0) {
                            bestll <- ll; bestres <- res1; bestprobs <- probs
                        }
                        llold <- ll
                        probsold <- probs
                    }
                    time0 <- FALSE
                    res <- list()
                    postprobs <- getAffinity(probs, affinity = affinity, norm = TRUE, logtype = logtype, mw = mw, data = data)
                    edgechange <- 0
                    thetachange <- 0
                    for (i in 1:k) {
                        if (!is.null(res1) & !("modules" %in% search)) {
                            start0 <- res1[[i]]$adj
                        } else {
                            start0 <- NULL
                        }
                        n <- getSgeneN(data)
                        dataF <- matrix(0, nrow(data), n)
                        colnames(dataF) <- 1:n
                        nozero <- which(postprobs[i, ] != 0)
                        if (length(nozero) != 0) {
                            if (length(nozero) == ncol(data)) {
                                dataR <- data
                                postprobsR <- postprobs[i, ]
                            } else {
                                dataR <- cbind(data[, nozero, drop = FALSE], dataF)
                                postprobsR <- c(postprobs[i, nozero], rep(0, n))
                            }
                        } else {
                            dataR <- dataF
                            postprobsR <- rep(0, n)
                        }
                        if (is.null(start0)) {
                            res[[i]] <- mynem(D = dataR, weights = postprobsR,
                                              search = search, start = start0,
                                              method = method,
                                              parallel = parallel2, reduce = reduce,
                                              runs = runs, verbose = verbose,
                                              redSpace = redSpace, ratio = ratio, domean = domean, modulesize = modulesize)#, subtopo = subtopoX)
                        } else {
                            test01 <- list()
                            test01scores <- numeric(3)
                            for (j in 1:3) {
                                if (j == 3) {
                                    start1 <- start0
                                } else {
                                    start1 <- start0*0 + (j - 1)
                                }
                                test01[[j]] <- mynem(D = dataR, weights = postprobsR,
                                              search = search, start = start1,
                                              method = method,
                                              parallel = parallel2, reduce = reduce,
                                              runs = runs, verbose = verbose,
                                              redSpace = redSpace, ratio = ratio, domean = domean, modulesize = modulesize)#, subtopo = subtopoX)
                                test01scores[j] <- max(test01[[j]]$scores)
                            }
                            res[[i]] <- test01[[which.max(test01scores)]]
                        }
                        edgechange <- edgechange + sum(abs(res[[i]]$adj - res1[[i]]$adj))
                        thetachange <- thetachange + sum(res[[i]]$subtopo != res1[[i]]$subtopo)
                        res[[i]]$D <- NULL
                        res[[i]]$subweights <- NULL
                        
                    }
                    evopen <- 0
                    if (evolution) {
                        ## order nems: (do i have to do this once or every damn time?)
                        res <- sortAdj(res)$res
                        evopen <- calcEvopen(res)
                    }
                    res1 <- res
                    ## e:
                    n <- ncol(res[[1]]$adj)
                    probs0 <- list()
                    probs0$probs <- matrix(0, k, ncol(data))
                    for (i in 1:3) { # i == 1 seems to win most of the time
                        if (i == 3) {
                            probs <- probsold
                        } else {
                            probs <- matrix(i - 1, k, ncol(data))
                        }
                        probs <- getProbs(probs, k, data, res, method, n, affinity, converged, subtopoX, ratio)
                        if (getLL(probs$probs, logtype = logtype, mw = mw, data = data) > getLL(probs0$probs, logtype = logtype, mw = mw, data = data)) {
                            probs0 <- probs
                        }
                    }
                    subtopoX <- probs0$subtopoX
                    probs <- probs0$probs
                    modelsize <- n*n*k
                    datasize <- nrow(data)*ncol(data)*k
                    ll <- getLL(probs, logtype = logtype, mw = mw, data = data) + evopen*datasize*(modelsize^-1)
                    mw <- apply(getAffinity(probs, affinity = affinity, norm = TRUE, logtype = logtype, mw = mw, data = data), 1, sum)
                    mw <- mw/sum(mw)
                    if(verbose) {
                        print(paste("ll: ", ll, sep = ""))
                        print(paste("changes in phi(s): ", edgechange, sep = ""))
                        print(paste("changes in theta(s): ", thetachange, sep = ""))
                        if (evolution) {
                            print(paste("evolution penalty: ", evopen, sep = ""))
                        }
                    }
                    lls <- c(lls, ll)
                    count <- count + 1
                }
                if (ll - bestll > 0) {
                    bestll <- ll; bestres <- res1; bestprobs <- probs
                }
                if (abs(ll - llold) > converged | llold > ll) { # llold > ll probably hardly happens but when it does it prob goes into a cycle for 100 iterations. not necesssarily. e.g. with modulesearch I do not give the last optimum as a starting network and cannot guarantee ll increase. improve modulesearch to add a starting network....
                    if (verbose & (increase | max_iter <= count)) {
                        print("no convergence")
                    }
                }
                ## ## one last M-step? I think this is actually wrong! If I do one last M-step I have to recalculate the weights, otherwise we loose consistency.
                limits <- list()
                limits$probs <- bestprobs # I have to set the old ones because if the max_iter hits, we lose the best otherwise!
                limits$res <- bestres
                limits$ll <- lls
                limits$k <- k
                limits$subtopo <- subtopoX
                return(limits)
            }
            if (!is.null(parallel)) {
                limits <- sfLapply(1:starts, do_inits)
            } else {
                limits <- lapply(1:starts, do_inits)
            }
            if (!is.null(parallel)) {
                sfStop()
            }
        }
    }
    best <- limits[[1]]
    if (length(limits) > 1) {
        for (i in 2:length(limits)) {
            if (max(best$ll) <= max(limits[[i]]$ll)) {
                best <- limits[[i]]
            }
        }
    }
    if (compress) {
        unique <- list()
        count <- 1
        added <- 1
        for (i in 1:length(best$res)) {
            if (i == 1) {
                unique[[1]] <- best$res[[1]]
                next()
            }
            count <- count + 1
            unique[[count]] <- best$res[[i]]
            added <- c(added, i)
            for (j in 1:(length(unique)-1)) {
                if (all(unique[[count]]$adj - unique[[j]]$adj  == 0)) {
                    unique[[count]] <- NULL
                    count <- count - 1
                    added <- added[-length(added)]
                    break()
                }
            }
        }
        ## I try to merge components based on the fact, that some components can be completely included in the others and are obsolete:
        dead <- NULL
        for (i in 1:length(unique)) {
            dups <- NULL
            a <- transitive.closure(unique[[i]]$adj, mat = TRUE)
            for (j in 1:length(unique)) {
                if (i %in% j | j %in% dead) { next() }
                b <- transitive.closure(unique[[j]]$adj, mat = TRUE)
                dups <- c(dups, which(apply(abs(a - b), 1, sum) == 0))
            }
            if (all(1:nrow(unique[[i]]$adj) %in% dups)) {
                dead <- c(dead, i)
            }
        }
        if (is.null(dead)) {
            ulength <- 1:length(unique)
        } else {
            added <- added[-dead]
            ulength <- (1:length(unique))[-dead]
        }
        count <- 0
        unique2 <- list()
        for (i in ulength) {
            count <- count + 1
            unique2[[count]] <- unique[[i]]    }
        unique <- unique2
        probs <- best$probs[added, , drop = FALSE]
        colnames(probs) <- colnames(D.backup)
        postprobs <- getAffinity(probs, affinity = affinity, norm = TRUE, logtype = logtype, mw = mw, data = data)
        if (!is.null(dim(postprobs))) {
            lambda <- apply(postprobs, 1, sum)
            lambda <- lambda/sum(lambda)
        } else {
            lambda <- 1
        }
        comp <- list()
        for (i in 1:length(unique)) {
            comp[[i]] <- list()
            comp[[i]]$phi <- unique[[i]]$adj
            comp[[i]]$theta <- unique[[i]]$subtopo
        }
    } else {
        probs <- best$probs[, , drop = FALSE]
        colnames(probs) <- colnames(D.backup)
        postprobs <- getAffinity(probs, affinity = affinity, norm = TRUE, logtype = logtype, mw = mw, data = data)
        if (!is.null(dim(postprobs))) {
            lambda <- apply(postprobs, 1, sum)
            lambda <- lambda/sum(lambda)
        } else {
            lambda <- 1
        }
        comp <- list()
        for (i in 1:length(best$res)) {
            comp[[i]] <- list()
            comp[[i]]$phi <- best$res[[i]]$adj
            comp[[i]]$theta <- best$res[[i]]$subtopo
        }
    }
    res <- list(limits = limits, comp = comp, data = D.backup, mw = lambda, probs = probs, ll = getLL(probs, logtype = logtype, mw = lambda, data = data))
    class(res) <- "mnem"
    return(res)
}
#' @noRd
#' @export
sortAdj <- function(res, list = FALSE) {
    resmat <- NULL
    for (i in 1:length(res)) {
        if (list) {
            resmat <- rbind(resmat, as.vector(transitive.closure(res[[i]], mat = TRUE))) # trans close?
        } else {
            resmat <- rbind(resmat, as.vector(transitive.closure(res[[i]]$adj, mat = TRUE))) # trans close?
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
    for (i in 1:length(res)) {
        res2[[i]] <- res[[resorder[i]]]
    }
    return(list(res = res2, order = resorder))
}
#' @noRd
#' @export
calcEvopen <- function(res) {
    evopen <- 0
    for (i in 1:(length(res)-1)) {
        evopen <- evopen + sum(abs(res[[i]]$adj - res[[(i+1)]]$adj))/length(res[[i]]$adj)
    }
    evopen <- -evopen#/(k-1)
    return(evopen)
}
#' @noRd
#' @export
llrScore <- function(data, adj, weights = NULL, ratio = TRUE) {
    if (is.null(weights)) {
        weights <- rep(1, ncol(data))
    }
    if (ratio) {
        score <- data%*%(adj*weights)
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
#' @export
modAdj <- function(adj, D) {
    Sgenes <- naturalsort(unique(colnames(D)))
    SgeneN <- getSgeneN(D)
    for (i in 1:SgeneN) {
        colnames(adj) <- rownames(adj) <- gsub(i, Sgenes[i], colnames(adj))
    }
    return(adj)
}
#' @noRd
#' @export
bootstrap <- function(x) {
    ## bootstrap on the components to get frequencies 
}

#' plot mnem result
#' @param x mnem object
#' @param oma outer margin
#' @param main main text
#' @param anno annotate cells
#' @param cexAnno text size
#' @param scale scale cells
#' @param global global clustering
#' @param egenes show egene attachments
#' @param sep seperate clusters
#' @param tsne use tsne
#' @param affinity logical use hard clustering if true
#' @param logtype log type of the data
#' @param cells show ccell attachments
#' @param pch cell symbol
#' @param legend show legend
#' @param showdata show data if true
#' @param bestCell show probability of best fitting cell
#' @param showprobs show probabilities
#' @param shownull show null node
#' @param ratio use log ratios (if true) or foldchanges
#' @param method "llr" for ratios
#' @param showweights shwo weights
#' @param ... additional parameters
#' @author Martin Pirkl
#' @return visualization of mnem result Rgraphviz
#' @export
#' @import
#' epiNEM
#' cluster
#' nem
#' graph
#' Rgraphviz
#' tsne
#' @examples
#' sim <- simData(Sgenes = 5, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' result <- mnem(data, k = 2, starts = 10)
#' plot(result)
plot.mnem <- function(x, oma = c(3,1,1,3), main = "M&NEM", anno = TRUE, cexAnno = 1, scale = NULL, global = TRUE, egenes = TRUE, sep = FALSE, tsne = FALSE, affinity = 0, logtype = 2, cells = TRUE, pch = ".", legend = FALSE, showdata = FALSE, bestCell = TRUE, showprobs = FALSE, shownull = TRUE, ratio = TRUE, method = "llr", showweights = TRUE, ...) {

    x2 <- x

    data <- x$data
    
    laymat <- rbind(1:(length(x$comp)+1), c(length(x$comp)+2, rep(length(x$comp)+3, length(x$comp))))

    if (legend & !showdata) {
        laymat <- matrix(1:(length(x$comp)+1), nrow = 1)
    }
    if (!legend & !showdata) {
        laymat <- matrix(1:(length(x$comp)), nrow = 1)
    }
    if (!legend & showdata) {
        laymat <- rbind(1:(length(x$comp)), c(length(x$comp)+1, rep(length(x$comp)+2, length(x$comp)-1)))
    }
    layout(laymat)
    
    par(oma=oma)
    if (legend) {
    plotDnf(c("Sgenes=Egenes", "Egenes=Cells", "Cells=Fit"), edgecol = rep("transparent", 3),
            nodeshape = list(Sgenes = "circle", Egenes = "box", Cells = "diamond", Fit = "circle"),
            nodelabel = list(Sgenes = "signaling\ngenes", Egenes = "effect\nreporters", Cells = "single\ncells", Fit = "highest\nresponsibility"),
            nodeheight = list(Sgenes = 2, Egenes = 0.5, Cells = 0.5, Fit = 0.5), nodewidth = list(Sgenes = 2, Egenes = 0.5, Cells = 0.5, Fit = 0.5),
            layout = "circo")
    }
    full <- x$comp[[1]]$phi
    mixnorm <- getAffinity(x$probs, affinity = affinity, norm = TRUE, logtype = logtype, mw = x$mw, data = data)
    mixnorm <- apply(mixnorm, 2, function(x) {
        xmax <- max(x)
        x[which(x != xmax)] <- 0
        x[which(x == xmax)] <- 1
        return(x)
    })
    if (is.null(dim(mixnorm))) {
        mixnorm <- matrix(mixnorm, 1, length(mixnorm))
    }
    realpct <- rowSums(mixnorm)
    realpct <- round(realpct/ncol(mixnorm), 3)*100
    unipct <- mixnorm
    unipct[, which(colSums(mixnorm) > 1)] <- 0
    unipct <- round(rowSums(unipct)/ncol(mixnorm), 3)*100
    Sgenes <- naturalsort(unique(colnames(x$data)))
    SgeneN <- length(Sgenes)
    for (i in 1:length(x$comp)) {

        shared <- unique(colnames(mixnorm)[which(apply(mixnorm, 2, function(x) return(sum(x != 0))) != 1 & mixnorm[i, ] != 0)])
        net <- x$comp[[i]]$phi
        for (j in 1:SgeneN) {
            colnames(net)[which(colnames(net) %in% j)] <- rownames(net)[which(rownames(net) %in% j)] <- Sgenes[j]
            shared[which(shared %in% j)] <- Sgenes[j]
        }
        graph <- adj2dnf(net)
        pathedges <- length(graph)
        if (egenes) {
            enodes <- list()
            enodeshape <- list()
            enodeheight <- list()
            enodewidth <- list()
            probs <- x$probs
            if (is.null(x$comp[[i]]$theta)) {
                weights <- getAffinity(x$probs, affinity = affinity, norm = TRUE, logtype = logtype, mw = x$mw, data = data)
                subtopo <- scoreAdj(modData(x$data), x$comp[[i]]$phi, method = method, weights = weights[i, ],
                                    ratio = ratio, ...)$subtopo
            } else {
                subtopo <- x$comp[[i]]$theta
            }
            for (j in 1:SgeneN) {
                tmpN <- paste("_9247E", j, sep = "_")
                enodes[[tmpN]] <- sum(subtopo == j)
                enodeshape[[tmpN]] <- "box"
                enodewidth[[tmpN]] <- 0.5
                enodeheight[[tmpN]] <- 0.5
                if (tmpN != 0) {
                    graph <- c(graph, paste(Sgenes[j], tmpN, sep = "="))
                }
            }
        } else {
            enodes <- NULL
            enodeshape <- NULL
            enodeheight <- NULL
            enodewidth <- NULL
        }
        if (cells) {
            datanorm <- modData(x$data)
            pnorm <- getAffinity(x$probs, affinity = affinity, norm = TRUE, logtype = logtype, mw = x$mw, data = data)
            pnorm <- apply(pnorm, 2, function(x) {
                xmax <- max(x)
                x[which(x != xmax)] <- 0
                x[which(x == xmax)] <- 1
                return(x)
            })
            if (is.null(dim(pnorm))) {
                pnorm <- matrix(pnorm, 1, length(pnorm))
            }
            cnodes <- list()
            cnodeshape <- list()
            cnodeheight <- list()
            cnodewidth <- list()
            for (j in 1:SgeneN) {
                tmpN <- paste("__9247C", j, sep = "_")
                cnodes[[tmpN]] <- sum(colnames(datanorm)[which(pnorm[i, ] == 1)] == j)
                cnodeshape[[tmpN]] <- "diamond"
                cnodewidth[[tmpN]] <- 0.5
                cnodeheight[[tmpN]] <- 0.5
                if (tmpN != 0) {
                    graph <- c(graph, paste(Sgenes[j], tmpN, sep = "="))
                }
            }
        } else {
            cnodes <- NULL
            cnodeshape <- NULL
            cnodeheight <- NULL
            cnodewidth <- NULL
        }
        if (bestCell) {
            bnodes <- bnodeshape <- bnodeheight <- bnodewidth <- list()
            if (nrow(x$probs) > 1) {
                gam <- getAffinity(x$probs, affinity = affinity, norm = TRUE, logtype = logtype, mw = x$mw, data = data)
            } else {
                gam <- (logtype^x$probs)*x$mw
                gam <- gam/gam
            }
            for (bnode in naturalsort(unique(colnames(gam)))) {
                tmpN <- paste0("_9247bnode", bnode)
                graph <- c(graph, paste0(bnode, "=_9247bnode", bnode))
                bnodes[[tmpN]] <- paste0(round(max(gam[i, which(colnames(gam) %in% bnode)]), 2)*100, "%")
                bnodeheight[[tmpN]] <- 0.5
                bnodewidth[[tmpN]] <- 0.5
                bnodeshape[[tmpN]] <- "circle"
            }
        } else {
            bnodes <- bnodeshape <- bnodeheight <- bnodewidth <- NULL
        }
        edgecol <- c(rep("black", pathedges), rep("grey", length(graph) - pathedges))
        if (egenes) {
            enodes[["_9247Null"]] <- "NULL"
            enodeshape[["_9247Null"]] <- "circle"
            enodewidth[["_9247Null"]] <- 1
            enodeheight[["_9247Null"]] <- 1
            nulltargets <- sum(subtopo == (SgeneN+1))
            enodes[["_9247Nulltargets"]] <- nulltargets
            enodeshape[["_9247Nulltargets"]] <- "box"
            enodewidth[["_9247Nulltargets"]] <- 0.5
            enodeheight[["_9247Nulltargets"]] <- 0.5
            graph <- c(graph, "_9247Null=_9247Nulltargets")
        }
        edgecol <- c(rep("black", pathedges), rep("grey", length(graph) - pathedges))
        if (showweights) {
            mainweights = paste("Cells: ", realpct[i], "% (unique: ", unipct[i], "%)\n
Mixture weight: ", round(x$mw[i], 3)*100, "%", sep = "")
        } else {
            mainweights <- NULL
        }
        plotDnf(graph, main = mainweights, bordercol = i+1, width = 1, connected = FALSE,
# signals = shared, inhibitors = Sgenes,
nodelabel = c(cnodes, enodes, bnodes),
nodeshape = c(cnodeshape, enodeshape, bnodeshape),
nodewidth = c(cnodewidth, enodewidth, bnodewidth),
nodeheight = c(cnodeheight, enodeheight, bnodeheight),
edgecol = edgecol)
        full <- net + full
    }

    if (showdata) {
        full[which(full > 1)] <- 1
        if (length(x$comp) > 1) {
            plot.adj(full)
        }
        if (nrow(pnorm) > 1) {
            pnorm <- apply(logtype^x$probs, 2, function(x) return(x/sum(x)))
            if (is.null(dim(pnorm))) {
                pnorm <- matrix(pnorm, 1, length(pnorm))
            }
            pcols <- rep(1, ncol(pnorm))
            pcols[which(apply(mixnorm, 2, max) == 1)] <- apply(mixnorm[, which(apply(mixnorm, 2, max) == 1)]*((matrix(1:nrow(mixnorm), nrow(mixnorm), ncol(mixnorm)))[, which(apply(mixnorm, 2, max) == 1), drop = FALSE]+1), 2, max)
            pcols[which(apply(mixnorm, 2, function(x) return(sum(x != 0))) != 1)] <- 1
            if (tsne) {
                pres <- list()
                pres$rotation <- tsne(t(pnorm))
            } else {
                pres <- prcomp(pnorm)
            }
            if (nrow(x$probs) == 2) {
                pres <- list()
                pres$rotation <- t(pnorm)
            }
            for (i in 1:SgeneN) {
                rownames(pres$rotation)[which(rownames(pres$rotation) %in% i)] <- Sgenes[i]
            }
            jittered <- pres$rotation[, 1:2]
            if (!sep) {
                global <- TRUE
            }
            if (global) {
                if (tsne) {
                    prtmp <- list()
                    prtmp$rotation <- tsne(t(x$data))
                } else {
                    prtmp <- prcomp(x$data)
                }
                jittered <- jittered*0
            }
            if (is.null(scale) ) {
                if (global) {
                    scale <- 1
                } else {
                    scale <- (max(jittered) - min(jittered))*0.5
                }
            }
            for (i in naturalsort(unique(pcols))) {
                maxind0 <- which(pcols == i)
                if (!global) {
                    if (tsne) {
                        prtmp <- list()
                        prtmp$rotation <- tsne(t(x$data[, which(pcols == i)]))
                    } else {
                        prtmp <- prcomp(x$data[, which(pcols == i)])
                    }
                    maxind <- 1:nrow(prtmp$rotation)
                } else {
                    maxind <- maxind0
                }
                if (sep) {
                    jittered[which(pcols == i), 1] <- (prtmp$rotation[maxind, 1] - mean(prtmp$rotation[maxind, 1]))*scale + jittered[maxind0, 1]
                    jittered[which(pcols == i), 2] <- (prtmp$rotation[maxind, 2] - mean(prtmp$rotation[maxind, 2]))*scale + jittered[maxind0, 2]
                } else {
                    if (!global) {
                        jittered[which(pcols == i), 1] <- (prtmp$rotation[maxind, 1] - mean(prtmp$rotation[maxind, 1]))
                        jittered[which(pcols == i), 2] <- (prtmp$rotation[maxind, 2] - mean(prtmp$rotation[maxind, 2]))
                    } else {
                        jittered[which(pcols == i), 1] <- prtmp$rotation[maxind, 1]
                        jittered[which(pcols == i), 2] <- prtmp$rotation[maxind, 2]
                    }
                }
            }
            if (anno) {
                plot(jittered, col = pcols, pch = "", main = main)
                text(jittered[, 1],
                     jittered[, 2],
                     labels = rownames(jittered),
                     srt = 45, pos = 1,
                     offset = 0, cex = cexAnno, col = pcols)
            } else {
                plot(jittered, col = pcols, main = main, pch = pch)
            }
            unique <- unique(pres$rotation[which(pcols != 1), 1:2])
            if (all(dim(unique) != 0) & !global) {
                for (i in 1:nrow(unique)) {
                    for (j in 1:nrow(unique)) {
                        if (i == j) { next() }
                        lines(unique[c(i,j), ], lty = 3, col = "grey")
                    }
                }
            }

        } else {

            prtmp <- prcomp(x$data)

            plot(prtmp$rotation[, 1:2], pch = ".")

        }
    }
}
