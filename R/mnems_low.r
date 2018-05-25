nemEst <- function(data, maxiter = 100, start = "full",
                   sumf = mean, alpha = 1, cut = 0,
                   kernel = "cosim", monoton = FALSE,
                   useCut = TRUE, useF = TRUE, ...) { # kernels can be cosim or cor # fun can be add, mult or sign
    if (sum(duplicated(colnames(data)) == TRUE) > 0) {
        data2 <- data[, -which(duplicated(colnames(data)) == TRUE)]
        for (j in unique(colnames(data))) {
            data2[, j] <- apply(data[, which(colnames(data) %in% j), drop = FALSE], 1, sumf)
        }
    } else {
        data2 <- data
    }
    R <- data2[, naturalsort(colnames(data2))] 
    n <- ncol(R)
    if (kernel %in% "cosim") {
        R2 <- t(R)%*%R
    }
    if (kernel %in% "cor") {
        R2 <- cor(R)
    }
    if (!(kernel %in% c("cosim", "cor"))) { stop("kernel neither set to 'cosim' nor 'cor'.") }
    if (alpha < 1) {
        C <- cor(R)
        C <- C2 <- solve(C, ...)
        for (r in 1:nrow(C)) {
            C[r, ] <- C[r, ]/(C2[r, r]^0.5)
        }
        for (c in 1:ncol(C)) {
            C[, c] <- C[, c]/(C2[c, c]^0.5)
        }
        diag(C) <- 1
        Cz <- apply(C, c(1,2), function(x) return(0.5*log((1+x)/(1-x)))) # conditional independence test with fisher-transform
        diag(Cz) <- 0
        Cz <- pnorm(((nrow(R) - n - 2 - 3)^0.5)*Cz)
        idx <- which(Cz >= alpha)
        Cp <- Cz
        Cp[-idx] <- 1
        Cp[idx] <- 0
    } else {
        Cp <- 1
        Cz <- 0
    }
    phibest <- phi <- matrix(0, n, n)
    rownames(phi) <- colnames(phi) <- colnames(R)
    E0 <- apply(R, 2, sum)
    phi <- phi[order(E0, decreasing = TRUE), order(E0, decreasing = TRUE)]
    phi[upper.tri(phi)] <- 1
    phi <- phi[naturalsort(rownames(phi)), naturalsort(colnames(phi))]
    phi <- transitive.closure(phi, mat = TRUE)
    E <- phi
    E <- E*Cp
    if (any(start %in% "full")) {
        phi <- phi
    } else if (any(start %in% "rand")) {
        phi <- phi*0
        diag(phi) <- 1
        phi[1:length(phi)] <- sample(c(0,1), length(phi), replace = TRUE)
    } else if (any(start %in% "null")) {
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
        P <- R%*%cbind(phi, 0)
        P[, grep("_", colnames(phi))] <- min(P)
        subtopo <- as.numeric(gsub(ncol(phi)+1, 0, apply(P, 1, function(x) return(which.max(x)))))
        theta <- t(R)*0
        theta[cbind(subtopo, 1:ncol(theta))] <- 1
        Oold <- O
        ll <- sum(diag(theta%*%P))
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
        ## we have to fill up the theta, however, this seems to worsen convergence and also places incorrect ones, needs improvement:
        nogenes <- which(apply(theta, 1, sum) == 0)
        nozeros <- which(t(P) > 0, arr.ind = TRUE)
        nozeros <- nozeros[which(nozeros[, 1] %in% nogenes), ]
        theta[nozeros] <- 1
        theta[grep("_", colnames(phi)), ] <- 0
        if (useF) {
            O <- t(R)%*%t(phi%*%theta)
        } else {
            O <- t(R)%*%t(theta)
        }
        if (useCut) {
            cutoff <- cut*max(abs(O))
            phi[which(O > cutoff & E == 1)] <- 1
            phi[which(O <= cutoff | E == 0)] <- 0
            phi <- transitive.closure(phi, mat = TRUE)
        } else {
            O <- O*E
            supertopo <- as.numeric(gsub(ncol(phi)+1, 0, apply(O, 1, function(x) return(which.max(x)))))
            phi <- phi*0
            phi[cbind(supertopo, 1:ncol(phi))] <- 1
            phi <- transitive.closure(phi, mat = TRUE)
        }
    }
    nem <- list(phi = phibest, theta = thetabest, iter = iter, ll = llbest, lls = lls, num = numbest, C = Cz, O = Obest, E = E0)
    class(nem) <- "nemEst" 
    return(nem)
}

lnem <- function(data, inference = "nemEst", control = NULL, donot = TRUE, ...) {
    kds <- unlist(lapply(colnames(data), function(x) return(length(unlist(strsplit(x, "_"))))))
    lnems <- list()
    data2 <- data[, which(kds <= 1)]
    if (donot) {
        for (i in 1:ncol(data2)) {
            data3 <- data2
            not <- paste0("!", colnames(data2)[i])
            data3[, i] <- -data3[, i]
            if (inference %in% "nemEst") {
                tmp <- nemEst(data3, ...)
            } else {
                colnames(data3)[ncol(data3)] <- "not"
                if (is.null(control)) {
                    if (all(data %in% c(0,1))) {
                        control <- set.default.parameters(setdiff(unique(colnames(data3)),"time"))
                    } else {
                        control <- set.default.parameters(setdiff(unique(colnames(data3)),"time"))
                        control$type <- "CONTmLLBayes"
                    }
                }
                tmp2 <- nem(data3, inference = inference, control = control, ...)
                tmp <- list()
                tmp$phi <- graph2adj(tmp2$graph)
                rownames(tmp$phi)[which(rownames(tmp$phi) %in% "not")] <- colnames(tmp$phi)[which(rownames(tmp$phi) %in% "not")] <- not
            }
            diag(tmp$phi) <- 0
            if (any(tmp$phi[i, ] == 1)) {
                diag(tmp$phi) <- 1
                lnems[[not]] <- tmp
            }
        }
    }
    for (i in 2:max(kds)) {
        for (j in which(kds == i)) {
            ## positive
            data3 <- cbind(data2, data[, j, drop = FALSE])
            multi <- colnames(data3)[ncol(data3)]
            if (inference %in% "nemEst") {
                tmp <- nemEst(data3, ...)
            } else {
                colnames(data3)[ncol(data3)] <- "multi"
                if (is.null(control)) {
                    if (all(data %in% c(0,1))) {
                        control <- set.default.parameters(setdiff(unique(colnames(data3)),"time"))
                    } else {
                        control <- set.default.parameters(setdiff(unique(colnames(data3)),"time"))
                        control$type <- "CONTmLLBayes"
                    }
                }
                tmp2 <- nem(data3, inference = inference, control = control, ...)
                tmp <- list()
                tmp$phi <- graph2adj(tmp2$graph)
                rownames(tmp$phi)[which(rownames(tmp$phi) %in% "multi")] <- colnames(tmp$phi)[which(rownames(tmp$phi) %in% "multi")] <- multi
            }
            singles <- unlist(strsplit(multi, "_"))
            direct <- colnames(tmp$phi)[which(apply(tmp$phi[which(colnames(tmp$phi) %in% singles), ], 2, sum) >= 1)]
            if (any(tmp$phi[which(rownames(tmp$phi) %in% c(multi)),
                            which(!(colnames(tmp$phi) %in% c(multi, singles, direct)))] == 1)) {
                lnems[[multi]] <- tmp
            }
            if (donot) {
                ## negative
                data3 <- cbind(data2, -data[, j, drop = FALSE])
                multi2 <- paste0("!", colnames(data3)[ncol(data3)])
                if (inference %in% "nemEst") {
                    tmp <- nemEst(data3, ...)
                } else {
                    colnames(data3)[ncol(data3)] <- "multi"
                    if (is.null(control)) {
                        if (all(data %in% c(0,1))) {
                            control <- set.default.parameters(setdiff(unique(colnames(data3)),"time"))
                        } else {
                            control <- set.default.parameters(setdiff(unique(colnames(data3)),"time"))
                            control$type <- "CONTmLLBayes"
                        }
                    }
                    tmp2 <- nem(data3, inference = inference, control = control, ...)
                    tmp <- list()
                    tmp$phi <- graph2adj(tmp2$graph)
                    rownames(tmp$phi)[which(rownames(tmp$phi) %in% "multi")] <- colnames(tmp$phi)[which(rownames(tmp$phi) %in% "multi")] <- multi
                }
                singles <- unlist(strsplit(multi, "_"))
                direct <- colnames(tmp$phi)[which(apply(tmp$phi[which(colnames(tmp$phi) %in% singles), ], 2, sum) >= 1)]
                if (any(tmp$phi[which(rownames(tmp$phi) %in% c(multi)),
                                which(!(colnames(tmp$phi) %in% c(multi, singles, direct)))] == 1)) {
                    lnems[[multi2]] <- tmp
                }
            }
        }
    }
    lnems <<- lnems
    dnf <- NULL
    for (i in 1:length(lnems)) {
        if (length(grep("!", names(lnems)[i])) > 0) {
            input <- gsub("!", "", unlist(strsplit(names(lnems)[i], "_")))
            dnf <- c(dnf, paste0(paste0("!", input, "="), colnames(lnems[[i]]$phi)[which(lnems[[i]]$phi[grep(gsub("!", "", names(lnems)[i]), rownames(lnems[[i]]$phi)), ] == 1)]))
        } else {
            tmp <- ladj(lnems[[i]]$phi)
            dnf <- c(dnf, tmp)
        }
    }
    dnf <- absorption3(unique(dnf))
    dnf <- transRed(dnf)
    return(list(lnems = lnems, dnf = dnf))
}

ladj <- function(adj) {
    for (i in notgrep("_", colnames(adj))) {
        adj[i, grep(rownames(adj)[i], rownames(adj))] <- 1
        adj[grep(rownames(adj)[i], rownames(adj)), i] <- 0
        for (j in which(adj[, i] == 1)) {
            if (adj[j, i] == 1) {
                adj[grep(paste(unlist(strsplit(rownames(adj)[j], "_")), collapse = ".*"), rownames(adj)), i] <- 0
                adj[j, i] <- 1
            }
        }
    }
    adj[grep("_", rownames(adj)), grep("_", colnames(adj))] <- 0
    adj[, intersect(which(apply(adj, 1, sum) == 0), grep("_", rownames(adj)))] <- 0
    rownames(adj)[grep("_", rownames(adj))] <- colnames(adj)[grep("_", colnames(adj))] <- paste0("AND", 1:length(grep("_", colnames(adj))))
    dnf <- NULL
    for (i in 1:nrow(adj)) {
        targets <- which(adj[i, ] == 1)
        for (j in targets) {
            if (length(grep("AND", colnames(adj)[j])) > 0) {
                dnf <- c(dnf, paste0(paste(sort(unique(c(rownames(adj)[i], rownames(adj)[which(adj[, j] == 1)]))), collapse = "+"), "=", rownames(adj)[which(adj[j, ] == 1)]))
            } else {
                dnf <- c(dnf, paste(c(rownames(adj)[i], "=", rownames(adj)[j]), collapse = ""))
            }
        }
    }
    dnf <- unique(dnf)
    dnf <- dnf[notgrep("AND", dnf)]
    dnf <- dnf[grep("=.*$", dnf)]
    return(dnf)
}

absorption3 <- function(bString, model = NULL) {
    graph <- bString
    degree <- unlist(lapply(graph, function(x) return(length(unlist(strsplit(x, "\\+"))))))
    identical <- NULL
    for (i in graph) {
        targets <- grep(paste("(?=.*", gsub("\\+", ")(?=.*", 
                                            gsub("=", ")(?=.*=", i)), ")", sep = ""), graph, 
                        perl = TRUE)
        toomuch <- grep(paste("!", gsub("\\+", "|!", gsub("=.*", 
                                                          "", i)), "", sep = ""), graph[targets])
        if (length(toomuch) > 0) {
            targets <- targets[-grep(paste("!", gsub("\\+", "|!", 
                                                     gsub("=.*", "", i)), "", sep = ""), graph[targets])]
        }
        if (length(targets) > 1) {
            identical <- c(identical, intersect(targets, which(degree == degree[which(graph %in% i)]))[-1])
            targets <- targets[-which(targets == which(graph %in% 
                                                       i))]
            targets <- intersect(targets, which(degree > degree[which(graph %in% i)]))
            if (is.null(model)) {
                if (sum(bString %in% graph[targets]) > 0) {
                    bString <- bString[-which(bString %in% graph[targets])]
                }
            }
            else {
                bString[which(model$reacID %in% graph[targets])] <- 0
            }
        }
    }
    if (is.null(model) & length(identical) > 0) {
        bString <- bString[-identical]
    }
    return(bString)
}

andgrep <- function(patterns, targets, fun = grep) {
    found <- as.numeric(names(which(table(unlist(lapply(patterns, fun, targets))) == length(patterns))))
    return(found)
}

notgrep <- function(pattern, targets) {
    found <- (1:length(targets))[-grep(pattern, targets)]
    return(found)
}

kernelnem <- function(R) {
    K <- t(R)%*%R*0
    for (i in 1:ncol(R)) {
        for (j in 1:ncol(R)) {
            Ri <- R[, i]
            Rj <- R[, j]
            Rj[which(Ri > 0 & Rj < 0)] <- 1
            K[i, j] <- t(Ri)%*%Rj
        }
    }
    K <- (K + t(K))/2
    return(K)
}

print.nemEst <- function(nem) {
}

clustNEM <- function(data, k = 2:5, ...) {
    smax <- 0
    K <- 1
    res <- NULL
    for (i in k) {
        d <- (1 - cor(data))/2
        d <- as.dist(d)
        kres <- kmeans(d, i)
        sres <- silhouette(kres$cluster, d)
        print(sum(sres[, 3]))
        if (sum(sres[, 3]) > smax) {
            Kres <- kres
            K <- i
            smax <- sum(sres[, 3])
        }
    }
    res <- list()
    for (i in 1:K) {
        if (sum(Kres$cluster == i) > 1) {
            res[[i]] <- mynem(data[, which(Kres$cluster == i)], ...)
            rownames(res[[i]]$adj) <- colnames(res[[i]]$adj) <- unique(naturalsort(names(which(Kres$cluster == i))))
        } else {
            res[[i]] <- list()
            res[[i]]$adj <- matrix(1, 1, 1)
            rownames(res[[i]]$adj) <- colnames(res[[i]]$adj) <- names(which(Kres$cluster == i))
        }
    }
    res$comp <- list()
    res$mw <- numeric(K)
    Sgenes <- length(unique(colnames(data)))
    for (i in 1:K) {
        res$comp[[i]] <- list()
        tmp <- res[[i]]$adj
        ## the next if biases towards clustNEM
        if (nrow(tmp) < Sgenes) {
            print("test")
            smiss <- unique(colnames(data)[-which(colnames(data) %in% colnames(tmp))])
            tmp <- rbind(cbind(tmp, matrix(0, nrow = nrow(tmp), ncol = length(smiss))), matrix(0, nrow = length(smiss), ncol = ncol(tmp) + length(smiss)))
            colnames(tmp)[(dim(res[[i]]$adj)+1):nrow(tmp)] <- rownames(tmp)[(dim(res[[i]]$adj)+1):nrow(tmp)] <- smiss
            tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
        }
        res$comp[[i]]$phi <- tmp
        res$mw[i] <- sum(Kres$cluster == i)/ncol(data)
    }
    res$probs <- matrix(res$mw, K, ncol(data))
    return(res)
}

## ideas:
## just use the time data for one static network, done!
## assume an evolving network:
tinem <- function(D, search = "greedy", start = NULL, method = "llr",
                  parallel = NULL, reduce = FALSE, weights = NULL, runs = 1, verbose = FALSE, redSpace = NULL,
                  trans.close = TRUE, subtopo = NULL, prior = NULL, ratio = TRUE, timeseries = NULL, max_iter = 100, evolution = TRUE) {
    if (is.null(timeseries)) {
        res <- mynem(D = data, search = search, start = start, method = method,
                     parallel = parallel, reduce = reduce, runs = runs,
                     verbose = verbose, redSpace = redSpace, ratio = ratio)
    } else {
        if (evolution) {
            res <- list()
            for (i in 1:ncol(timeseries)) { # maybe join some together?
                if (verbose) {
                    print(paste("timestep: ", i, sep = ""))
                }
                res[[i]] <- tinem(D, search = search, start = start, method = method,
                                  parallel = parallel, reduce = reduce, weights = weights,
                                  runs = runs, verbose = verbose, redSpace = redSpace,
                                  trans.close = trans.close, subtopo = subtopo, prior = prior,
                                  ratio = ratio, timeseries = matrix(timeseries[, i], ncol = 1), max_iter = max_iter,
                                  evolution = FALSE)
            }
        } else {
            timeseries2 <- timeseries
            score <- -sum(abs(data)) -sum(abs(timeseries))
            oldscore <- -Inf
            iter <- 0
            ## intelligent prelabelling:
            Sgenes <- length(unique(colnames(data)))
            datacor <- cor(data, timeseries)
            timelab <- apply(datacor, 2, which.max)
            timelab[which(timelab > Sgenes)] <- timelab[which(timelab > Sgenes)] %% Sgenes
            timelab[which(timelab == 0)] <- Sgenes
            ## timelab <- sample(unique(colnames(data)), ncol(timeseries), replace = TRUE)
            while(score > oldscore & iter < max_iter) {
                iter <- iter + 1
                oldscore <- score
                colnames(timeseries2) <- timelab
                fulldata <- cbind(data, timeseries2)
                res <- mynem(D = fulldata, search = search, start = start, method = method,
                             parallel = parallel, reduce = reduce, runs = runs,
                             verbose = verbose, redSpace = redSpace, ratio = ratio)
                subtopo <- res$subtopo
                adj1 <- transitive.closure(res$adj, mat = TRUE)
                adj2 <- adj1[, subtopo]
                tmp <- llrScore(t(timeseries), t(adj2), ratio = ratio)
                timelab <- apply(tmp, 1, which.max)
                if (verbose) {
                    print(paste("iteration: ", iter, sep = ""))
                    print(paste("score: ", score, sep = ""))
                }
                score <- res$score
            }
            res$timelab <- timelab
        }
    }

    return(res)
    
}
#' @noRd
#' @export
modules <- function(D, method = "llr", weights = NULL, reduce = FALSE, start = NULL,
                    verbose = FALSE, trans.close = TRUE, redSpace = NULL,
                    subtopo = NULL, ratio = TRUE, parallel = NULL, prior = NULL,
                    modulesize = 4, search = "exhaustive", domean = TRUE) {
    D <- data <- modData(D)
    n <- getSgeneN(D)
    Sgenes <- getSgenes(D)
    if (domean) {
        mD <- matrix(0, nrow(D), length(Sgenes))
        if (!is.null(weights)) {
            D <- t(t(D)*weights)
            weights <- rep(1, ncol(mD))
        }
        for (i in 1:length(Sgenes)) {
            mD[, i] <- apply(D[, which(colnames(D) %in% i), drop = FALSE], 1, mean)
        }
        D <- mD
        colnames(D) <- 1:length(Sgenes)
        sumdata <- data <- D
    } else {
        sumdata <- matrix(0, nrow(data), n)
        if (!is.null(weights)) {
            D <- t(t(D)*weights)
            weights <- rep(1, ncol(sumdata))
        }
        for (i in 1:n) {
            sumdata[, i] <- apply(D[, which(colnames(D) %in% i), drop = FALSE], 1, sum)
        }
        colnames(sumdata) <- 1:n
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
    for (i in 1:length(groups)) {
        subset <- which(hcut == i)
        if (verbose) { print(paste(c("calculating module", subset), collapse = " ")) }
        if (length(subset) > 1) {
            subdata <- data[, which(colnames(data) %in% subset)]
            if (is.null(start)) {
                start2 <- start
            } else {
                start2 <- start[which(rownames(start) %in% subset), which(colnames(start) %in% subset)]
            }
            tmp <- mynem(subdata, search = search, method = method, start = start2,
                         parallel = parallel, reduce = reduce,
                         weights = weights[which(colnames(data) %in% subset)], verbose = verbose,
                         redSpace = redSpace, trans.close = trans.close,
                         subtopo = subtopo, prior = prior, ratio = ratio,
                         domean = FALSE)
            if (is.null(fullnet)) {
                fullnet <- tmp$adj
            } else {
                tmpnames <- c(colnames(fullnet), colnames(tmp$adj))
                fullnet <- rbind(cbind(fullnet, matrix(0, nrow(fullnet), ncol(tmp$adj))),
                                 cbind(matrix(0, nrow(tmp$adj), ncol(fullnet)), tmp$adj))
                colnames(fullnet) <- rownames(fullnet) <- as.numeric(tmpnames)
            }
        } else {
            if (is.null(dim(fullnet))) {
                fullnet <- matrix(1, 1, 1)
                colnames(fullnet) <- rownames(fullnet) <- subset
            } else {
                fullnet <- rbind(cbind(fullnet, 0), 0)
                colnames(fullnet)[ncol(fullnet)] <- rownames(fullnet)[nrow(fullnet)] <- subset
            }
        }
    }
    fullnet <- transitive.reduction(fullnet)
    fullnet <- fullnet[order(as.numeric(rownames(fullnet))), order(as.numeric(colnames(fullnet)))]
    return(fullnet)
}
#' @noRd
#' @export
getSgeneN <- function(data) {
    Sgenes <- length(unique(unlist(strsplit(colnames(data), ","))))
    return(Sgenes)
}
#' @noRd
#' @export
getSgenes <- function(data) {
    Sgenes <- sort(as.numeric(unique(unlist(strsplit(colnames(data), ",")))))
    return(Sgenes)
}
#' @noRd
#' @export
mynem <- function(D, search = "greedy", start = NULL, method = "llr",
                  parallel = NULL, reduce = FALSE, weights = NULL, runs = 1,
                  verbose = FALSE, redSpace = NULL,
                  trans.close = TRUE, subtopo = NULL, prior = NULL,
                  ratio = TRUE, domean = TRUE, modulesize = 5, ...) { # reduce might not work as expected
    if ("modules" %in% search) {
        if (length(search) > 1) {
            search <- search[-which(search %in% "modules")]
        } else {
            search <- "greedy"
        }
        if (length(unique(colnames(D))) > modulesize) {
            start <- modules(D, method = method, weights = weights, reduce = reduce, verbose = verbose, start = start,
                             trans.close = trans.close, redSpace = redSpace, subtopo = subtopo,
                             ratio = ratio, parallel = parallel, prior = prior,
                             modulesize = modulesize, search = search, domean = domean)
        }
        if (search %in% "exhaustive") {
            search <- "greedy"
        }
    }
    D.backup <- D
    D <- modData(D)
    colnames(D) <- gsub("\\..*", "", colnames(D))
    Sgenes <- getSgenes(D)
    if (domean) {
        mD <- matrix(0, nrow(D), length(Sgenes))
        if (!is.null(weights)) {
            D <- t(t(D)*weights)
            weights <- rep(1, ncol(mD))
        }
        for (i in 1:length(Sgenes)) {
            mD[, i] <- apply(D[, which(colnames(D) %in% i), drop = FALSE], 1, mean)
        }
        D <- mD
        colnames(D) <- 1:length(Sgenes)
    }
    if (is.null(start)) {
        start2 <- "full"
        start <- better <- matrix(0, length(Sgenes), length(Sgenes))
    } else {
        if (length(start) == 1) {
            start2 <- start
            start <- better <- matrix(0, length(Sgenes), length(Sgenes))
        } else {
            better <- start
        }
    }
    diag(start) <- diag(better) <- 1
    colnames(better) <- rownames(better) <- colnames(start) <- rownames(start) <- Sgenes
    score <- scoreAdj(D, better, method = method, weights = weights, subtopo = subtopo,
                      prior = prior, ratio = ratio)
    score <- score$score
    oldscore <- score
    allscores <- score
    traClo <- nem::transitive.closure
    
    if (!is.null(parallel)) {
        transitive.closure <- nem::transitive.closure
        transitive.reduction <- nem::transitive.reduction
        sfInit(parallel = TRUE, cpus = parallel)
        sfExport("modules", "D", "start", "better", "traClo", "method", "scoreAdj",
                 "weights", "transitive.closure", "llrScore",
                 "transitive.reduction")
    }

    if (search %in% "greedy") {
        for (iter in 1:runs) {
            if (iter > 1) {
                better <- matrix(sample(c(0,1), nrow(better)*ncol(better), replace = TRUE), nrow(better), ncol(better))
                colnames(better) <- rownames(better) <- sample(Sgenes, length(Sgenes))
                better[lower.tri(better)] <- 0
                diag(better) <- 1
                better <- better[order(as.numeric(rownames(better))), order(as.numeric(colnames(better)))]
                score <- scoreAdj(D, better, method = method, weights = weights,
                                  subtopo = subtopo, prior = prior, ratio = ratio)
                subtopo <- score$subtopo
                score <- score$score
                oldscore <- score
                allscores <- score
            }
            stop <- F
            while(!stop) {
                doScores <- function(i) {
                    ## new <- better
                    ## new[i] <- 1 - new[i]
                    ## new <- traClo(new, mat = TRUE, loops = TRUE)
                    new <- models[[i]]
                    score <- scoreAdj(D, new, method = method, weights = weights,
                                      subtopo = subtopo, prior = prior,
                                      ratio = ratio)
                    subtopo <- score$subtopo
                    score <- score$score
                    return(score)
                }
                models <- unique(c(get.insertions(better), get.reversions(better), get.deletions(better)))
                if (is.null(parallel)) {
                    scores <- unlist(lapply((1:length(models)), doScores))
                } else {
                    scores <- unlist(sfLapply((1:length(models)), doScores))
                }
                scores[is.na(scores)] <- 0
                best <- models[[which.max(scores)]]
                best <- traClo(best, mat = TRUE)
                if (max(scores, na.rm = TRUE) > oldscore | (max(scores, na.rm = TRUE) == oldscore & sum(better == 1) > sum(best == 1))) {
                    better <- best
                    better <- traClo(better, mat = TRUE)
                    oldscore <- max(scores)
                    allscores <- c(allscores, oldscore)
                } else {
                    stop <- T
                }
            }
            if (iter > 1) {
                if (oldscore > oldscore2) {
                    better2 <- better
                    allscores2 <- allscores
                    oldscore2 <- oldscore
                }
            } else {
                better2 <- better
                allscores2 <- allscores
                oldscore2 <- oldscore
            }
        }
        better <- better2
        allscores <- allscores2
        oldscore <- oldscore2
    }

    if (search %in% "exhaustive") {
        models <- enumerate.models(length(Sgenes), Sgenes, trans.close = trans.close, verbose = verbose)
        doScores <- function(i) {
            adj <- models[[i]]
            score <- scoreAdj(D, adj, method = method, weights = weights,
                              subtopo = subtopo, prior = prior,
                              ratio = ratio)
            subtopo <- score$subtopo
            score <- score$score
            return(score)
        }
        if (is.null(parallel)) {
            scores <- unlist(lapply(1:length(models), doScores))
        } else {
            scores <- unlist(sfLapply(1:length(models), doScores))
        }
        best <- which.max(scores)
        better <- traClo(models[[best]], mat = TRUE)
        diag(better) <- 1
    }

    if (search %in% "exhaustive2") {
        dec2bin <- function(x) {
            if (x == 0) {
                y <- 0
            } else {
                tmp <- rev(as.integer(intToBits(x)))
                y <- tmp[which(tmp != 0)[1]:length(tmp)]
            }
            return(y)
        }
        n <- length(Sgenes)
        spaceExp <- 2^(n*n)
        ## if (verbose) { print(paste("different network structures: ", spaceExp, sep = "")) }
        if (!is.null(redSpace)) {
            spaceUnique <- redSpace
            ## if (verbose) { print(paste("different equivalence classes: ", length(spaceUnique), sep = "")) }
        } else {
            if (reduce) {
                createFull <- function(i) {
                    bivec <- dec2bin(i-1)
                    bivec <- c(rep(0, n*n - length(bivec)), bivec)
                    adj <- traClo(matrix(bivec, n, n), mat = TRUE)
                    diag(adj) <- 1
                    return(paste(as.vector(adj), collapse = ""))
                }
                if (is.null(parallel)) {
                    spaceFull <- do.call("c", lapply(1:spaceExp, createFull))
                } else {
                    spaceFull <- do.call("c", sfLapply(1:spaceExp, createFull))
                }
                spaceUnique <- (1:spaceExp)[-which(duplicated(spaceFull) == TRUE)]
                redSpace <- spaceUnique
                ## if (verbose) { print(paste("different equivalence classes: ", length(spaceUnique), sep = "")) }
            } else {
                spaceUnique <- 1:spaceExp
            }
        }
        scores <- numeric(spaceExp)
        doScores <- function(i) {
            bivec <- dec2bin(i-1)
            bivec <- c(rep(0, n*n - length(bivec)), bivec)
            adj <- traClo(matrix(bivec, n, n), mat = TRUE)
            diag(adj) <- 1
            colnames(adj) <- rownames(adj) <- Sgenes
            score <- scoreAdj(D, adj, method = method, weights = weights,
                              subtopo = subtopo, prior = prior,
                              ratio = ratio)
            subtopo <- score$subtopo
            score <- score$score
            return(score)
        }
        if (is.null(parallel)) {
            scores[spaceUnique] <- unlist(lapply(spaceUnique, doScores))
        } else {
            scores[spaceUnique] <- unlist(sfLapply(spaceUnique, doScores))
        }
        bestbivec <- dec2bin(which.max(scores)-1)
        bestbivec <- c(rep(0, n*n - length(bestbivec)), bestbivec)
        better <- traClo(matrix(bestbivec, n, n), mat = TRUE)
        colnames(better) <- rownames(better) <- Sgenes
        diag(better) <- 1
    }

    if (search %in% "estimate") {
        if (!is.null(weights)) {
            Dw <- t(t(D)*weights)
        } else {
            Dw <- D
        }
        tmp <- nemEst(Dw, start = start2, ...)
        subtopo <- apply(tmp$theta, 2, function(x) return(which(x == 1)))
        better <- tmp$phi
        oldscore <- tmp$ll
        allscores <- tmp$lls
        subweights <- t(Dw%*%cbind(tmp$phi, 0))
    }

    if (!is.null(parallel)) {
        sfStop()
    }

    if (is.null(subtopo)) {
        subtopo <- scoreAdj(D, better, method = method, weights = weights,
                            subtopo = subtopo, prior = prior,
                            ratio = ratio)
        subweights <- subtopo$subweights
        subtopo <- subtopo$subtopo
    }

    better <- nem::transitive.reduction(better)
    better <- better[order(as.numeric(rownames(better))), order(as.numeric(colnames(better)))]
    nem <- list(adj = better, score = oldscore, scores = allscores, redSpace = redSpace, subtopo = subtopo, D = D.backup, subweights = subweights)
    class(nem) <- "mynem"
    return(nem)
}
#' @noRd
#' @export
plot.mynem <- function(x, ...) {
    translate <- sort(unique(colnames(x$D)))

    adj <- x$adj

    for (i in 1:length(translate)) {
        colnames(adj)[which(colnames(adj) %in% i)] <- translate[i]
        rownames(adj)[which(rownames(adj) %in% i)] <- translate[i]
    }

    plot.adj(adj)
    
}

#' simulate single cell data from a mixture of networks
#' @param Sgenes number of Sgenes
#' @param Egenes number of Egenes
#' @param subsample range to subsample data. 1 means the full simulated data is used.
#' @param Nems numberof components
#' @param reps number of relicates, if set (not realistice for cells)
#' @param mw mixture weights
#' @param evolution evovling network, if set to true
#' @param nCells number of cells
#' @param uninform number of uninformative Egenes
#' @param unitheta uniform theta, if true
#' @author Martin Pirkl
#' @return simulation object with meta information and data
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
simData <- function(Sgenes = 5, Egenes = 1, subsample = 1,
                    Nems = 2, reps = NULL, mw = NULL, evolution = FALSE,
                    nCells = 1000, uninform = 0, unitheta = FALSE, edgeprob = 0.5) {
    Nem <- list()
    data <- NULL
    index <- NULL
    theta <- list()
    for (i in 1:Nems) {
        if (i == 1 | !evolution) {
            adj <- matrix(sample(c(0,1), Sgenes^2, replace = TRUE, prob = c(1-edgeprob, edgeprob)), Sgenes, Sgenes)
            adj <- adj[order(apply(adj, 1, sum), decreasing = TRUE), order(apply(adj, 2, sum), decreasing = FALSE)]
            adj[lower.tri(adj)] <- 0
            diag(adj) <- 1
            adj <- nem::transitive.closure(adj, mat = TRUE)
            colnames(adj) <- rownames(adj) <- sample(1:Sgenes, Sgenes)
            adj <- adj[order(as.numeric(rownames(adj))), order(as.numeric(colnames(adj)))]
        } else {
            children <- which(apply(Nem[[(i-1)]], 1, sum) == 0)
            parents <- which(apply(Nem[[(i-1)]], 2, sum) == 0)
            if (length(children) > 1) {
                child <- sample(children, 1)
            } else {
                child <- children
            }
            parents <- c(which(apply(Nem[[(i-1)]], 2, sum) == 0), child)
            if (length(parents) > 1) {
                if (child %in% parents) {
                    parent <- sample(parents[-which(parents %in% child)], 1)
                } else {
                    parent <- sample(parents, 1)
                }
            } else {
                parent <- parents
            }
            adj <- Nem[[(i-1)]]
            adj[, child] <- 0
            adj[child, parent] <- 1
            adj <- transitive.closure(adj, mat = TRUE)
        }
        if (is.null(reps)) {
            reps2 <- ceiling(nCells/Sgenes)
        } else {
            reps2 <- reps
        }
        Nem[[i]] <- nem::transitive.reduction(adj)
        data_tmp <- t(adj)
        colntmp <- rep(1:ncol(data_tmp), reps2)
        data_tmp <- data_tmp[, rep(1:ncol(data_tmp), reps2)]
        colnames(data_tmp) <- colntmp
        if (!is.null(mw)) {
            tmpsamp <- sample(1:ncol(data_tmp), ceiling(mw[i]*ncol(data_tmp)))
            data_tmp <- data_tmp[, tmpsamp, drop = FALSE] # since this is a uniform sampling, it's ok
        } else {
            if (is.null(reps)) {
                data_tmp <- data_tmp[, 1:ceiling(nCells/Nems)]
            }
        }
        index <- c(index, rep(i, ncol(data_tmp)))
        data_tmp <- data_tmp[rep(1:Sgenes, each = Egenes), , drop = FALSE]
        if (!unitheta) {
            eorder <- sample(1:nrow(data_tmp), nrow(data_tmp))
            data_tmp <- data_tmp[eorder, ]
            theta[[i]] <- rownames(data_tmp)[eorder]
        }
        data <- cbind(data, data_tmp)
    }
    if (subsample < 1) { # exchange subsample with mixture weights
        subsample <- sample(1:ncol(data), ceiling(ncol(data)*subsample))
        data <- data[, subsample] # set reps high and subsample low/high to make it more realistic
        index <- rep(1:Nems, each = Sgenes*reps)[subsample]
    }
    if (uninform > 0) {
        data <- rbind(data, matrix(sample(c(0,1), ncol(data)*uninform, replace = TRUE), uninform, ncol(data)))
    }
    return(list(Nem = Nem, theta = theta, data = data, index = index))
}
#' @noRd
#' @export
hamSim <- function(a, b, diag = 1, symmetric = TRUE) {
    Sgenes <- unique(colnames(a))
    ham <- numeric(ncol(b))
    for (i in 1:ncol(b)) {
        c <- b[, i]
        tmp <- numeric(sum(colnames(a) %in% colnames(b)[i]))
        count <- 0
        for (j in which(colnames(a) %in% colnames(b)[i])) {
            d <- a[, j]
            count <- count + 1
            tmp[count] <- 1 - sum(abs(c - d))/(nrow(b) - 1*diag)
        }
        ham[i] <- max(tmp)
    }
    ham2 <- sum(ham)/length(ham)
    if (symmetric) {
        for (i in 1:ncol(a)) {
            c <- a[, i]
            tmp <- numeric(sum(colnames(b) %in% colnames(a)[i]))
            count <- 0
            for (j in which(colnames(b) %in% colnames(a)[i])) {
                d <- b[, j]
                count <- count + 1
                tmp[count] <- 1 - sum(abs(c - d))/(nrow(a) - 1*diag)
            }
            ham[i] <- max(tmp)
        }
        ham3 <- sum(ham)/length(ham)
    } else {
        ham3 <- 1
    }
    ham <- min(c(ham2, ham3))
    return(ham)
}

## functions from the nem package. update to use the native nem package functions!
#' @noRd
#' @export
get.insertions = function(Phi, trans.close=TRUE){
    idx = which(Phi == 0)
    models = list()
    if(length(idx) > 0){
        for(i in 1:length(idx)){ # test all possible new edges
            Phinew = Phi
            Phinew[idx[i]] = 1
            if(trans.close)
                Phinew = transitive.closure(Phinew, mat=TRUE,loops=TRUE) 
            models[[i]] <- Phinew
        }
    } 
    models       
}
#' @noRd
#' @export
get.deletions = function(Phi){
    Phi = Phi - diag(ncol(Phi))
    idx = which(Phi == 1)
    models = list()
    if(length(idx) > 0){
        for(i in 1:length(idx)){ # test all possible edge deletions
            Phinew = Phi
            Phinew[idx[i]] = 0
            diag(Phinew) = 1
            models[[i]] <- Phinew
        }
    } 
    models       
}
#' @noRd
#' @export
get.reversions = function(Phi){
    idx = which(Phi + t(Phi) == 1, arr.ind=TRUE)
    models = list()
    if(NROW(idx) > 0){
        for(i in 1:NROW(idx)){ # test all possible edge reversions
            Phinew = Phi
            Phinew[idx[i,1],idx[i,2]] = 0
            Phinew[idx[i,2],idx[i,1]] = 1
            diag(Phinew) = 1
            models[[i]] <- Phinew
        }
    } 
    models       
}
#' @noRd
#' @export
enumerate.models <- function(x,name=NULL, trans.close=TRUE, verbose=TRUE) {

if (length(x) == 1) {
            n <- as.numeric(x)
	    if(is.null(name))
            	name <- letters[1:n]
        } else {
           n <- length(x)
           name <- x
        }

#------------------
# Sanity checks    

 if (n==1) stop("nem> choose n>1!")
 if (n>5)  stop("nem> exhaustive enumeration not feasible with more than 5 perturbed genes")            
 if (n==5) cat ("nem> this will take a while ... \n") 

#------------------

  bc <- bincombinations(n*(n-1))
  fkt1 <- function(x,n,name) {
    M <- diag(n)
    M[which(M==0)]<-x
    dimnames(M) <- list(name,name)	
    if(trans.close)    
    	M <- transitive.closure(M,mat=TRUE,loops=TRUE)    
    return(list(M))
  }
  
  models <- apply(bc,1,fkt1,n,name) 
  models <- unique(matrix(unlist(models),ncol=n*n,byrow=TRUE))
  
  fkt2 <- function(x,n,name){
     M <- matrix(x,n)
     dimnames(M) <- list(name,name)
     return(list(M))
  }
  models <- unlist(apply(models,1,fkt2,n,name),recursive=FALSE)
    
  if (verbose) cat("Generated",length(models),"unique models ( out of", 2^(n*(n-1)), ")\n")

  return(models)
}

#' @noRd
#' @export
adj2dnf <- function(A) {

  dnf <- NULL
  
  for (i in 1:ncol(A)) {
    for (j in 1:nrow(A)) {
      if (i %in% j) { next() }
      if (A[i, j] == 1) {
        dnf <- c(dnf, paste(colnames(A)[i], rownames(A)[j], sep = "="))
      }
      if (A[i, j] == -1) {
        dnf <- c(dnf, paste("!", colnames(A)[i], "=", rownames(A)[j], sep = ""))
      }
    }
  }

  dnf <- unique(dnf)
  
  return(dnf)

}
#' @noRd
#' @export
plot.adj <- function(x) {
    adj2graph <- function(adj.matrix) {
        V   <- rownames(adj.matrix)
        edL <- vector("list", length=nrow(adj.matrix))
        names(edL) <- V
        for (i in 1:nrow(adj.matrix)) {
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
#' @export
graph2adj <- function(gR) {
    adj.matrix <- matrix(0,
                         length(nodes(gR)),
                         length(nodes(gR))
                         )
    rownames(adj.matrix) <- nodes(gR)
    colnames(adj.matrix) <- nodes(gR)
    for (i in 1:length(nodes(gR))) {
        adj.matrix[nodes(gR)[i],adj(gR,nodes(gR)[i])[[1]]] <- 1
    }
    return(adj.matrix)
}
#' @noRd
#' @export
scoreAdj <- function(D, adj, method = "llr", weights = NULL,
                     trans.close = TRUE, subtopo = NULL,
                     prior = NULL, ratio = TRUE) {
    adj <- transitive.closure(adj, mat = TRUE)
    adj1 <- cbind(adj[colnames(D), ], "0" = 0)
    
    if (method %in% "llr") {
        score <- llrScore(D, adj1, weights = weights, ratio = ratio)
    }
    if (is.null(subtopo)) {
        subtopo <- apply(score, 1, which.max)
    }
    subweights <- score
    ll <- "max"
    if (ll %in% "max") {
        score <- sum(score[cbind(1:nrow(score), subtopo)])#/nrow(D) # I could normalize the effects matrix to 0 1 to make this a normalized score between 0 1
    }
    if (ll %in% "marg") {
        score <- log(sum(apply(exp(score), 1, sum))/length(score))
    }
    if (!is.null(prior)) {
        prior <- transitive.reduction(prior)
        adj <- transitive.reduction(adj)
    } else {
        prior <- adj
    }
    score <- score - sum(abs(prior - adj))/length(prior)
    return(list(score = score, subtopo = subtopo, subweights = subweights))
}

#' @noRd
#' @export
adj2dnf <- function(A) {

  dnf <- NULL
  
    for (i in 1:ncol(A)) {
        dnf <- c(dnf, rownames(A))
        for (j in 1:nrow(A)) {
            ## if (i %in% j) { next() }
            if (A[i, j] == 1) {
                dnf <- c(dnf, paste(colnames(A)[i], rownames(A)[j], sep = "="))
            }
            if (A[i, j] == -1) {
                dnf <- c(dnf, paste("!", colnames(A)[i], "=", rownames(A)[j], sep = ""))
            }
        }
    }
    
    dnf <- unique(dnf)
    
    return(dnf)
    
}

#' function for visualizing graphs in normal form
#' @param dnf Hyper-graph in disjunctive normal form.
#' @param freq Frequency of hyper-edges.
#' @param stimuli Vertices which can be stimulated.
#' @param signals Vertices which regulate E-genes.
#' @param inhibitors Vertices which can be inhibited.
#' @param connected If TRUE, only includes vertices which are connected to other vertices.
#' @param CNOlist CNOlist object. Optional instead of stimuli, inhibitors or signals.
#' @param cex Global font size.
#' @param fontsize Vertice label size.
#' @param labelsize Edge label size.
#' @param type Different plot types.
#' @param lwd Line width.
#' @param edgelwd Edgeline width.
#' @param legend 0 shows no legend. 1 shows legend as a graph. 2 shows legend in a box.
#' @param x x coordinate of box legend.
#' @param y y coordinate of box legend.
#' @param xjust Justification of legend box left, right or center (-1,1,0).
#' @param yjust Justification of legend box top, bottom or middle (-1,1,0).
#' @param width Vertice width.
#' @param height Vertice height.
#' @param rankdir not used
#' @param rank not used
#' @param layout Graph layout. See graphvizCapabilities()$layoutTypes.
#' @param main Main title.
#' @param sub Subtitle.
#' @param cex.main Main title font size.
#' @param cex.sub Subtitle font size.
#' @param col.sub Font color of subtitle.
#' @param fontcolor Global font color.
#' @param nodestates Binary state of each vertice.
#' @param simulate Simulate stimulation and inhibition of a list of vertices. E.g. simulate = list(stimuli = c("A", "B"), inhibitors = c("C", "D")).
#' @param andcolor not used
#' @param edgecol Vector with colors for every edge of the graph (not hyper-graph). E.g. an AND gate consists of three distinct edges.
#' @param labels Vector with labels for the edges.
#' @param labelcol Vector with label colors for the edges.
#' @param nodelabel List of vertices with labels as input. E.g. labels = list(A="test", B="label for B").
#' @param nodecol List of vertices with colors as input.
#' @param bordercol List of vertices with colors as input.
#' @param nodeshape List of vertices with shapes (diamond, box, square,...).
#' @param verbose Verbose output.
#' @param edgestyle not used
#' @param nodeheight List of vertices with height as input.
#' @param nodewidth List of vertices with width as input.
#' @param edgewidth Vector with edge widths.
#' @param lty Vector with edge styles (line, dotted,...).
#' @param hierarchy List with the hierarchy of the vertices. E.g. list(top = c("A", "B"), bottom = c("C", "D")).
#' @param showall See "connected" above.
#' @param nodefontsize not used
#' @param edgehead Vector with edge heads.
#' @param edgelabel Vector with edge labels.
#' @param edgetail Vector with edge tails.
#' @param bool If TRUE, only shows normal graph and no AND gates.
#' @param draw Do not plot the graph and only output the graphNEL object.
#' @param \dots additional arguments
#' @author Martin Pirkl
#' @return Rgraphviz object
#' @export
#' @import
#' Rgraphviz
#' @examples
#' g <- c("!A+B=G", "C=G", "!D=G", "C+D+E=G")
#' plotDnf(g)
plotDnf <- function(dnf = NULL, freq = NULL, stimuli = c(), signals = c(), inhibitors = c(), connected = TRUE,  CNOlist = NULL, cex = NULL, fontsize = NULL, labelsize = NULL, type = 2, lwd = 1, edgelwd = 1, legend = 0, x = 0, y = 0, xjust = 0, yjust = 0, width = 1, height = 1, rankdir = "TB", rank = "same", layout = "dot", main = "", sub = "", cex.main = 1.5, cex.sub = 1, col.sub = "grey", fontcolor = NULL, nodestates = NULL, simulate = NULL, andcolor = "transparent", edgecol = NULL, labels = NULL, labelcol = "blue", nodelabel = NULL, nodecol = NULL, bordercol = NULL, nodeshape = NULL, verbose = FALSE, edgestyle = NULL, nodeheight = NULL, nodewidth = NULL, edgewidth = NULL, lty = NULL, hierarchy = NULL, showall = FALSE, nodefontsize = NULL, edgehead = NULL, edgelabel = NULL, edgetail = NULL, bool = TRUE, draw = TRUE, ...) {
    if (is.matrix(dnf)) {
        dnf <- adj2dnf(transitive.reduction(dnf))
    }

    if (!bool & length(grep("\\+", dnf)) > 0) {
        dnf <- dnf[-grep("\\+", dnf)]
    }

    graph <- dnf

    if (!is.null(hierarchy)) {
        if (!showall) {
            nodes <- unique(unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+")))
            for (i in 1:length(hierarchy)) {
                hierarchy[[i]] <- intersect(hierarchy[[i]], nodes)
            }
        }
        graph2 <- NULL
        for (i in graph) {
            input <- unlist(strsplit(i, "="))
            output <- input[2]
            input <- gsub("!", "", unlist(strsplit(input[1], "\\+")))
            for (j in input) {
                graph2 <- c(graph2, paste(j, output, sep = "="))
            }
            graph2 <- unique(graph2)
        }
        hgraph <- NULL
        hgraph2 <- NULL
        for (i in 1:(length(hierarchy)-1)) {
            for (j in hierarchy[[i]]) {
                for (k in hierarchy[[i+1]]) {
                    hgraph <- c(hgraph, paste(j, k, sep = "="))
                    hgraph2 <- c(hgraph2, paste(k, j, sep = "="))
                }
            }
        }
        if (sum(hgraph %in% graph2) > 0) {
            hgraph <- hgraph[-which(hgraph %in% graph2)]
        }
        if (sum(hgraph2 %in% graph2) > 0) {
            hgraph <- hgraph[-which(hgraph2 %in% graph2)]
        }
        dnf <- c(graph, hgraph)
        ## update all the parameters e.g. edgestyle...
        if (is.null(edgecol)) {
            edgecol <- c(rep("black", length(grep("\\+", graph))+length(unlist(strsplit(graph, "\\+")))), rep("transparent", length(hgraph)))
            dnf2 <- dnf
            if (length(grep("\\+", dnf2)) > 0) {
                dnf2[-grep("\\+", dnf2)] <- gsub("=", "", dnf2[-grep("\\+", dnf2)])
            } else {
                dnf2 <- gsub("=", "", dnf2)
            }
            edgecol[grep("!", unlist(strsplit(unlist(strsplit(dnf2, "\\+")), "=")))] <- "red"
        } else {
            if (length(edgecol) == 1) {
                edgecol <- c(rep(edgecol, length(grep("\\+", graph))+length(unlist(strsplit(graph, "\\+")))), rep("transparent", length(hgraph)))
            } else {
                edgecol <- c(edgecol, rep("transparent", length(hgraph)))
            }
        }
    } else {
        if (is.null(lty) & !is.null(dnf)) {
            lty <- c(rep("solid", length(grep("\\+", graph))+length(unlist(strsplit(graph, "\\+")))))
        }
    }

    graph <- dnf

    dolegend <- FALSE
    if (is.null(dnf)) {
        dnf <- c("A=B")
        dolegend <- TRUE
    }
    
    if (!is.null(simulate)) {
        nodestates <- simulateDnf(graph, stimuli = simulate$stimuli, inhibitors = simulate$inhibitors)
    }
    
    if (is.null(freq)) {
        use.freq = FALSE
    } else {
        use.freq = TRUE
        if (is.null(labels)) {
            labels <- as.character(round(freq, 2)*100)
        }
    }
    if (is.null(labels)) {
        labels <- rep("", length(dnf))
    }
    
    if (is.null(fontsize)) {
        fontsize <- ""
    }
    if (is.null(labelsize)) {
        labelsize <- fontsize
    }

    if (!is.null(CNOlist)) {
        if (length(stimuli) == 0) {
            stimuli <- colnames(CNOlist@stimuli)
        }
        if (length(signals) == 0) {
            signals <- colnames(CNOlist@signals[[1]])
        }
        if(length(inhibitors) == 0) {
            inhibitors <- colnames(CNOlist@inhibitors)
        }
    }

    if (connected) {
        Vneg <- unique(c(c(unique(unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+"))))))
        if (length(grep("\\+", dnf)) > 0) {
            Vneg <- c(Vneg, paste("and", 1:length(grep("\\+", dnf)), sep = ""))
        }
        V <- unique(gsub("!", "", Vneg))
        stimuli <- intersect(stimuli, V)
        signals <- intersect(signals, V)
        inhibitors <- intersect(inhibitors, V)
    } else {
        Vneg <- unique(c(c(unique(unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+")))), stimuli, signals, inhibitors))
        if (length(grep("\\+", dnf)) > 0) {
            Vneg <- c(Vneg, paste("and", 1:length(grep("\\+", dnf)), sep = ""))
        }
        V <- unique(gsub("!", "", Vneg))
    }

    V <- sort(V)
    
    Vneg <- c(V, Vneg[grep("!", Vneg)])

    if (!is.null(nodecol)) {
        if (length(nodecol) == 1 & !is.list(nodecol)) {
            nodecol.tmp <- nodecol
            nodecol <- list()
            for (i in V) {
                if (length(grep("and", i)) == 0) {
                    nodecol[[i]] <- nodecol.tmp
                }
            }
        }
    }

    if (!is.null(bordercol)) {
        if (length(bordercol) == 1 & !is.list(bordercol)) {
            bordercol.tmp <- bordercol
            bordercol <- list()
            for (i in V) {
                if (length(grep("and", i)) == 0) {
                    bordercol[[i]] <- bordercol.tmp
                }
            }
        }
    }

    E <- list()

    for (i in V) {
        E[[i]] <- list()
    }

    Eneg <- list()

    for (i in Vneg) {
        Eneg[[i]] <- list()
    }

    count <- 0
    
    for (i in dnf) {
        tmp <- unlist(strsplit(i, "="))
        if (length(tmp)==1) {
            Eneg[[tmp]][["edges"]] <- c(Eneg[[tmp]][["edges"]], NULL)
            tmp <- gsub("!", "", tmp)
            E[[tmp]][["edges"]] <- c(E[[tmp]][["edges"]], NULL)
        } else {
            tmp2 <- unlist(strsplit(tmp[1], "\\+"))
            if (length(tmp2) > 1) {
                count <- count + 1
                Eneg[[paste("and", count, sep = "")]][["edges"]] <- c(Eneg[[paste("and", count, sep = "")]][["edges"]], which(Vneg %in% tmp[2]))
                E[[paste("and", count, sep = "")]][["edges"]] <- c(E[[paste("and", count, sep = "")]][["edges"]], which(V %in% tmp[2]))
                for (j in tmp2) {
                    Eneg[[j]][["edges"]] <- c(Eneg[[j]][["edges"]], which(Vneg %in% paste("and", count, sep = "")))
                    j <- gsub("!", "", j)
                    E[[j]][["edges"]] <- c(E[[j]][["edges"]], which(V %in% paste("and", count, sep = "")))
                }
            } else {
                Eneg[[tmp2]][["edges"]] <- c(Eneg[[tmp2]][["edges"]], which(Vneg %in% tmp[2]))
                tmp2 <- gsub("!", "", tmp2)
                E[[tmp2]][["edges"]] <- c(E[[tmp2]][["edges"]], which(V %in% tmp[2]))
            }
        }
    }

    g <- new("graphNEL",nodes=V,edgeL=E,edgemode="directed")

    gneg <- new("graphNEL",nodes=Vneg,edgeL=Eneg,edgemode="directed")

    nodes <- buildNodeList(g)

    edges <- buildEdgeList(g)
    
    nodesneg <- buildNodeList(gneg)

    edgesneg <- buildEdgeList(gneg)

    edgesnew <- list()
    
    for (i in sort(names(edges))) {
        edgesnew <- c(edgesnew, edges[[i]])
    }

    names(edgesnew) <- sort(names(edges))

    edges <- edgesnew

    edgesnegnew <- list()
    
    for (i in names(edgesneg)) {#[order(gsub("!", "", names(edgesneg)))]) {
        edgesnegnew <- c(edgesnegnew, edgesneg[[i]])
    }

    names(edgesnegnew) <- names(edgesneg)#[order(gsub("!", "", names(edgesneg)))]

    edgesneg <- edgesnegnew

    if (verbose) {
        print(paste("order of nodes: ", paste(names(nodes), collapse = ", "), sep = ""))
        print(paste("order of edges: ", paste(names(edges), collapse = ", "), sep = ""))
    }

    edges <- edgesneg
    names(edges) <- gsub("!", "", names(edges))

    for (i in 1:length(edges)) {
        edges[[i]]@from <- gsub("!", "", edges[[i]]@from)
    }
    
    nodeshape2 <- nodeshape
    nodeshape <- list()
    if (length(nodeshape2) == 1 & !(is.list(nodeshape2))) {
        for (i in 1:length(nodes)) {
            nodeshape[[nodes[[i]]@name]] <- nodeshape2
        }
    } else {
        nodeshape <- nodeshape2
    }
    
    for (i in 1:length(nodes)) {
        nodes[[i]]@attrs$height <- height
        nodes[[i]]@attrs$width <- width
        if (!is.null(nodelabel[[nodes[[i]]@name]])) {
            nodes[[i]]@attrs$name <- nodelabel[[nodes[[i]]@name]]
        }
        if (!is.null(nodeheight[[nodes[[i]]@name]])) {
            nodes[[i]]@attrs$height <- nodeheight[[nodes[[i]]@name]]
        }
        if (!is.null(nodewidth[[nodes[[i]]@name]])) {
            nodes[[i]]@attrs$width <- nodewidth[[nodes[[i]]@name]]
        }
        if (length(grep("and", nodes[[i]]@name)) > 0) {
            if (is.null(nodelabel)) {
                nodelabel <- list()
            }
            nodelabel[[nodes[[i]]@name]] <- "AND"
            nodes[[i]]@attrs$label <- ""
            nodes[[i]]@attrs$fontcolor <- andcolor
            if (is.null(bordercol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$color <- "grey"
            } else {
                nodes[[i]]@attrs$color <- bordercol[[nodes[[i]]@name]]
            }
            if (is.null(nodeshape[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$shape <- "box"
            } else {
                nodes[[i]]@attrs$shape <- nodeshape[[nodes[[i]]@name]]
            }
            if (is.null(nodecol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$fillcolor <- "grey"
            } else {
                nodes[[i]]@attrs$fillcolor <- nodecol[[nodes[[i]]@name]]
            }
            nodes[[i]]@attrs$width <- "0.5"
            nodes[[i]]@attrs$height <- "0.5"
            if (type == 2) {
                nodes[[i]]@attrs$fontsize <- "0"
            } else {
                nodes[[i]]@attrs$fontsize <- "0"
            }
        } else {
            nodes[[i]]@attrs$fontsize <- as.character(fontsize)
            if (is.null(nodecol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$fillcolor <- "white"
            } else {
                nodes[[i]]@attrs$fillcolor <- nodecol[[nodes[[i]]@name]]
            }
            if (is.null(nodeshape[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$shape <- "ellipse"
            } else {
                nodes[[i]]@attrs$shape <- nodeshape[[nodes[[i]]@name]]
            }
            if (is.null(bordercol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$color <- "black"
            } else {
                nodes[[i]]@attrs$color <- bordercol[[nodes[[i]]@name]]
            }
            if (names(nodes)[i] %in% stimuli & is.null(nodeshape[[nodes[[i]]@name]])) {
                if (type == 2) {
                    nodes[[i]]@attrs$shape <- "diamond"
                } else {
                    nodes[[i]]@attrs$shape <- "box"
                }
            }
            if (names(nodes)[i] %in% signals & is.null(nodecol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$fillcolor <- "lightblue"
            }
            if (names(nodes)[i] %in% inhibitors & is.null(bordercol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$color <- "red"
            }
        }
        if (!is.null(nodestates)) {
            if (sum(names(nodestates) %in% nodes[[i]]@name) == 1) {
                if (nodestates[which(names(nodestates) %in% nodes[[i]]@name)] == 0) {
                    if (is.null(nodecol[[nodes[[i]]@name]])) {
                        nodes[[i]]@attrs$fillcolor <- "white"
                    }
                    if (is.null(bordercol[[nodes[[i]]@name]])) {
                        nodes[[i]]@attrs$color <- "black"
                    }
                }
                if (nodestates[which(names(nodestates) %in% nodes[[i]]@name)] == 1) {
                    if (is.null(nodecol[[nodes[[i]]@name]])) {
                        nodes[[i]]@attrs$fillcolor <- "green"
                    }
                    if (is.null(bordercol[[nodes[[i]]@name]])) {
                        nodes[[i]]@attrs$color <- "black"
                    }
                }
            }
        }
    }
    if (length(edges) > 0) {
        for (i in 1:length(edges)) {
            edges[[i]]@attrs$fontsize <- as.character(labelsize)
            if (length(grep("and", names(edges)[i])) > 0) {
                tmp <- unlist(strsplit(names(edges)[i], "~"))
                k <- grep("and", tmp)
                inputN <- length(grep(tmp[k], edges))
                k <- as.numeric(gsub("and", "", tmp[k]))
                ## try to get the index of the correct edgecol:
                k2 <- grep("\\+", dnf)[k]
                if (grep("and", tmp) == 2) {
                    inputN2 <- which(gsub("!", "", unlist(strsplit(gsub("=.*", "", dnf[k2]), "\\+"))) %in% tmp[1])
                } else {
                    inputN2 <- length(unlist(strsplit(gsub("=.*", "", dnf[k2]), "\\+"))) + 1
                }
                if (k2 == 1) {
                    edgecolindex <- inputN2
                } else {
                    if (length(grep("\\+", graph[1:(k2-1)])) == 0) {
                        edgecolindex <- length(graph[1:(k2-1)]) + inputN2
                    } else {
                        edgecolindex <- length(unlist(strsplit(dnf[1:(k2-1)], "\\+"))) + length(grep("\\+", dnf[1:(k2-1)])) + inputN2
                    }
                }
                ## end
                inputN2 <- grep(tmp[1], unlist(strsplit(dnf[grep("\\+", dnf)[k]], "\\+")))-1
                edges[[i]]@attrs$style <- lty[grep("\\+", dnf)[k]]
                edges[[i]]@attrs$label <- labels[grep("\\+", dnf)[k]]
                if (use.freq) {
                    edges[[i]]@attrs$weight <- freq[grep("\\+", dnf)[k]]
                    edges[[i]]@attrs$fontcolor <- "blue"
                }
                if (!is.null(edgewidth)) {
                    edges[[i]]@attrs$weight <- edgewidth[grep("\\+", dnf)[k]]
                }
                if (!is.null(edgestyle)) {
                    if (!is.na(edgestyle[grep("\\+", dnf)[k]])) {
                        edges[[i]]@attrs$style <- edgestyle[edgecolindex] # [grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
                    }
                }
                if (!is.null(edgelabel)) {
                    if (!is.na(edgelabel[grep("\\+", dnf)[k]])) {
                        edges[[i]]@attrs$label <- edgelabel[edgecolindex] # [grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
                    }
                }
                if (length(grep("!", names(edgesneg)[i])) > 0) {
                    edges[[i]]@attrs$arrowhead <- "tee"
                    edges[[i]]@attrs$color <- "red"
                    if (!is.null(edgecol)) {
                        if (!is.na(edgecol[grep("\\+", dnf)[k]])) {
                            edges[[i]]@attrs$color <- edgecol[edgecolindex] # [grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
                        }
                    }
                    if (!is.null(edgehead)) {
                        if (!is.na(edgehead[grep("\\+", dnf)[k]])) {
                            edges[[i]]@attrs$arrowhead <- edgehead[edgecolindex] # [grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
                        }
                    }
                    if (!is.null(edgetail)) {
                        if (!is.na(edgetail[grep("\\+", dnf)[k]])) {
                            edges[[i]]@attrs$arrowtail <- edgetail[edgecolindex] # [grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
                        }
                    }
                } else {
                    edges[[i]]@attrs$arrowhead <- "open"
                    edges[[i]]@attrs$color <- "black"
                    if (gsub("and.*", "and", tmp[1]) %in% "and") {
                        if (!is.null(edgecol)) {
                            if (!is.na(edgecol[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$color <- edgecol[edgecolindex] # [grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k])])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k])])) > 1)]
                            }
                        }
                        if (!is.null(edgehead)) {
                            if (!is.na(edgehead[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$arrowhead <- edgehead[edgecolindex] # [grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k])])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k])])) > 1)]
                            }
                        }
                        if (!is.null(edgetail)) {
                            if (!is.na(edgetail[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$arrowtail <- edgetail[edgecolindex] # [grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k])])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k])])) > 1)]
                            }
                        }
                    } else {
                        if (!is.null(edgecol)) {
                            if (!is.na(edgecol[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$color <- edgecol[edgecolindex] # [grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
                            }
                        }
                        if (!is.null(edgehead)) {
                            if (!is.na(edgehead[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$arrowhead <- edgehead[edgecolindex] # [grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
                            }
                        }
                        if (!is.null(edgetail)) {
                            if (!is.na(edgetail[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$arrowtail <- edgetail[edgecolindex] # [grep("\\+", dnf)[k]+inputN+sum(unlist(gregexpr("\\+", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep("\\+", dnf)[k]-1)])) > 1)+inputN2]
                            }
                        }
                    }
                }
            } else {
                tmp <- unlist(strsplit(names(edges)[i], "~"))
                ## try to get the index of the correct edgecol:
                if (length(grep("!", names(edgesneg)[i])) == 0) {
                    k2 <- grep(paste("^", tmp[1], "=", tmp[2], sep = ""), dnf)
                } else {
                    k2 <- grep(paste("^!", tmp[1], "=", tmp[2], sep = ""), dnf)
                }
                if (k2 == 1) {
                    edgecolindex <- k2
                } else {
                    if (length(grep("\\+", dnf[1:(k2-1)])) == 0) {
                        edgecolindex <- k2
                    } else {
                        edgecolindex <- length(unlist(strsplit(dnf[1:(k2-1)], "\\+"))) + length(grep("\\+", dnf[1:(k2-1)])) + 1
                    }
                }
                ## end
                if (length(grep("!", names(edgesneg)[i])) > 0) {
                    edges[[i]]@attrs$style <- lty[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                    edges[[i]]@attrs$label <- labels[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                    if (use.freq) {
                        edges[[i]]@attrs$weight <- freq[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                        edges[[i]]@attrs$fontcolor <- "blue"
                    }
                    if (!is.null(edgewidth)) {
                        edges[[i]]@attrs$weight <- edgewidth[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                    }
                    if (!is.null(edgestyle)) {
                        if (!is.na(edgestyle[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$style <- edgestyle[edgecolindex] # [grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
                        }
                    }
                    if (!is.null(edgelabel)) {
                        if (!is.na(edgelabel[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$label <- edgelabel[edgecolindex] # [grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
                        }
                    }
                    edges[[i]]@attrs$arrowhead <- "tee"
                    edges[[i]]@attrs$color <- "red"
                    if (!is.null(edgecol)) {
                        if (!is.na(edgecol[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$color <- edgecol[edgecolindex] # [grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
                        }
                    }
                    if (!is.null(edgehead)) {
                        if (!is.na(edgehead[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$arrowhead <- edgehead[edgecolindex] # [grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
                        }
                    }
                    if (!is.null(edgetail)) {
                        if (!is.na(edgetail[grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$arrowtail <- edgetail[edgecolindex] # [grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
                        }
                    }
                } else {
                    edges[[i]]@attrs$style <- lty[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                    edges[[i]]@attrs$label <- labels[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                    if (use.freq) {
                        edges[[i]]@attrs$weight <- freq[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                        edges[[i]]@attrs$fontcolor <- "blue"
                    }
                    if (!is.null(edgewidth)) {
                        edges[[i]]@attrs$weight <- edgewidth[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)]
                    }
                    if (!is.null(edgestyle)) {
                        if (!is.na(edgestyle[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$style <- edgestyle[edgecolindex] # [grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
                        }
                    }
                    if (!is.null(edgelabel)) {
                        if (!is.na(edgelabel[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$label <- edgelabel[edgecolindex] # [grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
                        }
                    }
                    edges[[i]]@attrs$arrowhead <- "open"
                    edges[[i]]@attrs$color <- "black"
                    if (!is.null(edgecol)) {
                        if (!is.na(edgecol[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$color <- edgecol[edgecolindex] # [grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
                        }
                    }
                    if (!is.null(edgehead)) {
                        if (!is.na(edgehead[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$arrowhead <- edgehead[edgecolindex] # [grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
                        }
                    }
                    if (!is.null(edgetail)) {
                        if (!is.na(edgetail[grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)])) {
                            edges[[i]]@attrs$arrowtail <- edgetail[edgecolindex] # [grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)+sum(unlist(gregexpr("\\+", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)+sum(unlist(gregexpr("\\+.*=.*", graph[1:(grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""), dnf)-1)])) > 1)]
                        }
                    }
                }
            }
        }
    }
    if (type == 1) {
        
        g2 <- agopen(name="boolean", nodes=nodes, recipEdges = "distinct", edges=edges, edgeMode="undirected", attrs=list(edge = list(), graph = list(lwd = lwd, rankdir = rankdir), node=list(lwd = lwd, fixedsize=FALSE)))

        plot(g2, "dot", lwd = lwd, ...)

    } else {

        arrowheads <- character()
        arrowtails <- character()
        arrowcolors <- character()
        arrowlabels <- character()
        arrowlwd <- character()
        arrowlty <- character()
        arrowfontsize <- character()
        if (length(edges) > 0) {
            for (i in 1:length(edges)) {
                if (length(edges[[i]]@attrs$style) == 0) { edges[[i]]@attrs$style <- "solid" }
                arrowlty <- c(arrowlty, edges[[i]]@attrs$style)
                arrowheads <- c(arrowheads, edges[[i]]@attrs$arrowhead)
                if (!is.null(edgetail)) {
                    arrowtails <- c(arrowtails, edges[[i]]@attrs$arrowtail)
                } else {
                    arrowtails <- c(arrowtails, "none")
                }
                arrowcolors <- c(arrowcolors, edges[[i]]@attrs$color)
                arrowfontsize <- c(arrowfontsize, edges[[i]]@attrs$fontsize)
                arrowlwd <- c(arrowlwd, edges[[i]]@attrs$weight)
                arrowlabels <- c(arrowlabels, edges[[i]]@attrs$label)
            }
        }
        arrowlwd <- as.numeric(arrowlwd)
        
        graph.trans <- NULL
        and.count <- 0
        for (i in dnf) {
            if (length(grep("\\+", i)) > 0) {
                and.count <- and.count + 1
                output <- unlist(strsplit(i, "="))
                input <- unlist(strsplit(output[1], "\\+"))
                output <- output[2]
                for (i in input) {
                    graph.trans <- c(graph.trans, paste(gsub("!", "", i), "~", paste("and", and.count, sep = ""), sep = ""))
                }
                graph.trans <- c(graph.trans, paste(paste("and", and.count, sep = ""), "~", output, sep = ""))
            } else {
                put <- unlist(strsplit(i, "="))
                graph.trans <- c(graph.trans, paste(gsub("!", "", put[1]), "~", put[2], sep = ""))
            }
        }

        if (length(edgecol) == length(arrowcolors)) {
            edgecol <- edgecol[order(match(graph.trans, names(edges)))]
            arrowcolors <- edgecol
        }
        
        nodeshapes <- character()
        nodecolors <- character()
        nodeheight <- character()
        nodewidth <- character()
        nodecolor <- character()
        
        for (i in 1:length(nodes)) {
            nodeshapes <- c(nodeshapes, nodes[[i]]@attrs$shape)
            nodecolors <- c(nodecolors, nodes[[i]]@attrs$fillcolor)
            nodeheight <- c(nodeheight, nodes[[i]]@attrs$height)
            nodewidth <- c(nodewidth, nodes[[i]]@attrs$width)
            nodecolor <- c(nodecolor, nodes[[i]]@attrs$color)
        }

        nodeheight[which(nodeheight == "0.4")] <- "0.2"

        if (is.null(lty) & is.null(edgestyle)) {
            arrowlty <- rep("solid", length(edges))
        }
        
        if (use.freq) {
            if (is.null(lty)) {
                if (edgestyle) {
                    arrowlty[which(as.numeric(arrowlabels) < 66)] <- "dashed"
                    arrowlty[which(as.numeric(arrowlabels) < 33)] <- "dotted"
                }
            }
            arrowlwd <- arrowlwd - min(arrowlwd)
            arrowlwd <- as.character((arrowlwd/max(arrowlwd)+0.1)*2*edgelwd)
        } else {
            if (is.null(edgewidth)) {
                arrowlwd <- rep(edgelwd, length(edges))
            }
        }

        if (is.null(edgewidth) & is.null(edgelabel)) {
            arrowlabels <- rep("", length(edges))
        }
        
        if (length(arrowlty) == 0) {
            arrowlty <- rep("solid", length(edges))
        }
        if (length(arrowlwd) == 0) {
            arrowlwd <- rep(lwd, length(edges))
        }
        
        names(arrowfontsize) <- names(arrowheads) <- names(arrowtails) <- names(arrowcolors) <- names(arrowlwd) <- names(arrowlty) <- names(arrowlabels) <- names(edges)

        names(nodecolor) <- names(nodewidth) <- names(nodeheight) <- names(nodeshapes) <- names(nodecolors) <- names(nodes)

        if (length(unique(names(edges))) < length(names(edges))) {
            for (i in names(edges)[-which(duplicated(names(edges)) == TRUE)]) {
                getpos <- grep(paste("^", i, "$", sep = ""), names(edges))
                if (length(getpos) > 1) {
                    if (use.freq) {
                        if (arrowheads[getpos[1]] %in% "tee") {
                            arrowlabels[getpos[1]] <- paste(paste(c("-", "+"), arrowlabels[getpos], sep = ""), collapse = "\n")
                        } else {
                            arrowlabels[getpos[1]] <- paste(paste(c("+", "-"), arrowlabels[getpos], sep = ""), collapse = "\n")
                        }
                    } else {
                        if (is.null(edgelabel)) {
                            arrowlabels[getpos[1]] <- ""
                        }
                    }
                    arrowheads[getpos[1]] <- "odiamond"
                    if (is.null(edgecol)) {
                        arrowcolors[getpos[1]] <- "black"
                    } else {
                        if (is.na(edgecol[getpos[1]])) {
                            arrowcolors[getpos[1]] <- "black"
                        }
                    }
                    arrowlwd[getpos[1]] <- as.character(mean(as.numeric(arrowlwd[getpos])))
                }
            }
        }
        if (length(labels) == length(arrowlabels) & is.null(edgelabel)) {
            arrowlabels[!is.na(labels)] <- labels[!is.na(labels)]
        }
        if (length(edgecol) == 1) {
            arrowcolors <- rep(edgecol, length(arrowcolors))
            names(arrowcolors) <- names(arrowlabels)
        }
        
        if (legend == 1 | legend == 3) {
            if (dolegend) {
                start <- 1
                g@nodes <- c("LEGEND:", "STIMULUS", "INHIBITOR", "SIGNAL", "NOTHING", "active", "inactive")
                g@edgeL <- list()
                g@edgeData@data <- list()
            } else {
                start <- length(g@nodes) + 1
                g@nodes <- c(g@nodes, "LEGEND:", "STIMULUS", "INHIBITOR", "SIGNAL", "NOTHING", "active", "inactive")
            }
            g@edgeL[["LEGEND:"]] <- list()
            g@edgeL[["STIMULUS"]] <- list()
            g@edgeL[["INHIBITOR"]] <- list()
            g@edgeL[["SIGNAL"]] <- list()
            g@edgeL[["NOTHING"]] <- list()
            g@edgeL[["active"]] <- list()
            g@edgeL[["inactive"]] <- list()
            g@edgeL[["LEGEND:"]][["edges"]] <- as.integer(start+1)
            g@edgeL[["STIMULUS"]][["edges"]] <- as.integer(start+2)
            g@edgeL[["INHIBITOR"]][["edges"]] <- as.integer(start+3)
            g@edgeL[["SIGNAL"]][["edges"]] <- as.integer(start+4)
            g@edgeL[["NOTHING"]][["edges"]] <- as.integer(start+5)
            g@edgeL[["active"]][["edges"]] <- c(as.integer(start+6), as.integer(start+4))
            g@edgeL[["inactive"]][["edges"]] <- as.integer(start+5)
            g@edgeData@data[["LEGEND:|STIMULUS"]] <- list()
            g@edgeData@data[["STIMULUS|INHIBITOR"]] <- list()
            g@edgeData@data[["INHIBITOR|SIGNAL"]] <- list()
            g@edgeData@data[["SIGNAL|NOTHING"]] <- list()
            g@edgeData@data[["NOTHING|active"]] <- list()
            g@edgeData@data[["active|inactive"]] <- list()
            g@edgeData@data[["active|NOTHING"]] <- list()
            g@edgeData@data[["inactive|active"]] <- list()
            g@edgeData@data[["LEGEND:|STIMULUS"]][["weight"]] <- 1
            g@edgeData@data[["STIMULUS|INHIBITOR"]][["weight"]] <- 1
            g@edgeData@data[["INHIBITOR|SIGNAL"]][["weight"]] <- 1
            g@edgeData@data[["SIGNAL|NOTHING"]][["weight"]] <- 1
            g@edgeData@data[["NOTHING|active"]][["weight"]] <- 1
            g@edgeData@data[["active|inactive"]][["weight"]] <- 1
            g@edgeData@data[["active|NOTHING"]][["weight"]] <- 1
            g@edgeData@data[["inactive|active"]][["weight"]] <- 1
            arrowheads <- c(arrowheads, "LEGEND:~STIMULUS" = "none", "STIMULUS~INHIBITOR" = "open", "INHIBITOR~SIGNAL" = "tee", "SIGNAL~NOTHING" = "odiamond", "NOTHING~active" = "open", "active~inactive" = "tee", "active~NOTHING" = "tee", "inactive~active" = "open")
            arrowcolors <- c(arrowcolors, "LEGEND:~STIMULUS" = "transparent", "STIMULUS~INHIBITOR" = "black", "INHIBITOR~SIGNAL" = "red", "SIGNAL~NOTHING" = "blue", "NOTHING~active" = "black", "active~inactive" = "red", "active~NOTHING" = "red", "inactive~active" = "black")
            arrowlabels <- c(arrowlabels, "LEGEND:~STIMULUS" = "", "STIMULUS~INHIBITOR" = "    positive", "INHIBITOR~SIGNAL" = "    negative", "SIGNAL~NOTHING" = "    ambiguous\npositive\nnegative", "NOTHING~active" = "    bidirectional\ndifferent", "active~inactive" = "    bidirectional\ndifferent", "active~NOTHING" = "", "inactive~active" = "")
            nodecolors <- c(nodecolors, "LEGEND:" = "white", "STIMULUS" = "white", "INHIBITOR" = "white", "SIGNAL" = "lightblue", "NOTHING" = "white", "active" = "green", "inactive" = "white")
            nodeheight <- c(nodeheight, "LEGEND:" = 0, "STIMULUS" = as.character(max(nodeheight)), "INHIBITOR" = as.character(max(nodeheight)), "SIGNAL" = as.character(max(nodeheight)), "NOTHING" = as.character(max(nodeheight)), "active" = as.character(max(nodeheight)), "inactive" = as.character(max(nodeheight)))
            nodewidth <- c(nodewidth, "LEGEND:" = as.character(max(nodewidth)), "STIMULUS" = as.character(max(nodewidth)), "INHIBITOR" = as.character(max(nodewidth)), "SIGNAL" = as.character(max(nodewidth)), "NOTHING" = as.character(max(nodewidth)), "active" = as.character(max(nodewidth)), "inactive" = as.character(max(nodewidth)))
            if (type == 2) {
                nodeshapes <- c(nodeshapes, "LEGEND:" = "box", "STIMULUS" = "diamond", "INHIBITOR" = "ellipse", "SIGNAL" = "ellipse", "NOTHING" = "ellipse", "active" = "ellipse", "inactive" = "ellipse")
            } else {
                nodeshapes <- c(nodeshapes, "LEGEND:" = "box", "STIMULUS" = "box", "INHIBITOR" = "ellipse", "SIGNAL" = "ellipse", "NOTHING" = "ellipse", "active" = "ellipse", "inactive" = "ellipse")
            }
            nodecolor <- c(nodecolor, "LEGEND:" = "white", "STIMULUS" = "black", "INHIBITOR" = "red", "SIGNAL" = "black", "NOTHING" = "black", "active" = "black", "inactive" = "black")
            dnf <- c(dnf, "NOTHING=active", "!active=NOTHING", "!active=inactive", "inactive=active")
        }
        nodelabels <- names(nodecolor)
        names(nodelabels) <- nodelabels
        for (i in 1:length(nodelabel)) {
            nodelabels[which(names(nodelabels) %in% names(nodelabel)[i])] <- nodelabel[i]
        }
        nodefontsizes <- NULL
        if (!is.null(nodefontsize)) {
            nodefontsizes <- rep(14, length(nodelabels))
            names(nodefontsizes) <- names(nodelabels)
            for (i in 1:length(nodefontsize)) {
                nodefontsizes[which(names(nodefontsizes) %in% names(nodefontsize)[i])] <- nodefontsize[[i]]
            }
        }
        g <- layoutGraph(g, edgeAttrs = list(arrowhead = arrowheads, color = arrowcolors, label = arrowlabels, arrowtail = arrowtails), nodeAttrs = list(labels = nodelabels, color = nodecolor, height = nodeheight, width = nodewidth, shape = nodeshapes, fillcolor = nodecolors), layoutType=layout)
        graph.par(list(graph=list(main = main, sub = sub, cex.main = cex.main, cex.sub = cex.sub, col.sub = col.sub), edges=list(textCol = labelcol, lwd = edgelwd, fontsize = labelsize), nodes=list(lwd = lwd, fontsize = fontsize, cex = cex)))
        edgeRenderInfo(g) <- list(lty = arrowlty, lwd = arrowlwd, label = arrowlabels)
        if (length(edges) > 0) {
            for (i in names(g@renderInfo@edges$direction)) {
                input <- unlist(strsplit(i, "~"))
                output <- input[2]
                input <- input[1]
                ambig <- FALSE
                if (paste("!", input, "=", output, sep = "") %in% dnf & paste("", input, "=", output, sep = "") %in% dnf) {
                    ambig <- TRUE
                }
                if ((length(grep("and", i)) == 0 & g@renderInfo@edges$direction[[i]] == "both") | ambig) {
                    pos <- which(names(g@renderInfo@edges$arrowhead) %in% i)
                    if (is.null(edgehead)) {
                        if (paste("!", input, "=", output, sep = "") %in% dnf) {
                            g@renderInfo@edges$arrowhead[pos] <- "tee"
                        }
                        if (paste(input, "=", output, sep = "") %in% dnf) {
                            g@renderInfo@edges$arrowhead[pos] <- "open"
                        }
                        if (paste("!", output, "=", input, sep = "") %in% dnf) {
                            g@renderInfo@edges$arrowtail[pos] <- "tee"
                        }
                        if (paste(output, "=", input, sep = "") %in% dnf) {
                            g@renderInfo@edges$arrowtail[pos] <- "open"
                        }
                        if (paste("!", output, "=", input, sep = "") %in% dnf & paste("", output, "=", input, sep = "") %in% dnf) {
                            g@renderInfo@edges$arrowtail[pos] <- "odiamond"
                        }
                        if (paste("!", input, "=", output, sep = "") %in% dnf & paste("", input, "=", output, sep = "") %in% dnf) {
                            ## for (f in 1:length(g@renderInfo@edges)) {
                            ##   g@renderInfo@edges[[f]] <- c(g@renderInfo@edges[[f]], g@renderInfo@edges[[f]][pos])
                            ##   names(g@renderInfo@edges[[f]])[length(g@renderInfo@edges[[f]])] <- paste(input, "recip=", output, "recip", sep = "")
                            ## }
                            ## g@renderInfo@edges$splines[[length(g@renderInfo@edges$enamesFrom)]][[1]]@cPoints[[1]]@x <- as.integer(g@renderInfo@edges$splines[[length(g@renderInfo@edges$enamesFrom)]][[1]]@cPoints[[1]]@x + 20)
                            ## g@renderInfo@edges$splines[[length(g@renderInfo@edges$enamesFrom)]][[1]]@cPoints[[2]]@x <- as.integer(g@renderInfo@edges$splines[[length(g@renderInfo@edges$enamesFrom)]][[1]]@cPoints[[2]]@x + 40)
                            ## g@renderInfo@edges$splines[[length(g@renderInfo@edges$enamesFrom)]][[1]]@cPoints[[3]]@x <- as.integer(g@renderInfo@edges$splines[[length(g@renderInfo@edges$enamesFrom)]][[1]]@cPoints[[3]]@x + 40)
                            ## g@renderInfo@edges$splines[[length(g@renderInfo@edges$enamesFrom)]][[1]]@cPoints[[4]]@x <- as.integer(g@renderInfo@edges$splines[[length(g@renderInfo@edges$enamesFrom)]][[1]]@cPoints[[4]]@x + 20)
                            ## g@renderInfo@edges$enamesFrom[length(g@renderInfo@edges$enamesFrom)] <- paste(input, "recip", sep = "")
                            ## g@renderInfo@edges$enamesTo[length(g@renderInfo@edges$enamesTo)] <- paste(output, "recip", sep = "")
                            ## g@renderInfo@edges$arrowhead[length(g@renderInfo@edges$arrowhead)] <- "open"
                            ## g@renderInfo@edges$arrowhead[pos] <- "tee"
                            ## g@nodes <- c(g@nodes, paste(input, "recip", sep = ""), paste(output, "recip", sep = ""))
                            ## for (f in 1:length(g@renderInfo@nodes)) {
                            ##   g@renderInfo@nodes[[f]] <- c(g@renderInfo@nodes[[f]], g@renderInfo@nodes[[f]][1], g@renderInfo@nodes[[f]][1])
                            ##   names(g@renderInfo@nodes[[f]])[(length(g@renderInfo@nodes[[f]])-1):length(g@renderInfo@nodes[[f]])] <- c(paste(input, "recip", sep = ""), paste(output, "recip", sep = ""))
                            ## }
                            ## g@renderInfo@nodes$col[(length(g@renderInfo@nodes$col)-1):length(g@renderInfo@nodes$col)] <- "transparent"
                            ## g@renderInfo@nodes$fill[(length(g@renderInfo@nodes$fill)-1):length(g@renderInfo@nodes$fill)] <- "transparent"
                            ## g@renderInfo@nodes$label[(length(g@renderInfo@nodes$label)-1):length(g@renderInfo@nodes$label)] <- ""
                            ## g@edgeL[[paste(input, "recip", sep = "")]]$edges <- length(g@edgeL)+1
                            ## g@edgeL[[paste(output, "recip", sep = "")]]$edges <- numeric()
                            ## g@edgeData@data[[paste(input, "recip|", output, "recip", sep = "")]]$weight <- 1
                            g@renderInfo@edges$arrowhead[pos] <- "odiamond"
                        }
                    }
                    if (is.null(edgecol)) {
                        if (g@renderInfo@edges$arrowtail[pos] == "open" & g@renderInfo@edges$arrowhead[pos] == "open") {
                            g@renderInfo@edges$col[pos] <- "black"
                        }
                        if (g@renderInfo@edges$arrowtail[pos] == "tee" & g@renderInfo@edges$arrowhead[pos] == "tee") {
                            g@renderInfo@edges$col[pos] <- "red"
                        }
                        if (g@renderInfo@edges$arrowtail[pos] != g@renderInfo@edges$arrowhead[pos]) {
                            g@renderInfo@edges$col[pos] <- "brown"
                        }
                        if (g@renderInfo@edges$arrowtail[pos] == "odiamond" | g@renderInfo@edges$arrowhead[pos] == "odiamond") {
                            g@renderInfo@edges$col[pos] <- "blue"
                        }
                    } else {
                        if (is.null(edgecol)) { # is.na(edgecol[pos])
                            if (g@renderInfo@edges$arrowtail[pos] == "open" & g@renderInfo@edges$arrowhead[pos] == "open") {
                                g@renderInfo@edges$col[pos] <- "black"
                            }
                            if (g@renderInfo@edges$arrowtail[pos] == "tee" & g@renderInfo@edges$arrowhead[pos] == "tee") {
                                g@renderInfo@edges$col[pos] <- "red"
                            }
                            if (g@renderInfo@edges$arrowtail[pos] != g@renderInfo@edges$arrowhead[pos]) {
                                g@renderInfo@edges$col[pos] <- "brown"
                            }
                            if (g@renderInfo@edges$arrowtail[pos] == "odiamond" | g@renderInfo@edges$arrowhead[pos] == "odiamond") {
                                g@renderInfo@edges$col[pos] <- "blue"
                            }
                        }
                    }
                }
            }
        }
                                        #g@renderInfo@nodes$labelX[grep("and", names(g@renderInfo@nodes$labelX))] <- -1000
                                        #g@renderInfo@nodes$labelY[grep("and", names(g@renderInfo@nodes$labelY))] <- -1000
        if (!is.null(simulate$draw)) {
            for (i in simulate$inhibitors) {
                ## add the inhibiting node
                g@nodes <- c(g@nodes, paste(i, "_inhibited", sep = ""))
                g@renderInfo@nodes$nodeX <- c(g@renderInfo@nodes$nodeX, g@renderInfo@nodes$nodeX[which(names(g@renderInfo@nodes$nodeX) %in% i)])
                g@renderInfo@nodes$nodeY <- c(g@renderInfo@nodes$nodeY, g@renderInfo@nodes$nodeY[which(names(g@renderInfo@nodes$nodeY) %in% i)])
                names(g@renderInfo@nodes$nodeX)[length(g@renderInfo@nodes$nodeX)] <- paste(i, "_inhibited", sep = "")
                names(g@renderInfo@nodes$nodeY)[length(g@renderInfo@nodes$nodeY)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$labelX <- c(g@renderInfo@nodes$labelX, g@renderInfo@nodes$labelX[which(names(g@renderInfo@nodes$labelX) %in% i)] + 200)
                g@renderInfo@nodes$labelY <- c(g@renderInfo@nodes$labelY, g@renderInfo@nodes$labelY[which(names(g@renderInfo@nodes$labelY) %in% i)])
                names(g@renderInfo@nodes$labelX)[length(g@renderInfo@nodes$labelX)] <- paste(i, "_inhibited", sep = "")
                names(g@renderInfo@nodes$labelY)[length(g@renderInfo@nodes$labelY)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$labelJust <- c(g@renderInfo@nodes$labelJust, "n")
                names(g@renderInfo@nodes$labelJust)[length(g@renderInfo@nodes$labelJust)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$col <- c(g@renderInfo@nodes$col, "transparent")
                names(g@renderInfo@nodes$col)[length(g@renderInfo@nodes$col)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$fill <- c(g@renderInfo@nodes$fill, "transparent")
                names(g@renderInfo@nodes$fill)[length(g@renderInfo@nodes$fill)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$shape <- c(g@renderInfo@nodes$shape, "box")
                names(g@renderInfo@nodes$shape)[length(g@renderInfo@nodes$shape)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$style <- c(g@renderInfo@nodes$style, "")
                names(g@renderInfo@nodes$style)[length(g@renderInfo@nodes$style)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$height <- c(g@renderInfo@nodes$height, 2)
                names(g@renderInfo@nodes$height)[length(g@renderInfo@nodes$height)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$rWidth <- c(g@renderInfo@nodes$rWidth, 1)
                names(g@renderInfo@nodes$rWidth)[length(g@renderInfo@nodes$rWidth)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$lWidth <- c(g@renderInfo@nodes$lWidth, 1)
                names(g@renderInfo@nodes$lWidth)[length(g@renderInfo@nodes$lWidth)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$label[[paste(i, "_inhibited", sep = "")]] <- ""#"X            X\n X       X\n  X  X\n   X\n  X  X\n X       X\nX            X"
                g@renderInfo@nodes$labelWidth <- c(g@renderInfo@nodes$labelWidth, g@renderInfo@nodes$labelWidth[which(names(g@renderInfo@nodes$labelWidth) %in% i)])
                names(g@renderInfo@nodes$labelWidth)[length(g@renderInfo@nodes$labelWidth)] <- paste(i, "_inhibited", sep = "")

                ## add the inhibiting edge
                tmp.name <- paste(i, "_inhibited", sep = "")
                g@edgeL[[tmp.name]] <- list()
                g@edgeL[[tmp.name]][["edges"]] <- which(g@nodes %in% i)
                g@edgeData@data[[paste(tmp.name, "|", i, sep = "")]] <- list()
                g@edgeData@data[[paste(tmp.name, "|", i, sep = "")]]$weight <- 1
                g@renderInfo@edges$enamesFrom <- c(g@renderInfo@edges$enamesFrom, tmp.name)
                names(g@renderInfo@edges$enamesFrom)[length(g@renderInfo@edges$enamesFrom)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$enamesTo <- c(g@renderInfo@edges$enamesTo, i)
                names(g@renderInfo@edges$enamesTo)[length(g@renderInfo@edges$enamesTo)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelJust <- c(g@renderInfo@edges$labelJust, NA)
                names(g@renderInfo@edges$labelJust)[length(g@renderInfo@edges$labelJust)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelX <- c(g@renderInfo@edges$labelX, NA)
                names(g@renderInfo@edges$labelX)[length(g@renderInfo@edges$labelX)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelY <- c(g@renderInfo@edges$labelY, NA)
                names(g@renderInfo@edges$labelY)[length(g@renderInfo@edges$labelY)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelWidth <- c(g@renderInfo@edges$labelWidth, NA)
                names(g@renderInfo@edges$labelWidth)[length(g@renderInfo@edges$labelWidth)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$label <- c(g@renderInfo@edges$label, "")
                names(g@renderInfo@edges$label)[length(g@renderInfo@edges$label)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$arrowhead <- c(g@renderInfo@edges$arrowhead, "tee")
                names(g@renderInfo@edges$arrowhead)[length(g@renderInfo@edges$arrowhead)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$arrowtail <- c(g@renderInfo@edges$arrowtail, "odot")
                names(g@renderInfo@edges$arrowtail)[length(g@renderInfo@edges$arrowtail)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$col <- c(g@renderInfo@edges$col, "firebrick")
                names(g@renderInfo@edges$col)[length(g@renderInfo@edges$col)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$lwd <- c(g@renderInfo@edges$lwd, lwd[1]*1)
                names(g@renderInfo@edges$lwd)[length(g@renderInfo@edges$lwd)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$lty <- c(g@renderInfo@edges$lty, "solid")
                names(g@renderInfo@edges$lty)[length(g@renderInfo@edges$lty)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$direction <- c(g@renderInfo@edges$direction, "forward")
                names(g@renderInfo@edges$direction)[length(g@renderInfo@edges$direction)] <- paste(tmp.name, "~", i, sep = "")
                ## calculate splines
                tmp.splines <- rep(g@renderInfo@nodes$labelY[which(names(g@renderInfo@nodes$labelY) %in% i)], 8)
                tmp.splines[c(7,5,3,1)] <- round(seq(g@renderInfo@nodes$labelX[which(names(g@renderInfo@nodes$labelX) %in% i)] + g@renderInfo@nodes$rWidth[[i]] + 10,
                                                     g@renderInfo@nodes$labelX[which(names(g@renderInfo@nodes$labelX) %in% i)] + 200 - g@renderInfo@nodes$rWidth[[i]] - 10,
                                                     length.out = 4))
                ## tmp.splines[c(2,4,6,8)] <- seq(tmp.splines[2] + 50, tmp.splines[8], length.out = 4)
                tmp.splines[2] <- tmp.splines[2] + 50
                tmp.splines <- as.integer(tmp.splines)
                g@renderInfo@edges$splines[[paste(tmp.name, "~", i, sep = "")]] <- g@renderInfo@edges$splines[[1]]
                for (j in 1:4) {
                    g@renderInfo@edges$splines[[paste(tmp.name, "~", i, sep = "")]][[1]]@cPoints[[j]]@x <- tmp.splines[c(1,3,5,7)][j]
                    g@renderInfo@edges$splines[[paste(tmp.name, "~", i, sep = "")]][[1]]@cPoints[[j]]@y <- tmp.splines[c(2,4,6,8)][j]
                }
            }
            for (i in simulate$stimuli) {
                ## add the stimulating node
                g@nodes <- c(g@nodes, paste(i, "_inhibited", sep = ""))
                g@renderInfo@nodes$nodeX <- c(g@renderInfo@nodes$nodeX, g@renderInfo@nodes$nodeX[which(names(g@renderInfo@nodes$nodeX) %in% i)])
                g@renderInfo@nodes$nodeY <- c(g@renderInfo@nodes$nodeY, g@renderInfo@nodes$nodeY[which(names(g@renderInfo@nodes$nodeY) %in% i)])
                names(g@renderInfo@nodes$nodeX)[length(g@renderInfo@nodes$nodeX)] <- paste(i, "_inhibited", sep = "")
                names(g@renderInfo@nodes$nodeY)[length(g@renderInfo@nodes$nodeY)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$labelX <- c(g@renderInfo@nodes$labelX, g@renderInfo@nodes$labelX[which(names(g@renderInfo@nodes$labelX) %in% i)] + 200)
                g@renderInfo@nodes$labelY <- c(g@renderInfo@nodes$labelY, g@renderInfo@nodes$labelY[which(names(g@renderInfo@nodes$labelY) %in% i)])
                names(g@renderInfo@nodes$labelX)[length(g@renderInfo@nodes$labelX)] <- paste(i, "_inhibited", sep = "")
                names(g@renderInfo@nodes$labelY)[length(g@renderInfo@nodes$labelY)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$labelJust <- c(g@renderInfo@nodes$labelJust, "n")
                names(g@renderInfo@nodes$labelJust)[length(g@renderInfo@nodes$labelJust)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$col <- c(g@renderInfo@nodes$col, "transparent")
                names(g@renderInfo@nodes$col)[length(g@renderInfo@nodes$col)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$fill <- c(g@renderInfo@nodes$fill, "transparent")
                names(g@renderInfo@nodes$fill)[length(g@renderInfo@nodes$fill)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$shape <- c(g@renderInfo@nodes$shape, "box")
                names(g@renderInfo@nodes$shape)[length(g@renderInfo@nodes$shape)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$style <- c(g@renderInfo@nodes$style, "")
                names(g@renderInfo@nodes$style)[length(g@renderInfo@nodes$style)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$height <- c(g@renderInfo@nodes$height, 2)
                names(g@renderInfo@nodes$height)[length(g@renderInfo@nodes$height)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$rWidth <- c(g@renderInfo@nodes$rWidth, 1)
                names(g@renderInfo@nodes$rWidth)[length(g@renderInfo@nodes$rWidth)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$lWidth <- c(g@renderInfo@nodes$lWidth, 1)
                names(g@renderInfo@nodes$lWidth)[length(g@renderInfo@nodes$lWidth)] <- paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$label[[paste(i, "_inhibited", sep = "")]] <- ""#"X            X\n X       X\n  X  X\n   X\n  X  X\n X       X\nX            X"
                g@renderInfo@nodes$labelWidth <- c(g@renderInfo@nodes$labelWidth, g@renderInfo@nodes$labelWidth[which(names(g@renderInfo@nodes$labelWidth) %in% i)])
                names(g@renderInfo@nodes$labelWidth)[length(g@renderInfo@nodes$labelWidth)] <- paste(i, "_inhibited", sep = "")

                ## add the stimulating edge
                tmp.name <- paste(i, "_inhibited", sep = "")
                g@edgeL[[tmp.name]] <- list()
                g@edgeL[[tmp.name]][["edges"]] <- which(g@nodes %in% i)
                g@edgeData@data[[paste(tmp.name, "|", i, sep = "")]] <- list()
                g@edgeData@data[[paste(tmp.name, "|", i, sep = "")]]$weight <- 1
                g@renderInfo@edges$enamesFrom <- c(g@renderInfo@edges$enamesFrom, tmp.name)
                names(g@renderInfo@edges$enamesFrom)[length(g@renderInfo@edges$enamesFrom)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$enamesTo <- c(g@renderInfo@edges$enamesTo, i)
                names(g@renderInfo@edges$enamesTo)[length(g@renderInfo@edges$enamesTo)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelJust <- c(g@renderInfo@edges$labelJust, NA)
                names(g@renderInfo@edges$labelJust)[length(g@renderInfo@edges$labelJust)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelX <- c(g@renderInfo@edges$labelX, NA)
                names(g@renderInfo@edges$labelX)[length(g@renderInfo@edges$labelX)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelY <- c(g@renderInfo@edges$labelY, NA)
                names(g@renderInfo@edges$labelY)[length(g@renderInfo@edges$labelY)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelWidth <- c(g@renderInfo@edges$labelWidth, NA)
                names(g@renderInfo@edges$labelWidth)[length(g@renderInfo@edges$labelWidth)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$label <- c(g@renderInfo@edges$label, "")
                names(g@renderInfo@edges$label)[length(g@renderInfo@edges$label)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$arrowhead <- c(g@renderInfo@edges$arrowhead, "open")
                names(g@renderInfo@edges$arrowhead)[length(g@renderInfo@edges$arrowhead)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$arrowtail <- c(g@renderInfo@edges$arrowtail, "odot")
                names(g@renderInfo@edges$arrowtail)[length(g@renderInfo@edges$arrowtail)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$col <- c(g@renderInfo@edges$col, "limegreen")
                names(g@renderInfo@edges$col)[length(g@renderInfo@edges$col)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$lwd <- c(g@renderInfo@edges$lwd, lwd[1]*1)
                names(g@renderInfo@edges$lwd)[length(g@renderInfo@edges$lwd)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$lty <- c(g@renderInfo@edges$lty, "solid")
                names(g@renderInfo@edges$lty)[length(g@renderInfo@edges$lty)] <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$direction <- c(g@renderInfo@edges$direction, "forward")
                names(g@renderInfo@edges$direction)[length(g@renderInfo@edges$direction)] <- paste(tmp.name, "~", i, sep = "")
                ## calculate splines
                tmp.splines <- rep(g@renderInfo@nodes$labelY[which(names(g@renderInfo@nodes$labelY) %in% i)], 8)
                tmp.splines[c(7,5,3,1)] <- round(seq(g@renderInfo@nodes$labelX[which(names(g@renderInfo@nodes$labelX) %in% i)] + g@renderInfo@nodes$rWidth[[i]] + 10,
                                                     g@renderInfo@nodes$labelX[which(names(g@renderInfo@nodes$labelX) %in% i)] + 200 - g@renderInfo@nodes$rWidth[[i]] - 10,
                                                     length.out = 4))
                ## tmp.splines[c(2,4,6,8)] <- seq(tmp.splines[2] + 50, tmp.splines[8], length.out = 4)
                tmp.splines[2] <- tmp.splines[2] + 50
                tmp.splines <- as.integer(tmp.splines)
                g@renderInfo@edges$splines[[paste(tmp.name, "~", i, sep = "")]] <- g@renderInfo@edges$splines[[1]]
                for (j in 1:4) {
                    g@renderInfo@edges$splines[[paste(tmp.name, "~", i, sep = "")]][[1]]@cPoints[[j]]@x <- tmp.splines[c(1,3,5,7)][j]
                    g@renderInfo@edges$splines[[paste(tmp.name, "~", i, sep = "")]][[1]]@cPoints[[j]]@y <- tmp.splines[c(2,4,6,8)][j]
                }
            }
            g@renderInfo@graph$bbox[2,1] <- g@renderInfo@graph$bbox[2,1] + 150
            g@renderInfo@graph$bbox[2,2] <- g@renderInfo@graph$bbox[2,2] + 25
        }
        g <- g
        if (draw) {
            renderGraph(g, lwd = lwd, recipEdges = "distinct", ...)
        }
    }
    if (legend == 2 | legend == 3) {
        legend(x = x, y = y, legend = c("signals are blue", "stimuli are diamonds/boxes", "inhibitors have a red border", "positive regulation is green ->", "negative regulation is red -|", "ambiguous regulation is black -o"), fill = c("lightblue", "white", "red", "green", "red", "black"), col = c("lightblue", "white", "red", "green", "red", "black"), yjust = yjust, xjust = xjust)
    }
    return(g)
}
