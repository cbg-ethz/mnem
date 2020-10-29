#' Implementation of the original NEM
#'
#' Infers a signalling pathway from peerturbation experiments.
#' @param D data matrix with observed genes as rows and knock-down
#' experiments as columns
#' @param search either "greedy", "modules" or "exhaustive" (not
#' recommended for more than five S-genes)
#' @param start either NULL ("null") or a specific network to start
#' the greedy
#' @param method "llr" for log odds or p-values densities or "disc"
#' for binary data
#' @param parallel NULL for no parallel optimization or an integer
#' for the number of threads
#' @param reduce reduce search space (TRUE) for exhaustive search
#' @param weights a numeric vector of weights for the columns of D
#' @param runs the number of runs for the greedy search
#' @param verbose for verbose output (TRUE)
#' @param redSpace reduced search space for exhaustive search;
#' see result of exhaustive
#' search with reduce = TRUE
#' @param trans.close if TRUE uses the transitive closure of adj
#' @param subtopo optional matrix with the subtopology theta as
#' adjacency matrix
#' @param prior a prior network matrix for adj
#' @param ratio if FALSE uses alternative distance for the model score
#' @param domean if TRUE summarizes duplicate columns
#' @param modulesize the max number of S-genes included in one module
#' for search = "modules"
#' @param fpfn numeric vector of length two with false positive and
#' false negative rates
#' @param Rho optional perturbation matrix
#' @param logtype log base of the log odds
#' @param modified if TRUE, assumes a prepocessed data matrix
#' @param ... optional parameters for future search methods
#' @return transitively closed matrix or graphNEL
#' @author Martin Pirkl
#' @export
#' @import graph
#' @importFrom matrixStats rowMaxs
#' @examples
#' D <- matrix(rnorm(100*3), 100, 3)
#' colnames(D) <- 1:3
#' rownames(D) <- 1:100
#' adj <- diag(3)
#' colnames(adj) <- rownames(adj) <- 1:3
#' scoreAdj(D, adj)
#' @importFrom utils getFromNamespace
nem <- function(D, search = "greedy", start = NULL, method = "llr",
                  parallel = NULL, reduce = FALSE, weights = NULL, runs = 1,
                  verbose = FALSE, redSpace = NULL,
                  trans.close = TRUE, subtopo = NULL, prior = NULL,
                  ratio = TRUE, domean = TRUE, modulesize = 5,
                  fpfn = c(0.1, 0.1), Rho = NULL, logtype = 2,
                  modified = FALSE, ...) {
    if (method %in% "disc") {
        D[which(D == 1)] <- log((1-fpfn[2])/fpfn[1])/log(logtype)
        D[which(D == 0)] <- log(fpfn[2]/(1-fpfn[1]))/log(logtype)
        method <- "llr"
    }
    get.insertions <- get.ins.fast
    get.reversions <- get.rev.tc
    get.deletions <- get.del.tc
    D.backup <- D
    if (!modified) {
        D <- modData(D)
        colnames(D) <- gsub("\\..*", "", colnames(D))
    }
    if (is.null(Rho)) {
        Rho <- getRho(D)
    }
    if ("modules" %in% search) {
        if (length(search) > 1) {
            search <- search[-which(search %in% "modules")]
        } else {
            search <- "greedy"
        }
        if (length(unique(colnames(D))) > modulesize) {
            start <- modules(D, method = method, weights = weights,
                             reduce = reduce, verbose = verbose, start = start,
                             trans.close = trans.close, redSpace = redSpace,
                             subtopo = subtopo, fpfn = fpfn,
                             ratio = ratio, parallel = parallel, prior = prior,
                             modulesize = modulesize, search = search,
                             domean = domean, Rho = Rho, logtype = logtype)
        }
        if (search %in% "exhaustive") {
            search <- "greedy"
        }
    }
    if (length(unique(colnames(D))) == ncol(D)) {
        domean <- FALSE
    }
    Sgenes <- NULL
    if (domean) {
        D <- doMean(D, weights = weights, Rho = Rho, logtype = logtype)
        weights <- NULL
        Rho <- getRho(D)
        Sgenes <- getSgenes(D)
        domean <- FALSE
    }
    if (is.null(Sgenes)) {
        Sgenes <- getSgenes(D)
    }
    if (is.null(start)) {
        start2 <- "null"
        start <- better <- matrix(0, length(Sgenes), length(Sgenes))
        diag(start) <- 1
        colnames(start) <- rownames(start) <- Sgenes
    } else {
        if (length(start) == 1) {
            if (!(search %in% "estimate")) {
                start2 <- start
                start <- better <- matrix(0, length(Sgenes), length(Sgenes))
                diag(start) <- 1
                colnames(start) <- rownames(start) <- Sgenes
            } else {
                better <- matrix(0, length(Sgenes), length(Sgenes))
                start2 <- start
            }
        } else {
            better <- start2 <- start
            diag(start) <- 1
            colnames(start) <- rownames(start) <- Sgenes
        }
    }
    diag(better) <- 1
    colnames(better) <- rownames(better) <- Sgenes
    score <- scoreAdj(D, better, method = method, weights = weights,
                      subtopo = subtopo,
                      prior = prior, ratio = ratio, fpfn = fpfn,
                      Rho = Rho)
    score <- score$score
    oldscore <- score
    allscores <- score

    if (!is.null(parallel)) {
        get.insertions <- get.ins.fast
        get.reversions <- get.rev.tc
        get.deletions <- get.del.tc
        eigenMapMatMult <- mnem:::eigenMapMatMult
        naturalsort <- naturalsort::naturalsort
        sfInit(parallel = TRUE, cpus = parallel)
        sfExport("modules", "D", "start", "better", "transitive.reduction",
                 "method", "scoreAdj", "weights", "mytc",
                 "llrScore", "get.deletions", "get.insertions",
                 "get.reversions", "getRho", "doMean","eigenMapMatMult")
    }

    if (search %in% "small") {
        search <- "greedy"
        max_iter <- 1
    } else {
        max_iter <- Inf
    }

    guess <- FALSE
    if (search %in% "fast") {
        search <- "greedy"
        guess <- TRUE
    }

    if (search %in% "greedy") {
        P <- oldadj <- NULL
        if (guess) {
            better <- nem(D, search = "estimate", start = start,
                            method = method, parallel = parallel,
                            reduce = reduce, weights = weights, runs = runs,
                            verbose = verbose, redSpace = redSpace,
                            trans.close = trans.close, subtopo = subtopo,
                            prior = prior, ratio = ratio, domean = domean,
                            fpfn = fpfn, Rho = Rho, logtype = logtype,
                            modified = modified, Sgenes = Sgenes)$adj
            better <- mytc(better)
        }
        for (iter in seq_len(runs)) {
            if (iter > 1) {
                better <- matrix(sample(c(0,1),nrow(better)*ncol(better),
                                       replace = TRUE, prob = c(0.9,0.1)),
                                nrow(better),
                                ncol(better))
                colnames(better) <- rownames(better) <-
                    sample(Sgenes, length(Sgenes))
                better[lower.tri(better)] <- 0
                diag(better) <- 1
                better <- better[order(as.numeric(rownames(better))),
                                 order(as.numeric(colnames(better)))]
                score <- scoreAdj(D, better, method = method,
                                  weights = weights,
                                  subtopo = subtopo, prior = prior,
                                  ratio = ratio, fpfn = fpfn,
                                  Rho = Rho)
                P <- score$subweights
                oldadj <- better
                score <- score$score
                oldscore <- score
                allscores <- score
            }
            stop <- FALSE
            count <- 0
            while(!stop) {
                oldadj <- better
                doScores <- function(i) {
                    new <- models[[i]]
                    score <- scoreAdj(D, new, method = method,
                                      weights = weights,
                                      subtopo = subtopo, prior = prior,
                                      ratio = ratio, fpfn = fpfn,
                                      Rho = Rho, P = P, oldadj = oldadj,
                                      trans.close = FALSE)
                    P <- score$subweights
                    score <- score$score
                    return(list(score, P))
                }
                models <- unique(c(get.insertions(better),
                                   get.reversions(better),
                                   get.deletions(better)))
                if (is.null(parallel)) {
                    scores <- unlist(lapply((seq_len(length(models))),
                                            doScores), recursive = FALSE)
                } else {
                    scores <- unlist(sfLapply((seq_len(length(models))),
                                              doScores), recursive = FALSE)
                }
                Ps <- scores[seq_len(length(models))*2]
                scores <- unlist(scores[seq_len(length(models))*2 - 1])
                scores[is.na(scores)] <- 0
                bestidx <- which.max(scores)
                best <- models[[bestidx]]
                if ((max(scores, na.rm = TRUE) > oldscore |
                    ((max(scores, na.rm = TRUE) == oldscore &
                     sum(better == 1) > sum(best == 1))))
                    & count < max_iter) {
                    better <- best
                    oldscore <- max(scores)
                    allscores <- c(allscores, oldscore)
                    P <- Ps[[which.max(scores)]]
                } else {
                    stop <- TRUE
                }
                count <- count + 1
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
        models <- enumerate.models(length(Sgenes), Sgenes,
                                   trans.close = trans.close,
                                   verbose = verbose)
        doScores <- function(i) {
            adj <- models[[i]]
            score <- scoreAdj(D, adj, method = method, weights = weights,
                              subtopo = subtopo, prior = prior,
                              ratio = ratio, fpfn = fpfn,
                              Rho = Rho)
            score <- score$score
            return(score)
        }
        if (is.null(parallel)) {
            scores <- unlist(lapply(seq_len(length(models)), doScores))
        } else {
            scores <- unlist(sfLapply(seq_len(length(models)), doScores))
        }
        best <- which.max(scores)
        better <- mytc(models[[best]])
        allscores <- scores
        oldscore <- max(scores)
        diag(better) <- 1
    }

    if (search %in% "estimate") {
        if (!is.null(weights)) {
            Dw <- D*rep(weights, rep(nrow(D), ncol(D)))
            weights <- NULL
        } else {
            Dw <- D
        }
        tmp <- nemEst(Dw, start = "null", method = method, fpfn = fpfn,
                      Rho = Rho, domean = domean, modified = TRUE, ...)
        tmp1 <- nemEst(Dw, start = "full", method = method, fpfn = fpfn,
                       Rho = Rho, domean = domean, modified = TRUE, ...)
        if (tmp1$ll > tmp$ll) { tmp <- tmp1 }
        if (is.matrix(start2)) {
            tmp2 <- nemEst(Dw, start = start2, method = method, fpfn = fpfn,
                           Rho = Rho, domean = domean, modified = TRUE, ...)
            if (tmp2$ll > tmp$ll) { tmp <- tmp2 }
        }
        better <- tmp$phi
        oldscore <- tmp$ll
        allscores <- tmp$lls
        subweights <- Dw%*%cbind(t(Rho)%*%tmp$phi, 0)
    }

    if (!is.null(parallel)) {
        sfStop()
    }

    if (is.null(subtopo)) {
        subtopo <- scoreAdj(D, better, method = method, weights = weights,
                            prior = prior,
                            ratio = ratio, fpfn = fpfn,
                            Rho = Rho, dotopo = TRUE)
        subweights <- subtopo$subweights
        subtopo <- subtopo$subtopo
    } else {
        subweights <- P
    }

    better <- transitive.reduction(better)
    better <- better[order(as.numeric(rownames(better))),
                     order(as.numeric(colnames(better)))]
    nem <- list(adj = better, score = oldscore, scores = allscores,
                redSpace = redSpace, subtopo = subtopo, D = D.backup,
                subweights = subweights)
    return(nem)
}
#' Network score
#'
#' Computes the fit (score of a network) of the data given a network matrix
#' @param D data matrix; use modified = FALSE
#' @param adj adjacency matrix of the network phi
#' @param method either llr if D consists of log odds or disc,
#' if D is binary
#' @param logtype log base of the log odds
#' @param weights a numeric vector of weights for the columns of D
#' @param trans.close if TRUE uses the transitive closure of adj
#' @param subtopo optional matrix with the subtopology theta as
#' adjacency matrix
#' @param prior a prior network matrix for adj
#' @param ratio if FALSE uses alternative distance for the model score
#' @param fpfn numeric vector of length two with false positive and
#' false negative rates
#' @param Rho optional perturbation matrix
#' @param dotopo if TRUE computes and returns the subtopology theta (optional)
#' @param P previous score matrix (only used internally)
#' @param oldadj previous adjacency matrix (only used internally)
#' @param modified if TRUE, assumes a prepocessed data matrix
#' @return transitively closed matrix or graphNEL
#' @author Martin Pirkl
#' @export
#' @import graph
#' @importFrom matrixStats rowMaxs
#' @examples
#' D <- matrix(rnorm(100*3), 100, 3)
#' colnames(D) <- 1:3
#' rownames(D) <- 1:100
#' adj <- diag(3)
#' colnames(adj) <- rownames(adj) <- 1:3
#' scoreAdj(D, adj)
scoreAdj <- function(D, adj, method = "llr", logtype = 2, weights = NULL,
                     trans.close = TRUE, subtopo = NULL,
                     prior = NULL, ratio = TRUE, fpfn = c(0.1, 0.1),
                     Rho = NULL, dotopo = FALSE,
                     P = NULL, oldadj = NULL, modified = TRUE) {
    if (!modified) { D <- modData(D) }
    if (method %in% "disc") {                                             
        D[which(D == 1)] <- log((1 - fpfn[2])/fpfn[1])/log(logtype)       
        D[which(D == 0)] <- log(fpfn[2]/(1 - fpfn[1]))/log(logtype)       
        method <- "llr"                                                   
    }   
    if (trans.close) {
        adj <- mytc(adj)
    }
    if (is.null(Rho)) {
        adj1 <- adj[colnames(D), ]
    } else {
        adj1 <- t(Rho)%*%adj
        adj1[which(adj1 > 1)] <- 1
    }
    if (method %in% "llr") {
        ll <- "max"
        if (is.null(P) | is.null(oldadj)) {
            score <- llrScore(D, adj1, weights = weights, ratio = ratio)
        } else {
            score <- P
            changeidx <- unique(floor(which(adj != oldadj)/nrow(adj)+0.99))
            score[, changeidx] <- llrScore(D, adj1[, changeidx],
                                           weights = weights, ratio = ratio)
        }
    }
    if (is.null(subtopo) & dotopo) {
        subtopo <- maxCol_row(cbind(0, score))
        subtopo <- subtopo - 1
    }
    subweights <- score
    if (ll %in% "max") {
        if (is.null(subtopo)) {
            score[which(score < 0)] <- 0
            score <- sum(matrixStats::rowMaxs(score))
        } else {
            score <- sum(score[cbind(seq_len(nrow(score)),
                                     as.numeric(subtopo))])
        }
    }
    if (ll %in% "marg") {
        score <- sum(score)
    }
    if (!is.null(prior)) {
        prior <- transitive.reduction(prior)
        adj <- transitive.reduction(adj)
        score <- score - sum(abs(prior - adj))/length(prior)
    }
    return(list(score = score, subtopo = subtopo, subweights = subweights))
}
#' Transitive reduction
#'
#' Computes the transitive reduction of an adjacency
#' matrix or graphNEL object.
#' Originally imported from the package 'nem'.
#' @param g adjacency matrix or graphNEL object
#' @author Holger Froehlich
#' @references R. Sedgewick, Algorithms, Pearson, 2002.
#' @return transitively reduced adjacency matrix
#' @export
#' @importFrom methods as
#' @examples
#' g <- matrix(c(0,0,0,1,0,0,0,1,0), 3)
#' rownames(g) <- colnames(g) <- seq_len(3)
#' g.tr <- transitive.reduction(g)
transitive.reduction <- function (g) {
    if (!(is(g, "matrix") | is(g, "graphNEL"))) 
        stop("Input must be an adjacency matrix or graphNEL object")
    if (is(g, "graphNEL")) {
        g = as(g, "matrix")
    }
    g = transitive.closure(g)
    g = g - diag(diag(g))
    type = (g > 1) * 1 - (g < 0) * 1
    for (y in seq_len(nrow(g))) {
        for (x in seq_len(nrow(g))) {
            if (g[x, y] != 0) {
                for (j in seq_len(nrow(g))) {
                  if ((g[y, j] != 0) & sign(type[x, j]) * sign(type[x, 
                    y]) * sign(type[y, j]) != -1) {
                    g[x, j] = 0
                  }
                }
            }
        }
    }
    g
}
#' Transitive closure of a directed acyclic graph (dag)
#'
#' Computes the transitive closure of a dag or only of a
#' deletion/addition of an edge
#' @param g graph as matrix or graphNEL object
#' @param u index of the parent of an edge (optional)
#' @param v index of the child of an edge (optional)
#' @return transitively closed matrix or graphNEL
#' @author Martin Pirkl
#' @export
#' @import graph
#' @importFrom methods as
#' @examples
#' g <- matrix(c(0,0,0,1,0,0,0,1,0), 3)
#' transitive.closure(g)
transitive.closure <- function(g, u = NULL, v = NULL) {
    if (is(g, "graphNEL")) {
        a <- as(g, "matrix")
    } else if (is(g, "matrix")) {
        a <- g
    } else {
        stop("g has to be a matrix or graphNEL object")
    }
    gtc <- mytc(a, u = u, v = v)
    if (is(g, "graphNEL")) {
        gtc <- as(gtc, "graphNEL")
    }
    return(gtc)
}
#' Simulation accuracy.
#'
#' Computes the accuracy of the fit between simulated and
#' inferred mixture.
#' @param x mnem object
#' @param y simulation object or another mnem object
#' @param strict if TRUE, accounts for over/underfitting, i.e.
#' the number of components
#' @param unique if TRUE, phis of x and y are made unique each
#' (FALSE if strict is TRUE)
#' @param type type of accuracy. "ham" for hamming, "sens" for
#' sensitivity and "spec for Specificity"
#' @return plot of EM convergence
#' @author Martin Pirkl
#' @export
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' result <- mnem(data, k = 2, starts = 1)
#' fitacc(result, sim)
#' fitacc(result, sim, type = "sens")
#' fitacc(result, sim, type = "spec")
#' fitacc(result, sim, strict = TRUE, type = "sens")
#' fitacc(result, sim, strict = TRUE, type = "spec")
fitacc <- function(x, y, strict = FALSE, unique = TRUE,
                   type = "ham") {
    for (i in seq_len(length(x$comp))) {
        x$comp[[i]]$theta <- NULL
    }
    if (strict) {
        unique <- FALSE
    }
    if (is.null(y$Nem)) {
        y2 <- y
        y$Nem <- list()
        for (i in seq_len(length(y$comp))) {
            y$Nem[[i]] <- y2$comp[[i]]$phi
        }
    }
    if (unique) {
        x <- unique(x$comp)
        y <- unique(y$Nem)
    } else {
        x <- x$comp
        y <- y$Nem
    }
    xn <- length(x)
    yn <- length(y)
    n <- nrow(x[[1]]$phi)
    if (strict) {
        score <- 0
        while(length(x) > 0 & length(y) > 0) {
            couple <- numeric(2)
            best <- -Inf
            for (i in seq_len(length(x))) {
                for (j in seq_len(length(y))) {
                    A <- mytc(x[[i]]$phi)
                    B <- mytc(y[[j]])
                    tmp <- (n*(n-1) - sum(abs(A - B)))/(n*(n-1))
                    comp <- tmp >= best
                    comp2 <- TRUE
                    if (type %in% "sens") {
                        tp <- sum(A == 1 & B == 1)
                        fn <- sum(A == 0 & B == 1)
                        tmp <- tp/(tp+fn)
                        comp2 <- tmp >= best
                    } else if (type %in% "spec") {
                        tn <- sum(A == 0 & B == 0)
                        fp <- sum(A == 1 & B == 0)
                        tmp <- tn/(tn+fp)
                        comp2 <- tmp >= best
                    }
                    if (comp & comp2) {
                        couple <- c(i,j)
                        best <- tmp
                    }
                }
            }
            score <- best + score
            x[[couple[1]]] <- NULL
            y[[couple[2]]] <- NULL
        }
        if (length(x) != 0) {
            for (i in seq_len(length(x))) {
                A <- mytc(x[[i]]$phi)
                B <- diag(n)*0
                if (type %in% "ham") {
                    tmp <- (n*(n-1) - sum(abs(A - B)))/(n*(n-1))
                } else if (type %in% "sens") {
                    tp <- sum(A == 1 & B == 1)
                    fn <- sum(A == 0 & B == 1)
                    tmp <- tp/(tp+fn)
                    tmp <- 0
                } else if (type %in% "spec") {
                    tn <- sum(A == 0 & B == 0)
                    fp <- sum(A == 1 & B == 0)
                    tmp <- tn/(tn+fp)
                }
                score <- score + tmp
            }
        }
        if (length(y) != 0) {
            for (i in seq_len(length(y))) {
                A <- diag(n)*0
                B <- mytc(y[[i]])
                if (type %in% "ham") {
                    tmp <- (n*(n-1) - sum(abs(A - B)))/(n*(n-1))
                } else if (type %in% "sens") {
                    tp <- sum(A == 1 & B == 1)
                    fn <- sum(A == 0 & B == 1)
                    tmp <- tp/(tp+fn)
                } else if (type %in% "spec") {
                    tn <- sum(A == 0 & B == 0)
                    fp <- sum(A == 1 & B == 0)
                    tmp <- tn/(tn+fp)
                }
                score <- score + tmp
            }
        }
        score <- score/max(c(xn,yn))
    } else {
        xmax <- numeric(xn)
        ymax <- numeric(yn)
        best <- -Inf
        for (i in seq_len(xn)) {
            for (j in seq_len(yn)) {
                A <- mytc(x[[i]]$phi)
                B <- mytc(y[[j]])
                if (type %in% "ham") {
                    tmp <- (n*(n-1) - sum(abs(A - B)))/(n*(n-1))
                } else if (type %in% "sens") {
                    tp <- sum(A == 1 & B == 1)
                    fn <- sum(A == 0 & B == 1)
                    tmp <- tp/(tp+fn)
                } else if (type %in% "spec") {
                    tn <- sum(A == 0 & B == 0)
                    fp <- sum(A == 1 & B == 0)
                    tmp <- tn/(tn+fp)
                }
                if (tmp >= best) {
                    couple <- c(i,j)
                    best <- tmp
                }
                if (tmp > xmax[i]) {
                    xmax[i] <- tmp
                }
                if (tmp > ymax[j]) {
                    ymax[j] <- tmp
                }
            }
        }
        score <- sum(c(xmax,ymax))/(xn+yn)
    }
    return(score)
}
#' Plot convergence of EM
#'
#' This function plots the convergence of the different EM iterations (four
#' figures, e.g. par(mfrow=(2,2))).
#' @param x mnem object
#' @param col vector of colors for the iterations
#' @param type see ?plot.default
#' @param convergence difference of when two log likelihoods
#' are considered equal; see also convergence for the function
#' mnem()
#' @param ... additional parameters for the plots/lines functions
#' @author Martin Pirkl
#' @return plot of EM convergence
#' @export
#' @importFrom grDevices rgb
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' result <- mnem(data, k = 2, starts = 1)
#' par(mfrow=c(2,2))
#' plotConvergence(result)
plotConvergence <- function(x, col = NULL, type = "b",
                                 convergence = 0.1, ...) {
    runs <- length(x$limits)
    if (is.null(col)) {
        col <- rgb(runif(runs),runif(runs),runif(runs),0.75)
    }
    ymin <- Inf
    xmax <- 0
    maxll <- -Inf
    maxllcount <- 0
    minlen <- 0
    for (i in seq_len(runs)) {
        if (max(x$limits[[i]]$ll) - maxll > convergence) {
            maxll <- max(x$limits[[i]]$ll)
            maxllcount <- 0
        }
        if (abs(max(x$limits[[i]]$ll) - maxll) <= convergence) {
            maxllcount <- maxllcount + 1
        }
        if (length(x$limits[[i]]$ll) == 2) {
            minlen <- minlen + 1
        }
        ymin <- min(c(ymin, x$limits[[i]]$ll))
        xmax <- max(c(xmax, length(x$limits[[i]]$ll)))
    }
    for (i in seq_len(runs)) {
        if (i == 1) {
            plot(x$limits[[i]]$ll, col = col[i], ylim = c(ymin, max(x$lls)),
                 xlim = c(1, xmax), xlab = "EM iterations", ylab = "log odds",
                 type = type,
                 main = paste0(maxllcount,
                               " (",
                               round(maxllcount/runs*100, 2),
                               "%) run(s) with maximum log odds at ",
                               convergence, " distance\n",
                               minlen, " (", round(minlen/runs*100, 2),
                               "%) run(s) started in local optimum"),
                 ...)
        } else {
            lines(x$limits[[i]]$ll, col = col[i], type = type, ...)
        }
    }
    if (!is.null(x$limits[[1]]$phievo)) {
        ymin <- Inf
        ymax <- -Inf
        xmax <- 0
        maxll <- -Inf
        maxllcount <- 0
        minlen <- 0
        for (i in seq_len(runs)) {
            if (max(x$limits[[i]]$phievo) - maxll > convergence) {
                maxll <- max(x$limits[[i]]$phievo)
                maxllcount <- 0
            }
            if (abs(max(x$limits[[i]]$phievo) - maxll) <= convergence) {
                maxllcount <- maxllcount + 1
            }
            if (length(x$limits[[i]]$phievo) == 2) {
                minlen <- minlen + 1
            }
            ymin <- min(c(ymin, x$limits[[i]]$phievo))
            ymax <- max(c(ymax, x$limits[[i]]$phievo))
            xmax <- max(c(xmax, length(x$limits[[i]]$phievo)))
        }
        for (i in seq_len(runs)) {
            if (i == 1) {
                plot(x$limits[[i]]$phievo, col = col[i], ylim = c(ymin, ymax),
                     xlim = c(1, xmax), xlab = "EM iterations",
                     main = expression(evolution ~ of ~ phi),
                     ylab = "edge changes",
                     type = type,
                     ...)
            } else {
                lines(x$limits[[i]]$phievo, col = col[i], type = type, ...)
            }
        }
        ymin <- Inf
        ymax <- -Inf
        xmax <- 0
        maxll <- -Inf
        maxllcount <- 0
        minlen <- 0
        for (i in seq_len(runs)) {
            if (max(x$limits[[i]]$thetaevo) - maxll > convergence) {
                maxll <- max(x$limits[[i]]$thetaevo)
                maxllcount <- 0
            }
            if (abs(max(x$limits[[i]]$thetaevo) - maxll) <= convergence) {
                maxllcount <- maxllcount + 1
            }
            if (length(x$limits[[i]]$thetaevo) == 2) {
                minlen <- minlen + 1
            }
            ymin <- min(c(ymin, x$limits[[i]]$thetaevo))
            ymax <- max(c(ymax, x$limits[[i]]$thetaevo))
            xmax <- max(c(xmax, length(x$limits[[i]]$thetaevo)))
        }
        for (i in seq_len(runs)) {
            if (i == 1) {
                plot(x$limits[[i]]$thetaevo, col = col[i], ylim = c(ymin, ymax),
                     xlim = c(1, xmax), xlab = "EM iterations",
                     main = expression(evolution ~ of ~ theta),
                     ylab = "edge changes",
                     type = type,
                     ...)
            } else {
                lines(x$limits[[i]]$thetaevo, col = col[i], type = type, ...)
            }
        }
        ymin <- Inf
        ymax <- -Inf
        xmax <- 0
        maxll <- -Inf
        maxllcount <- 0
        minlen <- 0
        for (i in seq_len(runs)) {
            if (max(x$limits[[i]]$mwevo) - maxll > convergence) {
                maxll <- max(x$limits[[i]]$mwevo)
                maxllcount <- 0
            }
            if (abs(max(x$limits[[i]]$mwevo) - maxll) <= convergence) {
                maxllcount <- maxllcount + 1
            }
            if (length(x$limits[[i]]$mwevo) == 2) {
                minlen <- minlen + 1
            }
            ymin <- min(c(ymin, x$limits[[i]]$mwevo))
            ymax <- max(c(ymax, x$limits[[i]]$mwevo))
            xmax <- max(c(xmax, length(x$limits[[i]]$mwevo)))
        }
        for (i in seq_len(runs)) {
            if (i == 1) {
                plot(x$limits[[i]]$mwevo, col = col[i], ylim = c(ymin, ymax),
                     xlim = c(1, xmax), xlab = "EM iterations",
                     main = expression(evolution ~ of ~ pi),
                     ylab = "sum of absolute distances",
                     type = type,
                     ...)
            } else {
                lines(x$limits[[i]]$mwevo, col = col[i], type = type, ...)
            }
        }
    }
}
#' Creating app data.
#'
#' This function is for the reproduction of the application
#' results in the vignette and publication. See the publication
#' Pirkl & Beerenwinkel (2018) on how to download the data files:
#' GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv
#' k562_both_filt.txt
#' GSM2396861_k562_ccycle_cbc_gbc_dict.csv
#' GSM2396858_k562_tfs_7_cbc_gbc_dict.csv
#' @param sets numeric vector with the data sets: 1 (CROPseq),
#' 2, 3 (both PERTURBseq); default is all three
#' @param m number of Sgenes (for testing)
#' @param n number of most variable E-genes (for testing)
#' @param o number of samples per S-gene (for testing)
#' @param maxk maximum number of component in mnem inference (default: 5)
#' @param parallel number of threads for parallelisation
#' @param path path to the data files path/file.csv: "path/"
#' @param types types of data/analysis; "data" creates the gene expression
#' matrix, "lods" includes the log odds, "mnem" additionally performes the
#' mixture nem analysis; default c("data", "lods", "mnem")
#' @param allcrop if TRUE, does not restrict and uses the full CROPseq dataset
#' @param multi if TRUE, includes cells with more than one perturbed gene
#' @param file path and filename of the rda file with the raw data from the
#' command "data <- createApp(..., types = "data")"
#' @param ... additional parameters for the mixture nem function
#' @return app data object
#' @author Martin Pirkl
#' @export
#' @importFrom data.table fread
#' @importFrom utils read.csv read.table
#' @import snowfall Linnorm
#' @examples
#' ## recreate the app data object (takes very long, i.e. days)
#'\dontrun{
#' createApp()
#' }
#' data(app)
createApp <- function(sets = seq_len(3), m = NULL, n = NULL, o = NULL,
                      maxk = 5, parallel = NULL, path = "",
                      types = c("data", "lods", "mnem"),
                      allcrop = FALSE, multi = FALSE, file = NULL, ...) {
    if ("mnem" %in% types) {
        types <- c(types, "lods")
    }
    if ("lods" %in% types & is.null(file)) {
        types <- c(types, "data")
    }
    data <- lods <- NULL
    ## load datasets
    barcodes <- list()
    if ("data" %in% types) {
        datas <- list()
        if (1 %in% sets) {
            datafile <- "GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv"
            data <- fread(paste0(path, datafile), data.table = FALSE)
            data.backup <- data
            counts <- data[, grep("condition|^stim", colnames(data))]
            counts <- counts[-(seq_len(5)), -1]
            counts <- matrix(as.numeric(as.character(unlist(counts))),
                             nrow(counts))
            rownames(counts) <- data[6:nrow(data), 1]
            colnames(counts) <- data[4, grep("^stim", colnames(data))]
            colnames(counts)[which(colnames(counts) %in% "CTRL")] <- ""
            counts <- counts[, -grep("DHODH|MVD|TUBB", colnames(counts))]
            data <- counts
            barcodes[[1]] <- colnames(data)
            if (!is.null(n)) {
                var <- apply(data, 1, var)
                data <- data[order(var, decreasing = TRUE)[seq_len(n)], ,
                             drop = FALSE]
            }
            if (!is.null(m)) {
                Sgenes <- unique(c("", naturalsort(unique(colnames(data)))))
                data <- data[, which(colnames(data) %in% Sgenes[seq_len(m+1)]),
                             drop = FALSE]
            }
            if (!is.null(o)) {
                Sgenes <- naturalsort(unique(colnames(data)))
                idx <- NULL
                for (i in Sgenes) {
                    idx <- c(idx, which(colnames(data) %in% i)[seq_len(o)])
                }
                data <- data[, idx]
            }
            datas[[1]] <- data
        }
        if (2 %in% sets | 3 %in% sets) {
            data <- fread(paste0(path, "k562_both_filt.txt"))
            rns <- data$GENE
            tmp <- grep("cc7d", colnames(data))
            data1 <- data[, ..tmp]
            tmp <- grep("p7d", colnames(data))
            data2 <- data[, ..tmp]
            ## alternative:
            ## data <- read.table(paste0(path, "k562_both_filt.txt"),
            ## header = TRUE)
            ## rns <- data$GENE
            ## tmp <- grep("cc7d", colnames(data))
            ## data1 <- data[, tmp]
            ## tmp <- grep("p7d", colnames(data))
            ## data2 <- data[, tmp]
            rm(data)
            gc()
            data <- data1
            rm(data1)
            gc()
            data <- as.matrix(data)
            rownames(data) <- rns
            if (!is.null(n)) {
                var <- apply(data, 1, var)
                data <- data[order(var, decreasing = TRUE)[seq_len(n)], ,
                             drop = FALSE]
            }
            datafile <- "GSM2396861_k562_ccycle_cbc_gbc_dict.csv"
            c2g <- read.table(paste0(path, datafile), fill = TRUE, sep = ",")
            cn <- character(ncol(data))
            for (i in seq_len(length(c2g[[1]]))) {
                cells <- unlist(strsplit(as.character(c2g[[2]][[i]]), ", "))
                cn[which(colnames(data) %in% cells)] <- paste(cn[
                    which(colnames(data) %in% cells)],
                    gsub("^c_|^c_sg|^m_|_[0-9]$|_10$", "",
                         as.character(c2g[[1]][[i]])), sep = "_")
            }
            barcodes[[2]] <- colnames(data)
            colnames(data) <- gsub("^_|^_p_sg|^_p_", "", cn)
            colnames(data)[grep("INTER", colnames(data))] <- ""
            if (length(grep("_", colnames(data))) > 0 & !multi) {
                data <- data[, -grep("_", colnames(data))]
            }
            if (!is.null(m)) {
                Sgenes <- unique(c("", naturalsort(unique(colnames(data)))))
                data <- data[, which(colnames(data) %in% Sgenes[seq_len(m+1)]),
                             drop = FALSE]
            }
            if (!is.null(o)) {
                Sgenes <- naturalsort(unique(colnames(data)))
                idx <- NULL
                for (i in Sgenes) {
                    idx <- c(idx, which(colnames(data) %in% i)[seq_len(o)])
                }
                data <- data[, idx]
            }
            datas[[2]] <- data
            rm(data)
            gc()
            data <- data2
            rm(data2)
            gc()
            data <- as.matrix(data)
            rownames(data) <- rns
            if (!is.null(n)) {
                var <- apply(data, 1, var)
                data <- data[order(var, decreasing = TRUE)[seq_len(n)], ,
                             drop = FALSE]
            }
            datafile <- "GSM2396858_k562_tfs_7_cbc_gbc_dict.csv"
            c2g <- read.table(paste0(path, datafile), fill = TRUE, sep = ",")
            cn <- character(ncol(data))
            for (i in seq_len(length(c2g[[1]]))) {
                cells <- unlist(strsplit(as.character(c2g[[2]][[i]]), ", "))
                cn[which(colnames(data) %in% cells)] <- paste(cn[
                    which(colnames(data) %in% cells)],
                    gsub("^c_|^c_sg|^m_|_[0-9]$|_10$", "",
                         as.character(c2g[[1]][[i]])), sep = "_")
            }
            barcodes[[3]] <- colnames(data)
            colnames(data) <- gsub("^_|^_p_sg|^_p_", "", cn)
            colnames(data)[grep("INTER", colnames(data))] <- ""
            if (length(grep("_", colnames(data))) > 0 & !multi) {
                data <- data[, -grep("_", colnames(data))]
            }
            if (!is.null(m)) {
                Sgenes <- unique(c("", naturalsort(unique(colnames(data)))))
                data <- data[, which(colnames(data) %in% Sgenes[seq_len(m+1)]),
                             drop = FALSE]
            }
            if (!is.null(o)) {
                Sgenes <- naturalsort(unique(colnames(data)))
                idx <- NULL
                for (i in Sgenes) {
                    idx <- c(idx, which(colnames(data) %in% i)[seq_len(o)])
                }
                data <- data[, idx]
            }
            datas[[3]] <- data
        }
        if (4 %in% sets | 5 %in% sets) {
            data <- fread(paste0(path, "dc_both_filt_fix_tp10k.txt"))
            data.backup <- data
            rownames(data) <- data.backup$GENE
            if (!is.null(n)) {
                var <- apply(data, 1, var)
                data <- data[order(var, decreasing = TRUE)[seq_len(n)], ,
                             drop = FALSE]
            }
            data <- data[, -1]
            data <- as.matrix(data)
            data1 <- data[, grep("0h", colnames(data))]
            data2 <- data[, grep("3h", colnames(data))]
            data <- data1
            datafile <- "GSM2396857_dc_0hr_cbc_gbc_dict.csv"
            c2g <- read.table(paste0(path, datafile), fill = TRUE, sep = ",")
            cn <- character(ncol(data))
            for (i in seq_len(length(c2g[[1]]))) {
                cells <- unlist(strsplit(as.character(c2g[[2]][[i]]), ", "))
                cn[which(colnames(data) %in% cells)] <- paste(cn[
                    which(colnames(data) %in% cells)],
                    gsub("^c_|^c_sg|^m_|_[0-9]$|_10$", "",
                         as.character(c2g[[1]][[i]])), sep = "_")
            }
            barcodes[[4]] <- colnames(data)
            colnames(data) <- gsub("^_|^_p_sg|^_p_", "", cn)
            colnames(data)[grep("INTER", colnames(data))] <- ""
            if (length(grep("_", colnames(data))) > 0 & !multi) {
                data <- data[, -grep("_", colnames(data))]
            }
            if (!is.null(m)) {
                Sgenes <- unique(c("", naturalsort(unique(colnames(data)))))
                data <- data[, which(colnames(data) %in% Sgenes[seq_len(m+1)]),
                             drop = FALSE]
            }
            if (!is.null(o)) {
                Sgenes <- naturalsort(unique(colnames(data)))
                idx <- NULL
                for (i in Sgenes) {
                    idx <- c(idx, which(colnames(data) %in% i)[seq_len(o)])
                }
                data <- data[, idx]
            }
            datas[[4]] <- data
            data <- data2
            datafile <- "GSM2396856_dc_3hr_cbc_gbc_dict_strict.csv"
            c2g <- read.table(paste0(path, datafile), fill = TRUE, sep = ",")
            cn <- character(ncol(data))
            for (i in seq_len(length(c2g[[1]]))) {
                cells <- unlist(strsplit(as.character(c2g[[2]][[i]]), ", "))
                cn[which(colnames(data) %in% cells)] <- paste(cn[
                    which(colnames(data) %in% cells)],
                    gsub("^c_|^c_sg|^m_|_[0-9]$|_10$", "",
                         as.character(c2g[[1]][[i]])), sep = "_")
            }
            barcodes[[5]] <- colnames(data)
            colnames(data) <- gsub("^_|^_p_sg|^_p_", "", cn)
            colnames(data)[grep("INTER", colnames(data))] <- ""
            if (length(grep("_", colnames(data))) > 0 & !multi) {
                data <- data[, -grep("_", colnames(data))]
            }
            if (!is.null(m)) {
                Sgenes <- unique(c("", naturalsort(unique(colnames(data)))))
                data <- data[, which(colnames(data) %in% Sgenes[seq_len(m+1)]),
                             drop = FALSE]
            }
            if (!is.null(o)) {
                Sgenes <- naturalsort(unique(colnames(data)))
                idx <- NULL
                for (i in Sgenes) {
                    idx <- c(idx, which(colnames(data) %in% i)[seq_len(o)])
                }
                data <- data[, idx]
            }
            datas[[5]] <- data
            data <- cbind(datas[[4]], datas[[5]])
            data <- exp(data) - 1
            exprslvl <- apply(data, 1, median)
            data <- data[which(exprslvl > 0), ]
            data <- Linnorm(data)
            datas[[4]] <- data[, seq_len(ncol(datas[[4]]))]
            datas[[5]] <- data[, -seq_len(ncol(datas[[4]]))]
            print("data normalized")
        }
    } else {
        load(file)
        datas <- list()
        for (i in sets) {
            if (!is.null(data[[i]])) {
                datas[[i]] <- data[[i]]$data
            }
        }
    }
    print("data loaded")
    ## dataset normalization, log odds computation and mnem inference
    app <- list()
    for (i in sets) {
        data <- datas[[i]]
        if ("lods" %in% types) {
            print(paste0("dataset ", i))
            ## data normalization
            if (i == 1) {
                exprslvl <- apply(data, 1, median)
                data <- data[which(exprslvl > 0), ]
                data <- t(t(data)/(colSums(data)/10000))
                data <- Linnorm(data)
                print("data normalized")
            } else if (i %in% c(2,3)) {
                data <- exp(data) - 1
                exprslvl <- apply(data, 1, median)
                data <- data[which(exprslvl > 0), ]
                data <- Linnorm(data)
                print("data normalized")
            }
            ## log ratio calculation
            llr <- data*0
            C <- which(colnames(data) %in% "")
            distrPar <- function(i, data, C) {
                cdistr <- ecdf(data[i, C])
                llrcol <- numeric(ncol(data))
                for (j in which(!(colnames(data) %in% ""))) {
                    gene <- colnames(data)[j]
                    D <- which(colnames(data) %in% gene)
                    ddistr <- ecdf(data[i, D])
                    llrcol[j] <-
                        log2(min(ddistr(data[i, j]),
                                 1 -
                                 ddistr(data[i, j]))/min(cdistr(data[i,j]),
                                                         1 -
                                                         cdistr(data[i,j])))
                }
                return(llrcol)
            }
            if (is.null(parallel)) {
                llr <- lapply(seq_len(nrow(data)), distrPar, data, C)
            } else {
                sfInit(parallel = TRUE, cpus = parallel)
                llr <- sfLapply(seq_len(nrow(data)), distrPar, data, C)
                sfStop()
            }
            llr <- do.call("rbind", llr)
            llr[is.na(llr)] <- 0
            llr[is.infinite(llr)] <- max(llr[!is.infinite(llr)])
            colnames(llr) <- colnames(data)
            llr <- llr[, which(!(colnames(data) %in% ""))]
            rownames(llr) <- rownames(data)
            print("log ratios computed")
            ## mixture nem
            colnames(llr) <- toupper(colnames(llr))
            if (i == 1 & !allcrop) {
                cropgenes <- c("LCK", "ZAP70", "PTPN6", "DOK2",
                               "PTPN11", "EGR3", "LAT")
                lods <- llr[, which(colnames(llr) %in% cropgenes)]
            } else {
                lods <- llr
            }
            badgenes <- "Tcrlibrary"
            badgenes <- grep(badgenes, rownames(lods))
            if (length(badgenes) > 0) {
                lods <- lods[-badgenes, ]
            }
            sdev <- apply(lods, 1, sd)
            lods <- lods[which(sdev > sd(lods)), ]
            n <- length(unique(colnames(lods)))
            if ("mnem" %in% types) {
                bics <- rep(Inf, maxk)
                res <- list()
                for (k in seq_len(maxk)) {
                    res[[k]] <- mnem(lods, k = k, parallel = parallel, ...)
                }
                app[[i]] <- res
                print("mixture nested effects model learned")
            } else {
                app[[i]] <- list(data = data, lods = lods,
                                 barcodes = barcodes[[i]])
            }
        } else {
            app[[i]] <- list(data = data, barcodes = barcodes[[i]])
        }
    }
    return(app)
}
#' Hierarchical mixture.
#'
#' This function does a hierarchical mixture. That means it uses the
#' approximate BIC to check, if there are more than one component. It
#' recursively splits the data if there is evidence for k > 1 components.
#' @param data data matrix either binary or log odds
#' @param k number of maximal components for each hierarchy leaf
#' @param logtype log type of the data
#' @param getprobspars list of parameters for the getProbs function
#' @param ... additional parameters for the mnem function
#' @return object of class mnem
#' @author Martin Pirkl
#' @export
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' result <- mnemh(data, starts = 1, k = 1)
mnemh <- function(data, k = 2, logtype = 2, getprobspars = list(), ...) {
    D <- data
    data <- modData(data)
    n <- getSgeneN(data)
    Sgenes <- naturalsort(getSgenes(data))
    tmp <- mnemh.rec(data, k=k, logtype=logtype, ...)
    K <- length(tmp)/11
    comp <- res <- limits <- list()
    lls <- numeric(K)
    for (i in seq_len(K)) {
        lls[i] <- tmp[[(11+11*(i-1))]]
        comp[[i]] <- res[[i]] <- list()
        tmp2 <- tmp[[(2+11*(i-1))]]
        tmp3 <- tmp2[[1]]$phi
        thetatmp <- tmp[[(2+11*(i-1))]][[1]]$theta
        if (nrow(tmp2[[1]]$phi) < n) {
            tmpdata <- tmp[[(3+11*(i-1))]]
            inGenes <- naturalsort(getSgenes(tmpdata))
            tmpidx <- which(Sgenes %in% inGenes)
            if (max(thetatmp) > length(tmpidx)) {
                thetatmp[which(thetatmp == max(thetatmp))] <- length(Sgenes)+1
            }
            for (j in seq_len(length(tmpidx))) {
                thetatmp <- gsub(j, tmpidx[j], thetatmp)
            }
            outGenes <- Sgenes[which(!(Sgenes %in% inGenes))]
            tmpidx2 <- which(Sgenes %in% outGenes)
            colnames(tmp3) <- rownames(tmp3) <- tmpidx
            tmp3 <- cbind(matrix(0, nrow(tmp3)+length(outGenes),
                                 length(outGenes)),
                          rbind(matrix(0, length(outGenes), ncol(tmp3)), tmp3))
            colnames(tmp3)[seq_len(length(outGenes))] <-
                rownames(tmp3)[seq_len(length(outGenes))] <- tmpidx2
            tmp3 <- tmp3[naturalorder(rownames(tmp3)),
                         naturalorder(colnames(tmp3))]
        }
        res[[i]]$adj <- comp[[i]]$phi <- tmp3
        comp[[i]]$theta <- tmp[[(2+11*(i-1))]][[1]]$theta
        limits[[i]] <- list()
        limits[[i]]$ll <- tmp[[(6+11*(i-1))]]
    }
    probs <- matrix(0, K, ncol(data))
    probs <- do.call(getProbs, c(list(probs=probs, k=K, data=data,
                                      res=res, n=n), getprobspars))
    colnames(probs$probs) <- colnames(D)
    res <- list(limits = limits, comp = comp, data = D, mw = probs$mw,
                probs = probs$probs, lls = lls, ll = probs$ll)
    class(res) <- "mnem"
    return(res)
}
#' Processed scRNAseq from pooled CRISPR screens
#'
#' Example data: mnem results for
#' the Dixit et al., 2016 and Datlinger et al., pooled CRISPR screens.
#' For details see the vignette or function createApp().
#' @name app
#' @docType data
#' @usage app
#' @references Datlinger, P., Rendeiro, A., Schmidl, C., Krausgruber, T.,
#' Traxler, P., Klughammer, J., Schuster, L. C., Kuchler, A., Alpar, D.,
#' and Bock, C. (2017). Pooled crispr screening with single-cell transcriptome
#' readout. Nature Methods, 14, 297-301.
#' @references Dixit, A., Parnas, O., Li, B., Chen, J., Fulco, C. P.,
#' Jerby-Arnon, L., Marjanovic, N. D., Dionne, D., Burks, T., Raychowdhury, R.,
#' Adamson, B., Norman, T. M., Lander, E. S., Weissman, J. S., Friedman, N., and
#' Regev, A. (2016). Perturb-seq: Dissecting molecular circuits with scalable
#' single-cell rna profiling of pooled genetic screens. Cell, 167(7),
#' 1853-1866.e17.
#' @examples
#' data(app)
NA
#' Calculate responsibilities.
#'
#' This function calculates the responsibilities
#' of each component for all cells from the expected log distribution of the
#' hidden data.
#' @param x log odds for l cells and k components as a kxl matrix
#' @param affinity 0 for standard soft clustering, 1 for hard clustering
#' during inference (not recommended)
#' @param norm if TRUE normalises to probabilities (recommended)
#' @param logtype logarithm type of the data (e.g. 2 for log2 data or exp(1)
#' for natural)
#' @param mw mixture weights of the components
#' @param data data in log odds
#' @param complete if TRUE, complete data log likelihood is considered (for
#' very large data sets, e.g. 1000 cells and 1000 E-genes)
#' @return responsibilities as a kxl matrix (k components, l cells)
#' @author Martin Pirkl
#' @export
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' result <- mnem(data, k = 2, starts = 1)
#' resp <- getAffinity(result$probs, mw = result$mw, data = data)
getAffinity <- function(x, affinity = 0, norm = TRUE, logtype = 2, mw = NULL,
                        data = matrix(0, 2, ncol(x)), complete = FALSE) {
    if (is.null(mw)) { mw <- rep(1, nrow(x))/nrow(x) }
    if (affinity == 1) {
        if (complete) {
            y <- apply(y, 2, function(x) {
                xmax <- max(x)
                maxnum <- 2^1023
                shrinkfac <- log(maxnum)/log(logtype)
                x <- x - (xmax - shrinkfac)
                return(x)
            })
        }
        y <- logtype^x
        y <- y*mw
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
                if (complete & any(y > 0)) {
                    y <- apply(y, 2, function(x) {
                        xmax <- max(x)
                        maxnum <- 2^1023
                        shrinkfac <- log(maxnum)/log(logtype)
                        x <- x - (xmax - shrinkfac/length(x))
                        return(x)
                    })
                }
                y <- logtype^y
                y <- y*mw
                y <- y/colSums(y)[col(y)]
            }
        } else {
            y <- matrix(1, nrow(x), ncol(x))
        }
    }
    y[which(is.na(y) == TRUE)] <- 0
    return(y)
}
#' Calculate negative penalized log likelihood.
#'
#' This function calculates
#' a negative penalized log likelihood given a object of class mnem. This
#' penalized likelihood is based on the normal likelihood and penalizes
#' complexity of the mixture components (i.e. the networks).
#' @param x mnem object
#' @param man logical. manual data penalty, e.g. man=TRUE and pen=2 for an
#' approximation of the Akaike Information Criterion
#' @param degree different degree of penalty for complexity: positive entries
#' of transitively reduced phis or phi^r (degree=0), phi^r and mixture
#' components minus one k-1 (1), phi^r, k-1 and positive entries of thetas (2),
#' positive entries of transitively closed phis or phi^t, k-1 (3), phi^t, theta,
#' k-1 (4, default), all entries of phis, thetas and k-1 (5)
#' @param logtype logarithm type of the data (e.g. 2 for log2 data or exp(1)
#' for natural)
#' @param pen penalty weight for the data (e.g. pen=2 for approximate Akaike
#' Information Criterion)
#' @param useF use F (see publication) as complexity instead of phi and theta
#' @param Fnorm normalize complexity of F, i.e. if two components have the
#' same entry in F, it is only counted once
#' @author Martin Pirkl
#' @return penalized log likelihood
#' @export
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' pen <- numeric(3)
#' result <- list()
#' for (k in seq_len(2)) {
#'     result[[k]] <- mnem(data, k = k, starts = 1)
#'     pen[k] <- getIC(result[[k]])
#' }
#' print(pen)
getIC <- function(x, man = FALSE, degree = 4, logtype = 2, pen = 2,
                  useF = FALSE, Fnorm = FALSE) {
    n <- ncol(x$data)
    if (useF) {
        for (i in seq_len(length(x$comp))) {
            tmp <- mytc(x$comp[[i]]$phi)
            tmp2 <- matrix(0, nrow = nrow(tmp),
                           ncol = length(x$comp[[i]]$theta))
            tmp3 <- x$comp[[i]]$theta
            tmp3[which(tmp3 > nrow(tmp2))] <- 0
            tmp2[cbind(tmp3, seq_len(ncol(tmp2)))] <- 1
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
            for (i in seq_len(length(x$comp))) {
                tmp <- transitive.reduction(x$comp[[i]]$phi)
                if (degree > 2) {
                    tmp <- mytc(x$comp[[i]]$phi)
                }
                diag(tmp) <- 0
                fpar <- fpar + sum(tmp != 0)
            }
            if (degree > 1 & (degree != 3)) {
                fpar <- fpar + length(x$comp)*length(x$comp[[1]]$theta)
            }
        } else {
            fpar <-(length(x$comp[[1]]$phi)+
                    length(x$comp[[1]]$theta))*length(x$comp)
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
#' Learn the number of components K and optimize the mixture.
#'
#' High level function for learning the number of components k, if unknown.
#' @param D data with cells indexing the columns and features (E-genes)
#' indexing the rows
#' @param ks vector of number of components k to test
#' @param man logical. manual data penalty, e.g. man=TRUE and pen=2 for an
#' approximation of the Akaike Information Criterion
#' @param degree different degree of penalty for complexity: positive entries
#' of transitively reduced phis or phi^r (degree=0), phi^r and mixture
#' components minus one k-1 (1), phi^r, k-1 and positive entries of thetas (2),
#' positive entries of transitively closed phis or phi^t, k-1 (3), phi^t, theta,
#' k-1 (4, default), all entries of phis, thetas and k-1 (5)
#' @param logtype logarithm type of the data (e.g. 2 for log2 data or exp(1)
#' for natural)
#' @param pen penalty weight for the data (e.g. pen=2 for approximate Akaike
#' Information Criterion)
#' @param useF use F (see publication) as complexity instead of phi and theta
#' @param Fnorm normalize complexity of F, i.e. if two components have the
#' same entry in F, it is only counted once
#' @param ... additional parameters for the mnem main function
#' @author Martin Pirkl
#' @return list containing the result of the best k as an mnem object and the
#' raw and penalized log likelihoods
#' @export
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' result <- mnemk(data, ks = seq_len(2), starts = 1)
mnemk <- function(D, ks = seq_len(5), man = FALSE, degree = 4, logtype = 2,
                  pen = 2, useF = FALSE, Fnorm = FALSE, ...) {
    kres <- list()
    kic <- numeric(length(ks))
    kll <- numeric(length(ks))
    for (i in ks) {
        kres[[i]] <- mnem(D, k = i, ...)
        kic[i] <- getIC(kres[[i]], man = man, degree = degree,
                        logtype = logtype, pen = pen, useF = useF,
                        Fnorm = Fnorm)
        kll[i] <- kres[[i]]$ll
    }
    res <- list(best = kres[[which.min(kic)]], ics = kic, lls = kll)
    return(res)
}
#' Mixture NEMs - main function.
#'
#' This function simultaneously learns a mixture
#' of causal networks and clusters of a cell population from single cell
#' perturbation data (e.g. log odds of fold change) with a multi-trait
#' readout. E.g. Pooled CRISPR scRNA-Seq data (Perturb-Seq. Dixit et al., 2016,
#' Crop-Seq. Datlinger et al., 2017).
#' @param D data with cells indexing the columns and features (E-genes)
#' indexing the rows
#' @param inference inference method "em" for expectation maximization
#' @param search search method for single network inference "greedy",
#' "exhaustive" or "modules" (also possible: "small", which is greedy with
#' only one edge change per M-step to make for a smooth convergence)
#' @param phi a list of n lists of k networks for n starts of the EM and
#' k components
#' @param theta a list of n lists of k attachment vector for the E-genes
#' for n starts of the EM and k components
#' @param mw mixture weights; if NULL estimated or uniform
#' @param method "llr" for log ratios or foldchanges as input (see ratio)
#' @param parallel number of threads for parallelization of the number of
#' em runs
#' @param reduce logical - reduce search space for exhaustive search to
#' unique networks
#' @param runs number of runs for greedy search
#' @param starts number of starts for the em
#' @param type initialize with responsibilities either by "random",
#' "cluster" (each S-gene is clustered and the different S-gene clustered
#' differently combined for several starts),
#' "cluster2" (clustNEM is used to infer reasonable phis, which are then
#' used as a start for one EM run), "cluster3" (global clustering as a start),
#' or "networks" (initialize with random phis)
#' @param complete if TRUE, optimizes the expected complete log likelihood
#' of the model, otherwise the log likelihood of the observed data
#' @param p initial probabilities as a k (components) times l (cells) matrix
#' @param k number of components
#' @param kmax maximum number of components when k=NULL is inferred
#' @param verbose verbose output
#' @param max_iter maximum iteration, if likelihood does not converge
#' @param parallel2 if parallel=NULL, number of threads for single component
#' optimization
#' @param converged absolute distance for convergence between new and old log
#' likelihood; if set to -Inf, the EM stops if neither the phis nor thetas were
#' changed in the most recent iteration
#' @param redSpace space for "exhaustive" search
#' @param affinity 0 is default for soft clustering, 1 is for hard clustering
#' @param evolution logical. If TRUE components are penelized for being
#' different from each other.
#' @param lambda smoothness value for the prior put on the components, if
#' evolution set to TRUE
#' @param subtopoX hard prior on theta as a vector with entry i equal to j,
#' if E-gene i is attached to S-gene j
#' @param ratio logical, if true data is log ratios, if false foldchanges
#' @param logtype logarithm type of the data (e.g. 2 for log2 data or exp(1)
#' for natural)
#' @param domean average the data, when calculating a single NEM (speed
#' improvment)
#' @param modulesize max number of S-genes per module in module search
#' @param compress compress networks after search (warning: penelized
#' likelihood not interpretable)
#' @param increase if set to FALSE, the algorithm will not stop if the
#' likelihood decreases
#' @param fpfn numeric vector of length two with false positive and false
#' negative rates for discrete data
#' @param Rho perturbation matrix with dimensions nxl with n S-genes and
#' l samples; either as probabilities with the sum of probabilities for a
#' sample less or equal to 1 or discrete with 1s and 0s
#' @param ksel character vector of methods for the inference of k; can combine
#' "hc" (hierarchical clustering) or "kmeans" with "silhouette", "BIC" or "AIC";
#' can also include "cor" for correlation distance (preferred)
#' instead of euclidean
#' @author Martin Pirkl
#' @return object of class mnem
#' \item{comp}{list of the component with each component being
#' a list of the causal network phi and the E-gene attachment theta}
#' \item{data}{input data matrix}
#' \item{limits}{list of results for all indpendent searches}
#' \item{ll}{log likelihood of the best model}
#' \item{lls}{log likelihood ascent of the best model search}
#' \item{mw}{vector with mixture weights}
#' \item{probs}{kxl matrix containing the cell log likelihoods
#' of the model}
#' @export
#' @import
#' cluster
#' graph
#' Rgraphviz
#' naturalsort
#' snowfall
#' lattice
#' stats4
#' stats
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' result <- mnem(data, k = 2, starts = 1)
mnem <- function(D, inference = "em", search = "greedy", phi = NULL,
                 theta = NULL, mw = NULL, method = "llr",
                 parallel = NULL, reduce = FALSE, runs = 1, starts = 3,
                 type = "networks", complete = FALSE,
                 p = NULL, k = NULL, kmax = 10, verbose = FALSE,
                 max_iter = 100, parallel2 = NULL, converged = -Inf,
                 redSpace = NULL, affinity = 0, evolution = FALSE, lambda = 1,
                 subtopoX = NULL, ratio = TRUE, logtype = 2,
                 domean = TRUE, modulesize = 5, compress = FALSE,
                 increase = TRUE, fpfn = c(0.1, 0.1), Rho = NULL,
                 ksel = c("kmeans", "silhouette", "cor")) {
    if (length(grep("_", colnames(D))) > 0 & is.null(Rho)) {
        Rho <- getRho(D)
    }
    if (!is.null(k)) {
        if (k == 1) {
            type <- "random"
        }
    }
    if (nrow(D) > 10^3 & !complete) {
        print(paste0("The input data set consists of ", nrow(D), " E-genes ",
                     "and overall ", length(D), " data points. ",
                     "Consider option 'complete = TRUE', if the",
                     "observed log-likelihood reaches infinity, i.e. error."))
    }
    if (all(D %in% c(0,1))) { method <- "disc" }
    if (method %in% "disc") {
        D[which(D == 1)] <- log((1-fpfn[2])/fpfn[1])/log(logtype)
        D[which(D == 0)] <- log(fpfn[2]/(1-fpfn[1]))/log(logtype)
        method <- "llr"
    }
    if (reduce & search %in% "exhaustive" & is.null(redSpace)) {
        colnames(D) <- gsub("_.*", "", colnames(D))
        if (any(duplicated(colnames(D)) == TRUE)) {
            D <- D[, -which(duplicated(colnames(D)) == TRUE)]
        }
        redSpace <- nem(D,
                          search = "exhaustive", reduce = TRUE,
                          verbose = verbose, parallel = c(parallel, parallel2),
                          subtopo = subtopoX, ratio = ratio, domean = FALSE,
                          modulesize = modulesize, logtype = logtype,
                          modified = TRUE, Sgenes = Sgenes, Rho = Rho)$redSpace
    }
    if (!is.null(parallel)) { if (parallel == 1) { parallel <- NULL } }
    D.backup <- D
    D <- modData(D)
    Sgenes <- getSgenes(D)
    data <- D
    D <- NULL
    ## learn k:
    if (is.null(k) & is.null(phi) & is.null(p)) {
        if (length(grep("_", colnames(data))) > 0) {
            datak <- data
            colnames(datak) <- gsub("_", "", colnames(datak))
            datak <- modData(datak)
        } else {
            datak <- data
        }
        tmp <- learnk(datak, kmax = kmax, ksel = ksel, starts = starts,
                      verbose = verbose)
        k <- tmp$k
        if (verbose) {
            print(paste("components detected: ", k, sep = ""))
        }
        ks <- tmp$ks
    } else {
        if (!is.null(p)) {
            k <- nrow(p)
        } else if (!is.null(phi)) {
            k <- length(phi[[1]])
            starts <- length(phi)
        }
        ks <- rep(k, length(Sgenes))
    }
    if (is.null(mw)) {
        mw <- rep(1, k)/k
    }
    if (!is.null(k)) {
        if ((type %in% "networks" | type %in% "networks2") & is.null(phi)) {
            meanet <- nem(D = data, search = search, start = phi,
                            method = method,
                            parallel = parallel2, reduce = reduce,
                            runs = runs,
                            verbose = verbose, redSpace = redSpace,
                            ratio = ratio, domean = domean,
                            modulesize = modulesize, Rho = Rho,
                            logtype = logtype, modified = TRUE, Sgenes = Sgenes)
            if (type %in% "networks2") {
                linets <- TRUE
            } else {
                linets <- FALSE
            }
            init <- initComps(data, k, starts, verbose, meanet, linets = linets)
            probscl <- 0
        } else if (!is.null(phi)) {
            init <- phi
        } else {
            if (type %in% "cluster") {
                probscl <- initps(data, ks, k, starts = starts)
                init <- NULL
            } else if (type %in% c("cluster2", "cluster3")) {
                nem <- TRUE
                if (type %in% "cluster3") { nem <- FALSE }
                init <- clustNEM(data, nem = nem, k = k, starts = starts,
                                 search = search, method = method,
                                 parallel = parallel,
                                 runs = runs, verbose = verbose, ratio = ratio,
                                 domean = domean, modulesize = modulesize,
                                 Rho = Rho, logtype = logtype, modified = TRUE)
                if (type %in% "cluster2") {
                    mw <- init$mw
                    init <- list(A = unlist(init$comp, recursive = FALSE))
                } else {
                    p <- matrix(0.5, k, ncol(data))
                    for (i in seq_len(k)) {
                        p[i, which(init$cluster == i)] <- 1.5
                    }
                    init <- NULL
                }
                starts <- 1
            } else {
                probscl <- 0
                init <- NULL
            }
        }
    } else {
        probscl <- 0
        mw <- 1
    }
    if (!is.null(parallel)) { parallel2 <- NULL }
    res1 <- NULL
    if (starts <= 1) { parallel2 <- parallel; parallel <- NULL }
    if (!is.null(parallel)) { parallel2 <- NULL }
    if (k == 1) {
        if (!is.null(init)) {
            start <- phi[[1]][[1]]
        } else {
            start <- NULL
        }
        if (!is.null(theta)) {
            subtopoX <- theta[[1]][[1]]
        }
        if (!is.null(parallel) & is.null(parallel2)) { parallel2 <- parallel }
        if (verbose) {
            print("no evidence for sub populations or k set to 1")
        }
        limits <- list()
        limits[[1]] <- list()
        limits[[1]]$res <- list()
        limits[[1]]$res[[1]] <- nem(D = data, search = search, start = start,
                                      method = method,
                                      parallel = parallel2, reduce = reduce,
                                      runs = runs,
                                      verbose = verbose, redSpace = redSpace,
                                      ratio = ratio, domean = domean,
                                      modulesize = modulesize, Rho = Rho,
                                      logtype = logtype, modified = TRUE,
                                      Sgenes = Sgenes, subtopo = subtopoX)
        limits[[1]]$res[[1]]$D <- NULL
        res <- list()
        res[[1]] <- list()
        res[[1]]$adj <- limits[[1]]$res[[1]]$adj
        res[[1]]$subtopo <- limits[[1]]$res[[1]]$subtopo
        n <- ncol(res[[1]]$adj)
        mw <- 1
        probs <- matrix(0, k, ncol(data))
        probs <- getProbs(probs, k, data, res, method, n, affinity,
                          converged, subtopoX, ratio, mw = mw, fpfn = fpfn,
                          Rho = Rho, complete = complete)
        subtopoX <- probs$subtopoX
        limits[[1]]$ll <- probs$ll
        limits[[1]]$probs <- probs$probs
        limits[[1]]$mw <- 1
    } else {
        if (inference %in% "em") {
            if (!is.null(parallel)) {
                doMean <- doMean
                rowRanges <- matrixStats::rowRanges
                get.insertions <- get.ins.fast
                get.reversions <- get.rev.tc
                get.deletions <- get.del.tc
                naturalsort <- naturalsort::naturalsort
                eigenMapMatMult <- mnem:::eigenMapMatMult
                transClose_W <- mnem:::transClose_W
                maxCol_row <- mnem:::maxCol_row
                transClose_Ins <- mnem:::transClose_Ins
                sfInit(parallel = TRUE, cpus = parallel)
                sfLibrary(naturalsort)
                sfExport("modules", "mw", "ratio", "getSgeneN", "modData",
                         "sortAdj", "calcEvopen", "evolution",
                         "getSgenes",
                         "getLL", "getAffinity", "D", "nem",
                         "scoreAdj", "max_iter", "random_probs",
                         "verbose", "llrScore", "search", "redSpace",
                         "affinity", "getProbs", "probscl", "method",
                         "getRho", "doMean","transitive.closure",
                         "transitive.reduction", "get.insertions",
                         "mytc", "get.deletions",
                         "get.reversions", "nemEst", "rowRanges",
                         "eigenMapMatMult", "get.ins.fast",
                         "get.rev.tc", "get.del.tc",
                         "transClose_W","maxCol_row","transClose_Ins")
            }
            do_inits <- function(s) {
                if (!is.null(init)) {
                    res1 <- list()
                    for (i in seq_len(length(init[[s]]))) {
                        res1[[i]] <- list()
                        res1[[i]][["adj"]] <- init[[s]][[i]]
                        if (!is.null(theta)) {
                            res1[[i]][["subtopo"]] <- theta[[s]][[i]]
                        }
                    }
                    k <- length(res1)
                    n <- ncol(res1[[1]]$adj)
                    if (is.null(p)) {
                        probs <- matrix(0, k, ncol(data))
                        probs0 <- getProbs(probs, k, data, res1, method, n,
                                          affinity, converged, subtopoX,
                                          ratio, mw = mw, fpfn = fpfn,
                                          Rho = Rho, complete = complete)
                        subtopoX <- probs0$subtopoX
                        mw <- probs0$mw
                        ll <- probs0$ll
                        probs <- probs0$probs
                    } else {
                        p <- p*t(t(p)*apply(abs(data), 2, sum))
                        probs <- log(p)/log(logtype)
                    }
                }
                if (is.null(init)) {
                    if (is.null(p)) {
                        if (length(probscl) >= s & type %in% "cluster") {
                            probs <- probscl[[s]]
                        } else {
                            probs <- random_probs(k, data, logtype = logtype)
                        }
                    } else {
                        k <- nrow(p)
                        probs <- log(p)/log(logtype)
                    }
                }
                mw <- apply(getAffinity(probs, affinity = affinity, norm = TRUE,
                                        logtype = logtype, mw = mw, data = data,
                                        complete = complete), 1,sum)
                mw <- mw/sum(mw)
                if (any(is.na(mw))) { mw <- rep(1, k)/k }
                if (verbose) {
                    print(paste("start...", s))
                }
                limits <- list()
                ll <- getLL(probs, logtype = logtype, mw = mw, data = data,
                            complete = complete)
                llold <- bestll <- -Inf
                bestmw <- mw
                lls <- phievo <- thetaevo <- mwevo <- NULL
                count <- 0
                time0 <- TRUE
                probsold <- probs
                stop <- FALSE
                while (((((ll - llold > converged & abs(ll - llold) > 0 &
                           increase) |
                          (abs(ll - llold) > converged & !increase) &
                          count < max_iter)) | time0) & !stop) {
                              if (!time0) {
                                  if (ll - bestll > 0) {
                                      bestll <- ll; bestres <- res1
                                      bestprobs <- probs; bestmw <- mw
                                  }
                                  llold <- ll
                                  probsold <- probs
                              }
                              res <- list()
                              postprobs <- getAffinity(probs, affinity =
                                                                  affinity,
                                                       norm = TRUE, logtype =
                                                                        logtype,
                                                       mw = mw, data = data,
                                                       complete = complete)
                              edgechange <- 0
                              thetachange <- 0
                              for (i in seq_len(k)) {
                                  if (!is.null(res1) &
                                      !("modules" %in% search)) {
                                      start0 <- res1[[i]]$adj
                                  } else {
                                      start0 <- NULL
                                  }
                                  n <- getSgeneN(data)
                                  dataF <- matrix(0, nrow(data), n)
                                  colnames(dataF) <- seq_len(n)
                                  nozero <- which(postprobs[i, ] != 0)
                                  if (length(nozero) != 0) {
                                      if (length(nozero) == ncol(data)) {
                                          dataR <- data
                                          postprobsR <- postprobs[i, ]
                                          RhoR <- Rho
                                      } else {
                                          dataR <- cbind(data[, nozero,
                                                              drop = FALSE],
                                                         dataF)
                                          postprobsR <- c(postprobs[i, nozero],
                                                          rep(0, n))
                                          RhoR <- cbind(
                                              Rho[, nozero, drop = FALSE],
                                              diag(n))
                                      }
                                  } else {
                                      dataR <- dataF
                                      postprobsR <- rep(0, n)
                                      RhoR <- Rho
                                  }
                                  if (!is.null(theta)) {
                                      thetaX <- theta[[s]][[i]]
                                  } else {
                                      thetaX <- NULL
                                  }
                                  if (is.null(start0)) {
                                      res[[i]] <- nem(D = dataR,
                                                        weights = postprobsR,
                                                        search = search,
                                                        start = start0,
                                                        method = method,
                                                        parallel = parallel2,
                                                        reduce = reduce,
                                                        runs = runs,
                                                        verbose = verbose,
                                                        redSpace = redSpace,
                                                        ratio = ratio,
                                                        domean = domean,
                                                        modulesize = modulesize,
                                                        Rho = RhoR,
                                                        logtype = logtype,
                                                        modified = TRUE,
                                                        Sgenes = Sgenes,
                                                        subtopo = thetaX)
                                  } else {
                                      test01 <- list()
                                      test01scores <- numeric(2)
                                      for (j in seq_len(2)) {
                                          if (j == 2) {
                                              start1 <- start0
                                          } else {
                                              start1 <- start0*0
                                          }
                                          test01[[j]] <-
                                              nem(D = dataR,
                                                    weights =
                                                        postprobsR,
                                                    search = search,
                                                    start = start1,
                                                    method = method,
                                                    parallel =
                                                        parallel2,
                                                    reduce = reduce,
                                                    runs = runs,
                                                    verbose = verbose,
                                                    redSpace =
                                                        redSpace,
                                                    ratio = ratio,
                                                    domean = domean,
                                                    odulesize =
                                                        modulesize,
                                                    Rho = RhoR,
                                                    logtype = logtype,
                                                    modified = TRUE,
                                                    Sgenes = Sgenes,
                                                    subtopo = thetaX)
                                          test01scores[j] <-
                                              max(test01[[j]]$scores)
                                      }
                                      res[[i]] <-
                                          test01[[which.max(test01scores)]]
                                  }
                                  edgechange <- edgechange +
                                      sum(abs(res[[i]]$adj - res1[[i]]$adj))
                                  thetachange <- thetachange +
                                      sum(res[[i]]$subtopo != res1[[i]]$subtopo)
                                  res[[i]]$D <- NULL
                                  res[[i]]$subweights <- NULL

                              }
                              evopen <- 0
                              if (evolution) {
                                  res <- sortAdj(res)$res
                                  evopen <- calcEvopen(res)
                              }
                              res1 <- res
                              ## E-step:
                              n <- ncol(res[[1]]$adj)
                              probs0 <- list()
                              probs0$probs <- probsold
                              probs <- probsold
                              probs <- getProbs(probs, k, data, res, method,
                                                n,
                                                affinity, converged,
                                                subtopoX, ratio,
                                                mw = mw, fpfn = fpfn,
                                                Rho = Rho, complete = complete)
                              if (getLL(probs$probs, logtype = logtype,
                                        mw = mw, data = data,
                                        complete = complete) >
                                  getLL(probs0$probs, logtype = logtype,
                                        mw = mw, data = data,
                                        complete = complete)) {
                                  probs0 <- probs
                              }
                              subtopoX <- probs0$subtopoX
                              probs <- probs0$probs
                              modelsize <- n*n*k
                              datasize <- nrow(data)*ncol(data)*k
                              ll <- getLL(probs, logtype = logtype, mw = mw,
                                          data = data, complete = complete) +
                                  evopen*datasize*(modelsize^-1)
                              mwold <- mw
                              mw <- apply(getAffinity(probs,
                                                      affinity = affinity,
                                                      norm = TRUE,
                                                      logtype = logtype,
                                                      mw = mw,
                                                      data = data,
                                                      complete = complete),
                                          1, sum)
                              mw <- mw/sum(mw)
                              if(verbose) {
                                  print(paste("ll: ", ll, sep = ""))
                                  print(paste("changes in phi(s): ",
                                              edgechange, sep = ""))
                                  print(paste("changes in theta(s): ",
                                              thetachange, sep = ""))
                                  if (evolution) {
                                      print(paste("evolution penalty: ",
                                                  evopen, sep = ""))
                                  }
                              }
                              lls <- c(lls, ll)
                              phievo <- c(phievo, edgechange)
                              thetaevo <- c(thetaevo, thetachange)
                              mwevo <- c(mwevo, sum(abs(mw-mwold)))
                              count <- count + 1
                              if (is.infinite(converged) & phievo[count] == 0
                                  & thetaevo[count] == 0 & !time0) {
                                  stop <- TRUE
                              }
                              time0 <- FALSE
                          }
                if (ll - bestll > 0) {
                    bestll <- ll; bestres <- res1
                    bestprobs <- probs; bestmw <- mw
                }
                if (abs(ll - llold) > converged | llold > ll) {
                    if (verbose & (increase | max_iter <= count)) {
                        print("no convergence")
                    }
                }
                limits <- list()
                limits$probs <- bestprobs
                limits$res <- bestres
                limits$ll <- lls
                limits$k <- k
                limits$subtopo <- subtopoX
                limits$mw <- bestmw
                limits$phievo <- phievo
                limits$thetaevo <- thetaevo
                limits$mwevo <- mwevo
                return(limits)
            }
            if (!is.null(parallel)) {
                limits <- sfLapply(seq_len(starts), do_inits)
            } else {
                limits <- lapply(seq_len(starts), do_inits)
            }
            if (!is.null(parallel)) {
                sfStop()
            }
        }
        ## if (inference %in% "greedy") {
        ##     lls <- NULL
        ##     stop <- FALSE
        ##     probsold <- matrix(0, k, ncol(data))
        ##     while(!stop) {
        ##         evopen <- 0
        ##         for (i in seq_len(k)) {
        ##             probs0 <- list()
        ##             probs0$probs <- probsold
        ##             probs <- probsold
        ##             probs <- getProbs(probs, k, data, res, method,
        ##                               n,
        ##                               affinity, converged,
        ##                               subtopoX, ratio,
        ##                               mw = mw, fpfn = fpfn,
        ##                               Rho = Rho, complete = complete)
        ##             if (getLL(probs$probs, logtype = logtype,
        ##                       mw = mw, data = data, complete = complete) >
        ##                 getLL(probs0$probs, logtype = logtype, mw = mw,
        ##                       data = data, complete = complete)) {
        ##                 probs0 <- probs
        ##             }
        ##             subtopoX <- probs0$subtopoX
        ##             probs <- probs0$probs
        ##             modelsize <- n*n*k
        ##             datasize <- nrow(data)*ncol(data)*k
        ##             ll <- getLL(probs, logtype = logtype, mw = mw,
        ##                         data = data, complete = complete) +
        ##                 evopen*datasize*(modelsize^-1)
        ##         }
        ##     }
        ## }
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
        for (i in seq_len(length(best$res))) {
            if (i == 1) {
                unique[[1]] <- best$res[[1]]
                next()
            }
            count <- count + 1
            unique[[count]] <- best$res[[i]]
            added <- c(added, i)
            for (j in seq_len(length(unique)-1)) {
                if (all(unique[[count]]$adj - unique[[j]]$adj  == 0)) {
                    unique[[count]] <- NULL
                    count <- count - 1
                    added <- added[-length(added)]
                    break()
                }
            }
        }
        dead <- NULL
        for (i in seq_len(length(unique))) {
            dups <- NULL
            a <- mytc(unique[[i]]$adj)
            for (j in seq_len(length(unique))) {
                if (i %in% j | j %in% dead) { next() }
                b <- mytc(unique[[j]]$adj)
                dups <- c(dups, which(apply(abs(a - b), 1, sum) == 0))
            }
            if (all(seq_len(nrow(unique[[i]]$adj)) %in% dups)) {
                dead <- c(dead, i)
            }
        }
        if (is.null(dead)) {
            ulength <- seq_len(length(unique))
        } else {
            added <- added[-dead]
            ulength <- (seq_len(length(unique)))[-dead]
        }
        count <- 0
        unique2 <- list()
        for (i in ulength) {
            count <- count + 1
            unique2[[count]] <- unique[[i]]    }
        unique <- unique2
        probs <- best$probs[added, , drop = FALSE]
        colnames(probs) <- colnames(D.backup)
        postprobs <- getAffinity(probs, affinity = affinity, norm = TRUE,
                                 logtype = logtype, mw = mw, data = data,
                                 complete = complete)
        if (!is.null(dim(postprobs))) {
            lambda <- apply(postprobs, 1, sum)
            lambda0 <- lambda/sum(lambda)
            lambda <- best$mw
        } else {
            lambda <- 1
        }
        comp <- list()
        for (i in seq_len(length(unique))) {
            comp[[i]] <- list()
            comp[[i]]$phi <- unique[[i]]$adj
            comp[[i]]$theta <- unique[[i]]$subtopo
        }
    } else {
        probs <- best$probs[, , drop = FALSE]
        colnames(probs) <- colnames(D.backup)
        postprobs <- getAffinity(probs, affinity = affinity, norm = TRUE,
                                 logtype = logtype, mw = best$mw, data = data,
                                 complete = complete)
        if (!is.null(dim(postprobs))) {
            lambda <- apply(postprobs, 1, sum)
            lambda0 <- lambda/sum(lambda)
            lambda <- best$mw
            if (k == 1) { lambda <- 1 }
        } else {
            lambda <- 1
        }
        comp <- list()
        for (i in seq_len(length(best$res))) {
            comp[[i]] <- list()
            comp[[i]]$phi <- best$res[[i]]$adj
            comp[[i]]$theta <- best$res[[i]]$subtopo
        }
    }
    res <- list(limits = limits, comp = comp, data = D.backup, mw = lambda0,
                probs = probs, lls = best$ll, phievo = best$phievo,
                thetaevo = best$thetaevo, mwevo = best$mwevo,
                ll = getLL(probs, logtype = logtype, mw = lambda, data = data,
                           complete = complete), complete = complete, Rho = Rho)
    class(res) <- "mnem"
    return(res)
}
#' Bootstrap.
#'
#' Run bootstrap simulations on the components (phi)  of an object of
#' class mnem.
#' @param x mnem object
#' @param size size of the booststrap simulations
#' @param p percentage of samples (e.g. for 100 E-genes p=0.5 means
#' sampling 50)
#' @param logtype logarithm type of the data (e.g. 2 for log2 data or exp(1)
#' for natural)
#' @param complete if TRUE, complete data log likelihood is considered (for
#' very large data sets, e.g. 1000 cells and 1000 E-genes)
#' @param ... additional parameters for hte nem function
#' @author Martin Pirkl
#' @return returns bootstrap support for each edge in each component (phi); list
#' of adjacency matrices
#' @export
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' result <- mnem(data, k = 2, starts = 1)
#' boot <- bootstrap(result, size = 2)
bootstrap <- function(x, size = 1000, p = 1, logtype = 2, complete = FALSE,
                      ...) {
    data <- x$data
    post <- getAffinity(x$probs, logtype = logtype, mw = x$mw, data = data,
                        complete = complete)
    bootres <- list()
    for (i in seq_len(length(x$comp))) {
        dataR <- t(t(data)*post[i, ])
        bootres[[i]] <- x$comp[[i]]$phi*0
        for (j in seq_len(size)) {
            dataRR <- dataR[sample(seq_len(nrow(dataR)),
                                   ceiling(p*nrow(dataR)), replace = TRUE), ]
            res <- nem(dataRR, ...)
            bootres[[i]] <- bootres[[i]]+mytc(res$adj)
        }
        bootres[[i]] <- bootres[[i]]/size
        colnames(bootres[[i]]) <- rownames(bootres[[i]]) <-
            naturalsort(unique(colnames(data)))
    }
    class(bootres) <- "bootmnem"
    return(bootres)
}
#' Plot bootstrap mnem result.
#'
#' @param x bootmnem object
#' @param reduce if TRUE transitively reduces the graphs
#' @param ... additional parameters for the plotting function plotDNF
#' @author Martin Pirkl
#' @return visualization of bootstrap mnem result with Rgraphviz
#' @export
#' @method plot bootmnem
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' result <- mnem(data, k = 2, starts = 1)
#' boot <- bootstrap(result, size = 2)
#' plot(boot)
plot.bootmnem <- function(x, reduce = TRUE, ...) {
    dnfs <- freqs <- NULL
    for (i in seq_len(length(x))) {
        adj <- adj2 <- x[[i]]
        rownames(adj) <- colnames(adj) <- paste0(rownames(adj), "_", i)
        adj[which(adj <= 0.5)] <- 0
        if (reduce) {
            adj <- transitive.reduction(adj)*x[[i]]
        } else {
            adj <- x[[i]]
        }
        adj2 <- adj
        adj2[which(adj > 0.5)] <- 0
        bidi <- which(adj2+t(adj2) > 0.5)
        adj[which(adj <= 0.5)] <- 0
        adj[bidi] <- (adj2+t(adj2))[bidi]
        diag(adj) <- 0
        dnf <- adj2dnf(apply(adj, c(1,2), ceiling))
        dnf <- dnf[-seq_len(nrow(adj))]
        freq <- as.vector(t(adj))
        freq <- freq[which(freq != 0)]
        dnfs <- c(dnfs, dnf)
        freqs <- c(freqs, freq)
    }
    plotDnf(dnfs, edgelabel = freqs, ...)
}
#' Plot mnem result.
#'
#' @param x mnem object
#' @param oma outer margin
#' @param main main text
#' @param anno annotate cells by their perturbed gene
#' @param cexAnno text size of the cell annotations
#' @param scale scale cells to show relative and not absolute distances
#' @param global if TRUE clusters all cells, if FALSE clusters cells within
#' a component
#' @param egenes show egene attachments, i.e. number of E-genes
#' assigned to each S-gene
#' @param sep seperate clusters and not put them on top of each other
#' for better visualization
#' @param tsne if TRUE use tsne instead of pca
#' @param affinity use hard clustering if TRUE
#' @param logtype logarithm type of the data (e.g. 2 for log2 data or exp(1)
#' for natural)
#' @param cells show cell attachments, .i.e how many cells are assigned
#' to each S-gene
#' @param pch cell symbol
#' @param legend show legend
#' @param showdata show data if TRUE
#' @param bestCell show probability of best fitting cell for each S-gene
#' @param showprobs if TRUE, shows responsibilities for all cells and components
#' @param shownull if TRUE, shows the null node
#' @param ratio use log ratios (TRUE) or foldchanges (FALSE)
#' @param method "llr" for ratios
#' @param showweights if TRUE, shows mixture weights for all components
#' @param ... additional parameters
#' @author Martin Pirkl
#' @return visualization of mnem result with Rgraphviz
#' @export
#' @method plot mnem
#' @import
#' cluster
#' graph
#' Rgraphviz
#' tsne
#' stats
#' @importFrom graphics layout par text
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' result <- mnem(data, k = 2, starts = 1)
#' plot(result)
plot.mnem <- function(x, oma = c(3,1,1,3), main = "M&NEM", anno = TRUE,
                      cexAnno = 1, scale = NULL, global = TRUE, egenes = TRUE,
                      sep = FALSE, tsne = FALSE, affinity = 0, logtype = 2,
                      cells = TRUE, pch = ".", legend = FALSE, showdata = FALSE,
                      bestCell = TRUE, showprobs = FALSE, shownull = TRUE,
                      ratio = TRUE, method = "llr", showweights = TRUE, ...) {
    complete <- x$complete
    x2 <- x
    data <- x$data
    Rho <- x$Rho
    if (is.null(Rho)) {
        Rho <- getRho(data)
    }

    laymat <- rbind(seq_len(length(x$comp)+1),
                    c(length(x$comp)+2, rep(length(x$comp)+3,
                                            length(x$comp))))

    if (legend & !showdata) {
        laymat <- matrix(seq_len(length(x$comp)+1), nrow = 1)
    }
    if (!legend & !showdata) {
        laymat <- matrix(seq_len(length(x$comp)), nrow = 1)
    }
    if (!legend & showdata) {
        laymat <- rbind(seq_len(length(x$comp)),
                        c(length(x$comp)+1, rep(length(x$comp)+2,
                                                length(x$comp)-1)))
    }
    layout(laymat)

    par(oma=oma)
    if (legend) {
        plotDnf(c("Sgenes=Egenes", "Egenes=Cells", "Cells=Fit"),
                edgecol = rep("transparent", 3),
                nodeshape = list(Sgenes = "circle", Egenes = "box",
                                 Cells = "diamond", Fit = "circle"),
                nodelabel = list(Sgenes = "signaling\ngenes",
                                 Egenes = "effect\nreporters",
                                 Cells = "single\ncells",
                                 Fit = "highest\nresponsibility"),
                nodeheight = list(Sgenes = 2, Egenes = 0.5,
                                  Cells = 0.5, Fit = 0.5),
                nodewidth = list(Sgenes = 2, Egenes = 0.5,
                                 Cells = 0.5, Fit = 0.5),
                layout = "circo")
    }
    full <- x$comp[[1]]$phi
    mixnorm <- getAffinity(x$probs, affinity = affinity, norm = TRUE,
                           logtype = logtype, mw = x$mw, data = data,
                           complete = complete)
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
    Sgenes <- getSgenes(data)
    SgeneN <- getSgeneN(data)
    for (i in seq_len(length(x$comp))) {
        shared <- shared0 <- unique(colnames(mixnorm)[
            which(apply(mixnorm, 2,function(x)
                return(sum(x != 0))) != 1 & mixnorm[i, ] != 0)])
        net <- x$comp[[i]]$phi
        for (j in seq_len(SgeneN)) {
            colnames(net)[which(colnames(x$comp[[i]]$phi) %in% j)] <-
                rownames(net)[which(rownames(x$comp[[i]]$phi) %in% j)] <-
                Sgenes[j]
            shared[which(shared0 %in% j)] <- Sgenes[j]
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
                weights <- getAffinity(x$probs, affinity = affinity,
                                       norm = TRUE, logtype = logtype,
                                       mw = x$mw, data = data,
                                       complete = complete)
                subtopo <- scoreAdj(modData(data), x$comp[[i]]$phi,
                                    method = method, weights = weights[i, ],
                                    ratio = ratio, ...)$subtopo
            } else {
                subtopo <- x$comp[[i]]$theta
            }
            for (j in seq_len(SgeneN)) {
                tmpN <- paste("_9247E", j, sep = "_")
                enodes[[tmpN]] <- sum(subtopo == j)
                enodeshape[[tmpN]] <- "box"
                enodewidth[[tmpN]] <- 0.5
                enodeheight[[tmpN]] <- 0.5
                if (tmpN != 0) {
                    graph <- c(graph, paste(Sgenes[j], tmpN, sep = "="))
                }
            }
            if (shownull) {
                enull <- sum(subtopo == 0)
                enodes[["_9247Null"]] <- "NULL"
                enodeshape[["_9247Null"]] <- "circle"
                enodewidth[["_9247Null"]] <- 1
                enodeheight[["_9247Null"]] <- 1
                if (enull > 0) {
                    nulltargets <- enull
                } else {
                    nulltargets <- sum(subtopo == (SgeneN+1))
                }
                enodes[["_9247Nulltargets"]] <- nulltargets
                enodeshape[["_9247Nulltargets"]] <- "box"
                enodewidth[["_9247Nulltargets"]] <- 0.5
                enodeheight[["_9247Nulltargets"]] <- 0.5
                graph <- c(graph, "_9247Null=_9247Nulltargets")
            }
        } else {
            enodes <- NULL
            enodeshape <- NULL
            enodeheight <- NULL
            enodewidth <- NULL
        }
        if (cells) {
            datanorm <- modData(data)
            pnorm <- mixnorm
            if (is.null(dim(pnorm))) {
                pnorm <- matrix(pnorm, 1, length(pnorm))
            }
            cnodes <- list()
            cnodeshape <- list()
            cnodeheight <- list()
            cnodewidth <- list()
            for (j in seq_len(SgeneN)) {
                tmpN <- paste("__9247C", j, sep = "_")
                cnodes[[tmpN]] <-
                    sum(Rho[j, which(pnorm[i, ] == 1)] > 0)
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
                gam <- getAffinity(x$probs, affinity = affinity, norm = TRUE,
                                   logtype = logtype, mw = x$mw, data = data,
                                   complete = complete)
            } else {
                gam <- (logtype^x$probs)*x$mw
                gam <- gam/gam
            }
            for (bnode in Sgenes) {
                tmpN <- paste0("_9247bnode", bnode)
                graph <- c(graph, paste0(bnode, "=_9247bnode", bnode))
                bnodes[[tmpN]] <-
                    paste0(round(max(gam[
                        i, grep(bnode, colnames(gam))]), 2)*100, "%")
                bnodeheight[[tmpN]] <- 0.5
                bnodewidth[[tmpN]] <- 0.5
                bnodeshape[[tmpN]] <- "circle"
            }
        } else {
            bnodes <- bnodeshape <- bnodeheight <- bnodewidth <- NULL
        }
        edgecol <- c(rep("black", pathedges),
                     rep("grey", length(graph) - pathedges))
        if (showweights) {
            mainweights = paste("Cells: ",
                                realpct[i], "% (unique: ", unipct[i], "%)\n
Mixture weight: ", round(x$mw[i], 3)*100, "%", sep = "")
        } else {
            mainweights <- NULL
        }
        plotDnf(graph, main = mainweights, bordercol = i+1, width = 1,
                connected = FALSE,
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
            plotDnf(full)
        }
        if (nrow(mixnorm) > 1) {
            pnorm <- getAffinity(x$probs, affinity = affinity, norm = TRUE,
                                 logtype = logtype, mw = x$mw, data = data,
                                 complete = complete)
            if (is.null(dim(pnorm))) {
                pnorm <- matrix(pnorm, 1, length(pnorm))
            }
            pcols <- rep(1, ncol(pnorm))
            pcols[which(apply(mixnorm, 2, max) == 1)] <-
                apply(mixnorm[, which(apply(mixnorm, 2, max) == 1)]*
                      ((matrix(seq_len(nrow(mixnorm)), nrow(mixnorm),
                               ncol(mixnorm)))[,which(apply(mixnorm,
                                                            2, max) == 1),
                                               drop = FALSE]+1), 2, max)
            pcols[which(apply(mixnorm, 2, function(x)
                return(sum(x != 0))) != 1)] <- 1
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
            for (i in seq_len(SgeneN)) {
                rownames(pres$rotation)[
                    which(rownames(pres$rotation) %in% i)] <- Sgenes[i]
            }
            jittered <- pres$rotation[, seq_len(2)]
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
                    maxind <- seq_len(nrow(prtmp$rotation))
                } else {
                    maxind <- maxind0
                }
                if (sep) {
                    jittered[which(pcols == i), 1] <-
                        (prtmp$rotation[maxind, 1] -
                         mean(prtmp$rotation[maxind, 1]))*scale +
                        jittered[maxind0, 1]
                    jittered[which(pcols == i), 2] <-
                        (prtmp$rotation[maxind, 2] -
                         mean(prtmp$rotation[maxind, 2]))*scale +
                        jittered[maxind0, 2]
                } else {
                    if (!global) {
                        jittered[which(pcols == i), 1] <-
                            (prtmp$rotation[maxind, 1] -
                             mean(prtmp$rotation[maxind, 1]))
                        jittered[which(pcols == i), 2] <-
                            (prtmp$rotation[maxind, 2] -
                             mean(prtmp$rotation[maxind, 2]))
                    } else {
                        jittered[which(pcols == i), 1] <-
                            prtmp$rotation[maxind, 1]
                        jittered[which(pcols == i), 2] <-
                            prtmp$rotation[maxind, 2]
                    }
                }
            }
            if (anno) {
                cellnames <- apply(Rho, 2, function(x) {
                    y <- which(x > 0)
                    y <- paste(Sgenes[y], collapse = "_")
                    return(y)
                })
                plot(jittered, col = pcols, pch = "", main = main,
                     xlab = "pca1", ylab = "pca2")
                text(jittered[, 1],
                     jittered[, 2],
                     labels = cellnames,
                     srt = 45, pos = 1,
                     offset = 0, cex = cexAnno, col = pcols)
            } else {
                plot(jittered, col = pcols, main = main, pch = pch,
                     xlab = "pca1", ylab = "pca2")
            }
            unique <- unique(pres$rotation[which(pcols != 1), seq_len(2)])
            if (all(dim(unique) != 0) & !global) {
                for (i in seq_len(nrow(unique))) {
                    for (j in seq_len(nrow(unique))) {
                        if (i == j) { next() }
                        lines(unique[c(i,j), ], lty = 3, col = "grey")
                    }
                }
            }

        } else {

            prtmp <- prcomp(x$data)

            plot(prtmp$rotation[, seq_len(2)], pch = ".")

        }
    }
}
#' Cluster NEM.
#'
#' This function clusters the data and performs standard nem on each cluster.
#' @param data data of log ratios with cells in columns and features in rows
#' @param k number of clusters to check
#' @param cluster given clustering has to correspond to the columns of data
#' @param starts number of random starts for the kmeans algorithm
#' @param logtype logarithm type of the data
#' @param nem if FALSE only clusters the data
#' @param getprobspars list of parameters for the getProbs function
#' @param getaffinitypars list of parameters for the getAffinity function
#' @param Rho perturbation matrix with dimensions nxl with n S-genes and
#' l samples; either as probabilities with the sum of probabilities for a
#' sample less or equal to 1 or discrete with 1s and 0s
#' @param ... additional arguments for standard nem function
#' @author Martin Pirkl
#' @return family of nems; the first k list entries hold full information of
#' the standard nem search
#' \item{comp}{list of all adjacency matrices phi}
#' \item{mw}{vector of mixture weights}
#' \item{probs}{fake cell probabilities (see mw: mixture weights)}
#' @export
#' @import
#' stats
#' cluster
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- (sim$data - 0.5)/0.5
#' data <- data + rnorm(length(data), 0, 1)
#' resulst <- clustNEM(data, k = 2:3)
clustNEM <- function(data, k = 2:10, cluster = NULL, starts = 1, logtype = 2,
                     nem = TRUE, getprobspars = list(),
                     getaffinitypars = list(), Rho = NULL, ...) {
    data <- modData(data)
    if (is.null(cluster)) {
        smax <- 0
        K <- 1
        res <- NULL
        for (i in k) {
            d <- (1 - cor(data))/2
            d <- as.dist(d)
            kres <- kmeans(d, i, nstart = starts)
            sres <- silhouette(kres$cluster, d)
            if (mean(sres[, 3]) > smax) {
                Kres <- kres
                K <- i
                smax <- mean(sres[, 3])
            } else {
                break()
            }
        }
    } else {
        Kres <- list()
        Kres$cluster <- cluster
        K <- max(cluster)
    }
    res <- list()
    if (nem) {
        for (i in seq_len(K)) {
            if (sum(Kres$cluster == i) > 1 &
                length(unique(
                    colnames(data[,which(Kres$cluster == i)]))) > 1) {
                datar <- data[, which(Kres$cluster == i)]
                Rhor <- Rho[, which(Kres$cluster == i)]
                res[[i]] <- nem(datar, Rho = Rhor, ...)
            } else {
                res[[i]] <- list()
                res[[i]]$adj <- matrix(1, 1, 1)
                rownames(res[[i]]$adj) <- colnames(res[[i]]$adj) <-
                    unique(colnames(data[, which(Kres$cluster == i),
                                         drop = FALSE]))
            }
        }
        res$comp <- list()
        res$mw <- numeric(K)
        Sgenes <- getSgeneN(data)
        for (i in seq_len(K)) {
            res$comp[[i]] <- list()
            tmp <- res[[i]]$adj
            if (nrow(tmp) < Sgenes) {
                smiss <- unique(colnames(data)[-which(colnames(data)
                                                      %in% colnames(tmp))])
                tmp <-
                    rbind(cbind(tmp, matrix(0, nrow = nrow(tmp),
                                            ncol = length(smiss))),
                          matrix(0, nrow = length(smiss),
                                 ncol =ncol(tmp) + length(smiss)))
                colnames(tmp)[(dim(res[[i]]$adj)[1]+1):nrow(tmp)] <-
                    rownames(tmp)[(dim(res[[i]]$adj)[1]+1):nrow(tmp)] <- smiss
                tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
            }
            res$comp[[i]]$phi <- tmp
            res$mw[i] <- sum(Kres$cluster == i)/ncol(data)
            res[[i]]$adj <- tmp
        }
        n <- getSgeneN(data)
        probs <- matrix(0, K, ncol(data))
        probs <- do.call(getProbs, c(list(probs=probs, k=K, data=data,
                                          res=res, n=n, mw=res$mw, Rho=Rho),
                                     getprobspars))
        colnames(probs$probs) <- colnames(data)
        res$ll <- probs$ll
        res$probs <- probs$probs
        res$mw <- probs$mw
    }
    res$cluster <- Kres$cluster
    return(res)
}
#' Simulate data.
#'
#' This function simulates single cell data from a random
#' mixture of networks.
#' @param Sgenes number of Sgenes
#' @param Egenes number of Egenes
#' @param Nems number of components
#' @param reps number of replicates, if set (not realistic for cells)
#' @param mw mixture weights (has to be vector of length Nems)
#' @param evolution evolving and not purely random network, if set to TRUE
#' @param nCells number of cells
#' @param uninform number of uninformative Egenes
#' @param unitheta uniform theta, if TRUE
#' @param edgeprob edge probability, value between 0 and 1 for sparse or
#' dense networks
#' @param multi a vector with the percentages of cell with multiple
#' perturbations, e.g. c(0.2,0.1,0) for 20% double and 10% triple and
#' no quadruple knock-downs
#' @param subsample range to subsample data. 1 means the full simulated data is
#' used
#' @param scalefree if TRUE, graph is scale free
#' @param badCells number of cells, which are just noise and not connected
#' to the ground truth network
#' @param exactProb logical; if TRUE generates random network with exact
#' fraction of edges
#' @param ... additional parameters for the scale free network
#' sampler (see 'nem' package)
#' @author Martin Pirkl
#' @return simulation object with meta information and data
#' \item{Nem}{list of adjacency matrixes generatign the data}
#' \item{theta}{E-gene attachaments}
#' \item{data}{data matrix}
#' \item{index}{index for which Nem generated which cell
#' (data column)}
#' \item{mw}{vector of input mixture weights}
#' @export
#' @import
#' cluster
#' graph
#' Rgraphviz
#' tsne
#' @importFrom utils combn
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
simData <- function(Sgenes = 5, Egenes = 1,
                    Nems = 2, reps = NULL, mw = NULL, evolution = FALSE,
                    nCells = 1000, uninform = 0, unitheta = FALSE,
                    edgeprob = 0.25, multi = FALSE, subsample = 1,
                    scalefree = FALSE, badCells = 0, exactProb = TRUE, ...) {
    if (!is.null(mw) & Nems != length(mw)) {
        print(paste0("Vector of mixture weights 'mw' must be the length of the",
                     " number of komponents 'Nems'. Input 'Nems=", Nems,
                     " is overridden by the length ", length(mw),
                     " of 'mw'."))
        Nems <- length(mw)
    }
    if (multi[1] == FALSE) { multi <- rep(0, Sgenes-1) }
    if (length(multi) < Sgenes-1) {
        multi <- c(multi, rep(0, Sgenes - length(multi) - 1))
    }
    Nem <- list()
    data <- NULL
    index <- NULL
    theta <- list()
    for (i in seq_len(Nems)) {
        if (i == 1 | !evolution) {
            if (scalefree) {
                adj <- sampleRndNetwork(seq_len(Sgenes), DAG = TRUE, ...)
            } else {
                adj <- matrix(sample(c(0,1), Sgenes^2, replace = TRUE,
                                     prob = c(1-edgeprob, edgeprob)),
                              Sgenes, Sgenes)
                adj <- adj[order(apply(adj, 1, sum), decreasing = TRUE),
                           order(apply(adj, 2, sum), decreasing = FALSE)]
                adj[lower.tri(adj)] <- 0
                diag(adj) <- 1
                adj <- mytc(adj)
                while(sum(adj)-Sgenes >
                      floor((Sgenes*(Sgenes/2)-Sgenes/2)*edgeprob)
                      & exactProb) {
                          over <- sum(adj) - Sgenes - 
                              floor((Sgenes*(Sgenes/2) - Sgenes/2)*edgeprob)
                          adjtr <- transitive.reduction(adj)
                          adj[sample(which(adjtr == 1),
                                     min(over, sum(adjtr)))] <- 0
                      }
            }
            colnames(adj) <- rownames(adj) <- sample(seq_len(Sgenes),
                                                     Sgenes)
            adj <- adj[order(as.numeric(rownames(adj))),
                       order(as.numeric(colnames(adj)))]
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
            adj <- adj
        }
        Nem[[i]] <- transitive.reduction(adj)
    }
    for (i in seq_len(Nems)) {
        if (is.null(reps)) {
            reps2 <- ceiling(nCells/Sgenes)
        } else {
            reps2 <- reps
        }
        adj <- mytc(Nem[[i]])
        data_tmp <- t(adj)
        colntmp <- rep(seq_len(ncol(data_tmp)), reps2)
        data_tmp <- data_tmp[, rep(seq_len(ncol(data_tmp)), reps2)]
        colnames(data_tmp) <- colntmp
        if (!is.null(mw)) {
            tmpsamp <- sample(seq_len(ncol(data_tmp)),
                              ceiling(mw[i]*ncol(data_tmp)))
            data_tmp <- data_tmp[, tmpsamp, drop = FALSE]
        } else {
            mw <- rep(1/Nems, Nems)
            if (is.null(reps)) {
                data_tmp <- data_tmp[, seq_len(ceiling(nCells/Nems))]
            }
        }
        if (multi[1] == TRUE) {
            for (j in seq_len(Sgenes-1)) {
                data_fake <- matrix(0, 1, choose(Sgenes, j+1))
                colnames(data_fake) <- apply(combn(seq_len(Sgenes), j+1), 2,
                                             paste, collapse = "_")
                Rho <- getRho(data_fake)
                adj1 <- t(adj)%*%Rho
                adj1[which(adj1 > 1)] <- 1
                colnames(adj1) <- colnames(data_fake)
                data_tmp <- cbind(data_tmp, adj1)
                colnames(data_tmp)[-seq_len(ncol(data_tmp)-ncol(adj1))] <-
                    colnames(adj1)
            }
        } else {
            for (j in which(multi != 0)) {
                data_fake <- matrix(0, 1, choose(Sgenes, j+1))
                colnames(data_fake) <- apply(combn(seq_len(Sgenes), j+1), 2,
                                             paste, collapse = "_")
                Rho <- getRho(data_fake)
                adj1 <- t(adj)%*%Rho
                adj1[which(adj1 > 1)] <- 1
                colnames(adj1) <- colnames(data_fake)
                cells <- sample(which(!(seq_len(ncol(data_tmp)) %in%
                                        grep("_", colnames(data_tmp)))),
                                ceiling(multi[j]*ncol(data_tmp)),
                                replace = TRUE)
                knockdowns <- sample(seq_len(ncol(adj1)),
                                     ceiling(multi[j]*ncol(data_tmp)),
                                     replace = TRUE)
                data_tmp[, cells] <- adj1[, knockdowns]
                colnames(data_tmp)[cells] <- colnames(adj1)[knockdowns]
            }
        }
        index <- c(index, rep(i, ncol(data_tmp)))
        data_tmp <- data_tmp[rep(seq_len(Sgenes), each = Egenes), ,
                             drop = FALSE]
        if (!unitheta) {
            eorder <- sample(seq_len(nrow(data_tmp)), nrow(data_tmp))
            data_tmp <- data_tmp[eorder, ]
            theta[[i]] <- as.numeric(c(rownames(data_tmp), rep(0, uninform)))
        }
        data <- cbind(data, data_tmp)
    }
    if (subsample < 1) {
        subsample <- sample(seq_len(ncol(data)), ceiling(ncol(data)*subsample))
        data <- data[, subsample]
        index <- rep(seq_len(Nems), each = Sgenes*reps)[subsample]
    }
    if (uninform > 0) {
        data <- rbind(data, matrix(sample(c(0,1),
                                          ncol(data)*uninform, replace = TRUE),
                                   uninform, ncol(data)))
    }
    if (badCells > 0) {
        data <- cbind(data, matrix(0,nrow(data),badCells))
    }
    sim <- list(Nem = Nem, theta = theta, data = data, index = index, mw = mw)
    class(sim) <- "mnemsim"
    return(sim)
}
#' Plot simulated mixture.
#'
#' @param x mnemsim object
#' @param data noisy data matrix (optional)
#' @param logtype logarithm type of the data
#' @param fuzzypars list of parameters for the function fuzzyindex
#' @param ... additional parameters for the plotting function plotDNF
#' @author Martin Pirkl
#' @return visualization of simulated mixture with Rgraphviz
#' @export
#' @method plot mnemsim
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' plot(sim)
plot.mnemsim <- function(x, data = NULL, logtype = 2, fuzzypars = list(), ...) {
    noisymix <- TRUE
    if (is.null(data)) {
        data <- x$data
        data[which(data == 0)] <- -1
        noisymix2 <- ""
        noisymix <- FALSE
    }
    if (length(grep("_", colnames(data))) > 0) {
        Rho <- getRho(data)
    } else {
        Rho <- NULL
    }
    par(mfrow=c(1,length(x$mw)))
    fuzzypars <- c(list(x=x, data=data, logtype=logtype), fuzzypars)
    probs <- do.call(fuzzyindex, fuzzypars)
    mw <- probs$mw
    probs <- probs$probs
    mw2 <- unlist(apply(probs, 2, function(x) { return(which(x == max(x))) }))
    mw2 <- table(mw2)/ncol(probs)
    mw3 <- unlist(apply(probs, 2, function(x) {
        y <- which(x == max(x))
        if (length(y) > 1) {
            return(NULL)
        } else {
            return(y)
        }
    }))
    mw3 <- table(mw3)/ncol(probs)
    for (i in seq(length(x$mw))) {
        if (noisymix) {
            noisymix2 <- paste0("\n\nNoisy Mixture weight: ",
                               round(mw[i]*100, 2), "%")
        }
        tmp <- x$Nem[[i]]
        colnames(tmp) <- rownames(tmp) <- naturalsort(unique(colnames(data)))
        plotDnf(tmp, bordercol = i+1,
                main = paste0("Cells: ",
                              round(mw2[i]*100),
                              "% (unique ", round(mw3[i]*100), "%)",
                              "\n\nInput Mixture weight: ",
                              round(x$mw[i]*100), "%",
                              noisymix2), ...)
    }
}
#' Calculate fuzzy ground truth.
#'
#' Calculates responsibilities and mixture weights based on
#' the ground truth and noisy data.
#' @param x mnemsim object
#' @param data noisy data matrix
#' @param logtype logarithm type of the data
#' @param complete if TRUE, complete data log likelihood is considered (for
#' very large data sets, e.g. 1000 cells and 1000 E-genes)
#' @param ... additional parameters for the function getAffinity
#' @author Martin Pirkl
#' @return list with cell log odds mixture weights and log likelihood
#' @export
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' data <- sim$data
#' data[which(sim$data == 1)] <- rnorm(sum(sim$data == 1), 1, 1)
#' data[which(sim$data == 0)] <- rnorm(sum(sim$data == 0), -1, 1)
#' fuzzy <- fuzzyindex(sim, data)
fuzzyindex <- function(x, data, logtype = 2, complete = FALSE, ...) {
    data <- modData(data)
    k <- length(x$Nem)
    res <- list()
    for (i in seq_len(k)) {
        res[[i]] <- list()
        res[[i]]$adj <- x$Nem[[i]]
        res[[i]]$subtopo <- x$theta[[i]]
    }
    if (length(grep("_", colnames(data))) > 0) {
        Rho <- getRho(data)
    } else {
        Rho <- NULL
    }
    n <- getSgeneN(data)
    probs <- matrix(0, k, ncol(x$data))
    probs <- getProbs(probs, k, data, res, n = n, mw = x$mw, logtype = logtype,
                      Rho=Rho, complete = complete)
    ll <- probs$ll
    mw <- probs$mw
    colnames(probs$probs) <- colnames(data)
    probs2 <- do.call(getAffinity, c(list(x=probs$probs, logtype=logtype,
                                          mw=mw), ...))
    return(list(probs = probs2, mw = mw, ll = ll))
}
#' Accuracy for two phis.
#'
#' This function uses the hamming distance to calculate
#' an accuracy for two networks (phi).
#' @param a adjacency matrix (phi)
#' @param b adjacency matrix (phi)
#' @param diag if 1 includes diagonal in distance, if 0 not
#' @param symmetric comparing a to b is asymmetrical, if TRUE includes
#' comparison b to a
#' @author Martin Pirkl
#' @return normalized hamming accuracy for a and b
#' @export
#' @import
#' cluster
#' graph
#' Rgraphviz
#' tsne
#' @examples
#' sim <- simData(Sgenes = 3, Egenes = 2, Nems = 2, mw = c(0.4,0.6))
#' similarity <- hamSim(sim$Nem[[1]], sim$Nem[[2]])
hamSim <- function(a, b, diag = 1, symmetric = TRUE) {
    Sgenes <- unique(colnames(a))
    ham <- numeric(ncol(b))
    for (i in seq_len(ncol(b))) {
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
        for (i in seq_len(ncol(a))) {
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
#' Plot disjunctive normal form.
#'
#' This function visualizes a graph encoded as a disjunctive nromal form.
#' @param dnf Hyper-graph in disjunctive normal form,
#' e.g. c("A=B", "A=C+D", "E=!B") with the child on the left and the parents
#' on the right of the equation with "A=C+D" for A = C AND D. Alternatively,
#' dnf can be an adjacency matrix, which is converted on the fly to a
#' disjunctive normal form.
#' @param freq Frequency of hyper-edges which are placed on the edges.
#' @param stimuli Highlights vertices which can be stimulated.
#' @param signals Highlights vertices which regulate E-genes.
#' @param inhibitors Highlights vertices which can be inhibited.
#' @param connected If TRUE, only includes vertices which are connected to other
#' vertices.
#' @param CNOlist CNOlist object. Optional instead of stimuli, inhibitors or
#' signals. See package CellNOptR.
#' @param cex Global font size.
#' @param fontsize Vertice label size.
#' @param labelsize Edge label size.
#' @param type Different plot types. 2 for Rgraphviz and 1 for graph.
#' @param lwd Line width.
#' @param edgelwd Edgeline width.
#' @param legend 0 shows no legend. 1 shows legend as a graph. 2 shows legend
#' in a standard box.
#' @param x x coordinate of box legend.
#' @param y y coordinate of box legend.
#' @param xjust Justification of legend box left, right or center (-1,1,0).
#' @param yjust Justification of legend box top, bottom or middle (-1,1,0).
#' @param width Vertice width.
#' @param height Vertice height.
#' @param layout Graph layout. See graphvizCapabilities()$layoutTypes.
#' @param main Main title.
#' @param sub Subtitle.
#' @param cex.main Main title font size.
#' @param cex.sub Subtitle font size.
#' @param col.sub Font color of subtitle.
#' @param fontcolor Global font color.
#' @param nodestates Binary state of each vertice.
#' @param simulate Simulate stimulation and inhibition of a list of vertices.
#' E.g. simulate = list(stimuli = c("A", "B"), inhibitors = c("C", "D")).
#' @param edgecol Vector with colors for every edge of the graph
#' (not hyper-graph). E.g. an AND gate consists of three distinct edges.
#' @param labels Vector with labels for the edges.
#' @param labelcol Vector with label colors for the edges.
#' @param nodelabel List of vertices with labels as input.
#' E.g. labels = list(A="test", B="label for B").
#' @param nodecol List of vertices with colors as input.
#' @param bordercol List of vertices with colors as input.
#' @param nodeshape List of vertices with shapes (diamond, box, square,...).
#' @param verbose Verbose output.
#' @param edgestyle set the edge style like dashed, can be numerical
#' @param nodeheight List of vertices with height as input.
#' @param nodewidth List of vertices with width as input.
#' @param edgewidth Vector with edge widths.
#' @param lty Vector with edge styles (line, dotted,...).
#' @param hierarchy List with the hierarchy of the vertices.
#' E.g. list(top = c("A", "B"), bottom = c("C", "D")).
#' @param showall See "connected" above.
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
#' g <- c("!A+B+C=G", "C=G", "!D=G")
#' plotDnf(g)
plotDnf <- function(dnf = NULL, freq = NULL, stimuli = c(), signals = c(),
                    inhibitors = c(), connected = TRUE,  CNOlist = NULL,
                    cex = NULL, fontsize = NULL, labelsize = NULL, type = 2,
                    lwd = 1, edgelwd = 1, legend = 0, x = 0, y = 0, xjust = 0,
                    yjust = 0, width = 1, height = 1, layout = "dot", main = "",
                    sub = "",
                    cex.main = 1.5, cex.sub = 1, col.sub = "grey",
                    fontcolor = NULL, nodestates = NULL, simulate = NULL,
                    edgecol = NULL, labels = NULL,
                    labelcol = "blue", nodelabel = NULL, nodecol = NULL,
                    bordercol = NULL, nodeshape = NULL, verbose = FALSE,
                    edgestyle = NULL, nodeheight = NULL, nodewidth = NULL,
                    edgewidth = NULL, lty = NULL, hierarchy = NULL,
                    showall = FALSE, edgehead = NULL,
                    edgelabel = NULL, edgetail = NULL, bool = TRUE,
                    draw = TRUE, ...) {
    if (is.matrix(dnf)) {
        if (all(dnf == 0)) {
            diag(dnf) <- 1
        }
        edgelwd <- edgelabel <- dnf[which(dnf != 0)]
        edgelabel <- round(edgelabel, 3)
        edgelwd <- edgelwd - min(edgelwd) + 1
        dnf[which(dnf != 0)] <- 1
        dnf <- dnf2 <- adj2dnf(dnf)
        dnf <- dnf[grep("=", dnf)]
        if (length(dnf) == 0) {
            dnf <- dnf2
        }
    }

    if (!bool & length(grep("\\+", dnf)) > 0) {
        dnf <- dnf[-grep("\\+", dnf)]
    }

    graph <- dnf

    if (!is.null(hierarchy)) {
        if (!showall) {
            nodes <-
                unique(unlist(strsplit(unlist(strsplit(dnf, "=")), "\\+")))
            for (i in seq_len(length(hierarchy))) {
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
        for (i in seq_len(length(hierarchy)-1)) {
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
            edgecol <- c(rep("black", length(grep("\\+", graph))+
                                      length(unlist(strsplit(graph, "\\+")))),
                         rep("transparent", length(hgraph)))
            dnf2 <- dnf
            if (length(grep("\\+", dnf2)) > 0) {
                dnf2[-grep("\\+", dnf2)] <- gsub("=", "",
                                                 dnf2[-grep("\\+", dnf2)])
            } else {
                dnf2 <- gsub("=", "", dnf2)
            }
            edgecol[grep("!", unlist(strsplit(unlist(
                                  strsplit(dnf2,"\\+")), "=")))] <- "red"
        } else {
            if (length(edgecol) == 1) {
                edgecol <- c(rep(edgecol, length(grep("\\+", graph))+
                                          length(unlist(
                                              strsplit(graph, "\\+")))),
                             rep("transparent", length(hgraph)))
            } else {
                edgecol <- c(edgecol, rep("transparent", length(hgraph)))
            }
        }
    } else {
        if (is.null(lty) & !is.null(dnf)) {
            lty <- c(rep("solid", length(grep("\\+", graph))+
                                  length(unlist(strsplit(graph, "\\+")))))
        }
    }

    graph <- dnf

    dolegend <- FALSE
    if (is.null(dnf)) {
        dnf <- c("A=B")
        dolegend <- TRUE
    }

    if (!is.null(simulate)) {
        nodestates <- simulateDnf(graph, stimuli = simulate$stimuli,
                                  inhibitors = simulate$inhibitors)
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
        Vneg <-unique(c(c(unique(unlist(strsplit(
            unlist(strsplit(dnf, "=")), "\\+"))))))
        if (length(grep("\\+", dnf)) > 0) {
            Vneg <- c(Vneg, paste("and", seq_len(length(grep("\\+", dnf))),
                                  sep = ""))
        }
        V <- unique(gsub("!", "", Vneg))
        stimuli <- intersect(stimuli, V)
        signals <- intersect(signals, V)
        inhibitors <- intersect(inhibitors, V)
    } else {
        Vneg <- unique(c(c(unique(unlist(strsplit(
            unlist(strsplit(dnf, "=")), "\\+")))), stimuli,
            signals, inhibitors))
        if (length(grep("\\+", dnf)) > 0) {
            Vneg <- c(Vneg, paste("and", seq_len(length(grep("\\+", dnf))),
                                  sep = ""))
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
                Eneg[[paste("and", count, sep = "")]][["edges"]] <-
                    c(Eneg[[paste("and", count, sep = "")]][["edges"]],
                      which(Vneg %in% tmp[2]))
                E[[paste("and", count, sep = "")]][["edges"]] <-
                    c(E[[paste("and", count, sep = "")]][["edges"]],
                      which(V %in% tmp[2]))
                for (j in tmp2) {
                    Eneg[[j]][["edges"]] <-
                        c(Eneg[[j]][["edges"]],
                          which(Vneg %in%paste("and", count, sep = "")))
                    j <- gsub("!", "", j)
                    E[[j]][["edges"]] <-
                        c(E[[j]][["edges"]],
                          which(V %in% paste("and", count, sep = "")))
                }
            } else {
                Eneg[[tmp2]][["edges"]] <-
                    c(Eneg[[tmp2]][["edges"]], which(Vneg %in% tmp[2]))
                tmp2 <- gsub("!", "", tmp2)
                E[[tmp2]][["edges"]] <-
                    c(E[[tmp2]][["edges"]], which(V %in% tmp[2]))
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

    for (i in names(edgesneg)) {
        edgesnegnew <- c(edgesnegnew, edgesneg[[i]])
    }

    names(edgesnegnew) <- names(edgesneg)

    edgesneg <- edgesnegnew

    if (verbose) {
        print(paste("order of nodes: ", paste(names(nodes),
                                              collapse = ", "), sep = ""))
        print(paste("order of edges: ", paste(names(edges),
                                              collapse = ", "), sep = ""))
    }

    edges <- edgesneg
    names(edges) <- gsub("!", "", names(edges))

    for (i in seq_len(length(edges))) {
        edges[[i]]@from <- gsub("!", "", edges[[i]]@from)
    }
    nodeshape2 <- nodeshape
    nodeshape <- list()
    if (length(nodeshape2) == 1 & !(is.list(nodeshape2))) {
        for (i in seq_len(length(nodes))) {
            nodeshape[[nodes[[i]]@name]] <- nodeshape2
        }
    } else {
        nodeshape <- nodeshape2
    }

    for (i in seq_len(length(nodes))) {
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
            if (names(nodes)[i] %in% stimuli &
                is.null(nodeshape[[nodes[[i]]@name]])) {
                if (type == 2) {
                    nodes[[i]]@attrs$shape <- "diamond"
                } else {
                    nodes[[i]]@attrs$shape <- "box"
                }
            }
            if (names(nodes)[i] %in% signals &
                is.null(nodecol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$fillcolor <- "lightblue"
            }
            if (names(nodes)[i] %in% inhibitors &
                is.null(bordercol[[nodes[[i]]@name]])) {
                nodes[[i]]@attrs$color <- "red"
            }
        }
        if (!is.null(nodestates)) {
            if (sum(names(nodestates) %in% nodes[[i]]@name) == 1) {
                if (nodestates[which(names(nodestates) %in%
                                     nodes[[i]]@name)] == 0) {
                    if (is.null(nodecol[[nodes[[i]]@name]])) {
                        nodes[[i]]@attrs$fillcolor <- "white"
                    }
                    if (is.null(bordercol[[nodes[[i]]@name]])) {
                        nodes[[i]]@attrs$color <- "black"
                    }
                }
                if (nodestates[which(names(nodestates) %in%
                                     nodes[[i]]@name)] == 1) {
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
        for (i in seq_len(length(edges))) {
            edges[[i]]@attrs$fontsize <- as.character(labelsize)
            if (length(grep("and", names(edges)[i])) > 0) {
                tmp <- unlist(strsplit(names(edges)[i], "~"))
                k <- grep("and", tmp)
                inputN <- length(grep(tmp[k], edges))
                k <- as.numeric(gsub("and", "", tmp[k]))
                ## try to get the index of the correct edgecol:
                k2 <- grep("\\+", dnf)[k]
                if (grep("and", tmp) == 2) {
                    inputN2 <-
                        which(gsub("!", "",
                                   unlist(strsplit(gsub("=.*", "",dnf[k2]),
                                                   "\\+")))
                              %in% tmp[1])
                } else {
                    inputN2 <-
                        length(unlist(strsplit(gsub("=.*", "",
                                                    dnf[k2]), "\\+"))) + 1
                }
                if (k2 == 1) {
                    edgecolindex <- inputN2
                } else {
                    if (length(grep("\\+", graph[seq_len(k2-1)])) == 0) {
                        edgecolindex <- length(graph[seq_len(k2-1)]) + inputN2
                    } else {
                        edgecolindex <-
                            length(unlist(strsplit(dnf[seq_len(k2-1)],
                                                   "\\+"))) +
                            length(grep("\\+", dnf[seq_len(k2-1)])) + inputN2
                    }
                }
                ## end
                inputN2 <- grep(tmp[1],
                                unlist(strsplit(dnf[grep("\\+", dnf)[k]],
                                                "\\+")))-1
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
                        edges[[i]]@attrs$style <- edgestyle[edgecolindex]
                    }
                }
                if (!is.null(edgelabel)) {
                    if (!is.na(edgelabel[grep("\\+", dnf)[k]])) {
                        edges[[i]]@attrs$label <- edgelabel[edgecolindex]
                    }
                }
                if (length(grep("!", names(edgesneg)[i])) > 0) {
                    edges[[i]]@attrs$arrowhead <- "tee"
                    edges[[i]]@attrs$color <- "red"
                    if (!is.null(edgecol)) {
                        if (!is.na(edgecol[grep("\\+", dnf)[k]])) {
                            edges[[i]]@attrs$color <- edgecol[edgecolindex]
                        }
                    }
                    if (!is.null(edgehead)) {
                        if (!is.na(edgehead[grep("\\+", dnf)[k]])) {
                            edges[[i]]@attrs$arrowhead <- edgehead[edgecolindex]
                        }
                    }
                    if (!is.null(edgetail)) {
                        if (!is.na(edgetail[grep("\\+", dnf)[k]])) {
                            edges[[i]]@attrs$arrowtail <- edgetail[edgecolindex]
                        }
                    }
                } else {
                    edges[[i]]@attrs$arrowhead <- "open"
                    edges[[i]]@attrs$color <- "black"
                    if (gsub("and.*", "and", tmp[1]) %in% "and") {
                        if (!is.null(edgecol)) {
                            if (!is.na(edgecol[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$color <- edgecol[edgecolindex]
                            }
                        }
                        if (!is.null(edgehead)) {
                            if (!is.na(edgehead[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$arrowhead <-
                                    edgehead[edgecolindex]
                            }
                        }
                        if (!is.null(edgetail)) {
                            if (!is.na(edgetail[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$arrowtail <-
                                    edgetail[edgecolindex]
                            }
                        }
                    } else {
                        if (!is.null(edgecol)) {
                            if (!is.na(edgecol[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$color <- edgecol[edgecolindex]
                            }
                        }
                        if (!is.null(edgehead)) {
                            if (!is.na(edgehead[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$arrowhead <-
                                    edgehead[edgecolindex]
                            }
                        }
                        if (!is.null(edgetail)) {
                            if (!is.na(edgetail[grep("\\+", dnf)[k]])) {
                                edges[[i]]@attrs$arrowtail <-
                                    edgetail[edgecolindex]
                            }
                        }
                    }
                }
            } else {
                tmp <- unlist(strsplit(names(edges)[i], "~"))
                if (length(grep("!", names(edgesneg)[i])) == 0) {
                    k2 <- grep(paste("^", tmp[1], "=", tmp[2], "$", sep = ""),
                               dnf)
                } else {
                    k2 <- grep(paste("^!", tmp[1], "=", tmp[2], "$", sep = ""),
                               dnf)
                }
                if (k2 == 1) {
                    edgecolindex <- k2
                } else {
                    if (length(grep("\\+", dnf[seq_len(k2-1)])) == 0) {
                        edgecolindex <- k2
                    } else {
                        edgecolindex <-
                            length(unlist(strsplit(dnf[seq_len(k2-1)],
                                                   "\\+"))) +
                            length(grep("\\+", dnf[seq_len(k2-1)])) + 1
                    }
                }
                ## end
                if (length(grep("!", names(edgesneg)[i])) > 0) {
                    edges[[i]]@attrs$style <-
                        lty[grep(paste("^!",tmp[1], "=", tmp[2], "$", sep = ""),
                                 dnf)]
                    edges[[i]]@attrs$label <-
                        labels[grep(paste("^!", tmp[1], "=",
                                          tmp[2], "$", sep = ""), dnf)]
                    if (use.freq) {
                        edges[[i]]@attrs$weight <-
                            freq[grep(paste("^!", tmp[1], "=",
                                            tmp[2], "$", sep = ""), dnf)]
                        edges[[i]]@attrs$fontcolor <- "blue"
                    }
                    if (!is.null(edgewidth)) {
                        edges[[i]]@attrs$weight <-
                            edgewidth[grep(paste("^!", tmp[1], "=",
                                                 tmp[2], "$", sep = ""), dnf)]
                    }
                    if (!is.null(edgestyle)) {
                        if (!is.na(edgestyle[grep(paste("^!",
                                                        tmp[1], "=",
                                                        tmp[2], "$", sep = ""),
                                                  dnf)])) {
                            edges[[i]]@attrs$style <- edgestyle[edgecolindex]
                        }
                    }
                    if (!is.null(edgelabel)) {
                        if (!is.na(edgelabel[grep(paste("^!",
                                                        tmp[1], "=",
                                                        tmp[2], "$", sep = ""),
                                                  dnf)])) {
                            edges[[i]]@attrs$label <- edgelabel[edgecolindex]
                        }
                    }
                    edges[[i]]@attrs$arrowhead <- "tee"
                    edges[[i]]@attrs$color <- "red"
                    if (!is.null(edgecol)) {
                        if (!is.na(edgecol[grep(paste("^!",
                                                      tmp[1], "=",
                                                      tmp[2], "$", sep = ""),
                                                dnf)])) {
                            edges[[i]]@attrs$color <- edgecol[edgecolindex]
                        }
                    }
                    if (!is.null(edgehead)) {
                        if (!is.na(edgehead[grep(paste("^!",
                                                       tmp[1], "=",
                                                       tmp[2], "$", sep = ""),
                                                 dnf)])) {
                            edges[[i]]@attrs$arrowhead <- edgehead[edgecolindex]
                        }
                    }
                    if (!is.null(edgetail)) {
                        if (!is.na(edgetail[grep(paste("^!",
                                                       tmp[1], "=",
                                                       tmp[2], "$", sep = ""),
                                                 dnf)])) {
                            edges[[i]]@attrs$arrowtail <- edgetail[edgecolindex]
                        }
                    }
                } else {
                    edges[[i]]@attrs$style <-
                        lty[grep(paste("^", tmp[1], "=",
                                       tmp[2], "$", sep = ""), dnf)]
                    edges[[i]]@attrs$label <-
                        labels[grep(paste("^", tmp[1], "=",
                                          tmp[2], "$", sep = ""), dnf)]
                    if (use.freq) {
                        edges[[i]]@attrs$weight <-
                            freq[grep(paste("^",
                                            tmp[1], "=",
                                            tmp[2], "$", sep = ""), dnf)]
                        edges[[i]]@attrs$fontcolor <- "blue"
                    }
                    if (!is.null(edgewidth)) {
                        edges[[i]]@attrs$weight <-
                            edgewidth[grep(paste("^", tmp[1], "=",
                                                 tmp[2], "$", sep = ""), dnf)]
                    }
                    if (!is.null(edgestyle)) {
                        if (!is.na(edgestyle[grep(paste("^",
                                                        tmp[1], "=",
                                                        tmp[2], "$", sep = ""),
                                                  dnf)])) {
                            edges[[i]]@attrs$style <- edgestyle[edgecolindex]
                        }
                    }
                    if (!is.null(edgelabel)) {
                        if (!is.na(edgelabel[grep(paste("^",
                                                        tmp[1], "=",
                                                        tmp[2], "$", sep = ""),
                                                  dnf)])) {
                            edges[[i]]@attrs$label <- edgelabel[edgecolindex]
                        }
                    }
                    edges[[i]]@attrs$arrowhead <- "open"
                    edges[[i]]@attrs$color <- "black"
                    if (!is.null(edgecol)) {
                        if (!is.na(edgecol[grep(paste("^",
                                                      tmp[1], "=",
                                                      tmp[2], "$", sep = ""),
                                                dnf)])) {
                            edges[[i]]@attrs$color <- edgecol[edgecolindex]
                        }
                    }
                    if (!is.null(edgehead)) {
                        if (!is.na(edgehead[grep(paste("^",
                                                       tmp[1], "=",
                                                       tmp[2], "$", sep = ""),
                                                 dnf)])) {
                            edges[[i]]@attrs$arrowhead <- edgehead[edgecolindex]
                        }
                    }
                    if (!is.null(edgetail)) {
                        if (!is.na(edgetail[grep(paste("^",
                                                       tmp[1], "=",
                                                       tmp[2], "$", sep = ""),
                                                 dnf)])) {
                            edges[[i]]@attrs$arrowtail <- edgetail[edgecolindex]
                        }
                    }
                }
            }
        }
    }
    if (type == 1) {

        g2 <- agopen(name="boolean", nodes=nodes, recipEdges = "distinct",
                     edges=edges, edgeMode="undirected",
                     attrs=list(edge = list(),
                                graph = list(lwd = lwd),
                                node=list(lwd = lwd, fixedsize=FALSE)))

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
            for (i in seq_len(length(edges))) {
                if (length(edges[[i]]@attrs$style) == 0) {
                    edges[[i]]@attrs$style <- "solid"
                }
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
                    graph.trans <-
                        c(graph.trans,
                          paste(gsub("!", "", i), "~",
                                paste("and", and.count, sep = ""),
                                sep = ""))
                }
                graph.trans <-
                    c(graph.trans, paste(paste("and", and.count, sep = ""),
                                         "~", output, sep = ""))
            } else {
                put <- unlist(strsplit(i, "="))
                graph.trans <-
                    c(graph.trans, paste(gsub("!", "", put[1]), "~",
                                         put[2], sep = ""))
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

        for (i in seq_len(length(nodes))) {
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

        names(arrowfontsize) <- names(arrowheads) <- names(arrowtails) <-
            names(arrowcolors) <- names(arrowlwd) <- names(arrowlty) <-
            names(arrowlabels) <- names(edges)

        names(nodecolor) <- names(nodewidth) <- names(nodeheight) <-
            names(nodeshapes) <- names(nodecolors) <- names(nodes)

        if (length(unique(names(edges))) < length(names(edges))) {
            for (i in names(edges)[-which(duplicated(names(edges)) == TRUE)]) {
                getpos <- grep(paste("^", i, "$", sep = ""), names(edges))
                if (length(getpos) > 1) {
                    if (use.freq) {
                        if (arrowheads[getpos[1]] %in% "tee") {
                            arrowlabels[getpos[1]] <-
                                paste(paste(c("-", "+"),
                                            arrowlabels[getpos], sep = ""),
                                      collapse = "\n")
                        } else {
                            arrowlabels[getpos[1]] <-
                                paste(paste(c("+", "-"),
                                            arrowlabels[getpos], sep = ""),
                                      collapse = "\n")
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
                    arrowlwd[getpos[1]] <-
                        as.character(mean(as.numeric(arrowlwd[getpos])))
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
                g@nodes <- c("LEGEND:", "STIMULUS", "INHIBITOR", "SIGNAL",
                             "NOTHING", "active", "inactive")
                g@edgeL <- list()
                g@edgeData@data <- list()
            } else {
                start <- length(g@nodes) + 1
                g@nodes <- c(g@nodes, "LEGEND:", "STIMULUS", "INHIBITOR",
                             "SIGNAL", "NOTHING", "active", "inactive")
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
            g@edgeL[["active"]][["edges"]] <- c(as.integer(start+6),
                                                as.integer(start+4))
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
            arrowheads <- c(arrowheads, "LEGEND:~STIMULUS" = "none",
                            "STIMULUS~INHIBITOR" = "open",
                            "INHIBITOR~SIGNAL" = "tee",
                            "SIGNAL~NOTHING" = "odiamond",
                            "NOTHING~active" = "open",
                            "active~inactive" = "tee", "active~NOTHING" = "tee",
                            "inactive~active" = "open")
            arrowcolors <- c(arrowcolors, "LEGEND:~STIMULUS" = "transparent",
                             "STIMULUS~INHIBITOR" = "black",
                             "INHIBITOR~SIGNAL" = "red",
                             "SIGNAL~NOTHING" = "blue",
                             "NOTHING~active" = "black",
                             "active~inactive" = "red",
                             "active~NOTHING" = "red",
                             "inactive~active" = "black")
            arrowlabels <- c(arrowlabels, "LEGEND:~STIMULUS" = "",
                             "STIMULUS~INHIBITOR" = "    positive",
                             "INHIBITOR~SIGNAL" = "    negative",
                             "SIGNAL~NOTHING" =
                                 "    ambiguous\npositive\nnegative",
                             "NOTHING~active" = "    bidirectional\ndifferent",
                             "active~inactive" = "    bidirectional\ndifferent",
                             "active~NOTHING" = "", "inactive~active" = "")
            nodecolors <- c(nodecolors, "LEGEND:" = "white",
                            "STIMULUS" = "white", "INHIBITOR" = "white",
                            "SIGNAL" = "lightblue", "NOTHING" = "white",
                            "active" = "green", "inactive" = "white")
            nodeheight <- c(nodeheight, "LEGEND:" = 0,
                            "STIMULUS" = as.character(max(nodeheight)),
                            "INHIBITOR" = as.character(max(nodeheight)),
                            "SIGNAL" = as.character(max(nodeheight)),
                            "NOTHING" = as.character(max(nodeheight)),
                            "active" = as.character(max(nodeheight)),
                            "inactive" = as.character(max(nodeheight)))
            nodewidth <- c(nodewidth, "LEGEND:" = as.character(max(nodewidth)),
                           "STIMULUS" = as.character(max(nodewidth)),
                           "INHIBITOR" = as.character(max(nodewidth)),
                           "SIGNAL" = as.character(max(nodewidth)),
                           "NOTHING" = as.character(max(nodewidth)),
                           "active" = as.character(max(nodewidth)),
                           "inactive" = as.character(max(nodewidth)))
            if (type == 2) {
                nodeshapes <- c(nodeshapes, "LEGEND:" = "box",
                                "STIMULUS" = "diamond", "INHIBITOR" = "ellipse",
                                "SIGNAL" = "ellipse", "NOTHING" = "ellipse",
                                "active" = "ellipse", "inactive" = "ellipse")
            } else {
                nodeshapes <- c(nodeshapes, "LEGEND:" = "box",
                                "STIMULUS" = "box", "INHIBITOR" = "ellipse",
                                "SIGNAL" = "ellipse", "NOTHING" = "ellipse",
                                "active" = "ellipse", "inactive" = "ellipse")
            }
            nodecolor <- c(nodecolor, "LEGEND:" = "white", "STIMULUS" = "black",
                           "INHIBITOR" = "red", "SIGNAL" = "black",
                           "NOTHING" = "black", "active" = "black",
                           "inactive" = "black")
            dnf <- c(dnf, "NOTHING=active", "!active=NOTHING",
                     "!active=inactive", "inactive=active")
        }
        nodelabels <- names(nodecolor)
        names(nodelabels) <- nodelabels
        for (i in seq_len(length(nodelabel))) {
            nodelabels[which(names(nodelabels) %in% names(nodelabel)[i])] <-
                nodelabel[i]
        }
        g <- layoutGraph(g, edgeAttrs = list(arrowhead = arrowheads,
                                             color = arrowcolors,
                                             label = arrowlabels,
                                             arrowtail = arrowtails),
                         nodeAttrs = list(labels = nodelabels,
                                          color = nodecolor,
                                          height = nodeheight,
                                          width = nodewidth, shape = nodeshapes,
                                          fillcolor = nodecolors),
                         layoutType=layout)
        graph.par(list(graph=list(main = main, sub = sub, cex.main = cex.main,
                                  cex.sub = cex.sub, col.sub = col.sub),
                       edges=list(textCol = labelcol, lwd = edgelwd,
                                  fontsize = labelsize),
                       nodes=list(lwd = lwd,fontsize = fontsize, cex = cex)))
        edgeRenderInfo(g) <- list(lty = arrowlty, lwd = arrowlwd,
                                  label = arrowlabels)
        if (length(edges) > 0) {
            for (i in names(g@renderInfo@edges$direction)) {
                input <- unlist(strsplit(i, "~"))
                output <- input[2]
                input <- input[1]
                ambig <- FALSE
                if (paste("!", input, "=", output, sep = "") %in% dnf &
                    paste("", input, "=", output, sep = "") %in% dnf) {
                    ambig <- TRUE
                }
                if ((length(grep("and", i)) == 0 &
                     g@renderInfo@edges$direction[[i]] == "both") | ambig) {
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
                        if (paste("!", output, "=", input, sep = "") %in% dnf &
                            paste("", output, "=", input, sep = "") %in% dnf) {
                            g@renderInfo@edges$arrowtail[pos] <- "odiamond"
                        }
                        if (paste("!", input, "=", output, sep = "") %in% dnf &
                            paste("", input, "=", output, sep = "") %in% dnf) {
                            g@renderInfo@edges$arrowhead[pos] <- "odiamond"
                        }
                    }
                    if (is.null(edgecol)) {
                        if (g@renderInfo@edges$arrowtail[pos] == "open" &
                            g@renderInfo@edges$arrowhead[pos] == "open") {
                            g@renderInfo@edges$col[pos] <- "black"
                        }
                        if (g@renderInfo@edges$arrowtail[pos] == "tee" &
                            g@renderInfo@edges$arrowhead[pos] == "tee") {
                            g@renderInfo@edges$col[pos] <- "red"
                        }
                        if (g@renderInfo@edges$arrowtail[pos] !=
                            g@renderInfo@edges$arrowhead[pos]) {
                            g@renderInfo@edges$col[pos] <- "brown"
                        }
                        if (g@renderInfo@edges$arrowtail[pos] == "odiamond" |
                            g@renderInfo@edges$arrowhead[pos] == "odiamond") {
                            g@renderInfo@edges$col[pos] <- "blue"
                        }
                    } else {
                        if (is.null(edgecol)) { # is.na(edgecol[pos])
                            if (g@renderInfo@edges$arrowtail[pos] == "open" &
                                g@renderInfo@edges$arrowhead[pos] == "open") {
                                g@renderInfo@edges$col[pos] <- "black"
                            }
                            if (g@renderInfo@edges$arrowtail[pos] == "tee" &
                                g@renderInfo@edges$arrowhead[pos] == "tee") {
                                g@renderInfo@edges$col[pos] <- "red"
                            }
                            if (g@renderInfo@edges$arrowtail[pos] !=
                                g@renderInfo@edges$arrowhead[pos]) {
                                g@renderInfo@edges$col[pos] <- "brown"
                            }
                            if (g@renderInfo@edges$arrowtail[pos] ==
                                "odiamond" |
                                g@renderInfo@edges$arrowhead[pos] ==
                                "odiamond") {
                                g@renderInfo@edges$col[pos] <- "blue"
                            }
                        }
                    }
                }
            }
        }
        if (!is.null(simulate$draw)) {
            for (i in simulate$inhibitors) {
                g@nodes <- c(g@nodes, paste(i, "_inhibited", sep = ""))
                g@renderInfo@nodes$nodeX <-
                    c(g@renderInfo@nodes$nodeX,
                      g@renderInfo@nodes$nodeX[which(
                          names(g@renderInfo@nodes$nodeX) %in% i)])
                g@renderInfo@nodes$nodeY <-
                    c(g@renderInfo@nodes$nodeY,
                      g@renderInfo@nodes$nodeY[which(
                          names(g@renderInfo@nodes$nodeY) %in% i)])
                names(g@renderInfo@nodes$nodeX)[
                    length(g@renderInfo@nodes$nodeX)] <-
                    paste(i, "_inhibited", sep = "")
                names(g@renderInfo@nodes$nodeY)[
                    length(g@renderInfo@nodes$nodeY)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$labelX <-
                    c(g@renderInfo@nodes$labelX,
                      g@renderInfo@nodes$labelX[
                          which(names(g@renderInfo@nodes$labelX) %in% i)] + 200)
                g@renderInfo@nodes$labelY <-
                    c(g@renderInfo@nodes$labelY,
                      g@renderInfo@nodes$labelY[
                          which(names(g@renderInfo@nodes$labelY) %in% i)])
                names(g@renderInfo@nodes$labelX)[
                    length(g@renderInfo@nodes$labelX)] <-
                    paste(i, "_inhibited", sep = "")
                names(g@renderInfo@nodes$labelY)[
                    length(g@renderInfo@nodes$labelY)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$labelJust <-
                    c(g@renderInfo@nodes$labelJust, "n")
                names(g@renderInfo@nodes$labelJust)[
                    length(g@renderInfo@nodes$labelJust)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$col <-
                    c(g@renderInfo@nodes$col, "transparent")
                names(g@renderInfo@nodes$col)[
                    length(g@renderInfo@nodes$col)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$fill <-
                    c(g@renderInfo@nodes$fill, "transparent")
                names(g@renderInfo@nodes$fill)[
                    length(g@renderInfo@nodes$fill)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$shape <-
                    c(g@renderInfo@nodes$shape, "box")
                names(g@renderInfo@nodes$shape)[
                    length(g@renderInfo@nodes$shape)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$style <- c(g@renderInfo@nodes$style, "")
                names(g@renderInfo@nodes$style)[
                    length(g@renderInfo@nodes$style)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$height <- c(g@renderInfo@nodes$height, 2)
                names(g@renderInfo@nodes$height)[
                    length(g@renderInfo@nodes$height)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$rWidth <- c(g@renderInfo@nodes$rWidth, 1)
                names(g@renderInfo@nodes$rWidth)[
                    length(g@renderInfo@nodes$rWidth)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$lWidth <- c(g@renderInfo@nodes$lWidth, 1)
                names(g@renderInfo@nodes$lWidth)[
                    length(g@renderInfo@nodes$lWidth)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$label[[paste(i, "_inhibited", sep = "")]] <-
                    ""
                g@renderInfo@nodes$labelWidth <-
                    c(g@renderInfo@nodes$labelWidth,
                      g@renderInfo@nodes$labelWidth[
                          which(names(g@renderInfo@nodes$labelWidth) %in% i)])
                names(g@renderInfo@nodes$labelWidth)[
                    length(g@renderInfo@nodes$labelWidth)] <-
                    paste(i, "_inhibited", sep = "")
                tmp.name <- paste(i, "_inhibited", sep = "")
                g@edgeL[[tmp.name]] <- list()
                g@edgeL[[tmp.name]][["edges"]] <- which(g@nodes %in% i)
                g@edgeData@data[[paste(tmp.name, "|", i, sep = "")]] <- list()
                g@edgeData@data[[paste(tmp.name, "|", i, sep = "")]]$weight <- 1
                g@renderInfo@edges$enamesFrom <-
                    c(g@renderInfo@edges$enamesFrom, tmp.name)
                names(g@renderInfo@edges$enamesFrom)[
                    length(g@renderInfo@edges$enamesFrom)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$enamesTo <- c(g@renderInfo@edges$enamesTo, i)
                names(g@renderInfo@edges$enamesTo)[
                    length(g@renderInfo@edges$enamesTo)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelJust <-
                    c(g@renderInfo@edges$labelJust, NA)
                names(g@renderInfo@edges$labelJust)[
                    length(g@renderInfo@edges$labelJust)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelX <- c(g@renderInfo@edges$labelX, NA)
                names(g@renderInfo@edges$labelX)[
                    length(g@renderInfo@edges$labelX)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelY <- c(g@renderInfo@edges$labelY, NA)
                names(g@renderInfo@edges$labelY)[
                    length(g@renderInfo@edges$labelY)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelWidth <-
                    c(g@renderInfo@edges$labelWidth, NA)
                names(g@renderInfo@edges$labelWidth)[
                    length(g@renderInfo@edges$labelWidth)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$label <- c(g@renderInfo@edges$label, "")
                names(g@renderInfo@edges$label)[
                    length(g@renderInfo@edges$label)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$arrowhead <-
                    c(g@renderInfo@edges$arrowhead, "tee")
                names(g@renderInfo@edges$arrowhead)[
                    length(g@renderInfo@edges$arrowhead)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$arrowtail <-
                    c(g@renderInfo@edges$arrowtail, "odot")
                names(g@renderInfo@edges$arrowtail)[
                    length(g@renderInfo@edges$arrowtail)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$col <-
                    c(g@renderInfo@edges$col, "firebrick")
                names(g@renderInfo@edges$col)[
                    length(g@renderInfo@edges$col)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$lwd <-
                    c(g@renderInfo@edges$lwd, lwd[1]*1)
                names(g@renderInfo@edges$lwd)[
                    length(g@renderInfo@edges$lwd)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$lty <-
                    c(g@renderInfo@edges$lty, "solid")
                names(g@renderInfo@edges$lty)[
                    length(g@renderInfo@edges$lty)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$direction <-
                    c(g@renderInfo@edges$direction, "forward")
                names(g@renderInfo@edges$direction)[
                    length(g@renderInfo@edges$direction)] <-
                    paste(tmp.name, "~", i, sep = "")
                ## calculate splines
                tmp.splines <-
                    rep(g@renderInfo@nodes$labelY[
                        which(names(g@renderInfo@nodes$labelY) %in% i)], 8)
                tmp.splines[c(7,5,3,1)] <-
                    round(seq(g@renderInfo@nodes$labelX[
                        which(names(g@renderInfo@nodes$labelX) %in% i)] +
                        g@renderInfo@nodes$rWidth[[i]] + 10,
                        g@renderInfo@nodes$labelX[
                            which(names(g@renderInfo@nodes$labelX) %in% i)] +
                        200 - g@renderInfo@nodes$rWidth[[i]] - 10,
                        length.out = 4))
                tmp.splines[2] <- tmp.splines[2] + 50
                tmp.splines <- as.integer(tmp.splines)
                g@renderInfo@edges$splines[[paste(tmp.name, "~", i,
                                                  sep = "")]] <-
                    g@renderInfo@edges$splines[[1]]
                for (j in seq_len(4)) {
                    tmpname <- paste(tmp.name, "~", i, sep = "")
                    g@renderInfo@edges$splines[[tmpname]][[1]]@cPoints[[j]]@x <-
                        tmp.splines[c(1,3,5,7)][j]
                    tmpname <- paste(tmp.name, "~", i, sep = "")
                    g@renderInfo@edges$splines[[tmpname]][[1]]@cPoints[[j]]@y <-
                        tmp.splines[c(2,4,6,8)][j]
                }
            }
            for (i in simulate$stimuli) {
                ## add the stimulating node
                g@nodes <- c(g@nodes, paste(i, "_inhibited", sep = ""))
                g@renderInfo@nodes$nodeX <-
                    c(g@renderInfo@nodes$nodeX,
                      g@renderInfo@nodes$nodeX[
                          which(names(g@renderInfo@nodes$nodeX) %in% i)])
                g@renderInfo@nodes$nodeY <-
                    c(g@renderInfo@nodes$nodeY,
                      g@renderInfo@nodes$nodeY[
                          which(names(g@renderInfo@nodes$nodeY) %in% i)])
                names(g@renderInfo@nodes$nodeX)[
                    length(g@renderInfo@nodes$nodeX)] <-
                    paste(i, "_inhibited", sep = "")
                names(g@renderInfo@nodes$nodeY)[
                    length(g@renderInfo@nodes$nodeY)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$labelX <-
                    c(g@renderInfo@nodes$labelX,
                      g@renderInfo@nodes$labelX[
                          which(names(g@renderInfo@nodes$labelX) %in% i)] + 200)
                g@renderInfo@nodes$labelY <-
                    c(g@renderInfo@nodes$labelY,
                      g@renderInfo@nodes$labelY[
                          which(names(g@renderInfo@nodes$labelY) %in% i)])
                names(g@renderInfo@nodes$labelX)[
                    length(g@renderInfo@nodes$labelX)] <-
                    paste(i, "_inhibited", sep = "")
                names(g@renderInfo@nodes$labelY)[
                    length(g@renderInfo@nodes$labelY)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$labelJust <-
                    c(g@renderInfo@nodes$labelJust, "n")
                names(g@renderInfo@nodes$labelJust)[
                    length(g@renderInfo@nodes$labelJust)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$col <-
                    c(g@renderInfo@nodes$col, "transparent")
                names(g@renderInfo@nodes$col)[
                    length(g@renderInfo@nodes$col)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$fill <-
                    c(g@renderInfo@nodes$fill, "transparent")
                names(g@renderInfo@nodes$fill)[
                    length(g@renderInfo@nodes$fill)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$shape <-
                    c(g@renderInfo@nodes$shape, "box")
                names(g@renderInfo@nodes$shape)[
                    length(g@renderInfo@nodes$shape)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$style <-
                    c(g@renderInfo@nodes$style, "")
                names(g@renderInfo@nodes$style)[
                    length(g@renderInfo@nodes$style)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$height <-
                    c(g@renderInfo@nodes$height, 2)
                names(g@renderInfo@nodes$height)[
                    length(g@renderInfo@nodes$height)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$rWidth <- c(g@renderInfo@nodes$rWidth, 1)
                names(g@renderInfo@nodes$rWidth)[
                    length(g@renderInfo@nodes$rWidth)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$lWidth <- c(g@renderInfo@nodes$lWidth, 1)
                names(g@renderInfo@nodes$lWidth)[
                    length(g@renderInfo@nodes$lWidth)] <-
                    paste(i, "_inhibited", sep = "")
                g@renderInfo@nodes$label[[paste(i, "_inhibited", sep = "")]] <-
                    ""
                g@renderInfo@nodes$labelWidth <-
                    c(g@renderInfo@nodes$labelWidth,
                      g@renderInfo@nodes$labelWidth[
                          which(names(g@renderInfo@nodes$labelWidth) %in% i)])
                names(g@renderInfo@nodes$labelWidth)[
                    length(g@renderInfo@nodes$labelWidth)] <-
                    paste(i, "_inhibited", sep = "")
                tmp.name <- paste(i, "_inhibited", sep = "")
                g@edgeL[[tmp.name]] <- list()
                g@edgeL[[tmp.name]][["edges"]] <- which(g@nodes %in% i)
                g@edgeData@data[[paste(tmp.name, "|", i, sep = "")]] <- list()
                g@edgeData@data[[paste(tmp.name, "|", i, sep = "")]]$weight <- 1
                g@renderInfo@edges$enamesFrom <-
                    c(g@renderInfo@edges$enamesFrom, tmp.name)
                names(g@renderInfo@edges$enamesFrom)[
                    length(g@renderInfo@edges$enamesFrom)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$enamesTo <- c(g@renderInfo@edges$enamesTo, i)
                names(g@renderInfo@edges$enamesTo)[
                    length(g@renderInfo@edges$enamesTo)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelJust <-
                    c(g@renderInfo@edges$labelJust, NA)
                names(g@renderInfo@edges$labelJust)[
                    length(g@renderInfo@edges$labelJust)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelX <- c(g@renderInfo@edges$labelX, NA)
                names(g@renderInfo@edges$labelX)[
                    length(g@renderInfo@edges$labelX)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelY <- c(g@renderInfo@edges$labelY, NA)
                names(g@renderInfo@edges$labelY)[
                    length(g@renderInfo@edges$labelY)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$labelWidth <-
                    c(g@renderInfo@edges$labelWidth, NA)
                names(g@renderInfo@edges$labelWidth)[
                    length(g@renderInfo@edges$labelWidth)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$label <- c(g@renderInfo@edges$label, "")
                names(g@renderInfo@edges$label)[
                    length(g@renderInfo@edges$label)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$arrowhead <-
                    c(g@renderInfo@edges$arrowhead, "open")
                names(g@renderInfo@edges$arrowhead)[
                    length(g@renderInfo@edges$arrowhead)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$arrowtail <-
                    c(g@renderInfo@edges$arrowtail, "odot")
                names(g@renderInfo@edges$arrowtail)[
                    length(g@renderInfo@edges$arrowtail)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$col <-
                    c(g@renderInfo@edges$col, "limegreen")
                names(g@renderInfo@edges$col)[
                    length(g@renderInfo@edges$col)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$lwd <-
                    c(g@renderInfo@edges$lwd, lwd[1]*1)
                names(g@renderInfo@edges$lwd)[
                    length(g@renderInfo@edges$lwd)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$lty <-
                    c(g@renderInfo@edges$lty, "solid")
                names(g@renderInfo@edges$lty)[
                    length(g@renderInfo@edges$lty)] <-
                    paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$direction <-
                    c(g@renderInfo@edges$direction, "forward")
                names(g@renderInfo@edges$direction)[
                    length(g@renderInfo@edges$direction)] <-
                    paste(tmp.name, "~", i, sep = "")
                tmp.splines <-
                    rep(g@renderInfo@nodes$labelY[
                        which(names(g@renderInfo@nodes$labelY) %in% i)], 8)
                tmp.splines[c(7,5,3,1)] <-
                    round(seq(g@renderInfo@nodes$labelX[
                        which(names(g@renderInfo@nodes$labelX) %in% i)] +
                        g@renderInfo@nodes$rWidth[[i]] + 10,
                        g@renderInfo@nodes$labelX[
                            which(names(g@renderInfo@nodes$labelX) %in% i)] +
                        200 - g@renderInfo@nodes$rWidth[[i]] - 10,
                        length.out = 4))
                tmp.splines[2] <- tmp.splines[2] + 50
                tmp.splines <- as.integer(tmp.splines)
                tmpname <- paste(tmp.name, "~", i, sep = "")
                g@renderInfo@edges$splines[[tmpname]] <-
                    g@renderInfo@edges$splines[[1]]
                for (j in seq_len(4)) {
                    tmpname <- paste(tmp.name, "~", i, sep = "")
                    g@renderInfo@edges$splines[[tmpname]][[1]]@cPoints[[j]]@x <-
                        tmp.splines[c(1,3,5,7)][j]
                    tmpname <- paste(tmp.name, "~", i, sep = "")
                    g@renderInfo@edges$splines[[tmpname]][[1]]@cPoints[[j]]@y <-
                        tmp.splines[c(2,4,6,8)][j]
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
        legend(x = x, y = y,
               legend = c("signals are blue",
                          "stimuli are diamonds/boxes",
                          "inhibitors have a red border",
                          "positive regulation is green ->",
                          "negative regulation is red -|",
                          "ambiguous regulation is black -o"),
               fill = c("lightblue", "white", "red", "green", "red", "black"),
               col = c("lightblue", "white", "red", "green", "red", "black"),
               yjust = yjust, xjust = xjust)
    }
    return(g)
}
