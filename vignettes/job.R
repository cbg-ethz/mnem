## ## tsting mpi01 49

## out <- 1 + 1

## save(out, file = "temp.rda")

## stop("test")

## ```r

## define dataset

#dataset <- "cropseq"

#dataset <- "perturbseq_cc7d"

#dataset <- "perturbseq_p7d"

dataset <- "perturbseq_dc3"

#dataset <- "perturbseq_dc0"

cores <- 78

library(snowfall)

maxk <- 10

dollr <- TRUE

getBIC <- function(x, AIC = FALSE, degree = 0) {
    fpar <- 0
    for (i in 1:length(x$best$res)) {
        tmp <- transitive.reduction(x$best$res[[i]]$adj)
        if (degree > 2) {
            ## tmp <- transitive.closure(x$best$res[[i]]$adj, mat = TRUE) # even though ffloops are determined by the trans reduction?
        }
        diag(tmp) <- 0
        fpar <- fpar + sum(tmp != 0) # edges
    }
    if (degree > 0) {
        fpar <- fpar + length(x$best$res) - 1 # mixture weights, which I don't set, but analytically estimate?
    }
    if (degree > 1 & (degree < 3 | degree > 3)) {
        fpar <- fpar + length(x$best$res)*nrow(x$data) # theta edges as free paras even though I use map?
    }
    n <- ncol(x$data)
    if (AIC) {
        bic <- 2*fpar - 2*max(x$best$ll)
    } else {
        bic <- log(n)*fpar - 2*max(x$best$ll)
    }
    return(bic)
}

####################################### perturb-seq and crop-seq log likelihood ratio:

if (dollr) {

    load(paste0(dataset, "_data.rda"))

    if (length(grep("perturbseq", dataset)) == 0) {

        exprslvl <- apply(data, 1, median)
        
        data <- t(t(data)/(colSums(data)/10000))

        data <- log2(data[which(exprslvl > 0), ] + 0.5)

        ## coldata <- matrix(colnames(data))

        ## colnames(coldata) <- "KO"

        ## coldata[which(coldata %in% "")] <- "CONTROL"

        ## colnames(data) <- paste0(colnames(data), "_", 1:ncol(data))
        
        ## dds <- DESeqDataSetFromMatrix(countData = data,
        ##                       colData = coldata,
        ##                       design = ~ KO)
        
        ## rld <- rlog(dds, blind=FALSE)
        
        ## exprslvl <- apply(data, 1, median)
        
        ## keep <- which(exprslvl > 0)

        ## data <- data[keep, ]

        ## design <- matrix(0, ncol(data), length(unique(colnames(data))))

        ## colnames(design) <- sort(unique(colnames(data)))

        ## rownames(design) <- colnames(data)

        ## for (i in 1:ncol(design)) {

        ##     design[which(rownames(design) %in% colnames(design)[i]), i] <- 1

        ## }

        ## data.voom <- voom(data, design, plot = TRUE)
        
    } else {

        data <- exp(data) - 1

        exprslvl <- apply(data, 1, median)
        
        data <- data[which(exprslvl > 0), ]

        if (length(grep("p7d|cc7d", dataset)) > 0) {

            data <- data[, -which(colnames(data) %in% "")]

            colnames(data)[grep("INTER", colnames(data))] <- ""

        }

        data <- log2(data + 0.5)

    }

    llr <- data*0

    C <- which(colnames(data) %in% "")

    densPar <- function(i, data, C) {
        llrcol <- numeric(ncol(data))
        for (j in which(!(colnames(data) %in% ""))) { # 1:ncol(data)) {
            gene <- colnames(data)[j]
            D <- which(colnames(data) %in% gene)
            cdens <- density(data[i, C], n = 1024, from = min(c(data[i, C], data[i, D])), to = max(c(data[i, C], data[i, D])))
            cdist <- which.min(abs(data[i, j] - cdens$x))
            ddens <- density(data[i, D], n = 1024, from = min(c(data[i, C], data[i, D])), to = max(c(data[i, C], data[i, D])))
            ddist <- which.min(abs(data[i, j] - ddens$x))
            llrcol[j] <- log2(ddens$y[ddist]/cdens$y[cdist])
        }
        return(llrcol)
    }
    
    ##library(ks) # why ?
    distrPar <- function(i, data, C) {
        llrcol <- numeric(ncol(data))
        for (j in which(!(colnames(data) %in% ""))) { # 1:ncol(data)) {
            gene <- colnames(data)[j]
            D <- which(colnames(data) %in% gene)
            cdistr <- ecdf(data[i, C])
            ##cdistr2 <- kcde(data[i, C])
            ##cdist <- which.min(abs(data[i, j] - cdistr$eval.points))
            ##cdist <- order(abs(data[i, j] - cdistr$eval.points))[1:2]
            ddistr <- ecdf(data[i, D])
            ##ddistr <- kcde(data[i, D])
            ##ddist <- which.min(abs(data[i, j] - ddistr$eval.points))
            ##ddist <- order(abs(data[i, j] - ddistr$eval.points))[1:2]
            llrcol[j] <- log2(min(ddistr(data[i, j]), 1 - ddistr(data[i, j]))/min(cdistr(data[i,j]), 1 - cdistr(data[i,j])))
            ##llrcol[j] <- log2(min(ddistr$estimate[ddist], 1 - ddistr$estimate[ddist])/min(cdistr$estimate[cdist], 1 - cdistr$estimate[cdist]))
            ##llrcol[j] <- log2(dist(ddistr$estimate[ddist])/dist(cdistr$estimate[cdist]))
        }
        return(llrcol)
    }

    sfInit(parallel = TRUE, cpus = cores)
    llr <- sfLapply(1:nrow(data), distrPar, data, C)
    llr <- do.call("rbind", llr)
    sfStop()

    llr[is.na(llr)] <- 0

    llr[is.infinite(llr)] <- max(llr[!is.infinite(llr)])

    colnames(llr) <- colnames(data)

    llr <- llr[, which(!(colnames(data) %in% ""))] # !!!

    rownames(llr) <- rownames(data)

    save(llr, file = paste0(dataset, "_llr.rda"))

    print("llr done")

}

## stop()

############### do on the aggregated graph:

source("mnem/R/mnems.r")
source("mnem/R/mnems_low.r")
library(cluster)
library(nem)

load(paste0(dataset, "_kegg.rda"))
load(paste0(dataset, "_llr.rda"))

if (length(grep("perturbseq", dataset)) == 0) {

    llr <- llr[, -grep("DHODH|MVD|TUBB", colnames(llr))]

}

llr <- t(apply(llr, 1, function(x) {
    x[is.infinite(x)] <- max(x[!is.infinite(x)])
    return(x)
}))

colnames(llr) <- toupper(colnames(llr))

if (length(grep("perturbseq", dataset)) == 0 | length(grep("dc", dataset)) > 0) {
    kegg <- kadjagg
    for (j in 1:nrow(genes)) {
        colnames(kegg) <- gsub(paste("^", genes[j, 2], "$", sep = ""), genes[j, 1], colnames(kegg))
        rownames(kegg) <- gsub(paste("^", genes[j, 2], "$", sep = ""), genes[j, 1], rownames(kegg))
    }
    lods <- llr[, which(colnames(llr) %in% colnames(kegg))]
    lods2 <- llr[, -which(colnames(llr) %in% colnames(kegg))]
} else {
    lods <- llr
}

badgenes <- "Tcrlibrary|^RP"

if (length(grep(badgenes, rownames(lods))) > 0) {
    lods <- lods[-grep(badgenes, rownames(lods)), ]
}

n <- length(unique(colnames(lods)))

lods <- lods/max(abs(lods))

bics <- rep(Inf, maxk)

res <- list()

for (k in 1:maxk) {
    
    ## k <- 2
    
    res[[k]] <- mnem(lods, starts = min(n*k*10, cores), type = "random", parallel = cores, k = k, verbose = TRUE, converged = 10^-1, search = "greedy")
    
    bics[k] <- getBIC(res[[k]])
    
    if (k > 2) {
        if (bics[k] > bics[(k-1)]) { break() }
    }

    save(res, file = paste0(dataset, "_mnem_agg.rda"))

}

res$bics <- bics

## no kegg genes:

if (length(grep(badgenes, rownames(lods2))) > 0) {
    lods2 <- lods2[-grep(badgenes, rownames(lods2)), ]
}

n <- length(unique(colnames(lods2)))

lods2 <- lods2/max(abs(lods2))

bics <- rep(Inf, maxk)

res2 <- list()

for (k in 1:maxk) {

    ## k <- 1
    
    res2[[k]] <- mnem(lods2, starts = min(n*k*10, cores), type = "random", parallel = cores, k = k, verbose = TRUE, converged = 10^-1, search = "greedy")
    
    bics[k] <- getBIC(res2[[k]])
    
    if (k > 2) {
        if (bics[k] > bics[(k-1)]) { break() }
    }

    save(res, res2, file = paste0(dataset, "_mnem_agg.rda"))

}

res2$bics <- bics

save(res, res2, file = paste0(dataset, "_mnem_agg.rda"))

stop("aggregate kegg done")

############################################# do mnem:

## source("mnem/R/mnems.r")
## source("mnem/R/mnems_low.r")
## library(cluster)
## library(nem)

## load(paste0(dataset, "_kegg.rda"))
## load(paste0(dataset, "_llr.rda"))

## if (length(grep("perturbseq", dataset)) == 0) {

##     llr <- llr[, -grep("DHODH|MVD|TUBB", colnames(llr))]

## }

## llr <- t(apply(llr, 1, function(x) {
##     x[is.infinite(x)] <- max(x[!is.infinite(x)])
##     return(x)
## }))

## colnames(llr) <- toupper(colnames(llr))

## done.genes <- NULL

## allgenes <- sort(unique(colnames(llr)))

## if (length(kadj) == 0) {

##     kadj <- list()

##     tmp <- matrix(0, length(unique(colnames(llr))), length(unique(colnames(llr))))

##     colnames(tmp) <- rownames(tmp) <- sort(unique(colnames(llr)))

##     kadj[[1]] <- tmp

## }

## for (i in 1:length(kadj)) {

##     kegg <- kadj[[i]]

##     for (j in 1:nrow(genes)) {
##         colnames(kegg) <- gsub(paste("^", genes[j, 2], "$", sep = ""), genes[j, 1], colnames(kegg))
##         rownames(kegg) <- gsub(paste("^", genes[j, 2], "$", sep = ""), genes[j, 1], rownames(kegg))
##     }
    
##     do.genes <- genes[which(genes[, 1] %in% colnames(kegg)), 1]

##     done.genes <- c(done.genes, do.genes)

## }

## if (!all(allgenes %in% done.genes)) {

##     do.genes <- allgenes[-which(allgenes %in% done.genes)]

##     tmp <- matrix(0, length(do.genes), length(do.genes))

##     colnames(tmp) <- rownames(tmp) <- sort(do.genes)

##     kadj[[(length(kadj)+1)]] <- tmp

## }

## allres <- list()

## for (i in 1:length(kadj)) {

##     ## i <- 1

##     kegg <- kadj[[i]]

##     for (j in 1:nrow(genes)) {
##         colnames(kegg) <- gsub(paste("^", genes[j, 2], "$", sep = ""), genes[j, 1], colnames(kegg))
##         rownames(kegg) <- gsub(paste("^", genes[j, 2], "$", sep = ""), genes[j, 1], rownames(kegg))
##     }
    
##     do.genes <- genes[which(genes[, 1] %in% colnames(kegg)), 1]

##     done.genes <- c(done.genes, do.genes)

##     data <- llr[, which(colnames(llr) %in% do.genes)]

##     badgenes <- "Tcrlibrary|^RP"
    
##     if (length(grep(badgenes, rownames(data))) > 0) {
##         data <- data[-grep(badgenes, rownames(data)), ]
##     }
    
##     data <- exp(data)

##     data <- log2(data)

##     n <- length(unique(colnames(data)))
    
##     ## this is important: furhter reduce egenes?

##     ## shrink to prevent infinite likelihood:
    
##     data <- data/max(abs(data))

##     bics <- rep(Inf, maxk)
    
##     res <- list()

##     for (k in 1:maxk) {

##         ## k <- 2
        
##         res[[k]] <- mnem(data, starts = min(n*k*10, cores), type = "random", parallel = cores, k = k, verbose = TRUE, converged = 10^-1, search = "modules")
        
##         bics[k] <- log(nrow(data)*ncol(data))*k - 2*max(res[[k]]$best$ll)

##         if (k > 2) {
##             if (bics[k] > bics[(k-1)]) { break() }
##         }
##     }

##     res$bics <- bics

##     allres[[i]] <- res

## }

## save(allres, file = paste0(dataset, "_mnem.rda"))

## stop("mnem done")

###########################################

stop()

rm *.sh.*

##

module load repo/grid
module load grid/grid

file=HolyshitItstheattackofEddieMunster.sh
echo "module load repo/grid" >> $file
echo "module load grid/grid" >> $file
echo "module load R/3.4.0" >> $file
echo "R --slave '--args arg1=1' < mnem/vignettes/job.R" >> $file

##qsub -q mpi01.q@bs-dsvr50.ethz.ch -pe make 1 $file

qsub -q mpi04-ht.q -pe make 78 $file

rm $file

##

## ```

bics <- rep(Inf, maxk)

res <- list()

for (k in 1:maxk) {

    ## k <- 2
    
    res[[k]] <- mnem(lods, starts = 15, type = "random", parallel = cores, k = k, verbose = TRUE, converged = 10^-1, search = "modules")
    
    bics[k] <- getBIC(lods, k, max(res[[k]]$best$ll))
    
    if (k > 2) {
        if (bics[k] > bics[(k-1)]) { break() }
    }

    save(res, file = "temp.rda")

}

res$bics <- bics

## get bics:

bics <- rep(Inf, maxk)

ll <- rep(Inf, maxk)

for (i in 1:length(res)) {

    bics[i] <- getBIC(res[[i]], degree = 4, AIC = TRUE)

    ll[i] <- max(res[[i]]$best$ll)

}

ll2 <- ll

ll <- ll - min(ll) + min(bics)

ll3 <- seq(min(bics), max(bics[!is.infinite(bics)]), length.out = 5)

par(mar=c(5,5,1,5))
plot(bics, type = "b", ylab = "AIC (red)", col = "red", xlab = "AIC and log likelihood as a function of number of components", yaxt = "n")
lines(ll, type = "b", col = "blue")
axis(4, ll3, round(ll3 + min(ll2) - min(bics)))
axis(2, ll3, round(ll3))
mtext("unnormalized log likelihood (blue)", side=4, line=3)

###### diff exprs:

