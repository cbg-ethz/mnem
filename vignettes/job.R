
## define dataset

#dataset <- "cropseq"

#dataset <- "perturbseq_cc7d"

dataset <- "perturbseq_p7d"

#dataset <- "perturbseq_dc3"

#dataset <- "perturbseq_dc0"

cores <- 78

library(snowfall)

maxk <- 5

####################################### perturb-seq and crop-seq log likelihood ratio:

load(paste0(dataset, "_data.rda"))

if (length(grep("perturbseq", dataset)) == 0) {
    
    data <- t(t(data)/(colSums(data)/1000000))

    exprslvl <- apply(data, 1, median)
   
    data <- data[which(exprslvl > 0), ]

} else {

    data <- exp(data) - 1

    exprslvl <- apply(data, 1, median)
    
    data <- data[which(exprslvl > 0), ]

}

## data <- log(data + 1) # I don't know

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
        llrcol[j] <- log(ddens$y[ddist]/cdens$y[cdist])
    }
    return(llrcol)
}

sfInit(parallel = TRUE, cpus = cores)
llr <- sfLapply(1:nrow(data), densPar, data, C)
llr <- do.call("rbind", llr)
sfStop()

colnames(llr) <- colnames(data)

llr <- llr[, which(!(colnames(data) %in% ""))] # !!!

rownames(llr) <- rownames(data)

save(llr, file = paste0(dataset, "_llr.rda"))

print("llr done")

## stop()

############################################# do mnem:

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

done.genes <- NULL

allgenes <- sort(unique(colnames(llr)))

if (length(kadj) == 0) {

    kadj <- list()

    tmp <- matrix(0, length(unique(colnames(llr))), length(unique(colnames(llr))))

    colnames(tmp) <- rownames(tmp) <- sort(unique(colnames(llr)))

    kadj[[1]] <- tmp

}

for (i in 1:length(kadj)) {

    kegg <- kadj[[i]]

    for (j in 1:nrow(genes)) {
        colnames(kegg) <- gsub(paste("^", genes[j, 2], "$", sep = ""), genes[j, 1], colnames(kegg))
        rownames(kegg) <- gsub(paste("^", genes[j, 2], "$", sep = ""), genes[j, 1], rownames(kegg))
    }
    
    do.genes <- genes[which(genes[, 1] %in% colnames(kegg)), 1]

    done.genes <- c(done.genes, do.genes)

}

if (!all(allgenes %in% done.genes)) {

    do.genes <- allgenes[-which(allgenes %in% done.genes)]

    tmp <- matrix(0, length(do.genes), length(do.genes))

    colnames(tmp) <- rownames(tmp) <- sort(do.genes)

    kadj[[(length(kadj)+1)]] <- tmp

}

allres <- list()

for (i in 1:length(kadj)) {

    ## i <- 1

    kegg <- kadj[[i]]

    for (j in 1:nrow(genes)) {
        colnames(kegg) <- gsub(paste("^", genes[j, 2], "$", sep = ""), genes[j, 1], colnames(kegg))
        rownames(kegg) <- gsub(paste("^", genes[j, 2], "$", sep = ""), genes[j, 1], rownames(kegg))
    }
    
    do.genes <- genes[which(genes[, 1] %in% colnames(kegg)), 1]

    done.genes <- c(done.genes, do.genes)

    data <- llr[, which(colnames(llr) %in% do.genes)]

    badgenes <- "Tcrlibrary|^RP"
    
    if (length(grep(badgenes, rownames(data))) > 0) {
        data <- data[-grep(badgenes, rownames(data)), ]
    }
    
    data <- exp(data)

    data <- log2(data)

    n <- length(unique(colnames(data)))
    
    ## this is important: furhter reduce egenes?

    ## shrink to prevent infinite likelihood:
    
    data <- data/max(abs(data))

    bics <- rep(Inf, maxk)
    
    res <- list()

    for (k in 1:maxk) {

        ## k <- 2
        
        res[[k]] <- mnem(data, starts = min(n*k*10, cores), type = "random", parallel = cores, k = k, verbose = TRUE, converged = 10^-1)
        
        bics[k] <- log(nrow(data)*ncol(data))*(k*n*n + k*nrow(data)) - 2*max(res[[k]]$best$ll)

        if (k > 2) {
            if (bics[k] > bics[(k-1)]) { break() }
        }
    }

    res$bics <- bics

    allres[[i]] <- res

}

save(allres, file = paste0(dataset, "_mnem.rda"))

stop("mnem done")

###########################################

stop()

module load repo/grid
module load grid/grid

file=Yougottagetmeouttahereyoureamathguy.sh
echo "module load repo/grid" >> $file
echo "module load grid/grid" >> $file
echo "module load R/3.4.0" >> $file
echo "R --slave '--args arg1=1' < mnem/vignettes/job.R" >> $file

qsub -q mpi04-ht.q -pe make 78 $file

rm $file

##

## pdf("temp.pdf", width = 8, height = 4)

## par(mfrow=c(1,2))

## plotDnf(c("A=C", "B=C"), width = 1)

## plotDnf(c("A+B=C"), width = 1)

## dev.off()

##

library(nem)

fullsample <- numeric()

n <- 3

for (i in 1:1000) {

    print(i)

    ## tmp <- graph2adj(unifDAG(n, weighted=FALSE))

    ## tmp <- transitive.closure(graph2adj(unifDAG(n, weighted=FALSE)), mat = TRUE)

    ## gamma <- sample(c(0.1, 10), 1)
    
    ## tmp <- transitive.closure(sampleRndNetwork(1:n, scaleFree=TRUE, gamma=gamma, trans.close=TRUE, DAG=TRUE), mat = TRUE)
    
    ## tmp <- sampleRndNetwork(1:n, scaleFree=TRUE, gamma=gamma, trans.close=FALSE, DAG=TRUE)
    
    ## tmp <- matrix(sample(c(0,1), n*n, replace = TRUE), n, n)

    ## tmp[lower.tri(tmp)] <- 0

    ## tmp <- transitive.closure(tmp, mat = TRUE)

    ## rownames(tmp) <- colnames(tmp) <- sample(1:n, n)

    ## tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]

    fullsample <- rbind(fullsample, as.vector(tmp))

}

fullsample <- apply(fullsample, 1, paste, collapse = "")

par(mfrow=c(1,2))

hist(table(fullsample), breaks = 10)

hist(table(fullsample2), breaks = 10)

hist(as.numeric(fullsample), breaks = 10)

hist(as.numeric(fullsample2), breaks = 10)

##
