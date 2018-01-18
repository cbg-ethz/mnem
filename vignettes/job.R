## ## tsting mpi01 49

## out <- 1 + 1

## save(out, file = "temp.rda")

## stop("test")

## ```r

source("mnem/R/mnems.r")
source("mnem/R/mnems_low.r")
library(cluster)
library(nem)
library(Rgraphviz)

## define dataset

#dataset <- "cropseq"

#dataset <- "perturbseq_cc7d"

#dataset <- "perturbseq_p7d"

#dataset <- "perturbseq_dc3"

#dataset <- "perturbseq_dc0"

args <- commandArgs()

dataset <- gsub("dataset=", "", args[grep("dataset=", args)])

parallel <- gsub("cores=", "", args[grep("cores=", args)])

starts <- gsub("starts=", "", args[grep("starts=", args)])

run <- gsub("run=", "", args[grep("run=", args)])

dollr <- as.numeric(gsub("dollr=", "", args[grep("dollr=", args)]))

dobig <- as.numeric(gsub("dobig=", "", args[grep("dobig=", args)]))

dosmall <- as.numeric(gsub("dosmall=", "", args[grep("dosmall=", args)]))

addendum <- paste0("_run", run)

library(snowfall)

maxk <- 5

print(dataset)
print(dobig)
print(dosmall)
print(dollr)

####################################### perturb-seq and crop-seq log likelihood ratio:

if (dollr) {

    load(paste0(dataset, "_data.rda"))

    if (length(grep("perturbseq", dataset)) == 0) {

        exprslvl <- apply(data, 1, median)
        
        data <- data[which(exprslvl > 0), ]

        data <- t(t(data)/(colSums(data)/10000))

        data <- log2(data + 0.5)
    
    } else {

        data <- exp(data) - 1

        exprslvl <- apply(data, 1, median)
        
        data <- data[which(exprslvl > 0), ]

        if (length(grep("p7d|cc7d", dataset)) > 0) {

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

    sfInit(parallel = TRUE, cpus = parallel)
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

########### do small set of genes:

if (dosmall) {

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

    if (length(grep("perturbseq", dataset)) == 0) {

        cropgenes <- c("LCK", "ZAP70", "PTPN6", "DOK2", "PTPN11", "EGR3", "LAT")
        
        lods <- llr[, which(colnames(llr) %in% cropgenes)]
        
    } else {
        
        lods <- llr
        
    }

    badgenes <- "Tcrlibrary"#|^RP" # really remove ribosomal protein genes?
    badgenes <- grep(badgenes, rownames(lods))

    if (length(badgenes) > 0) {
        lods <- lods[-badgenes, ]
    }

    sdev <- apply(lods, 1, sd)

    lods <- lods[which(sdev > sd(lods)), ] ## this really?

    n <- length(unique(colnames(lods)))

    lods <- lods#/(sum(lods)/1023) # max(abs(lods))

    bics <- rep(Inf, maxk)

    res <- list()
    
    for (k in 1:maxk) {

        ## k <- 2
        
        res[[k]] <- mnem(lods, starts = starts, type = "random", parallel = parallel, k = k, verbose = TRUE, converged = 10^-1, search = "modules")
        
        ## bics[k] <- getIC(res[[k]], AIC = TRUE, degree = 0) 
        
        ## if (k > 2) {
        ##     if (bics[k] > bics[(k-1)]) { break() }
        ## }

        res[[k]]$data <- NULL
        
        save(res, file = paste0("cropseq/", dataset, "_mnem_small", addendum, ".rda"))
        
    }
    
    save(res, file = paste0("cropseq/", dataset, "_mnem_small", addendum, ".rda"))

    stop("small set done")

}

############### do on the aggregated graph:

if (dobig) {

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

    badgenes <- "Tcrlibrary"#|^RP""
    badgenes <- grep(badgenes, rownames(lods))

    if (length(badgenes) > 0) {
        lods <- lods[-badgenes, ]
    }

    n <- length(unique(colnames(lods)))

    lods <- lods # max(abs(lods)) normalising the max to one does not seem to be enough... weird

    bics <- rep(Inf, maxk)

    res2 <- list()

    res3 <- list()

    for (i in 1:2) {

        ## i <- 1

        if (i == 2 & length(grep("perturbseq", dataset)) > 0) { break() }

        res <- list()

        if (i == 2) {

            if (length(grep(badgenes, rownames(lods2))) > 0) {
                lods2 <- lods2[-grep(badgenes, rownames(lods2)), ]
            }

            n <- length(unique(colnames(lods2)))

            lods2 <- lods2/(sum(lods2)/1023) # max(abs(lods2))

            lods <- lods2
        }
        
        for (k in 1:maxk) {
            
            ## k <- 1
            
            res[[k]] <- mnem(lods, starts = starts, type = "random", parallel = parallel, k = k, verbose = TRUE, converged = 10^-1, search = "modules", modulesize = 5)

            ## bics[k] <- getIC(res[[k]], AIC = TRUE, degree = 0) 
            
            ## if (k > 2) {
            ##     if (bics[k] > bics[(k-1)]) { break() }
            ## }

            res[[k]]$data <- NULL

            save(res, res3, file = paste0("cropseq/", dataset, "_mnem_big", addendum, ".rda"))

        }

        if (i == 1) {
            res3 <- res
        } else {
            res2 <- res
            res <- res3
        }
        
    }

    save(res, res2, file = paste0("cropseq/", dataset, "_mnem_big", addendum, ".rda"))

    stop("big kegg done")

}

###########################################

stop()

rm *.sh.*
    
##
    
##qsub -q mpi01.q@bs-dsvr50.ethz.ch -pe make 1 $file

## on the grid:
    
## several starts on one core each:

rm *.sh.*

##

module load repo/grid
module load grid/grid

for i in `seq 1 100`; do
    cores=1
    file=${i}.sh
    echo "module load repo/grid" >> $file
    echo "module load grid/grid" >> $file
    echo "module load R/3.4.0" >> $file
    echo "R --slave --args 'dataset=perturbseq_p7d' 'cores=$cores' 'run=$i' 'starts=1' 'dollr=0' 'dobig=1' 'dosmall=0' < mnem/vignettes/job.R" >> $file
    qsub -q mpi01.q -pe make $cores $file
    rm $file
done

## on euler:

## several starts on one core each:

R
q("no")
module load r/3.4.0

for i in `seq 2 100`; do
    cores=1
    file=Withsomeluck${i}.sh
    echo "R --slave --args 'dataset=cropseq' 'cores=$cores' 'run=$i' 'starts=1' 'dollr=0' 'dobig=1' 'dosmall=0' < mnem/vignettes/job.R" >> $file
bsub -M 10000 -q normal.24h -n 1 -e logs/${file}_${i}_error.txt -o logs/${file}_${i}_output.txt < $file
    rm $file
done

## get the best for several on one core each:

setwd("~/Mount/Euler")

maxk <- 5

starts <- 100

dataset <- "cropseq"

#dataset <- "perturbseq_cc7d"

#dataset <- "perturbseq_p7d"

bigorsmall <- "small"

lls <- matrix(0, 5, starts)

llmins <- matrix(0, 5, starts)

resMax <- list()

for (i in 1:starts) {
    print(i)
    if (file.exists(paste0("cropseq/", dataset, "_mnem_", bigorsmall, "_run", i, ".rda"))) {
        load(paste0("cropseq/", dataset, "_mnem_", bigorsmall, "_run", i, ".rda"))
    } else {
        next()
    }
    lls[1, i] <- res[[1]]$ll
    for (j in 1:min(maxk, length(res))) {
        lls[j, i] <- res[[j]]$ll
        llmins[j, i] <- min(res[[j]]$limits[[1]]$ll)
        if (i == 1 | length(resMax) < j) {
            resMax[[j]] <- res[[j]]
        } else {
            if (resMax[[j]]$ll < res[[j]]$ll) {
                resMax[[j]] <- res[[j]]
            }
        }
    }
}

## load log odds before next step

res <- resMax

res[[5]]$data <- res[[4]]$data <- res[[3]]$data <- res[[2]]$data <- res[[1]]$data <- lods

## save(res, file = paste0(dataset, "_mnem_", bigorsmall, "_final.rda"))

load(paste0(dataset, "_mnem_", bigorsmall, "_final.rda"))

## i <- 2

## llstmp <- c(lls[i, ], llmins[i, ])

## llstmp <- llstmp[which(llstmp != 0)]

## plot(llstmp)

## get bics:

maxk <- length(res)

source("mnem/R/mnems.r")
source("mnem/R/mnems_low.r")
library(graph)
library(nem)
library(Rgraphviz)
library(epiNEM)
library(naturalsort)

bics <- rep(0, maxk)

ll <- rep(0, maxk)

for (i in 1:maxk) {
    
    bics[i] <- getIC(res[[i]], useF = F, Fnorm = T)

    ll[i] <- max(res[[i]]$ll)

}

ll2 <- ll

ll <- (ll/(max(ll)-min(ll)))*(max(bics)-min(bics))

ll <- ll - min(ll) + min(bics)

ll3 <- seq(min(bics), max(bics[!is.infinite(bics)]), length.out = 5)

if (length(grep("perturbseq", dataset)) != 0) {
    bigorsmall <- "" # for perturbseq
} else {
    bigorsmall <- paste0("_", bigorsmall)
}

## pdf(paste0(gsub("perturbseq_", "", dataset), bigorsmall, "_bic.pdf"), width = 7, height = 6)

par(mar=c(5,5,2,5))
plot(bics, type = "b", ylab = "modified log likelihood ratio", col = "red", xlab = "modified (red) and raw log likelihood (blue) as a function of number of components", yaxt = "n", ylim = c(min(min(bics,ll)), max(max(bics,ll))), xaxt = "n")
lines(ll, type = "b", col = "blue")
axis(4, ll3, round(seq(min(ll2), max(ll2), length.out = 5)))
axis(2, ll3, round(ll3))
axis(1, 1:maxk, 1:maxk)
mtext("raw log likelihood ratio", side=4, line=3)

dev.off()

i <- which.min(bics)

gamma <- getAffinity(res[[i]]$probs, mw = res[[i]]$mw)

## HeatmapOP(gamma, breaks = 100, bordercol = "transparent")

pdf(paste0(gsub("perturbseq_", "", dataset), bigorsmall, "_hist.pdf"), width = 5, height = 5)

hist(gamma, main = "Histogram of responsibilities", xlab = "responsibilities")

dev.off()

i <- which.min(bics)

pdf(paste0(gsub("perturbseq_", "", dataset), bigorsmall, "_gamma.pdf"), width = 30, height = 10)

gamma <- getAffinity(res[[i]]$probs, mw = res[[i]]$mw)

HeatmapOP(gamma, breaks = 100, bordercol="transparent", cexCol = 0.4, xrot = 45, Rowv = F)

dev.off()

pdf(paste0(gsub("perturbseq_", "", dataset), bigorsmall, "_network.pdf"), width = 30, height = 10)

plot(res[[i]])

dev.off()

pdf(paste0(gsub("perturbseq_", "", dataset), bigorsmall, "_network_overfit.pdf"), width = 30, height = 10)

plot(res[[(i+1)]])

dev.off()

if (length(grep("cc7d", dataset)) == 0) {
    pdf(paste0(gsub("perturbseq_", "", dataset), bigorsmall, "_network_wo.pdf"), width = 10, height = 7)
} else {
    pdf(paste0(gsub("perturbseq_", "", dataset), bigorsmall, "_network_wo.pdf"), width = 14, height = 14) # cc7d
}

plot(res[[i]], egenes = F, bestCell=F, cells = F, showweights = F)

dev.off()

if (length(grep("p7d", dataset)) == 1) {
    pdf(paste0(gsub("perturbseq_", "", dataset), bigorsmall, "_network_overfit_wo.pdf"), width = 8, height = 5)
    plot(res[[(i+1)]], egenes = F, bestCell=F, cells = F, showweights = F)
    dev.off()
}

if (length(grep("cc7d", dataset)) == 1) {
    pdf(paste0(gsub("perturbseq_", "", dataset), bigorsmall, "_network_overfit_wo.pdf"), width = 11.5, height = 8.5)
    plot(res[[(i+1)]], egenes = F, bestCell=F, cells = F, showweights = F)
    dev.off()
}

pdf(paste0(gsub("perturbseq_", "", dataset), bigorsmall, "_kegg.pdf"), width = 10, height = 10)

plot.adj(kadjagg)

dev.off()

## get distance for kegg and combined result:

library(bnstruct)
#kadjagg2 <- transitive.closure(kadjagg[order(rownames(kadjagg)), order(colnames(kadjagg))], mat = T)
kadjagg2 <- kadjagg[order(rownames(kadjagg)), order(colnames(kadjagg))]
diag(kadjagg2) <- 0
#resnem <- transitive.closure(annotAdj(res[[1]]$comp[[1]]$phi, lods), mat = TRUE)
resnem <- annotAdj(res[[1]]$comp[[1]]$phi, lods)
diag(resnem) <- 0
shd(resnem, kadjagg2)

resagg <- annotAdj(res[[i]]$comp[[1]]$phi, lods)
for (j in 1:i) {
    #resagg <- resagg + transitive.closure(res[[i]]$comp[[j]]$phi, mat = TRUE)
    resagg <- resagg + res[[i]]$comp[[j]]$phi
    print(shd(kadjagg2, res[[i]]$comp[[j]]$phi))
}
diag(resagg) <- 0
resagg[which(resagg > 0)] <- 1

shd(resagg, kadjagg2)

par(mfrow=c(1,3))
plot.adj(kadjagg2)
plot.adj(resagg)
plot.adj(resnem)

annotAdj <- function(adj, data) {
    Sgenes <- sort(unique(colnames(data)))
    colnames(adj) <- rownames(adj) <- sort(Sgenes)
    return(adj)
}

## ```

