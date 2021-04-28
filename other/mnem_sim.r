
system("rm .RData")

rm(list=ls())

## i <- as.numeric(commandArgs(TRUE)[1])

## source("mnems.r")
## source("mnems_low.r")
## library(data.table)
## library(naturalsort)
## library(nem)
## library(Linnorm)

## app <- createApp(sets = i, dataonly = TRUE, multi = TRUE)

## save(app, file = paste0("app", i, ".rda"))

## stop()

##

load("sim_edgeprobs.rda")

doSgenes <- c(5,10,25,50,100,200)

source("mnems.r")
source("mnems_low.r")
library(Rcpp)
sourceCpp(code=readChar("mm.cpp", file.info("mm.cpp")$size))
source("nempi_main.r")
source("nempi_low.r")

library(nem)
library(naturalsort)
library(Rgraphviz)
library(cluster)
library(matrixStats)
library(e1071)

Sgenes <- as.numeric(commandArgs(TRUE)[1])
Egenes <- as.numeric(commandArgs(TRUE)[2])
uninform <- as.numeric(commandArgs(TRUE)[3])
starts <- as.numeric(commandArgs(TRUE)[4])
dorun <- commandArgs(TRUE)[5]
if (length(grep(":", dorun)) > 0) {
    dorun <- as.numeric(gsub(":.*", "", dorun)):as.numeric(gsub(".*:", "", dorun))
} else {
    dorun <- as.numeric(dorun)
}
donoise <- commandArgs(TRUE)[6]
if (length(grep(":", donoise)) > 0) {
    donoise <- as.numeric(gsub(":.*", "", donoise)):as.numeric(gsub(".*:", "", donoise))
} else {
    donoise <- as.numeric(donoise)
}
doNem <- commandArgs(TRUE)[7]
if (length(grep(":", doNem)) > 0) {
    doNem <- as.numeric(gsub(":.*", "", doNem)):as.numeric(gsub(".*:", "", doNem))
} else if (length(grep(",", doNem)) > 0) {
    doNem <- as.numeric(unlist(strsplit(doNEM, ",")))
} else{
    doNem <- as.numeric(doNem)
}
nCells <- as.numeric(commandArgs(TRUE)[8])
scalefree <- as.numeric(commandArgs(TRUE)[9])
prob <- as.numeric(commandArgs(TRUE)[10])
subn <- as.numeric(commandArgs(TRUE)[11])
subruns <- as.numeric(commandArgs(TRUE)[12])
donempi <- as.numeric(commandArgs(TRUE)[13])
subtype <- commandArgs(TRUE)[14]
if (is.na(subtype)) { subtype <- "networks" }

save <- TRUE

## Sgenes <- 10; Egenes <- 2; uninform <- floor(Sgenes*Egenes*0.1); starts <- 10; dorun <- 1; donoise <- 1; doNem <- 2; nCells <- 1000; scalefree <- 0; prob <- 0.2; save <- FALSE; subn <- 2; subruns <- 100; donempi <- 1; parallel <- 4; subtype <- "networks"

print(Sgenes)
print(Egenes)
print(uninform)
print(starts)
print(dorun)
print(donoise)
print(doNem)

runs <- dorun
Nems <- doNems <- 1:100
noises <- 1:3

if (prob == 0) {
    bestprob <- min(edgeprobs[which.min(abs(doSgenes - Sgenes))])
} else {
    bestprob <- prob
}

results <- array(0, c(length(runs), length(noises), length(Nems), 9, 9), dimnames = list(paste0("run", seq_len(length(runs))), paste0("noise", noises), paste0("k=", Nems), c("nem", "networks", "networks2", "cluster", "cluster2", "cluster3", "random", "kmeans", "mergeN"), c("adjusted rand index", "rand index", "component ratio", "hamming accuracy", "pairwise max accuracy strict", "strict sens", "strict spec", "time", "ll")))

gtnfeatures <- array(0, c(length(runs), length(noises), length(Nems), 5), dimnames = list(paste0("run", seq_len(length(runs))), paste0("noise", noises), paste0("k=", Nems), c("network penalty", "median mw", "network probability", "min/max mw", "mw sd")))

## if (Sgenes <= 50) {
##     type <- "random"
##     search <- "greedy"
## } else {
##     type <- "random"
##     search <- "estimate"
## }

search <- "greedy"
parallel <- 0

for (run in dorun) {
    cat(paste0(run, "."))
    for (j in donoise) {
        for (i in doNem) {
            ## run <- 1; j <- 1; i <- 2
            mw <- runif(i, 0.1, 1)
            mw <- mw/sum(mw)
            data <- simData(Sgenes=Sgenes, Nems=i, nCells=nCells, uninform=uninform, Egenes=Egenes, edgeprob = bestprob, mw = mw, scalefree = scalefree)

            lods <- data$data
            lods[which(data$data == 1)] <- rnorm(sum(data$data == 1), 1, j)
            lods[which(data$data == 0)] <- rnorm(sum(data$data == 0), -1, j)

            cordist <- as.dist((1-cor(lods))/2)

            simfull <- bigphi(data)

            gtnfeatures[run, j, i, 1] <- calcEvopen(sortAdj(data$Nem, list = TRUE)$res, list = TRUE)
            gtnfeatures[run, j, i, 2] <- median(mw)
            gtnfeatures[run, j, i, 3] <- fullTransProb(data$Nem)
            gtnfeatures[run, j, i, 4] <- min(mw)/max(mw)
            gtnfeatures[run, j, i, 5] <- sd(mw)

            if (Sgenes <= 20) {
                start <- as.numeric(format(Sys.time(), "%s"))
                mnemres <- mnem(lods, search = search, k = 1)
                results[run, j, i, 1, 8] <- as.numeric(format(Sys.time(), "%s")) - start
                resfull <- bigphi(mnemres)
                results[run, j, i, 1, 4] <- hamSim(simfull, resfull)
                results[run, j, i, 1, 5] <- fitacc(mnemres, data, strict = TRUE)
                results[run, j, i, 1, 6] <- fitacc(mnemres, data, strict = TRUE, type = "sens")
                results[run, j, i, 1, 7] <- fitacc(mnemres, data, strict = TRUE, type = "spec")
                cells <- rep(1, length(data$index))
                results[run, j, i, 1, 1] <- classAgreement(table(data$index, cells))$crand
                results[run, j, i, 1, 2] <- classAgreement(table(data$index, cells))$rand
                results[run, j, i, 1, 3] <- nrow(mnemres$probs)/i
                results[run, j, i, 1, 9] <- mnemres$ll

                start <- as.numeric(format(Sys.time(), "%s"))
                mnemres <- mnem(lods, search = search, starts = starts, type = "networks")
                results[run, j, i, 2, 8] <- as.numeric(format(Sys.time(), "%s")) - start
                resfull <- bigphi(mnemres)
                results[run, j, i, 2, 4] <- hamSim(simfull, resfull)
                results[run, j, i, 2, 5] <- fitacc(mnemres, data, strict = TRUE)
                results[run, j, i, 2, 6] <- fitacc(mnemres, data, strict = TRUE, type = "sens")
                results[run, j, i, 2, 7] <- fitacc(mnemres, data, strict = TRUE, type = "spec")
                cells <- apply(getAffinity(mnemres$probs, mw = mnemres$mw), 2, which.max)
                results[run, j, i, 2, 1] <- classAgreement(table(data$index, cells))$crand
                results[run, j, i, 2, 2] <- classAgreement(table(data$index, cells))$rand
                results[run, j, i, 2, 3] <- nrow(mnemres$probs)/i
                results[run, j, i, 2, 9] <- mnemres$ll

                start <- as.numeric(format(Sys.time(), "%s"))
                mnemres <- mnem(lods, search = search, starts = starts, type = "networks2")
                results[run, j, i, 3, 8] <- as.numeric(format(Sys.time(), "%s")) - start
                resfull <- bigphi(mnemres)
                results[run, j, i, 3, 4] <- hamSim(simfull, resfull)
                results[run, j, i, 3, 5] <- fitacc(mnemres, data, strict = TRUE)
                results[run, j, i, 3, 6] <- fitacc(mnemres, data, strict = TRUE, type = "sens")
                results[run, j, i, 3, 7] <- fitacc(mnemres, data, strict = TRUE, type = "spec")
                cells <- apply(getAffinity(mnemres$probs, mw = mnemres$mw), 2, which.max)
                results[run, j, i, 3, 1] <- classAgreement(table(data$index, cells))$crand
                results[run, j, i, 3, 2] <- classAgreement(table(data$index, cells))$rand
                results[run, j, i, 3, 3] <- nrow(mnemres$probs)/i
                results[run, j, i, 3, 9] <- mnemres$ll

                start <- as.numeric(format(Sys.time(), "%s"))
                mnemres <- mnem(lods, search = search, starts = starts, type = "cluster")
                results[run, j, i, 4, 8] <- as.numeric(format(Sys.time(), "%s")) - start
                resfull <- bigphi(mnemres)
                results[run, j, i, 4, 4] <- hamSim(simfull, resfull)
                results[run, j, i, 4, 5] <- fitacc(mnemres, data, strict = TRUE)
                results[run, j, i, 4, 6] <- fitacc(mnemres, data, strict = TRUE, type = "sens")
                results[run, j, i, 4, 7] <- fitacc(mnemres, data, strict = TRUE, type = "spec")
                cells <- apply(getAffinity(mnemres$probs, mw = mnemres$mw), 2, which.max)
                results[run, j, i, 4, 1] <- classAgreement(table(data$index, cells))$crand
                results[run, j, i, 4, 2] <- classAgreement(table(data$index, cells))$rand
                results[run, j, i, 4, 3] <- nrow(mnemres$probs)/i
                results[run, j, i, 4, 9] <- mnemres$ll

                start <- as.numeric(format(Sys.time(), "%s"))
                mnemres <- mnem(lods, search = search, starts = starts, type= "cluster2")
                results[run, j, i, 5, 8] <- as.numeric(format(Sys.time(), "%s")) - start
                resfull <- bigphi(mnemres)
                results[run, j, i, 5, 4] <- hamSim(simfull, resfull)
                results[run, j, i, 5, 5] <- fitacc(mnemres, data, strict = TRUE)
                results[run, j, i, 5, 6] <- fitacc(mnemres, data, strict = TRUE, type = "sens")
                results[run, j, i, 5, 7] <- fitacc(mnemres, data, strict = TRUE, type = "spec")
                cells <- apply(getAffinity(mnemres$probs, mw = mnemres$mw), 2, which.max)
                results[run, j, i, 5, 1] <- classAgreement(table(data$index, cells))$crand
                results[run, j, i, 5, 2] <- classAgreement(table(data$index, cells))$rand
                results[run, j, i, 5, 3] <- nrow(mnemres$probs)/i
                results[run, j, i, 5, 9] <- mnemres$ll

                start <- as.numeric(format(Sys.time(), "%s"))
                mnemres <- mnem(lods, search = search, starts = starts, type= "cluster3")
                results[run, j, i, 6, 8] <- as.numeric(format(Sys.time(), "%s")) - start
                resfull <- bigphi(mnemres)
                results[run, j, i, 6, 4] <- hamSim(simfull, resfull)
                results[run, j, i, 6, 5] <- fitacc(mnemres, data, strict = TRUE)
                results[run, j, i, 6, 6] <- fitacc(mnemres, data, strict = TRUE, type = "sens")
                results[run, j, i, 6, 7] <- fitacc(mnemres, data, strict = TRUE, type = "spec")
                cells <- apply(getAffinity(mnemres$probs, mw = mnemres$mw), 2, which.max)
                results[run, j, i, 6, 1] <- classAgreement(table(data$index, cells))$crand
                results[run, j, i, 6, 2] <- classAgreement(table(data$index, cells))$rand
                results[run, j, i, 6, 3] <- nrow(mnemres$probs)/i
                results[run, j, i, 6, 9] <- mnemres$ll

                as.numeric(format(Sys.time(), "%s"))
                mnemres <- mnem(lods, search = search, starts = starts, type= "random")
                results[run, j, i, 7, 8] <- as.numeric(format(Sys.time(), "%s")) - start
                resfull <- bigphi(mnemres)
                results[run, j, i, 7, 4] <- hamSim(simfull, resfull)
                results[run, j, i, 7, 5] <- fitacc(mnemres, data, strict = TRUE)
                results[run, j, i, 7, 6] <- fitacc(mnemres, data, strict = TRUE, type = "sens")
                results[run, j, i, 7, 7] <- fitacc(mnemres, data, strict = TRUE, type = "spec")
                cells <- apply(getAffinity(mnemres$probs, mw = mnemres$mw), 2, which.max)
                results[run, j, i, 7, 1] <- classAgreement(table(data$index, cells))$crand
                results[run, j, i, 7, 2] <- classAgreement(table(data$index, cells))$rand
                results[run, j, i, 7, 3] <- nrow(mnemres$probs)/i
                results[run, j, i, 7, 9] <- mnemres$ll
            }

            if (Sgenes <= 20) { clustnem <- TRUE } else { clustnem <- FALSE }
            start <- as.numeric(format(Sys.time(), "%s"))
            cres <- clustNEM(lods, search = search, starts = starts, nem = clustnem)
            results[run, j, i, 8, 8] <- as.numeric(format(Sys.time(), "%s")) - start
            if (Sgenes <= 20) {
                resfull <- bigphi(cres)
                results[run, j, i, 8, 4] <- hamSim(simfull, resfull)
                results[run, j, i, 8, 5] <- fitacc(cres, data, strict = TRUE)
                results[run, j, i, 8, 6] <- fitacc(cres, data, strict = TRUE, type = "sens")
                results[run, j, i, 8, 7] <- fitacc(cres, data, strict = TRUE, type = "spec")
            }
            results[run, j, i, 8, 1] <- classAgreement(table(data$index, cres$cluster))$crand
            results[run, j, i, 8, 2] <- classAgreement(table(data$index, cres$cluster))$rand
            if (Sgenes <= 20) {
                results[run, j, i, 8, 3] <- nrow(cres$probs)/i
                results[run, j, i, 8, 9] <- cres$ll
            } else {
                results[run, j, i, 8, 3] <- length(table(cres$cluster))/i
            }

            start <- as.numeric(format(Sys.time(), "%s"))
            subcl <- matrix(0, subruns, ncol(lods))
            Sgenes2 <- unique(unlist(strsplit(colnames(lods), "_")))
            for (subrun in 1:subruns) {
                sub <- sample(Sgenes2, subn)
                subs <- grep(paste(paste0("^", sub, "$"), collapse = "|"), colnames(lods))
                lodsub <- lods
                colnames(lodsub)[-subs] <- ""
                ## this would be super advanced:
                if (donempi) {
                    np <- nempi(lodsub, converged = 0.1, full = 0)
                    Rho <- np$Gamma
                    if (nrow(Rho) == subn+1) { stop() }
                    Rho[which(Rho > 0.5)] <- 1
                    Rho[which(Rho < 1)] <- 0
                    clidx <- which(apply(Rho, 2, sum) > 0)
                    lodsub <- lodsub[, clidx]
                    Rho <- Rho[, clidx]
                    msub <- mnem(lodsub, search = search, starts = starts, type = subtype, Rho = Rho)
                } else {
                    clidx <- which(colnames(lodsub) != "")
                    lodsub <- lodsub[, clidx]
                    msub <- mnem(lodsub, search = search, starts = starts, type = subtype)
                }
                cl <- getAffinity(msub$probs, mw = msub$mw, complete = TRUE)
                cl <- apply(cl, 2, function(x) {
                    y <- which.max(x);
                    return(y)
                })
                if (length(cl) == 0) { cl <- rep(1, ncol(lodsub)) }
                subcl[subrun, clidx] <- cl # somehow add the uncertainty? clnumber.uncertainty...
                cat(paste0(subrun, ","))
            }
            cooc <- matrix(0, ncol(lods), ncol(lods))
            coocfun <- function(x, subcl) {
                acl <- x
                bcl <- abs(subcl - acl)
                bcl[which(acl == 0), ] <- 1
                bcl[which(abs(bcl) >= 1)] <- 1
                return(matrix(1, 1, nrow(bcl))%*%(1-bcl))
            }
            if (parallel > 0) {
                sfInit(parallel=1, cpus=parallel)
                cooc <- sfApply(subcl, 2, coocfun, subcl)
                sfStop()
            } else {
                cooc <- apply(subcl, 2, coocfun, subcl)
            }
            d <- 1 - cooc/max(cooc)
            diag(d) <- 0
            d <- as.dist(d)
            K <- median(apply(subcl, 1, max))
            kres <- kmeans(d, centers=K, nstart = starts)
            fullcl <- kres
            cells <- fullcl$cluster
            results[run, j, i, 9, 8] <- as.numeric(format(Sys.time(), "%s")) - start
            results[run, j, i, 9, 1] <- classAgreement(table(data$index, cells))$crand
            results[run, j, i, 9, 2] <- classAgreement(table(data$index, cells))$rand
            if (Sgenes <= 20) {
                mnemres <- mnem(lods, search = search, starts = starts, type = "networks")
            }
            if (Sgenes <= 20) {
                resfull <- bigphi(mnemres)
                results[run, j, i, 9, 4] <- hamSim(simfull, resfull)
                results[run, j, i, 9, 5] <- fitacc(mnemres, data, strict = TRUE)
                results[run, j, i, 9, 6] <- fitacc(mnemres, data, strict = TRUE, type = "sens")
                results[run, j, i, 9, 7] <- fitacc(mnemres, data, strict = TRUE, type = "spec")
            }
            if (Sgenes <= 20) {
                results[run, j, i, 9, 3] <- nrow(mnemres$probs)/i
                results[run, j, i, 9, 9] <- mnemres$ll
            } else {
                results[run, j, i, 9, 3] <- K/i
            }

            ## results[run, j, i, , ]

            if (save) {
                path <- "/cluster/work/bewi/members/mpirkl/mnem_sim_results/"
                for (filen in 1:100) {
                    filename <- paste0(path, Sgenes, "_", Egenes, "_", uninform, "_", i, "_", j, "_", donempi, "_", filen, ".rda")
                    if (!file.exists(filename)) {
                        save(results, gtnfeatures, file = filename)
                        break()
                    }
                }
            }
        }
    }
}

stop("success")

## hpc commands:

module load bioconductor/3.6

module load curl/7.49.1

module load gmp/5.1.3

## general: bsub -M 100000 -q normal.24h -n 1 -e error.txt -o output.txt -R "rusage[mem=100000]" "R --silent --no-save --args '2' < testing/vignettes/mnem_sim.r"

ram=1500;

##

Sgenes=50 # Sgenes
Egenes=2 # Egenes
uninform=$(( $(($Sgenes * $Egenes)) / 10)) # uninformative genes
starts=10 # em starts
##
runs=1 # runs
noise=1 # noises (max 3)
comps=2 # components
cells=1000 # number of cells
sf=0
prob=0.1 # check edgemeans
subn=10
subruns=100
nempi=0
mergetype=networks

rm .RData

rm error.txt

rm output.txt

## Sgenes Egenes uninformative emstarts greedyruns noises comps ncells
bsub -M ${ram} -q normal.4h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${Sgenes}' '${Egenes}' '${uninform}' '${starts}' '${runs}' '${noise}' '${comps}' '${cells}' '${sf}' '${prob}' '${subn}' '${subruns}' '${nempi}' '${mergetype}' < mnem_sim.r"

for z in {2..100}; do
bsub -M ${ram} -q normal.4h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${Sgenes}' '${Egenes}' '${uninform}' '${starts}' '${runs}' '${noise}' '${comps}' '${cells}' '${sf}' '${prob}' '${subn}' '${subruns}' '${nempi}' '${mergetype}' < mnem_sim.r"
done

## Sgenes >= 50

for comps in {3..10}; do
for z in {2..100}; do
bsub -M ${ram} -q normal.4h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${Sgenes}' '${Egenes}' '${uninform}' '${starts}' '${runs}' '${noise}' '${comps}' '${cells}' '${sf}' '${prob}' '${subn}' '${subruns}' '${nempi}' < mnem_sim.r"
done
done

## subsetting to >= 10 Sgenes

for comps in {3..10}; do
for noise in {1..3}; do
for z in {2..100}; do
bsub -M ${ram} -q normal.4h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${Sgenes}' '${Egenes}' '${uninform}' '${starts}' '${runs}' '${noise}' '${comps}' '${cells}' '${sf}' '${prob}' '${subn}' '${subruns}' '${nempi}' < mnem_sim.r"
done
done

## read sim results:

Sgenes <- 25
Egenes <- 2
uninform <- floor((Sgenes*Egenes)/10)
donempi <- 1
#path <- "~/Mount/Leoshare/cluster_sim_res/"
path <- "~/Mount/Eulershare/mnem_sim_results/"

full <- list()
library(abind)
for (Nem in 1:10) {
    full[[Nem]] <- list()
    noisear <- noisefit <- NULL
    for (noise in 1:3) {
        acc <- NULL
        for (filen in 1:100) {
            filename <- paste0(path, Sgenes, "_", Egenes, "_", uninform, "_", Nem, "_", noise, "_", donempi, "_", filen, ".rda")
            if (file.exists(filename)) {
                load(filename)
                acc <- abind(acc, results, along = 1)
                cat(paste0(filen, "."))
            } else {
                if (filen == 0) {
                    break()
                }
            }
        }
        full[[Nem]][[noise]] <- acc
    }
}

## save(full, file = paste("mnem", "clust", Sgenes, Egenes, uninform, ".rda", sep = "_"))

load(paste("mnem", "clust", Sgenes, Egenes, uninform, ".rda", sep = "_"))

methods <- c(1,2,3,4,5,6,7,8,9)
feat <- c(1,3,8)
pdf("temp.pdf", height = 5, width = length(feat)*5)
par(mfrow=c(1,length(feat)))
library(abind)
for (Nem in 1:length(full)) {
    noisear <- noisefit <- NULL
    if (length(full) < Nem) { break() }
    if (length(full[[Nem]]) == 0) { next() }
    for (noise in seq_len(length(full[[Nem]]))) {
        if (length(full[[Nem]]) < noise) { break() }
        if (!is.null(full[[Nem]][[noise]])) {
            minrun <- min(c(dim(noisear)[1], dim(full[[Nem]][[noise]])[1]))
            feats <- list()
            for (f in 1:length(feat)) {
                if (length(feats) < f) {
                    feats[[f]] <- full[[Nem]][[noise]][1:minrun,noise,Nem,methods,feat[f]]
                } else {
                    feats[[f]] <- cbind(feats[[f]][1:minrun,], full[[Nem]][[noise]][1:minrun,noise,Nem,methods,feat[f]])
                }
            }
        }
    }
    for (f in 1:length(feat)) {
        if (dimnames(full[[Nem]][[noise]])[[5]][feat[f]] %in% c("time", "component ratio")) {
            ylim <- NULL
        } else {
            ylim <- c(0,1)
        }
        boxplot(feats[[f]], main = dimnames(full[[Nem]][[noise]])[[5]][feat[f]], ylim = ylim)
        abline(v=c(0.5+length(methods),0.5+length(methods)*2), lty = 2)
    }
}
dev.off()

## plot:

types <- dim(results)[4]
pdf("temp.pdf", height = 22, width = 22)
par(mfrow=c(dim(results)[5],length(Nemsdone)))
for (j in 1:dim(results)[5]) {
    for (i in 1:length(Nemsdone)) {
        ymax <- max(as.vector(cbind(results[,1,i,1:types,j],results[,2,i,1:types,j],results[,3,i,1:types,j])))
        ymin <- min(as.vector(cbind(results[,1,i,1:types,j],results[,2,i,1:types,j],results[,3,i,1:types,j])))
        if (is.na(ymax)) { ymax <- 1 }
        if (is.na(ymin)) { ymin <- 0 }
        if (j %in% c(4,5,7)) { ymin <- 0.5 }
        boxplot(cbind(results[,1,i,1:types,j],results[,2,i,1:types,j],results[,3,i,1:types,j]), col = 2:7, ylim = c(ymin,max(ymax, 1)), main = paste0("K = ", i), xaxt = "n", xlab = expression(sigma), ylab = dimnames(results)[[5]][j])
        axis(1, 1:(types*3), c(rep(1, types), rep(2, types), rep(3, types)))
        abline(h=0,lty=3)
        abline(v=types*1:2+0.5, col = rgb(0,0,0,0.75), lty = 2)
    }
}
dev.off()

## gtn comparison

pdf("temp.pdf", height = 20, width = 6)
par(mfrow=c(7,2))
K <- 3
i <- 2
j <- 1
results1 <- results0[which(gtnfeatures0[,i,K,j] < median(gtnfeatures0[,i,K,j])),,,,]
results2 <- results0[which(gtnfeatures0[,i,K,j] > median(gtnfeatures0[,i,K,j])),,,,]
boxplot(results1[,i,K,1:types,1], ylim = c(0,1), col = 2:6, main = paste0("K = ", K), xaxt = "n", xlab = expression(sigma), ylab = "Adjusted")
abline(h=seq(0,1,length.out = 11), lty = 2, col = rgb(0,0,0,0.25))

## plot:

Nemsdone <- 1:5

index <- c(1,4:5)

types <- dim(ari)[4]
pdf("temp.pdf", height = 6, width = 14)
par(mfrow=c(2,max(Nemsdone)))
for (i in Nemsdone) {
    boxplot(cbind(results[,1,i,index]), col = 2:6, ylim = c(0.5,1), main = paste0("K = ", i), xaxt = "n", xlab = "", ylab = "Mixture Network Accuracy")
    axis(1, 1:3, c("ks", "kmeans", "hier"))
    abline(h=0,lty=3)
    abline(v=types*1:2+0.5, col = rgb(0,0,0,0.75), lty = 2)
}
for (i in Nemsdone) {
    boxplot(cbind(results[,1,i,index]+1), col = 2:6, main = paste0("K = ", i), xaxt = "n", xlab = "", ylab = "Running Time")#, log = "y")
    axis(1, 1:(types*3), c(rep(1, types), rep(2, types), rep(3, types)))
    abline(h=0,lty=3)
    axis(1, 1:3, c("ks", "kmeans", "hier"))
}
dev.off()

## play time:

source("mnem/R/mnems.r")
source("mnem/R/mnems_low.r")
library(nem)
library(naturalsort)
library(snowfall)
library(Rgraphviz)
library(cluster)
library(matrixStats)
library(microbenchmark)
library(e1071)
library(Rcpp)
sourceCpp("mnem/src/mm.cpp")

library(MASS)

matrixnem <- function(D, nsample = NULL, ...) {
    if (is.null(nsample)) {
        nsample <- nrow(D)
    }
    D <- modData(D)
    D <- doMean(D)
    F <- t(D)
    F[which(F > 0)] <- 1
    F[which(F < 0)] <- 0
    n <- getSgeneN(D)
    Sgenes <- getSgenes(D)
    phi <- matrix(0, n, n)
    diag(phi) <- 1
    colnames(phi) <- rownames(phi) <- Sgenes
    phis <- list()
    scores <- numeric(nsample)
    for (i in seq_len(nsample)) {
        theta <- t(F)*0
        theta <- apply(F, 2, function(x) {
            y <- x*0
            y[sample(which(x == 1), 1)] <- 1
            return(y)
        })
        thetainv <- ginv(theta)
        phi <- F%*%thetainv
        phis[[i]] <- mytc(phi)
        scores[i] <- sum(rowMaxs(llrScore(D, phis[[i]], ...)))
    }
    phi <- phis[[which.max(scores)]]
    score <- max(scores)
    return(return(list(adj = phi, score = score, scores = scores)))
}

## find perfect probs for a 50% dense network:

freqs <- seq(0.05,0.25,length.out=5)
ngenes <- c(3,5,10,20,30,40,50,100,200)
edgemeds <- matrix(0, length(ngenes), length(freqs))
rownames(edgemeds) <- ngenes
colnames(edgemeds) <- freqs
edgemeans <- edgemeds
size <- 100

for (i in seq_len(length(ngenes))) {
    cat(i)
    n <- ngenes[i]
    for (j in seq_len(length(freqs))) {
        dens <- numeric(size)
        for (k in seq_len(size)) {
            sim <- simData(Sgenes=ngenes[i], Nems=1, nCells=ngenes[i], uninform=0, Egenes=1, edgeprob = freqs[j], scalefree = 0)
            tmp <- mytc(sim$Nem[[1]])
            diag(tmp) <- 0
            dens[k] <- sum(tmp == 1)/(n*(n-1)*0.5)
        }
        edgemeds[i, j] <- median(dens)
        edgemeans[i, j] <- mean(dens)
    }
}

save(edgemeans, edgemeds, file = "sim_edge_densitiy.rda")

## bench mynem:

load("sim_edgeprobs.rda")

Sgenes <- 3
Egenes <- 2
uninform <- 0#floor(Sgenes*Egenes*0.1)

doSgenes <- c(5,10,25,50,100,200)

bestprob <- min(edgeprobs[which.min(abs(doSgenes - Sgenes))+1])

sim <- simData(Sgenes=Sgenes, Nems=1, nCells=Sgenes, uninform=uninform, Egenes=Egenes, edgeprob = bestprob*0.1, scalefree = TRUE)

lods <- sim$data
lods[which(sim$data == 1)] <- rnorm(sum(sim$data == 1), 1, 0.1)
lods[which(sim$data == 0)] <- rnorm(sum(sim$data == 0), -1, 0.1)

phi <- mytc(sim$Nem[[1]])
subtopo <- sim$theta[[1]]
theta <- matrix(0, ncol(phi), length(subtopo))
theta[cbind(subtopo, 1:length(subtopo))] <- 1
F <- phi%*%theta

jacobi <- function(F, s = 5) {

    ranks <- numeric(s)

    for (i in 1:s) {

        P0 <- matrix(runif(ncol(F)*nrow(F)), ncol(F))

        P1 <- matrix(runif(ncol(F)*nrow(F)), ncol(F))

        J <- list()

        count <- 1

        x <- 0

        ## count4 <- 0

        ## count3 <- 1

        ## J[[count3]] <- matrix(0, length(F)*2, length(F)*2)

        J <- matrix(0, 2^ncol(F)-1, 2*ncol(F)*nrow(F))

        X <- matrix(0, 2^ncol(F)-1, ncol(F))

        while(any(x == 0)) {

            ## if (count4 == length(F)*2) {

            ##     J[[count3]] <- diag(eigen(J[[count3]])$values)[1:Matrix::rankMatrix(J[[count3]]), ]

            ##     count3 <- count3 + 1
            ##     J[[count3]] <- matrix(0, length(F)*2, length(F)*2)
            ##     count4 <- 0

            ## }

            x <- numeric(ncol(F))
            dec <- floor(count/2)
            x[1] <- count - dec*2
            count2 <- 2
            while(dec != 0) {
                tmp <- floor(dec/2)
                x[count2] <- dec - tmp*2
                dec <- tmp
                count2 <- count2 + 1
            }
            x <- c(x, rep(0, ncol(F)-length(x)))
            X[count, ] <- x

            P0X <- P0

            P0X[which(x == 0), ] <- 1 - P0X[which(x == 0), ]

            P1X <- P1

            P1X[which(x == 0), ] <- 1 - P1X[which(x == 0), ]

            LL <- log(P1X/P0X)%*%F

            L <- exp(sum(diag(LL)))

            z <- rep(x, nrow(F)*2)

            z <- (z - 0.5)/0.5

            DL <- L/c(as.vector(P1X), as.vector(P0X))*z

            DL <- DL*c(rep(1, length(F)), rep(-1, length(F)))

            DL <- DL*c(as.vector(t(F)), as.vector(t(F)))

            ## J[[count3]][count, ] <- DL

            J[count, ] <- DL

            ## count4 <- count4 + 1

            count <- count + 1

        }

        ranks[i] <- Matrix::rankMatrix(J)

    }

}

search <- "greedy" # search <- "estimate" # search <- "fast" # search <- "greedy.rank"

## source("mnem/R/mnems_low.r"); sourceCpp("mnem/src/mm.cpp")

Rprof("temp.txt", line.profiling=TRUE)
nem2 <- mynem(lods, search = search)#, start = sim$Nem[[1]])
Rprof(NULL)
summaryRprof("temp.txt", lines = "show")$sampling.time
head(summaryRprof("temp.txt", lines = "show")$by.self, 10)

sum(abs(transitive.closure(sim$Nem[[1]], mat = TRUE) - transitive.closure(nem2$adj, mat = TRUE)) != 0)

Rprof("temp.txt", line.profiling=TRUE)
nem2 <- matrixnem(lods)
Rprof(NULL)
summaryRprof("temp.txt", lines = "show")$sampling.time
head(summaryRprof("temp.txt", lines = "show")$by.self, 10)

sum(abs(transitive.closure(sim$Nem[[1]], mat = TRUE) - transitive.closure(nem2$adj, mat = TRUE)) != 0)

par(mfrow=c(1,2))
plotDnf(sim$Nem[[1]])
plotDnf(nem2$adj)

runs <- 10
ham <- ll <- time <- matrix(0, runs, 3)
for (i in 1:runs) {
    cat(i)
    sim <- simData(Sgenes=Sgenes, Nems=1, nCells=Sgenes, uninform=uninform, Egenes=Egenes, edgeprob = bestprob*0.1, scalefree = TRUE)
    lods <- sim$data
    lods[which(sim$data == 1)] <- rnorm(sum(sim$data == 1), 1, 0.1)
    lods[which(sim$data == 0)] <- rnorm(sum(sim$data == 0), -1, 0.1)

    start <- as.numeric(Sys.time())
    tmp <- mynem(lods, search = "fast")
    time[i, 1] <- as.numeric(Sys.time()) - start
    start <- as.numeric(Sys.time())
    tmp2 <- mynem(lods, search = "estimate")
    time[i, 2] <- as.numeric(Sys.time()) - start
    start <- as.numeric(Sys.time())
    ## control <- set.default.parameters(setdiff(unique(colnames(lods)),"time"))
    ## control$type <- "CONTmLLRatio"
    ## tmp3 <- nem(lods, control=control, inference = search2)
    tmp3 <- matrixnem(lods, nsample = 10)
    time[i, 3] <- as.numeric(Sys.time()) - start
    ## tmp3$adj <- transitive.closure(as(tmp3$graph, "matrix"), mat = TRUE)
    ##tmp3$score <- tmp3$mLL

    ham[i, 1] <- sum(abs(transitive.closure(sim$Nem[[1]], mat = TRUE) - transitive.closure(tmp$adj, mat = TRUE)) != 0)
    ham[i, 2] <- sum(abs(transitive.closure(sim$Nem[[1]], mat = TRUE) - transitive.closure(tmp2$adj, mat = TRUE)) != 0)
    ham[i, 3] <- sum(abs(transitive.closure(sim$Nem[[1]], mat = TRUE) - transitive.closure(tmp3$adj, mat = TRUE)) != 0)

    ll[i, 1] <- tmp$score
    ll[i, 2] <- tmp2$score
    ll[i, 3] <- tmp3$score

}

par(mfrow=c(1,3))
boxplot(ham)
boxplot(ll)
boxplot(time)

##

search2 <- "nem.greedy"
## search2 <- "ModuleNetwork"

Rprof("temp.txt", line.profiling=TRUE)
control <- set.default.parameters(setdiff(unique(colnames(lods)),"time"))
control$type <- "CONTmLLRatio"
nem3 <- nem(lods, control=control, inference = search2)
Rprof(NULL)
summaryRprof("temp.txt", lines = "show")$sampling.time
head(summaryRprof("temp.txt", lines = "show")$by.self, 10)

sum(abs(transitive.closure(sim$Nem[[1]], mat = TRUE) - transitive.closure(as(nem3$graph, "matrix"), mat = TRUE)) != 0)

par(mfrow=c(1,2))
plotDnf(sim$Nem[[1]])
plotDnf(transitive.closure(as(nem3$graph, "matrix"), mat = TRUE))

Sgenes <- 10
microbenchmark(
    transitive.closure(matrix(sample(c(0,1), Sgenes^2, replace = TRUE), Sgenes, Sgenes), mat = TRUE),
    mytc(matrix(sample(c(0,1), Sgenes^2, replace = TRUE), Sgenes, Sgenes)),
    times = 10000)

microbenchmark(
    sort(unique(which(adj1 - oldadj != 0, arr.ind = TRUE)[, 2])),
    unique(which(adj1 - oldadj != 0, arr.ind = TRUE)[, 2]),
    which(adj1 - oldadj != 0, arr.ind = TRUE)[, 2],
    times = 100)

## bench mnem:

source("mnem/R/mnems.r")
source("mnem/R/mnems_low.r")
library(nem)
library(naturalsort)
library(snowfall)
library(Rgraphviz)
library(cluster)
library(matrixStats)
library(microbenchmark)
library(e1071)
library(Rcpp)
sourceCpp("mnem/src/mm.cpp")

load("sim_edgeprobs.rda")

Sgenes <- 20
Egenes <- 10
uninform <- floor(Sgenes*Egenes*0.1)

doSgenes <- c(5,10,25,50,100,200)

bestprob <- 0.25 # min(edgeprobs[which.min(abs(doSgenes - Sgenes))])

mw <- runif(2, 0.1, 1)
mw <- mw/sum(mw)
sim <- simData(Sgenes=Sgenes, nCells=1000, uninform=uninform, Egenes=Egenes, edgeprob = bestprob, mw = mw, scalefree = 0)

lods <- sim$data
lods[which(sim$data == 1)] <- rnorm(sum(sim$data == 1), 1, 0.1)
lods[which(sim$data == 0)] <- rnorm(sum(sim$data == 0), -1, 0.1)

search <- "greedy"

res <- mnem(lods, search = search, phi = list(A = sim$Nem), theta = list(A = sim$theta), mw = mw)

probs <- getAffinity(res$probs, mw = res$mw)
cells <- apply(probs, 2, which.max)

classAgreement(table(sim$index, cells))
fitacc(res, sim, strict = TRUE)
fitacc(res, sim, strict = TRUE, type = "sens")
fitacc(res, sim, strict = TRUE, type = "spec")

par(mfrow=c(2,3))
pcadata <- prcomp(lods)
plot(pcadata$rotation[, 1:2], col = cells+1)
plot(pcadata$rotation[, 1:2], col = sim$index+1)
plotConvergence(res)

##

res <- mnem(lods, starts = 10, search = search, type = "random")

## Rprof("temp.txt", line.profiling=TRUE)
## res <- mnemh(lods, starts = 10, search = search, converged = -Inf, type = "random")
## Rprof(NULL)
## summaryRprof("temp.txt", lines = "show")$sampling.time
## head(summaryRprof("temp.txt", lines = "show")$by.self)

search <- "estimate"

res <- mnem(lods, starts = 10, search = search, type = "cluster3")

res <- mnem(lods, starts = 10, search = search, type = "cluster2")

probs <- getAffinity(res$probs, mw = res$mw)
cells <- apply(probs, 2, which.max)

classAgreement(table(sim$index, cells))
fitacc(res, sim, strict = TRUE)
fitacc(res, sim, strict = TRUE, type = "sens")
fitacc(res, sim, strict = TRUE, type = "spec")

par(mfrow=c(2,3))
pcadata <- prcomp(lods)
plot(pcadata$rotation[, 1:2], col = cells+1)
plot(pcadata$rotation[, 1:2], col = sim$index+1)
plotConvergence(res)

source("mnem/R/mnems.r")
source("mnem/R/mnems_low.r")
library(cluster)

cres <- clustNEM(lods, cluster = NULL, search = search, starts = 10, k = length(res$comp))

classAgreement(table(sim$index, cres$cluster))
fitacc(cres, sim, strict = TRUE)
fitacc(cres, sim, strict = TRUE, type = "sens")
fitacc(cres, sim, strict = TRUE, type = "spec")

par(mfrow=c(1,2))
pcadata <- prcomp(lods)
plot(pcadata$rotation[, 1:2], col = cres$cluster+1)
plot(pcadata$rotation[, 1:2], col = sim$index+1)

##

source("mnem/R/mnems.r")
source("mnem/R/mnems_low.r")
library(Rcpp)
sourceCpp("mnem/src/mm.cpp")

library(nem)
data(Ivanova2006RNAiTimeSeries)

data <- matrix(0, dim(dat)[2], dim(dat)[1]*dim(dat)[3])
for (i in seq_len(dim(dat)[1])) {
    data[, seq_len(dim(dat)[3]) + (i-1)*dim(dat)[3]] <- dat[i,,]
}

rownames(data) <- dimnames(dat)[[2]]
colnames(data) <- rep(dimnames(dat)[[3]], dim(dat)[1])

res <- mnemk(data, ks=1:8, starts = 10, type = "cluster")

par(mfrow=c(1,4))
res <- mnem(data, k = 8, starts = 10)
plotConvergence(res)
abline(h=c(max(res$lls), 4500), lty = 2)

plot(res)

cells <- apply(getAffinity(res$probs, mw = res$best$mw), 2, which.max)

plot(cells, rep(1:8, each = 6), main = cor(cells, rep(1:8, each = 6)))

pcares <- prcomp(data)
plot(pcares$rotation[, 1:2], col = cells)

##

source("mnem/R/mnems.r")
source("mnem/R/mnems_low.r")
library(nem)
library(naturalsort)
library(snowfall)
library(Rgraphviz)
library(cluster)
library(matrixStats)
library(microbenchmark)
library(e1071)
library(kernlab)
library(Rcpp)
sourceCpp("mnem/src/mm.cpp")

Sgenes <- 5
Egenes <- 10

sim <- simData(Sgenes=Sgenes, Nems=1, nCells=100, uninform=0, Egenes=Egenes, edgeprob = 0.25)

Phi <- mytc(sim$Nem[[1]])
Theta <- matrix(0, ncol(Phi), length(sim$theta[[1]]))
Theta[cbind(as.numeric(sim$theta[[1]]), 1:ncol(Theta))] <- 1

test <- function(R, Phi, Theta) {
    F <- Phi%*%Theta
    RF <- R*t(F[colnames(R), ])
    S <- t(RF)%*%RF
    return(S)
}

lods <- sim$data
lods[which(lods == 1)] <- rnorm(sum(lods == 1), 1, 0.5)
lods[which(lods == 0)] <- rnorm(sum(lods == 0), -1, 0.5)

R <- test(lods, Phi, Theta)

class(R) <- "kernelMatrix"

nemca <- kpca(R, features = Sgenes)

dims <- c(1,2)
par(mfrow=c(2,2))
plotDnf(sim$Nem[[1]])
plot(nemca@rotated[, dims], pch = as.character(1:5))

sim2 <- simData(Sgenes=Sgenes, Nems=1, nCells=100, uninform=0, Egenes=Egenes, edgeprob = 0.25)

Phi2 <- mytc(sim2$Nem[[1]])
Theta2 <- matrix(0, ncol(Phi), length(sim2$theta[[1]]))
Theta2[cbind(as.numeric(sim2$theta[[1]]), 1:ncol(Theta2))] <- 1

R2 <- test(lods, Phi2, Theta2)

class(R2) <- "kernelMatrix"

nemca2 <- kpca(R2, features = Sgenes)

plotDnf(sim2$Nem[[1]])
plot(nemca2@rotated[, dims], pch = as.character(1:5))

##

nemsim <- function(x, y, method = "disc") {
    if (method %in% "disc") {
        z <- sum(x == 1 & y == 1)
        u <- sum(x == 1 & y == 0)
        v <- sum(x == 0 & y == 1)
        z <- z - min(u,v)
        return(z)
    }
}

data <- sim$data

flip <- sample(1:length(data), floor(0.1*length(data)))

data[flip] <- 1 - data[flip]

R <- matrix(0, ncol(data), ncol(data))

for (i in 1:ncol(data)) {
    for (j in 1:ncol(data)) {
        R[i,j] <- nemsim(data[, i], data[, j])
    }
}

R <- R - min(R)

dims <- c(1,2)
par(mfrow=c(1,2))
class(R) <- "kernelMatrix"

nemca <- kpca(R, features = Sgenes)

plotDnf(sim$Nem[[1]])
plot(nemca@rotated[, dims], pch = as.character(1:5))

R2 <- t(data)%*%data

class(R2) <- "kernelMatrix"

nemca2 <- kpca(R2, features = Sgenes)

plotDnf(sim$Nem[[1]])
plot(nemca2@rotated[, dims], pch = as.character(1:5))
