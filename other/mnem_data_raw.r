### job script to turn normalized count data in log likelihood ratios and apply mnem.

##library(mnem)
library(cluster)
library(nem)
library(Rgraphviz)
library(naturalsort)
library(snowfall)
library(SCnorm)
library(Linnorm)
source("mnem/R/mnems.r")
source("mnem/R/mnems_low.r")

## define dataset

#dataset <- "cropseq"

#dataset <- "perturbseq_cc7d"

#dataset <- "perturbseq_p7d"

args <- commandArgs()

dataset <- gsub("dataset=", "", args[grep("dataset=", args)])

parallel <- gsub("cores=", "", args[grep("cores=", args)])

starts <- gsub("starts=", "", args[grep("starts=", args)])

run <- gsub("run=", "", args[grep("run=", args)])

dollr <- as.numeric(gsub("dollr=", "", args[grep("dollr=", args)]))

donorm <- as.numeric(gsub("donorm=", "", args[grep("donorm=", args)]))

dosmall <- as.numeric(gsub("dosmall=", "", args[grep("dosmall=", args)]))

addendum <- paste0("_run", run)

maxk <- 5

print(dataset)
print(dosmall)
print(dollr)
print(donorm)

####################################### perturb-seq and crop-seq log likelihood ratio:

if (donorm) {

    load(paste0(dataset, "_data.rda"))

    if (length(grep("perturbseq", dataset)) == 0) {

        exprslvl <- apply(data, 1, median)
        
        data <- data[which(exprslvl > 0), ]

        data <- t(t(data)/(colSums(data)/10000))

        data <- Linnorm(data)

        ## data <- log2(data + 0.5)
    
    } else {

        data <- exp(data) - 1

        exprslvl <- apply(data, 1, median)

        data <- data[which(exprslvl > 0), ]

        data <- Linnorm(data)

        colnames(data)[grep("INTER", colnames(data))] <- ""

        ## data <- log2(data + 0.5)

    }
    
    save(data, file = paste0(dataset, "_data_norm.rda"))

    print("norm done")

}

if (dollr) {

    load(paste0(dataset, "_data_norm.rda"))

    llr <- data*0

    C <- which(colnames(data) %in% "")

    distrPar <- function(i, data, C) {
        llrcol <- numeric(ncol(data))
        for (j in which(!(colnames(data) %in% ""))) {
            gene <- colnames(data)[j]
            D <- which(colnames(data) %in% gene)
            cdistr <- ecdf(data[i, C])
            ddistr <- ecdf(data[i, D])
            llrcol[j] <- log2(min(ddistr(data[i, j]), 1 - ddistr(data[i, j]))/min(cdistr(data[i,j]), 1 - cdistr(data[i,j])))
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

    ## if (length(grep("perturbseq", dataset)) == 0) {

    ##     llr <- llr[, -grep("DHODH|MVD|TUBB", colnames(llr))]


    ## }

    llr <- t(apply(llr, 1, function(x) {
        x[is.infinite(x)] <- max(x[!is.infinite(x)])
        return(x)
    }))

    colnames(llr) <- toupper(colnames(llr))

    if (length(grep("perturbseq", dataset)) == 0) {

        search <- "greedy"

        cropgenes <- c("LCK", "ZAP70", "PTPN6", "DOK2", "PTPN11", "EGR3", "LAT")
        
        lods <- llr[, which(colnames(llr) %in% cropgenes)]
        
    } else {

        search <- "greedy"
        
        lods <- llr
        
    }

    badgenes <- "Tcrlibrary"
    badgenes <- grep(badgenes, rownames(lods))

    if (length(badgenes) > 0) {
        lods <- lods[-badgenes, ]
    }

    sdev <- apply(lods, 1, sd)

    lods <- lods[which(sdev > sd(lods)), ]

    print(dim(lods))

    n <- length(unique(colnames(lods)))

    lods <- lods

    bics <- rep(Inf, maxk)

    res <- list()
    
    for (k in 1:maxk) {

        res[[k]] <- mnem(lods, starts = starts, parallel = parallel, k = k, verbose = TRUE, converged = 10^-1, search = search)
        
        res[[k]]$data <- NULL # important to save space
        
        save(res, file = paste0("cropseq/", dataset, "_mnem_small", addendum, ".rda"))
        
    }
    
    save(res, file = paste0("cropseq/", dataset, "_mnem_small", addendum, ".rda"))

    stop("small set done")

}

###########################################

stop()

#### analyze:

maxk <- 5

starts <- 100

dataset <- "cropseq"

#dataset <- "perturbseq_cc7d"

#dataset <- "perturbseq_p7d"

## load data:

load(paste0(dataset, "_llr.rda"))

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

badgenes <- "Tcrlibrary"
badgenes <- grep(badgenes, rownames(lods))

if (length(badgenes) > 0) {
    lods <- lods[-badgenes, ]
}

sdev <- apply(lods, 1, sd)

lods <- lods[which(sdev > sd(lods)), ]

print(dim(lods))

n <- length(unique(colnames(lods)))

lods <- lods

## get results:

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

cropres <- res

cc7dres <- res

p7dres <- res

res2 <- list(cropres, cc7dres, p7dres)

app <- res2

## make app smaller:

for (i in 1:3) {
    for (j in 1:5) {
        app[[i]][[j]]$data <- app[[i]][[j]]$data[1, , drop = FALSE]
        app[[i]][[j]]$limits <- NULL
    }
}

## save(app, file = "app.rda")

## some analysis:

load("app.rda")

res <- app[[2]]

ics <- lls <- numeric(maxk)

for (i in 1:maxk) {
    ics[i] <- getIC(res[[i]])
    lls[i] <- res[[i]]$ll
    print(getIC(res[[i]]))
}

par(mfrow=c(1,2))
plot(ics, type = "b")
plot(lls, type = "b")

paltmp <- palette()

paltmp[3] <- "blue"

paltmp[4] <- "brown"

palette(paltmp)

pdf("temp.pdf", width = 12, height = 10)#width = 16, height = 14)#width = 10, height = 8)
plot(res[[which.min(ics)]], showweights = FALSE, showprobs = FALSE, bestCell = FALSE, shownull = FALSE, cells = FALSE, egenes = FALSE)
dev.off()

sets <- c("cropseq_small", "cc7d", "p7d")

for (i in 1:3) {
    tmplls <- numeric()
    for (j in 1:5) {
        tmplls <- c(tmplls, getIC(app[[i]][[j]]))
    }
    for (j in 1:5) {
        if (j == 3) { tmp <- "_overfit" }
        if (j == 1) { tmp <- "_underfit" }
        if (j == 2) { tmp <- "" }
        if (j %in% c(4,5)) { tmp <- j }
        tmpwidth <- j*8
        if (j == 1 & i == 2) { tmpwidth <- 2*tmpwidth }
        tmpwidth <- tmpwidth + 4
        pdf(paste0(sets[i], "_network", tmp, ".pdf"), height = 14, width = tmpwidth)
        plot(app[[i]][[j]])
        dev.off()

    }
}

##

graph <- c("TOP=PTGER2", "TOP=CIT", "PTGER2=AURKA", "CIT=AURKA")

pdf("Fig9.pdf", width = 8, height = 4)
par(mfrow=c(1,2))
plotDnf(graph, width = 1, nodelabel = list(TOP = ""),
        bordercol = list(TOP = "transparent", PTGER2 = "blue", CIT = "blue", AURKA = "blue"))
plotDnf(graph, width = 1, nodelabel = list(TOP = "", AURKA = "Hidden\n Player"),
        bordercol = list(TOP = "transparent", PTGER2 = "blue", CIT = "blue", AURKA = "blue"))
dev.off()



transitions_0 <- matrix(c(4,1,0), nrow=3, ncol=1)
transitions_1 <- matrix(c(3,1,0,0,1,0,0,0,0), nrow = 3, ncol = 3)
transitions_2 <- matrix(c(2,0,1,0,0,2,2,1,3), nrow =3, ncol =3)
transitions_3 <- matrix(c(4,1,0), nrow = 1, ncol = 3)
transitions <- list(transitions_0, transitions_1, transitions_2, transitions_3)

trans_prob_0 <- matrix(rep(0,6), nrow = 3, ncol = 3)
trans_prob_1 <- matrix(rep(0,9), nrow = 3, ncol = 3)
trans_prob_2 <- matrix(rep(0,9), nrow = 3, ncol = 3)
trans_prob_3 <- matrix(rep(0,3), nrow = 3, ncol = 3)
transition_probabilities <- list(trans_prob_0, trans_prob_1, trans_prob_2, trans_prob_3)

for (i in 1:4){
    trans <- transitions[[i]]
    for (k in 1:ncol(trans)){
        for (l in 1:nrow(trans)){
            transition_probabilities[[i]][l,k] <- (trans[l,k] + 1)/(sum(trans[,k])+nrow(trans)) }
    }
    print( transition_probabilities[[i]])
}

