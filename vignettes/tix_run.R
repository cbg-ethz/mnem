
## ## test:

## system("rnai-query compose --gene hspa5 --pathogen salmonella --study infectx --design p --sample 1000 /cluster/work/bewi/members/simondi/data/tix/database/tix_index.db temp.tsv")

## pw <- as.numeric(commandArgs(TRUE)[1])

## print(pw)

## stop()

pw <- as.numeric(commandArgs(TRUE)[1])

if (is.na(pw)) { pw <- NULL }

## setwd("~/Mount/Leo/"); pw <- 1

## stop(pw)

## loading packages:

library(nem)
library(graph)
library(cluster)
library(snow)
library(snowfall)

library(KEGG.db)
library(KEGGgraph)
library(Rgraphviz)
library(dplyr)

source("mnem/R/mnems.r")
source("mnem/R/mnems_low.r")
  
## here comes the R code for leo jobs:

tixdata <- function(x) {
    
    tmp <- t(read.delim(x))
    
    meta <- tmp[1:which(rownames(tmp) %in% "object_idx"), ]
    
    tmp <- apply(tmp[-(1:which(rownames(tmp) %in% "object_idx")), ], c(1,2), as.numeric)
    
    colnames(tmp) <- meta[which(rownames(meta) %in% "gene"), ]
    
    plates <- sort(unique(meta[which(rownames(meta) %in% "plate"), ]))
    
    return(list(data = tmp, meta = meta, plates = plates))

}

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


load("tix/kadj.rda")

pathogen <- "salmonella"

load("tix/genes.rda")

cores <- 1

## first get the tsv:

if (is.null(pw)) {
    dopws <- 1:length(kadj)
} else {
    dopws <- pw
}

for (j in dopws) {

    print(names(kadj)[j])

    kadjtmp <- kadj[[j]]

    print(dim(kadjtmp))
    
    kadjtmp <- kadjtmp[which(rownames(kadjtmp) %in% genes[, 1]), which(colnames(kadjtmp) %in% genes[, 1]), drop = FALSE]

    print(dim(kadjtmp))
    
    if (dim(kadjtmp)[1] < 2) { next() }
    
    if (dim(kadjtmp)[1] > 5) {

        sample <- 0
        
        ## how to compress pathway: or sample randomly

        if (sample) {
            randgenes <- sample(1:nrow(kadjtmp), 5)
            kadjtmp <- kadjtmp[randgenes, randgenes]
        } else {
            edgen <- order(apply(kadjtmp, 1, sum)+apply(kadjtmp, 2, sum), decreasing = TRUE)
            
            kadjtmp <- kadjtmp[edgen[1:min(nrow(kadjtmp), 5)], edgen[1:min(nrow(kadjtmp), 5)]]
        }
    }
    
    cells <- 1000
    
    for (i in rownames(kadjtmp)) {

        print(i)

        if (!file.exists(paste0("tix/data/", pathogen, "_", i, ".tsv"))) {
            system(paste("rnai-query compose --gene ",
                         tolower(i),
                         " --pathogen ", pathogen, " --study infectx --design p --sample ", cells, " /cluster/work/bewi/members/simondi/data/tix/database/tix_index.db tix/data/", pathogen, "_", i, ".tsv", sep = ""))
        } #else { next() }

        if (!file.exists(paste0("tix/data/", pathogen, "_", i, ".tsv"))) { next() }

        tmp <- tixdata(paste0("tix/data/", pathogen, "_", i, ".tsv"))

        meta <- tmp$meta

        plates <- tmp$plates

        data <- tmp$data

        for (k in plates) {

            print(k)

            if (!file.exists(paste0("tix/data/", pathogen, "_", k, ".tsv"))) {
                system(paste("rnai-query compose --gene ",
                             "mock",
                             " --plate ",
                             k,
                             " --pathogen ", pathogen, " --study infectx --design p --sample ", cells, " /cluster/work/bewi/members/simondi/data/tix/database/tix_index.db tix/data/", pathogen, "_", k, ".tsv", sep = ""))
            }

            if (!file.exists(paste0("tix/data/", pathogen, "_", i, ".tsv"))) { next() }

            tmp <- tixdata(paste0("tix/data/", pathogen, "_", k, ".tsv"))

            colnames(tmp$data) <- rep("", ncol(tmp$data))
            
            datap <- cbind(tmp$data, data[, which(meta[which(rownames(meta) %in% "plate"), ] %in% k)])

            C <- which(colnames(datap) %in% "")

            start <- Sys.time()
            if (cores > 1) {
                sfInit(parallel = TRUE, cpus = cores)
                tmp <- sfLapply(1:nrow(datap), distrPar, datap, C)
                sfStop()
            } else {
                tmp <- lapply(1:nrow(datap), distrPar, datap, C)
            }
            tmp <- do.call("rbind", tmp)
            colnames(tmp) <- rep(i, ncol(tmp))
            if (exists("llr")) {
                llr <- cbind(llr, tmp)
            } else {
                llr <- tmp
            }
            
            end <- Sys.time()
            print(end - start)
            save(llr, file = paste0("tix/llr/pw_", j, ".rda"))
            
        }
        system(paste0("rm tix/data/", pathogen, "_", plates, ".tsv"))
    }
    system(paste0("rm tix/data/", pathogen, "_", rownames(kadjtmp), ".tsv"))
}

stop("done")
