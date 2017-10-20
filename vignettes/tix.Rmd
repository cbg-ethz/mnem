## getting arguments:

system("rnai-query query --db /cluster/work/bewi/members/simondi/data/tix/database/tix_index.db --gene hspa5 --pathogen salmonella --study infectx --design p --sample 1000 temp.tsv")

stop()

pw <- as.numeric(commandArgs(TRUE)[1])

if (is.na(pw)) { pw <- NULL }

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

## ## testing 1,2,3

## data <- simData(Nems=1, Egenes = 100, reps = 5)$data

## res <- mynem(data, parallel = 10)

## save(res, file = "temp.rda")

## stop()

## pathogen <- "shigella"

## k <- "dz05-1e"

## cells <- 1000

## system(paste("rnai-query query --db /cluster/work/bewi/members/simondi/data/tix/database/tix_index.db --gene mock --plate ",
##              k,
##              " --pathogen ", pathogen, " --study infectx --design p --sample ", cells, " tix/data/", pathogen, "_", k, ".tsv", sep = ""))

## stop()
            
## here comes the R code for leo jobs:

tixdata <- function(x) {
    
    tmp <- t(read.delim(x))
    
    meta <- tmp[1:which(rownames(tmp) %in% "object_idx"), ]
    
    tmp <- apply(tmp[-(1:which(rownames(tmp) %in% "object_idx")), ], c(1,2), as.numeric)
    
    colnames(tmp) <- meta[which(rownames(meta) %in% "gene"), ]
    
    plates <- sort(unique(meta[which(rownames(meta) %in% "plate"), ]))
    
    return(list(data = tmp, meta = meta, plates = plates))

}

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
 
load("tix/kadj.rda")

pathogen <- "salmonella"

## system(paste0("rnai-query select gene --db /cluster/work/bewi/members/simondi/data/tix/database/tix_index.db --pathogen ", pathogen, " tix/genes_tmp.tsv"))

## genes <- read.delim("tix/genes_tmp.rda")

load("tix/genes.rda")

cores <- 64

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
            system(paste("rnai-query query --db /cluster/work/bewi/members/simondi/data/tix/database/tix_index.db --gene ",
                         tolower(i),
                         " --pathogen ", pathogen, " --study infectx --design p --sample ", cells, " tix/data/", pathogen, "_", i, ".tsv", sep = ""))
        } else { next() }

        if (!file.exists(paste0("tix/data/", pathogen, "_", i, ".tsv"))) { next() }

        tmp <- tixdata(paste0("tix/data/", pathogen, "_", i, ".tsv"))

        meta <- tmp$meta

        plates <- tmp$plates

        data <- tmp$data

        for (k in plates) {

            if (!file.exists(paste0("tix/data/", pathogen, "_", k, ".tsv"))) {
                system(paste("rnai-query query --db /cluster/work/bewi/members/simondi/data/tix/database/tix_index.db --gene ",
                             "mock",
                             " --plate ",
                             k,
                             " --pathogen ", pathogen, " --study infectx --design p --sample ", cells, " tix/data/", pathogen, "_", k, ".tsv", sep = ""))
            }

            ## if (!file.exists(paste0("tix/data/", pathogen, "_", i, ".tsv"))) { next() }

            ## tmp <- tixdata(paste0("tix/data/", pathogen, "_", k, ".tsv"))

            ## colnames(tmp$data) <- rep("", ncol(tmp$data))
            
            ## datap <- cbind(tmp$data, data[, which(meta[which(rownames(meta) %in% "plate"), ] %in% k)])

            ## C <- which(colnames(datap) %in% "")

            ## start <- Sys.time()
            ## sfInit(parallel = TRUE, cpus = cores)
            ## ##tmp <- sfLapply(1:4, densPar, datap, C) # for local test
            ## tmp <- sfLapply(1:nrow(datap), densPar, datap, C)
            ## sfStop()
            ## llr <- cbind(llr, do.call("rbind", tmp))
            ## end <- Sys.time()
            ## print(end - start)
            
        }
    }

    ## save(llr, file = paste0("tix/pw_", j, ".rda"))

    ## system("rm *.tsv")
    
}

stop("done")

## run job

stop()

source activate tix

module load r/3.3.3

## old:

bsub -M 10000 -q normal.4h -n 1 -R "rusage[mem=10000]" "R --slave < mnem/vignettes/tix.r > output.txt"

bsub -M 100000 -W 10 -n 1 -R "rusage[mem=100000]" "R --slave < mnem/vignettes/tix.r > output.txt"

## batch:

for i in {1..225}
do

bsub -M 10000 -q normal.4h -n 1 -e logs/error$i.txt -o logs/output$i.txt -R "rusage[mem=10000]" "R --silent --no-save --args $i < mnem/vignettes/tix.r"

sleep 1
done

## single:

rm error.txt

rm output.txt

bsub -M 10000 -q normal.4h -n 1 -e error.txt -o output.txt -R "rusage[mem=10000]" "R --silent --no-save < mnem/vignettes/tix.r"

##

rm error.txt

rm output.txt

bsub -M 10000 -q normal.4h -n 1 -e error.txt -o output.txt -R "rusage[mem=10000]" "R --silent --no-save --args 1 < mnem/vignettes/tix.r"

stop()

## here comes addtional tix R code, which is not run on any grid:

source("https://bioconductor.org/biocLite.R")

biocLite(c("nem", "graph", "cluster", "snow", "snowfall", "epiNEM", "KEGG.db", "KEGGgraph", "Rgraphviz", "dplyr", "annotables"))

##

source("mnem/R/mnems.r")
source("mnem/R/mnems_low.r")

library(nem)
library(graph)
library(cluster)
library(bnem)
library(snow)
library(snowfall)
library(epiNEM)

library(KEGG.db)
library(KEGGgraph)
library(Rgraphviz)
library(dplyr)
library(annotables)

## get all pathways:

## system("rnai-query select --db /net/bs-filesvr05/export/group/beerenwinkel/simon/data/target_infect_x/data_base/tix.db pathogen > pathogen.tsv")

setwd("~/Mount/Leo/")

## system("rnai-query select --db /cluster/work/bewi/members/simondi/data/tix/database/tix_index.db pathogen > pathogen.tsv")

pathogen <- read.delim("pathogen.tsv", header = FALSE)

pathogen <- as.character(pathogen[, 1])

p <- pathogen[1]

## system("rnai-query select --db /cluster/work/bewi/members/simondi/data/tix/database/tix_index.db gene > genes.tsv")

genes <- read.delim("genes.tsv", header = TRUE)

genes <- toupper(as.character(genes[, 1]))

genes2 <- character(length(genes))

ids <- numeric(length(genes))

count <- 0

for (i in 1:length(genes)) {

    idx <- which(grch38$symbol %in% genes[i])

    if (length(idx) == 0) { next() }

    if (length(idx) > 1) {

        idx <- (idx[!is.na(idx)])[1]

    }

    id <- grch38$entrez[idx]

    ids[count+1] <- id

    genes2[count+1] <- genes[i]

    count <- count + 1

}

genes <- cbind(genes2, ids)

genes <- genes[order(genes[, 1]), ]

genes <- genes[-which(genes[, 2] %in% "0"), ]

genes <- genes[!is.na(genes[, 2]), ]

genes <- genes[!duplicated(genes[, 2]), ]

## save(genes, file = "tix/genes.rda")

load("tix/genes.rda")

kegg <- list()

for (i in genes[, 2]) {
    
    try(kegg[[i]] <- mget(i, KEGGEXTID2PATHID))

}

pathways <- sort(unique(unlist(kegg)))

for (i in pathways) {

    tryres <- try(retrieveKGML(pathwayid=gsub("hsa", "", i), organism="hsa", destfile=paste("tix/pathways/", i, ".xml", sep = ""), method="internal", quiet=TRUE), silent = TRUE)

}

kadj <- list()

count <- 0

for (i in pathways) {

    ## i <- pathways[1]

    if (!file.exists(paste("tix/pathways/", i, ".xml", sep = ""))) { next() }

    keggpw <- parseKGML(paste("tix/pathways/", i, ".xml", sep = ""))
    
    kgraph <- KEGGpathway2Graph(keggpw)
    
    kadjtmp <- graph2adj(kgraph)

    kadjtmp <- kadjtmp[which(gsub("hsa:", "", rownames(kadjtmp)) %in% genes[, 2]), which(gsub("hsa:", "", rownames(kadjtmp)) %in% genes[, 2]), drop = FALSE]

    print(dim(kadjtmp))

    if (dim(kadjtmp)[1] == 1) { next() }

    if (any(kadjtmp == 1)) {
        kadjtmp <- transitive.closure(kadjtmp, mat = TRUE)
    }
    
    diag(kadjtmp) <- 0

    ## I'll do this later:

    ## edgen <- order(apply(kadjtmp, 1, sum)+apply(kadjtmp, 2, sum), decreasing = TRUE)

    ## kadjtmp <- kadjtmp[edgen[1:min(nrow(kadjtmp), 5)], edgen[1:min(nrow(kadjtmp), 5)]]

    kadjtmp <- kadjtmp[order(rownames(kadjtmp)), order(rownames(kadjtmp))]

    rownames(kadjtmp) <- colnames(kadjtmp) <- genes[which(genes[, 2] %in% gsub("hsa:", "", rownames(kadjtmp))), 1]

    count <- count + 1
    
    kadj[[count]] <- kadjtmp

    names(kadj)[count] <- i

    save(kadj, file = "tix/kadj.rda")

}

p <- pathogen[2]

## old but maybe usefull:

## gene ids part II:

genes <- toupper(colnames(scorec)[order(scorec[3, ]*(-scorec[1, ]), decreasing = TRUE)[1:ngenes]])

library("biomaRt")
ensembl=useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

genesBM <- 
getBM(attributes=c('ensembl_gene_id', 'entrezgene', 'hgnc_symbol', 'hgnc_id'), 
      filters = 'hgnc_symbol', 
      values = genes, 
      mart = ensembl)

genesBM <- genesBM[-which(is.na(genesBM[, 2]) == TRUE), ]

## reactome:

library(org.Hs.eg.db)
library(ReactomePA)
x <- enrichPathway(gene=genes[, 2], pvalueCutoff=0.05, readable=T)
x <- as.data.frame(x)
head(x)

pathways <- unique(sort(x$Description))

viewPathway(pathways[1])
