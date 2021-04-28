Head <- function(x, n=10) {
    y <- x[1:n, 1:n]
    return(y)
}

invertRNA <- function(x) {
    y <- gsub("A", "B", gsub("C", "D", gsub("G", "H", gsub("T", "U", x))))
    y <- gsub("B", "T", gsub("D", "G", gsub("H", "C", gsub("U", "A", y))))
    return(y)
}

dc0 <- read.delim("data/perturbseq/GSE90063_dc0hr_umi_wt.txt")

dc0c2g <- read.table("data/perturbseq/GSE90063_RAW/GSM2396857_dc_0hr_cbc_gbc_dict.csv", fill = T, sep = ",")

tmp <- read.csv("data/perturbseq/GSE90063_RAW/GSM2396857_dc_0hr_cbc_gbc_dict_trunk.csv", header = F)

tmp <- as.vector(tmp)

dc0cells <- read.table("data/perturbseq/GSE90063_RAW/GSM2396857_dc_0hr_cellnames.csv")

dc0genes <- read.delim("data/perturbseq/GSE90063_RAW/GSM2396857_dc_0hr_genenames.csv")

cells <- gsub("_.*", "", gsub(".*,", "", as.character(dc0cells[[1]])))

genes <- gsub("_.*", "", gsub(".*,", "", as.character(dc0genes[[1]])))

tmp <- unlist(strsplit(as.character(dc0c2g[[2]][1]), ","))

gene <- unlist(strsplit(tmp[1], "_"))[2]

tmp <- gsub(" ", "", gsub("_.*", "", tmp))[-1]

which(tmp %in% colnames(dc0))

which(colnames(dc0) %in% tmp)

sum(tmp %in% cells)

data <- dc0

colnames(data) <- 1:ncol(data)

for (i in 1:length(dc0c2g[[1]])) {

    tmp <- unlist(strsplit(as.character(dc0c2g[[1]][i]), ","))

    gene <- unlist(strsplit(tmp[1], "_"))[2]
    
    tmp <- gsub(" ", "", gsub("_.*", "", tmp))[-1]
    
    colnames(data)[which(colnames(dc0) %in% tmp)] <- paste(colnames(data)[which(colnames(dc0) %in% tmp)], gene, sep = ".")

}

colnames(data)[grep("\\.", colnames(data))]

sum(colnames(dc0) %in% cells)

