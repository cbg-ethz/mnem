
args <- commandArgs()
Sgenes <- as.numeric(gsub("Sgenes=", "", args[grep("Sgenes=", args)]))
run <- gsub("run=", "", args[grep("run=", args)])
noise <- gsub("noise=", "", args[grep("noise=", args)])
nem <- gsub("nem=", "", args[grep("nem=", args)])

# nem <- 1; noise <- run <- 1; Sgenes <- 20

if (length(grep(":", nem)) > 0) {
    nem <- as.numeric(gsub(":.*", "", nem)):as.numeric(gsub(".*:", "", nem))
} else {
    nem <- as.numeric(nem)
}

if (length(grep(":", noise)) > 0) {
    noise <- as.numeric(gsub(":.*", "", noise)):as.numeric(gsub(".*:", "", noise))
} else {
    noise <- as.numeric(noise)
}

if (length(grep(":", run)) > 0) {
    run <- as.numeric(gsub(":.*", "", run)):as.numeric(gsub(".*:", "", run))
} else {
    run <- as.numeric(run)
}

library(cluster)
library(nem)
## library(epiNEM)
library(snowfall)
library(Rgraphviz)
source("mnem/R/mnems.r")
source("mnem/R/mnems_low.r")

start1 <- as.numeric(format(Sys.time(), "%s"))

## variable parameters:
runs <- 100
noises <- c(1, 2.5, 5, 10)
nems <- 1:10
maxk <- 5

## fixed parameters:
starts <- 10 # make it dependent on the learnt k
search <- "modules"
verbose <- FALSE

## mixing parameters:
Egenes <- 2
nCells <- 1000

simres <- array(0, c(runs, length(noises), length(nems), 4, 6), list(paste("run_", 1:runs, sep = ""), paste("noise_", noises, sep = ""), paste("components_", nems, sep = ""), c("mnem", "nem", "random", "random2"), c("time", "overfit", "accuracy", "sensitivity", "specificity", "mixing")))

for (i in nem) {
    for (j in noise) {    
        if (file.exists(paste("temp/simres_mnem_", Sgenes, "_", run, "_", j, "_", i, ".rda", sep = ""))) {
            print(paste("temp/simres_mnem_", Sgenes, "_", run, "_", j, "_", i, ".rda", sep = ""))
            stop("simulation result already exists")
        }
    }
}

for (donoise in noise) {
    for (donem in nem) {

        ## donoise <- noise[1]; donem <- nem[1]

        print(paste(rep("_", 100), collapse = ""))
        print(paste("run", run))
        print(paste("noise", noises[donoise]))
        print(paste("nems", nems[donem]))
 
        mw <- runif(nems[donem], 0.1, 1)
        mw <- mw/sum(mw)
        sim <- simData(Sgenes = Sgenes, Egenes = Egenes, nCells = nCells,
                       Nems = nems[donem],
                       mw = mw, uninform = floor(Sgenes*Egenes*0.1))
        simfull <- NULL
        fullsim <- sim$Nem[[1]]*0
        for (i in 1:length(sim$Nem)) {
            tmp <- transitive.closure(sim$Nem[[i]], mat = TRUE)
            simfull <- cbind(simfull, t(tmp))
            fullsim <- fullsim + transitive.closure(sim$Nem[[i]], mat = TRUE)
        }
        fullsim[which(fullsim > 1)] <- 1
        diag(fullsim) <- 1
        
        ## mnem:
        data <- sim$data
        data <- (data - 0.5)/0.5
        data <- data + rnorm(length(sim$data), 0, noises[donoise])
        ## learn k:
        ## random p:
        p <- NULL
        start <- as.numeric(format(Sys.time(), "%s"))
        res <- list()
        if (maxk == 1) {
            res <- mnem(data, starts = starts, search = search, verbose = verbose)
        } else {
            bics <- rep(Inf, maxk)
            for (k in 1:maxk) {
                res[[k]] <- mnem(data, starts = starts, search = search, k = k, verbose = verbose)
                bics[k] <- getIC(res[[k]], AIC = FALSE)
            }
            res1 <- res[[1]]
            res <- res[[which.min(bics)]]
        }
        simres[run, donoise, donem, 1, 1] <- as.numeric(format(Sys.time(), "%s")) - start
        resfull <- NULL
        fullres <- res$comp[[1]]$phi*0
        for (i in 1:length(res$comp)) {
            tmp <- transitive.closure(res$comp[[i]]$phi, mat = TRUE)
            resfull <- cbind(resfull, t(tmp))
            fullres <- fullres + transitive.closure(res$comp[[i]]$phi, mat = TRUE)
        }
        fullres[which(fullres > 1)] <- 1
        diag(fullres) <- 1
        
        simres[run, donoise, donem, 1, 2] <- length(res$comp)/nems[donem]
        tp <- sum(fullres == 1 & fullsim == 1) - ncol(fullres)
        tn <- sum(fullres == 0 & fullsim == 0)
        fp <- sum(fullres == 1 & fullsim == 0)
        fn <- sum(fullres == 0 & fullsim == 1)
        simres[run, donoise, donem, 1, 4] <- tp/(tp+fn)
        simres[run, donoise, donem, 1, 5] <- tn/(tn+fp)
        simres[run, donoise, donem, 1, 3] <- hamSim(simfull, resfull)
        if (length(mw) == length(res$mw)) {
            simres[run, donoise, donem, 1, 6] <- sum(dist(rbind(sort(res$mw), sort(mw))))
        }
        if (length(mw) > length(res$mw)) {
            simres[run, donoise, donem, 1, 6] <- sum(dist(rbind(sort(c(res$mw, rep(0, length(mw) - length(res$mw)))), sort(mw))))
        }
        if (length(mw) < length(res$mw)) {
            simres[run, donoise, donem, 1, 6] <- sum(dist(rbind(sort(c(mw, rep(0, length(res$mw) - length(mw)))), sort(res$mw))))
        }
        
        ## nem:
        start <- as.numeric(format(Sys.time(), "%s"))
        nemres <- mynem(data, search = search)
        simres[run, donoise, donem, 2, 1] <- as.numeric(format(Sys.time(), "%s")) - start
        if (maxk == 1) {
            fullnem <- transitive.closure(nemres$adj, mat = TRUE)
        } else {
            fullnem <- transitive.closure(res1$comp[[1]]$phi, mat = TRUE)
        }
        diag(fullnem) <- 1

        simres[run, donoise, donem, 2, 2] <- 1/nems[donem]
        tp <- sum(fullnem == 1 & fullsim == 1) - ncol(fullnem)
        tn <- sum(fullnem == 0 & fullsim == 0)
        fp <- sum(fullnem == 1 & fullsim == 0)
        fn <- sum(fullnem == 0 & fullsim == 1)
        simres[run, donoise, donem, 2, 4] <- tp/(tp+fn)
        simres[run, donoise, donem, 2, 5] <- tn/(tn+fp)
        simres[run, donoise, donem, 2, 3] <- hamSim(simfull, t(fullnem))
        simres[run, donoise, donem, 2, 6] <- sum(dist(rbind(sort(c(1, rep(0, length(mw) - 1))), sort(mw))))
        
        ## random:
        fullrand <- matrix(sample(c(0,1), Sgenes^2, replace = TRUE), Sgenes, Sgenes)
        diag(fullrand) <- 1
        rownames(fullrand) <- colnames(fullrand) <- 1:Sgenes
        
        simres[run, donoise, donem, 4, 2] <- 1/nems[donem]
        tp <- sum(fullrand == 1 & fullsim == 1) - ncol(fullrand)
        tn <- sum(fullrand == 0 & fullsim == 0)
        fp <- sum(fullrand == 1 & fullsim == 0)
        fn <- sum(fullrand == 0 & fullsim == 1)
        simres[run, donoise, donem, 4, 4] <- tp/(tp+fn)
        simres[run, donoise, donem, 4, 5] <- tn/(tn+fp)
        simres[run, donoise, donem, 4, 3] <- hamSim(simfull, fullrand)
        simres[run, donoise, donem, 4, 6] <-  sum(dist(rbind(sort(c(1, rep(0, length(mw) - 1))), sort(mw))))
        
    }
}

end1 <- as.numeric(format(Sys.time(), "%s"))

print(end1 - start1)

## simres[1,,,,]

for (i in nem) {
    for (j in noise) {    
        save(simres, noises, nems, file = paste("temp/simres_mnem_", Sgenes, "_", run, "_", j, "_", i, ".rda", sep = ""))
    }
}

## save(simres, noises, file = paste("temp/simres_mnem_", Sgenes, "_", run, "_", noise, "_", nem, "_", as.numeric(format(Sys.time(), "%s")), ".rda", sep = ""))

stop("simulations are done")
