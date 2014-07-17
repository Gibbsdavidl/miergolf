

# epistasis testing #
#The generalized R² was originally proposed by Cox & Snell (1989):
#1 - (L(0) / L(t)) ^ {2/n}

epinet <- function(net, dat, pheno, grp, whichgrp) {

    newtab <- data.frame()  
    idx <- grp == whichgrp
    datdat <- dat[,6:ncol(dat)]
    datdat <- datdat[,idx]
    dat <- cbind(dat[,1:5], datdat)
    nc <- ncol(dat)
    pheno <- pheno[idx]
  
    for (i in 1:nrow(net)) {
      print(i)
      try({
          geneA <- net[i,1]
          geneB <- net[i,2]
          
          if ( (geneA %in% dat$Gene) && (geneB %in% dat$Gene) ) {
              # as long as it's in there ... then search for good pair
              g1 <- which(dat$Gene == geneA)           # index to the variants
              g2 <- which(dat$Gene == geneB)
              sampP1 <- (dat$MAF[g1])/sum(dat$MAF[g1]) # sample according to MAF
              sampP2 <- (dat$MAF[g2])/sum(dat$MAF[g2])
              maxcod <- 0
              hillClimbingSteps <- 20000               # number of hill climbing steps
              at1 <- sample(g1,1)
              at2 <- sample(g2,1)                      # where we start
              g1dat <- as.numeric(dat[at1,6:nc])       # get the data out
              g2dat <- as.numeric(dat[at2,6:nc])
              count = 0                                # stop if no improvement
              d <- data.frame(Y=pheno,
                              G1=g1dat,
                              G2=g2dat)
              d$Y <- as.factor(d$Y)
                 
              for (step in 1:hillClimbingSteps) {
                  if (sample(x=c(0,1),1) == 1){ # probabilisically try one or the other
                      try1 <- sample(x=g1,size=1,prob=sampP1)
                      try2 <- at2
                      d$G1 <- as.numeric(dat[try1,6:nc])
                      d$G2 <- as.numeric(dat[try2,6:nc])
                  } else {
                      try1 <- at1
                      try2 <- sample(x=g2,size=1,prob=sampP2)
                      d$G1 <- as.numeric(dat[try1,6:nc])
                      d$G2 <- as.numeric(dat[try2,6:nc])
                  }
                  modI  <- glm(formula=Y~G1:G2, data=d, family=binomial(link=logit))
                  modN  <- glm(formula=Y~NULL, data=d, family=binomial(link=logit))
                  codI  <- 1- (exp(logLik(modN)[1]) / exp(logLik(modI)[1])) ^ (2/nrow(d))

                  if (codI > maxcod) {   # up the hill!
                      maxcod <- codI
                      at1 <- try1
                      at2 <- try2
                      count <- 0
                      cat (at1, at2, maxcod, step, "\n", sep="\t")
                  } else {               # no change
                      try1 <- at1
                      try2 <- at2
                      count <- count+1   # not getting better
                  }
                  if (count > 2000) {
                      print ("not getting better")
                      break
                  }
              }
              step
              at1
              at2
              maxcod

              g1dat <- as.numeric(dat[at1,6:nc])
              g2dat <- as.numeric(dat[at2,6:nc])
  
              d <- data.frame(Y=pheno,
                              G1=g1dat,
                              G2=g2dat)
              d$Y <- as.factor(d$Y)
    
              modI  <- glm(formula=Y~G1:G2, data=d, family=binomial(link=logit))
              modG1 <- glm(formula=Y~G1,    data=d, family=binomial(link=logit)) 
              modG2 <- glm(formula=Y~G2,    data=d, family=binomial(link=logit)) 
              modN  <- glm(formula=Y~NULL,  data=d, family=binomial(link=logit))
        
              codI  <- 1- (exp(logLik(modN)[1]) / exp(logLik(modI)[1]))  ^ (2/nrow(d))
              codG1 <- 1- (exp(logLik(modN)[1]) / exp(logLik(modG1)[1])) ^ (2/nrow(d))
              codG2 <- 1- (exp(logLik(modN)[1]) / exp(logLik(modG2)[1])) ^ (2/nrow(d))
        
              #if (codI > codG1 && codI > codG2) {
              #    print ("super sweet bro")
              newtab <- rbind(newtab, data.frame(at1, at2, geneA, geneB, codI, codG1, codG2,
                                                     modI$coefficients[1], modI$coefficients[2],
                                                     logLik(modI)[1]
                                                     ))
              #}
              #else {
              #    print ("negatory, dude.")
              #}
          }
      })
  }
    
  newtab
}


#####################################################################

library(doParallel)
library(foreach)
library(iterators)

registerDoParallel(cores=12)

load("groups.rda")

dat <- read.table(gzfile("haepic_data.txt.gz"), header=T, stringsAsFactors=F)
dat <- dat[dat$MAF > 0.01,]
bigmaf <- which(dat$MAF > 0.5)
dat$MAF[bigmaf] <- 1-dat$MAF[bigmaf]

haepic <- read.table("haepic_genes-regulate-genes.txt", stringsAsFactors=F)

nhaepic <- nrow(haepic)
chunk <- floor(13286/12)
idx <- sapply(1:11, function(a) rep(a,chunk), simplify=F)
idx[[12]] <- rep(12,nhaepic-(11*chunk))
netl <- split(haepic, as.factor(unlist(idx)))

placenta <- read.table("B_MRGE_Placenta_Related_NB____.txt")
placenta <- placenta[,2]

x <- foreach(i=1:12, .combine='rbind') %dopar% epinet(netl[[i]], dat, placenta, grp, 1)

########################################################################

library(doParallel)
library(foreach)
library(iterators)

registerDoParallel(cores=12)

load("groups.rda")

dat <- read.table(gzfile("haepic_data.txt.gz"), header=T, stringsAsFactors=F)
dat <- dat[dat$MAF > 0.01,]
bigmaf <- which(dat$MAF > 0.5)
dat$MAF[bigmaf] <- 1-dat$MAF[bigmaf]

haepic <- read.table("haepic_genes-regulate-genes.txt", stringsAsFactors=F)

nhaepic <- nrow(haepic)
chunk <- floor(13286/12)
idx <- sapply(1:11, function(a) rep(a,chunk), simplify=F)
idx[[16]] <- rep(12,nhaepic-(11*chunk))
netl <- split(haepic, as.factor(unlist(idx)))

eclampsia <- read.table("B_CLIN_Preeclampsia_Eclampsia_M____.txt")
eclampsia <- eclampsia[,2]

x <- foreach(i=1:12, .combine='rbind') %dopar% epinet(netl[[i]], dat, eclampsia, grp, 1)


####################################
