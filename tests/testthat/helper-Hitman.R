library(covr)
library(limma) #for EList-class
library(ppcor)
library(Hitman)
library(testthat)

# example from limma::contrasts.fit
set.seed(0)
M <- matrix(rnorm(100*6, sd=0.3), nrow=100, ncol=6)
dimnames(M) <- list(paste0("gene", 1:nrow(M)), paste0("sample", 1:ncol(M)))
design <- cbind(First3Arrays=c(1,1,1,0,0,0), Last3Arrays=c(0,0,0,1,1,1))
grp <- rep(c("First3", "Last3"), each=3)
M[1, 1:3] <- M[1, 1:3] + 2

el <- new("EList")
el$E <- M
ww <- matrix(rexp(n=nrow(M)*ncol(M)), ncol=ncol(M), nrow=nrow(M))
el$weights <- ww

phenotype <- rnorm(ncol(M))
names(phenotype) <- colnames(M)
pheno.minus.one <- phenotype
pheno.minus.one[1:3] <- pheno.minus.one[1:3]-1
pheno2 <- pheno.minus.one
pheno2[1] <- NA
covar <- rnorm(length(pheno.minus.one))

grp.mat <- ezlimma:::batch2design(grp)

# create associated phenotype, to avoid hitman warning about weak assoc
set.seed(0)
pheno.v <- stats::setNames(rnorm(ncol(M)), nm=colnames(M))
pheno.v[1:3] <- pheno.v[1:3]-3
ee <- pheno.v + rnorm(length(pheno.v), sd=0.1)
grp2 <- ezlimma:::batch2design(grp)[,1]
names(grp2) <- colnames(M)
covar.tmp <- rnorm(length(pheno.v))

hm <- hm.tmp <- hitman(E=grp2, M=M, Y=M[1,], check.names = FALSE)

hm.tmp$MY_dir.p <- hm.tmp$MY.p
hm.tmp$EM_dir.p <- hm.tmp$EM.p
ey.sign <- sign(cor(M[1,], grp2))
p.cols <- c("EM_dir.p", "MY_dir.p")
# rerun this on hm
ret <- modify_hitman_pvalues(tab=hm.tmp, overall.sign = ey.sign, stat.cols=c("EM.z", "MY.z"), p.cols=p.cols)

# replication too
set.seed(0)
tab.tmp <- data.frame(matrix(NA, nrow=100, ncol=4, dimnames=list(paste0("r", 1:100), c("stat1", "p1", "stat2", "p2"))))
tab.tmp[, c(1, 3)] <- rnorm(n=200)
tab.tmp[, 2] <- 2*pnorm(-abs(tab.tmp[, 1]))
tab.tmp[, 4] <- 2*pnorm(-abs(tab.tmp[, 3]))
