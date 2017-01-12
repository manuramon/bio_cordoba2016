#
# Bio Cordoba practice session

library(ggplot2); theme_set(theme_bw())
library(dplyr)

# 1. Examine PLINK outputs
# missing data reports
imiss <- read.table("Output_12_13_2016/geno_cordoba_stats.imiss"
                    , header = TRUE, stringsAsFactors = TRUE)
lmiss <- read.table("Output_12_13_2016/geno_cordoba_stats.lmiss"
                    , header = TRUE, stringsAsFactors = TRUE)
lmiss_aDNA <- read.table("Output_12_13_2016/geno_cordoba_stats_aDNA.lmiss"
                    , header = TRUE, stringsAsFactors = TRUE)
lmiss_rDNA <- read.table("Output_12_13_2016/geno_cordoba_stats_rDNA.lmiss"
                    , header = TRUE, stringsAsFactors = TRUE)

barplot(1-imiss$F_MISS, names.arg = imiss$IID, 
        horiz = TRUE, las=1,
        col = ifelse(imiss$F_MISS>0.1, "tomato", "lightgray"))

quantile(lmiss_aDNA$F_MISS, probs = seq(0,1,.1), na.rm = TRUE)

# missigness per sample
table(lmiss$N_MISS)/nrow(lmiss)*100
cumsum(rev(table(lmiss$N_MISS))) # use to decided the --geno value in plink
table(lmiss_rDNA$N_MISS)/nrow(lmiss_rDNA)*100
table(lmiss_aDNA$N_MISS)/nrow(lmiss_aDNA)*100

op <- par(mfrow=c(3,1), mar=c(4,4,2,1))
hist(lmiss$F_MISS, main="ALL", xlab="miss_freq")
hist(lmiss_aDNA$F_MISS, main="aDNA", xlab="miss_freq")
hist(lmiss_rDNA$F_MISS, main="rDNA", xlab="miss_freq")
par(op)

op <- par(mfrow=c(4,1), mar=c(4,4,2,1))
plot(1-lmiss$F_MISS[lmiss$CHR==24], cex=.2, type="h", 
     xlab = "position (bp)", ylab="missigness (%)", main="ALL", 
     col="gray60", las=1)
plot(1-lmiss_aDNA$F_MISS[lmiss_aDNA$CHR==24], cex=.2, type="h",
     xlab = "position (bp)", ylab="missigness (%)", main="aDNA", 
     col="gray60", las=1)
plot(1-lmiss_rDNA$F_MISS[lmiss_rDNA$CHR==24], cex=.2, type="h",
     xlab = "position (bp)", ylab="missigness (%)", main="rDNA", 
     col="gray60", las=1)
miss_diff <- (1-lmiss_rDNA$F_MISS[lmiss_rDNA$CHR==24])-
    (1-lmiss_aDNA$F_MISS[lmiss_aDNA$CHR==24])
plot(abs(miss_diff), cex=.2, type="h",
     xlab = "position (bp)", ylab="missigness (%)", main="rDNA vs. aDNA", 
     col=ifelse(miss_diff<0, "tomato","gray60"), las=1)
par(op)


# allele freq
snpfrq <- read.table("Output_12_13_2016/geno_cordoba_stats.frq"
                     , header = TRUE, stringsAsFactors = TRUE)
snpfrq_aDNA <- read.table("Output_12_13_2016/geno_cordoba_stats_aDNA.frq"
                     , header = TRUE, stringsAsFactors = TRUE)
snpfrq_rDNA <- read.table("Output_12_13_2016/geno_cordoba_stats_rDNA.frq"
                          , header = TRUE, stringsAsFactors = TRUE)

op <- par(mfrow=c(3,1), mar=c(4,4,2,1))
hist(snpfrq$MAF, main="ALL", xlab="miss_freq")
hist(snpfrq_aDNA$MAF, main="aDNA", xlab="miss_freq")
hist(snpfrq_rDNA$MAF, main="rDNA", xlab="miss_freq")
par(op)

table(round(snpfrq_rDNA$MAF, 1))

op <- par(mfrow=c(3,1), mar=c(4,4,2,1))
barplot(table(round(snpfrq$MAF, 1)), main="ALL", space = 0.05)
barplot(table(round(snpfrq_aDNA$MAF, 1)), main="aDNA", space = 0.05)
barplot(table(round(snpfrq_rDNA$MAF, 1)), main="rDNA", space = 0.05)
par(op)

op <- par(mfrow=c(4,1), mar=c(4,4,2,1))
plot(snpfrq$MAF[snpfrq$CHR==24], cex=.2, type="h",
     xlab = "position (bp)", ylab="MAF", main="ALL")
plot(snpfrq_aDNA$MAF[snpfrq_aDNA$CHR==24], cex=.2, type="h",
     xlab = "position (bp)", ylab="MAF", main="aDNA")
plot(snpfrq_rDNA$MAF[snpfrq_rDNA$CHR==24], cex=.2, type="h",
     xlab = "position (bp)", ylab="MAF", main="rDNA")
plot(snpfrq_rDNA$MAF[snpfrq_rDNA$CHR==24] - snpfrq_aDNA$MAF[snpfrq_aDNA$CHR==24], 
     cex=.2, type="h",
     xlab = "position (bp)", ylab="MAF diff", main="rDNA vs. aDNA")
par(op)

# moving average
ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
op <- par(mfrow=c(2,1), mar=c(3,3,1,1))
kk <- ma(snpfrq_rDNA$MAF, n=6)
plot(kk, cex=.2, type="n")
points(kk, cex=.2, col=factor(snpfrq_rDNA$CHR))
kk <- ma(snpfrq_aDNA$MAF, n=6)
plot(kk, cex=.2, type="n")
points(kk, cex=.2, col=factor(snpfrq_aDNA$CHR))
par(op)


# HWE
snphwe <- read.table("Output_12_13_2016/geno_cordoba_stats.hwe",
                     header = TRUE, stringsAsFactors = FALSE)
names(snphwe)[7:8] <- c("o_het","e_het")
sum(as.numeric(sapply(strsplit(snphwe$GENO, "/"), "[[", 1))) # AA
sum(as.numeric(sapply(strsplit(snphwe$GENO, "/"), "[[", 2))) # Aa
sum(as.numeric(sapply(strsplit(snphwe$GENO, "/"), "[[", 3))) # aa

snphwe_aDNA <- read.table("Output_12_13_2016/geno_cordoba_stats_aDNA.hwe",
                     header = TRUE, stringsAsFactors = FALSE)
names(snphwe_aDNA)[7:8] <- c("o_het","e_het")
sum(as.numeric(sapply(strsplit(snphwe_aDNA$GENO, "/"), "[[", 1))) # AA
sum(as.numeric(sapply(strsplit(snphwe_aDNA$GENO, "/"), "[[", 2))) # Aa
sum(as.numeric(sapply(strsplit(snphwe_aDNA$GENO, "/"), "[[", 3))) # aa


## LD
library(MCMCpack)
con <- pipe("awk '{print $2}' Output_12_13_2016/geno_cordoba.fam")
ids <- scan(con, what = "char"); close(con)
reldat <- scan("Output_12_13_2016/geno_cordoba_LD.rel")
relmat <- xpnd(reldat)
dimnames(relmat) <- list(ids, ids)
relmat <- relmat[c(4,1:3,6,5,8,9:10,7),c(4,1:3,6,5,8,9:10,7)]

library(gplots)
library(reshape2)
heatmap.2(relmat/2)

ggplot(melt(relmat/2), aes(Var1, Var2)) + 
    geom_tile(aes(fill=value)) + 
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(x = "", y = "") 
relmat[lower.tri(relmat)] <- NA
ggplot(na.omit(melt(relmat/2)), aes(Var1, Var2)) + 
    geom_tile(aes(fill=value)) + 
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(x = "", y = "") +
    theme(axis.text.x=element_text(angle=90),
          axis.ticks=element_blank(),
          axis.line=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_line(color='#eeeeee'))

ldsnp <- read.table("Output_12_13_2016/geno_cordoba_LD2.ld",
                    header = TRUE, stringsAsFactors = FALSE)
plot(ldsnp$R2, cex=.4, pch=4, col=factor(ldsnp$CHR_A))
plot(ldsnp$R2[ldsnp$R2>.5], cex=.4, pch=4, 
     col=factor(ldsnp$CHR_A[ldsnp$R2>.5]))
ldsnp$dist <- round((ldsnp$BP_B - ldsnp$BP_A)/1000)
ldsnp$distc <- c(1:21)[cut(ldsnp$dist, 
                          c(seq(0,100,10),seq(150,500,50),750,1000), 
                          include.lowest = T)]
plot(aggregate(R2~distc, ldsnp, mean), 
     type="l", las=1, xlim=c(0,21),xaxt="n",
     xlab="Distance (Mb)", ylab="LD (r2)")
axis(1, at=c(1,6,11,14,19,21), labels = c(0,50,100,250,500,1000))


## PCA from plink
myeval <- read.table("Output_12_13_2016/geno_cordoba_pca.eigenval")
myevec <- read.table("Output_12_13_2016/geno_cordoba_pca.eigenvec",
                     stringsAsFactors = FALSE)
plot(myeval$V1, type="o", 
     ylab="eigenvalues", xlab="")
abline(h=1, col="orange")

pop_code <- c(rep("FOS",4),"OCP","FOS","ME","OCP","ME","ME")
plot(myevec[, 3:4],
     xlim=c(-.5, .5), ylim=c(-.5, .7),
     xlab="eigenvector 1", ylab="eigenvector 2",
     col=c("blue", "orange", "tomato")[factor(pop_code)],
     pch=c(17,16,16)[factor(pop_code)])
abline(h=0,v=0, lty=2)
legend("topright", legend = c("FOS","ME","OCP"), 
       col = c("blue", "tomato","orange"),
       pch=c(17,16,16), bty="n")

## SNPrelate
library(gdsfmt)
library(SNPRelate)

bed.fn <- "Output_12_13_2016/geno_cordoba.bed"
fam.fn <- "Output_12_13_2016/geno_cordoba.fam"
bim.fn <- "Output_12_13_2016/geno_cordoba.bim"
vcf.fn <- "Output_12_13_2016/geno_cordoba.vcf"

snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "geno_cordoba.gds")
# snpgdsVCF2GDS(vcf.fn, "geno_cordoba.gds", method = "biallelic.only")
snpgdsSummary("geno_cordoba.gds")

# Open the GDS file
genofile <- snpgdsOpen("geno_cordoba.gds")

# LD-based SNP pruning
set.seed(1234)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.1)
# Get all selected snp id
snpset.id <- unlist(snpset)

# PCA
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=3)
pop_code <- c(rep("FOS",4),"OCP","FOS","ME","OCP","ME","ME")
# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
tab$clase <- ifelse(substr(tab$sample.id, 1,5)=="0_FOS", 1, 0)
plot(tab$EV2, tab$EV1, 
     xlab="eigenvector 2", 
     ylab="eigenvector 1",
     col=c("blue", "orange", "tomato")[factor(pop_code)],
     pch=ifelse(tab$clase==1, 17, 16))
abline(h=0, v=0, lty=2)
legend("bottomleft", legend = c("FOS","ME","OCP"), 
       col = c("blue", "tomato","orange"),
       pch=c(17,16,16), bty="n")

# lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
# pairs(pca$eigenvect[,1:4], 
#       col=ifelse(tab$clase==1, "blue", "orange"),
#       upper.panel=NULL,
#       labels=lbls)

## genetic tree
dissMatrix <- snpgdsDiss(genofile , snp.id=snpset.id, num.thread=3, verbose=TRUE)
snpHCluster <-  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.25)
arbolG <- snpgdsCutTree(snpHCluster, 
                         z.threshold=15, outlier.n=5, 
                         n.perm = 5000, samp.group=NULL,
                         col.outlier="red", col.list=NULL, pch.outlier=4, 
                         pch.list=NULL,label.H=FALSE, label.Z=TRUE, verbose=TRUE)
# snpgdsDrawTree(arbolG, leaflab="perpendicular")
plot(arbolG$dendrogram, leaflab = "perpendicular",
     ylab="individual dissimilarity")
library("ape")
op <- par(mar=c(1,1,1,1), mgp=c(0,0,0), cex=1)
plot(as.phylo(snpHCluster$hclust), type = "fan",
     tip.color = c("blue", "orange", "tomato")[factor(pop_code)])
par(op)

# Parallel coordinates plot for the top principal components
# library(MASS)
# parcoord(pca$eigenvect, col=ifelse(tab$clase==1, "blue", "orange"))

# Fst estimation
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
snpgdsFst(genofile, sample.id=sample.id, 
          population=as.factor(pop_code),
          method="W&C84")



# ADMIXTURE
admix <- read.table("Output_12_13_2016/geno_cordoba.2.Q")
op <- par(mar=c(3,3,1,1), mgp=c(2,.7,0))
barplot(t(as.matrix(admix)), col=rainbow(2),
        space = 0.05,
        xlab="Individual #", ylab="Ancestry", border=NA)
par(op)

admix <- read.table("Output_12_13_2016/geno_cordoba_pruned.2.Q")
op <- par(mar=c(3,3,1,1), mgp=c(2,.7,0))
barplot(t(as.matrix(admix)), col=rainbow(2),
        names=ids, las=2,
        space = 0.05,
        xlab="Individual #", ylab="Ancestry", border=NA)
par(op)


## Venn diagram
library(VennDiagram)


## ROH
myroh <- read.table("Output_12_13_2016/geno_cordoba_roh.LROH", header = TRUE)
names(myroh)[6] <- "no_snp"
myroh$len <- round((myroh$MAX_END - myroh$MIN_START)/1000000)
table(myroh$INDV)
plot(no_snp ~ len, myroh, pch=20,
     xlab="Length of Segments (Mb)", ylab="Number of ROHs",
     col=ifelse(substr(myroh$INDV, 1,4)=="0_OC", "orange", "blue" ))
legend("bottomright", legend = c("ME","OCP"), 
       pch=20, col=c("blue","orange"), bty="n")


## CNVs
