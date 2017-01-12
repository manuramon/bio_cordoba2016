#' Theory behind association studies


#' Linkage Desequilibrium
#' example
ld_ex01 <- c("AB","AB","ab","aB","ab","ab","Ab","AB","Ab","AB")
library(stringr)
sum(str_count(ld_ex01, "A"))/10
sum(str_count(ld_ex01, "a"))/10
sum(str_count(ld_ex01, "B"))/10
sum(str_count(ld_ex01, "b"))/10

# observed freqs
sum(str_count(ld_ex01, "AB"))/10
sum(str_count(ld_ex01, "aB"))/10
sum(str_count(ld_ex01, "Ab"))/10
sum(str_count(ld_ex01, "ab"))/10

# Expected freqs: freq(AB) = freq(A) x freq(B)
sum(str_count(ld_ex01, "A"))/10 * sum(str_count(ld_ex01, "B"))/10
sum(str_count(ld_ex01, "a"))/10 * sum(str_count(ld_ex01, "B"))/10
sum(str_count(ld_ex01, "A"))/10 * sum(str_count(ld_ex01, "b"))/10
sum(str_count(ld_ex01, "a"))/10 * sum(str_count(ld_ex01, "b"))/10



#' PLINK files
#' MAP file
map <- data.frame(chromosome = 1,
                  snp = "snp1",
                  cm = 0,
                  bp = 1000
)
write.table(map, "data/assoc_01.map",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#' PED file
ped <- data.frame(FamilyID = "FAM001",
                  IndividualID = 1:250,
                  PaternalID = 0,  # unknown
                  MaternalID = 0,  # unknown
                  Sex = 0          # unknown
                  )
#' Example from Irizarry
disease <- factor(c(rep(0,180),rep(1,20),rep(0,40),rep(1,10)),
                  labels=c("control","cases"))
genotype <- factor(c(rep("AA/Aa",200), rep("aa",50)),
                   levels=c("AA/Aa","aa"))
dat <- data.frame(disease, genotype)
dat <- dat[sample(nrow(dat)),]     #shuffle them up head(dat)
ped$Phenotype <- as.numeric(dat$disease)
snp1 <- as.character(dat$genotype)
snp1[snp1=="aa"] <- "GG"
snp1[snp1=="AA/Aa"] <- sample(c("AA","AG"), prob = c(1/3, 2/3), 
                                  size = 200, replace = TRUE)
ped$snp1a <- substr(snp1, 1, 1)
ped$snp1b <- substr(snp1, 2, 2)
write.table(ped, "data/assoc_01.ped",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#' Chisq test
tab <- table(genotype, disease)
tab
# odds of having the disease of you are "aa":     10/40 = 0.25
# odds of having the disease of you are "AA/Aa": 20/180 = 0.1111
OR <- (tab[2,2]/tab[2,1])/(tab[1,2]/tab[1,1])
assoc_test <- chisq.test(tab, correct = FALSE)
assoc_test$observed
assoc_test$expected
assoc_test$p.value


#' With plink
system("sh assoc_01.sh")


#' With SNPassoc in R
library(SNPassoc)
mySNP <- snp(snp1, sep="")
summary(mySNP)
ped2 <- ped[, 1:6]
ped2$snp1 <- snp1
pos <- map[, c(2,1,4)]
names(pos) <- c("snp","chr","pos")
ped2 <- setupSNP(ped2, colSNPs = 7, sep="", info = pos)
tableHWE(ped2)
gwas01 <- WGassociation(Phenotype ~ snp1, data=ped2)
gwas01
association(Phenotype ~ snp1, data=ped2)

## Large samples, small p-values
tab2 <- tab*10
chisq.test(tab2, correct = FALSE)


## Manhattan plots
library(qqman)
head(gwasResults) # see data
str(gwasResults)  # what's in data?
as.data.frame(table(gwasResults$CHR)) # num. markers per chr

# wrong plot
plot(gwasResults$P, cex=.2, pch="+")

# plot -log10(P)
manhattan(gwasResults, 
          ylim = c(0, 10), cex = 0.6, cex.axis = 0.9, 
          suggestiveline = F, genomewideline = F, 
          chrlabs = as.character(c(1:22)))

# highlight SNPs of interest
str(snpsOfInterest)
manhattan(gwasResults, 
          ylim = c(0, 10), cex = 0.6, cex.axis = 0.9, 
          col = c("blue4", "orange3"), 
          suggestiveline = F, genomewideline = F, 
          highlight = snpsOfInterest,
          chrlabs = as.character(c(1:22)))

# plot a specific CHR
manhattan(subset(gwasResults, CHR==3),
          ylim = c(0, 10), #xlim = c(200,500),
          cex = 0.6, cex.axis = 0.9,
          col = c("blue4", "orange3"),
          suggestiveline = F, genomewideline = F,
          highlight = snpsOfInterest)

# QQ-plot
qq(gwasResults$P)


## plot other scores
gwasResults <- transform(gwasResults, 
                         zscore = qnorm(P/2, lower.tail = FALSE))
head(gwasResults)
manhattan(gwasResults, 
          p = "zscore", logp=FALSE, ylab="Z-score",
          ylim = c(0, 7), cex = 0.6, cex.axis = 0.9, 
          col = c("blue4", "orange3"), 
          suggestiveline = F, genomewideline = F, 
          highlight = snpsOfInterest,
          chrlabs = as.character(c(1:22)))

library(qvalue)
myqval <- qvalue(gwasResults$P, fdr.level = 0.01)
gwasResults <- transform(gwasResults, 
                         qvalue = myqval$qvalues)
manhattan(gwasResults, 
          p = "qvalue", logp=TRUE, 
          ylab=bquote(-log[10]~Q-value),
          ylim = c(0, 10), cex = 0.6, cex.axis = 0.9, 
          col = c("blue4", "orange3"), 
          suggestiveline = F, genomewideline = F, 
          highlight = snpsOfInterest,
          chrlabs = as.character(c(1:22)))

