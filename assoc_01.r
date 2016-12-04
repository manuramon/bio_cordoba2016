#' Theory behind association studies

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

table(rep(ped$Phenotype,2), c(ped$snp1a,ped$snp1b))
chisq.test(table(rep(ped$Phenotype,2), c(ped$snp1a,ped$snp1b)), correct = FALSE)

fisher.test(table(rep(ped$Phenotype,2), c(ped$snp1a,ped$snp1b)))

