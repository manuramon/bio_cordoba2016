#' ROH example
#' 

# Generate Ped and Map files

# Part 1: Generate Ped file
# Step 1: Generate Individual ID, Family ID, Paternal ID, Maternal ID, Sex, Phenotype columns

Fam.dat <- 1:500
Ind.dat <- 1:500
Pat.dat <- rep(0,500)
Mat.dat <- rep(0,500)
x <- 1:2
Sex.dat <- sample(x,500,replace=TRUE)
Phen.dat1 <- rep(1,250)
Phen.dat2 <- rep(2,250)
Phen.dat <- c(Phen.dat1,Phen.dat2)

ped.1t6<- cbind(Fam.dat,Ind.dat,Pat.dat,Mat.dat,Sex.dat,Phen.dat)

# Step 2: Generate genotype data
bases <- c("A","T","C","G")
nsnps <- 54609*2
nsnps <- 200000*2
base.nums <- sample(x,nsnps*500,replace=TRUE)
snp.matrix <- matrix(base.nums,nrow=500) 

i<- 1
snp.base.matrix <- matrix(nrow=nrow(snp.matrix),ncol=ncol(snp.matrix))
col.index <- matrix(1:(nsnps*4),nrow=nsnps,ncol=2,byrow=TRUE)
base.index <- matrix(nrow=nsnps,ncol=2)
replace.index <-which(snp.matrix[,col.index[i,]]==1)
replace.index2 <-which(snp.matrix[,col.index[i,]]==2)

for (i in 1:nsnps){
    base.index[i,] <- bases[sample(1:4,2,replace=FALSE)]
    snp.base.matrix[,col.index[i,]] <- replace(snp.matrix[,col.index[i,]],replace.index, base.index[i,][1])
    snp.base.matrix[,col.index[i,]] <- replace(snp.base.matrix[,col.index[i,]],replace.index2,base.index[i,][2])
}

# Step 3: Bind genotype data with first 6 columns

ped <- cbind(ped.1t6,snp.base.matrix)
write.table(ped,"data/artificial.ped",quote=FALSE,row.names=FALSE,col.names=FALSE)

# Part 2: Generate Map file

# Read in data from Map file to create current map:
map<-read.table("../Genotipos/Toros/gwas_toros.map", 
                head=FALSE, stringsAsFactors = FALSE)
names(map) <- c("chr","snp","cm","bp")
map <- subset(map, !chr %in% c("0","X","Y"))
map$chr <- as.numeric(map$chr)

library(dplyr)
setG <- map %>%  group_by(chr) %>% 
        summarise(n=n(), min=min(bp), max=max(bp)) %>% 
        data.frame()
setG$ng <- round(setG$n/sum(setG$n)*200000)
map <- c()
for (i in 1:nrow(setG)) {
    chr <- rep(i, length.out = setG$ng[i])
    pos <- sort(sample(setG$min[i]:setG$max[i], setG$ng[i]))
    print(c(setG$ng[i],length(chr), length(pos)))
    map <- rbind(map, cbind(chr, pos))
}
newmap <- cbind(chr=map[,"chr"],
                snp=paste0("snp", 1:nrow(map)),
                cm=0,
                bp=map[,"pos"])
newmap <- rbind(newmap, cbind(chr="29",snp="snp200000",cm="0",bp="51494790"))
write.table(newmap, "data/artificial.map", 
            quote=FALSE,row.names=FALSE,col.names=FALSE)


## ----------------------------

w.snp <- c(20,50,100,250,500)
for (i in w.snp){
    system(paste0("bin/plink_1.9_mac --cow --file data/artificial",
                  " --homozyg",
                  " --homozyg-window-snp ", i,
                  " --homozyg-snp 50 ",
                  " --homozyg-window-missing 5 --homozyg-window-threshold .05",
                  " --homozyg-window-het 5 --homozyg-density 50",
                  " --out data/artificial_", i))
}

# each file is read into R
roh20  <- read.table("data/artificial_20.hom.indiv",  header=TRUE)
roh50  <- read.table("data/artificial_50.hom.indiv",  header=TRUE)
roh100 <- read.table("data/artificial_100.hom.indiv", header=TRUE)
roh250 <- read.table("data/artificial_250.hom.indiv", header=TRUE)
roh500 <- read.table("data/artificial_500.hom.indiv", header=TRUE)

plot(NSEG ~ I(roh20$KB/1000), roh20,
     xlab="Length of Segments",
     ylab="Number of ROHs", 
     main="Runs of Homozygosity SNP window=20")



