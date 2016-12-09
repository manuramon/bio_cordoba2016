#
# Bio Cordoba practice session

library(ggplot2); theme_set(theme_bw())

# 1. Examine PLINK outputs
# missing data reports
imiss <- read.table("Output_12_09_2016/geno_cordoba_stats.imiss"
                    , header = TRUE, stringsAsFactors = TRUE)
lmiss <- read.table("Output_12_09_2016/geno_cordoba_stats.lmiss"
                    , header = TRUE, stringsAsFactors = TRUE)
lmiss_aDNA <- read.table("Output_12_09_2016/geno_cordoba_stats_aDNA.lmiss"
                    , header = TRUE, stringsAsFactors = TRUE)
lmiss_rDNA <- read.table("Output_12_09_2016/geno_cordoba_stats_rDNA.lmiss"
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
snpfrq <- read.table("Output_12_09_2016/geno_cordoba_stats.frq"
                     , header = TRUE, stringsAsFactors = TRUE)
snpfrq_aDNA <- read.table("Output_12_09_2016/geno_cordoba_stats_aDNA.frq"
                     , header = TRUE, stringsAsFactors = TRUE)
snpfrq_rDNA <- read.table("Output_12_09_2016/geno_cordoba_stats_rDNA.frq"
                          , header = TRUE, stringsAsFactors = TRUE)

op <- par(mfrow=c(3,1), mar=c(4,4,2,1))
hist(snpfrq$MAF, main="ALL", xlab="miss_freq")
hist(snpfrq_aDNA$MAF, main="aDNA", xlab="miss_freq")
hist(snpfrq_rDNA$MAF, main="rDNA", xlab="miss_freq")
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

