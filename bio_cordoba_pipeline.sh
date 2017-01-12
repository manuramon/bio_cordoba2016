#!/bin/sh

## Create files from Genome Studio Reports to PLINK
## Manuel Ramon, 18-abr-2016

## Initialize some vars
HEAP=1024 # BEAGLE Heap Size

outdir=Output_$(date +"%m_%d_%Y")
mkdir $outdir
infile="data/snp50K/Cordoba_FinalReport"
outfile="geno_cordoba"
breed="sheep"
beagle_dir="bin"
bin_dir="bin"


# HEADER
echo "\n\n  Running Genome Studio 2 PLINK script."
echo "  "$(date)"\n"
echo "  Breed:" $breed
echo "  Data file:" $infile.txt
echo "  Output files:" $outfile.map, $outfile.ped
echo ""

# SAMPLES FILE


# MAP FILE
  # Column 1 - Chr: Chr ($19)
  # Column 2 - SNP: SNP Name ($1)
  # Column 3 - Genetic distance: 0
  # Column 4 - Base-pair position: Position ($20)
awk 'BEGIN {FS="\t"; OFS="\t"}
     NR==11 {sample_ID=$2}; 
     NR>10 {if($2!=sample_ID) {exit 1}; 
            print $19,$1,0,$20}' $infile.txt  > $outdir/$outfile.map

# LGEN (PED) FILE
# SNP_name, Sample_ID, Allele1_AB, Allele2_AB and GC_Score 
  # Column 1 - Family ID: 0 or Sample Group ($7)
  # Column 2 - Individual ID: Sample ID ($2)
  # Column 3 - SNP: SNP Name ($1)
  # Column 4 - A1: Allele1-Top ($3)
  # Column 5 - A2: Allele2-Top ($4)
  # Column 6 - GC Score: GC.Score ($5) - NOT 
 

# FAM FILE
  # Column 1 - Family ID: 0 or Sample Group ($7)
  # Column 2 - Individual ID: Sample ID ($2)
  # Column 3 - Paternal ID:  0 (from pedigree)
  # Column 4 - Maternal ID:  0 (from pedigree)
  # Column 5 - Sex:  0 (from pedigree)
  # Column 6 - Phenotype:  0 (from data)
awk 'BEGIN {FS="\t"; OFS="\t"}
     NR>10 {print 0,$2,0,0,0,0}' $infile.txt  > $outdir/temp
uniq $outdir/temp > $outdir/$outfile.fam
rm $outdir/temp


# translate missing genotypes from "- -" to "0 0" 
awk -v OFS='\t' '{sub(/\t-\t-/, "\t0\t0"); print}' $outdir/$outfile.temp.lgen > $outdir/$outfile.lgen
rm $outdir/$outfile.temp.lgen


# PLINK. FIRST STEPS

# IMP!!! Use --keep-allele-order to avoid changes of ref alleles
# .bim file contains chromosome, marker name, genetic position, physical position, 
# and allele 1 (minor?) and allele 2 (major?).

# IMP: update map information
# Oar_v4.0 downloaded from NPchiMp (http://bioinformatics.tecnoparco.org/SNPchimp/)
awk  'BEGIN{FS=","; OFS="\t"}; NR>1{print $12,$11}' data/SNPchimp_oarv4_0.csv > $outdir/oar_v4_bp_positions.txt
$bin_dir/plink_1.9_mac --noweb --sheep --lfile  $outdir/$outfile \
    --allow-extra-chr 0 --biallelic-only \
    --update-map $outdir/oar_v4_bp_positions.txt\
    --keep-allele-order \
    --make-bed \
    --out $outdir/${outfile}_v4

$bin_dir/plink_1.9_mac --noweb --sheep --lfile  $outdir/$outfile \
    --allow-extra-chr 0 --biallelic-only \
    --keep-allele-order \
    --make-bed \
    --out $outdir/${outfile}_v4
    
$bin_dir/plink_1.9_mac --noweb --sheep --bfile  $outdir/${outfile}_v4 \
    --make-bed \
    --out $outdir/${outfile}

rm $outdir/${outfile}_v4*


# transpose data
$bin_dir/plink_1.9_mac --noweb --sheep --bfile  $outdir/${outfile} \
    --recode transpose \
    --out $outdir/${outfile}_trans



# 1st check data 
$bin_dir/plink_1.9_mac --noweb --sheep --bfile  $outdir/${outfile} \
    --allow-extra-chr 0 \
    --freq --missing --hardy --het small-sample --ibc --nonfounders\
    --out $outdir/${outfile}_stats
    
# do the same but spliting aDNA
awk '{if(substr($2,1,3)=="FOS") print $1,$2}' $outdir/${outfile}.fam > $outdir/aDNA_samples
$bin_dir/plink_1.9_mac --noweb --sheep --bfile  $outdir/$outfile \
    --allow-extra-chr 0 --biallelic-only \
    --keep $outdir/aDNA_samples \
    --freq --missing --hardy --het small-sample  --ibc --nonfounders\
    --out $outdir/${outfile}_stats_aDNA

awk '{if(substr($2,1,3)!="FOS") print $1,$2}' $outdir/${outfile}.fam > $outdir/rDNA_samples
$bin_dir/plink_1.9_mac --noweb --sheep --bfile  $outdir/$outfile \
    --allow-extra-chr 0 --biallelic-only \
    --keep $outdir/rDNA_samples \
    --freq --missing --hardy --het small-sample  --ibc --nonfounders\
    --out $outdir/${outfile}_stats_rDNA
    
rm ${outdir}/*.nosex
rm ${outdir}/*.log


# 2nd Quality Control (QC) 

# If we make QC using default thresholds with aDNA and rDNA data 
# we will lose many SNPS and also samples (see plink log)
$bin_dir/plink_1.9_mac --sheep --bfile  $outdir/$outfile \
    --allow-extra-chr --biallelic-only \
    --geno 0.1  --mind 0.1 --maf 0.001 --missing \
    --keep-allele-order \
    --make-bed \
    --out $outdir/${outfile}_temp

# IMP: check missing data reports!!!
# remove Chr 0, X, Y
awk '$1==0 || $1=="CONTIG" || $1=="OAR" {print $2}' $outdir/$outfile.map > $outdir/excludeSNP

# QC
$bin_dir/plink_1.9_mac --noweb --sheep --bfile  $outdir/$outfile \
    --biallelic-only \
    --geno 0.7 --mind 1 --maf 0.0001 --missing \
    --exclude $outdir/excludeSNP \
    --keep-allele-order \
    --make-bed \
    --out $outdir/$outfile

rm $outdir/${outfile}_temp*

# 3rd recode as VCF file

$bin_dir/plink_1.9_mac --noweb --sheep --bfile  $outdir/$outfile \
    --allow-extra-chr 0 \
    --recode vcf \
    --out $outdir/$outfile

# there is a row given and advise that must be deleted from VCF file
sed '/may not be based on real reference genome/d' $outdir/$outfile.vcf > $outdir/temp.vcf
mv $outdir/temp.vcf $outdir/$outfile.vcf
# use --real-ref-alleles

# Use VCFTOOLS to filter data
vcftools --vcf $outdir/${outfile}.vcf --chr 10 \
    --remove-filtered-all --remove-indels --max-alleles 2 \
    --recode --stdout | \
    gzip -c > $outdir/${outfile}_chr10.vcf.gz
    
vcftools --gzvcf $outdir/${outfile}_chr10.vcf.gz \
    --freq \
    --out $outdir/${outfile}_chr10


# BEAGLE. Phase genotypes
# create pedigree file to Beagle (if there are pairs ot trios)
# perform phasing and imputing by chromosome 

# $bin_dir/plink_1.9_mac --noweb --sheep \
#     --vcf  $outdir/$outfile.vcf \
#     --list-duplicate-vars \
#     --out $outdir/${outfile}_repe
    
for chrom in {1..27}
do
  vcftools --vcf $outdir/$outfile.vcf --chr ${chrom} \
    --remove-filtered-all --remove-indels --max-alleles 2 --recode --stdout | \
    gzip -c > $outdir/${outfile}_chr${chrom}.vcf.gz

  java -Xmx${HEAP}m -jar ~/libs/beagle/beagle_v4.0.jar \
    gt=$outdir/${outfile}_chr${chrom}.vcf.gz \
    nthreads=3 \
    out=$outdir/${outfile}_chr${chrom}_phase
  
  rm $outdir/${outfile}_chr${chrom}.vcf.gz
  #rm $outdir/${outfile}_chr${chrom}_phase.log
  
done

# concatenate phased and imputed VCF files
gunzip $outdir/${outfile}_chr*_phase.vcf.gz
cp $outdir/${outfile}_chr1_phase.vcf $outdir/${outfile}_phase.vcf

for chrom in {2..27}
do
 tail -n +11 -q $outdir/${outfile}_chr${chrom}_phase.vcf >> $outdir/${outfile}_phase.vcf
done

rm $outdir/${outfile}_chr*_phase.log
rm $outdir/${outfile}_chr*_phase.vcf
gzip $outdir/${outfile}_phase.vcf
gzip $outdir/${outfile}.vcf

# back to plink format. Use a2-allele from previous .bim file
awk '{print $2,$6}' $outdir/${outfile}.bim > ${outdir}/a2_reference.txt
$bin_dir/plink_1.9_mac --noweb --sheep --vcf  $outdir/${outfile}_phase.vcf.gz \
    --a2-allele ${outdir}/a2_reference.txt \
    --recode \
    --out $outdir/${outfile}_phased_kk



# Compare sequences

# brew info vcftools
# export PERL5LIB=/usr/local/lib/perl5/site_perl:${PERL5LIB}
awk '{printf "%s,", $1"_"$2}' $outdir/aDNA_samples
vcf-subset -c 0_FOS66,0_FOS82,0_FOS102,0_FOS59,0_FOS105 $outdir/${outfile}.vcf.gz | \
    bgzip -c > $outdir/${outfile}_aDNA.vcf.gz

awk '{printf "%s,", $1"_"$2}' $outdir/rDNA_samples
vcf-subset -c 0_OCP2,0_ME1076,0_OCP4,0_ME960,0_ME961 $outdir/${outfile}.vcf.gz | \
    bgzip -c > $outdir/${outfile}_rDNA.vcf.gz

tabix -p vcf $outdir/${outfile}_aDNA.vcf.gz
tabix -p vcf $outdir/${outfile}_rDNA.vcf.gz

vcf-compare  $outdir/${outfile}_aDNA.vcf.gz $outdir/${outfile}_rDNA.vcf.gz


## LD
# relationship matrix (not LD sensitive). Plot in R
$bin_dir/plink_1.9_mac --sheep \
    --bfile $outdir/$outfile \
    --make-rel \
    --out $outdir/${outfile}_LD
    
$bin_dir/plink_1.9_mac --sheep \
    --bfile $outdir/$outfile \
    --r2 --matrix \
    --out $outdir/${outfile}_LD

# LD between adjacent SNPs
$bin_dir/plink_1.9_mac --sheep \
    --bfile $outdir/$outfile \
    --r2 --ld-window 2 --ld-window-kb 1000 --ld-window-r2 0 \
    --out $outdir/${outfile}_LD2


## PCA analysis
# rscript
$bin_dir/plink_1.9_mac --sheep \
    --bfile $outdir/$outfile \
    --pca \
    --out $outdir/${outfile}_pca


## ADMIXTURE
$bin_dir/admixture $outdir/$outfile.bed 2
# $bin_dir/admixture -B $outdir/$outfile.bed 2 -j3
mv ${outfile}* $outdir/.

# choose optimal K
for k in {2..4}
do 
   $bin_dir/admixture --cv $outdir/$outfile.bed $k | tee $outdir/log_k$k
done
grep -h "CV" $outdir/log_k*
rm ${outfile}.* 

# consider LD
# thin markers according to the obs. sample correlation coeffs.
# output: a pruned subset of markers that are in approximate linkage 
# equilibrium with each other
# three parameters: a window size in variant count or kilobase, 
# a variant count to shift the window at the end of each step, and 
# a pairwise r2 threshold
$bin_dir/plink_1.9_mac --noweb --sheep --bfile $outdir/$outfile \
    --indep-pairwise 50 10 0.1 \
    --out $outdir/$outfile

$bin_dir/plink_1.9_mac --noweb --sheep --bfile $outdir/$outfile \
    --extract $outdir/$outfile.prune.in \
    --make-bed \
    --out $outdir/${outfile}_pruned

$bin_dir/admixture $outdir/${outfile}_pruned.bed 2
mv ${outfile}* $outdir/.



## ROH
$bin_dir/plink_1.9_mac --noweb --sheep \
    --vcf $outdir/${outfile}.vcf.gz \
    --homozyg \
    --homozyg-window-snp 10 \
    --homozyg-window-het 1 \
    --homozyg-kb 1000 \
    --homozyg-window-missing 2 \
    --out $outdir/${outfile}_roh

# with vcftools
for chrom in {1..27}
do
vcftools --gzvcf $outdir/${outfile}.vcf.gz \
    --chr $chrom \
    --LROH \
    --out $outdir/${outfile}_roh_$chrom
done
mv $outdir/${outfile}_roh_1.LROH $outdir/${outfile}_roh.LROH
for chrom in {2..27}
do
   tail -n +2 -q $outdir/${outfile}_roh_${chrom}.LROH >> $outdir/${outfile}_roh.LROH
done
rm $outdir/${outfile}_roh_*



## CNV



## SEQUENCES
# Another VCF example
zless data/VCF/all.fb.cf.gz
wc -l data/VCF/all.fb.vcf.gz
vcftools --gzvcf data/VCF/all.fb.vcf.gz \
    --chr 18 --chr 20 \
    --remove-filtered-all --remove-indels \
    --max-alleles 2 \
    --recode \
    --out data/VCF/onlysnps
    
gzip data/VCF/onlysnps.recode.vcf

$bin_dir/plink_1.9_mac --cow \
    --vcf data/VCF/onlysnps.recode.vcf.gz \
    --allow-extra-chr 0 --biallelic-only strict list \
    --vcf-filter \
    --recode \
    --out data/VCF/SNPS2plink

awk '{print $1}' data/VCF/SNPS2plink.map | sort | uniq -c
head data/VCF/SNPS2plink.ped | cut -d' ' -f1-18


echo ""
echo "\n  pipeline finished"
echo ""

### EOF