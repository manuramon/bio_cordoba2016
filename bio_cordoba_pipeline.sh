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
plink_dir="bin"


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
awk 'BEGIN {FS="\t"; OFS="\t"}
     NR>10 {print 0,$2,$1,$3,$4}' $infile.txt  > $outdir/$outfile.temp.lgen

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
awk -F',' 'NR>1{print $12,$11}' data/SNPchimp_oarv4_0.csv > $outdir/oar_v4_bp_positions.txt
$plink_dir/plink_1.9_mac --noweb --sheep --lfile  $outdir/$outfile \
    --allow-extra-chr 0 --biallelic-only \
    --update-map $outdir/oar_v4_bp_positions.txt\
    --keep-allele-order \
    --make-bed\
    --out $outdir/${outfile}

$plink_dir/plink_1.9_mac --noweb --sheep --bfile  $outdir/$outfile \
    --allow-extra-chr 0 \
    --keep-allele-order \
    --make-bed\
    --out $outdir/${outfile}


# 1st check data 
$plink_dir/plink_1.9_mac --noweb --sheep --bfile  $outdir/$outfile \
    --allow-extra-chr 0 --biallelic-only \
    --freq --missing  --het small-sample --ibc --nonfounders\
    --out $outdir/${outfile}_stats
    
# do the same but spliting aDNA
awk '{if(substr($2,1,3)=="FOS") print $1,$2}' $outdir/${outfile}.fam > $outdir/aDNA_samples
$plink_dir/plink_1.9_mac --noweb --sheep --bfile  $outdir/$outfile \
    --allow-extra-chr 0 --biallelic-only \
    --keep $outdir/aDNA_samples \
    --freq --missing  --het small-sample  --ibc --nonfounders\
    --out $outdir/${outfile}_stats_aDNA

awk '{if(substr($2,1,3)!="FOS") print $1,$2}' $outdir/${outfile}.fam > $outdir/rDNA_samples
$plink_dir/plink_1.9_mac --noweb --sheep --bfile  $outdir/$outfile \
    --allow-extra-chr 0 --biallelic-only \
    --keep $outdir/rDNA_samples \
    --freq --missing  --het small-sample  --ibc --nonfounders\
    --out $outdir/${outfile}_stats_rDNA
    
rm ${outdir}/*.nosex
rm ${outdir}/*.log


# 2nd Quality Control (QC) 

# If we make QC using default thresholds with aDNA and rDNA data we will lose many SNPS
$plink_dir/plink_1.9_mac --sheep --bfile  $outdir/$outfile \
    --allow-extra-chr --biallelic-only \
    --geno 0.1  --mind 0.1 --maf 0.001 --missing \
    --keep-allele-order \
    --make-bed \
    --out $outdir/${outfile}_temp

# IMNP: check missing data reports!!!
$plink_dir/plink_1.9_mac --noweb --sheep --bfile  $outdir/$outfile \
    --allow-extra-chr --biallelic-only \
    --geno 0.5 --mind 1 --maf 0.0001 --missing \
    --keep-allele-order \
    --make-bed \
    --out $outdir/$outfile

rm $outdir/${outfile}_temp*

# 3rd recode as VCF file

$plink_dir/plink_1.9_mac --noweb --sheep --bfile  $outdir/$outfile \
    --allow-extra-chr 0 \
    --recode vcf \
    --out $outdir/$outfile

# there is a row given and advise that must be deleted from VCF file
sed '/may not be based on real reference genome/d' $outdir/$outfile.vcf > $outdir/temp.vcf
mv $outdir/temp.vcf $outdir/$outfile.vcf


# Use VCFTOOLS to filter data
vcftools --vcf $outdir/${outfile}.vcf --not-chr 0 \
    --remove-filtered-all --remove-indels --max-alleles 2 --recode --stdout | \
    gzip -c > $outdir/${outfile}_NOT0.vcf.gz


# BEAGLE. Phase genotypes
# create pedigree file to Beagle (if there are pairs ot trios)
# awk 'BEGIN{FS=";"; OFS=" "}; NR>1{print 0,$1,$2,$3}' ../ped_25082015 > ped2bgl.txt

# get reference panel from:
# wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/sheep_9940/VCF/*.gz
#vcftools --gzvcf vcf_chr_1.vcf.gz --remove-filtered-all --remove-indels --max-alleles 2 \
#    --recode \
##   --recode-INFO-all \ # to keep info
#    --stdout | gzip -c > vcf_chr_1.recode.vcf.gz


# perform phasing and imputing by chromosome 
for chrom in {1..27}
do
  vcftools --vcf $outdir/$outfile.vcf --chr ${chrom} \
    --remove-filtered-all --remove-indels --max-alleles 2 --recode --stdout | \
    gzip -c > $outdir/${outfile}_chr${chrom}.vcf.gz

  java -Xmx${HEAP}m -jar ~/libs/beagle/beagle_v4.1.jar \
    gt=$outdir/${outfile}_chr${chrom}.vcf.gz \
    nthreads=3 \
    ne=200\
    out=$outdir/${outfile}_chr${chrom}_phase
  
  rm $outdir/${outfile}_chr${chrom}.vcf.gz
  #rm $outdir/${outfile}_chr${chrom}_phase.log
  
done

# concatenate phased and imputed VCF files
gunzip $outdir/${outfile}_chr*_phase.vcf.gz
cp $outdir/${outfile}_chr1_phase $outdir/${outfile}_phase.vcf

for chrom in {2..27}
do
 tail -n +11 -q $outdir/${outfile}_chr${chrom}_phase.vcf >> $outdir/${outfile}_phase.vcf
done

rm $outdir/${outfile}_chr*_phase.log
rm $outdir/${outfile}_chr*_phase.vcf
gzip $outdir/${outfile}_phase.vcf
gzip $outdir/${outfile}.vcf

# back to plink format. Use a2-allele from previous .bim file
awk '{print $2,$6}' $outdir/${outfile}.bim > $outdir/a2_reference.txt
$plink_dir/plink_1.9_mac --noweb --sheep --vcf  $outdir/${outfile}_phase.vcf.gz \
    --recode --a2-allele a2_reference.txt --make-bed \
    --out $outdir/${outfile}_phased

 

echo ""
echo "\n  Manchega_pipeline_V0.1.sf finished"
echo ""

### EOF