### ====================================================================
### PoolSeq Analysis
### Tribolium castaneum RNAi and 17-DMAG normal versus reduced-eye lines
### Author: Mohammed Errbii 2025
### ====================================================================

This pipline outlines the step to reproduce the bioinformatic analyses presented in the manuscript by [Sayed et al](https://doi.org/10.21203/rs.3.rs-6320655/v1).

⸻

# Genome Data Processing and Variant Calling Pipeline


#### Trimming Adapters and Low-Quality Bases

```bash
ls *_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz//' | \
parallel --jobs 8 --eta \
'java -jar trimmomatic-0.38.jar PE \
{=._R1_001.fastq.gz=}R1_001.fastq.gz {=._R1_001.fastq.gz=}R2_001.fastq.gz \
paired/{=._R1_001.fastq.gz=}_paired_R1.fastq.gz unpaired/{=._R1_001.fastq.gz=}_unpaired_R1.fastq.gz \
paired/{=._R1_001.fastq.gz=}_paired_R2.fastq.gz unpaired/{=._R1_001.fastq.gz=}_unpaired_R2.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads SLIDINGWINDOW:4:20 MINLEN:40'
```

⸻

#### Mapping Paired Reads to the Reference Genome

```bash
for infile in *_1.fq.gz; do
  base=$(basename -s _1.fq.gz ${infile})
  bwa mem -t 30 genome.fa ${infile} ${base}_2.fq.gz > bamfiles/${base}.sam
done
```

⸻

#### Clean SAM Files

```bash
ls bamfiles/*.sam | parallel --jobs 4 --eta \
'java -jar picard.jar CleanSam I={} O={.}C.sam'
```

⸻

#### Convert SAM to Sorted BAM

```bash
ls bamfiles/*C.sam | parallel --jobs 4 --eta \
'samtools view -Sb {} | samtools sort -o {.}S.bam'
```

⸻

#### Fix Mate Information

```bash
ls bamfiles/*CS.bam | parallel --jobs 4 --eta \
'java -jar picard.jar FixMateInformation I={} O={.}F.bam ADD_MATE_CIGAR=true ASSUME_SORTED=true'
```

⸻

#### Add Read Groups

Requires a tab-separated rgid.txt file with fields: input_bam, ID, PU, SM, LB

```bash
parallel -a rgid.txt --colsep '\t' --jobs 4 --eta \
'java -jar picard.jar AddOrReplaceReadGroups I={1} O={5}R.bam RGID={2} RGLB={4} RGPL=illumina RGPU={3} RGSM={4}'
```

⸻

#### Mark Duplicates

```bash
ls *R.bam | parallel --jobs 4 --eta \
'java -jar picard.jar MarkDuplicates I={} O={.}D.bam M={.}D.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 ASSUME_SORTED=true'
```

⸻

#### Index BAM Files

```bash
ls *D.bam | parallel --jobs 4 --eta 'samtools index {}'
```

⸻

# Pool-seq Analyses

#### Generate mpileup

```bash
samtools mpileup -B sample1.bam sample2.bam sample3.bam sample4.bam > all_samples.mpileup
```
#### Create Sync File

```bash
java -ea -Xmx12g -jar mpileup2sync.jar \
--input all_samples.mpileup \
--output ../popoolation2/all_samples.sync \
--fastq-type sanger \
--min-qual 20 \
--threads 20
```
#### Perform Fisher’s Exact Test (FET)

```bash
perl fisher-test.pl \
--input all_samples.sync \
--min-count 4 \
--min-coverage 20 \
--max-coverage 384,332,387,310 \
--min-covered-fraction 0 \
--output RNe_RRe_DNe_DRe.covFiltered.fet \
--suppress-noninformative \
--window-summary-method geometricmean
```
--max-coverage 2 % (equivalent to 384 for RNAi_normal-eye sample, 332 for RNAi reduced_eye, 387 for 17-DMAG normal_eye, and 310 for 17-DMAG reduced_eye)

⸻

# SNP Calling & Annotation

#### Call Variants Per Sample (HaplotypeCaller)

```bash
ls *D.bam | parallel --jobs 4 --eta \
'java -jar gatk.jar HaplotypeCaller -R genome.fa -I {} -O {.}.g.vcf.gz -ERC GVCF -ploidy 2'
```
#### Combine GVCFs

```bash
java -jar gatk.jar CombineGVCFs -R genome.fa \
--variant sample1.g.vcf.gz --variant sample2.g.vcf.gz \
--variant sample3.g.vcf.gz --variant sample4.g.vcf.gz\
-O combined.g.vcf.gz
```
#### Genotype GVCFs

```bash
java -jar gatk.jar GenotypeGVCFs \
-R genome.fa \
-V combined.g.vcf.gz \
-O genotyped.vcf.gz
```
#### Extract SNPs

```bash
java -jar gatk.jar SelectVariants \
-V genotyped.vcf.gz \
--select-type-to-include SNP \
-O snps.vcf.gz
```

#### Filter SNPs

```bash
vcftools \
--gzvcf snps.vcf.gz \
--recode-INFO-all \
--max-meanDP 231 \
--minDP 5 \
--min-meanDP 47.8 \
--max-missing 1 \
--out snps.filtered \
--recode
```

⸻

# SNP Annotation with SnpEff

```bash
java -jar snpEff.jar -v -stats snpeff_stats.html Tcas5.2 snps.filtered.vcf.gz | bgzip > snps.filtered.ann.vcf.gz
```