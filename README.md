# Applied Human Computational Genomics 2017
Author: Khushbu Patel
About: Second year Bioinformatics Master's student at Georgia Institute of Technology

## Background - Liquid Biopsy
## What is Liquid Biopsy?
Liquid Biopsies are non invasive tests performed to look for circulating cell free DNA from tumor cells in blood. A liquid biopsy has a potential to detect cancer at an early stage. 

### Liquid Biopsy Companies
* GRAIL
* Qiagen
* Biocept

### Alterations expected to be detected in Liquid Biopsy
* Somatic mutations-LOH events
* Somatic copy number alterations (CNAs) in tumor-normal exome data
* Epigenetic modifications 

### Examples of famous Variant calling pipelines to detect somatic and germline mutations
* MuTect2
* SomaticSnipper
* VarScan2
* MuSE
* LoFreq
* GDC DNA-seq 


### Sequence Coverage Analysis
[Illumina Sequence Coverage Calculator](https://support.illumina.com/downloads/sequencing_coverage_calculator.html)

Sequencing coverage describes the average number of reads that align to, or "cover," known reference bases. The next-generation sequencing (NGS) coverage level often determines whether variant discovery can be made with a certain degree of confidence at particular base positions.
Higher number of base coverage indicate that each base is covered by a greater number of aligned sequence reads, so base calls can be made with a higher degree of confidence.

Liquid biopsies require a much higher level of coverage than solid tumor biopsies
We need orders of magnitude higher than standard 900-1100X coverage

We need:-
To calculate coverage needed to accurately detect somatic mutations




# Applied Human Computational Genomics Pipeline
Variant calling pipeline for genomic data analysis

Aim: To detect somatic mutations in the pool of circulating cell free DNA by variant calling.

## Requirements

1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)

#Note: Java v.1.8 is required to run Picard

### Reference genome

Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

#### Create FASTA index file
```
samtools faidx reference.fa
```
#### Create FASTA dictionary file
```
java -jar picard.jar CreateSequenceDictionary R=reference O=dictionary
```

## Test data

Use the following protocol to download and prepare test dataset from NIST sample NA12878

```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```

## Help

To access help use the following command:

```{sh}
python3 ahcg_pipeline.py -h
```

## Installation
GATK version 3.4

## Virtual Box Setup
[Click here for instructions to setup](https://www.perkin.org.uk/posts/create-virtualbox-vm-from-the-command-line.html)

## Command to run pipeline
```
python ahcg_pipeline.py -i /path/to/FASTQ1 /path/to/FASTQ2 -o /path/to/output -p /path/to/picard.jar \
-g /path/to/GenomeAnalysisTK.jar -b /path/to/bowtie2 -w /path/to/genome.bt2 -r /path/to/genome.fa \
-d /path/to/genome.vcf -t /path/to/trimmomatic-0.36.jar -a /path/to/adapters.fa -o /path/to/output \
-s /path/to/guardant360.refGene_hg38.genes.bed
```

### Extract regions of interest from BAM file
```
samtools view input.bam "Chr10:18000-45500" > output.bam
```

### Determine per base sequence coverage of genome
```
samtools depth -r CHR1:1000-2000 input.bam > output.bed
```

### Obtain exon coordinates using UCSC genome browser
[Click here for the steps to Download the exon coordinates](https://github.com/kpatel427/ahcg2017_starterpipeline/blob/master/UCSC.pdf)

### Filter variants by quality, depth of coverage & type of mutation
[Filtering out variants to find "interesting / relevant variants" ](http://snpeff.sourceforge.net/SnpSift.html)

We choose the specific parameters here assuming-
* liquid biopsy samples
* we have exome sequencing data
* coverage is 200x to 250x 

```
java -jar GenomeAnalysisTK.jar -R FASTA -T SelectVariants -V input.vcf -o output.vcf -select "QUAL >= 30 && DP >= 25"

```

### Filtering CNVs 
Variants that are not naturally present since birth has to be distinguished from somatic variants. We need to incoorporate assessment of copy number variations and evaluation of contamination by normal cells in our pipelines.

[Control-FREEC (Control-FREE Copy number and allelic content caller)](http://boevalab.com/FREEC/tutorial.html) — a tool that annotates genotypes and discovers CNAs and LOH. 

Running first on test Data [(De Mattos-Arruda et al., 2015)](https://github.com/kpatel427/ahcg2017_starterpipeline/blob/master/Paper.pdf)

Design: Exome capture was performed using the Nextera Rapid Capture Exome kit; Illumina Hiseq 2000; WXS

1. Control sample - SRR2530741 (Germline)
```
/path/to/freec/freec -conf config_ctrl.txt
```
2. Test Data (50Kb window) - SRR2530742 (Tumor)
```
/PATH_TO_FREEC/freec -conf config_BL.txt
```

### Final pipeline version 
1.0.8
