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

### Features of Pipeline
Final version [1.0.8]()

* Filters variants by parameters QUAL>=30 and DP>=25
* Calculates measures of descriptive statistics for certain genes
* Retreves SRA samples automatically
* Detects copy number variants (CNVs) 
* Plots results for CNVs

## Requirements

1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)
6. [Control-FREEC - version 11.0](http://boevalab.com/FREEC/index.html#downloads)

#Note: Java v.1.8 is required to run Picard

### Reference genome

Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)
Fasta index file and Fasta dictionary file is required; steps to create them are described below.

##### Create FASTA index file
```
samtools faidx reference.fa
```
##### Create FASTA dictionary file
```
java -jar picard.jar CreateSequenceDictionary R=reference O=dictionary
```

### Help

To access help use the following command:

```{sh}
python3 ahcg_pipeline.py -h
```


## Virtual Box Setup
[Click here for instructions to setup](https://www.perkin.org.uk/posts/create-virtualbox-vm-from-the-command-line.html)

## Command to run pipeline
```
ahcg_pipeline_v1.0.1.py -c config_file.txt
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


