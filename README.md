# sars-cov-2 detection using sequencing and bioinformatics

The purpose here is to provide a hand-on method to analyze sequencing data that come from COVID suspicious samples.

## Sequencing reads pre-process

Software: Trimmomatic-0.39 (http://www.usadellab.org/cms/?page=trimmomatic); FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

```
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE $path-to-inputfolder/*_R1_001.fastq.gz $path-to-inputfolder/*_R2_001.fastq.gz \ 
$path-to-outputfolder/forward_paired.fq.gz $path-to-outputfolder/forward_unpaired.fq.gz \
$path-to-outputfolder/reverse_paired.fq.gz $path-to-outputfolder/reverse_unpaired.fq.gz \
ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:6:8:10:2:keepBothReads \
SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36
#Trim adapters, cut off quality < 15 regions
#reverse_paired.fq.gz forward_paired.fq.gz will be kept for following analysis

fastqc -o outputfolder $path-to-outputfolder/forward_paired.fq.gz
fastqc -o outputfolder $path-to-outputfolder/reverse_paired.fq.gz
#check the quality

```

## Mapping the sars-cov-2

Software: bwa (https://github.com/lh3/bwa); samtools (https://github.com/samtools/samtools)

```
bwa index $reference.fasta
bwa mem $reference.fasta forward_paired.fq.gz reverse_paired.fq.gz > $sample.sam
samtools sort $sample.sam -o $sample_sorted.bam
samtools index $sample_sorted.bam

```

## Metagenomic analysis

Software:

```
kraken2 --db minikraken_8GB_20200312 --threads 20 --gzip-compressed $sample.fq.gz > $sample.kraken
cut -f2,3 $sample.kraken > $sample.kraken_for_krona
ktImportTaxonomy -i -k -o krona_$sample.html $sample.kraken_for_krona

```
