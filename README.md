# SARS-CoV-2 detection from clinical samples using sequencing and bioinformatics method

The objective here is to provide a hand-on method to analyze sequencing data that come from COVID suspicious clinial samples. 
Questions we are trying to answer are:
1. If the clinical sample is SARS-CoV-2 positive/negative
2. The microbiome composition in the clinical sample

Pipeline input is Illumina short-reads sequencing data; Output includes number of reads mapped to SARS-CoV-2 reference, possible consensus genome generation if SARS-CoV-2 was detected, and pie chart of sample composition.

## Sequencing reads pre-process

This step filtered out adapters and low quality regions.
Software: Trimmomatic-0.39 (http://www.usadellab.org/cms/?page=trimmomatic); FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

```
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE $path-to-inputfolder/*_R1_001.fastq.gz $path-to-inputfolder/*_R2_001.fastq.gz \ 
$path-to-outputfolder/forward_paired.fq.gz $path-to-outputfolder/forward_unpaired.fq.gz \
$path-to-outputfolder/reverse_paired.fq.gz $path-to-outputfolder/reverse_unpaired.fq.gz \
ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-PE-2.fa:6:8:10:2:keepBothReads \
SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36
#Trim adapters, cut off quality < 15 regions
#reverse_paired.fq.gz forward_paired.fq.gz were kept for following analysis

fastqc -o outputfolder $path-to-outputfolder/forward_paired.fq.gz
fastqc -o outputfolder $path-to-outputfolder/reverse_paired.fq.gz
#confirm the reads quality

```

## Mapping sequencing reads to SARS-CoV-2 reference

Software: bwa (https://github.com/lh3/bwa); samtools (https://github.com/samtools/samtools)

```
bwa index $reference.fasta
bwa mem $reference.fasta forward_paired.fq.gz reverse_paired.fq.gz > $sample.sam

samtools sort $sample.sam -o $sample_sorted.bam
samtools index $sample_sorted.bam

```
The bam file can be visualized in IGV (http://software.broadinstitute.org/software/igv/) to visualize any mapped reads. 
To extract the mapped reads, use the following command.

```
samtools view -u -f 1 -F 12 $sample_sorted.bam > $sample_mapped.bam

samtools sort $sample_mapped.bam -o $sample_mapped_sorted.bam

bedtools bamtofastq -i $sample_mapped_sorted.bam -fq $sample_mapped.1.fastq -fq2 $sample_mapped.2.fastq

```
The mapped reads should be double checked by BLASTn to confirm they are SARS-CoV-2 reads. De novo assembly can be performed for the mapped reads.

## de novo assembly

Software: SPAdes-3.14.1 (https://github.com/ablab/spades); Quast (https://github.com/ablab/quast)

```
SPAdes-3.14.1-Linux/bin/spades.py --rna -1 $sample_mapped.1.fastq -2 $sample_mapped.2.fastq \
-o $assembly-output-folder
quast-master/quast-lg.py $assembly-output-folder/hard_filtered_transcripts.fasta -m 36 \
-o $assembly-quality-assess-folder

```

## Hypothesis-free taxonomic classification to reveal microbiome composition

Different from the approach above that has a suspected pathogen and hypothesis, the whole genome taxonomic classification can detect pathogens in a hypothesis-free way.

Software: Kraken2 (https://ccb.jhu.edu/software/kraken2/); Krona (https://github.com/marbl/Krona)

```
kraken2 --db minikraken_8GB_20200312 --threads 20 --gzip-compressed $sample.fq.gz > $sample.kraken
cut -f2,3 $sample.kraken > $sample.kraken_for_krona
ktImportTaxonomy -i -k -o krona_$sample.html $sample.kraken_for_krona

```
##Contact Me

Comments and discussions are welcome, email: stan@stjude.org
