# protist-genomics
Assembly and annotation of protist genomes

# Table of contents:

1. [Introduction](#introduction)
2. [Dependencies](#dependencies)
3. Genome assembly
4. Transcriptome assembly
5. Contamination removal
6. [Repeat masking](#repeat-finding-and-masking)
7. [Gene prediction](#gene-prediction)


# Introduction

This page describes my pipeline for assembling and annotating protist genomes. My goal is to have a quick pipeline that works with protists from all eukaryotic supergroups with no need to change default settings. Required input data: genome and RNA-seq reads from the same species.

Any suggestions and critiques are welcome!

# Dependencies

For detecting and masking repeats, I use two programs: [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/)
and [RepeatMasker](http://www.repeatmasker.org/). However, these programs can be sometimes difficult to install and rely on numerous dependencies. If you know of an alternative way how to do this, please let me know!

For gene prediction, I use [BRAKER2](https://github.com/Gaius-Augustus/BRAKER). Braker2 relies on Augustus and GeneMarkET and can take advantage of RNA-Seq data. For mapping RNA-Seq reads, I use the [STAR aligner](https://github.com/alexdobin/STAR).

# Repeat finding and masking
1) build a RepeatModeler database
```
/opt/RepeatModeler-open-1.0.11/BuildDatabase -name species_name -engine ncbi scaffolds_filtered.fasta
```
2) run RepeatModeler (if it fails with an out of bound error, rerun it until it works)
```
/opt/RepeatModeler-open-1.0.11/RepeatModeler -pa 24 -database species_name -engine ncbi >& run.out
```
3) run RepeatMasker
```
/opt/RepeatMasker/RepeatMasker -pa 10 -lib consensi.fa.classified scaffolds_filtered.fasta
```
# Gene prediction
1) run RNA-Seq mapping with STAR
```
/opt/STAR-2.7.0a/STAR --runMode genomeGenerate --runThreadN 24   --genomeDir genomeDir/ --genomeFastaFiles scaffolds_filtered.fasta --limitGenomeGenerateRAM 67543940821

/opt/STAR-2.7.0a/STAR --runThreadN 24 --genomeDir genomeDir/ --readFilesIn R1_001.fastq.gz R2_001.fastq.gz --readFilesCommand zcat
```
2) Convert and sort sam to sorted.bam with Samtools

3) run braker with hints from RNA-Seq data
```
braker.pl --species=species_name --genome=scaffolds_filtered.fasta --bam=../STAR_RNA-Seq_mapping/Aligned.out.sorted.bam
```
