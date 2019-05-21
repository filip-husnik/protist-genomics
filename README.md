# protist-genomics
Assembly and annotation of protist genomes

# Table of contents:

1. [Introduction](#introduction)
2. [Dependencies](#dependencies)
3. [Genome assembly](#genome-assembly)
4. [Transcriptome assembly](#transcriptome-assembly)
5. Contamination removal
6. [Repeat masking](#repeat-finding-and-masking)
7. [Gene prediction](#gene-prediction)


# Introduction

This page describes my pipeline for assembling and annotating protist genomes. My goal is to have a quick pipeline that works with protists from all eukaryotic supergroups with no need to change default settings. Required input data: genome and RNA-seq reads from the same species.

Any suggestions and critiques are welcome!

# Dependencies

- [FastP](https://github.com/OpenGene/fastp)
- [SPAdes](http://cab.spbu.ru/software/spades)
- [Megahit](https://github.com/voutcn/megahit)
- [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki)
- [Blobtools](https://blobtools.readme.io/docs)
- [Autometa](https://bitbucket.org/jason_c_kwan/autometa)
- [BUSCO](https://busco.ezlab.org/)
- [STAR](https://github.com/alexdobin/STAR)
- [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/)
- [RepeatMasker](http://www.repeatmasker.org/)
- [BRAKER2](https://github.com/Gaius-Augustus/BRAKER)

RepeatModeler and RepeatMaskers can be sometimes difficult to install and rely on numerous dependencies. If you know of an alternative way how to do this, please let me know!

# Genome assembly

Quality trimming

```
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz --unpaired1 out.RS.fq.gz --unpaired2 out.RS.fq.gz
```
Genome assembly (genomes <1000Mbp)
```
/opt/SPAdes-3.13.0-Linux/bin/spades.py -o default_spades --pe1-1 out.R1.fq.gz --pe1-2 out.R2.fq.gz --pe1-s out.RS.fq.gz --careful --threads 24
```
Contamination assessment for the SPAdes assembly
```
cd default_spades
blastn -task megablast -query scaffolds.fasta -db /scratch/NCBI_NT/nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -culling_limit 5 -num_threads 16 -evalue 1e-25 -max_target_seqs 5 > scaffolds_vs_nt.blastn
/opt/bin/diamond blastx --query scaffolds.fasta --max-target-seqs 1 --sensitive --threads 14 --db /scratch/uniprot/uniprot_ref_proteomes.diamond.dmnd --evalue 1e-25 --outfmt 6 --out scaffolds.fasta.vs.uniprot_ref.mts1.1e25.out
/opt/blobtools/blobtools taxify -f scaffolds.fasta.vs.uniprot_ref.mts1.1e25.out -m /scratch/uniprot/uniprot_ref_proteomes.taxids -s 0 -t 2 
/opt/blobtools/blobtools create -i scaffolds.fasta -y spades -t scaffolds_vs_nt.blastn -t scaffolds.fasta.vs.uniprot_ref.mts1.1e25.taxified.out
/opt/blobtools/blobtools blobplot -i blobDB.json -r superkingdom
/opt/blobtools/blobtools blobplot -i blobDB.json
/opt/blobtools/blobtools view -i blobDB.json --rank all
```
SSU contamination assessment (for reads <150bp, e.g. HiSeq)
```
/opt/phyloFlash-pf3.3b1/phyloFlash.pl -lib LIBNAME -CPUS 16 -read1 out.R1.fq.gz -read2 out.R2.fq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -everything -id 55 -trusted scaffolds.fasta
```
SSU contamination assessment (for reads of 250bp, e.g. MiSeq)
```
/opt/phyloFlash-pf3.3b1/phyloFlash.pl -lib LIBNAME -CPUS 16 -read1 out.R1.fq.gz -read2 out.R2.fq.gz -dbhome /Data/filip/phyloFlash_DB/128/ -everything -id 55 -trusted scaffolds.fasta -readlength 250
```
Completeness assessment
```
/opt/busco/scripts/run_BUSCO.py -i scaffolds.fasta -c 16 -m genome -l /opt/busco/databases/eukaryota_odb9/ --long
```
Genome assembly for genomes and metagenomes >1000Mbp

```
/opt/SPAdes-3.13.0-Linux/bin/metaspades.py -o default_metaspades --pe1-1 out.R1.fq.gz --pe1-2 out.R2.fq.gz --pe1-s out.RS.fq.gz --threads 24
```
```
megahit -t 24 -1 out.R1.fq.gz -2 out.R2.fq.gz -r out.RS.fq.gz
```
Metagenome binning with Autometa

```
run_autometa.py --assembly scaffolds.fasta --processors 16 --length_cutoff 500 --maketaxtable --ML_recruitment
```

Splitting bacterial contigs into genome bins (e.g. when interested in symbionts)

```
cluster_process.py --bin_table ML_recruitment_output.tab --column ML_expanded_clustering --fasta Bacteria.fasta --do_taxonomy --db_dir /home/filip/autometa/databases --output_dir cluster_process_output
```
Visualizing the bacterial bins

```
R
library(ggplot2)
data = read.table('ML_recruitment_output.tab', header=TRUE, sep='\t')

```
```
ggplot( data, aes( x = bh_tsne_x, y = bh_tsne_y, col = ML_expanded_clustering )) + geom_point( aes( alpha = 0.5, size = sqrt( data$length ) / 100 )) + guides( color = 'legend', size = 'none', alpha = 'none' ) + theme_classic() + xlab('BH-tSNE X') + ylab('BH-tSNE Y') + guides( color = guide_legend( title = 'Cluster/bin' ))
```
```
ggplot( data, aes( x = bh_tsne_x, y = bh_tsne_y, col = phylum )) + geom_point( aes( alpha = 0.5, size = sqrt( data$length ) / 100 )) + guides( color = 'legend', size = 'none', alpha = 'none' ) + theme_classic() + xlab('BH-tSNE X') + ylab('BH-tSNE Y') + guides( color = guide_legend( title = 'Phylum' ))
```
```
ggplot( data, aes( x = cov, y = gc, col = ML_expanded_clustering )) + geom_point( aes( alpha = 0.5, size = sqrt( data$length ) / 100 )) + guides( color = 'legend', size = 'none', alpha = 'none' ) + theme_classic() + xlab('Coverage') + ylab('GC (%)') + guides( color = guide_legend( title = 'Cluster/bin' )) + scale_x_continuous( limits = c( 200, 250 ))
```
Extracting a bacterial genome from a metagenome based on its assembly graph

1. Open assembly_graph_with_scaffolds.gfa in Bandage [https://github.com/rrwick/Bandage]
2. Identify your genome of interest and select its graph (interconnected nodes)
3. Output->Save selected nodes sequences to FASTA -> symbiont_selected_nodes.fasta
4. Find and extract scaffolds corresponding to these nodes
```
grep "NODE_" symbiont_selected_nodes.fasta | cut -f 2 -d '_' > symbiont_selected_nodes.txt
grep -f symbiont_selected_nodes.txt assembly_graph_with_scaffolds.gfa | cut -f 2 | sort | uniq > symbiont_selected_nodes_to_scaffolds.fasta
```

Annotating a bacterial genome
```
/opt/prokka/bin/prokka --outdir PROKKA_annotation --compliant --gram neg --rfam scaffolds.fasta
```

Dot plot alignment of two closely related bacterial genomes
```
nucmer reference.fasta query.fasta
mummerplot -s large --layout -t png out.delta
```
Dot plot alignment of two bacterial genomes
```
promer reference.fasta query.fasta
mummerplot -s large --layout -t png out.delta
```


# Transcriptome assembly

Quality trimming

Standard RNA extraction and library prep (e.g. from culture)

```
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz --unpaired1 out.RS.fq.gz --unpaired2 out.RS.fq.gz
```

Single-cell RNA-Seq (SMART-Seq2 protocol with IS primers)

```
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz --unpaired1 out.RS.fq.gz --unpaired2 out.RS.fq.gz --adapter_fasta /opt/Nextera_adapters_ISprimer.fa
```
Transcriptome assembly with RNA-SPAdes
```
/opt/SPAdes-3.13.0-Linux/bin/rnaspades.py -o default_rnaspades --pe1-1 out.R1.fq.gz --pe1-2 out.R2.fq.gz --pe1-s out.RS.fq.gz --threads 16 --memory 150
```
Transcriptome assembly with Trinity
```
/opt/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --left out.R1.fq.gz --right out.R2.fq.gz --CPU 16 --full_cleanup --output default_trinity --max_memory 150G
```

# Repeat finding and masking

Repeat masking will in most cases increase gene prediction accuracy. It reduces prediction of false genes in repetitive and low complexity genome regions. However, it's always useful to compare gene prediction between an unmasked and masked assembly to see any negative effects (for example, masking of real genes in gene-rich and repeat-poor genomes).

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

3) run Braker2 with hints from RNA-Seq data
```
braker.pl --species=species_name --genome=scaffolds_filtered.fasta --bam=../STAR_RNA-Seq_mapping/Aligned.out.sorted.bam
```
