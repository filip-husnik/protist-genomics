# protist-genomics
Assembly and annotation of protist genomes

# Table of contents:

1. [Introduction](#introduction)
2. [Dependencies](#dependencies)
3. [Genome assembly](#genome-assembly)
4. [Transcriptome assembly](#transcriptome-assembly)
5. Contamination removal
6. [Repeat finding and soft masking](#repeat-finding-and-soft-masking)
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

Quality trimming (+provide adapters with --adapter_fasta adapters.fasta depending on your library prep)

```
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz --unpaired1 out.RS.fq.gz --unpaired2 out.RS.fq.gz --detect_adapter_for_pe --detect_adapter_for_pe --html LIBNAME_fastp.html --json LIBNAME_fastp.json --thread 12
```

Genome assembly (genomes <1000Mbp)
```
/opt/SPAdes-3.13.0-Linux/bin/spades.py -o default_spades --pe1-1 out.R1.fq.gz --pe1-2 out.R2.fq.gz --pe1-s out.RS.fq.gz --careful --threads 24
```
Contamination assessment for the SPAdes assembly with Blobtools1 (every protist genome is a metagenome)
```
cd default_spades
blastn -task megablast -query scaffolds.fasta -db /scratch/NCBI_NT/nt -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -culling_limit 5 -num_threads 16 -evalue 1e-25 -max_target_seqs 5 > scaffolds_vs_nt.blastn
/opt/bin/diamond blastx --query scaffolds.fasta --max-target-seqs 1 --sensitive --threads 14 --db /scratch/uniprot/uniprot_ref_proteomes.diamond.dmnd --evalue 1e-25 --outfmt 6 --out scaffolds.fasta.vs.uniprot_ref.mts1.1e25.out
/opt/blobtools/blobtools taxify -f scaffolds.fasta.vs.uniprot_ref.mts1.1e25.out -m /scratch/uniprot/uniprot_ref_proteomes.taxids -s 0 -t 2 
/opt/blobtools/blobtools create -i scaffolds.fasta -y spades -t scaffolds_vs_nt.blastn -t scaffolds.fasta.vs.uniprot_ref.mts1.1e25.taxified.out
/opt/blobtools/blobtools blobplot -i blobDB.json -r superkingdom
/opt/blobtools/blobtools blobplot -i blobDB.json
/opt/blobtools/blobtools blobplot -i blobDB.json -r order
/opt/blobtools/blobtools view -i blobDB.json --rank all
```

Contamination assessment for the SPAdes assembly with Blobtools2 (every protist genome is a metagenome)
```
blastn -db nt -query assembly.fasta -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 10 -max_hsps 1 -evalue 1e-25 -num_threads 16 -out blast.out

diamond blastx --query assembly.fasta --db /path/to/uniprot.db.with.taxids --outfmt 6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore --sensitive --max-target-seqs 1 --evalue 1e-25 --threads 16 > diamond.out

minimap2 -ax sr -t 16 assembly.fasta reads_1.fastq.gz reads_2.fastq.gz | samtools sort -@16 -O BAM -o assembly.reads.bam -

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
/opt/busco/scripts/run_BUSCO.py -i scaffolds.fasta -c 16 -m genome -l /opt/busco/databases/eukaryota_odb9/ --long --out BUSCO_genome
```
Genome assembly for genomes and metagenomes >1000Mbp

```
/opt/SPAdes-3.13.0-Linux/bin/metaspades.py -o default_metaspades --pe1-1 out.R1.fq.gz --pe1-2 out.R2.fq.gz --pe1-s out.RS.fq.gz --threads 24
```
```
/opt/megahit/megahit -t 24 -1 out.R1.fq.gz -2 out.R2.fq.gz -r out.RS.fq.gz
/opt/megahit/megahit_toolkit contig2fastg 141 ./intermediate_contigs/k141.contigs.fa > k141.fastg
```
Getting sequencing coverage for the Megahit assembly (Bowtie2+map2cov or BBmap)
```
/opt/bbmap/bbwrap.sh ref=final.contigs.fa in=out.R1.fq.gz in2=out.R2.fq.gz out=aln.sam.gz kfilter=22 subfilter=15 maxindel=80
/opt/bbmap/pileup.sh in=aln.sam.gz out=cov.txt

bowtie2-build final.contigs.fa final.contigs.fa
bowtie2 -p 16 -q --mm -x final.contigs.fa -1 out.R1.fq.gz -2 out.R2.fq.gz > aligned.sam
/opt/blobtools/blobtools map2cov -i final.contigs.fa -s aligned.sam
grep "#" -v aligned.sam.cov | cut -f 1,3 > aligned.coverage.autometa.tab

```
Extracting reads mapping to a subset of contigs (map all reads to the complete assembly and use bamfilter to filter based on contig ids)
```
bowtie2-build contigs.fasta contigs.fasta
bowtie2 -x contigs.fasta --reorder --mm -p 24 -1 forward.fastq.gz -2 reverse.fastq.gz | samtools view -Sb - | samtools sort -m 5G -o sorted.bam -f
/opt/blobtools/blobtools bamfilter -b contigs.fasta.mapped.bam -i ids.txt --sort --keep --threads 24
cat *.InIn.fq *.InUn.fq > subset_mapped.fastq
cat subset_mapped.fastq | paste - - - - - - - - | tee >(cut -f 1-4 | tr "\t" "\n" > subset_mapped.F.fastq) | cut -f 5-8 | tr "\t" "\n" > subset_mapped.R.fastq
```

Metagenome binning with Autometa

Running all steps in sequence (change -k to archaea if needed)
```
run_autometa.py --assembly scaffolds.fasta --processors 16 --length_cutoff 3000 --maketaxtable --ML_recruitment --output_dir autometa_default -db /media/Data1/Filip_autometa/autometa/databases -k bacteria
```

Running step by step
```
make_taxonomy_table.py --assembly scaffolds.fasta --processors 16 --length_cutoff 3000
run_autometa.py --assembly Bacteria.fasta --processors 16 --length_cutoff 3000 --taxonomy_table taxonomy.tab
ML_recruitment.py --contig_tab recursive_dbscan_output.tab --k_mer_matrix k-mer_matrix --out_table ML_recruitment_output.tab
```

Splitting bacterial contigs into genome bins (e.g. when interested in symbionts)

```
cluster_process.py --bin_table ML_recruitment_output.tab --column ML_expanded_clustering --fasta Bacteria.fasta --do_taxonomy --db_dir /scratch/Filip_storage/autometa/databases --output_dir cluster_process_output
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
ggplot( data, aes( x = bh_tsne_x, y = bh_tsne_y, col = order )) + geom_point( aes( alpha = 0.5, size = sqrt( data$length ) / 100 )) + guides( color = 'legend', size = 'none', alpha = 'none' ) + theme_classic() + xlab('BH-tSNE X') + ylab('BH-tSNE Y') + guides( color = guide_legend( title = 'Order' ))
```
Note: change coverage limits based on your coverage c( 0, 5000 )
```
ggplot( data, aes( x = cov, y = gc, col = ML_expanded_clustering )) + geom_point( aes( alpha = 0.5, size = sqrt( data$length ) / 100 )) + guides( color = 'legend', size = 'none', alpha = 'none' ) + theme_classic() + xlab('Coverage') + ylab('GC (%)') + guides( color = guide_legend( title = 'Cluster/bin' )) + scale_x_continuous( limits = c( 0, 5000 ))
```
Extracting a bacterial genome from a metagenome based on its assembly graph (often works for symbionts)

1. Open assembly_graph_with_scaffolds.gfa in Bandage [https://github.com/rrwick/Bandage]
2. Identify your genome of interest and select its graph (interconnected nodes)
3. Output->Save selected nodes sequences to FASTA -> symbiont_selected_nodes.fasta
4. Find and extract scaffolds corresponding to these nodes

```
#for SPAdes assemblies
grep "NODE_" symbiont_selected_nodes.fasta | cut -f 2 -d '_' > symbiont_selected_nodes.txt
grep -f symbiont_selected_nodes.txt assembly_graph_with_scaffolds.gfa | cut -f 2 | cut -f 1,2,3,4,5,6 -d '_' | sort | uniq > symbiont_selected_nodes_SPAdes.txt

#for Megahit k141 assemblies
grep ">" symbiont_selected_nodes.fasta | cut -f 2 -d'_' | sed s"/^/k141_/" | sed s"/+/ /"g | sed s"/-/ /"g |  sort | uniq > symbiont_selected_nodes_Megahit.txt

perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' symbiont_selected_nodes_SPAdes.txt scaffolds.fasta > symbiont_selected_nodes_to_scaffolds.fasta
```

Annotating a bacterial genome
```
/opt/prokka/bin/prokka --outdir PROKKA_annotation --compliant --gram neg --rfam scaffolds.fasta
```
Inferring orthologs, gene alignments, and gene trees for a set of proteomes (one proteome per fasta file)

```
#First, rename your proteins to have consistent headers
awk '/^>/{print ">SpeciesXYZ_" ++i; next}{print}' < proteome.faa > proteome_renamed.faa
```
```
/opt/OrthoFinder-2.3.3/orthofinder -f directory_with_proteomes -t 24 -M msa -A mafft -S diamond -T fasttree
```

Trimming single-gene alignments with trimal
```
for file in *.fasta; do trimal -in "$file" -out trimal_"$file" -automated1; done
```
Concatenating single-gene protein alignments into a multi-gene alignment
```
#remove protein identifiers from fasta files (only species names needed for concatenation)
for file in trimal_*; do sed "s/_[^_]*$//" "$file" > for_concatenation_"$file"; done
```
```
java -jar /opt/phyutility/phyutility.jar -concat -in trimal_* -out concatenated_alignment.fasta -aa
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
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz --unpaired1 out.RS.fq.gz --unpaired2 out.RS.fq.gz --html LIBNAME_fastp.html --json LIBNAME_fastp.json --thread 12
```

Single-cell RNA-Seq (SMART-Seq2 protocol with IS primers)

```
fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz --unpaired1 out.RS.fq.gz --unpaired2 out.RS.fq.gz --adapter_fasta /opt/Nextera_adapters_ISprimer.fa --html LIBNAME_fastp.html --json LIBNAME_fastp.json --thread 12
```
Transcriptome assembly with RNA-SPAdes
```
/opt/SPAdes-3.13.0-Linux/bin/rnaspades.py -o default_rnaspades --pe1-1 out.R1.fq.gz --pe1-2 out.R2.fq.gz --pe1-s out.RS.fq.gz --threads 16 --memory 150
```
Transcriptome assembly with Trinity
```
/opt/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --left out.R1.fq.gz --right out.R2.fq.gz --CPU 16 --full_cleanup --output default_trinity --max_memory 150G
```
Proteome prediction with TransDecoder
```
/opt/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t transcripts.fasta
/opt/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t transcripts.fasta
```
Completeness assessment with BUSCO
```
/opt/busco/scripts/run_BUSCO.py -i transcripts.fasta -c 16 -m transcriptome -l /opt/busco/databases/eukaryota_odb9/ --long --out BUSCO_transcr
/opt/busco/scripts/run_BUSCO.py -i transcripts.fasta -c 16 -m proteins -l /opt/busco/databases/eukaryota_odb9/ --long --out BUSCO_prot 
```

# Repeat finding and soft masking

Repeat masking will in most cases increase gene prediction accuracy. It reduces prediction of false genes in repetitive and low complexity genome regions. However, it's always useful to compare gene prediction between an unmasked and masked assembly to see any negative effects (for example, masking of real genes in gene-rich and repeat-poor genomes).

1) build a RepeatModeler database
```
/opt/RepeatModeler-open-1.0.11/BuildDatabase -name species_name -engine ncbi scaffolds_filtered.fasta
```
2) run RepeatModeler (if it fails with an out of bound error, rerun it until it works)
```
/opt/RepeatModeler-open-1.0.11/RepeatModeler -pa 24 -database species_name -engine ncbi >& run.out
```
3) run RepeatMasker to hard mask the reference (replace repeat sequences with Ns) 
```
/opt/RepeatMasker/RepeatMasker -pa 24 -html -gff -s -no_is -small -lib consensi.fa.classified scaffolds_filtered.fasta
```
4) run ProcessRepats to soft mask the reference instead (i.e. change repeat sequences to lowercase, recommended)
```
/opt/RepeatMasker/ProcessRepeats -maskSource scaffolds_filtered.fasta -xsmall scaffolds_filtered.fasta.cat.gz
```

# Gene prediction
(A) using the soft masked assembly generated above as a reference
1) run RNA-Seq mapping with STAR
```
mkdir genomeDir
/opt/STAR-2.7.0a/STAR --runMode genomeGenerate --runThreadN 24   --genomeDir genomeDir --genomeFastaFiles scaffolds_filtered_softmasked.fasta --limitGenomeGenerateRAM 67543940821
/opt/STAR-2.7.0a/STAR --runThreadN 24 --genomeDir genomeDir --readFilesIn R1_001.fastq.gz R2_001.fastq.gz --readFilesCommand zcat
```
2) Convert and sort sam to sorted.bam with Samtools (recent version, such as 1.9)
```
samtools view -bS Aligned.out.sam --threads 16 | samtools sort - -o Aligned.out.sorted.bam --threads 16
samtools index Aligned.out.sorted.bam
```

3) run Braker2 with hints from RNA-Seq data

```
braker.pl --species=species_name --genome=scaffolds_filtered_softmasked.fasta --bam=../STAR_RNA-Seq_mapping/Aligned.out.sorted.bam --cores=24 --softmasking

/opt/augustus_3.3_works/scripts/getAnnoFasta.pl --seqfile=genome.fa augustus.hints.gff
```
If GeneMark fails, renew its license in your home directory [https://github.com/Gaius-Augustus/BRAKER#perl-pipeline-dependencies].

(B) using the unmasked assembly as a reference
1) run RNA-Seq mapping with STAR
```
mkdir genomeDir
/opt/STAR-2.7.0a/STAR --runMode genomeGenerate --runThreadN 24   --genomeDir genomeDir --genomeFastaFiles scaffolds_filtered.fasta --limitGenomeGenerateRAM 67543940821
/opt/STAR-2.7.0a/STAR --runThreadN 24 --genomeDir genomeDir --readFilesIn R1_001.fastq.gz R2_001.fastq.gz --readFilesCommand zcat
```
2) Convert and sort sam to sorted.bam with Samtools (recent version, such as 1.9)
```
samtools view -bS Aligned.out.sam --threads 16 | samtools sort - -o Aligned.out.sorted.bam --threads 16
samtools index Aligned.out.sorted.bam
```
3) run Braker2 with hints from RNA-Seq data
```
braker.pl --species=species_name --genome=scaffolds_filtered.fasta --bam=Aligned.out.sorted.bam --cores=24
/opt/augustus_3.3_works/scripts/getAnnoFasta.pl --seqfile=genome.fa augustus.hints.gff
```
