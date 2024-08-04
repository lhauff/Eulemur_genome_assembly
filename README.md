## Introduction 
Documentation of code used to conduct _de novo_ genome assembly for the endangered _Eulemur rufifrons_ from sequencing efforts conducted entirely with portable Nanopore sequencing in the host country of Madagascar. 

### Programs 
See programs.txt for full list of software programs and versions used for project

## Dorado (0.3.4) Rebasecalling 
Convert FAST5 files to pod5 (pod5 now available output for Dorado)
```shell
pod5 convert fast5 ./input/*.fast5  --output converted.pod5
pod5 merge e1.pod5 e2.pod5  --output merged.pod5
```
Download Dorado model
```shell
dorado download --model dna_r9.4.1_e8_sup@v3.6
```
Run Dorado
```shell
dorado basecaller -x "cuda:0" \
         dna_r9.4.1_e8_sup@v3.6 \
         pod5/eul_${SEQ}/ \
         --emit-fastq > fastqs/eul_${SEQ}_pod5_sup.fastq
```
Concatenate, rename, and gzip
```shell
cat *.fastq > eul_all.fq
mv eul_all.fq eul_all.fastq
gzip eul_all.fastq > eul_all.fastq.gz
```
## Nuclear Genome assembly 
### Trim adapters with PoreChop (0.2.4)
```shell
porechop -i fastqs/eul_all.fastq.gz -o fastqs/eul_trimmed.fastq -t 24
```
### Flye (2.9.1) assembly 
```shell
flye --nano-hq fastqs/eul_trimmed.fastq \
	--out-dir results/flye \
	--genome-size 3.0g \
	--asm-coverage 40 \
	-t 40
```

### Polish with Medaka (1.9.1)
Broke polishing into three distinct stages for memory purposes
```shell
mini_align -i fastqs/eul_all.fastq -r results/flye/eulemur_rufifrons_all_assembly.fasta -m \
    -p results/medaka/ -t 24

medaka consensus results/medaka/calls_to_draft.bam results/medaka/contigs_batch$BATCH.hdf \
    --model r941_min_sup_g507 --batch 200 --threads 2 \
    --regions $CONTIGS

medaka stitch results/medaka/*.hdf polished_eul_assembly.fasta
```
### Remove haplotigs with Purge_Haplotigs (1.1.2)
Map reads to assembly
```shell
minimap2 -t 20 -ax map-ont $FA $FQ | samtools sort -m 1G -o $OUT/aligned.bam -T tmp.ali
```
Generate read coverage histogram and manually determine cutoffs for low, medium, and high coverage
```shell
purge_haplotigs hist -b $BAM  -g $FA -t 16
```
Analyze coverage by contigs using cutoffs determined above
```shell
purge_haplotigs cov -i aligned.bam.gencov -l 5 -m 35 -h 75 \
	-o 75high/e_ruf_coverage_stats.csv \
	-j 80 -s 80
```
Run purging pipeline
```shell
purge_haplotigs purge -g $FA -c 75high/e_ruf_coverage_stats.csv \
	-t 16 \
	-o 75high/e_ruf_curated
```

### Remove contamination with Kraken2 (2.1.3)
Build standard database
```shell
./kraken2-build --standard --threads 24 --db database/
```
Run Kraken2
```shell
kraken/kraken2 --db kraken/standard_db/ --threads 24 --use-names --report results/kraken/eul_ruf_contam_kraken.report.txt /scratch/lrh85/nanopore/results/repeat/eulemur_rufifrons_all.fasta.masked --classified-out results/kraken/classified/eul_ruf_contam.kraken.classified.out.txt --unclassified-out results/kraken/unclassified/eul_ruf_contam.kraken.unclassified.out.txt 
```
Pull taxids from Kraken2 report and remove from genome
```shell
awk '/Bacteria/,/Viruses/' eul_ruf_all_kraken.report.txt | head -n -1 | cut -f 5 > bacteria_taxids.txt
grep -f bacteria_taxids.txt -Fw eul_ruf_all.kraken.classified.out.txt > bacteria_contigs.txt
cat bacteria_contigs.txt | cut -d " " -f1 > contigs_to_remove.txt
awk '(NR==FNR) { toRemove[$1]; next }
	/^>/ { p=1; for(h in toRemove) if ( h ~ $0) p=0 }
	p' contigs_to_remove.txt eulemur_rufifrons_all.fasta.masked > eulemur_rufirons_filtered_masked.fasta
```
Removing contigs smaller than 3k with seqkit (2.5.1)
```shell
cat eulemur_rufirons_filtered_masked.fasta | seqkit seq -m 3000 > eulruf_trim_filt_purged_3krm.fasta
```

### Masking repetitive regions with RepeatMasker (4.1.5)
```shell
FA=eulruf_trim_filt_purged_3krm.fasta
/cache/home/lrh85/RepeatMasker/RepeatMasker $FA -pa 54 -species Primate -xsmall
```

### Scaffold genome to reference with RagTag (2.1.0)
Scaffolding to _Lemur catta_ reference genome (mLemCat1.pri)
```shell
OUT=eulruf_scaffold_lcattaREF
REF=ref_genomes/lemur_catta/ncbi_dataset/data/GCF_020740605.2/GCF_020740605.2_mLemCat1.pri_genomic.fa
RES=results/repeat/eulemur_rufifrons_filtered_maksed.fasta

#initial scaffolding
ragtag.py scaffold \
	-o results/ragtag/$OUT \
	-w -r -t 20 \
	--mm2-params '-x map-ont' \
	$REF \
	$RES
```

### Additional assembly cleanup
Removing contigs which erroneously mapped to Y chromosome and unpolished mitochondrial genome
```shell
awk '(NR==FNR) { toRemove[$1]; next }
	/^>/ { p=1; for(h in toRemove) if ( h ~ $0) p=0 }
	p' rm.txt ragtag.scaffold.fasta > eulruf_filtered_purge_rmy.fasta
```
Adding polished mitochondrial genome produced with Flye meta mode (see below)
```shell
cat eulruf_filtered_purge_rmy.fasta polished_eulruf_mito_assembly.fasta > eulemur_rufifrons.fasta
```

### Evaulating genome completeness with QUAST (5.2.0)
Representative scripts to evaluate completeness of generated assembly by comparing to reference genome
```shell
quast $FA  \
	-r ref_genomes/lemur_catta/ncbi_dataset/data/RefSeq/GCF_020740605.2_mLemCat1.pri_genomic.fa \
        -g ref_genomes/lemur_catta/ncbi_dataset/data/RefSeq/genomic.gff \
        --nanopore fastqs/eul_trimmed.fasta \
        -o results/results/quast/eul_filt_purge_3krm_masked_scaff \
        --large --no-snps
```

### Evalutating genome completeness with BUSCO (5.7.0)
Representative scripts to evaluate completeness of generated assembly using vertebrate- and primate-specific lineage databases, respectively. 

```shell
busco -i $FASTA -m genome -l vertebrata_odb10 -c 24 -o $OUT_DIR -f --miniprot
busco -i $FASTA -m genome -l primates_odb10 -c 24 -o $OUT_DIR -f --miniprot
```

### Chaining to human genome for TOGA (1.1.5) annotation

```
Convert ref and target genome to 2bit
```shell
faToTwoBit genome.fa genome.2bit
```

Chain target (human) to query (_Eulemur rufifrons_ scaffolded to _Lemur_catta_)
```shell
/cache/home/lrh85/make_lastz_chains/make_chains.py \
        --pd results/chains/rufi_scaf_lcatta_chained_human -f \
        --cluster_executor slurm \
        --cluster_queue p_deenr_1 \
        human eulemur \
        2bit/hg38.2bit 2bit/eulruf_scaffold_lcattaRS.2bit
```

Rename chromosomes back
```shell
CHAIN=rufi_lcatta
RPATH=/home/lrh85/make_lastz_chains/standalone_scripts/

$RPATH/rename_chromosomes_back.py --rename_table_query eul_rufifrons_chrom_rename_table.tsv hg38.eul_rufifrons.final.chain
```

### Annotating with TOGA (1.1.5)
Convert target GTF to bed-12
```shell
gtfToGenePred ref_genomes/human/hg38.ncbiRefSeq.gtf ref_genomes/human/temp.genePred
genePredToBed ref_genomes/human/temp.genePred ref_genomes/human/hg38.bed
```
Annotating against set of human genes from <https://github.com/hillerlab/TOGA/tree/master/supply>
```shell
TOGA_PATH=/home/lrh85/TOGA
BED_PATH=ref_genomes/human/GRCh38/toga.transcripts.bed
ISO_PATH=ref_genomes/human/GRCh38/toga.isoforms.hg38.txt

$TOGA_PATH/toga.py results/chains/rufi_scaf_lcatta_chained_human/hg38.eulruf.renamed.final.chain $BED_PATH 2bit/hg38.2bit 2bit/eulruf_scaffold_lcattaRS.2bit \
	--pd results/toga/rufifrons --pn eul_rufifrons \
	-i $ISO_PATH --chn 20 --cb 8,16,32,64,128,256 \
	--cjn 100 \
	--nd $TOGA_PATH/nextflow_logs --nc $TOGA_PATH/nextflow_config_files 
```
Interpretting TOGA output by extracting Intact (I), Partially Intact (PI) or Uncertain Loss (UL) genes:
```shell
grep PROJECTION loss_summ_data.tsv | awk '{if ($3 == "I" || $3 == "PI" || $3 == "UL") print $2}' | sort > I_PI_UL.txt

sort -k4,4 query_annotation.bed -o query_annotation.sorted.bed

join -1 1 -2 4 I_PI_UL.txt query_annotation.sorted.bed -t $'\t'  > query_annotation.I_PI_UL.bed
```
To generate table of gene loss
```shell
grep GENE loss_summ_data.tsv | cut -f3 | sort | uniq -c
```

## Mitochondiral genome assembly

### Preparing reads
Map reads to _Lemur catta_ mitogenome
```shell
minimap2 -t 20 -ax map-ont ref_genomes/lemur_catta/ncbi_dataset/data/RefSeq/lcattaRS_mito.fasta fastqs/eul_all.fastq.gz > results/minimap/eulruf/eulruf_mito.sam
```
Subsample reads
```shell
seqkit sample -n 10000 eulruf_mito.fastq -o eulruf_mito_subsample.fastq
```
### Flye assembly in meta mode
```shell
flye --nano-hq fastqs/eulruf_mito_subsample.fastq \
        --out-dir results/flye/eulruf_mito_subsample \
        --genome-size 17k --meta \
        -t 40
```

### Polishing with Medaka
```shell
medaka_consensus -b 100 -i fastqs/eulruf_mito_subsample.fastq -d results/flye/eulruf_mito_subsample/eulruf_mito.fasta -o results/medaka/ -t 24 -m r941_min_sup_g507
```





