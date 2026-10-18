# Linux Command-Line Cheatsheet for Bioinformatics

A practical reference for everyday bioinformatics work on Linux: shell basics, text processing, sequence files (FASTA/FASTQ), alignments (SAM/BAM), variants (VCF), intervals (BED/GTF), common tools and HPC.

> Conventions: `ref.fa` = reference genome, `r1.fq.gz`/`r2.fq.gz` = paired-end reads, `-t`/`-p`/`-@` = threads.

---

## Table of Contents

1. [File Formats at a Glance](#1-file-formats-at-a-glance)
2. [Shell Basics](#2-shell-basics)
3. [Text Processing Power Tools](#3-text-processing-power-tools)
4. [Compression & Archives](#4-compression--archives)
5. [FASTA One-Liners](#5-fasta-one-liners)
6. [FASTQ One-Liners](#6-fastq-one-liners)
7. [seqkit](#7-seqkit)
8. [Getting Data (NCBI, SRA, Ensembl)](#8-getting-data-ncbi-sra-ensembl)
9. [Quality Control & Trimming](#9-quality-control--trimming)
10. [Read Alignment](#10-read-alignment)
11. [samtools (SAM/BAM/CRAM)](#11-samtools-sambamcram)
12. [Variant Calling & bcftools (VCF)](#12-variant-calling--bcftools-vcf)
13. [bedtools (Intervals)](#13-bedtools-intervals)
14. [GTF/GFF Annotation](#14-gtfgff-annotation)
15. [RNA-seq Quantification](#15-rna-seq-quantification)
16. [BLAST](#16-blast)
17. [Environments: conda / mamba](#17-environments-conda--mamba)
18. [Loops & Parallelism](#18-loops--parallelism)
19. [HPC: SLURM, tmux, Transfers](#19-hpc-slurm-tmux-transfers)
20. [Data Integrity & Good Practice](#20-data-integrity--good-practice)

---

## 1. File Formats at a Glance

| Format | Content | Key facts |
|---|---|---|
| FASTA (`.fa`, `.fasta`, `.fna`, `.faa`) | Sequences | Header line starts with `>` |
| FASTQ (`.fq`, `.fastq`) | Reads + qualities | 4 lines per read: `@id`, seq, `+`, quality |
| SAM / BAM / CRAM | Alignments | SAM = text, BAM = binary, CRAM = reference-compressed |
| VCF / BCF | Variants | Header lines start with `#`; BCF = binary |
| BED | Intervals | **0-based**, half-open; `chrom start end` |
| GTF / GFF3 | Annotation | **1-based**, closed; 9 columns |
| BigWig / bedGraph | Coverage tracks | bedGraph = text, BigWig = indexed binary |

---

## 2. Shell Basics

```bash
pwd                         # where am I
ls -lh                      # list with human-readable sizes
ls -lhS                     # sort by size
cd ~/projects/rnaseq        # change directory
mkdir -p data/{raw,trimmed,aligned}   # create nested dirs in one go
cp -r src/ dest/            # copy recursively
mv old.txt new.txt          # move / rename
rm -i file                  # remove (ask first)
ln -s /path/to/ref.fa .     # symlink large files instead of copying

less -S file.tsv            # view without line wrapping (great for wide tables)
head -n 20 file; tail -n 20 file
tail -f log.txt             # follow a growing log
wc -l file                  # count lines

find . -name "*.fq.gz"      # find files
find . -name "*.bam" -size +1G   # big BAMs
du -sh *                    # size of each item
df -h                       # disk space
chmod +x script.sh          # make executable
which samtools              # where is a tool
history | grep bwa          # find earlier commands
```

**Redirection & pipes**

```bash
cmd > out.txt               # overwrite
cmd >> out.txt              # append
cmd 2> err.log              # stderr to file
cmd > out.txt 2>&1          # stdout + stderr
cmd1 | cmd2                 # pipe
cmd | tee out.txt           # show AND save
```

---

## 3. Text Processing Power Tools

```bash
# grep
grep -c "pattern" file             # count matching lines
grep -v "^#" file.vcf              # drop header lines
grep -w "BRCA1" genes.tsv          # whole word
grep -i "kinase" annot.txt         # case-insensitive
grep -E "chr(1|2|X)\b" file        # extended regex
grep -f ids.txt file               # match a list of IDs

# cut / sort / uniq
cut -f1,3 file.tsv                 # columns 1 and 3 (tab)
cut -d',' -f2 file.csv             # comma-delimited
sort -k2,2n file                   # numeric sort on column 2
sort -k1,1 -k2,2n file.bed         # sort BED by chrom, then start
sort -u file                       # unique sorted
sort file | uniq -c | sort -nr     # frequency table

# awk
awk '{print $1, $3}' file                      # columns
awk -F'\t' '$5 > 30' file                      # filter rows
awk '{s+=$2} END {print s}' file               # column sum
awk '{s+=$2} END {print s/NR}' file            # column mean
awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2}' pos.txt   # 1-based pos → BED

# sed
sed 's/chr//' file                 # remove "chr" prefix (first per line)
sed 's/old/new/g' file             # replace all
sed -n '10,20p' file               # print lines 10–20
sed -i 's/\r$//' file              # fix Windows line endings in place

# other
paste a.txt b.txt                  # merge side by side
join -t $'\t' a.sorted b.sorted    # join on column 1 (both sorted)
tr 'acgt' 'ACGT'                   # translate characters
column -t file.tsv | less -S       # pretty-print table
comm -12 <(sort a) <(sort b)       # lines common to both files
diff a.txt b.txt
```

---

## 4. Compression & Archives

```bash
gzip file / gunzip file.gz
zcat file.gz | head                # read without decompressing
zless file.gz
zgrep "pattern" file.gz
pigz -p 8 file                     # parallel gzip (much faster)

bgzip file.vcf                     # block gzip (needed for indexing)
tabix -p vcf file.vcf.gz           # index VCF
tabix -p bed file.bed.gz           # index BED
tabix file.vcf.gz chr1:1000-2000   # query region

tar -czvf archive.tar.gz folder/   # create
tar -xzvf archive.tar.gz           # extract
tar -tzvf archive.tar.gz           # list contents
unzip file.zip
```

---

## 5. FASTA One-Liners

```bash
grep -c "^>" seqs.fa                         # number of sequences
grep "^>" seqs.fa | head                     # headers
grep "^>" seqs.fa | sed 's/>//' > ids.txt    # list IDs

# total bases
grep -v "^>" seqs.fa | tr -d '\n' | wc -c

# length of each sequence
awk '/^>/{if(id)print id"\t"len; id=substr($1,2); len=0; next}
     {len+=length($0)} END{print id"\t"len}' seqs.fa

# linearize (one line per sequence)
awk '/^>/{printf("%s%s\t",(N>0?"\n":""),$0);N++;next}
     {printf("%s",$0)} END{printf("\n")}' seqs.fa

# simplify headers (keep first word)
awk '/^>/{print $1; next}{print}' seqs.fa > clean.fa

# index & extract a region (samtools)
samtools faidx ref.fa
samtools faidx ref.fa chr1:10000-10500
cut -f1,2 ref.fa.fai > genome.sizes          # chrom sizes for bedtools
```

---

## 6. FASTQ One-Liners

```bash
# number of reads
echo $(( $(zcat r1.fq.gz | wc -l) / 4 ))

# first 3 reads
zcat r1.fq.gz | head -n 12

# FASTQ → FASTA
zcat r1.fq.gz | awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' > r1.fa

# read length distribution
zcat r1.fq.gz | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c

# average read length
zcat r1.fq.gz | awk 'NR%4==2{s+=length($0); n++} END{print s/n}'

# keep first 100k reads
zcat r1.fq.gz | head -n 400000 | gzip > r1.sub.fq.gz

# check pairs have the same read count
for f in r1.fq.gz r2.fq.gz; do echo "$f $(( $(zcat $f | wc -l)/4 ))"; done
```

---

## 7. seqkit

Fast, multithreaded FASTA/FASTQ toolkit — often replaces the one-liners above.

```bash
seqkit stats *.fq.gz                     # counts, lengths, N50 (-a for all stats)
seqkit seq -m 1000 seqs.fa               # keep seqs ≥ 1000 bp
seqkit seq -r -p seqs.fa                 # reverse complement
seqkit fx2tab -n -l -g seqs.fa           # name, length, GC%
seqkit grep -f ids.txt seqs.fa           # extract by ID list
seqkit grep -s -p GAATTC seqs.fa         # search sequence motif
seqkit sample -p 0.1 -s 42 r1.fq.gz      # subsample 10% (seed 42)
seqkit head -n 1000 r1.fq.gz             # first N records
seqkit rmdup -s seqs.fa                  # remove duplicate sequences
seqkit sort -l -r seqs.fa                # sort by length, longest first
seqkit split -p 4 seqs.fa                # split into 4 parts
seqkit fq2fa r1.fq.gz > r1.fa            # FASTQ → FASTA
seqkit translate cds.fa                  # DNA → protein
```

---

## 8. Getting Data (NCBI, SRA, Ensembl)

```bash
# plain downloads
wget -c https://example.org/file.fa.gz      # -c resumes
curl -L -O https://example.org/file.fa.gz

# SRA (sra-tools)
prefetch SRR1234567
fasterq-dump --split-files -e 8 SRR1234567  # -> SRR1234567_1.fastq, _2.fastq
pigz -p 8 SRR1234567_*.fastq

# NCBI datasets CLI
datasets download genome accession GCF_000001405.40 --include genome,gtf
unzip ncbi_dataset.zip
datasets summary genome taxon "Escherichia coli" --reference

# Entrez Direct
esearch -db nucleotide -query "NC_045512" | efetch -format fasta > sarscov2.fa
esearch -db protein -query "BRCA1 human" | efetch -format fasta > brca1.faa

# Ensembl (example: human GRCh38)
wget https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

---

## 9. Quality Control & Trimming

```bash
mkdir -p qc
fastqc -t 8 -o qc/ *.fq.gz                   # per-sample QC reports
multiqc qc/ -o qc/                           # aggregate into one report

# fastp (adapter + quality trimming, paired-end)
fastp -i r1.fq.gz -I r2.fq.gz \
      -o r1.trim.fq.gz -O r2.trim.fq.gz \
      -h fastp.html -j fastp.json -w 8

# Trimmomatic (paired-end)
trimmomatic PE -threads 8 r1.fq.gz r2.fq.gz \
  r1.P.fq.gz r1.U.fq.gz r2.P.fq.gz r2.U.fq.gz \
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:36

# cutadapt
cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -q 20 -m 30 \
  -o r1.trim.fq.gz -p r2.trim.fq.gz r1.fq.gz r2.fq.gz
```

---

## 10. Read Alignment

```bash
# BWA-MEM (DNA short reads)
bwa index ref.fa
bwa mem -t 8 -R "@RG\tID:s1\tSM:s1\tPL:ILLUMINA" ref.fa r1.fq.gz r2.fq.gz \
  | samtools sort -@ 4 -o s1.bam
samtools index s1.bam

# Bowtie2
bowtie2-build ref.fa ref_idx
bowtie2 -p 8 -x ref_idx -1 r1.fq.gz -2 r2.fq.gz | samtools sort -o s1.bam

# minimap2 (long reads / assemblies / short reads)
minimap2 -ax map-ont ref.fa ont.fq.gz | samtools sort -o ont.bam   # Nanopore
minimap2 -ax map-hifi ref.fa hifi.fq.gz | samtools sort -o hifi.bam # PacBio HiFi
minimap2 -ax sr ref.fa r1.fq.gz r2.fq.gz | samtools sort -o sr.bam  # short reads
minimap2 -x asm5 ref.fa asm.fa > asm.paf                           # assembly vs ref

# HISAT2 (RNA-seq, splice-aware)
hisat2-build ref.fa ref_idx
hisat2 -p 8 -x ref_idx -1 r1.fq.gz -2 r2.fq.gz | samtools sort -o s1.bam

# STAR (RNA-seq, splice-aware)
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir star_idx \
     --genomeFastaFiles ref.fa --sjdbGTFfile annot.gtf --sjdbOverhang 100
STAR --runThreadN 8 --genomeDir star_idx \
     --readFilesIn r1.fq.gz r2.fq.gz --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate --outFileNamePrefix s1_
```

---

## 11. samtools (SAM/BAM/CRAM)

```bash
samtools view -h in.bam | less -S            # view with header
samtools view -H in.bam                      # header only
samtools view -b -o out.bam in.sam           # SAM → BAM
samtools view -C -T ref.fa -o out.cram in.bam # BAM → CRAM
samtools sort -@ 8 -o sorted.bam in.bam      # coordinate sort
samtools sort -n -o namesorted.bam in.bam    # sort by read name
samtools index sorted.bam                    # creates .bai

# statistics
samtools flagstat in.bam                     # mapping summary
samtools idxstats in.bam                     # reads per chromosome
samtools stats in.bam > in.stats             # detailed (works with multiqc)
samtools coverage in.bam                     # coverage per contig
samtools depth -a in.bam > depth.txt         # per-base depth

# filtering
samtools view -b -F 4 in.bam > mapped.bam            # mapped only
samtools view -b -f 4 in.bam > unmapped.bam          # unmapped only
samtools view -b -q 30 in.bam > mapq30.bam           # MAPQ ≥ 30
samtools view -b -f 2 in.bam > proper_pairs.bam      # properly paired
samtools view -b in.bam chr1:1000000-2000000 > region.bam   # needs index
samtools view -c in.bam                              # count reads

# duplicates (needs fixmate on name-sorted input)
samtools sort -n -o ns.bam in.bam
samtools fixmate -m ns.bam fm.bam
samtools sort -o pos.bam fm.bam
samtools markdup -r pos.bam dedup.bam                # -r removes dups

# other
samtools merge merged.bam a.bam b.bam
samtools fastq -1 r1.fq.gz -2 r2.fq.gz in.bam        # BAM → FASTQ
```

**Common SAM flags** (decode any value at the [Broad flag explainer](https://broadinstitute.github.io/picard/explain-flags.html)):

| Flag | Meaning |
|---|---|
| 1 | read paired |
| 2 | properly paired |
| 4 | read unmapped |
| 8 | mate unmapped |
| 16 | read on reverse strand |
| 256 | secondary alignment |
| 1024 | PCR/optical duplicate |
| 2048 | supplementary alignment |

---

## 12. Variant Calling & bcftools (VCF)

```bash
# call variants
bcftools mpileup -Ou -f ref.fa s1.bam | bcftools call -mv -Oz -o s1.vcf.gz
bcftools index s1.vcf.gz

# inspect
bcftools view -h s1.vcf.gz | tail             # header
bcftools view -H s1.vcf.gz | head             # records only
bcftools view -H s1.vcf.gz | wc -l            # number of variants
bcftools stats s1.vcf.gz > s1.stats
bcftools query -l s1.vcf.gz                   # sample names

# filter
bcftools view -i 'QUAL>30 && INFO/DP>10' s1.vcf.gz -Oz -o filt.vcf.gz
bcftools view -v snps s1.vcf.gz               # SNPs only (indels: -v indels)
bcftools view -r chr1:1000-5000 s1.vcf.gz     # region (needs index)
bcftools view -s sampleA multi.vcf.gz         # one sample

# extract fields to a table
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n' s1.vcf.gz

# normalize, merge, consensus
bcftools norm -f ref.fa -m -both s1.vcf.gz -Oz -o norm.vcf.gz
bcftools merge a.vcf.gz b.vcf.gz -Oz -o merged.vcf.gz
bcftools consensus -f ref.fa s1.vcf.gz > s1.consensus.fa
```

Output types: `-Ov` VCF, `-Oz` gzipped VCF, `-Ob` BCF, `-Ou` uncompressed BCF (fast for pipes).

---

## 13. bedtools (Intervals)

Most commands want **sorted** input: `sort -k1,1 -k2,2n in.bed > in.sorted.bed`

```bash
bedtools intersect -a a.bed -b b.bed               # overlapping parts
bedtools intersect -a a.bed -b b.bed -u            # A features with any overlap
bedtools intersect -a a.bed -b b.bed -v            # A features with NO overlap
bedtools intersect -a a.bed -b b.bed -wa -wb       # report both
bedtools intersect -a reads.bam -b genes.bed > in_genes.bam

bedtools merge -i sorted.bed                       # merge overlapping
bedtools merge -i sorted.bed -d 100                # merge within 100 bp
bedtools closest -a a.sorted.bed -b b.sorted.bed   # nearest feature
bedtools slop -i a.bed -g genome.sizes -b 500      # extend ±500 bp
bedtools flank -i a.bed -g genome.sizes -l 2000 -r 0 -s   # upstream (promoters)
bedtools complement -i sorted.bed -g genome.sizes  # gaps
bedtools getfasta -fi ref.fa -bed a.bed -fo a.fa   # sequences (-s for strand)
bedtools coverage -a genes.bed -b reads.bam        # coverage per feature
bedtools genomecov -ibam in.bam -bg > cov.bedgraph # coverage track
bedtools bamtobed -i in.bam > reads.bed
bedtools makewindows -g genome.sizes -w 10000 > windows.bed
```

---

## 14. GTF/GFF Annotation

```bash
cut -f3 annot.gtf | grep -v "^#" | sort | uniq -c   # feature types
awk '$3=="gene"' annot.gtf | wc -l                 # number of genes

# gene IDs
grep -o 'gene_id "[^"]*"' annot.gtf | cut -d'"' -f2 | sort -u > gene_ids.txt

# genes → BED (GTF 1-based → BED 0-based)
awk 'BEGIN{OFS="\t"} $3=="gene" {
  match($0, /gene_id "[^"]+"/); id=substr($0, RSTART+9, RLENGTH-10)
  print $1, $4-1, $5, id, ".", $7}' annot.gtf > genes.bed

# GFF3 → GTF
gffread annot.gff3 -T -o annot.gtf

# transcript sequences from genome + annotation
gffread -w transcripts.fa -g ref.fa annot.gtf
```

---

## 15. RNA-seq Quantification

```bash
# featureCounts (Subread) — counts from aligned BAMs
featureCounts -T 8 -a annot.gtf -o counts.txt -p --countReadPairs *.bam
# -s 0/1/2 for unstranded/stranded/reverse-stranded libraries

# HTSeq
htseq-count -f bam -r pos -s no s1.bam annot.gtf > s1.counts

# Salmon (alignment-free)
salmon index -t transcripts.fa -i salmon_idx -p 8
salmon quant -i salmon_idx -l A -1 r1.fq.gz -2 r2.fq.gz -p 8 -o quant/s1

# kallisto
kallisto index -i kallisto.idx transcripts.fa
kallisto quant -i kallisto.idx -o kq/s1 -t 8 r1.fq.gz r2.fq.gz
```

Downstream differential expression is usually done in R (DESeq2, edgeR) or Python (PyDESeq2).

---

## 16. BLAST

```bash
makeblastdb -in db.fa -dbtype nucl -out mydb          # -dbtype prot for proteins

blastn -query q.fa -db mydb -outfmt 6 -evalue 1e-5 \
       -max_target_seqs 5 -num_threads 8 -out hits.tsv
blastp -query q.faa -db prot_db -outfmt 6 -out hits.tsv
blastx  # nucleotide query vs protein db
tblastn # protein query vs nucleotide db

# custom columns
blastn -query q.fa -db mydb \
  -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" -out hits.tsv

# remote search against NCBI (small queries only)
blastn -query q.fa -db nt -remote -outfmt 6 -out remote_hits.tsv
```

`-outfmt 6` columns: `qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore`

Faster protein alternative: **DIAMOND** (`diamond makedb --in db.faa -d db` → `diamond blastp -d db -q q.faa -o hits.tsv`).

---

## 17. Environments: conda / mamba

Install [Miniforge](https://github.com/conda-forge/miniforge) (includes `mamba`), then:

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

mamba create -n bio samtools bcftools bedtools seqkit fastp fastqc multiqc bwa minimap2
conda activate bio
conda deactivate

conda env list                         # list environments
conda list                             # packages in current env
mamba install -n bio salmon            # add a package
mamba search -c bioconda star          # available versions

conda env export --from-history > environment.yml   # share your env
mamba env create -f environment.yml                 # recreate it
conda env remove -n bio
```

Containers: `docker run --rm -v $PWD:/data biocontainers/samtools:v1.9-4-deb_cv1 samtools --version` · on HPC use `apptainer`/`singularity`.

---

## 18. Loops & Parallelism

```bash
# loop over paired samples
for r1 in data/*_R1.fq.gz; do
  s=$(basename "$r1" _R1.fq.gz)
  r2=data/${s}_R2.fq.gz
  bwa mem -t 8 ref.fa "$r1" "$r2" | samtools sort -o aligned/${s}.bam
  samtools index aligned/${s}.bam
done

# loop over a sample list
while read s; do echo "processing $s"; done < samples.txt

# GNU parallel
parallel -j 4 fastqc {} ::: *.fq.gz
parallel -j 4 'samtools index {}' ::: *.bam
ls *.bam | parallel -j 4 'samtools flagstat {} > {.}.flagstat'

# xargs
ls *.bam | xargs -P 4 -I {} samtools index {}

# strict, safer scripts
#!/usr/bin/env bash
set -euo pipefail
```

---

## 19. HPC: SLURM, tmux, Transfers

```bash
# SLURM
sbatch job.sh                 # submit
squeue -u $USER               # my jobs
scancel <jobid>               # cancel
sinfo                         # partitions / nodes
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS   # usage after finish
srun --pty -c 4 --mem 16G -t 2:00:00 bash              # interactive session
```

Minimal `job.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=align
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out

source ~/miniforge3/etc/profile.d/conda.sh
conda activate bio
bwa mem -t $SLURM_CPUS_PER_TASK ref.fa r1.fq.gz r2.fq.gz | samtools sort -o s1.bam
```

```bash
# keep sessions alive after disconnecting
tmux new -s work              # new session
# Ctrl-b d                    # detach
tmux ls                       # list sessions
tmux attach -t work           # reattach

nohup ./run.sh > run.log 2>&1 &    # run in background

# transfers
scp file.bam user@server:/path/
rsync -avP data/ user@server:/path/data/   # resumable, only copies changes
ssh user@server
htop                                        # monitor CPU/RAM
nproc; free -h                              # cores; memory
```

---

## 20. Data Integrity & Good Practice

```bash
md5sum *.fq.gz > checksums.md5     # create checksums
md5sum -c checksums.md5            # verify after download/transfer
gzip -t file.gz && echo OK         # test gzip integrity
samtools quickcheck *.bam && echo "BAMs OK"

samtools --version | head -1       # record tool versions
```

**Habits that save you later**

- Never edit raw data — keep `data/raw/` read-only (`chmod -R a-w data/raw`).
- One folder per project: `data/`, `results/`, `scripts/`, `logs/`, `envs/`.
- Script everything; keep a `README` with commands and versions.
- Check chromosome naming matches across files (`chr1` vs `1`).
- Remember: **BED is 0-based**, **GTF/VCF/SAM are 1-based**.
- For bigger pipelines look at **Snakemake** or **Nextflow** (nf-core).

---

## Useful Links

- [samtools / bcftools docs](https://www.htslib.org/doc/)
- [bedtools docs](https://bedtools.readthedocs.io/)
- [seqkit docs](https://bioinf.shenwei.me/seqkit/)
- [Bioconda](https://bioconda.github.io/)
- [nf-core pipelines](https://nf-co.re/)
- [Explain SAM flags](https://broadinstitute.github.io/picard/explain-flags.html)

---

*Pull requests and additions welcome.*
