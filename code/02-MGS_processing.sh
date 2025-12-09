#!/bin/bash

# Setup ----

cd ~/CSA/MGS
conda activate /home/robyn/anaconda3/envs/biobakery3


# Process reads ----

parallel -j 1 --link `# Run kneaddata` \
  'kneaddata -i {1} -i {2} -o kneaddata/out/ \
    -db /home/shared/bowtiedb/GRCh38_PhiX \
    --trimmomatic /usr/local/prg/Trimmomatic-0.36/ \
    -t 4 \
    --trimmomatic-options "SLIDINGWINDOW:4:20 MINLEN:50" \
    --bowtie2-options "--very-sensitive --dovetail" \
    --remove-intermediate-output' \
  ::: ~/CAIS-microbiome/data/raw/LangilleCSABedrestMGS_RunNS125/*_R1.fastq \
  ::: ~/CAIS-microbiome/data/raw/LangilleCSABedrestMGS_RunNS125/*_R2.fastq
    
concat_paired_end.pl `# Concatenate reads` \
  -p 4 \
  --no_R_match \
  -o cat_reads \
  kneaddata/out/*_paired_*.fastq


# Functional analysis ----

for file in kneaddata/cat_reads/* `# Run HUMAnN3 on all samples``
  do humann \
    --input ${file} \
    --output humann3/raw_output \
    --search-mode uniref90 \
    --threads 18
  done

humann_join_tables `# Join pathway abundance data into single file` \
  --input humann3/raw_output \
  --file_name pathabundance
  --output humann3/pathabundance.tsv
  done

humann_split_stratified_table `# Unstratify data` \
  --input pathabundance.tsv \
  --output .