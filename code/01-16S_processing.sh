#!/bin/bash

# Setup and activate QIIME2 environment ----

cd ~/CAIS-microbiome/data/16S
conda activate qiime2-2022.8 # Version 2022.8

# .qzv files were generated for all .qza files
mkdir ~/CAIS-microbiome/results/qzv_files/ 


# Import raw reads into QIIME2 ----

mkdir reads

qiime tools import \ 
  --type SampleData[PairedEndSequencesWithQuality] \ 
  --input-path ~/CAIS-microbiome/data/raw/LangilleCSABedrest_16S \
  --output-path reads/reads_raw.qza \ 
  --input-format CasavaOneEightSingleLanePerSampleDirFmt 

qiime demux summarize \ 
  --i-data reads/reads_raw.qza \
  --o-visualization ~/CAIS-microbiome/results/qzv_files/reads_raw.qzv


# Trim V4/V5 primers ----

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences reads/reads_raw.qza \
  --p-cores 4 \
  --p-front-f GTGYCAGCMGCCGCGGTAA \
  --p-front-r CCGYCAATTYMTTTRAGTTT \
  --p-discard-untrimmed \ 
  --p-no-indels \
  --o-trimmed-sequences reads/reads_trimmed.qza

qiime demux summarize `# Generate .qzv` \
  --i-data reads/reads_trimmed.qza \
  --o-visualization ~/CAIS-microbiome/results/qzv_files/reads_trimmed.qzv


# Denoise samples using DADA2 ----

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs reads/reads_trimmed.qza \
  --p-trunc-len-f 270 \
  --p-trunc-len-r 210 \
  --p-max-ee-f 3 \
  --p-max-ee-r 3 \
  --p-n-threads 30 \
  --output-dir denoising

qiime feature-table summarize `# Generate .qzv for table` \
  --i-table denoising/table.qza \
  --o-visualization ~/CAIS-microbiome/results/qzv_files/table.qzv

qiime feature-table tabulate-seqs `# Generate .qzv for repseqs` \
  --i-data denoising/representative_sequences.qza \
  --o-visualization ~/CAIS-microbiome/results/qzv_files/representative_sequences.qzv


# Remove possible bleed-through sequences ----

mkdir tables
mkdir repseqs

qiime feature-table filter-features \
  --i-table denoising/table.qza \
  --p-min-frequency 23 `# 0.1% of mean sample depth` \
  --p-min-samples 1 \
  --o-filtered-table tables/table_bleed.qza \
  
qiime feature-table summarize `# Generate .qzv` \
  --i-table tables/table_bleed.qza \
  --o-visualization ~/CAIS-microbiome/results/qzv_files/table_bleed.qzv

qiime feature-table filter-seqs `# Filter representative sequences` \
   --i-data denoising/representative_sequences.qza \
   --i-table tables/table_bleed.qza \
   --o-filtered-data repseqs/repseqs_bleed.qza

qiime feature-table tabulate-seqs `# Generate .qzv` \
  --i-data qiime2/filtering/repseqs_bleed.qza \
  --o-visualization ~/CAIS-microbiome/results/qzv_files/repseqs_bleed.qzv


# Taxonomic classification ----

qiime feature-classifier classify-sklearn \
  --i-reads denoising/representative_sequences.qza \
  --i-classifier /home/shared/taxa_classifiers/qiime2-2022.8_classifiers/silva-138-99-nb-classifier.qza \
  --p-n-jobs 18 \
  --output-dir taxonomy


# Remove contaminating sequences ----

qiime taxa filter-table \
  --i-table tables/table_bleed.qza \
  --i-taxonomy taxonomy/classification.qza \
  --p-exclude mitochondria,chloroplast `# Remove mitochondria, chloroplast seqs` \
  --p-include p__ `# Remove seqs unclassified at the phylum level` \
  --o-filtered-table tables/table_bleed_decon.qza

qiime feature-table summarize `# Generate .qzv` \
  --i-table tables/table_bleed_decon.qza \
  --o-visualization ~/CAIS-microbiome/results/qzv_files/table_bleed_decon.qzv
  
qiime feature-table filter-seqs `# Filter representative sequences` \
   --i-data repseqs/repseqs_bleed.qza \
   --i-table tables/table_bleed_decon.qza \ 
   --o-filtered-data repseqs/repseqs_bleed_decon.qza

qiime feature-table tabulate-seqs `# Generate .qzv` \
  --i-data repseqs/repseqs_bleed_decon.qza \
  --o-visualization ~/CAIS-microbiome/results/qzv_files/repseqs_bleed_decon.qzv


# Split into gut and oral samples ----

for type in "gut" "oral"
  do qiime feature-table filter-samples \
    --i-table tables/table_bleed_decon.qza \
    --m-metadata-file ~/CAIS-microbiome/data/sample_metadata.tsv \
    --p-where "sample_type=='${type}'" \
    --o-filtered-table tables/table_${type}.qza
  done
  
for type in "gut" "oral"
  do qiime feature-table summarize `# Generate .qzv` \
    --i-table tables/table_${type}.qza \
    --o-visualization ~/CAIS-microbiome/results/qzv_files/table_${type}.qzv
  done
  
for type in "gut" "oral" `# Split representative sequences`
  do qiime feature-table filter-seqs \
    --i-data repseqs/repseqs_bleed_decon.qza \
    --i-table tables/table_${type}.qza \
    --o-filtered-data repseqs/repseqs_${type}.qza 
  done
  
for type in "gut" "oral" `# Generate .qzv`
  do qiime feature-table tabulate-seqs \
    --i-data qiime2/filtering/repseqs_${type}.qza \
    --o-visualization ~/CAIS-microbiome/results/qzv_files/repseqs_${type}.qzv 
  done

# Generate phylogeny ----

qiime fragment-insertion sepp `# All samples` \
  --i-representative-sequences repseqs/repseqs_bleed_decon.qza \ 
  --i-reference-database /home/shared/rRNA_db/16S/sepp-refs-gg-13-8.qza  \
  --o-tree phylogeny/tree_all.qza  \
  --o-placements phylogeny/placements_all.qza \
  --p-threads 30

for type in "gut" "oral" `# Gut and oral samples`
  do qiime fragment-insertion sepp \
    --i-representative-sequences repseqs/repseqs_${type}.qza \ 
    --i-reference-database /home/shared/rRNA_db/16S/sepp-refs-gg-13-8.qza  \
    --o-tree phylogeny/tree_${type}.qza  \
    --o-placements phylogeny/placements_${type}.qza \
    --p-threads 30
  done 

# Collapse to genus and family level ----

mkdir collapse

for type in "gut" "oral"; do `# Collapse split tables to genus level`
  qiime taxa collapse \
    --i-table tables/table_${type}.qza  \
    --i-taxonomy taxonomy/classification.qza \
    --p-level 6 \
    --o-collapsed-table collapse/table_${type}_genus.qza; done

qiime taxa collapse \ `# Collapse combined table to family level`
  --i-table tables/table_${type}.qza  \
  --i-taxonomy taxonomy/classification.qza \
  --p-level 5 \
  --o-collapsed-table collapse/table_${type}_family.qza

# Calculate diversity ----

mkdir diversity
mkdir diversity/qzas

for type in "gut" "oral"; do
  qiime gemelli phylogenetic-rpca-with-taxonomy `# Calculate phyloegentic RPCA` \
    --i-table tables/table_${type}.qza \
    --i-phylogeny phylogeny/tree_${type}.qza  \
    --m-taxonomy-file taxonomy/classification.qza \
    --o-biplot diversity/beta_analysis/rpca_pcoa_${type}.qza \
    --o-distance-matrix diversity/beta_analysis/rpca_dm_${type}.qza \
    --o-counts-by-node-tree diversity/beta_analysis/rpca_phylo-tree_${type}.qza \
    --o-counts-by-node diversity/beta_analysis/rpca_phylo-table_${type}.qza \
    --o-t2t-taxonomy diversity/beta_analysis/rpca_phylo-taxonomy_${type}.qza; done

rm diversity/beta_analysis/rpca_phylo-tree_*.qza `# Remove unused results`
rm diversity/beta_analysis/rpca_phylo-table_*.qza
rm diversity/beta_analysis/rpca_phylo-taxonomy_*.qza

conda deactivate
conda activate q2-boots-amplicon-2025.7 `# Load q2-boots to perform rarefaction`

# # Calculate alpha diversity with rarefaction
for metric in "observed_features" "shannon" "faith_pd"; do 
  qiime boots alpha \
    --i-table tables/table_bleed_decon.qza \
    --i-phylogeny phylogeny/tree_all.qza \
    --p-sampling-depth 2500 \
    --p-metric ${metric} \
    --p-n 100 \
    --p-no-replacement \
    --p-average-method median \
    --o-average-alpha-diversity diversity/qzas/${metric}.qza; done

# Calculate Bray Curtis and Unifrac with rarefaction
for type in "gut" "oral"; do 
  for metric in "braycurtis" "weighted_unifrac"; do
    qiime boots beta \
      --i-table tables/table_${type}.qza \
      --i-phylogeny phylogeny/tree_${type}.qza \
      --p-metric ${metric} \
      --p-sampling-depth 2500 \
      --p-n 100 \
      --p-no-replacement \
      --p-average-method medoid \
      --o-average-distance-matrix diversity/qzas/${metric}_dm_${type}.qza; done; done

for type in "gut" "oral"; do `# Apply principal coordinates analysis`
  for metric in "braycurtis" "weighted_unifrac"; do
    qiime diversity pcoa \
      --i-distance-matrix diversity/beta_analysis/${metric}_dm_${type}.qza \
      --o-pcoa diversity/beta_analysis/${metric}_pcoa_${type}.qza; done; done

# Export .qza files for statistical analysis ----

for infile in collapse/*; do `# Export feature tables`
  base="$(basename "$infile" .qza)"
  qiime tools export \
      --input-path ${infile} \
      --output-path . 
  biom convert \
      -i feature-table.biom \
      -o collapse/${base}.tsv \
      --to-tsv 
  rm feature-table.biom; done

for infile in diversity/qza/; do `# Export alpha diversity`
  base="$(basename "$infile" .qza)"
  qiime tools export \
    --input-path ${infile} \
    --output-path diversity
  mv diversity/alpha-diversity.tsv diversity/${base}.tsv
  rm diversity/alpha-diversity.tsv; done

for infile in diversity/qzas/*; do `# Export diversity results`
  base="$(basename "$infile" .qza)"
  qiime tools export \
    --input-path ${infile} \
    --output-path diversity
  mv diversity/alpha-diversity.tsv diversity/${base}.tsv
  mv diversity/distance-matrix.tsv diversity/${base}.tsv
  mv diversity/ordination.txt diversity/${base}.txt; done