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


# Collapse and export tables to genus and family level ----

mkdir collapse

for type in "gut" "oral"
  qiime taxa collapse \
    --i-table tables/table_${type}.qza  \
    --i-taxonomy taxonomy/classification.qza \
    --p-level 6 \
    --o-collapsed-table collapse/table_${type}_genus.qza
    done
  done

for type in "gut" "oral" `# Export .qza to .tsv`
  do 
    qiime tools export \
      --input-path collapse/table_${type}_genus.qza \
      --output-path . \
    biom convert 
      -i feature-table.biom \
      -o collapse/table_${type}_genus.tsv \
      --to-tsv 
    rm feature-table.biom
  done

qiime taxa collapse `# Collapse combined gut-oral table to family level` \
  --i-table tables/table_bleed_decon.qza  \
  --i-taxonomy taxonomy/classification.qza \
  --p-level 5 \
  --o-collapsed-table collapse/table_family.qza

qiime tools export `# Export family feature table` \
  --input-path collapse/table_family.qza \
  --output-path collapse
biom convert \
  -i collapse/feature-table.biom \
  -o collapse/table_family.tsv \
  --to-tsv 
rm feature-table.biom


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

# Repeated rarefaction ----

mkdir rarefaction

for ((i=1; i<=100; i++))
  do qiime feature-table rarefy 
    --i-table tables/table_bleed_decon.qza \
    --p-sampling-depth 2500 `# Rarefy to 2,500 reads` \
    --o-rarefied-table rarefaction/table_${i}.qza
  done


# Calculate and export alpha diversity ----

mkdir diversity
mkdir diversity/alpha_vectors

for metric in "observed_features" "shannon" `# Non-phylogenetic metrics`
  do 
    for ((i=1; i<=100; i++)) # Calculate diversity for all subsampled tables
      do 
        qiime diversity alpha \
          --i-table rarefaction/table_${i}.qza \
          --p-metric ${metric} \
          --o-alpha-diversity diversity/alpha_vectors/${metric}_${i}.qza 
        # Export diversity for each subsampled table
        qiime tools export \ 
          --input-path diversity/alpha_vectors/${metric}_${i}.qza \
          --output-path diversity 
        # Copy sample diversity to temporary diversity table
        cat diversity/alpha-diversity.tsv >> diversity/${metric}_temp.tsv 
        # Calculate average diversity for new subsample with previous average
        awk '{
          sum[$1] += $2;
          count[$1]++;
        }
        END {
          for (sample_id in sum) {
            print sample_id, sum[sample_id] / count[sample_id];
          }
        }' diversity/${metric}_temp.tsv \
        | sort > diversity/${metric}.tsv
    done
  done
  
for ((i=1; i<=100; i++)) `# Calculate diversity for all subsampled tables`
  do 
    qiime diversity alpha-phylogenetic `# Non-phylogenetic metrics` \
      --i-table rarefaction/table_${i}.qza \
      --i-phylogeny phylogeny/tree_${type}.qza \
      --p-metric 'faith_pd' \
      --o-alpha-diversity diversity/alpha_vectors/faith_${i}.qza 
    # Export diversity for each subsampled table
    qiime tools export \
      --input-path diversity/alpha_vectors/${metric}_${i}.qza \
      --output-path diversity 
    # Copy sample diversity to temporary diversity table
    cat diversity/alpha-diversity.tsv >> diversity/${metric}_temp.tsv 
    # Calculate average diversity for new subsample with previous average
    awk '{
      sum[$1] += $2;
      count[$1]++;
    }
    END {
      for (sample_id in sum) {
        print sample_id, sum[sample_id] / count[sample_id];
      }
    }' diversity/${metric}_temp.tsv \
    | sort > diversity/${metric}.tsv
  done

# Calculate and export beta-diversity ----

mkdir gemelli 

for type in "gut" "oral"
  do
    qiime gemelli phylogenetic-rpca-with-taxonomy `# Calculate phyloegentic RPCA` \
      --i-table tables/table_bleed_decon.qza \
      --i-phylogeny phylogeny/tree_${type}.qza  \
      --m-taxonomy-file taxonomy/classification.qza \
      --o-biplot diversity/gemelli/biplot_${type}.qza \
      --o-distance-matrix diversity/gemelli/dm_${type}.qza \
      --o-counts-by-node-tree diversity/gemelli//phylo-tree_${type}.qza \
      --o-counts-by-node diversity/gemelli//phylo-table_${type}.qza \
      --o-t2t-taxonomy diversity/gemelli//phylo-taxonomy_${type}.qza \
  done

for type in "gut" "oral" # Export distance matrices
  do 
    qiime tools export \
      --input-path diversity/gemelli/dm_${type}.qza \
      --output-path diversity 
    mv diversity/distance-matrix.tsv diversity/dm_${type}.txt
  done
  
for type in "gut" "oral" # Export ordination plots
  do 
    qiime tools export \
      --input-path diversity/gemelli/biplot_${type}.qza \
      --output-path diversity 
    mv diversity/ordination.txt diversity/ordination_${type}.txt
  done