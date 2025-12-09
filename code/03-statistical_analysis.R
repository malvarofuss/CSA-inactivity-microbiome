# Setup ----

set.seed(0509)

library(tidyverse)
library(glue)
library(rstatix) 
library(vegan)
library(phyloseq)
library(microbiome)
library(lme4)
library(lmerTest)

parse_sample_metadata <- function(data) { 
  data %>%
    mutate(
      participant_id = str_sub(sample_id, 1, 4),
      sample_type = case_when(
        str_detect(sample_id, "G") ~  "gut",
        str_detect(sample_id, "O") ~  "oral"),
      phase = case_when(
        str_detect(sample_id, "BDC") ~ "baseline",
        str_detect(sample_id, "HDT") ~ "bedrest",
        str_detect(sample_id, "R") ~ "recovery",
        str_detect(sample_id, "W4") ~ "4w",
        str_detect(sample_id, "M4") ~ "4m"),
      day = case_when(
        str_detect(sample_id, "BDC") ~  NA,
        str_detect(sample_id, "HDT") ~ str_extract(sample_id, "\\d+$"),
        str_detect(sample_id, "R") ~ str_extract(sample_id, "\\d+$"),
        str_detect(sample_id, "W4") ~ NA,
        str_detect(sample_id, "M4") ~ NA) %>% as.numeric())
}


# Participant summary ----

participants <- # Load participant metadata
  read_tsv("data/participant_metadata.tsv", col_types = c("ffffDnnnDn")) %>%
  mutate(cohort = str_replace(cohort, "Group", "Cohort"))

participants %>% # Calculate participant BMI
  mutate(bmi = weight_kg / ((height_cm / 100)^2))

participants %>% # Calculate mean age across all participants
  summarise(mean_age = mean(age))

frailty <- # Load frailty data
  read_tsv("data/frailty_data.tsv", col_types = c("nfci"))

participants %>% # Participant summary by sex and exercise group
  inner_join( # Add baseline FI-36 measurements
    frailty %>% filter(phase == "baseline")) %>%
  group_by(exercise, sex) %>% # Group data
  summarise( # Summarize age, height, weight and BMI (mean and SD)
    age_mean = mean(age),
    age_sd = sd(age), 
    height_mean = mean(height_cm),
    height_sd = sd(height_cm),
    weight_mean = mean(weight_kg),
    weight_sd = sd(weight_kg),
    bmi_mean = mean(bmi),
    bmi_sd = sd(bmi),
    fi_mean = mean(fi36),
    fi_sd = sd(fi36)) %>%
  write_tsv("results/participant_summary.tsv")

participants %>% # Kruskal-Wallis test for BMI
    mutate(group = interaction(exercise, sex)) %>%
    kruskal_test(bmi ~ group) 

participants %>% # Calculate effect size for BMI KW test
    mutate(group = interaction(exercise, sex)) %>%
    kruskal_effsize(bmi ~ group)

participants %>% # Kruskal-Wallis test for baseline FI-36
    inner_join(frailty %>% filter(phase == "baseline")) %>%
    mutate(group = interaction(exercise, sex)) %>%
    kruskal_test(fi36 ~ group)


# Family relative abundance ----

feature_tables <- list() # More feature tables will be added for DA analysis

# Transform to relative abundance
feature_tables$family_relab <-
  read_tsv("data/16S/collapse/table_family.tsv", skip = 1) %>% # Load data
  rename(`feature_name` = `#OTU ID`) %>%
  column_to_rownames("feature_name") %>% 
  otu_table(taxa_are_rows = TRUE) %>% # Transform to relative abundance
  phyloseq() %>% 
  transform("compositional") %>% 
  as.data.frame.matrix() %>% 
  rownames_to_column(var = "feature_name") %>%
  pivot_longer( # Transform to long format
    cols = -1, 
    names_to = "sample_id", 
    values_to = "relab") %>% 
  parse_sample_metadata() %>% 
  inner_join(participants) # Append participant metadata

# Calculate top 6 most abundant families
feature_tables$family_relab %>% 
  group_by(sample_type, feature_name) %>%
  summarise( # Calculate total abundance per family
    `total_abundance` = sum(c(as.numeric(`relab`)), na.rm = TRUE)) %>%
  distinct(feature_name, total_abundance) %>%
  group_by(sample_type) %>%
  slice_max(order_by = total_abundance, n = 6) %>%
  write_tsv("results/top_6_families.tsv")


# Beta diversity analysis ----

# Load beta diversity results
beta_diversity <- list( 
  dm_gut = # Gut distance matrix
    read_tsv("data/16S/diversity/dm_gut.tsv") %>%
      rename("sample_id" = 1),
  dm_oral = # Oral distance matrix
    read_tsv("data/16S/diversity/dm_oral.tsv") %>%
      rename("sample_id" = 1),
  pcoa_gut = # Gut ordination plot
    read_tsv(
      "data/16S/diversity/ordination_gut.txt",
      skip = 7, col_names = c("sample_id", "PC1", "PC2")) %>%
      filter(str_detect(sample_id, "PT")) %>%
      mutate(across(2:3, as.numeric)), 
  pcoa_oral = # Oral ordination plot
    read_tsv(
      "data/16S/diversity/ordination_oral.txt",
      skip = 7, , col_names = c("sample_id", "PC1", "PC2"))  %>%
      filter(str_detect(sample_id, "PT")) %>%
      mutate(across(2:3, as.numeric))) %>%
  map(parse_sample_metadata) %>% # Parse sample ID for all results
  map(~ left_join(.x, participants)) # Join participant metadata

# Format beta-diversity matrices as distances for PERMANOVA 
beta_diversity$dist_gut <- # Gut samples
  beta_diversity$dm_gut %>% 
  select(all_of(.[["sample_id"]])) %>% 
  as.dist()
beta_diversity$dist_oral <- # Oral samples
  beta_diversity$dm_oral %>% 
  select(all_of(.[["sample_id"]])) %>% 
  as.dist()

# Calculate association with participant ID, sex, cohort and exercise group
calculate_permanova <- function(sample_type, study_variable) {
  set.seed(0509)
  dist_mat <- beta_diversity[[paste0("dist_", sample_type)]]
  dm_mat   <- beta_diversity[[paste0("dm_", sample_type)]] 
  
  formula <- as.formula(paste("dist_mat ~", study_variable))
  
  adonis2(formula, data = dm_mat, method = "robust.aitchison") %>%
    rownames_to_column(var = "rowname") %>%
    write_tsv(
      glue("results/permanovas/{sample_type}_{study_variable}.tsv"))
}

expand_grid( # Run all PERMANOVA tests
  c("gut", "oral"), 
  c("participant_id", "cohort", "sex", "exercise")) %>%
pmap(~ calculate_permanova(..1, ..2))

# Association with FI-36 (uncorrected model)
beta_diversity$dist_gut_frailty <- # Format distances (gut)
  beta_diversity$dm_gut %>%
  right_join(frailty) %>%
  filter(!is.na(`fi36`), !is.na(sample_id)) %>%
  select(all_of(.[["sample_id"]])) %>%
  as.dist()
beta_diversity$dist_oral_frailty <- # Format distances (oral)
  beta_diversity$dm_oral %>%
  right_join(frailty) %>%
  filter(!is.na(`fi36`), !is.na(sample_id)) %>%
  select(all_of(.[["sample_id"]])) %>%
  as.dist()
set.seed(0509)
adonis2( # Calculate PERMANOVA (gut)
  beta_diversity$dist_gut_frailty ~ fi36, 
  by = "margin",
  method = "robust.aitchison",
  data = beta_diversity$dm_gut %>% 
    right_join(frailty) %>%
    filter(!is.na(`fi36`), !is.na(`sample_id`)) %>%
    select(all_of(.[["sample_id"]]), fi36)) %>%
  rownames_to_column(var = "rowname") %>%
  write_tsv("results/permanovas/gut_frailty_uncorrected.tsv")
set.seed(0509)
adonis2( # Calculate PERMANOVA (oral)
  beta_diversity$dist_oral_frailty ~ fi36, 
  by = "margin",
  method = "robust.aitchison",
  data = beta_diversity$dm_oral %>% 
    right_join(frailty) %>%
    filter(!is.na(`fi36`), !is.na(`sample_id`)) %>%
    select(all_of(.[["sample_id"]]), fi36)) %>%
  rownames_to_column(var = "rowname") %>%
  write_tsv("results/permanovas/oral_frailty_uncorrected.tsv")

# Association with FI-36 (model adjusted for participant ID)
set.seed(0509)
adonis2( # Calculate PERMANOVA (gut)
  beta_diversity$dist_gut_frailty ~ fi36 + participant_id, 
  by = "margin",
  method = "robust.aitchison",
  data = beta_diversity$dm_gut %>% 
    right_join(frailty) %>%
    filter(!is.na(`fi36`), !is.na(`sample_id`)) %>%
    select(all_of(.[["sample_id"]]), fi36, participant_id)) %>%
  rownames_to_column(var = "rowname") %>%
  write_tsv("results/permanovas/gut_frailty_corrected.tsv")
set.seed(0509)
adonis2( # Calculate PERMANOVA (oral)
  beta_diversity$dist_oral_frailty ~ fi36 + participant_id, 
  by = "margin",
  method = "robust.aitchison",
  data = beta_diversity$dm_oral %>% 
    right_join(frailty) %>%
    filter(!is.na(`fi36`), !is.na(`sample_id`)) %>%
    select(all_of(.[["sample_id"]]), fi36, , participant_id)) %>%
  rownames_to_column(var = "rowname") %>%
  write_tsv("results/permanovas/oral_frailty_corrected.tsv")

# Association with FI-36 (outlier removed, model adjusted for participant ID)
beta_diversity$dist_gut_frailty_noout <- 
  beta_diversity$dm_gut %>% 
  right_join(frailty) %>%
  filter(!is.na(`fi36`), !is.na(sample_id)) %>%
  filter(!is.na(sample_id)) %>%
  filter(fi36 < 0.6) %>% # Remove high-frailty outlier
  select(all_of(.[["sample_id"]])) %>% 
  as.dist()
beta_diversity$dist_oral_frailty_noout <- 
  beta_diversity$dm_oral %>% 
  right_join(frailty) %>%
  filter(!is.na(`fi36`), !is.na(sample_id)) %>% 
  filter(fi36 < 0.6) %>% # Remove high-frailty outlier
  select(all_of(.[["sample_id"]])) %>% 
  as.dist()
set.seed(0509)
adonis2( # Gut samples
  beta_diversity$dist_gut_frailty_noout ~ fi36 + participant_id, 
  by = "margin",
  method = "robust.aitchison",
  data = beta_diversity$dm_gut %>% 
    right_join(frailty) %>%
    filter(!is.na(`fi36`), !is.na(sample_id)) %>% 
    filter(fi36 < 0.6) %>% # Remove high-frailty outlier
    select(all_of(.[["sample_id"]]), fi36, participant_id)) %>%
  rownames_to_column(var = "rowname") %>%
  write_tsv("results/permanovas/gut_frailty_corrected_noout.tsv")
set.seed(0509)
adonis2( # Gut samples
  beta_diversity$dist_oral_frailty_noout ~ fi36 + participant_id, 
  by = "margin",
  method = "robust.aitchison",
  data = beta_diversity$dm_oral %>% 
    right_join(frailty) %>%
    filter(!is.na(`fi36`), !is.na(sample_id)) %>% 
    filter(fi36 < 0.6) %>% # Remove high-frailty outlier
    select(all_of(.[["sample_id"]]), fi36, , participant_id)) %>%
  rownames_to_column(var = "rowname") %>%
  write_tsv("results/permanovas/oral_frailty_corrected_noout.tsv")

# Alpha diversity analysis

# Load alpha diversity data
alpha_diversity <- 
  c("observed_features", "shannon", "faith") %>% # Load all metrics 
  map_dfr(~ read_tsv(
    glue("data/16S/diversity/{.x}.tsv"),
    skip = 1, col_types = "cd",
    col_names = c("sample_id", "value")) %>%
  mutate(metric = .x)) %>%
  parse_sample_metadata() %>% 
  left_join(participants) # Append participant metadata

# Calculate effect of bedrest using linear mixed models
alpha_diversity_analysis <- function(formula, outfile) {

  set.seed(0509)
  
  alpha_diversity %>%
    filter(phase == "bedrest") %>%
    group_by(sample_type, exercise, metric) %>%
    summarise(
      model = list(
        tryCatch({
          lmer(formula, data = cur_data())}, error = function(e) "error"))) %>%
    mutate(
      .keep = "unused",
      coefficients = map(
        model, ~ .x %>% summary() %>% coefficients() %>% 
          as_tibble() %>% slice(n()))) %>%
    unnest_wider("coefficients") %>%
    rownames_to_column(var = "rowname") %>%
    write_tsv(outfile)
}

alpha_diversity_analysis( # Run models for effect of day of HDBR
  formula = value ~ cohort + sex + day + (1 + day | participant_id),
  outfile = "results/alpha_diversity_day_bedrest.tsv")

alpha_diversity_analysis( # Run models for interaction effect with sex
  formula = value ~ cohort + sex * day + (1 + day | participant_id),
  outfile = "results/alpha_diversity_sex_interaction.tsv")

# Differential abundance analysis ----

# Load feature count data
feature_tables <- append(feature_tables, list(
  gut = # Gut genus count feature table
    read_tsv("data/16S/collapse/table_genus_gut.tsv", skip = 1),
  oral = # Oral genus count feature table
    read_tsv("data/16S/collapse/table_genus_oral.tsv", skip = 1),
  pathway = read_tsv( # Pathway count feature table
    "data/MGS/humann3/unstratify/pathabundance_unstratified.tsv") %>%
    rename_with(~ str_remove(.x, "_Abundance-RPKs")) %>% 
    rename_with(~ str_remove(.x, "_Abundance")) %>%
    rename_with(~ str_replace(.x, "HDT([1-9])$", "HDT0\\1")), 
  metabolite = # Metabolite abundance feature table
    read_tsv("data/metabolites/metabolites.txt") %>%
    filter( # Keep microbiota-related, tier 1 metabolites only
      `Category` == "Microbiota-related",
      `ID Level` == "Tier 1") %>%
    rename_with( 
     ~ str_replace(.x, "HDT([1-9])$", "HDT0\\1")) %>%
    select(-1, -2) %>%
    rename(`feature_name` = 1)))

# Transform count data to CLR abundance
transform_count_to_CLR <- function(count_table) {
  count_table %>%
  rename(`feature_name` = 1) %>% # Reformat feature table
  column_to_rownames("feature_name") %>% 
  otu_table(taxa_are_rows = TRUE) %>% # Import to phyloseq object
  phyloseq() %>% 
  transform("clr") %>%  # Apply CLR transformation
  as.data.frame.matrix() %>% 
  rownames_to_column(var = "feature_name") # Return to original format
}

feature_tables$gut_clr <- # 16S gut samples
  feature_tables$gut %>% 
  transform_count_to_CLR()

feature_tables$oral_clr <- # 16S oral samples
  feature_tables$oral %>% 
  transform_count_to_CLR()

feature_tables$pathway_clr <- # Metegenomic samples
  feature_tables$pathway %>% 
  transform_count_to_CLR()

# Remove low prevalence features
remove_low_prevalence_features <- function(clr_table, count_table, cut_off) {
  
  low_prevalence_features <- # Identify features with prevalence > 0.1
    count_table %>%
    rename(`feature_name` = 1) %>% 
    pivot_longer( # Pivot to long format
      cols = -1, 
      names_to = "sample_id", 
      values_to = "count") %>%
    group_by(feature_name) %>% 
    mutate( # Calculate prevalence of each feature across all samples
      `prevalence` = sum(c(as.numeric(`count`)) > 0, na.rm = TRUE)) %>%
    filter(prevalence > {{ cut_off }}) %>% # Filter under cut-off
    pull(feature_name) %>% 
    unique()
    
  # Filter CLR table
  clr_table %>% filter(`feature_name` %in% low_prevalence_features)
}

feature_tables$gut_clr_prev10 <- # Gut-associated genera
  remove_low_prevalence_features(
    feature_tables$gut_clr, feature_tables$gut, 34)

feature_tables$oral_clr_prev10 <- # Oral-associated genera
  remove_low_prevalence_features(
    feature_tables$oral_clr, feature_tables$oral, 34)

feature_tables$pathway_clr_prev10 <- # MetaCyc pathways
  remove_low_prevalence_features(
    feature_tables$pathway_clr, feature_tables$pathway, 9)

# Calculate differential abundance
calculate_differential_abundance <- function(feature_table) {
  feature_table %>%
  pivot_longer( # Pivot to long format
    cols = -1, 
    names_to = "sample_id", 
    values_to = "abundance") %>%
  parse_sample_metadata() %>% # Add sample and participant metadata
  left_join(participants) %>%
  filter(phase == "bedrest") %>% # Filter bedrest timepoints
  group_by(exercise, feature_name) %>% 
  summarise( # Run LMMs for all features
    model = list( 
      tryCatch({
        lmer(
          abundance ~ cohort + sex + day + (1 | participant_id))}))) %>%
  mutate( # Format results
    .keep = "unused",
    coefficients = map(
      model, ~ .x %>% 
        summary() %>% coefficients() %>% as_tibble() %>% slice(n()))) %>%
  unnest_wider("coefficients") %>%
  group_by(exercise) %>% # Correct p-values
  adjust_pvalue(p.col = "Pr(>|t|)", output.col = "p_BH", method = "BH") 
}

calculate_differential_abundance(feature_tables$gut_clr_prev10) %>% 
  filter(`p_BH` < 0.05) %>% # Save statistically significant results
  write_tsv("results/differential_abundance/genus_gut.tsv")

calculate_differential_abundance(feature_tables$oral_clr_prev10) %>%
  filter(`p_BH` < 0.05) %>% # Save statistically significant results
  write_tsv("results/differential_abundance/genus_oral.tsv")

calculate_differential_abundance( # MetaCyc pathways
  feature_tables$pathway_clr_prev10 %>% 
    filter(!feature_name %in% c("UNMAPPED", "UNINTEGRATED"))) %>% 
  filter(`Pr(>|t|)` < 0.05) %>% # Save statistically significant results
  write_tsv("results/differential_abundance/pathway_gut.tsv")

calculate_differential_abundance(feature_tables$metabolite) %>%
  filter(`Pr(>|t|)` < 0.05) %>% # Save statistically significant results
  write_tsv("results/differential_abundance/metabolite_gut.tsv")