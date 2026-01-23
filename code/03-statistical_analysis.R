# 03 - Statistical analysis

## Setup ----

setwd(".")
set.seed(0509)

library(tidyverse) # Load libraries
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

## Section 3.1 (Results) ----

participants <- # Load participant metadata
  read_tsv("data/participant_metadata.tsv", col_types = c("ffffDnnnDn")) %>%
  mutate(cohort = str_replace(cohort, "Group", "Cohort"))

frailty <- list(# Load frailty data
  index = read_tsv("data/frailty_data.tsv", col_types = c("nfci")),
  frailty_components =
    read_csv("data/frailty_components.csv") %>%
    mutate(.keep = "unused",
           participant_id = `Participant`,
           phase = case_when(
             `Timeline` == "Baseline" ~ "baseline",
             str_detect(`Timeline`, "HDT") ~ "bedrest",
             str_detect(`Timeline`, "Recovery") ~ "recovery",
             `Timeline` == "4 Weeks" ~ "4w",
             `Timeline` == "4 Months" ~ "4m"),
           day = case_when(
             str_detect(`Timeline`, "HDT") ~ str_extract(`Timeline`, "\\d+$"),
             str_detect(`Timeline`, "R") ~ str_extract(`Timeline`, "\\d+$")) %>%
             as.numeric()) %>%
    mutate(across(1:33, ~ case_when(. == 2 ~ 0, TRUE ~ .)))) %>%
  reduce(full_join)

### Participant and frailty summary ----

participants %>% # Calculate mean age across all participants
  summarise(mean_age = mean(age), sd_age = sd(age))

participants %>% # Calculate participant BMI
  mutate(bmi = weight_kg / ((height_cm / 100)^2))

participants %>% # Grouped summary for age, height, weight, and BMI
  group_by(exercise, sex) %>% 
  summarise( # Calculate mean and SD
    mean(age), sd(age), mean(height_cm), sd(height_cm),
    mean(weight_kg), sd(weight_kg), mean(bmi), sd(bmi))

inner_join(participants, frailty) %>% # Participant frailty at baseline
  filter(phase == "baseline") %>%
  group_by(exercise, sex) %>%
  summarise(mean(fi36), sd(fi36))

inner_join(participants, frailty) %>% # Change in frailty during bedrest
  filter(phase == "bedrest", day == 14) %>%  
  summarise(median(fi36))

participants %>% # Kruskal-Wallis test for BMI
    mutate(group = interaction(exercise, sex)) %>%
    kruskal_test(bmi ~ group) 

participants %>% # Calculate effect size for BMI KW test
    mutate(group = interaction(exercise, sex)) %>%
    kruskal_effsize(bmi ~ group)

participants %>% # Kruskal-Wallis test for baseline frailty
    inner_join(frailty %>% filter(phase == "baseline")) %>%
    mutate(group = interaction(exercise, sex)) %>%
    kruskal_test(fi36 ~ group)

### Sample missingness ----

samples <- list( # Load sample metadata
  gut = read_csv("data/16S/sample_depth_gut.csv"),
  oral = read_csv("data/16S/sample_depth_oral.csv")) 

samples$missingness_bedrest <- samples %>% # Calculate missingess 
  map(~ .x %>% 
    filter(str_detect(sample_id, "HDT")) %>%
    parse_sample_metadata() %>%
    complete(day, participant_id, sample_type) %>%
    mutate(present = as.integer(!is.na(`seq_depth`)))) %>%
  reduce(rbind)

samples$missingness_bedrest %>% # Missingness distribution by participant
  group_by(participant_id, sample_type) %>%
  summarise(perc_missing = mean(present == 0)) %>%
  group_by(sample_type) %>%
  summarise(
    median_missing = median(perc_missing) * 100,
    IQR_missing = IQR(perc_missing) * 100,
    min_missing = min(perc_missing) * 100,
    max_missing = max(perc_missing) * 100)

map(c("day", "sex", "exercise"), function(var) { # Differences in missingness
  samples$missingness_bedrest %>%
    inner_join(participants, by = "participant_id") %>%
    glmer(
      formula = as.formula(glue("present ~ {var} + (1 | participant_id)")),
      data = ., family = binomial) %>%
    summary()
})

### Family relative abundance ----

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


## Section 3.2. Beta diversity analysis ----

beta_diversity <- list()

beta_diversity$matrices <- # Load distance matrices
  list.files(
    "data/16S/diversity", pattern = "dm.*\\.tsv$", full.names = TRUE) %>%
  set_names(basename(.) %>% 
    str_remove("\\.tsv$") %>% str_remove("_dm")) %>%
  map(~ .x %>% read_tsv %>% rename("sample_id" = 1))

# Calculate association with participant ID, sex, cohort and exercise group
expand_grid( 
  names(beta_diversity$matrices), 
  c("participant_id", "cohort", "sex", "exercise")) %>%
  pmap(~ {
    matrix <- ..1
    variable <- ..2
    
    data <- beta_diversity$matrices[[matrix]] %>% # Append metadata
      parse_sample_metadata() %>% inner_join(participants)
    
    samples <- data$sample_id # Extract sample IDs
    
    dist_mat <- data %>% # Format matrix as distances
      select(all_of(samples)) %>%
      as.dist()
    
    vars <- data %>% select(variable, participant_id)
    
    formula <- as.formula(paste("dist_mat ~", variable))
    
    set.seed(509)
    adonis2(formula, data = vars, by = "margin") %>%
      slice(1) %>%
      mutate(matrix = {{matrix}}, variable = {{variable}})
      
  }) %>% reduce(rbind) %>%
  write_tsv("results/permanovas_study_variables.tsv")

# Calculate association with FI-36
map(beta_diversity$matrices, ~ {

  data <- .x %>%
    parse_sample_metadata() %>%
    inner_join(frailty %>% filter(!is.na(fi36)))
    # filter(fi36 < 0.6) # Uncomment to remove frailty outlier for analysis

  samples <- data$sample_id
  
  dist_mat <- data %>%
    select(all_of(samples)) %>%
    as.dist()
    
  vars <- data %>% select(fi36, participant_id)
  
  set.seed(0509)
  adonis2(
    dist_mat ~ fi36 + participant_id, 
    data = vars, by = "margin")
})
  
# Frailty component analysis
expand_grid( 
  names(beta_diversity$matrices), 
  colnames(frailty %>% select(contains("Bother_")))) %>%
  pmap(~ {
    matrix <- ..1
    component <- ..2
    
    data = beta_diversity$matrices[[matrix]] %>%
      parse_sample_metadata() %>%
      inner_join(frailty %>% filter(!is.na(!!sym(component))))
    
    dist_mat <- data %>%
      select(all_of(.$sample_id)) %>%
      as.dist()
    
    vars <- data %>% select(all_of(component)) 
    
    n_component <- sum(vars[[component]])
    
    formula <- as.formula(paste("dist_mat ~ ", component))
    set.seed(0509)
    
    adonis2(formula, data = vars) %>%
      slice(1) %>% 
      mutate(
        matrix = {{matrix}},
        component = str_remove(component, "Bother_"),
        n = n_component)
  }) %>% reduce(rbind) %>%
  group_by(matrix) %>%
  adjust_pvalue(p.col = "Pr(>F)", output.col = "p_BH", method = "BH") %>%
  write_tsv("results/permanovas_frailty_components.tsv")

# Alpha diversity analysis ----

# Load alpha diversity data
alpha_diversity <- 
  c("observed_features", "shannon", "faith_pd") %>% # Load all metrics 
  map_dfr(~ read_table(
    glue("data/16S/diversity/{.x}.tsv"),
    skip = 1, 
    col_types = "cd",
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
    read_tsv("data/16S/collapse/table_gut_genus.tsv", skip = 1),
  oral = # Oral genus count feature table
    read_tsv("data/16S/collapse/table_oral_genus.tsv", skip = 1),
  pathway = read_tsv( # Pathway count feature table
    "data/MGS/humann3/pathabundance_unstratified.tsv") %>%
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
