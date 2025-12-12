# Setup ----

# Load packages
library(tidyverse)
library(scales)
library(knitr)
library(ggpubr)
library(ggrepel)
library(phyloseq)
library(microbiome)
library(patchwork)
library(cowplot)
instlibrary(ggh4x)

# Initialize lists
figures <- list()
figures$fig1 <- list()
figures$fig2 <- list()
figures$fig3 <- list()
figures$fig4 <- list()
figures$fig5 <- list()
figures$figS4 <- list()
figures$figS5 <- list()

# Define functions
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

# Assign unique timepoint to visualize data across all study phases
assign_unique_time <- function(data_table) { 
  data_table %>%
    mutate( # Assign unique timepoint for plotting
      timepoint = case_when(           
        phase == "baseline" ~ 0,           
        phase == "bedrest" ~ day,           
        phase == "recovery" ~ 14 + day,           
        phase == "4w" ~ 23,           
        phase == "4m" ~ 25))
}

# Add vertical lines to separate study phases
add_phase_lines <- function() { 
  
  list(
  geom_vline(
    xintercept = c(1, 14, 21), # Start of HDBR, end of HDBR, end of REC
    linetype = "dotted", colour = "gray20", linewidth = 0.3),
  scale_x_continuous(
    breaks = c(1, 14, 21),
    labels = c("HDT-1", "HDT-14", "REC-7")))
}

# Figures S1, S2 and S3 ----

plot_sample_summary <- function(data) {
  
  data %>%
    parse_sample_metadata() %>%
    mutate(time_id = str_extract(sample_id, "(?<=-).*")) %>% 
    complete(time_id, participant_id) %>%
    mutate(
      phase = case_when(
        str_detect(time_id, "BDC") ~ "baseline",
        str_detect(time_id, "HDT") ~ "bedrest",
        str_detect(time_id, "R") ~ "recovery",
        str_detect(time_id, "W4") ~ "4w",
        str_detect(time_id, "M4") ~ "4m")) %>%
    filter(!is.na(participant_id)) %>%
    left_join(participants) %>%
    
    ggplot(aes(x = participant_id, y = time_id)) + 
    facet_grid(
      fct_relevel(phase, "baseline", "bedrest", "recovery", "4w", "4m") ~ ., 
      scales = "free", space = "free") +
    xlab(NULL) + ylab(NULL) + 
    scale_x_discrete(position = "top") +
    scale_y_discrete(limits = rev) + 
    geom_tile(aes(fill = seq_depth), colour = "white", width = 0.9) +
    scale_fill_viridis_c(
      option = "viridis", na.value = "grey70",
      trans = "log10", labels = label_number(),
      name = "Sequencing depth (total reads)",
      guide = guide_colorbar(
        title.position = "top", 
        title.hjust = 0.5)) +
    
    theme_bw() + 
    theme(
      strip.text = element_blank(), 
      strip.background = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = -0.1),
      legend.position = "bottom",
      legend.key.width = unit(1.8, "cm"),
      legend.key.height = unit(0.3, "cm"))
}

# Figure S1
figures$figS1 <- # Plot summary of 16S gut samples
  read_csv("data/16S/sample_depth_gut.csv") %>%
  plot_sample_summary()

ggsave(
  "results/figures/FigS1.png", figures$figS1,
  width = 6.5, height = 7)

# Figure S2
figures$figS2 <- # Plot summary of 16S gut samples
  read_csv("data/16S/sample_depth_oral.csv") %>%
  plot_sample_summary()

ggsave(
  "results/figures/FigS2.png", figures$figS1,
  width = 6.5, height = 7)

# Figure S3
figures$figS3 <- # Plot summary of fecal metagenomic samples
  read_tsv("data/MGS/kneaddata/kneaddata_read_counts.txt") %>%
  rename(`sample_id` = `Sample`) %>%
  rename(`seq_depth` = `trimmed pair1`) %>%
  mutate(sample_id = sub("_.*", "", sample_id)) %>%
  mutate(sample_id = str_replace(sample_id, "HDT([1-9])$", "HDT0\\1")) %>%
  select(sample_id, seq_depth) %>%
  plot_sample_summary() + 
  geom_point( # Add extra geom for metabolomic samples 
    data = feature_tables$metabolite %>%
      select(-1) %>%
      pivot_longer(everything(), names_to = "sample_id") %>%
      parse_sample_metadata() %>%
      mutate(time_id = str_extract(sample_id, "(?<=-).*")),
    aes(x = participant_id, y = time_id),
    colour = "white",
    size = 2)

ggsave(
  "results/figures/FigS3.png", figures$figS3,
  width = 6.5, height = 4)


# Figure 1 ----

figures$fig1$a <-
  
  feature_tables$family_relab %>%
    right_join(read_tsv("results/top_6_families.tsv")) %>%
    arrange(total_abundance) %>%
    mutate(
      sample_id = sample_id %>% str_remove("[GO]"), # No blank rows in figure
      sample_type = case_when( # For facet grid labels
        sample_type == "gut" ~ "Gut",
        sample_type == "oral" ~ "Oral"),
      feature_name = # Format family strings
      str_extract(feature_name, "f__.+") %>% str_remove("f__"),
      participant_id = case_when( # Add exercise group to participant ID
       exercise == 0 ~ paste0("C - ", participant_id),
       exercise == 1 ~ paste0("E - ", participant_id)),
      participant_id = case_when( # Add sex to participant ID
        sex == 0 ~ paste0(participant_id, "F"),
        sex == 1 ~ paste0(participant_id, "M"))) %>%
      assign_unique_time() %>%
  
  ggplot(
    mapping = aes(x = relab, y = sample_id)) +
    facet_grid2(
      participant_id ~ sample_type, 
      scales = "free_y", switch = "y") +
    xlab(NULL) + ylab(NULL) +
    geom_bar(
      aes(fill = feature_name), stat = "identity", 
      width = 1, alpha = 0.8) +
    scale_fill_brewer(palette = "Paired") +
    guides(fill = guide_legend(nrow = 6)) +
    labs(tag = "a.") + 
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.3, "cm"),
      strip.placement = "inside",
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.spacing = unit(0, "lines"),
      panel.border = element_rect(fill = NA, colour = "black"),
      # strip.placement = "outside",
      # legend.margin = margin(t=0, unit = "pt"),
      strip.background = element_rect(fill = "gray95"),
      strip.text.x = element_text(size = 11),
      strip.text.y.left = element_text(size = 8, angle = 0),
      strip.clip = "on",
      axis.text.x = element_blank())

figures$fig1$b <- 
  ggplot(
    data = 
      read_tsv(
        "data/16S/diversity/ordination_gut.txt", 
        skip = 7, col_names = c("sample_id", "PC1", "PC2")) %>%
      filter(str_detect(sample_id, "PT")) %>%
      mutate(`PC2` = str_remove(`PC2`, "\t.*") %>% as.numeric()) %>%
      left_join(samples),
    mapping = aes(x = PC1, y = PC2)) +
  xlab("PC1 (68.1%)") + ylab("PC2 (27.0%)") +
  geom_point(aes(colour = participant_id), size = 0.6) + 
  geom_label_repel(
    data = . %>% group_by(participant_id) %>% slice(1),
    aes(colour = participant_id, group = 1, label = participant_id),
    size = 3, segment.color = "transparent", fontface = "bold") +
  guides(colour = "none") +
  annotate(
    "text", label = "R\u00B2 = 0.953\np < 0.001", 
    x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, 
    size = 4, fontface = "bold") +
  labs(tag = "b.") + 
  theme_bw() + 
  theme(
    legend.position = "none",
    aspect.ratio = 1)

figures$fig1$c <- 
  ggplot(
    data = 
      read_tsv(
        "data/16S/diversity/ordination_oral.txt", 
        skip = 7, col_names = c("sample_id", "PC1", "PC2")) %>%
      filter(str_detect(sample_id, "PT")) %>%
      mutate(`PC2` = str_remove(`PC2`, "\t.*") %>% as.numeric()) %>%
      left_join(samples),
    mapping = aes(x = PC1, y = PC2)) +
  
  xlab("PC1 (54.3%)") + ylab("PC2 (33.8%)") +
  geom_point(aes(colour = participant_id), size = 0.6) + 
  geom_label_repel(
    data = . %>% group_by(participant_id) %>% slice(1),
    aes(colour = participant_id, group = 1, label = participant_id),
    size = 3, segment.color = "transparent", fontface = "bold") +
  guides(colour = "none") +
  annotate(
    "text", label = "R\u00B2 = 0.845\np < 0.001", 
    x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, 
    size = 4, fontface = "bold") +
  labs(tag = "c.") + 
  theme_bw() + 
  theme(
    legend.position = "none",
    aspect.ratio = 1)

figures$fig1$final <-
  (((figures$fig1$a / guide_area()) + 
     plot_layout(heights = c(1000, 1))) |
  ((figures$fig1$b / figures$fig1$c) +
     plot_layout(heights = c(1, 1)))) +
  theme(
    text = element_text(family = "Arial"),
    plot.margin = margin(0, 0, 0, 0),
    plot.tag.position = c(0, 1),  
    # plot.tag.location = "panel",
    plot.tag = element_text(margin = margin(b = -5)))  

ggsave( 
  filename = "results/figures/Fig1.png", plot = figures$fig1$final, 
  width = 6.8, height = 7.5, units = "in", dpi = 600) 

# Figure 2 ----

figures$fig2$ab <- 
  list(
    beta_diversity$pcoa_gut %>%
      ggplot(mapping = aes(x = PC1, y = PC2)) + 
      xlab(NULL) + ylab(NULL) +
      labs(tag = "a.") + 
      geom_point(aes(colour = cohort), size = 0.1, alpha = 0.6) +
      stat_ellipse(aes(colour = cohort), alpha = 0.8) +
      guides(
        colour = guide_legend(override.aes = list(
          shape = 16, size = 1.8, alpha = 1))) +
      annotate(
        "text", label = "R\u00B2 = 0.124\np <  0.001", 
        x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, 
        size = 2.8, fontface = "bold"),
    beta_diversity$pcoa_gut %>%
      ggplot(mapping = aes(x = PC1, y = PC2)) + 
      xlab(NULL) + ylab(NULL) +
      geom_point(aes(colour = sex), size = 0.1, alpha = 0.6) + 
      stat_ellipse(aes(colour = sex), alpha = 0.8) +
      scale_colour_manual(
        labels = c("Female", "Male"), 
        values = c("#b07aa1", "#59a14f")) + 
      guides(
        colour = guide_legend(override.aes = list(
          shape = 16, size = 1.8, alpha = 1))) +
      annotate(
        "text", label = "R\u00B2 = 0.0957\np <  0.001", 
        x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, 
        size = 2.8, fontface = "bold"),
    beta_diversity$pcoa_gut %>%
      ggplot(mapping = aes(x = PC1, y = PC2)) + 
      xlab(NULL) + ylab(NULL) +
      geom_point(aes(colour = exercise), size = 0.1, alpha = 0.6) +
      stat_ellipse(aes(colour = exercise), alpha = 0.8) +
      scale_colour_manual(
        labels = c("Control", "Exercise"), 
        values = c("#f28e2b", "#4e79a7")) + 
      guides(
        colour = guide_legend(override.aes = list(
          shape = 16, size = 1.8, alpha = 1))) +
      annotate(
        "text", label = "R\u00B2 = 0.0667\np <  0.001", 
        x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, 
        size = 2.8, fontface = "bold"),
    beta_diversity$pcoa_oral %>%
      ggplot(mapping = aes(x = PC1, y = PC2)) + 
      xlab(NULL) + ylab(NULL) +
      labs(tag = "b.") + 
      geom_point(aes(colour = cohort), size = 0.1, alpha = 0.6) + 
      stat_ellipse(aes(colour = cohort), alpha = 0.8) + 
      guides(
        colour = guide_legend(override.aes = list(
          shape = 16, size = 1.8, alpha = 1))) +
      annotate(
        "text", label = "R\u00B2 = 0.155\np <  0.001", 
        x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, 
        size = 2.8, fontface = "bold"),
    beta_diversity$pcoa_oral %>%
      ggplot(mapping = aes(x = PC1, y = PC2)) + 
      xlab(NULL) + ylab(NULL) +
      geom_point(aes(colour = sex), size = 0.1, alpha = 0.6) + 
      stat_ellipse(aes(colour = sex), alpha = 0.8) +
      scale_colour_manual(
        labels = c("Female", "Male"), 
        values = c("#b07aa1", "#59a14f")) +
      guides(
        colour = guide_legend(override.aes = list(
          shape = 16, size = 1.8, alpha = 1))) +
      annotate(
        "text", label = "R\u00B2 = 0.124\np <  0.001", 
        x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, 
        size = 2.8, fontface = "bold"),
    beta_diversity$pcoa_oral %>%
      ggplot(mapping = aes(x = PC1, y = PC2)) + 
      xlab(NULL) + ylab(NULL) +
      geom_point(aes(colour = exercise), size = 0.1, alpha = 0.6) +
      stat_ellipse(aes(colour = exercise), alpha = 0.8) +
      scale_colour_manual(
        labels = c("Control", "Exercise"), 
        values = c("#f28e2b", "#4e79a7")) +
      guides(
        colour = guide_legend(override.aes = list(
          shape = 16, size = 1.8, alpha = 1))) +
      annotate(
        "text", label = "R\u00B2 = 0.0185\np =  0.005", 
        x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, 
        size = 2.8, fontface = "bold"))

figures$fig2$final <- 
  wrap_plots(figures$fig2$ab, nrow = 2, ncol = 3) +
  plot_layout(guides = "collect") &
  theme_bw() &
  theme(
    aspect.ratio = 1,
    legend.position = "right",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank())

ggsave(
  filename = "results/figures/Fig2.png", 
  plot = figures$fig2$final,
  width = 5.2, height = 3, units = "in", dpi = 600) 


# Figure 3 ----

figures$fig3$a <- 
  
  # Figure data
  beta_diversity$pcoa_gut %>% # Load PCoA results
    
    inner_join(frailty) %>% # Select FI-36
    mutate( # Format data
      exercise = case_when(
        exercise == 0 ~ "C", exercise == 1 ~ "E"),
      sex = case_when(
        sex == 0 ~ "F", sex == 1 ~ "M"),
      subject_id = str_remove(participant_id, "PT")) %>%
  
  # Generate plot
  ggplot(mapping = aes(x = PC1, y = PC2)) + 
    xlab("PC1 (68.1%)") + ylab("PC2 (27.0%)") +
    labs(tag = "a.", colour = "FI-36") +
    # Normal points (FI36 < 0.45)
    geom_point(
      data = . %>% filter(fi36 <= 0.45),
      aes(colour = fi36), size = 1) + 
    # Outliers (FI36 > 0.6)
    geom_point(
      data = . %>% filter(fi36 > 0.6),
      color = "firebrick4", size = 2) +
    # PERMANOVA results
    annotate(
      "text", label = "R\u00B2 = 0.000228\np =  0.588", 
      x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 3.2, fontface = "bold")

figures$fig3$b <- 
  
  beta_diversity$pcoa_oral %>%
    inner_join(frailty) %>%
    mutate(
      exercise = case_when(
        exercise == 0 ~ "C", exercise == 1 ~ "E"),
      sex = case_when(
        sex == 0 ~ "F", sex == 1 ~ "M"),
      subject_id = str_remove(participant_id, "PT")) %>%
  
  ggplot(mapping = aes(x = PC1, y = PC2)) + 
    xlab("PC1 (54.3%)") + ylab("PC2 (33.8%)") +
    labs(tag = "b.", colour = "FI-36") +
    
    # Normal points (FI36 < 0.45)
    geom_point(
      data = . %>% filter(fi36 <= 0.45),
      aes(colour = fi36), size = 1) + 
  
    # Outliers (FI36 > 0.6)
  geom_point(
    data = . %>% filter(fi36 > 0.6),
    color = "firebrick4", size = 2) +
    
    annotate(
      "text", label = "R\u00B2 = 0.00329\np = 0.086", 
      x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 3.2, fontface = "bold") +
    
    theme_bw() + 
    theme(
    )

figures$fig3$final <-
  
  # Define plot layout
  figures$fig3$a + figures$fig3$b + 
  plot_layout(guides = "collect") &
  
  # Format legend
  scale_colour_gradient(
    low = "gray80", high = "red", 
    limits = c(0, 0.45), breaks = c(0, 0.2, 0.4)) &
  # Format axes
  scale_x_continuous(breaks = c(-0.1, 0, 0.1, 0.2)) &
  scale_y_continuous(breaks = c(-0.1, 0, 0.1, 0.2)) &
  # Set figure theme
  theme_bw() &
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    legend.position = "right",
    legend.title = element_text("FI-36", size = 11),
    legend.justification = "center",
    legend.key.width = unit(0.5, "cm"),
    legend.key.height = unit(0.65, "cm"),
    plot.tag = element_text(size = 14),
    aspect.ratio = 1
  )

# Save figure
ggsave( 
  filename = "results/figures/Fig3.png", plot = figures$fig3$final, 
  width = 6, height = 2.35, units = "in", dpi = 600) 


# Figure 4 ----

figures$fig4$a <- 
  
  alpha_diversity %>% 
    # filter(sample_id)
    parse_sample_metadata() %>%
  mutate(
    sample_type = case_when(
        str_detect(sample_id, "G") ~  "gut",
        str_detect(sample_id, "O") ~  "oral")) %>%
  dplyr::filter(metric == "faith", sample_type == "gut") %>%
    dplyr::filter(phase %in% c("bedrest", "recovery")) %>%
  
    # filter(metric == "faith", phase %in% c("bedrest", "recovery")) %>%
    assign_unique_time() %>%
  
  left_join(participants) %>%
  
  ggplot(aes(x = timepoint, y = value)) + 
  # facet_grid(. ~ sample_type) + 
    xlab(NULL) + ylab("Phylogenetic diversity") +
  labs(tag = "a.") +
    add_phase_lines() +
  geom_line(
      aes(group = `participant_id`, colour = as.factor(exercise)),
      linewidth = 0.4, alpha = 0.3,  key_glyph = "blank") +
    geom_smooth(
      # data = raw_data %>% 
      #   dplyr::filter(phase != "baseline"),
      aes(colour = as.factor(exercise)),
      method = "lm", se = FALSE, linewidth = 1,
      formula = 
        y ~ 
        x + I((x - 14) * ( x > 14)) + 
        I((x - 21) * (x > 21)) + 
        I((x - 23) * (x > 23))) +
  geom_point( # Dummy geom for legend
    aes(colour = as.factor(exercise)),
    shape = 15, size = 4, alpha = 0) +
  scale_colour_manual(
    values = c("#f28e2b", "#4e79a7"),
    # breaks = c("C", "E"),
    labels = c("Control", "Exercise")) +
  guides(
    colour = guide_legend(override.aes = list(
    shape = 15, size = 3, linewidth = 0, alpha = 1))) +
  annotate(
      "text", label = "p (control) =  0.028; p (exercise) = 0.86", 
      x = -Inf, y = Inf, hjust = -0.08, vjust = 2, size = 3, fontface = "bold") +
  theme_bw() + 
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      strip.background = element_rect(fill = "gray95"),
      strip.text.x = element_text(size = 12),
      strip.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 8))

figures$fig4$b <- 
  
  alpha_diversity %>% 
    # filter(sample_id)
    parse_sample_metadata() %>%
  mutate(
    sample_type = case_when(
        str_detect(sample_id, "G") ~  "gut",
        str_detect(sample_id, "O") ~  "oral")) %>%
  dplyr::filter(metric == "faith", sample_type == "oral") %>%
    dplyr::filter(phase %in% c("bedrest", "recovery")) %>%
  
    # filter(metric == "faith", phase %in% c("bedrest", "recovery")) %>%
    assign_unique_time() %>%
  
  left_join(participants) %>%
  
  ggplot(aes(x = timepoint, y = value)) + 
  # facet_grid(. ~ sample_type) + 
    xlab(NULL) + ylab("Phylogenetic diversity") +
  labs(tag = "b.") +
    add_phase_lines() +
  geom_line(
      aes(group = `participant_id`, colour = as.factor(exercise)),
      linewidth = 0.4, alpha = 0.3) +
    geom_smooth(
      # data = raw_data %>% 
      #   dplyr::filter(phase != "baseline"),
      aes(colour = as.factor(exercise)),
      method = "lm", se = FALSE, linewidth = 1,
      formula = 
        y ~ 
        x + I((x - 14) * ( x > 14)) + 
        I((x - 21) * (x > 21)) + 
        I((x - 23) * (x > 23))) +
  geom_point( # Dummy geom for legend
    aes(colour = as.factor(exercise)),
    shape = 15, size = 4, alpha = 0) +
  scale_colour_manual(
    values = c("#f28e2b", "#4e79a7"),
    # breaks = c("C", "E"),
    labels = c("Control", "Exercise")) +
  guides(
    colour = guide_legend(override.aes = list(
      shape = 15, size = 3, linewidth = 0, alpha = 1))) +
   annotate(
      "text", label = "p (control) =  0.54; p (exercise) = 0.33", 
      x = -Inf, y = Inf, hjust = -0.08, vjust = 2, size = 3, fontface = "bold") +
  theme_bw() + 
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      strip.background = element_rect(fill = "gray95"),
      strip.text.x = element_text(size = 12),
      strip.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 8))

figures$fig4$final <-
  
  # Define plot layout
  figures$fig4$a + figures$fig4$b &
  plot_layout(guides = "collect") +
  
  # Format legend
  
  # Format axes

  # Set figure theme
  theme_bw() &
  theme(
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11),
    legend.position = "bottom",
    # legend.title = element_text("FI-36", size = 12),
    legend.justification = "center",
    plot.tag = element_text(size = 14)
  )

# Save figure
ggsave( 
  filename = "results/figures/Fig4.png", plot = figures$fig4$final, 
  width = 6.8, height = 3, units = "in", dpi = 600) 


# Figures 5, S4 and S5 ----

# Define functions for generating all differential abundance figures
plot_da_estimates <- function(da_results, p_col) { # Left hand side of plot
  
  ggplot(
    data =
      read_tsv({{ da_results }}) %>%
      left_join(read_csv("~/CAIS-microbiome/data/da_results_recoding.csv")) %>%
      group_by(exercise) %>%
      arrange(Estimate) %>%
      mutate(
        result_number = row_number(),
        exercise = case_when(
          exercise == 0 ~ "Control",
          exercise == 1 ~ "Exercise"),
        direction = case_when(
          as.numeric(Estimate) < 0 ~ "Decrease",
          as.numeric(Estimate) > 0 ~ "Increase"),
        feature_name_number = paste0(feature_name_new, " (", result_number, ")"),
        feature_name_number = reorder(feature_name_number, -Estimate)
      ),
    
    mapping = aes(x = 0, y = feature_name_number)
  ) +
    
    facet_grid(
      exercise ~ ., scales = "free", space = "free", switch = "y"
    ) +
    
    geom_point(
      data = . %>% arrange(Estimate),
      aes(size = abs(Estimate), colour = direction)
    ) +
    
    scale_colour_manual(
      values = c("#76b7b2", "#edc948"),
      labels = c("Decrease", "Increase")
    ) +
    guides(size = "none", colour = guide_legend(override.aes = list(size = 6))) +
    
    geom_text(
      aes(
        x = 0,
        label = case_when(
          as.numeric(.data[[p_col]]) < 0.001 ~ "***",
          as.numeric(.data[[p_col]]) < 0.01  ~ "**",
          as.numeric(.data[[p_col]]) < 0.05  ~ "*"
        )
      ),
      size = 3, colour = "gray20"
    ) +
    
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 20, face = "italic"),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.position = "bottom",
      panel.border = element_blank(),
      strip.placement = "outside",
      strip.background = element_rect(fill = "gray95")
    )
}

plot_feature_trajectories <- function( # Right hand side of plot
    da_results, feature_table, features_control, features_exercise) {
  
  ggplot(

    data = # Load CLR-transformed
      read_tsv({{ da_results }}) %>%
      group_by(exercise) %>%
      arrange(Estimate) %>%
      mutate(result_number = row_number()) %>%
        left_join(
          {{ feature_table }} %>%
            pivot_longer(
              cols = -1, 
              names_to = "sample_id", 
              values_to = "abundance")) %>%
      parse_sample_metadata() %>%
      dplyr::filter(phase %in% c("bedrest", "recovery")) %>%
      assign_unique_time(),
    
    mapping = aes(x = timepoint, y = as.numeric(abundance))) +
    facet_grid(exercise ~ ., scales = "free") +
    force_panelsizes(rows = c(features_control, features_exercise)) +
    xlab(NULL) + ylab("Abundance") +

    scale_y_continuous(
      position = "right",
      breaks = c(-6, -4, -2, 0, 2, 4, 6, 8)) +

    add_phase_lines() +

    geom_smooth(
      aes(colour = `feature_name`),
      method = "lm", se = FALSE, linewidth = 0.2, alpha = 0.7,
      formula =
        y ~
        x +
        I((x - 1) * (x > 1)) +
        I((x - 14) * (x > 14)) +
        I((x - 21) * (x > 21))) +
    geom_text_repel(
      data = . %>%
        group_by(exercise, participant_id, feature_name, result_number) %>%
        slice_tail(n = 1) %>%
        group_by(exercise, Estimate, feature_name, result_number) %>%
        summarise(abundance = mean(as.numeric(abundance))) %>%
        mutate(timepoint = 21),
      aes(x = 21, label = `result_number`, colour = `feature_name`), 
      stat = "unique",
      size = 3, fontface = "bold", alpha = 1) +
    scale_fill_brewer(palette = "Paired") +
    guides(colour = "none") + 

    theme_bw() +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank(),
      axis.text.x = element_text(size= 8),
      axis.title.y = element_text(size = 10, colour = "gray20"),
      axis.text.y = element_text(size = 8))
}

ggsave( 
  filename = "results/figures/Fig5.png",
  plot = 
    (plot_da_estimates(
      "results/differential_abundance/genus_gut.tsv", "p_BH") + # 16S gut samples
       labs(tag = "a.") +  
    plot_feature_trajectories(
      "results/differential_abundance/genus_gut.tsv", 
      feature_tables$gut_clr, 
      14, 4) + 
    plot_layout(widths = c(1, 2)))
    /
    (plot_da_estimates(
      "results/differential_abundance/genus_oral.tsv", "p_BH") + # 16S oral samples
       labs(tag = "b.") + 
    plot_feature_trajectories(
      "results/differential_abundance/genus_oral.tsv", 
      feature_tables$oral_clr, 
      3, 2)  +
    plot_layout(widths = c(1, 2))) +
    plot_layout(
      guides = "collect", 
      heights = c(3, 1)) 
    &
    theme(
      strip.text = element_text(size = 8),
      legend.position = "bottom",
      legend.justification = "center",
      axis.text.y = element_text(size = 8)),
    
    width = 6, height = 8.5, 
    units = "in", dpi = 600) 

ggsave( 
  filename = "results/figures/FigS4.png",
  plot = 
    plot_da_estimates(
      "results/differential_abundance/pathway_gut.tsv", "Pr(>|t|)") + # 16S gut samples
           plot_feature_trajectories(
             "results/differential_abundance/pathway_gut.tsv", 
             feature_tables$pathway_clr, 
             14, 
             6) &
     plot_layout(
       widths = c(1, 2)) &
           theme(
             axis.text.y = element_text(size = 7),
             plot.margin = margin(0, 0, 0, 0, "cm")),
  width = 6.5, 
  height = 7, 
  units = "in", 
  dpi = 600) 

ggsave( 
  filename = "results/figures/FigS5.png",
  plot = 
    plot_da_estimates("results/differential_abundance/metabolite_gut.tsv", "Pr(>|t|)") + # 16S gut samples
           plot_feature_trajectories(
             "results/differential_abundance/metabolite_gut.tsv", 
             feature_tables$metabolite, 
             16, 
             13) &
     plot_layout(
       widths = c(1, 2)) &
           theme(
             axis.text.y = element_text(size = 7),
             plot.margin = margin(0, 0, 0, 0, "cm")),
  width = 6.5, 
  height = 7, 
  units = "in", 
  dpi = 600)