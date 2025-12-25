# Load required libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(gridExtra)
library(reshape2)
library(tidyverse)
library(tibble)
# library(fgsea)
# library(org.Mm.eg.db)
# getwd()
setwd("your directory")
# Load data exactly as is
df <- read_csv("AD_GSEA_p1-3_PDpathways_Notes_with_brain_regions.csv")

# Store original dimensions for comparison
original_rows <- nrow(df)
original_cols <- ncol(df)

# Apply your filter
df_clean <- df %>% 
  filter(!is.na(pval), !is.na(padj), !is.na(NES))

# COMPREHENSIVE INTEGRITY CHECK
cat("DATA INTEGRITY CHECK AFTER FILTERING\n")
cat("=====================================\n")

# 1. Dimension Analysis
cat("DIMENSION ANALYSIS:\n")
cat("  Original dimensions:", original_rows, "x", original_cols, "\n")
cat("  Filtered dimensions:", nrow(df_clean), "x", ncol(df_clean), "\n")
cat("  Rows removed:", original_rows - nrow(df_clean), "\n")

# 2. What was filtered out
cat("FILTERING IMPACT:\n")
na_pval <- sum(is.na(df$pval))
na_padj <- sum(is.na(df$padj))
na_nes <- sum(is.na(df$NES))
cat("  Original rows with NA pval:", na_pval, "\n")
cat("  Original rows with NA padj:", na_padj, "\n")
cat("  Original rows with NA NES:", na_nes, "\n")

# Count rows with ANY NA in these columns
any_na <- sum(is.na(df$pval) | is.na(df$padj) | is.na(df$NES))
cat("  Total rows removed:", any_na, "\n\n")

# 3. Verify filtering worked
cat("FILTERING VERIFICATION:\n")
cat("  No NA in pval:", sum(is.na(df_clean$pval)) == 0, "\n")
cat("  No NA in padj:", sum(is.na(df_clean$padj)) == 0, "\n")
cat("  No NA in NES:", sum(is.na(df_clean$NES)) == 0, "\n\n")

# 4. Data quality checks
cat("DATA QUALITY:\n")
cat("  pval range:", round(min(df_clean$pval), 8), "to", round(max(df_clean$pval), 3), "\n")
cat("  padj range:", round(min(df_clean$padj), 8), "to", round(max(df_clean$padj), 3), "\n")
cat("  NES range:", round(min(df_clean$NES), 3), "to", round(max(df_clean$NES), 3), "\n")
cat("  Unique pathways:", length(unique(df_clean$pathway)), "\n")
cat("  Unique comparisons:", length(unique(df_clean$Comparison)), "\n\n")

# 5. Column completeness check
cat("COLUMN COMPLETENESS:\n")
for(col in colnames(df_clean)) {
  na_count <- sum(is.na(df_clean[[col]]))
  completeness <- round(100 * (nrow(df_clean) - na_count) / nrow(df_clean), 1)
  cat(sprintf("  %s: %d NA (%.1f%% complete)\n", col, na_count, completeness))
}

# 6. Data validation
cat("\nDATA VALIDATION:\n")
if(any(df_clean$pval < 0 | df_clean$pval > 1)) {
  cat("  ⚠️  WARNING: Invalid pval values\n")
} else {
  cat("  ✓ pval values valid [0,1]\n")
}

if(any(df_clean$padj < 0 | df_clean$padj > 1)) {
  cat("  ⚠️  WARNING: Invalid padj values\n")
} else {
  cat("  ✓ padj values valid [0,1]\n")
}

if(any(is.infinite(df_clean$NES))) {
  cat("  ⚠️  WARNING: Infinite NES values\n")
} else {
  cat("  ✓ No infinite NES values\n")
}

# Display basic information about the dataset
cat("Dataset dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
cat("Column names:\n")
print(colnames(df))

# Check pathway column
cat("\nPathway analysis:\n")
cat("Total pathways:", length(df$pathway), "\n")
cat("Unique pathways:", length(unique(df$pathway)), "\n")
cat("Duplicate pathways:", sum(duplicated(df$pathway)), "\n")

# Create clean dataset
key_columns <- c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge", "Comparison", "Matched_Brain_Region")
existing_columns <- key_columns[key_columns %in% colnames(df_clean)]
df_clean <- df_clean[, existing_columns]


cat("DIRECT PATHWAY ENRICHMENT PLOTS\n")
cat("Using Original NES, Pathway, and Comparison Data\n")
cat("================================================================================\n")

cat("Dataset Overview:\n")
cat("• Total entries:", nrow(df_clean), "\n")
cat("• Unique pathways:", length(unique(df_clean$pathway)), "\n")
cat("• Unique comparisons:", length(unique(df_clean$Comparison)), "\n")
cat("• Unique matched brain regions:", length(unique(df_clean$Matched_Brain_Region)), "\n")
cat("• NES range:", round(min(df_clean$NES), 3), "to", round(max(df_clean$NES), 3), "\n\n")

# GsEA plot for brain regions
# Calculate total regions and pathways 
total_regions <- n_distinct(df_clean$Matched_Brain_Region)

# Combined analysis: regional statistics + frequency-based 
pathway_combined <- df_clean %>%
  filter(padj < 0.05) %>%
  
  # First: Calculate pathway statistics by region
  group_by(Matched_Brain_Region, pathway) %>%
  summarise(
    avg_NES = mean(NES, na.rm = TRUE),
    sd_NES = sd(NES, na.rm = TRUE),
    avg_padj = mean(padj, na.rm = TRUE),
    n_comparisons = n(), # How many times this pathway appears in this region (across different time points/conditions)
    .groups = 'drop'
  ) %>%
  
  # Second: Calculate frequency statistics per pathway
  group_by(pathway) %>%
  mutate(
    n_regions = n_distinct(Matched_Brain_Region), # How many different brain regions show this pathway
    n_occurrences = sum(n_comparisons), #Total number of times this pathway appears across all regions
    avg_NES = avg_NES
  ) %>%
  ungroup()

# Get top 20 up and top 20 down per region (No fisher, just GSEA statistics)
top_pathways_by_region <- pathway_combined %>%
  filter(avg_padj < 0.05) %>%
  group_by(Matched_Brain_Region) %>%
  mutate(
    rank_up = ifelse(avg_NES > 1.5, rank(-avg_NES), NA),
    rank_down = ifelse(avg_NES < -1.5, rank(avg_NES), NA)
  ) %>%
  filter((rank_up <= 20 & !is.na(rank_up)) | (rank_down <= 20 & !is.na(rank_down))) %>%
  ungroup() %>%
  mutate(
    pathway_short = str_trunc(pathway, width = 100),
    direction = ifelse(avg_NES > 0, "Positive", "Negative"),
    significance = case_when(
      avg_padj < 0.001 ~ "***",
      avg_padj < 0.01 ~ "**",
      avg_padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Create output directory for individual plots
dir.create("pathway_plots_by_Brain_region_GSEA", showWarnings = FALSE)

# Get list of unique comparisons (NOT regions)
unique_comparisons <- unique(top_pathways_by_region$Matched_Brain_Region)

# Loop through each comparison and create individual plots
for (current_comparison in unique_comparisons) {
  
  # Filter data for this comparison
  comparison_data <- top_pathways_by_region %>%
    filter(Matched_Brain_Region == current_comparison)
  
  # Skip if no data
  if (nrow(comparison_data) == 0) next
  
  # Create plot for this comparison
  p <- ggplot(comparison_data, 
              aes(x = NES, y = reorder(pathway_short, avg_NES))) +
    geom_segment(aes(x = 0, xend = avg_NES, 
                     y = pathway_short, yend = pathway_short,
                     color = direction),
                 linewidth = 1.5, alpha = 0.7) +
    geom_text(aes(label = significance, 
                  x = avg_NES + ifelse(avg_NES > 0, 0.15, -0.15)),
              size = 4, fontface = "bold", 
              hjust = ifelse(comparison_data$avg_NES > 0, 0, 1)) +
    geom_vline(xintercept = 0, linetype = "solid", 
               color = "gray30", linewidth = 1) +
    scale_color_manual(values = c("Positive" = "indianred", 
                                  "Negative" = "steelblue")) +
    labs(
      title = paste("Top Pathways in", current_comparison),
      subtitle = "avg_padj: *** p < 0.001, ** p < 0.01, * p < 0.05",
      x = "Average Normalized Enrichment Score (NES)",
      y = "Pathway",
      color = "Direction"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Create safe filename (remove special characters)
  safe_filename <- gsub("[^A-Za-z0-9_]", "_", current_comparison)
  filename <- paste0("pathway_plots_by_Brain_region_GSEA/", safe_filename, "_top_pathways.png")
  
  # Save plot
  ggsave(filename, p, width = 12, height = 10, dpi = 300)
  
  cat("Saved:", filename, "\n")
}

# GSEA Plot for Comparison
# Get top 20 up and top 20 down per comparison (No fisher, just GSEA statistics)
top_pathways_by_comparison <- df_clean %>%
  filter(padj < 0.05) %>%
  group_by(Comparison) %>%
  mutate(
    rank_up = ifelse(NES > 1.5, rank(-NES), NA),
    rank_down = ifelse(NES < -1.5, rank(NES), NA)
  ) %>%
  filter((rank_up <= 20 & !is.na(rank_up)) | (rank_down <= 20 & !is.na(rank_down))) %>%
  ungroup() %>%
  mutate(
    pathway_short = str_trunc(pathway, width = 100),
    direction = ifelse(NES > 0, "Positive", "Negative"),
    significance = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

# Create output directory for individual plots
dir.create("pathway_plots_by_comparison_GSEA", showWarnings = FALSE)

# Get list of unique comparisons (NOT regions)
unique_comparisons <- unique(top_pathways_by_comparison$Comparison)

# Loop through each comparison and create individual plots
for (current_comparison in unique_comparisons) {
  
  # Filter data for this comparison
  comparison_data <- top_pathways_by_comparison %>%
    filter(Comparison == current_comparison)
  
  # Skip if no data
  if (nrow(comparison_data) == 0) next
  
  # Create plot for this comparison
  p <- ggplot(comparison_data, 
              aes(x = NES, y = reorder(pathway_short, NES))) +
    geom_segment(aes(x = 0, xend = NES, 
                     y = pathway_short, yend = pathway_short,
                     color = direction),
                 linewidth = 1.5, alpha = 0.7) +
    geom_text(aes(label = significance, 
                  x = NES + ifelse(NES > 0, 0.15, -0.15)),
              size = 4, fontface = "bold", 
              hjust = ifelse(comparison_data$NES > 0, 0, 1)) +
    geom_vline(xintercept = 0, linetype = "solid", 
               color = "gray30", linewidth = 1) +
    scale_color_manual(values = c("Positive" = "indianred", 
                                  "Negative" = "steelblue")) +
    labs(
      title = paste("Top Pathways in", current_comparison),
      subtitle = "padj: *** p < 0.001, ** p < 0.01, * p < 0.05",
      x = "Normalized Enrichment Score (NES)",
      y = "Pathway",
      color = "Direction"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  # Create safe filename (remove special characters)
  safe_filename <- gsub("[^A-Za-z0-9_]", "_", current_comparison)
  filename <- paste0("pathway_plots_by_comparison_GSEA/", safe_filename, "_top_pathways.png")
  
  # Save plot
  ggsave(filename, p, width = 12, height = 10, dpi = 300)
  
  cat("Saved:", filename, "\n")
}


