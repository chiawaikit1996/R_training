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
library(fgsea)
library(org.Mm.eg.db)
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


# # Show some examples
# cat("\nExamples (rowname vs pathway):")
# for(i in 1:5) {
#   match_status <- ifelse(rownames(df)[i] == df$pathway[i], "✓ EXACT MATCH", "✗ NO MATCH")
#   cat(sprintf("\nRow %d: %s", i, match_status))
#   cat(sprintf("\n  Rowname: %s", rownames(df)[i]))
#   cat(sprintf("\n  Pathway: %s", df$pathway[i]))
# }

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

# Method 1: Fisher's Exact Test 
# Tests if a pathway appears in more regions than expected by chance.

library(dplyr)

# Calculate total regions and unique pathways
total_regions <- n_distinct(df_clean$Matched_Brain_Region)
total_unique_pathways <- df_clean %>% 
  filter(padj < 0.01) %>% 
  distinct(pathway) %>% 
  nrow()

# Combined analysis with CORRECTED Fisher's test
pathway_combined_fisher <- df_clean %>%
  filter(padj < 0.01) %>%
  
  # First: Calculate pathway statistics by region
  group_by(Matched_Brain_Region, pathway) %>%
  summarise(
    avg_NES = mean(NES, na.rm = TRUE),
    sd_NES = sd(NES, na.rm = TRUE),
    avg_padj = mean(padj, na.rm = TRUE),
    n_comparisons = n(),
    .groups = 'drop'
  ) %>%
  
  # Second: Calculate frequency statistics per pathway
  group_by(pathway) %>%
  mutate(
    n_regions = n_distinct(Matched_Brain_Region),
    n_occurrences = sum(n_comparisons)
  ) %>%
  ungroup() %>%
  
  # Third: CORRECTED Fisher's exact test
  group_by(pathway) %>%
  mutate(
    fisher_pval = {
      # Contingency table:
      #                In pathway    Not in pathway
      # In regions          a              b
      # Not in regions      c              d
      
      a <- dplyr::first(n_regions)  # regions with this pathway
      b <- total_regions - a         # regions without this pathway
      c <- total_unique_pathways - 1 # other pathways
      d <- max(1, c * total_regions - a)  # ensure positive
      
      fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
    }
  ) %>%
  ungroup() %>%
  
  mutate(
    fisher_padj = p.adjust(fisher_pval, method = "BH")
  ) %>%
  
  arrange(Matched_Brain_Region, desc(avg_NES))


# View results
head(pathway_combined_fisher, 20)

library(dplyr)
library(ggplot2)
library(stringr)

# Step 1: Identify top 20 pathways by cross-region prevalence
top_20_pathways <- pathway_combined_fisher %>%
  group_by(pathway) %>%
  summarise(
    n_regions = dplyr::first(n_regions),
    fisher_padj = dplyr::first(fisher_padj),
    avg_abs_NES = mean(abs(avg_NES)),
    .groups = 'drop'
  ) %>%
  filter(fisher_padj < 0.05) %>%
  arrange(desc(n_regions), fisher_padj, desc(avg_abs_NES)) %>%
  slice_head(n = 20) %>%
  pull(pathway)

# Step 2: Prepare heatmap data with only top 20 pathways
heatmap_data <- pathway_combined_fisher %>%
  filter(pathway %in% top_20_pathways) %>%
  filter(avg_padj < 0.05) %>%  # Only show region-specific significant results
  mutate(
    pathway_short = str_trunc(pathway, width = 100),
    NES_capped = pmax(pmin(avg_NES, 3), -3)  # Cap for visualization
  )

# Step 3: Create ordered factor for proper sorting
pathway_order <- pathway_combined_fisher %>%
  filter(pathway %in% top_20_pathways) %>%
  group_by(pathway) %>%
  summarise(
    avg_NES_global = mean(avg_NES),
    n_regions = dplyr::first(n_regions),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_regions), desc(avg_NES_global)) %>%
  mutate(pathway_short = str_trunc(pathway, width = 100))

heatmap_data <- heatmap_data %>%
  mutate(
    pathway_short = factor(pathway_short, levels = pathway_order$pathway_short)
  )

# Step 4: Create heatmap
ggplot(heatmap_data, 
       aes(x = Matched_Brain_Region, y = pathway_short)) +
  geom_tile(aes(fill = NES_capped), color = "white", size = 0.5) +
  geom_text(aes(label = ifelse(fisher_padj < 0.001, "***",
                               ifelse(fisher_padj < 0.01, "**",
                                      ifelse(fisher_padj < 0.05, "*", "")))),
            size = 3, fontface = "bold") +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "indianred",
    midpoint = 0,
    name = "Avg NES",
    limits = c(-3, 3)
  ) +
  labs(
    title = "Top 20 Cross-Region Pathways: Heatmap View",
    subtitle = "Pathways ranked by cross-region prevalence | *** Fisher padj < 0.001, ** < 0.01, * < 0.05",
    x = "Brain Region",
    y = "Pathway"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave("pathway_heatmap_top20_cross_region_fisher.png", width = 14, height = 11, dpi = 300, bg = "white")

# Print summary
cat("\n=== TOP 20 PATHWAYS SUMMARY ===\n")
print(pathway_order %>% dplyr::select(pathway_short, n_regions, avg_NES_global))



# # Get top 20 up and top 20 down per region
# top_pathways_by_region <- pathway_combined_fisher %>%
#   filter(fisher_padj < 0.05) %>%
#   group_by(Matched_Brain_Region) %>%
#   mutate(
#     rank_up = ifelse(avg_NES > 1.5, rank(-avg_NES), NA),
#     rank_down = ifelse(avg_NES < -1.5, rank(avg_NES), NA)
#   ) %>%
#   filter((rank_up <= 20 & !is.na(rank_up)) | (rank_down <= 20 & !is.na(rank_down))) %>%
#   ungroup() %>%
#   mutate(
#     pathway_short = str_trunc(pathway, width = 100),
#     direction = ifelse(avg_NES > 0, "Positive", "Negative"),
#     significance = case_when(
#       fisher_pval < 0.001 ~ "***",
#       fisher_pval < 0.01 ~ "**",
#       fisher_pval < 0.05 ~ "*",
#       TRUE ~ "ns"
#     )
#     
#   )
# 
# # Create output directory for individual plots
# dir.create("pathway_plots_by_region_fisher", showWarnings = FALSE)
# 
# # Get list of unique brain regions
# regions <- unique(top_pathways_by_region$Matched_Brain_Region)
# 
# # Loop through each region and create individual plots
# for (region in regions) {
#   
#   # Filter data for this region
#   region_data <- top_pathways_by_region %>%
#     filter(Matched_Brain_Region == region)
#   
#   # Skip if no data
#   if (nrow(region_data) == 0) next
#   
#   # Create plot for this region
#   p <- ggplot(region_data, 
#               aes(x = avg_NES, y = reorder(pathway_short, avg_NES))) +
#     geom_segment(aes(x = 0, xend = avg_NES, 
#                      y = pathway_short, yend = pathway_short,
#                      color = direction),
#                  linewidth = 1.5, alpha = 0.7) +
#     # geom_point(aes(size = -log10(fisher_padj), color = direction),
#     #            alpha = 0.8) +
#     geom_text(aes(label = significance, 
#                   x = avg_NES + ifelse(avg_NES > 0, 0.15, -0.15)),
#               size = 4, fontface = "bold", 
#               hjust = ifelse(region_data$avg_NES > 0, 0, 1)) +
#     geom_vline(xintercept = 0, linetype = "solid", 
#                color = "gray30", linewidth = 1) +
#     scale_color_manual(values = c("Positive" = "steelblue", 
#                                   "Negative" = "indianred")) +
#     # scale_size_continuous(name = "-log10(Fisher padj)", range = c(3, 10)) +
#     labs(
#       title = paste("Top Pathways in", region),
#       subtitle = "Fisher padj *** p < 0.001, ** p < 0.01, * p < 0.05",
#       x = "Average Normalized Enrichment Score (NES)",
#       y = "Pathway",
#       color = "Direction"
#     ) +
#     theme_bw(base_size = 12) +
#     theme(
#       plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
#       plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5),
#       axis.text.y = element_text(size = 9),
#       legend.position = "right",
#       panel.grid.major.y = element_line(color = "gray90"),
#       panel.grid.minor = element_blank(),
#       panel.background = element_rect(fill = "white", color = NA),
#       plot.background = element_rect(fill = "white", color = NA)
#     )
#   
#   # Create safe filename (remove special characters)
#   safe_filename <- gsub("[^A-Za-z0-9_]", "_", region)
#   filename <- paste0("pathway_plots_by_region_fisher/", safe_filename, "_top_pathways.png")
#   
#   # Save plot
#   ggsave(filename, p, width = 12, height = 10, dpi = 300)
#   
#   cat("Saved:", filename, "\n")
# }
# 
# cat("\n=== SUMMARY ===\n")
# cat("Total regions plotted:", length(regions), "\n")
# cat("Files saved in: pathway_plots_by_region_fisher/\n")

# Method 2: Binomial Test 
# Tests if pathway frequency across regions exceeds random expectation.

total_regions <- n_distinct(df_clean$Matched_Brain_Region)

pathway_combined_BM <- df_clean %>%
  filter(padj < 0.05) %>%
  
  group_by(Matched_Brain_Region, pathway) %>%
  summarise(
    avg_NES = mean(NES, na.rm = TRUE),
    sd_NES = sd(NES, na.rm = TRUE),
    avg_padj = mean(padj, na.rm = TRUE),
    n_comparisons = n(),
    .groups = 'drop'
  ) %>%
  
  group_by(pathway) %>%
  mutate(
    n_regions = n_distinct(Matched_Brain_Region),
    n_occurrences = sum(n_comparisons),
    total_pathways = n_distinct(df_clean %>% filter(padj < 0.01) %>% pull(pathway))
  ) %>%
  ungroup() %>%
  
  # Binomial test: Is pathway more widespread than expected?
  rowwise() %>%
  mutate(
    expected_prob = 1 / total_pathways,  # null hypothesis probability
    binom_pval = binom.test(
      n_regions,           # observed successes
      total_regions,       # trials
      expected_prob,       # expected probability
      alternative = "greater"
    )$p.value
  ) %>%
  ungroup() %>%
  
  mutate(
    binom_padj = p.adjust(binom_pval, method = "BH")
  ) %>%
  # select(-total_pathways, -expected_prob) %>%
  
  arrange(Matched_Brain_Region, desc(abs(avg_NES)))
# Step 1: Identify top 20 pathways by cross-region prevalence
top_20_pathways <- pathway_combined_BM %>%
  group_by(pathway) %>%
  summarise(
    n_regions = dplyr::first(n_regions),
    binom_padj = dplyr::first(binom_padj),
    avg_abs_NES = mean(abs(avg_NES)),
    .groups = 'drop'
  ) %>%
  filter(binom_padj < 0.05) %>%
  arrange(desc(n_regions), binom_padj, desc(avg_abs_NES)) %>%
  slice_head(n = 20) %>%
  pull(pathway)

# Step 2: Prepare heatmap data with only top 20 pathways
heatmap_data <- pathway_combined_BM %>%
  filter(pathway %in% top_20_pathways) %>%
  filter(avg_padj < 0.05) %>%  # Only show region-specific significant results
  mutate(
    pathway_short = str_trunc(pathway, width = 100),
    NES_capped = pmax(pmin(avg_NES, 3), -3)  # Cap for visualization
  )

# Step 3: Create ordered factor for proper sorting
pathway_order <- pathway_combined_BM %>%
  filter(pathway %in% top_20_pathways) %>%
  group_by(pathway) %>%
  summarise(
    avg_NES_global = mean(avg_NES),
    n_regions = dplyr::first(n_regions),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_regions), desc(avg_NES_global)) %>%
  mutate(pathway_short = str_trunc(pathway, width = 100))

heatmap_data <- heatmap_data %>%
  mutate(
    pathway_short = factor(pathway_short, levels = pathway_order$pathway_short)
  )

# Step 4: Create heatmap
ggplot(heatmap_data, 
       aes(x = Matched_Brain_Region, y = pathway_short)) +
  geom_tile(aes(fill = NES_capped), color = "white", size = 0.5) +
  geom_text(aes(label = ifelse(binom_padj < 0.001, "***",
                               ifelse(binom_padj < 0.01, "**",
                                      ifelse(binom_padj < 0.05, "*", "")))),
            size = 3, fontface = "bold") +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "indianred",
    midpoint = 0,
    name = "Avg NES",
    limits = c(-3, 3)
  ) +
  labs(
    title = "Top 20 Cross-Region Pathways: Heatmap View",
    subtitle = "Pathways ranked by cross-region prevalence | *** Binomial padj < 0.001, ** < 0.01, * < 0.05",
    x = "Brain Region",
    y = "Pathway"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave("pathway_heatmap_top20_cross_region_BM_stats.png", width = 14, height = 11, dpi = 300, bg = "white")

# Print summary
cat("\n=== TOP 20 PATHWAYS SUMMARY ===\n")
print(pathway_order %>% select(pathway_short, n_regions, avg_NES_global))

# Get top 20 up and top 20 down per region
top_pathways_by_region <- pathway_combined_BM %>%
  filter(binom_padj < 0.05) %>%
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
      binom_padj < 0.001 ~ "***",
      binom_padj < 0.01 ~ "**",
      binom_padj < 0.05 ~ "*",
      TRUE ~ "ns"
    )
    
  )

# Create output directory for individual plots
dir.create("pathway_plots_by_region_BM", showWarnings = FALSE)

# Get list of unique brain regions
regions <- unique(top_pathways_by_region$Matched_Brain_Region)

# Loop through each region and create individual plots
for (region in regions) {
  
  # Filter data for this region
  region_data <- top_pathways_by_region %>%
    filter(Matched_Brain_Region == region)
  
  # Skip if no data
  if (nrow(region_data) == 0) next
  
  # Create plot for this region
  p <- ggplot(region_data, 
              aes(x = avg_NES, y = reorder(pathway_short, avg_NES))) +
    geom_segment(aes(x = 0, xend = avg_NES, 
                     y = pathway_short, yend = pathway_short,
                     color = direction),
                 linewidth = 1.5, alpha = 0.7) +
    # geom_point(aes(size = -log10(fisher_padj), color = direction),
    #            alpha = 0.8) +
    geom_text(aes(label = significance, 
                  x = avg_NES + ifelse(avg_NES > 0, 0.15, -0.15)),
              size = 4, fontface = "bold", 
              hjust = ifelse(region_data$avg_NES > 0, 0, 1)) +
    geom_vline(xintercept = 0, linetype = "solid", 
               color = "gray30", linewidth = 1) +
    scale_color_manual(values = c("Positive" = "steelblue", 
                                  "Negative" = "indianred")) +
    # scale_size_continuous(name = "-log10(Fisher padj)", range = c(3, 10)) +
    labs(
      title = paste("Top Pathways in", region),
      subtitle = "Binomial padj *** p < 0.001, ** p < 0.01, * p < 0.05",
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
  safe_filename <- gsub("[^A-Za-z0-9_]", "_", region)
  filename <- paste0("pathway_plots_by_region_BM/", safe_filename, "_top_pathways.png")
  
  # Save plot
  ggsave(filename, p, width = 12, height = 10, dpi = 300)
  
  cat("Saved:", filename, "\n")
}

cat("\n=== SUMMARY ===\n")
cat("Total regions plotted:", length(regions), "\n")
cat("Files saved in: pathway_plots_by_region_BM/\n")

