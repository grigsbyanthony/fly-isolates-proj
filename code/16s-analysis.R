# Blow fly isolate exploration and data visualization
# 2025-08-29

# ==================================================================================
# Load required libraries
# ==================================================================================
library(phyloseq)
library(qiime2R)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggthemes)
library(vegan)
library(RColorBrewer)
library(patchwork)
library(reshape2)
library(ggsignif)  # For adding significance bars
library(rstatix)   # For statistical testing

# ==================================================================================
# Load custom ggplot2 theme
# ==================================================================================
theme_pub <- function(base_size = 14, base_family = "CMU Sans Serif", legend_pos = "bottom") {
  ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      panel.background = element_blank(),
      plot.background = element_blank(),
      panel.border = element_blank(),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, margin = margin(r = 10)),
      axis.title.x = element_text(margin = margin(t = 5)),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_blank(),
      legend.position = legend_pos,
      legend.direction = "horizontal",
      legend.key.size = grid::unit(0.4, "cm"),
      legend.title = element_text(face = "bold"),
      plot.margin = grid::unit(c(10, 5, 5, 5), "mm"),
      strip.background = element_rect(fill = "#f0f0f0", colour = "#f0f0f0"),
      strip.text = element_text(face = "bold")
    )
} 

# ==================================================================================
# Load qiime2 artifacts as consolidated phyloseq object
# ==================================================================================
ps <- qiime2R::qza_to_phyloseq(
  features = "table-filtered-all.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

# Print phyloseq object summary
cat("Phyloseq object summary:\n")
print(ps)

# ==================================================================================
# Alpha Diversity Analysis
# ==================================================================================

# Calculate alpha diversity metrics
alpha_div <- estimate_richness(ps, measures = c("Observed", "Shannon"))

# Add sample metadata
sample_data_df <- data.frame(sample_data(ps))
alpha_div$SampleID <- rownames(alpha_div)

# Fix sample name formatting if needed
alpha_div$SampleID <- gsub("\\.", "-", alpha_div$SampleID)

# Merge with metadata
alpha_div <- merge(alpha_div, sample_data_df, by.x = "SampleID", by.y = "row.names")

# Reshape data to long format for faceting by alpha diversity metric
alpha_long <- alpha_div %>%
  select(SampleID, Sex, Observed, Shannon) %>%
  pivot_longer(
    cols = c(Observed, Shannon),
    names_to = "Alpha_Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Alpha_Metric = factor(Alpha_Metric, 
                         levels = c("Observed", "Shannon"),
                         labels = c("Observed ASVs", "Shannon Diversity"))
  )

# Create alpha diversity bar plot with points (lollipop style)
alpha_plot <- ggplot(alpha_long, aes(x = SampleID, y = Value, fill = SampleID, color = SampleID)) +
  geom_bar(stat = "identity", alpha = 0.7, width = 0.7) +
  geom_point(size = 5, position = position_dodge(width = 0.7)) +
  facet_wrap(~ Alpha_Metric, scales = "free_y", ncol = 1) +
  labs(
    x = "Sample ID",
    y = "Alpha Diversity Value"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

# Display the plot
print(alpha_plot)

# Save the plot
ggsave("alpha_diversity_barplot.png", 
       plot = alpha_plot,
       width = 10, height = 8, 
       dpi = 300, 
       bg = "white")

# Perform statistical tests to compare alpha diversity between sexes for each metric
stat_tests <- alpha_long %>%
  group_by(Alpha_Metric) %>%
  t_test(Value ~ Sex) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "Sex")

# Print the statistical test results
cat("\nStatistical test results for alpha diversity by Sex:\n")
print(stat_tests)

# Alternative visualization: Boxplot grouped by Sex with significance bars if applicable
alpha_grouped <- ggplot(alpha_long, aes(x = Sex, y = Value)) +  # Removed fill=Sex from here
  geom_boxplot(aes(fill = Sex), alpha = 0.7, outlier.shape = NA, show.legend = FALSE) +  # Keep fill but hide legend
  geom_jitter(width = 0.2, height = 0, alpha = 0.7, size = 4, aes(color = SampleID)) +
  facet_wrap(~ Alpha_Metric, scales = "free_y", ncol = 1) +
  labs(
    x = "Sex",
    y = "Alpha Diversity Value",
    color = "Sample ID"  # Rename the legend title
  ) +
  theme_pub(legend_pos = "bottom") +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  guides(fill = "none")  # Explicitly remove the fill guide

# Add significance bars if any tests are significant
if (any(stat_tests$p.adj < 0.05)) {
  alpha_grouped <- alpha_grouped +
    geom_signif(
      data = stat_tests %>% filter(p.adj < 0.05),
      aes(xmin = xmin, xmax = xmax, annotations = p.signif, y_position = y.position),
      manual = TRUE,
      tip_length = 0.01,
      vjust = 0.5
    )
}

# Display the grouped plot
print(alpha_grouped)

# Save the grouped plot
ggsave("alpha_diversity_by_sex.png", 
       plot = alpha_grouped,
       width = 10, height = 8, 
       dpi = 300, 
       bg = "white")

# ==================================================================================
# Beta Diversity Analysis - Bray-Curtis PCoA
# ==================================================================================

# Get OTU table for distance calculation
otu_table_data <- as.matrix(otu_table(ps))
if(taxa_are_rows(ps)) {
  otu_table_data <- t(otu_table_data)
}

# Calculate Bray-Curtis distance matrix
bray_dist <- vegdist(otu_table_data, method = "bray")

# Perform PCoA (Principal Coordinates Analysis)
pcoa_result <- cmdscale(bray_dist, eig = TRUE, k = nrow(otu_table_data) - 1)

# Extract coordinates and eigenvalues
pcoa_coords <- pcoa_result$points
eigenvalues <- pcoa_result$eig

# Calculate percentage of variance explained
percent_explained <- round(eigenvalues / sum(eigenvalues) * 100, 2)

# Create data frame for plotting
pcoa_df <- data.frame(
  Sample = rownames(pcoa_coords),
  PC1 = pcoa_coords[, 1],
  PC2 = pcoa_coords[, 2],
  stringsAsFactors = FALSE
)

# Add metadata
pcoa_df <- merge(pcoa_df, sample_data_df, by.x = "Sample", by.y = "row.names")

# Create PCoA plot colored by Sex
pcoa_plot <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Sex)) +
  geom_point(size = 5, alpha = 0.8) +
  stat_ellipse(aes(group = Sex), type = "norm", level = 0.95, linetype = 2) +
  scale_color_brewer(palette = "Set1") +
  labs(
    x = paste0("PC1 (", percent_explained[1], "% variance)"),
    y = paste0("PC2 (", percent_explained[2], "% variance)")
  ) +
  theme_pub(legend_pos = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 4)))

# Display the plot
print(pcoa_plot)

# Save the plot
ggsave("beta_diversity_bray_curtis_pcoa.png", 
       plot = pcoa_plot,
       width = 10, height = 8, 
       dpi = 300, 
       bg = "white")

# Perform PERMANOVA tests
# 1. Test for Sex effect
permanova_sex <- adonis2(bray_dist ~ Sex, 
                        data = sample_data_df, 
                        permutations = 999)

cat("\nPERMANOVA Results (testing for Sex effect):\n")
print(permanova_sex)

# 2. Test for Species effect
permanova_species <- adonis2(bray_dist ~ Species, 
                            data = sample_data_df, 
                            permutations = 999)

cat("\nPERMANOVA Results (testing for Species effect):\n")
print(permanova_species)

# 3. Test for interaction between Species and Sex
permanova_interaction <- adonis2(bray_dist ~ Species * Sex, 
                                data = sample_data_df, 
                                permutations = 999)

cat("\nPERMANOVA Results (testing for Species * Sex interaction):\n")
print(permanova_interaction)

# Export PERMANOVA results to a consolidated CSV file
# Convert each PERMANOVA result to a data frame and add a column for the test type
permanova_sex_df <- as.data.frame(permanova_sex)
permanova_sex_df$Test <- "Sex Effect"
permanova_sex_df$Factor <- rownames(permanova_sex_df)

permanova_species_df <- as.data.frame(permanova_species)
permanova_species_df$Test <- "Species Effect"
permanova_species_df$Factor <- rownames(permanova_species_df)

permanova_interaction_df <- as.data.frame(permanova_interaction)
permanova_interaction_df$Test <- "Species * Sex Interaction"
permanova_interaction_df$Factor <- rownames(permanova_interaction_df)

# Combine all results into a single data frame
all_permanova_results <- rbind(
  permanova_sex_df,
  permanova_species_df,
  permanova_interaction_df
)

# Reorder columns for better readability
all_permanova_results <- all_permanova_results[, c("Test", "Factor", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]

# Write to CSV
write.csv(all_permanova_results, "permanova_results.csv", row.names = FALSE)

cat("\nPERMANOVA results exported to 'permanova_results.csv'\n")

# ==================================================================================
# Relative Abundance Analysis
# ==================================================================================

# Transform to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# Function to prepare data for visualization
prepare_abundance_data <- function(ps_obj, tax_level, top_n = 10) {
  # Extract data from phyloseq object
  otu_table <- as.data.frame(otu_table(ps_obj))
  tax_table <- as.data.frame(tax_table(ps_obj))
  sample_data <- as.data.frame(sample_data(ps_obj))
  
  # Ensure OTUs are rows
  if (!taxa_are_rows(ps_obj)) {
    otu_table <- t(otu_table)
  }
  
  # Merge taxonomy with OTU table
  otu_tax <- merge(otu_table, tax_table, by = "row.names")
  rownames(otu_tax) <- otu_tax$Row.names
  otu_tax$Row.names <- NULL
  
  # Aggregate by taxonomic level
  tax_level_data <- otu_tax %>%
    group_by(!!sym(tax_level)) %>%
    summarise(across(where(is.numeric), sum)) %>%
    filter(!is.na(!!sym(tax_level)) & !!sym(tax_level) != "")
  
  # Get top taxa based on total abundance across all samples
  total_abundance <- rowSums(select(tax_level_data, where(is.numeric)))
  tax_level_data$TotalAbundance <- total_abundance
  
  top_taxa <- tax_level_data %>%
    arrange(desc(TotalAbundance)) %>%
    head(top_n) %>%
    pull(!!sym(tax_level))
  
  # Create a copy for reshaping
  tax_data_for_reshape <- tax_level_data %>%
    select(-TotalAbundance) %>%
    # Convert to long format
    pivot_longer(cols = -!!sym(tax_level), names_to = "SampleID", values_to = "Abundance")
  
  # Group low-abundance taxa into "Other" category PER SAMPLE
  plot_data <- tax_data_for_reshape %>%
    mutate(TaxaGroup = ifelse(!!sym(tax_level) %in% top_taxa, !!sym(tax_level), "Other")) %>%
    group_by(SampleID, TaxaGroup) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  # Add sample metadata
  plot_data <- merge(plot_data, sample_data, by.x = "SampleID", by.y = "row.names")
  
  return(plot_data)
}

# Create faceted barplot for Phylum and Genus
create_faceted_barplot <- function(ps_obj, tax_level = "Phylum", top_n = 10, facet_var = NULL) {
  # Prepare data
  plot_data <- prepare_abundance_data(ps_obj, tax_level, top_n)
  
  # Order taxa by overall abundance
  taxa_order <- plot_data %>%
    group_by(TaxaGroup) %>%
    summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
    arrange(desc(TotalAbundance)) %>%
    pull(TaxaGroup)
  
  # Replace underscores with spaces in genus names if needed
  if(tax_level == "Genus") {
    plot_data$TaxaGroup <- gsub("_", " ", plot_data$TaxaGroup)
    taxa_order <- gsub("_", " ", taxa_order)
  }
  
  # Set factor levels for ordering
  plot_data$TaxaGroup <- factor(plot_data$TaxaGroup, levels = rev(taxa_order))
  
  # Create color palette
  n_taxa <- length(unique(plot_data$TaxaGroup))
  if(n_taxa <= 10) {
    colors <- brewer.pal(max(3, n_taxa), "Set3")
  } else {
    colors <- colorRampPalette(brewer.pal(8, "Set3"))(n_taxa)
  }
  
  # Create barplot with faceting if facet_var is provided
  p <- ggplot(plot_data, aes(x = SampleID, y = Abundance, fill = TaxaGroup)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = colors, name = tax_level) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1), expand = c(0, 0)) +
    labs(
      x = "Sample ID",
      y = "Relative Abundance"
    ) +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right",
      legend.direction = "vertical"
    )
  
  # Add faceting if facet_var is provided
  if(!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = "free_x")
  }
  
  return(p)
}

# Create individual barplots for Phylum and Genus
phylum_barplot <- create_faceted_barplot(ps_rel, "Phylum", top_n = 10)
genus_barplot <- create_faceted_barplot(ps_rel, "Genus", top_n = 10)

# Create a combined plot with facets for Phylum and Genus
# First, prepare data for both taxonomic levels
phylum_data <- prepare_abundance_data(ps_rel, "Phylum", top_n = 3)  # Limit to top 3 phyla
genus_data <- prepare_abundance_data(ps_rel, "Genus", top_n = 10)

# Ensure we're keeping the top 10 genera and not grouping them into "Other"
cat("\nTop 10 genera by abundance:\n")
top_genera <- genus_data %>%
  group_by(TaxaGroup) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  filter(TaxaGroup != "Other") %>%
  head(10) %>%
  pull(TaxaGroup)
print(top_genera)

# Print the unique genera to debug
cat("\nUnique genera found:\n")
print(unique(genus_data$TaxaGroup))

# Add a column to identify the taxonomic level
phylum_data$TaxLevel <- "Phylum"
genus_data$TaxLevel <- "Genus"

# Combine the data
combined_data <- rbind(phylum_data, genus_data)

# Order taxa by overall abundance within each taxonomic level
phylum_order <- phylum_data %>%
  group_by(TaxaGroup) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  pull(TaxaGroup)

genus_order <- genus_data %>%
  group_by(TaxaGroup) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  pull(TaxaGroup)

# Replace underscores with spaces in genus names
genus_data$TaxaGroup <- gsub("_", " ", genus_data$TaxaGroup)
combined_data$TaxaGroup <- ifelse(combined_data$TaxLevel == "Genus", 
                                 gsub("_", " ", combined_data$TaxaGroup),
                                 combined_data$TaxaGroup)

# Make sure we're not losing genera in the normalization process
cat("\nUnique genera after underscore replacement:\n")
print(unique(genus_data$TaxaGroup))

# Create color palettes for each taxonomic level
phylum_colors <- brewer.pal(min(length(unique(phylum_data$TaxaGroup)), 9), "Set1")
genus_colors <- brewer.pal(min(length(unique(genus_data$TaxaGroup)), 9), "Set3")

# Ensure all samples sum to 100% for each taxonomic level
combined_data_normalized <- combined_data %>%
  group_by(TaxLevel, SampleID) %>%
  mutate(Total = sum(Abundance)) %>%
  mutate(Abundance = Abundance / Total) %>%
  select(-Total)

# Create separate color palettes for each taxonomic level
phylum_data_subset <- combined_data_normalized %>% filter(TaxLevel == "Phylum")
genus_data_subset <- combined_data_normalized %>% filter(TaxLevel == "Genus")

phylum_unique_taxa <- unique(phylum_data_subset$TaxaGroup)
genus_unique_taxa <- unique(genus_data_subset$TaxaGroup)

n_phylum <- length(phylum_unique_taxa)
n_genus <- length(genus_unique_taxa)

# Custom color palette for phylum
# Reorder phylum taxa to put "Other" at the bottom of the legend
if("Other" %in% phylum_unique_taxa) {
  phylum_unique_taxa_ordered <- c(phylum_unique_taxa[phylum_unique_taxa != "Other"], "Other")
} else {
  phylum_unique_taxa_ordered <- phylum_unique_taxa
}

phylum_colors <- colorRampPalette(brewer.pal(min(9, max(3, n_phylum)), "Set1"))(n_phylum)
names(phylum_colors) <- phylum_unique_taxa_ordered

# Custom color palette for genus with "Other" as grey
genus_custom_palette <- c("#ea5545", "#f46a9b", "#ef9b20", "#edbf33", "#ede15b", "#bdcf32", "#87bc45", "#27aeef", "#b33dc6")

# Make sure we have non-Other taxa
if(length(genus_unique_taxa) > 1) {
  # Assign colors to taxa (excluding "Other")
  non_other_taxa <- genus_unique_taxa[genus_unique_taxa != "Other"]
  
  # Make sure we don't exceed the palette length
  genus_colors <- genus_custom_palette[1:min(length(genus_custom_palette), length(non_other_taxa))]
  names(genus_colors) <- non_other_taxa
  
  # Add grey for "Other" category
  if("Other" %in% genus_unique_taxa) {
    genus_colors <- c(genus_colors, "Other" = "#808080")  # Grey for "Other"
  }
} else {
  # If we only have "Other", just use grey
  genus_colors <- c("Other" = "#808080")
}

# Print the color assignments to debug
cat("\nGenus color assignments:\n")
print(genus_colors)

# Combine the color palettes
all_colors <- c(phylum_colors, genus_colors)

# Create the faceted barplot with separate legends
taxa_barplot <- ggplot(combined_data_normalized, aes(x = SampleID, y = Abundance, fill = TaxaGroup)) +
  geom_col(position = "stack") +
  facet_grid(TaxLevel ~ ., scales = "free_y", switch = "y") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = all_colors, breaks = c(phylum_unique_taxa, genus_unique_taxa)) +
  labs(
    x = "Sample ID",
    y = "Relative Abundance"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.direction = "vertical",
    strip.text = element_text(size = 12, face = "bold")
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE))

# Reorder genus data to put "Other" at the bottom of the stack
if(length(genus_unique_taxa) > 1) {
  genus_data_subset <- genus_data_subset %>%
    mutate(TaxaGroup = factor(TaxaGroup, 
                             levels = c(non_other_taxa[order(non_other_taxa)], "Other")))
}

# Reorder phylum data to put "Other" at the bottom of the stack
if("Other" %in% phylum_unique_taxa) {
  phylum_data_subset <- phylum_data_subset %>%
    mutate(TaxaGroup = factor(TaxaGroup, 
                             levels = c(phylum_unique_taxa[phylum_unique_taxa != "Other" & phylum_unique_taxa != ""], "Other")))
}

# Create a direct approach for the genus plot using the original genus_data
# This bypasses any potential issues with the combined_data_normalized
genus_data_direct <- genus_data %>%
  # Replace underscores with spaces
  mutate(TaxaGroup = gsub("_", " ", TaxaGroup)) %>%
  # Normalize to ensure each sample sums to 100%
  group_by(SampleID) %>%
  mutate(Total = sum(Abundance)) %>%
  mutate(Abundance = Abundance / Total) %>%
  select(-Total) %>%
  ungroup()

# Order taxa by overall abundance
genus_order_direct <- genus_data_direct %>%
  group_by(TaxaGroup) %>%
  summarise(TotalAbundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(TotalAbundance)) %>%
  pull(TaxaGroup)

# Reorder with "Other" at the bottom
if("Other" %in% genus_order_direct) {
  genus_order_direct <- c(genus_order_direct[genus_order_direct != "Other"], "Other")
}

# Set factor levels for ordering
genus_data_direct$TaxaGroup <- factor(genus_data_direct$TaxaGroup, levels = genus_order_direct)

# Create custom color palette for genus
genus_custom_palette <- c("#ea5545", "#f46a9b", "#ef9b20", "#edbf33", "#ede15b", "#bdcf32", "#87bc45", "#27aeef", "#b33dc6")
genus_colors_direct <- c()

# Assign colors to non-Other taxa
non_other_taxa_direct <- genus_order_direct[genus_order_direct != "Other"]
for(i in 1:length(non_other_taxa_direct)) {
  if(i <= length(genus_custom_palette)) {
    genus_colors_direct[non_other_taxa_direct[i]] <- genus_custom_palette[i]
  } else {
    # If we run out of colors in the palette, use a default color
    genus_colors_direct[non_other_taxa_direct[i]] <- "#000000"
  }
}

# Add grey for "Other" category
if("Other" %in% genus_order_direct) {
  genus_colors_direct["Other"] <- "#808080"  # Grey for "Other"
}

# Create separate plots for Phylum and Genus with their own legends
phylum_plot <- ggplot(phylum_data_subset, aes(x = SampleID, y = Abundance, fill = TaxaGroup)) +
  geom_col(position = "stack") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = phylum_colors, name = "Phylum") +
  labs(
    x = "Sample ID",
    y = "Relative Abundance"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.direction = "vertical",
    strip.text = element_text(size = 12, face = "bold")
  )

# Use the direct approach for the genus plot
genus_plot <- ggplot(genus_data_direct, aes(x = SampleID, y = Abundance, fill = TaxaGroup)) +
  geom_col(position = "stack") +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = genus_colors_direct, name = "Genus") +
  labs(
    x = "Sample ID",
    y = "Relative Abundance"
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.direction = "vertical",
    strip.text = element_text(size = 12, face = "bold")
  )

# Combine the plots using patchwork
combined_taxa_plot <- phylum_plot / genus_plot +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

# Display the plots
print(taxa_barplot)
print(combined_taxa_plot)

# Save the plots
ggsave("taxonomic_relative_abundance_faceted.png", 
       plot = taxa_barplot,
       width = 14, height = 10, 
       dpi = 300, 
       bg = "white")

ggsave("taxonomic_relative_abundance_separate_legends.png", 
       plot = combined_taxa_plot,
       width = 14, height = 12, 
       dpi = 300, 
       bg = "white")

# ==================================================================================
# Create a table of top 10 most abundant features by taxonomic level
# ==================================================================================

# Extract ASV table and taxonomy table from phyloseq object
asv_table_df <- as.data.frame(otu_table(ps_rel))
tax_table_df <- as.data.frame(tax_table(ps_rel))

# Ensure ASVs are rows
if (!taxa_are_rows(ps_rel)) {
  asv_table_df <- t(asv_table_df)
}

# Calculate total abundance for each ASV across all samples
asv_table_df$TotalAbundance <- rowSums(asv_table_df)

# Add row names as a column
asv_table_df$ASV <- rownames(asv_table_df)

# Merge with taxonomy table
merged_data <- merge(asv_table_df, tax_table_df, by.x = "ASV", by.y = "row.names")

# Select relevant columns and sort by total abundance
top_features <- merged_data %>%
  select(ASV, TotalAbundance, Kingdom, Phylum, Class, Order, Family, Genus) %>%
  arrange(desc(TotalAbundance)) %>%
  head(10)  # Get top 10 most abundant features

# Clean up taxonomy names (remove prefixes like "k__", "p__", etc.)
top_features <- top_features %>%
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus), 
                ~ gsub("^[kpcofgs]__", "", .))) %>%
  # Replace empty strings with NA
  mutate(across(c(Kingdom, Phylum, Class, Order, Family, Genus), 
                ~ ifelse(. == "", NA, .)))

# Create a more readable table without relative abundance
top_features_table <- top_features %>%
  select(ASV, Phylum, Order, Family, Genus)

# Print the table
cat("\nTop 10 most abundant features across all samples:\n")
print(top_features_table)

# Save the table to CSV
write.csv(top_features_table, "top_10_abundant_features.csv", row.names = FALSE)
cat("\nTable of top 10 most abundant features saved to 'top_10_abundant_features.csv'\n")

# Create a more compact table for visualization
compact_table <- top_features %>%
  select(Phylum, Order, Family, Genus)

# Save the compact table to CSV
write.csv(compact_table, "top_10_abundant_features_compact.csv", row.names = FALSE)
cat("\nCompact table of top 10 most abundant features saved to 'top_10_abundant_features_compact.csv'\n")

# Create a formatted table for publication
library(kableExtra)
if(requireNamespace("kableExtra", quietly = TRUE)) {
  formatted_table <- top_features %>%
    select(Phylum, Order, Family, Genus)
  
  # Create a nicely formatted HTML table
  html_table <- kable(formatted_table, format = "html", caption = "Top 10 Most Abundant Features") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)
  
  # Save the HTML table
  writeLines(as.character(html_table), "top_10_abundant_features_table.html")
  cat("\nFormatted HTML table saved to 'top_10_abundant_features_table.html'\n")
} else {
  cat("\nkableExtra package not available. Skipping formatted table creation.\n")
}

# ==================================================================================
# Create an alluvial plot for taxonomic hierarchy of top 20 features
# ==================================================================================

# Check if ggalluvial package is installed, if not, install it
if (!requireNamespace("ggalluvial", quietly = TRUE)) {
  install.packages("ggalluvial")
}
library(ggalluvial)

# Prepare data for alluvial plot
# First, extract the taxonomic data for the top 20 features
alluvial_data <- top_features %>%
  select(ASV, TotalAbundance, Phylum, Order, Family, Genus) %>%
  # Replace NA values with "Unknown" for visualization
  mutate(across(c(Phylum, Order, Family, Genus), ~ifelse(is.na(.), "Unknown", .))) %>%
  # Filter out features where Genus is "Unknown" or "Streptococcus"
  filter(Genus != "Unknown" & Genus != "Streptococcus") %>%
  # Consolidate Enterobacteriaceae families and simplify order names
  mutate(
    Family = case_when(
      Family %in% c("Enterobacteriaceae_A_732334", "Enterobacteriaceae_A_725029") ~ "Enterobacteriaceae A",
      TRUE ~ Family
    ),
    Order = case_when(
      Order == "Enterobacterales_737866" ~ "Enterobacterales",
      TRUE ~ Order
    ),
    Genus = case_when(
      Genus == "Vagococcus_B" ~ "Vagococcus B",
      TRUE ~ Genus
    ),
    Phylum = case_when(
      Phylum == "Bacillota_I" ~ "Bacillota I",
      TRUE ~ Phylum
    )
  )

# Print the number of features after filtering
cat("\nNumber of features after removing Unknown genera and Streptococcus:", nrow(alluvial_data), "\n")

# Create a long format dataset for the alluvial plot
alluvial_long <- alluvial_data %>%
  pivot_longer(
    cols = c(Phylum, Order, Family, Genus),
    names_to = "TaxonomicLevel",
    values_to = "Taxon"
  ) %>%
  # Set the taxonomic level as a factor with the correct order
  mutate(TaxonomicLevel = factor(TaxonomicLevel, levels = c("Phylum", "Order", "Family", "Genus")))

# Create a frequency table for the alluvial plot
alluvial_freq <- alluvial_long %>%
  group_by(TaxonomicLevel, Taxon) %>%
  summarise(Frequency = n(), .groups = "drop") %>%
  # Add a column for the total abundance
  left_join(
    alluvial_long %>%
      group_by(TaxonomicLevel, Taxon) %>%
      summarise(TotalAbundance = sum(TotalAbundance), .groups = "drop"),
    by = c("TaxonomicLevel", "Taxon")
  )

# Create a color palette for the taxa
# Get unique taxa across all levels
unique_taxa <- unique(alluvial_freq$Taxon)
n_taxa <- length(unique_taxa)

# Generate a color palette
if (n_taxa <= 20) {
  # Use a combination of Set1, Set2, and Set3 for up to 20 colors
  colors <- c(
    brewer.pal(min(9, n_taxa), "Set1"),
    brewer.pal(min(8, max(0, n_taxa - 9)), "Set2"),
    brewer.pal(min(12, max(0, n_taxa - 17)), "Set3")
  )[1:n_taxa]
} else {
  # For more than 20 taxa, use a colorRampPalette
  colors <- colorRampPalette(
    c(
      brewer.pal(9, "Set1"),
      brewer.pal(8, "Set2"),
      brewer.pal(12, "Set3")
    )
  )(n_taxa)
}

# Assign colors to taxa
names(colors) <- unique_taxa

# Create the alluvial plot
# First, prepare the data for the alluvial plot
alluvial_plot_data <- alluvial_data %>%
  # Add a unique ID for each ASV
  mutate(id = row_number())

# Calculate abundance by taxonomic level for ordering
phylum_abundance <- alluvial_plot_data %>%
  group_by(Phylum) %>%
  summarise(TotalAbundance = sum(TotalAbundance)) %>%
  arrange(desc(TotalAbundance))

order_abundance <- alluvial_plot_data %>%
  group_by(Order) %>%
  summarise(TotalAbundance = sum(TotalAbundance)) %>%
  arrange(desc(TotalAbundance))

family_abundance <- alluvial_plot_data %>%
  group_by(Family) %>%
  summarise(TotalAbundance = sum(TotalAbundance)) %>%
  arrange(desc(TotalAbundance))

genus_abundance <- alluvial_plot_data %>%
  group_by(Genus) %>%
  summarise(TotalAbundance = sum(TotalAbundance)) %>%
  arrange(desc(TotalAbundance))

# Set factor levels for ordering (most abundant at the top)
alluvial_plot_data$Phylum <- factor(alluvial_plot_data$Phylum, levels = phylum_abundance$Phylum)
alluvial_plot_data$Order <- factor(alluvial_plot_data$Order, levels = order_abundance$Order)
alluvial_plot_data$Family <- factor(alluvial_plot_data$Family, levels = family_abundance$Family)
alluvial_plot_data$Genus <- factor(alluvial_plot_data$Genus, levels = genus_abundance$Genus)

# Create the alluvial plot with ordered strata
alluvial_plot <- ggplot(
  alluvial_plot_data,
  aes(
    y = TotalAbundance,
    axis1 = Phylum,
    axis2 = Order,
    axis3 = Family,
    axis4 = Genus
  )
) +
  geom_alluvium(aes(fill = Phylum), width = 0.4, alpha = 0.8) +  # Adjusted width to 0.4
  geom_stratum(width = 0.4, fill = "white", color = "grey") +    # Adjusted width to 0.4
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, family = "CMU Sans Serif") +
  scale_x_discrete(
    limits = c("Phylum", "Order", "Family", "Genus"),
    expand = c(0.05, 0.05)
  ) +
  scale_y_continuous(expand = c(0, 0)) +  # Remove padding around y-axis
  scale_fill_manual(values = colors) +
  labs(
    # Title removed
    y = "Relative Abundance",
    fill = "Phylum"
  ) +
  theme_pub(legend_pos = "bottom") +  # Changed legend position to bottom
  theme(
    legend.position = "bottom",       # Ensure legend is at the bottom
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.key.width = unit(1, "cm"),  # Wider legend keys
    legend.key.height = unit(0.5, "cm")
  ) +
  guides(fill = guide_legend(ncol = 3))  # Display legend in 3 columns for better layout

# Save the alluvial plot
ggsave("taxonomic_hierarchy_alluvial_top10.png", 
       plot = alluvial_plot,
       width = 14, height = 10, 
       dpi = 300, 
       bg = "white")

cat("\nAlluvial plot of taxonomic hierarchy saved to 'taxonomic_hierarchy_alluvial_top10.png'\n")

# Create a more detailed alluvial plot with lodes
# This version shows each individual ASV as a separate flow
alluvial_lodes <- to_lodes_form(
  alluvial_plot_data,
  axes = c("Phylum", "Order", "Family", "Genus"),
  id = "id"
)

# Add total abundance to the lodes data
alluvial_lodes <- alluvial_lodes %>%
  left_join(
    alluvial_plot_data %>% select(id, TotalAbundance),
    by = "id"
  )

# Order the strata in the lodes data
# First, create a mapping of taxonomic level to ordered taxa
taxa_order_map <- list(
  Phylum = phylum_abundance$Phylum,
  Order = order_abundance$Order,
  Family = family_abundance$Family,
  Genus = genus_abundance$Genus
)

# Apply the ordering to the lodes data
alluvial_lodes <- alluvial_lodes %>%
  mutate(
    stratum = factor(
      stratum,
      levels = taxa_order_map[[as.character(x)]]
    )
  )

# Create the detailed alluvial plot with ordered strata
detailed_alluvial_plot <- ggplot(
  alluvial_lodes,
  aes(x = x, y = TotalAbundance, stratum = stratum, alluvium = id, fill = stratum, label = stratum)
) +
  geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "darkgray", alpha = 0.7) +
  geom_stratum(width = 0.13, alpha = 0.8) +  # Adjusted width to 0.13 (smaller value = wider bars)
  geom_text(stat = "stratum", size = 3, family = "CMU Sans Serif") +
  scale_x_discrete(
    limits = c("Phylum", "Order", "Family", "Genus"),
    expand = c(0.05, 0.05)
  ) +
  scale_y_continuous(expand = c(0, 0)) +  # Remove padding around y-axis
  labs(
    # Title removed
    y = "Relative Abundance",
    fill = "Taxon"
  ) +
  theme_pub(legend_pos = "bottom") +  # Changed legend position to bottom
  theme(
    legend.position = "bottom",       # Ensure legend is at the bottom
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.key.width = unit(1, "cm"),  # Wider legend keys
    legend.key.height = unit(0.5, "cm")
  ) +
  guides(fill = guide_legend(ncol = 4))  # Display legend in 4 columns for better layout

# Save the detailed alluvial plot
ggsave("taxonomic_hierarchy_detailed_alluvial_top10.png", 
       plot = detailed_alluvial_plot,
       width = 16, height = 12, 
       dpi = 300, 
       bg = "white")

cat("\nDetailed alluvial plot of taxonomic hierarchy saved to 'taxonomic_hierarchy_detailed_alluvial_top10.png'\n")
