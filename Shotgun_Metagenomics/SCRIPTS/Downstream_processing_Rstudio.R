############ Shotgun metagenomicss
# Load required packages
library(qiime2R)
library(mia)
library(miaViz)
library(TreeSummarizedExperiment)
library(dplyr)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(ecodist)
library(data.table)
library(vegan)
library(ggforce)
library(ggplot2)
library(readr)
library(tidyr)
library(scater)
library(devtools)


##TSE from Metaphlan
tse_Metaphlan <- importMetaPhlAn("path/metaphlan_db_meta4_combined_reports.txt", "Z:/PhD CGTG/Experiments/Metagenomics/human_urine/Batch 4/Puhti/metadata.tsv", package = "mia")

# change MI_ID rownames of the tse, to our sample_ID
col_data_Metaphlan <- colData(tse_Metaphlan)
rownames(col_data_Metaphlan) <- col_data_Metaphlan$Patient_day
colData(tse_Metaphlan) <- col_data_Metaphlan
col_data_Metaphlan <- as.data.frame(colData(tse_Metaphlan))


# Use species as rownames instead of strains
# Check for NAs in species column
sum(is.na(rowData(tse_Metaphlan)$species))

# Check for duplicate species names
sum(duplicated(rowData(tse_Metaphlan)$species))

# If duplicates exist, make them unique
rownames(tse_Metaphlan) <- make.unique(rowData(tse_Metaphlan)$species)


# Set rownames to species
rownames(tse_Metaphlan) <- rowData(tse_Metaphlan)$species

# Verify
rownames(tse_Metaphlan) |> head()


# Check relative abundance, is in %
colSums(assay(tse_Metaphlan, 1))

assay(tse_Metaphlan, "relabundance") <- assay(tse_Metaphlan, "metaphlan") / 100

# Verify
colSums(assay(tse_Metaphlan, "relabundance")) |> head()  # should be ~1





################################################
## Alpha diversity - Shannon - Genus/species

########## Response #############
tse <- tse_Metaphlan

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Trial == "T563")]
tse <- tse[, which(colData(tse)$Response %in% c("No_DC", "DC"))]

tse <- mergeFeaturesByRank(tse, rank ="species", onRankOnly=TRUE)

# Estimate (observed) richness
tse <- addAlpha(
  tse, assay.type = "relabundance", index = "shannon", name = "observed",
  detection = 10)

# Check some of the first values in colData
head(tse$observed)

# Set the specific order for the Patient factor levels
colData(tse)$Patient <- factor(colData(tse)$Patient, 
                               levels = c("30104", "30105", "30207", "30209", 
                                          "30110", "30111", "30112", "30210", 
                                          "30211", "30212", "20216", "20217", 
                                          "20219", "20111", "10114", "10115", 
                                          "10116"))

# Plot with patients in the specified order
plotColData(
  tse,
  y = "observed",      
  x = "Patient",       
  colour_by = "Response",   
  point_size = 4 +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
        axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size = 18),  
        axis.title.y = element_text(size = 18),  
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  
        legend.text = element_text(size = 14),  
        legend.title = element_text(size = 16),  
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
  ) +
  labs(
    title = "Alpha diversity (Shannon)",
    x = "
    Patient", 
    y = expression(Richness[Observed]), 
    colour = "Response"  
  ) +
  scale_colour_manual(values = c("DC" = "blue", "No_DC" = "orange")) +
  guides(colour = guide_legend(override.aes = list(size = 6)))  


library(dplyr)
library(ggplot2)
library(ggsignif)

# Build plotting data from colData
df_alpha <- as.data.frame(colData(tse)) %>%
  dplyr::select(Response, observed) %>%
  dplyr::filter(!is.na(Response), !is.na(observed)) %>%
  dplyr::mutate(
    Response = factor(Response, levels = c("DC", "No_DC"))
  )

#  Colours 
fill_cols   <- c("DC" = "blue", "No_DC" = "orange")
border_cols <- c("DC" = "#6E6E6E", "No_DC" = "#C07D10")

#  Plot 
violin_alpha <- ggplot(df_alpha, aes(x = Response, y = observed)) +
  
  geom_violin(
    aes(fill = Response),
    color = "black",
    width = 0.6,
    trim  = FALSE,
    alpha = 0.7
  ) +
  
  geom_boxplot(
    aes(fill = Response),
    color         = "black",
    width         = 0.12,        
    outlier.shape = NA
  ) +
  
  geom_point(
    aes(color = Response),
    position = position_jitter(width = 0.08, height = 0, seed = 42),
    alpha    = 0.8,
    size     = 2.5
  ) +
  
  geom_signif(
    comparisons      = list(c("DC", "No_DC")),
    test             = "wilcox.test",
    map_signif_level = function(p) paste0("p = ", signif(p, 3)),
    step_increase    = 0.12,
    color            = "black",
    size             = 0.6,
    textsize         = 4
  ) +
  
  labs(
    x     = NULL,
    y     ="Richness",
    title = "Alpha Diversity (Shannon)"
  ) +
  
  scale_fill_manual(values  = fill_cols) +
  scale_color_manual(values = border_cols) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  
  theme_minimal() +
  theme(
    plot.title       = element_text(face = "bold", size = 14,
                                    hjust = 0.5, margin = margin(b = 10)),
    axis.title.y     = element_text(size = 15, margin = margin(r = 10)),
    axis.text        = element_text(size = 13),
    axis.text.x      = element_text(color = "black", size = 14),
    legend.position  = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black", linewidth = 0.6)
  )

print(violin_alpha)




##############################################
## Beta diversity - genus/species

  tse <- tse_Metaphlan
tse <- tse[, which(colData(tse)$Trial == "T563")]

# species level
tse <- mergeFeaturesByRank(tse, rank ="species", onRankOnly=TRUE)

tse <- tse[, which(colData(tse)$Response %in% c("DC", "No_DC"))]

# Extract day 1 samples
tse_rel <- tse
tse_day1 <- tse_rel[, rownames(subset(colData(tse_rel), Day == "1"))]

# Or Include only species present in at least 2 or more patients
species_to_keep <- rowSums(assay(tse_day1, "relabundance") > 0) >= 2
tse_filtered <- tse_day1[species_to_keep, ]

tse_day1_rel_abund_assay <- assays(tse_filtered)$relabundance
bray_curtis_dist <- vegan::vegdist(t(tse_day1_rel_abund_assay), method = "bray")
#install.packages("ecodist")
#library(ecodist)
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
#PCoA1 and PCoA2
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])
# plot
bray_curtis_pcoa_df <- cbind(bray_curtis_pcoa_df,
                             patient_status = colData(tse_filtered)$Response,
                             label = colData(tse_filtered)$Patient)


# For DC vs No_DC
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, 
                           aes(x = pcoa1, y = pcoa2, color = patient_status)) +
  geom_point() + 
  geom_text(aes(label = label), size = 3, vjust = 1.5, show.legend = FALSE) +  
  
  # Add circles around the groups with colored outlines 
  geom_mark_ellipse(aes(group = patient_status), 
                    fill = NA,  
                    linetype = "solid", 
                    show.legend = FALSE) +  
  
  labs(x = "PC1",
       y = "PC2", 
       colour = "Response",
       title = "Beta diversity between groups") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  
        axis.text.y = element_text(size = 14),  
        axis.title.x = element_text(size = 14),  
        axis.title.y = element_text(size = 14),  
        plot.title = element_text(size = 12, face = "bold"),  
        legend.text = element_text(size = 14),  
        legend.title = element_text(size = 16),  
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),  
        panel.background = element_rect(fill = "white", color = NA), 
        axis.line = element_line(color = "black", linewidth = 0.4)
        ) +
  scale_x_continuous(expand = expansion(mult = 0.3)) + 
  scale_y_continuous(expand = expansion(mult = 0.3)) +
  scale_color_manual(values = c("DC" = "blue", "No_DC" = "orange")) + 
  guides(colour = guide_legend(override.aes = list(size = 6)))  

# Display the plot
bray_curtis_plot



### Differential abundance analysis - LinDA 
#install.packages("MicrobiomeStat")
library(MicrobiomeStat)

tse <- tse_Metaphlan

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Trial == "T563")]
tse <- tse[, which(colData(tse)$Response %in% c("No_DC", "DC"))]

# Extract species level
tse_spec <- mergeFeaturesByRank(tse, rank ="species", onRankOnly=TRUE)

# Replace all spaces and hyphens with underscores in row names of tse_species
rownames(tse_spec) <- gsub("[ -]", "_", rownames(tse_spec))

# Run LinDA
linda_out <- linda(feature.dat = as.data.frame(assay(tse_spec)),
                   meta.dat = as.data.frame(colData(tse_spec)),
                   formula = "~ Response",
                   alpha = 0.01,
                   prev.filter = 0,
                   mean.abund.filter = 0,
                   n.cores = 1)

head(linda_out)

# Extract the results data frame from linda_out
results <- linda_out$output$ResponseNo_DC

# Convert row names to a column named 'taxon'
results$taxon <- rownames(results)

head(results)
# Filter for significance
significant_results <- results[results$pvalue < 0.3, ]
# Flip the direction if needed
significant_results$log2FoldChange <- -significant_results$log2FoldChange
# Create a new column to represent the Response
significant_results$Response <- ifelse(significant_results$log2FoldChange > 0, "DC", "No_DC")
print(significant_results)

# Remove "s__" prefix
significant_results$taxon <- gsub("^s__", "", significant_results$taxon)
# Remove rows where "taxon" is empty or NA
significant_results <- significant_results[significant_results$taxon != "" & !is.na(significant_results$taxon), ]

significant_results$highlight <- ifelse(significant_results$Response == "DC", "lightblue", "#FFC300")

# Border color column
significant_results$border <- ifelse(significant_results$taxon == "Prevotella_timonensis", "darkred", NA)

# Text color column
significant_results$text_color <- ifelse(significant_results$taxon == "Prevotella_timonensis", "darkred", "black")

# Create the bar plot
bar_plot <- ggplot(significant_results, aes(x = reorder(taxon, log2FoldChange), 
                                            y = log2FoldChange, 
                                            fill = highlight,
                                            color = border)) +
  geom_bar(stat = "identity", show.legend = TRUE, linewidth = 0.8) +
  scale_fill_manual(name = NULL,
                    values = c("lightblue" = "lightblue", "#FFC300" = "#FFC300"),
                    breaks = c("lightblue", "#FFC300"),
                    labels = c("lightblue" = "DC", 
                               "#FFC300"   = "No_DC")) +
  scale_color_identity() +
  coord_flip() +
  labs(title = "Differential abundant species 
      (shotgun metagenomics)", x = "Taxon
       ", y = "Log Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14, face = "italic",
                                   color = significant_results$text_color[order(significant_results$log2FoldChange)]),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),  
        plot.title = element_text(size = 20, face = "bold", hjust = 0),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )

print(bar_plot)



####### genus

tse <- tse_Metaphlan

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Trial == "T563")]
tse <- tse[, which(colData(tse)$Response %in% c("No_DC", "DC"))]

# Extract genus level
tse_spec <- mergeFeaturesByRank(tse, rank ="genus", onRankOnly=TRUE)

# Replace all spaces and hyphens with underscores in row names of tse_species
rownames(tse_spec) <- gsub("[ -]", "_", rownames(tse_spec))

# Run LinDA
linda_out <- linda(feature.dat = as.data.frame(assay(tse_spec)),
                   meta.dat = as.data.frame(colData(tse_spec)),
                   formula = "~ Response",
                   alpha = 0.05,
                   prev.filter = 0,
                   mean.abund.filter = 0)

head(linda_out)

# Extract the results data frame from linda_out
results <- linda_out$output$ResponseNo_DC

# Convert row names to a column named 'taxon'
results$taxon <- rownames(results)

head(results)

# Filter for significant results based on adjusted p-value threshold (adjust as needed)
significant_results <- results[results$pvalue < 0.3, ]
print(significant_results)

# Filter for significant results based on log2FoldChange thresholds
significant_results <- significant_results[significant_results$log2FoldChange < 0 | significant_results$log2FoldChange > 0,]
print(significant_results)

# Create a new column to represent the log2FoldChange with a negative sign for DC group
significant_results$log2FoldChange <- -significant_results$log2FoldChange

# Create a new column for color mapping based on log2FoldChange (blue for DC, orange for No_DC)
significant_results$highlight <- ifelse(significant_results$log2FoldChange > 0, "blue", "orange")
print(significant_results)

# Create a new column to represent the Response (DC = blue, No_DC = orange)
significant_results$Response <- ifelse(significant_results$log2FoldChange > 0, "DC", "No_DC")
print(significant_results)
# Remove "g__" prefix
significant_results$taxon <- gsub("^g__", "", significant_results$taxon)
# Remove rows where "taxon" is empty or NA
significant_results <- significant_results[significant_results$taxon != "" & !is.na(significant_results$taxon), ]

# Create the bar plot with the 'Response' column for color mapping
bar_plot <- ggplot(significant_results, aes(x = reorder(taxon, log2FoldChange), y = log2FoldChange, fill = Response)) +
  geom_bar(stat = "identity", show.legend = TRUE) +  
  scale_fill_manual(values = c("DC" = "blue", "No_DC" = "orange")) + 
  coord_flip() +  
  labs(title = "Differential abundant genus 
  (shotgun metagenomics)", x = "Taxon", y = "Log Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 14),  
        axis.text.y = element_text(size = 14, face = "italic"),  
        axis.title.x = element_text(size = 18),  
        axis.title.y = element_text(size = 18),  
        plot.title = element_text(size = 20, face = "bold", hjust = 0),  
        legend.text = element_text(size = 14), 
        legend.title = element_text(size = 16),  
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
  ) +
  guides(fill = guide_legend(title = "", 
                             labels = c("DC", "No_DC")))  

# Print the plot
print(bar_plot)





### Wilcoxon barplot
# Prevotella timonensis 
library(writexl)
library(ggsignif)
library(ggrepel)

tse <- tse_Metaphlan

tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Trial == "T563")]
tse <- tse[, which(colData(tse)$Response %in% c("No_DC", "DC"))]

grep("Prevotella", rowData(tse)$species, value = TRUE)

# Extract species level
tse_species <- mergeFeaturesByRank(tse, rank = "species", onRankOnly = TRUE)

# Replace all spaces and hyphens with underscores in row names of tse_species
rownames(tse_species) <- gsub("[ -]", "_", rownames(tse_species))

metadata_df <- as.data.frame(colData(tse))

# Convert colData to a data frame
metadata_df <- metadata_df[metadata_df$Day == "1" & metadata_df$Response %in% c("DC", "No_DC"), ]

# Transpose assay data
tse_assay <- t(assay(tse, "relabundance"))  

# Merge metadata with transposed assay data
merged_data <- merge(metadata_df, as.data.frame(tse_assay), by.x = "row.names", by.y = "row.names", all.x = TRUE)

# DC first in boxplot
merged_data$Response <- factor(merged_data$Response, levels = c("DC", "No_DC"))

# Initialize a list to store all plots
plot_list <- list()

# Use all columns from the merged_data (or merged_data)
colnames(merged_data) <- gsub("\\[|\\]", "", colnames(merged_data))
all_columns <- colnames(merged_data)[11:ncol(merged_data)]
print(all_columns)
## Change the 11:ncol from where the first species starts!!

#check relative abundance
row_sums <- rowSums(merged_data[, 16:ncol(merged_data)], na.rm = TRUE)  
print(row_sums)
summary(row_sums) 


# Remove columns where all values are 0 or NA
merged_data <- merged_data[, apply(merged_data, 2, function(x) any(!is.na(x)))]

# Ensure all column names are valid R variable names
colnames(merged_data) <- make.names(colnames(merged_data))

# Get the correctly formatted column names
all_columns <- colnames(merged_data)[11:ncol(merged_data)]

# Remove s__ prefix from species column names
colnames(merged_data)[11:ncol(merged_data)] <- gsub("^s__", "", colnames(merged_data)[11:ncol(merged_data)])
all_columns <- colnames(merged_data)[11:ncol(merged_data)]

merged_data <- na.omit(merged_data)

merged_data$Patient <- as.character(merged_data$Patient)

# Isolate data for Prevotella_timonensis 
df_prev <- merged_data %>%
  dplyr::select(Patient, Response, Prevotella_timonensis) %>%
  dplyr::filter(Response %in% c("DC", "No_DC"))

# Plot 
prevotella_plot <- ggplot(df_prev, aes(x = Response, y = Prevotella_timonensis, fill = Response)) +
  stat_summary(fun = "mean", geom = "bar", color = "black", width = 0.6) +
  geom_jitter(aes(color = Response), width = 0.15, size = 2, alpha = 0.8) +
  geom_text_repel(
    aes(label = Patient, color = Response),
    size = 3,
    box.padding = 0.1,
    point.padding = 0.2,
    segment.color = NA,
    show.legend = FALSE
  ) +
  geom_signif(
    comparisons      = list(c("DC", "No_DC")),
    test             = "wilcox.test",
    map_signif_level = function(p) paste0("p = ", signif(p, 3)),
    y_position = max(df_prev$Prevotella_timonensis, na.rm = TRUE) * 1.15,
    step_increase    = 0.12,
    color            = "black",
    size             = 0.6,
    textsize         = 4
  ) +
  scale_fill_manual(values = c("DC" = "lightblue", "No_DC" = "#FFC300")) +
  scale_color_manual(values = c("DC" = "blue", "No_DC" = "orange")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(
    title = "Prevotella timonensis",
    x = NULL,
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(size = 14, face = "bold.italic", hjust = 0.5),
    axis.title.y     = element_text(size = 14, margin = margin(r = 10)),
    axis.text        = element_text(size = 12),
    axis.text.x      = element_text(color = "black", size = 12),
    legend.position  = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line        = element_line(color = "black", linewidth = 0.6)
  )

print(prevotella_plot)


####### Baseline, Day 8 and Day 36 comparison of prevotella timonensis

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

tse <- tse_Metaphlan

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Day %in% c("1", "8", "36"))]
tse <- tse[, which(colData(tse)$Trial == "T563")]
tse <- tse[, which(colData(tse)$Response %in% c("No_DC", "DC"))]

tse <- mergeFeaturesByRank(tse, rank ="species", onRankOnly=TRUE)


# Remove "s__" prefix from species names
rowData(tse)$species <- gsub("^s__", "", rowData(tse)$species)


# Select only specific Prevotella species
#target_species <- c("Prevotella timonensis_A") #DADA2
target_species <- c("Prevotella_timonensis") #MetaPhlan
tse_prevotella <- tse[rowData(tse)$species %in% target_species, ]

# Aggregate species-level data
No_DC_rel <- tse_prevotella[, colData(tse_prevotella)$Response == "No_DC"]
tse_No_DC_species <- agglomerateByRank(No_DC_rel, rank = "species", onRankOnly = TRUE)

DC_rel <- tse_prevotella[, colData(tse_prevotella)$Response == "DC"]
tse_DC_species <- agglomerateByRank(DC_rel, rank = "species", onRankOnly = TRUE)

# Merge dataset at species level
combined_tse <- mergeFeaturesByRank(tse_prevotella, rank = "species", onRankOnly = TRUE)

# Update row names with Patient ID
col_data_tse <- colData(combined_tse)
rownames(col_data_tse) <- col_data_tse$Response_Patient_day
colData(combined_tse) <- col_data_tse

# Ensure numeric ordering of Patient_day for correct plotting sequence
colData(combined_tse)$Patient_day <- as.numeric(as.character(colData(combined_tse)$Patient_day))
combined_tse <- combined_tse[, order(colData(combined_tse)$Patient_day)]

# Convert Patient_day to factor for proper plotting order
colData(combined_tse)$Patient_day <- factor(colData(combined_tse)$Patient_day, 
                                            levels = sort(unique(colData(combined_tse)$Patient_day)))
                                   
# Order samples as specified
desired_order <- c("DC_30105_1", "DC_30105_8", "DC_30105_36", "DC_30111_1", "DC_30111_8", "DC_30111_36", "DC_30112_1", "DC_30112_8", "DC_30112_36",
                   "DC_30207_1", "DC_30207_8", "DC_30207_36", "DC_30210_1", "DC_30210_8", "DC_30210_36",
                   "No_DC_30104_1", "No_DC_30104_36", "No_DC_30110_1", "No_DC_30110_8", "No_DC_30110_36",
                   "No_DC_30209_1", "No_DC_30209_8", "No_DC_30209_36", "No_DC_30212_1", "No_DC_30212_8")
colData(combined_tse)$Response_Patient_day <- factor(colData(combined_tse)$Response_Patient_day, levels = desired_order)
combined_tse <- combined_tse[, order(colData(combined_tse)$Response_Patient_day)]

### Convert Data to ggplot-Compatible Format
df <- as.data.frame(assay(combined_tse, "relabundance"))
df$species <- rowData(combined_tse)$species

# Transform data to long format
df_long <- df %>%
  pivot_longer(cols = -species, names_to = "Sample", values_to = "Abundance")

# Merge with metadata
metadata <- as.data.frame(colData(combined_tse))
metadata$Sample <- rownames(metadata)

df_long <- dplyr::left_join(df_long, metadata, by = "Sample")

df_long$Sample <- factor(df_long$Sample, levels = desired_order)

grep("Prevotella", rowData(tse)$species, value = TRUE)

# DADA2
#species_colors <- c(
#  "Prevotella_timonensis" = "lightblue"
#)

# Metaphlan
species_colors <- c(
  "Prevotella_timonensis" = "lightblue"
)

patient_colors <- c(
  "30105" = "#c6dbef",
  "30111" = "#9ecae1",
  "30112" = "#6baed6",
  "30207" = "#3182bd",
  "30210" = "#08519c",
  "30104" = "#fdd0a2",
  "30110" = "#f16913",
  "30209" = "#cb181d",
  "30212" = "darkred"
)

# add colors to dataframe
df_long$Patient <- as.character(df_long$Patient)

df_long$color <- patient_colors[df_long$Patient]

## Plot
ggplot(df_long, aes(x = Sample, y = Abundance, fill = color)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_identity() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    legend.text = element_text(size = 14, face = "italic"),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    panel.background = element_blank(),
    panel.grid = element_blank()
  ) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = expression(bold(P.~timonensis ~ "abundance during therapy")),
    x = "",
    y = "Relative Abundance (%)",
    fill = "species"
  )




## Patient overview with fixed effect

library(readxl)
library(dplyr)
library(tidyr)

# Read data
fixed_effect_df <- read_excel(
  "Z:/PhD CGTG/Experiments/Metagenomics/human_urine/Fixed_effects.xlsx"
) %>%
  mutate(
    Patients = as.character(Patients),
    Age = as.character(Age)
  )


# Convert to long format
df_long <- fixed_effect_df %>%
  pivot_longer(
    cols = c(
      `Response`,
      `Age`,
      `Hospital`,
      `Primary Tumor`,
      `Histological subtype`,
      `Stage at diagnosis`
    ),
    names_to = "Variable",
    values_to = "Value"
  ) 

# Bin age groups
df_long <- df_long %>%
  mutate(
    Value = case_when(
      Variable == "Age" & Value >= 50 & Value < 60 ~ "50–60",
      Variable == "Age" & Value >= 60 & Value < 70 ~ "60–70",
      Variable == "Age" & Value >= 70 & Value <= 80 ~ "70–80",
      TRUE ~ as.character(Value)
    )
  )

# Annotation colors
annotation_colors <- c(
  # Response
  "No_DC" = "darkslategray2",
  "DC"  = "#FFDD44", 
  
  # Age 
  "50–60" = "#6baed6",
  "60–70" = "#3182bd",
  "70–80" = "#08519c",
  
  # Hospital
  "Mayo Clinic" = "darkorchid1",
  "Docrates Cancer Center" = "darkorchid4",
  
  # Primary Tumor
  "Epithelial Ovarian Cancer"        = "burlywood",
  "Fallopian Tube Cancer"  = "darkorange",
  "Primary Peritoneal Cancer"         = "darkorange3",
  
  # Histological subtype
  "Serous"        = "deeppink4",
  "High-grade serous"  = "deeppink3",
  "Mucinous carcinoma"         = "#df65b0",
  "Carcinosarcoma"         = "#d4b9da",
  
  # Stage at diagnosis
  "IIIA2"        = "#b2df8a",
  "IIIB"  = "#66c2a5",
  "IIIC"         = "#41ae76",
  "IVA"         = "#238b45",
  "IVB"         = "#013220"
)

# Order variables (controls plot & legend order)
df_long$Variable <- factor(
  df_long$Variable,
  levels = c(
    "Response",
    "Age",
    "Hospital",
    "Primary Tumor",
    "Histological subtype",
    "Stage at diagnosis"
  )
)

# Alphabetical row order
variable_order <- c("Response", "Age", "Hospital", "Primary Tumor", "Histological subtype", "Stage at diagnosis")
df_long$Variable <- factor(df_long$Variable, levels = variable_order)

# Sample order by increasing Microbiome present 
microbiome_order <- df_long %>%
  filter(Variable == "Response") %>%
  mutate(MicrobiomeNum = ifelse(Value == "DC", 1, 0)) %>%
  arrange(MicrobiomeNum) %>%
  pull(Patients)

df_long$Patients <- factor(df_long$Patients, levels = microbiome_order)



#  Plot 
ggplot() +
geom_tile(
  data = df_long %>% filter(Variable == "Response"),
  aes(x = Patients, y = Variable, fill = Value),
  color = "white", size = 0.3
) +
  scale_fill_manual(
    name = "Response",
    values = scales::alpha(
      c("No_DC" = "darkslategray2", "DC" = "#FFDD44"),
      0.60
    ),
    guide = guide_legend(order = 1, title.position = "top",
                         direction = "vertical")
  ) +
  
  ggnewscale::new_scale_fill() +
 
geom_tile(
  data = df_long %>% filter(Variable == "Age"),
  aes(x = Patients, y = Variable, fill = Value),
  color = "white", size = 0.3
) +
  scale_fill_manual(
    name = "Age",
    values = scales::alpha(
      c(
        "50–60" = "#6baed6",
        "60–70" = "#3182bd",
        "70–80" = "#08519c"
      ),
      0.60
    ),
    guide = guide_legend(order = 2, title.position = "top",
                         direction = "vertical")
  ) +
  ggnewscale::new_scale_fill() +
geom_tile(
  data = df_long %>% filter(Variable == "Hospital"),
  aes(x = Patients, y = Variable, fill = Value),
  color = "white", size = 0.3
) +
  scale_fill_manual(
    name = "Hospital",
    values = scales::alpha(
      c(
        "Mayo Clinic" = "darkorchid1",
        "Docrates Cancer Center" = "darkorchid4"
      ),
      0.60
    ),
    guide = guide_legend(order = 3, title.position = "top",
                         direction = "vertical")
  ) +
  
  ggnewscale::new_scale_fill() +
geom_tile(
  data = df_long %>% filter(Variable == "Primary Tumor"),
  aes(x = Patients, y = Variable, fill = Value),
  color = "white", size = 0.3
) +
  scale_fill_manual(
    name = "Primary Tumor",
    values = scales::alpha(
      c(
        "Epithelial Ovarian Cancer" = "burlywood",
        "Fallopian Tube Cancer" = "darkorange",
        "Primary Peritoneal Cancer" = "darkorange3"
      ),
      0.60
    ),
    guide = guide_legend(order = 4, title.position = "top",
                         direction = "vertical")
  ) +
  
  ggnewscale::new_scale_fill() +
geom_tile(
  data = df_long %>% filter(Variable == "Histological subtype"),
  aes(x = Patients, y = Variable, fill = Value),
  color = "white", size = 0.3
) +
  scale_fill_manual(
    name = "Histological subtype",
    values = scales::alpha(
      c(
        "Serous" = "deeppink4",
        "High-grade serous" = "deeppink3",
        "Mucinous carcinoma" = "#df65b0",
        "Carcinosarcoma" = "#d4b9da"
      ),
      0.60
    ),
    guide = guide_legend(order = 5, title.position = "top",
                         direction = "vertical")
  ) +
  
  ggnewscale::new_scale_fill() +
geom_tile(
  data = df_long %>% filter(Variable == "Stage at diagnosis"),
  aes(x = Patients, y = Variable, fill = Value),
  color = "white", size = 0.3
) +
  scale_fill_manual(
    name = "Stage at diagnosis",
    values = scales::alpha(
      c(
        "IIIA2" = "#b2df8a",
        "IIIB"  = "#66c2a5",
        "IIIC"  = "#41ae76",
        "IVA"   = "#238b45",
        "IVB"   = "#013220"
      ),
      0.60
    ),
    guide = guide_legend(order = 6, title.position = "top",
                         direction = "vertical")
  ) +
scale_y_discrete(limits = rev(levels(df_long$Variable))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1.1, hjust = 0.8, size = 11, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 11, face = "bold", angle = 0),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    
    legend.spacing.x = unit(0.6, "cm"),
    legend.box.spacing = unit(0.6, "cm")
  )
