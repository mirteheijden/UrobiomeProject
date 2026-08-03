############ Shotgun metagenomics 
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
tse_Metaphlan <- importMetaPhlAn("path/human_urine/Batch 4/MetaPhlan/metaphlan_db_meta4_combined_reports.txt", "path/human_urine/Batch 4/Puhti/metadata.tsv", package = "mia")

# change MI_ID rownames of the tse, to our sample_ID
col_data_Metaphlan <- colData(tse_Metaphlan)
rownames(col_data_Metaphlan) <- col_data_Metaphlan$Patient_day
colData(tse_Metaphlan) <- col_data_Metaphlan
col_data_Metaphlan <- as.data.frame(colData(tse_Metaphlan))


# Use species as rownames instead of strains
# Check for NAs in species column
sum(is.na(rowData(tse_Metaphlan)$species))

# Check for duplicate species names (can happen with SGBs)
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
# GG2
#tse <- tse_GG2
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
  y = "observed",      # y variable (numeric)
  x = "Patient",       # x variable (ordered categorical - patient names)
  colour_by = "Response",   # Colour by another variable
  point_size = 4
  #size_by = "observed" # Scale dot size by a numeric variable (e.g., "observed")
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Larger x-axis labels
        axis.text.y = element_text(size = 14),  # Larger y-axis labels
        axis.title.x = element_text(size = 18),  # Larger axis title
        axis.title.y = element_text(size = 18),  
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Bigger, bold title
        legend.text = element_text(size = 14),  # Bigger legend text
        legend.title = element_text(size = 16),  # Bigger legend title
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
  ) +
  labs(
    title = "Alpha diversity (Shannon)",
    x = "
    Patient", 
    y = expression(Richness[Observed]), 
    colour = "Response"  # Change legend title
  ) +
  scale_colour_manual(values = c("DC" = "blue", "No_DC" = "orange")) +
  guides(colour = guide_legend(override.aes = list(size = 6)))  # Adjust dot size in legend


library(dplyr)
library(ggplot2)
library(ggsignif)

# ── Build plotting data from colData ──────────────────────────────────────────
df_alpha <- as.data.frame(colData(tse)) %>%
  dplyr::select(Response, observed) %>%
  dplyr::filter(!is.na(Response), !is.na(observed)) %>%
  dplyr::mutate(
    Response = factor(Response, levels = c("DC", "No_DC"))
  )

# ── Colours ───────────────────────────────────────────────────────────────────
fill_cols   <- c("DC" = "blue", "No_DC" = "orange")
border_cols <- c("DC" = "#6E6E6E", "No_DC" = "#C07D10")

# ── Plot ──────────────────────────────────────────────────────────────────────
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
    width         = 0.12,        # narrow inner boxplot
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
## Response ###

# GG2
#tse <- tse_GG2
tse <- tse_Metaphlan
tse <- tse[, which(colData(tse)$Trial == "T563")]

# species level
tse <- mergeFeaturesByRank(tse, rank ="species", onRankOnly=TRUE)

tse <- tse[, which(colData(tse)$Response %in% c("DC", "No_DC"))]

# Metaphlan (check availability altExpNames(tse))
#tse_rel <- transformAssay(tse, assay.type = "relabundance", method = "relabundance")

# Extract day 1 samples
tse_rel <- tse
tse_day1 <- tse_rel[, rownames(subset(colData(tse_rel), Day == "1"))]

# Or Include only species present in at least 2 or more patients
species_to_keep <- rowSums(assay(tse_day1, "relabundance") > 0) >= 2
tse_filtered <- tse_day1[species_to_keep, ]

#assay_names <- assayNames(tse_filtered)
#print(assay_names)
tse_day1_rel_abund_assay <- assays(tse_filtered)$relabundance
bray_curtis_dist <- vegan::vegdist(t(tse_day1_rel_abund_assay), method = "bray")
#install.packages("ecodist")
#library(ecodist)
bray_curtis_pcoa <- ecodist::pco(bray_curtis_dist)
#bray_curtis_pcoa$vectors
#PCoA1 and PCoA2
bray_curtis_pcoa_df <- data.frame(pcoa1 = bray_curtis_pcoa$vectors[,1], 
                                  pcoa2 = bray_curtis_pcoa$vectors[,2])
# Create a plot
# Adds the variable you want to use for coloring to the data frame
bray_curtis_pcoa_df <- cbind(bray_curtis_pcoa_df,
                             patient_status = colData(tse_filtered)$Response,
                             label = colData(tse_filtered)$Patient)


# For DC vs No_DC
bray_curtis_plot <- ggplot(data = bray_curtis_pcoa_df, 
                           aes(x = pcoa1, y = pcoa2, color = patient_status)) +
  geom_point() +  # Points for each sample
  geom_text(aes(label = label), size = 3, vjust = 1.5, show.legend = FALSE) +  # Labels for each point (Patient names)
  
  # Add circles around the groups with colored outlines and no fill
  geom_mark_ellipse(aes(group = patient_status), 
                    fill = NA,  # No fill
                    linetype = "solid",  # Outline type
                    show.legend = FALSE) +  # Hide legend for the ellipses
  
  labs(x = "PC1",
       y = "PC2", 
       colour = "Response",
       title = "Beta diversity between groups") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),  # Larger x-axis labels
        axis.text.y = element_text(size = 14),  # Larger y-axis labels
        axis.title.x = element_text(size = 14),  # Larger axis title
        axis.title.y = element_text(size = 14),  
        plot.title = element_text(size = 12, face = "bold"),  # Bigger, bold title
        legend.text = element_text(size = 14),  # Bigger legend text
        legend.title = element_text(size = 16),  # Bigger legend title
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
        panel.background = element_rect(fill = "white", color = NA),  # White background
        axis.line = element_line(color = "black", linewidth = 0.4)
        ) +
  scale_x_continuous(expand = expansion(mult = 0.3)) + # Adjust horizontal space
  scale_y_continuous(expand = expansion(mult = 0.3)) + # Adjust vertical space
  scale_color_manual(values = c("DC" = "blue", "No_DC" = "orange")) + # Define custom colors
  guides(colour = guide_legend(override.aes = list(size = 6)))  # Adjust dot size in legend

# Display the plot
bray_curtis_plot




#########################################################################
### Differential abundance analysis
### ANCOM-BC - GG2
library(ANCOMBC)

tse <- tse_Metaphlan

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Trial == "T563")]
tse <- tse[, which(colData(tse)$Response %in% c("No_DC", "DC"))]

# Convert to phyloseq object
pseq <- convertToPhyloseq(tse, assay.type = "relabundance")
pseq_genus <- phyloseq::tax_glom(pseq, taxrank = "genus")

# ── Replace row names with genus names ────────────────────────────────────────
genus_names <- as.character(tax_table(pseq_genus)[, "genus"])
taxa_names(pseq_genus) <- genus_names

# Perform ANCOMBC analysis, comparing DC vs No_DC
out_gg2_genus = ancombc(
  data = pseq_genus, 
  formula = "Response", 
  p_adj_method = "fdr", 
  lib_cut = 0, 
  group = "Response", 
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.01, 
  global = TRUE
)

# Extract results for No_DC (which includes log fold changes for the comparison between DC and No_DC)
res <- out_gg2_genus$res

# Compute log fold change for DC vs No_DC
results <- data.frame(
  taxon = res$lfc$taxon,
  log_fold_change = -res$lfc$ResponseNo_DC,
  q_value = as.numeric(as.character(res$q_val$ResponseNo_DC)),
  differentially_abundant = res$diff_abn$ResponseNo_DC
)

# Filter for significant results
significant_results <- results[results$differentially_abundant == TRUE, ]
significant_results <- significant_results[significant_results$q_value < 0.05,]
significant_results <- significant_results[significant_results$log_fold_change < 0 | significant_results$log_fold_change > 0,]

# Print the significant results
print(significant_results)

## Bar plot

# Extract the row names and genus column from rowData(tse)
row_data <- rowData(tse)
row_names_to_genus <- data.frame(
  row_name = rownames(row_data),
  genus = row_data$genus
)

# Create a lookup table for easy matching
lookup_table <- setNames(row_names_to_genus$genus, row_names_to_genus$row_name)

# Replace values in the 'taxon' column of significant_results
significant_results <- significant_results %>%
  mutate(taxon = ifelse(taxon %in% names(lookup_table), 
                        lookup_table[taxon], 
                        taxon))

# Add Response column based on log_fold_change direction
significant_results$Response <- ifelse(significant_results$log_fold_change > 0, "DC", "No_DC")

# Remove rows where 'taxon' is empty or NA
significant_results <- significant_results[!is.na(significant_results$taxon) & significant_results$taxon != "", ]
# Remove "g__" prefix
significant_results$taxon <- gsub("^g__", "", significant_results$taxon)
# Remove rows where "taxon" is empty or NA
significant_results <- significant_results[significant_results$taxon != "" & !is.na(significant_results$taxon), ]

# Assign colors based on Response only
significant_results$highlight <- ifelse(significant_results$Response == "DC", "blue", "orange")

# Ensure 'highlight' is a factor
significant_results$highlight <- factor(significant_results$highlight, 
                                        levels = c("blue", "orange"))

# Create the bar plot with custom colors
bar_plot <- ggplot(significant_results, aes(x = reorder(taxon, log_fold_change), 
                                            y = log_fold_change, fill = highlight)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  scale_fill_manual(values = c("blue" = "blue", "orange" = "orange"),
                    labels = c("blue" = "DC", 
                               "orange" = "No_DC")) +
  coord_flip() +
  labs(title = "", x = "Taxon
       ", y = "Log Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),  
        plot.title = element_text(size = 16, face = "bold", hjust = 0),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(title = ""))

print(bar_plot)





#########################################################################
### Differential abundance analysis
### ANCOM-BC - GG2
library(ANCOMBC)

tse <- tse_Metaphlan

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Trial == "T563")]
tse <- tse[, which(colData(tse)$Response %in% c("No_DC", "DC"))]

# Convert to phyloseq object
pseq <- convertToPhyloseq(tse, assay.type = "relabundance")
pseq_species <- phyloseq::tax_glom(pseq, taxrank = "species")

# ── Replace row names with species names ────────────────────────────────────────
species_names <- as.character(tax_table(pseq_species)[, "species"])
taxa_names(pseq_species) <- species_names


# Perform ANCOMBC analysis, comparing DC vs No_DC
out_gg2_species = ancombc(
  data = pseq_species, 
  formula = "Response", 
  p_adj_method = "fdr", 
  lib_cut = 0, 
  group = "Response", 
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.01, 
  global = TRUE
)

# Extract results for No_DC (which includes log fold changes for the comparison between DC and No_DC)
res <- out_gg2_species$res

# Compute log fold change for DC vs No_DC
results <- data.frame(
  taxon = res$lfc$taxon,
  log_fold_change = -res$lfc$ResponseNo_DC,
  q_value = as.numeric(as.character(res$q_val$ResponseNo_DC)),
  differentially_abundant = res$diff_abn$ResponseNo_DC
)

# Filter for significant results
significant_results <- results[results$differentially_abundant == TRUE, ]
significant_results <- significant_results[significant_results$q_value < 0.05,]
significant_results <- significant_results[significant_results$log_fold_change < 0 | significant_results$log_fold_change > 0,]

# Print the significant results
print(significant_results)

## Bar plot

# Extract the row names and species column from rowData(tse)
row_data <- rowData(tse)
row_names_to_species <- data.frame(
  row_name = rownames(row_data),
  species = row_data$species
)

# Create a lookup table for easy matching
lookup_table <- setNames(row_names_to_species$species, row_names_to_species$row_name)

# Replace values in the 'taxon' column of significant_results
significant_results <- significant_results %>%
  mutate(taxon = ifelse(taxon %in% names(lookup_table), 
                        lookup_table[taxon], 
                        taxon))

# Add Response column based on log_fold_change direction
significant_results$Response <- ifelse(significant_results$log_fold_change > 0, "DC", "No_DC")

# Remove rows where 'taxon' is empty or NA
significant_results <- significant_results[!is.na(significant_results$taxon) & significant_results$taxon != "", ]
# Remove "s__" prefix
significant_results$taxon <- gsub("^s__", "", significant_results$taxon)
# Remove rows where "taxon" is empty or NA
significant_results <- significant_results[significant_results$taxon != "" & !is.na(significant_results$taxon), ]


# Assign colors: gray for all, blue for DC with Prevotella, #FFC300 for No_DC with Prevotella
significant_results$highlight <- "gray"
significant_results$highlight[grepl("^Prevotella", significant_results$taxon, ignore.case = TRUE) & significant_results$Response == "DC"] <- "lightblue"
significant_results$highlight[grepl("^Prevotella", significant_results$taxon, ignore.case = TRUE) & significant_results$Response == "No_DC"] <- "#FFC300"

# Ensure 'highlight' is a factor with the desired order
significant_results$highlight <- factor(significant_results$highlight, 
                                        levels = c("lightblue", "#FFC300", "gray"))


# Create the bar plot with custom colors
bar_plot <- ggplot(significant_results, aes(x = reorder(taxon, log_fold_change), 
                                            y = log_fold_change, fill = highlight)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  scale_fill_manual(values = c("lightblue" = "lightblue", "#FFC300" = "#FFC300", "gray" = "gray"),
                    labels = c("lightblue" = "Prevotella in DC", 
                               "#FFC300" = "Prevotella in No_DC", 
                               "gray" = "Other Taxa")) +
  coord_flip() +
  labs(title = "", x = "Taxon
       ", y = "Log Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 14),  # Larger x-axis labels
        axis.text.y = element_text(size = 14),  # Larger y-axis labels
        axis.title.x = element_text(size = 18),  # Larger axis title
        axis.title.y = element_text(size = 18),  
        plot.title = element_text(size = 16, face = "bold", hjust = 0),  # Bigger, bold title
        legend.text = element_text(size = 14),  # Bigger legend text
        legend.title = element_text(size = 16),  # Bigger legend title
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
  ) +
  guides(fill = guide_legend(title = ""))

# Print the plot
print(bar_plot)





### LinDA 
#install.packages("MicrobiomeStat")
library(MicrobiomeStat)
####### species

#tse <- tse_GG2
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
#significant_results$highlight <- factor(significant_results$highlight, 
#                                        levels = c("lightblue", "#FFC300"))

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
                    labels = c("lightblue" = "DC", 
                               "#FFC300"   = "No_DC")) +
  scale_color_identity() +
  coord_flip() +
  labs(title = "Differential abundant species 
      (shotgun metagenomics)", x = "Taxon
       ", y = "Log Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14, 
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

# Counts of P. timonensis
# colors "lightgreen", "deepskyblue2"
assay(tse_spec)["s__Prevotella_timonensis", ]



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


# color codes: G = b_l_ue, o_r_a_nge and S = l_i_g_htblue, #F_F_C_300
# Create the bar plot with the 'Response' column for color mapping
bar_plot <- ggplot(significant_results, aes(x = reorder(taxon, log2FoldChange), y = log2FoldChange, fill = Response)) +
  geom_bar(stat = "identity", show.legend = TRUE) +  # Show legend
  scale_fill_manual(values = c("DC" = "blue", "No_DC" = "orange")) +  # blue for DC, #DFFC300 for No_DC (positive log2FoldChange)
  coord_flip() +  # Flip coordinates to make taxa names readable
  labs(title = "Differential abundant genus 
  (shotgun metagenomics)", x = "Taxon", y = "Log Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 14),  # Larger x-axis labels
        axis.text.y = element_text(size = 14),  # Larger y-axis labels
        axis.title.x = element_text(size = 18),  # Larger axis title
        axis.title.y = element_text(size = 18),  
        plot.title = element_text(size = 20, face = "bold", hjust = 0),  # Bigger, bold title
        legend.text = element_text(size = 14),  # Bigger legend text
        legend.title = element_text(size = 16),  # Bigger legend title
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank(),  # Remove minor grid lines
  ) +
  guides(fill = guide_legend(title = "", 
                             labels = c("DC", "No_DC")))  # Custom labels for the legend

# Print the plot
print(bar_plot)





### Wilcoxon barplots
################## Response #########################
library(writexl)
library(ggsignif)

#tse <- tse_GG2
tse <- tse_Metaphlan

tse <- tse[, which(colData(tse)$Day == "1")]
#tse <- tse[, which(colData(tse)$Trial == "T563")]
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
tse_assay <- t(assay(tse, "relabundance"))  # Use the manually transformed relative abundance data

# Merge metadata with transposed assay data
merged_data <- merge(metadata_df, as.data.frame(tse_assay), by.x = "row.names", by.y = "row.names", all.x = TRUE)

# DC first in boxplot
merged_data$Response <- factor(merged_data$Response, levels = c("DC", "No_DC"))

# Directory to save the combined PDF
output_dir <- "path/human_urine/RESULTS"
output_pdf <- file.path(output_dir, "Wilcoxon_species_MetaPhlan4.pdf")


# Initialize a list to store all plots
plot_list <- list()

# Use all columns from the merged_data (or merged_data)
colnames(merged_data) <- gsub("\\[|\\]", "", colnames(merged_data))
all_columns <- colnames(merged_data)[11:ncol(merged_data)]
print(all_columns)
## Change the 11:ncol from where the first species starts!!

#check relative abundance
row_sums <- rowSums(merged_data[, 11:ncol(merged_data)], na.rm = TRUE)  
print(row_sums)
summary(row_sums)  # Check if values are around 1


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


library(ggplot2)
library(ggrepel)

### PDF
plot_list <- list()

library(writexl) 
results_df <- data.frame(Taxon = character(), P_value = numeric(), stringsAsFactors = FALSE)

for (col in all_columns) {
  
  wilcox_result <- wilcox.test(as.formula(paste(col, "~ Response")),
                               data = merged_data,  
                               subset = merged_data$Response %in% c("DC", "No_DC"),
                               exact = FALSE)
  
  results_df <- rbind(results_df, data.frame(Taxon = col, P_value = wilcox_result$p.value))
  
  y_min <- min(merged_data[[col]], na.rm = TRUE)
  y_max <- max(merged_data[[col]], na.rm = TRUE)
  y_range <- y_max - y_min
  y_limits <- c(y_min - 0.1 * y_range, y_max + 0.3 * y_range)
  
  p <- ggplot(merged_data, aes_string(x = "Response", y = col, fill = "Response")) +
    stat_summary(fun = "mean", geom = "bar", color = "black", width = 0.6) +
    geom_jitter(aes_string(color = "Response"), width = 0.15, size = 2, alpha = 0.8) +
    geom_text_repel(
      aes_string(label = "Patient", color = "Response"),
      size = 3,
      segment.size = 0,
      box.padding = 0.1,
      point.padding = 0.2,
      segment.color = NA,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = c("DC" = "lightblue", "No_DC" = "#FFC300")) +
    scale_color_manual(values = c("DC" = "blue", "No_DC" = "orange")) +
    labs(
      title = col,
      x = NULL,
      y = "Relative Abundance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      axis.line = element_line(color = "black", linewidth = 0.6)
    ) +
    coord_cartesian(ylim = y_limits, clip = "off") +
    geom_signif(
      test = "wilcox.test",
      comparisons = list(c("DC", "No_DC")),
      map_signif_level = function(p) paste0("p = ", signif(p, digits = 3)),
      step_increase = 0.1,
      color = "black",
      size = 0.7,
      textsize = 4
    )
  
  plot_list[[length(plot_list) + 1]] <- p
  
}  

# Print Prevotella_timonensis plot
#prevotella_plot <- plot_list[[which(all_columns == "Prevotella_timonensis")]]
print(prevotella_plot)


# Save all plots to PDF (optional)
pdf(output_pdf, width = 8, height = 6)
for (plt in plot_list) {
  print(plt)
}
dev.off()

# Optional: Adjust p-values for multiple testing
results_df$Adjusted_P <- p.adjust(results_df$P_value, method = "BH")

# Write results to Excel
write_xlsx(results_df, "path/human_urine/RESULTS/wilcoxon_barplots_species_MetaPhlan4.xlsx")



# Boxplot Prevotella timonensis only
library(writexl)
library(ggsignif)
library(ggrepel)

#tse <- tse_GG2
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
tse_assay <- t(assay(tse, "relabundance"))  # Use the manually transformed relative abundance data

# Merge metadata with transposed assay data
merged_data <- merge(metadata_df, as.data.frame(tse_assay), by.x = "row.names", by.y = "row.names", all.x = TRUE)

# DC first in boxplot
merged_data$Response <- factor(merged_data$Response, levels = c("DC", "No_DC"))

# Directory to save the combined PDF
output_dir <- "path/human_urine/RESULTS"
output_pdf <- file.path(output_dir, "Wilcoxon_species_MetaPhlan4.pdf")


# Initialize a list to store all plots
plot_list <- list()

# Use all columns from the merged_data (or merged_data)
colnames(merged_data) <- gsub("\\[|\\]", "", colnames(merged_data))
all_columns <- colnames(merged_data)[11:ncol(merged_data)]
print(all_columns)
## Change the 16:ncol from where the first species starts!!

#check relative abundance
row_sums <- rowSums(merged_data[, 16:ncol(merged_data)], na.rm = TRUE)  
print(row_sums)
summary(row_sums)  # Check if values are around 1


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

# ── Isolate data for Prevotella_timonensis ────────────────────────────────────
df_prev <- merged_data %>%
  dplyr::select(Patient, Response, Prevotella_timonensis) %>%
  dplyr::filter(Response %in% c("DC", "No_DC"))

# ── Plot ──────────────────────────────────────────────────────────────────────
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
    plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
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








####### Baseline comparison of prevotella timonensis

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# Filter data: Select Trial T563 and Day 1 samples
#tse <- tse_phylo # use capital S and G
tse <- tse_Metaphlan #use small s and g
tse <- tse[, colData(tse)$Day == "1"]
tse <- tse[, colData(tse)$Trial == "T563"]

# Remove "s__" prefix from species names
rowData(tse)$species <- gsub("^s__", "", rowData(tse)$species)

# Filter samples with known disease outcome (No_DC or DC)
tse <- tse[, colData(tse)$Response %in% c("No_DC", "DC")]

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

### Filter for DC Group
#combined_tse_DC <- combined_tse[, colData(combined_tse)$Response == "DC"]

# Remove "_1" suffix from sample names
colData(combined_tse)$Response_Patient_day <- gsub("_1$", "", colData(combined_tse)$Response_Patient_day)
rownames(colData(combined_tse)) <- colData(combined_tse)$Response_Patient_day

# Order samples as specified
#desired_order <- c("DC_30111", "DC_30112", "DC_30210", "DC_30105", "DC_30207", "No_DC_30104", "No_DC_30209", "No_DC_30110", "No_DC_30212")
#colData(combined_tse)$Response_Patient_day <- factor(colData(combined_tse)$Response_Patient_day, levels = desired_order)
#combined_tse <- combined_tse[, order(colData(combined_tse)$Response_Patient_day)]

### Convert Data to ggplot-Compatible Format
df <- as.data.frame(assay(combined_tse, "metaphlan"))
df$species <- rowData(combined_tse)$species

# Transform data to long format
df_long <- df %>%
  pivot_longer(cols = -species, names_to = "Sample", values_to = "Abundance")

# Merge with metadata
metadata <- as.data.frame(colData(combined_tse))
metadata$Sample <- rownames(metadata)
df_long <- merge(df_long, metadata, by = "Sample")

# Remove "_1" from sample names (if still present)
df_long$Sample <- gsub("_1$", "", df_long$Sample)

# Set sample order explicitly
#df_long$Sample <- factor(df_long$Sample, levels = desired_order)

grep("Prevotella", rowData(tse)$species, value = TRUE)

# DADA2
#species_colors <- c(
#  "Prevotella_timonensis" = "lightblue"
#)

# Metaphlan
species_colors <- c(
  "Prevotella_timonensis" = "lightblue"
)

### Generate ggplot2 Stacked Bar Plot
ggplot(df_long, aes(x = Sample, y = Abundance, fill = species)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = species_colors) +
  theme_minimal(base_size = 14) +  # Increase overall text size
  theme(
    axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),  # Add space to Y-axis title
    legend.text = element_text(size = 14, face = "italic"),  # Italic legend text
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    panel.background = element_blank(),  # Remove background
    panel.grid = element_blank()  # Remove gridlines
  ) +
  scale_y_continuous(labels = scales::percent) +  
  labs(title = "", x = "", y = "Relative Abundance (%)", fill = "species")






####### Baseline, Day 8 and Day 36 comparison of prevotella timonensis

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

tse <- tse_Metaphlan #already in relabundance style

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

### Filter for DC Group
#combined_tse_DC <- combined_tse[, colData(combined_tse)$Response == "DC"]

# Remove "_1" suffix from sample names
#colData(combined_tse)$Response_Patient_day <- gsub("_1$", "", colData(combined_tse)$Response_Patient_day)
#rownames(colData(combined_tse)) <- colData(combined_tse)$Response_Patient_day

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



# Remove "_1" from sample names (if still present)
#df_long$Sample <- gsub("_1$", "", df_long$Sample)

# Set sample order explicitly
#df_long$Sample <- factor(df_long$Sample, levels = desired_order)

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
  scale_fill_identity() +   # <-- IMPORTANT
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


# Write results to Excel
write_xlsx(df_long, "path/human_urine/RESULTS/T563_shotgun_metagenomics/df_long_metaphlan.xlsx")








## Line plot

library(dplyr)
library(ggplot2)

# Clean and prepare data
df_plot <- df_long %>%
  mutate(
    Day = as.factor(Day),
    Day = factor(Day, levels = c("1", "8", "36")),
    Patient_label = paste0(Response, "_", Patient)
  )

# Define colors (blue = DC, orange = No_DC)
patient_colors <- c(
  "DC" = "#1f77b4",       # blue
  "No_DC" = "#ff7f0e"     # orange
)

# Plot
ggplot(df_plot, aes(x = Day, y = Abundance, group = Patient_label)) +
  
  # Lines per patient
  geom_line(aes(color = Response), linewidth = 1, alpha = 0.7) +
  
  # Points per timepoint
  geom_point(aes(shape = Day, color = Response), size = 3, stroke = 1.2) +
  
  # Colors
  scale_color_manual(values = patient_colors, name = "Response") +
  
  # Shapes for days
  scale_shape_manual(
    values = c("1" = 1, "8" = 17, "36" = 16),
    name = "Day"
  ) +
  
  # Labels
  labs(
    title = "Prevotella timonensis during therapy",
    x = "Day",
    y = "Relative abundance (%)"
  ) +
  
  # Convert to %
  scale_y_continuous(labels = scales::percent) +
  
  # Style (match your previous plot)
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", linewidth = 0.8)
  )


# New ilne plot test
library(RColorBrewer)

# Creat color paletter per patient
# Get unique patients grouped by response
dc_patients <- unique(df_plot$Patient[df_plot$Response == "DC"])
nodc_patients <- unique(df_plot$Patient[df_plot$Response == "No_DC"])

# Blue shades for DC (5)
dc_colors <- colorRampPalette(brewer.pal(9, "Blues"))(length(dc_patients))

# Orange/red shades for No_DC (3)
nodc_colors <- colorRampPalette(brewer.pal(9, "OrRd"))(length(nodc_patients))

# Combine into named vector
patient_colors <- c(
  setNames(dc_colors, dc_patients),
  setNames(nodc_colors, nodc_patients)
)

# Assign unique shapes per patient
shape_values <- c(16, 17, 15, 3, 7, 8, 18, 4)  # 8 distinct symbols

patient_shapes <- setNames(
  shape_values[1:length(unique(df_plot$Patient))],
  unique(df_plot$Patient)
)

# Create group legend labels
df_plot <- df_plot %>%
  mutate(
    Patient = factor(Patient),
    Respons = factor(Response, levels = c("DC", "No_DC"))
    )

# plot
ggplot(df_plot, aes(x = Day, y = Abundance, group = Patient)) +
  
  geom_line(aes(color = Patient), linewidth = 1, alpha = 0.8) +
  
  geom_point(
    aes(color = Patient, shape = Patient),
    size = 3, stroke = 1.2
  ) +
  
  # Apply custom colors & shapes
  scale_color_manual(
    values = patient_colors,
    labels = function(x) {
      ifelse(
        x %in% dc_patients,
        paste0("Patient (DC): ", x),
        paste0("Patient (No_DC): ", x)
      )
    }
  ) +
  
  scale_shape_manual(
    values = patient_shapes,
    labels = function(x) {
      ifelse(
        x %in% dc_patients,
        paste0("Patient (DC): ", x),
        paste0("Patient (No_DC): ", x)
      )
    }
  ) +
  
  labs(
    title = "Prevotella timonensis over time",
    x = "Day",
    y = "Relative abundance (%)",
    color = "Patients",
    shape = "Patients"
  ) +
  
  scale_y_continuous(labels = scales::percent) +
  
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.title.x = element_text(face = "bold", size = 14),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 11),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "black", linewidth = 0.8)
  )








# HEATMAP Significant ANCOM-BC spcies

tse <- tse_GG2

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Trial == "T563")]
tse <- tse[, which(colData(tse)$Response %in% c("No_DC", "DC"))]
# Remove the "s__" prefix from the "species" column in rowData
rowData(tse)$species <- gsub("^s__", "", rowData(tse)$species)

# Replace all spaces and hyphens with underscores in row names of tse_species
rownames(tse_species) <- gsub("[ -]", "_", rownames(tse_species))

# Transform count assay to relative abundances
tse_rel <- transformAssay(tse,
                          assay.type = "counts",
                          method = "relabundance")

# Convert to phyloseq object
pseq <- makePhyloseqFromTreeSummarizedExperiment(tse_rel)
pseq_species <- phyloseq::tax_glom(pseq, taxrank = "species")

# Perform ANCOMBC analysis, comparing DC vs No_DC
out_gg2_species = ancombc(
  phyloseq = pseq_species, 
  formula = "Response", 
  p_adj_method = "fdr", 
  lib_cut = 0, 
  group = "Response", 
  struc_zero = TRUE, 
  neg_lb = TRUE, 
  tol = 1e-5, 
  max_iter = 100, 
  conserve = TRUE, 
  alpha = 0.05, 
  global = TRUE
)

# Extract results for No_DC (which includes log fold changes for the comparison between DC and No_DC)
res <- out_gg2_species$res

# Compute log fold change for DC vs No_DC
# Flip the log fold change if needed, since currently it is for No_DC as the numerator
results <- data.frame(
  taxon = res$lfc$taxon,
  # Flip the fold change calculation: we assume the log fold change from No_DC to DC
  log_fold_change = -res$lfc$ResponseNo_DC,  # Reverse the sign for the log fold change
  q_value = as.numeric(as.character(res$q_val$ResponseNo_DC)),
  differentially_abundant = res$diff_abn$ResponseNo_DC
)


# Filter for significant results
significant_results <- results[results$differentially_abundant == TRUE, ]

# Filter for significant results based on q-value threshold
significant_results <- significant_results[significant_results$q_value < 0.05,]

# Filter for significant results based on log_fold_change thresholds
significant_results <- significant_results[significant_results$log_fold_change < 0 | significant_results$log_fold_change > 0,]

# Print the significant results
print(significant_results)

## bar plot

# Extract the row names and species column from rowData(tse_rel)
row_data <- rowData(tse_rel)
row_names_to_species <- data.frame(
  row_name = rownames(row_data),
  species = row_data$species
)

# Create a lookup table for easy matching
lookup_table <- setNames(row_names_to_species$species, row_names_to_species$row_name)

# Replace values in the 'taxon' column of significant_results
significant_results <- significant_results %>%
  mutate(taxon = ifelse(taxon %in% names(lookup_table), 
                        lookup_table[taxon], 
                        taxon))

# Create a new column for color mapping based on log_fold_change (blue for DC, #DFFC300 for No_DC)
significant_results$highlight <- ifelse(significant_results$log_fold_change > 0, "blue", "#orange")

# Create a new column to represent the Response category (DC = blue, No_DC = #DFFC300)
significant_results$Response <- ifelse(significant_results$log_fold_change > 0, "DC", "No_DC")

# Remove rows where 'taxon' is empty or NA
significant_results <- significant_results[!is.na(significant_results$taxon) & significant_results$taxon != "", ]
# Remove "g__" prefix
significant_results$taxon <- gsub("^g__", "", significant_results$taxon)
# Remove rows where "taxon" is empty or NA
significant_results <- significant_results[significant_results$taxon != "" & !is.na(significant_results$taxon), ]

# Clean up taxon names in significant_results
significant_results$taxon <- gsub("\\s+", "_", significant_results$taxon)  # Replace spaces with underscores
significant_results$taxon <- gsub("__", "_", significant_results$taxon)  # Remove double underscores if present
significant_results$taxon <- trimws(significant_results$taxon)  # Trim leading/trailing spaces


## Merged_data fr
tse <- tse_GG2

# Filter samples based on metadata
tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Trial == "T563")]
tse <- tse[, which(colData(tse)$Response %in% c("No_DC", "DC"))]

# Remove the "s__" prefix from the "species" column in rowData
rowData(tse)$species <- gsub("^s__", "", rowData(tse)$species)

# Extract species level
tse_species <- mergeFeaturesByRank(tse, rank = "species", onRankOnly = TRUE)

# Replace spaces and hyphens with underscores in row names
rownames(tse_species) <- gsub("[ -]", "_", rownames(tse_species))

# --------- MANUAL RELATIVE ABUNDANCE TRANSFORMATION ---------

# Get the counts matrix
counts_matrix <- assay(tse_species, "counts")

# Normalize by dividing each value by the column (sample) sum
rel_abundance_matrix <- sweep(counts_matrix, 2, colSums(counts_matrix), FUN = "/")

# Store as a new assay
assay(tse_species, "relabundance") <- rel_abundance_matrix

# Assign tse_rel to the transformed tse_species
tse_rel <- tse_species

tse_prevotella <- tse_rel[grepl("Prevotella", rowData(tse_rel)$species, ignore.case = TRUE), ]

# Convert colData to a data frame
metadata_df <- as.data.frame(colData(tse_rel))

# Transpose assay data
tse_assay <- t(assay(tse_rel, "relabundance"))  # Use the relative abundance data

# Merge metadata with transposed assay data
merged_data <- merge(metadata_df, as.data.frame(tse_assay), by.x = "row.names", by.y = "row.names", all.x = TRUE)



## Combine and make plot 
# Extract only the significant taxa from ANCOM-BC results
significant_taxa <- significant_results$taxon

# Clean up formatting: remove leading/trailing spaces and ensure consistent case
significant_taxa <- trimws(significant_taxa)  # Remove any extra spaces
significant_taxa <- tolower(significant_taxa) # Convert to lowercase for matching

# Get all columns that contain "Prevotella" in merged_data
prevotella_columns <- grep("Prevotella", colnames(merged_data), value = TRUE, ignore.case = TRUE)

# Clean column names to match formatting of significant_taxa
cleaned_colnames <- tolower(trimws(prevotella_columns))  # Normalize formatting

# Find the intersection between cleaned column names and significant taxa
significant_prevotella <- prevotella_columns[cleaned_colnames %in% significant_taxa]

# Debugging step: Print matched taxa
print(significant_prevotella)

# If no significant Prevotella taxa exist, stop with a message
if (length(significant_prevotella) == 0) {
  stop("No significant Prevotella taxa found in ANCOM-BC results.")
}

# Convert merged_data to long format for ggplot, using only significant taxa
heatmap_data <- merged_data %>%
  select(Response_Patient_day, all_of(significant_prevotella)) %>%
  pivot_longer(cols = -Response_Patient_day, names_to = "Taxon", values_to = "Abundance")

# Get x-axis positions for the first 5 and last 4 samples
first_5_x_min <- 1
first_5_x_max <- 5

last_4_x_min <- (length(unique(heatmap_data$Response_Patient_day)) - 3)
last_4_x_max <- length(unique(heatmap_data$Response_Patient_day))

# Create the heatmap
heatmap_plot <- ggplot(heatmap_data, aes(x = Response_Patient_day, y = reorder(Taxon, -Abundance), fill = Abundance)) +
  geom_tile(color = "white") +  # White borders for better separation
  scale_fill_gradientn(colors = c("lightyellow", "orange", "red"), 
                       values = scales::rescale(c(0, 0.001, 0.0025, 0.005, max(heatmap_data$Abundance))),
                       limits = c(0, max(heatmap_data$Abundance)),  
                       na.value = "white") +  
  labs(title = "Prevotella species", 
       x = "", 
       y = "", 
       fill = "Rel. Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  
        axis.text.y = element_text(size = 10), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  # Add black rectangle around the first 5 samples
  annotate("rect", xmin = first_5_x_min - 0.5, xmax = first_5_x_max + 0.5, 
           ymin = 0.5, ymax = length(unique(heatmap_data$Taxon)) + 0.5, 
           color = "black", fill = NA, linewidth = 1.2) +
  
  # Add black rectangle around the last 4 samples
  annotate("rect", xmin = last_4_x_min - 0.5, xmax = last_4_x_max + 0.5, 
           ymin = 0.5, ymax = length(unique(heatmap_data$Taxon)) + 0.5, 
           color = "black", fill = NA, linewidth = 1.2)

# Print the heatmap
print(heatmap_plot)










# HEATMAP all species
##TSE from GG2
tse_GG2 <- importQIIME2("path/human_urine/Batch 2/Puhti/RESULTS_taxprofiler/Snakemake_GG2/counts.qza", taxonomy = "path/human_urine/Batch 2/Puhti/RESULTS_taxprofiler/Snakemake_GG2/taxonomy.qza", sampleMetaFile="path/human_urine/Batch 2/Puhti/config/metadata_GG2.tsv")
# change rownames of the tse, to our sample_ID
col_data_GG2 <- colData(tse_GG2)
rownames(col_data_GG2) <- col_data_GG2$Patient_day
colData(tse_GG2) <- col_data_GG2
col_data_GG2 <- as.data.frame(colData(tse_GG2))

################## Response #########################

filter_species <- function(tse_subset) {
  genus_counts <- rowSums(assay(tse_subset) > 0)  # Count non-zero occurrences
  tse_subset_filtered <- tse_subset[genus_counts >= 0, ]  # Keep genera in ≥1 samples
  return(tse_subset_filtered)
}

# Load data
tse <- tse_GG2

tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Trial == "T563")]
tse <- tse[, which(colData(tse)$Response %in% c("No_DC", "DC"))]

# Remove the "s__" prefix from the "species" column in rowData
rowData(tse)$species <- gsub("^s__", "", rowData(tse)$species)

# Extract species level
tse_species <- mergeFeaturesByRank(tse, rank = "species", onRankOnly = TRUE)

# Replace all spaces and hyphens with underscores in row names of tse_species
rownames(tse_species) <- gsub("[ -]", "_", rownames(tse_species))

# Apply genus filtering
tse_species <- filter_species(tse_species)

# --------- MANUAL RELATIVE ABUNDANCE TRANSFORMATION ---------
counts_matrix <- assay(tse_species, "counts")
rel_abundance_matrix <- sweep(counts_matrix, 2, colSums(counts_matrix), FUN = "/")
assay(tse_species, "relabundance") <- rel_abundance_matrix
tse_rel <- tse_species

# Extract all Prevotella species from tse_rel
all_prevotella <- tse_rel[grepl("Prevotella", rowData(tse_rel)$species, ignore.case = TRUE), ]

# Convert colData to a data frame
metadata_df <- as.data.frame(colData(tse_rel))

# Transpose assay data
tse_assay <- t(assay(tse_rel, "relabundance"))
tse_assay_df <- as.data.frame(tse_assay)

tse_assay_filtered <- tse_assay_df[, apply(tse_assay_df, 2, function(x) any(x > 0, na.rm = TRUE))]

merged_data_all_prevotella <- merge(metadata_df, tse_assay_filtered, by.x = "row.names", by.y = "row.names", all.x = TRUE)

all_prevotella_columns <- grep("^Prevotella", colnames(merged_data_all_prevotella), value = TRUE, ignore.case = TRUE)

heatmap_data_all_prevotella <- merged_data_all_prevotella %>%
  select(Response_Patient_day, all_of(all_prevotella_columns)) %>%
  pivot_longer(cols = -Response_Patient_day, names_to = "Taxon", values_to = "Abundance")

first_5_x_min <- 1
first_5_x_max <- 5
last_4_x_min <- (length(unique(heatmap_data_all_prevotella$Response_Patient_day)) - 3)
last_4_x_max <- length(unique(heatmap_data_all_prevotella$Response_Patient_day))

heatmap_data_all_prevotella$Response_Patient_day <- gsub("_1$", "", heatmap_data_all_prevotella$Response_Patient_day)

heatmap_plot <- ggplot(heatmap_data_all_prevotella, aes(x = Response_Patient_day, y = reorder(Taxon, -Abundance), fill = Abundance)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("lightyellow", "orange", "red"), 
                       values = scales::rescale(c(0, 0.001, 0.0025, 0.005, max(heatmap_data_all_prevotella$Abundance))),
                       limits = c(0, max(heatmap_data_all_prevotella$Abundance)),  
                       na.value = "white") +  
  labs(title = "Prevotella species", 
       x = "", 
       y = "", 
       fill = "Rel. Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  
        axis.text.y = element_text(size = 10), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  annotate("rect", xmin = first_5_x_min - 0.5, xmax = first_5_x_max + 0.5, 
           ymin = 0.5, ymax = length(unique(heatmap_data_all_prevotella$Taxon)) + 0.5, 
           color = "black", fill = NA, linewidth = 1.2) +
  
  annotate("rect", xmin = last_4_x_min - 0.5, xmax = last_4_x_max + 0.5, 
           ymin = 0.5, ymax = length(unique(heatmap_data_all_prevotella$Taxon)) + 0.5, 
           color = "black", fill = NA, linewidth = 1.2)

print(heatmap_plot)





###### Heatmap with signficance based on wilcoxon test

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
library(MicrobiomeStat)
library(tibble)
library(ANCOMBC)
library(foreach)
library(rngtools)
library(ggplot2)
library(tidyr)

filter_species <- function(tse_subset) {
  genus_counts <- rowSums(assay(tse_subset) > 0)  # Count non-zero occurrences
  tse_subset_filtered <- tse_subset[genus_counts >= 0, ]  # Keep genera in ≥1 samples
  return(tse_subset_filtered)
}

# Load data
tse <- tse_Metaphlan

tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Response %in% c("No_DC", "DC"))]

# Remove the "s__" prefix from the "species" column in rowData
rowData(tse)$species <- gsub("^s__", "", rowData(tse)$species)

# Extract species level
tse_species <- mergeFeaturesByRank(tse, rank = "species", onRankOnly = TRUE)

# Replace all spaces and hyphens with underscores in row names of tse_species
rownames(tse_species) <- gsub("[ -]", "_", rownames(tse_species))

# Apply genus filtering
tse_species <- filter_species(tse_species)

# --------- MANUAL RELATIVE ABUNDANCE TRANSFORMATION ---------
counts_matrix <- assay(tse_species, "counts")
rel_abundance_matrix <- sweep(counts_matrix, 2, colSums(counts_matrix), FUN = "/")
assay(tse_species, "relabundance") <- rel_abundance_matrix
tse_rel <- tse_species

# Extract all Prevotella species from tse_rel
all_prevotella <- tse_rel[grepl("Prevotella", rowData(tse_rel)$species, ignore.case = TRUE), ]

# Convert colData to a data frame
metadata_df <- as.data.frame(colData(tse_rel))

# Transpose assay data
tse_assay <- t(assay(tse_rel, "relabundance"))
tse_assay_df <- as.data.frame(tse_assay)

tse_assay_filtered <- tse_assay_df[, apply(tse_assay_df, 2, function(x) any(x > 0, na.rm = TRUE))]

merged_data_all_prevotella <- merge(metadata_df, tse_assay_filtered, by.x = "row.names", by.y = "row.names", all.x = TRUE)

all_prevotella_columns <- grep("^Prevotella", colnames(merged_data_all_prevotella), value = TRUE, ignore.case = TRUE)

heatmap_data_all_prevotella <- merged_data_all_prevotella %>%
  select(Response_Patient_day, all_of(all_prevotella_columns)) %>%
  pivot_longer(cols = -Response_Patient_day, names_to = "Taxon", values_to = "Abundance")

first_5_x_min <- 1
first_5_x_max <- 5
last_4_x_min <- (length(unique(heatmap_data_all_prevotella$Response_Patient_day)) - 3)
last_4_x_max <- length(unique(heatmap_data_all_prevotella$Response_Patient_day))

heatmap_data_all_prevotella$Response_Patient_day <- gsub("_1$", "", heatmap_data_all_prevotella$Response_Patient_day)

library(ggplot2)
library(dplyr)
library(ggsignif)

# Run Wilcoxon rank sum test for each species and store raw p-values
p_values <- sapply(all_prevotella_columns, function(species) {
  result <- wilcox.test(as.formula(paste(species, "~ Response")), 
                        data = merged_data_all_prevotella,
                        subset = merged_data_all_prevotella$Response %in% c("DC", "No_DC"),
                        exact = FALSE)
  return(result$p.value)
})

# Assign significance labels based on raw p-values
significance_labels <- ifelse(p_values < 0.001, "***", 
                              ifelse(p_values < 0.01, "**", 
                                     ifelse(p_values < 0.05, "*", "")))

# Create a data frame for annotation
significance_df <- data.frame(
  Taxon = names(p_values),
  Significance = significance_labels
)

# Merge significance labels with heatmap data
heatmap_data_all_prevotella <- left_join(heatmap_data_all_prevotella, significance_df, by = "Taxon")

# Plot heatmap with significance annotations
heatmap_plot <- ggplot(heatmap_data_all_prevotella, aes(x = Response_Patient_day, y = reorder(Taxon, -Abundance), fill = Abundance)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = c("lightyellow", "orange", "red"), 
                       values = scales::rescale(c(0, 0.001, 0.0025, 0.005, max(heatmap_data_all_prevotella$Abundance))),
                       limits = c(0, max(heatmap_data_all_prevotella$Abundance)),  
                       na.value = "white") +  
  labs(title = "Prevotella species", 
       x = "", 
       y = "", 
       fill = "Rel. Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),  
        axis.text.y = element_text(size = 10), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  
  annotate("rect", xmin = first_5_x_min - 0.5, xmax = first_5_x_max + 0.5, 
           ymin = 0.5, ymax = length(unique(heatmap_data_all_prevotella$Taxon)) + 0.5, 
           color = "black", fill = NA, linewidth = 1.2) +
  
  annotate("rect", xmin = last_4_x_min - 0.5, xmax = last_4_x_max + 0.5, 
           ymin = 0.5, ymax = length(unique(heatmap_data_all_prevotella$Taxon)) + 0.5, 
           color = "black", fill = NA, linewidth = 1.2) +
  
  # Add significance labels on the heatmap
  geom_text(data = heatmap_data_all_prevotella, 
            aes(x = max(Response_Patient_day), y = Taxon, label = Significance), 
            color = "black", size = 5, hjust = -0.5) 

print(heatmap_plot)






########################################################################
###### Functional analysis - Humann3

# Dot plot - KEGG
library(TreeSummarizedExperiment)
library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)

# Gene family data
# EC
#gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_EC_unstratified.txt", header = T)
# KO
#gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_KO_unstratified.txt", header = T)
# GO
gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_GO_unstratified.txt", header = T)


# --- Step 1: Set row names from first column and clean gene_family_data ---
rownames(gene_family_data) <- gene_family_data[, 1]
gene_family_data <- gene_family_data[, -1]

# Extract sample IDs from column names (e.g., S00TL.0001)
extract_sample_id <- function(x) str_extract(x, "S00TL\\.\\d{4}")
colnames(gene_family_data) <- sapply(colnames(gene_family_data), extract_sample_id)

# --- Step 2: Load metadata ---
metadata <- read.table(
  "path/human_urine/Batch 3/puhti/metadata.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)

# Extract and format sample ID from metadata row names
metadata$Sample_ID <- str_extract(rownames(metadata), "S00TL-\\d{4}")
metadata$Sample_ID <- gsub("-", ".", metadata$Sample_ID)  # match colnames format

# --- Step 3: Match and filter ---
common_ids <- intersect(colnames(gene_family_data), metadata$Sample_ID)
gene_family_data <- gene_family_data[, common_ids, drop = FALSE]
metadata <- metadata[metadata$Sample_ID %in% common_ids, ]
metadata <- metadata[match(common_ids, metadata$Sample_ID), ]
rownames(metadata) <- metadata$Sample_ID
stopifnot(all(colnames(gene_family_data) == rownames(metadata)))

# --- Step 4: Create TSE object ---
tse <- TreeSummarizedExperiment(
  assays = list(counts = gene_family_data),
  colData = metadata
)

# Keep only Day 1 and samples with known response
tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Response %in% c("DC", "No_DC"))]

# --- Step 5: Wilcoxon test per feature ---
counts <- assay(tse)
meta <- colData(tse)
group <- factor(meta$Response, levels = c("DC", "No_DC"))

wilcox_results <- apply(counts, 1, function(feature_counts) {
  test <- wilcox.test(feature_counts ~ group)
  return(c(statistic = test$statistic, pval = test$p.value))
})

wilcox_df <- as.data.frame(t(wilcox_results))
wilcox_df$Feature <- rownames(wilcox_df)
wilcox_df$padj <- p.adjust(wilcox_df$pval, method = "fdr")

# --- Step 6: Add GO categories ---
wilcox_df$GO_cat <- case_when(
  str_detect(wilcox_df$Feature, "\\[BP\\]") ~ "BP",
  str_detect(wilcox_df$Feature, "\\[MF\\]") ~ "MF",
  str_detect(wilcox_df$Feature, "\\[CC\\]") ~ "CC",
  TRUE ~ "Other"
)

# --- Step 7: Filter significant results and calculate GeneRatio ---
significant_results <- wilcox_df %>% filter(pval < 0.05)

abundance_matrix_subset <- counts[significant_results$Feature, , drop = FALSE]
significant_results$GeneRatio <- rowSums(abundance_matrix_subset > 0) / ncol(abundance_matrix_subset)

# --- Step 8: Prepare data for plotting ---
prep_plot_data <- function(df, go_type) {
  df %>%
    filter(GO_cat == go_type, pval < 0.05) %>%
    mutate(Feature = gsub("\\.", " ", Feature)) %>%
    arrange(pval) %>%
    mutate(Feature = factor(Feature, levels = rev(Feature)))
}

# --- Step 9: Plot function ---
plot_go <- function(plot_data, go_type) {
  if (nrow(plot_data) == 0) return(NULL)
  
  ggplot(plot_data, aes(x = GeneRatio, y = Feature)) +
    geom_point(aes(color = pval), size = 4) +
    scale_color_gradient(low = "#2E8B57", high = "skyblue") +
    theme_bw(base_size = 14) +
    labs(
      title = paste("Functional Enrichment -", go_type),
      x = "Gene Ratio",
      y = NULL,
      color = "p-value"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      legend.title = element_text(face = "bold")
    )
}

# --- Step 10: Generate plots per GO type ---
plot_bp <- plot_go(prep_plot_data(significant_results, "BP"), "BP")
plot_mf <- plot_go(prep_plot_data(significant_results, "MF"), "MF")
plot_cc <- plot_go(prep_plot_data(significant_results, "CC"), "CC")

# --- Step 11: Show results ---
print(plot_bp)
print(plot_mf)
print(plot_cc)

# Or combine vertically
(plot_bp / plot_mf / plot_cc)


write.xlsx(significant_results[, c("Feature", "pval")],
           file = "path/human_urine/RESULTS/T563_shotgun_metagenomics/wilcoxon_kegg_results.xlsx",
           row.names = FALSE)




### Maaslin2


# Gene family data
# EC
#gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_EC_unstratified.txt", header = T)
# KO
#gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_KO_unstratified.txt", header = T)
# GO
gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_GO_unstratified.txt", header = T)
# MetaCyc (ERROR!)
#gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_MetaCyc_unstratified.txt", header = T)


# Set row names from first column and remove that column
rownames(gene_family_data) <- gene_family_data[, 1]
gene_family_data <- gene_family_data[, -1]

# Extract sample IDs from column names (e.g., S00TL.0001)
extract_sample_id <- function(x) str_extract(x, "S00TL\\.\\d{4}")
colnames(gene_family_data) <- sapply(colnames(gene_family_data), extract_sample_id)

# --- Load metadata ---
metadata <- read.table(
  "path/human_urine/Batch 3/puhti/metadata.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)

# Extract and format sample ID from metadata row names
metadata$Sample_ID <- str_extract(rownames(metadata), "S00TL-\\d{4}")
metadata$Sample_ID <- gsub("-", ".", metadata$Sample_ID)  # Convert to match format in gene_family_data

# --- Match and filter based on common sample IDs ---
common_ids <- intersect(colnames(gene_family_data), metadata$Sample_ID)

# Subset and reorder
gene_family_data <- gene_family_data[, common_ids, drop = FALSE]
metadata <- metadata[metadata$Sample_ID %in% common_ids, ]
metadata <- metadata[match(common_ids, metadata$Sample_ID), ]

# --- Set metadata rownames to Sample_ID (not Response_Patient_day!) ---
rownames(metadata) <- metadata$Sample_ID

# --- Final check: column names of gene_family_data must match rownames of metadata ---
stopifnot(all(colnames(gene_family_data) == rownames(metadata)))

# --- Create TreeSummarizedExperiment object ---
tse <- TreeSummarizedExperiment(
  assays = list(counts = gene_family_data),
  colData = metadata
)

tse

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$'Response' %in% c("DC", "No_DC"))]

# Step 4: Check the TSE object
print(tse)


# Maaslin_DC
library(Maaslin2)

# Higher in DC, run this line before Maaslin2 run
colData(tse)$Response <- relevel(as.factor(colData(tse)$Response), ref = "No_DC")


maaslin2_out <- Maaslin2(input_data = as.data.frame(t(assay(tse))),
                         input_metadata = as.data.frame(colData(tse)),
                         output = "EC_1",
                         transform = "AST",
                         fixed_effects = c("Response", "DC", "No_DC"),
                         # you can also fit MLM by specifying random effects
                         # random_effects = c(...),
                         normalization = "TSS",
                         standardize = FALSE,
                         # filtering was previously performed
                         min_prevalence = 0)

library(dplyr)
library(knitr)

maaslin2_out$results %>%
  filter(qval < 0.5) %>%
  knitr::kable()


#heatmap

# Filter results with p-value < 0.01
significant_results <- maaslin2_out$results %>%
  filter(pval < 0.01)

# Set row names of significant_results to values in column 1
rownames(significant_results) <- significant_results[, 1]

# Remove the first column as it's now redundant as row names
significant_results <- significant_results[, -1]

# Extract the row names from the significant_results dataframe
original_rownames <- rownames(significant_results)

# Remove the leading "X" only if it is present at the start
cleaned_rownames <- sub("^X", "", original_rownames)
# Replace any non-alphanumeric character with a period
cleaned_rownames <- gsub("[^[:alnum:]]", ".", cleaned_rownames)


# Set the modified row names back to the significant_results dataframe
rownames(significant_results) <- cleaned_rownames

# Print the new row names to verify
print(rownames(significant_results))




# Extract row names from the TreeSummarizedExperiment object
rownames_tse <- rownames(tse)

# Remove the leading "X" if it is present at the start
cleaned_rownames_tse <- sub("^X", "", rownames_tse)
# Replace any non-alphanumeric character with a period
cleaned_rownames_tse <- gsub("[^[:alnum:]]", ".", cleaned_rownames_tse)

# Set the cleaned row names back to the tse dataframe
rownames(tse) <- cleaned_rownames_tse

# Print the new row names to verify
print(rownames(tse))






# Extract abundance values for significant features
abundance_matrix <- assay(tse)[rownames(significant_results), ]

# Print the new row names to verify
print(rownames(abundance_matrix))


### Heatmap
# Extract metadata2 for grouping
metadata_tse <- colData(tse)
group_column <- "Response"

# Filter metadata to include only data for which "Day" = "1"
metadata_filter1 <- metadata_tse[metadata_tse$Day == "1", ]

# Filter metadata2 for No_DC and DC groups
metadata_filter2 <- metadata_filter1[metadata_filter1[, group_column] %in% c("DC", "No_DC"), ]

# Extract row names from metadata_filter2
selected_patients <- rownames(metadata_filter2)

# Subset abundance_matrix to keep only columns present in selected_patients
abundance_matrix_subset <- abundance_matrix[, selected_patients, drop = FALSE]

# -----------------------------
# ✅ Custom name_mapping
# -----------------------------
name_mapping <- c(
  "S00TL.0001" = "No_DC_30104_1",
  "S00TL.0003" = "DC_30105_1",
  "S00TL.0006" = "DC_30207_1",
  "S00TL.0009" = "No_DC_30209_1",
  "S00TL.0012" = "No_DC_30110_1",
  "S00TL.0015" = "DC_30111_1",
  "S00TL.0018" = "DC_30112_1",
  "S00TL.0021" = "DC_30210_1",
  "S00TL.0026" = "No_DC_30212_1"
)

# Map old column names to new names using the name_mapping
old_col_names <- colnames(abundance_matrix_subset)
new_col_names <- name_mapping[old_col_names]

# Fallback to original names if not mapped
new_col_names[is.na(new_col_names)] <- old_col_names[is.na(new_col_names)]

# Apply new column names
colnames(abundance_matrix_subset) <- new_col_names

# Print to verify
print(colnames(abundance_matrix_subset))

# -----------------------------
# ✅ Plot heatmaphttp://127.0.0.1:11347/graphics/plot_zoom_png?width=1239&height=760
# -----------------------------
library(pheatmap)

pheatmap(as.matrix(abundance_matrix_subset),
         width = 5,
         height = 5,
         fontsize_row = 14,
         fontsize_col = 18)




## Different plot
#library(ggplot2)
#library(dplyr)

# Step 1: Prepare plot data
plot_data <- significant_results %>%
  mutate(
    Feature = rownames(significant_results),
    Count = rowSums(abundance_matrix_subset[rownames(significant_results), ] > 0),  # how many samples express it
    GeneRatio = Count / ncol(abundance_matrix_subset),  # ratio of samples
    Category = "KO"  # you can later split by "EC", "GO", etc.
  )

# Optional: clean up long feature names (e.g. KO terms)
plot_data$Feature <- gsub("_", " ", plot_data$Feature)

# Order features by effect size or p-value for better plotting
plot_data <- plot_data %>%
  arrange(pval) %>%
  mutate(Feature = factor(Feature, levels = rev(Feature)))

# Step 2: Dot Plot
ggplot(plot_data, aes(x = GeneRatio, y = Feature)) +
  geom_point(aes(size = Count, color = pval)) +
  scale_color_gradient(low = "red", high = "skyblue") +
  theme_bw(base_size = 14) +
  labs(
    title = "Functional Enrichment of KEGG Features",
    x = "Gene Ratio",
    y = NULL,
    color = "p-value",
    size = "Sample Count"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )

## Colors based on Effect size (coef)
ggplot(plot_data, aes(x = GeneRatio, y = Feature)) +
  geom_point(aes(size = Count, color = coef)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_bw(base_size = 14) +
  labs(
    title = "Functional Enrichment by Response (DC vs No_DC)",
    x = "Gene Ratio",
    y = NULL,
    color = "Effect Size (coef)",
    size = "Sample Count"
  ) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_text(face = "bold")
  )





################ barplots 
library(stringr)
# EC
#gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_EC_unstratified.txt", header = T)
# KO
#gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_KO_unstratified.txt", header = T)
# GO
gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_GO_unstratified.txt", header = T)
# MetaCyc (ERROR!)
#gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_MetaCyc_unstratified.txt", header = T)


# Set row names from first column and remove that column
rownames(gene_family_data) <- gene_family_data[, 1]
gene_family_data <- gene_family_data[, -1]

# Extract sample IDs from column names (e.g., S00TL.0001)
extract_sample_id <- function(x) str_extract(x, "S00TL\\.\\d{4}")
colnames(gene_family_data) <- sapply(colnames(gene_family_data), extract_sample_id)

# --- Load metadata ---
metadata <- read.table(
  "path/human_urine/Batch 3/puhti/metadata_GG2.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)

# Extract and format sample ID from metadata row names
metadata$Sample_ID <- str_extract(rownames(metadata), "S00TL-\\d{4}")
metadata$Sample_ID <- gsub("-", ".", metadata$Sample_ID)  # Convert to match format in gene_family_data

# --- Match and filter based on common sample IDs ---
common_ids <- intersect(colnames(gene_family_data), metadata$Sample_ID)

# Subset and reorder
gene_family_data <- gene_family_data[, common_ids, drop = FALSE]
metadata <- metadata[metadata$Sample_ID %in% common_ids, ]
metadata <- metadata[match(common_ids, metadata$Sample_ID), ]

# --- Set metadata rownames to Sample_ID (not Response_Patient_day!) ---
rownames(metadata) <- metadata$Sample_ID

# --- Final check: column names of gene_family_data must match rownames of metadata ---
stopifnot(all(colnames(gene_family_data) == rownames(metadata)))

# --- Create TreeSummarizedExperiment object ---
tse <- TreeSummarizedExperiment(
  assays = list(counts = gene_family_data),
  colData = metadata
)

tse

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$'Response' %in% c("DC", "No_DC"))]


# Step 4: Check the TSE object
print(tse)

# Convert colData to a data frame
metadata_df <- as.data.frame(colData(tse))

# Transpose assay data
tse_assay <- t(assay(tse))

# Merge metadata with transposed assay data
merged_data <- merge(metadata_df, as.data.frame(tse_assay), by.x = "row.names", by.y = "row.names", all.x = TRUE)

# Rename patient 20103 to 20103(PR)
#merged_data <- merged_data %>%
#  mutate(Patient = ifelse(Patient == 20103, "20103(PR)", as.character(Patient)))

# Ensure correct factor order
merged_data$Response <- factor(merged_data$Response, levels = c("DC", "No_DC"))

# Define colors for plot
fill_colors <- c("DC" = "skyblue", "No_DC" = "#9DC183")
point_colors <- c("DC" = "darkblue", "No_DC" = "darkgreen")

# ---- STEP 1: Identify glycan-related KO columns ----
#glycan_cols <- grep("glycan|mucin|sialidase|fucosidase|hexosaminidase|galactosidase|glycosidase|neuraminidase", 
#                   colnames(merged_data), ignore.case = TRUE, value = TRUE)
glycan_cols <- grep("galactosidase|fucosidase|sialidase|glycoside hydrolase|glycosyltransferase|polysaccharide lyase|carbohydrate esterase", 
                    colnames(merged_data), ignore.case = TRUE, value = TRUE)
#glycan_cols <- grep("peptidoglycan biosynthetic process", 
#                    colnames(merged_data), ignore.case = TRUE, value = TRUE)

# Optional: Filter out columns with all NA or no variation
glycan_cols <- glycan_cols[sapply(glycan_cols, function(col) {
  values <- merged_data[[col]]
  !all(is.na(values)) && length(unique(na.omit(values))) > 1
})]

library(stringr)
library(patchwork)

# ---- STEP 2: Assign GO category to each glycan column ----
go_map <- data.frame(
  col = glycan_cols,
  GO_cat = case_when(
    str_detect(glycan_cols, "\\[BP]") ~ "BP",
    str_detect(glycan_cols, "\\[MF]") ~ "MF",
    str_detect(glycan_cols, "\\[CC]") ~ "CC",
    TRUE ~ "Other"
  ),
  stringsAsFactors = FALSE
)

# ---- STEP 3: Make plots per GO category ----
make_plots_for_go <- function(go_type) {
  cols_for_go <- go_map$col[go_map$GO_cat == go_type]
  
  if (length(cols_for_go) == 0) return(NULL)  # skip if no matches
  
  glycan_plots <- lapply(cols_for_go, function(col) {
    plot_title <- paste0(col, "\n")
    make_gene_plot(col, plot_title)
  })
  
  wrap_plots(glycan_plots, ncol = 3) + 
    plot_annotation(title = paste("GO category:", go_type))
}

# ---- STEP 4: Generate separate plots ----
plot_bp <- make_plots_for_go("BP")
plot_mf <- make_plots_for_go("MF")
plot_cc <- make_plots_for_go("CC")

# Show them one by one
print(plot_bp)
print(plot_mf)
print(plot_cc)

# Or combine vertically
#(plot_bp / plot_mf / plot_cc)



# Load libraries
library(ggplot2)
library(ggsignif)
library(patchwork)
library(grid)        # For grid.newpage()
library(grDevices)   # For dev.new() if needed

# Number of plots per page
plots_per_page <- 3

# Divide glycan plots into chunks
plot_chunks <- split(glycan_plots, ceiling(seq_along(glycan_plots) / plots_per_page))

# Loop over each chunk and open a new window for it
for (i in seq_along(plot_chunks)) {
  # Optional: open a new graphics window (not needed in RStudio but useful in some IDEs)
  if (!interactive()) dev.new(width = 14, height = 5)  # only opens if not in RStudio
  
  # Start a new plotting page
  grid.newpage()
  
  # Combine and print the chunk
  chunk_plot <- wrap_plots(plot_chunks[[i]], ncol = 3)
  print(chunk_plot)
  
  # Optional: pause between pages
  if (interactive() && i < length(plot_chunks)) {
    readline(prompt = "Press [Enter] to see the next page of plots...")
  }
}





### TEST! Heatmap for functional analysis of 30111

# Dot plot - KEGG
library(TreeSummarizedExperiment)
library(stringr)
library(TreeSummarizedExperiment)
library(stringr)

# Gene family data
# EC
#gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_EC_unstratified.txt", header = T)
# KO
#gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_KO_unstratified.txt", header = T)
# GO
gene_family_data <- read.delim("path/human_urine/Batch 3/Humann3/RenormRename_genefamilies_Uniref90_GO_unstratified.txt", header = T)


# Set row names from first column and remove that column
rownames(gene_family_data) <- gene_family_data[, 1]
gene_family_data <- gene_family_data[, -1]

# Extract sample IDs from column names (e.g., S00TL.0001)
extract_sample_id <- function(x) str_extract(x, "S00TL\\.\\d{4}")
colnames(gene_family_data) <- sapply(colnames(gene_family_data), extract_sample_id)

# --- Load metadata ---
metadata <- read.table(
  "path/human_urine/Batch 3/puhti/metadata_GG2.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)

# Extract and format sample ID from metadata row names
metadata$Sample_ID <- str_extract(rownames(metadata), "S00TL-\\d{4}")
metadata$Sample_ID <- gsub("-", ".", metadata$Sample_ID)  # Convert to match format in gene_family_data

# --- Match and filter based on common sample IDs ---
common_ids <- intersect(colnames(gene_family_data), metadata$Sample_ID)

# Subset and reorder
gene_family_data <- gene_family_data[, common_ids, drop = FALSE]
metadata <- metadata[metadata$Sample_ID %in% common_ids, ]
metadata <- metadata[match(common_ids, metadata$Sample_ID), ]

# --- Set metadata rownames to Sample_ID (not Response_Patient_day!) ---
rownames(metadata) <- metadata$Sample_ID

# --- Final check: column names of gene_family_data must match rownames of metadata ---
stopifnot(all(colnames(gene_family_data) == rownames(metadata)))

# --- Create TreeSummarizedExperiment object ---
tse <- TreeSummarizedExperiment(
  assays = list(counts = gene_family_data),
  colData = metadata
)

tse

# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Patient == "30111")]




# Extract count data and metadata
counts <- assay(tse)
meta <- colData(tse)

# Subset to day 1 and patient 30111
tse <- tse[, which(colData(tse)$Day == "1" & colData(tse)$Patient == "30111")]

# Extract count data and metadata
counts <- assay(tse)
meta <- colData(tse)

# Convert counts into a data frame for plotting
expr_df <- as.data.frame(counts)
expr_df$Feature <- rownames(expr_df)

# If multiple samples exist for this patient, average them
expr_df$MeanExpr <- rowMeans(expr_df[, -ncol(expr_df), drop = FALSE])

# Add a GO category column based on rownames/Feature
expr_df$GO_cat <- case_when(
  str_detect(expr_df$Feature, "\\[BP\\]") ~ "BP",
  str_detect(expr_df$Feature, "\\[MF\\]") ~ "MF",
  str_detect(expr_df$Feature, "\\[CC\\]") ~ "CC",
  TRUE ~ "Other"
)

# Function to plot top 30 for a category
plot_top_go <- function(df, go_type) {
  top_expr <- df %>%
    filter(GO_cat == go_type) %>%
    arrange(desc(MeanExpr)) %>%
    head(30)
  
  ggplot(top_expr, aes(x = MeanExpr, y = reorder(Feature, MeanExpr))) +
    geom_point(color = "#2E8B57", size = 3) +
    theme_bw(base_size = 14) +
    labs(
      title = paste("Top 30", go_type, "GO terms in Patient 30111 (Day 1)"),
      x = "Mean Expression",
      y = "GO Term"
    ) +
    theme(
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
}

# Generate three plots
p_bp <- plot_top_go(expr_df, "BP")
p_mf <- plot_top_go(expr_df, "MF")
p_cc <- plot_top_go(expr_df, "CC")

# If you want them together
library(patchwork)
p_bp / p_mf / p_cc


#write.xlsx(significant_results[, c("Feature", "pval")],
#           file = "path/human_urine/RESULTS/T563_shotgun_metagenomics/wilcoxon_kegg_results.xlsx",
#           row.names = FALSE)





## Patient overview with fixed effect (legenda on bottom)

library(readxl)
library(dplyr)
library(tidyr)

# Read data
fixed_effect_df <- read_excel(
  "path/human_urine/Fixed_effects.xlsx"
)


fixed_effect_df <- read_excel(
  "path/human_urine/Fixed_effects.xlsx"
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
  
  # Age bins
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

# ---------------- Alphabetical row order (top → bottom) ----------------
variable_order <- c("Response", "Age", "Hospital", "Primary Tumor", "Histological subtype", "Stage at diagnosis")
df_long$Variable <- factor(df_long$Variable, levels = variable_order)

# ---------------- Sample order by increasing Microbiome present ----------------
microbiome_order <- df_long %>%
  filter(Variable == "Response") %>%
  mutate(MicrobiomeNum = ifelse(Value == "DC", 1, 0)) %>%
  arrange(MicrobiomeNum) %>%
  pull(Patients)

df_long$Patients <- factor(df_long$Patients, levels = microbiome_order)



# ---------------- Plot ----------------
ggplot() +
  
  # ---------------- Response ----------------
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
  
  # ---------------- Age ----------------
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
  
  # ---------------- Hospital ----------------
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
  
  # ---------------- Primary Tumor ----------------
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
  
  # ---------------- Histological subtype ----------------
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
  
  # ---------------- Stage ----------------
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
  
  # ---------------- Theme ----------------
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
