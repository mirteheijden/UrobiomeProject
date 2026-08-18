# load
library(TreeSummarizedExperiment)
library(taxa)
library(microbiome) 
library(dplyr)
library(readr)
library(stringr)
library(tidyverse)
library(mia)
library(qiime2R)
library(miaViz)
library(phyloseq)
library(ggplot2)
library(patchwork)
library(ecodist)
library(data.table)
library(vegan)
library(ggforce)
library(RColorBrewer)
library(readr)
library(tibble)
library(openxlsx)
library(ggrepel)
library(scater)
library(reshape2)
library(scales) 
library(viridis) 
library(forcats)  
library(ggsignif)

# Import data
kraken2_data <- read_csv("path/Kraken2_data_bacteria_only.csv")

# Assign Taxonomic Ranks using lvl_type
lvl_type_to_rank <- c(
  R = "root", D = "domain", D1 = "domain", D2 = "domain",
  K = "kingdom",
  P = "phylum", P1 = "phylum",
  C = "class", C1 = "class", C2 = "class",
  O = "order", O1 = "order", O2 = "order",
  F = "family", F1 = "family",
  G = "genus", G1 = "genus",
  S = "species",
  U = "unclassified"
)

kraken2_data <- kraken2_data %>%
  mutate(rank = lvl_type_to_rank[lvl_type])

# Prepare the Counts Matrix
counts <- kraken2_data[, grep("_all$", colnames(kraken2_data))]
counts <- counts[, !grepl("^tot_", colnames(counts))]  # Exclude 'tot_' columns
#colnames(counts) <- gsub("_all$", "", colnames(counts))  # Remove '_all' suffix

# Make unique rownames using taxonomic prefix
rownames(counts) <- paste0(
  kraken2_data$rank, "__",
  str_replace_all(kraken2_data$name, " ", "_"),
  "_", seq_len(nrow(kraken2_data))
)

# Prepare Taxonomy (Row Data)
taxonomy <- dplyr::select(kraken2_data, taxid, name, rank) %>%
  dplyr::rename(name = name)

# Include Metadata (Column Data)
metadata <- read_tsv("path/config/metadata.tsv")
colnames_ordered <- colnames(counts)
metadata <- metadata %>% arrange(match(kraken_id, colnames_ordered))

# Create the TreeSummarizedExperiment Object
tse <- TreeSummarizedExperiment(
  assays = list(counts = as.matrix(counts)),
  rowData = taxonomy,
  colData = DataFrame(metadata)
)


### Adjust the lay-out for rowData(tse), so that it is accepted by tse codes
# Extract the `rowData(tse)` as a dataframe
taxonomy_df <- as.data.frame(rowData(tse))

# Reshape the dataframe so that ranks are columns and taxid is the rowname
reshaped_taxonomy <- taxonomy_df %>%
  # Spread the data by rank
  pivot_wider(names_from = rank, values_from = name, values_fn = list(name = first)) %>%
  # Set taxid as rownames
  column_to_rownames("taxid")

# Assign the reshaped dataframe back to rowData(tse)
rowData(tse) <- DataFrame(reshaped_taxonomy)

# Check the result
head(rowData(tse))

### Venn Diagrams
#install.packages("ggVennDiagram")
library(ggVennDiagram)

### Adjust the lay-out for rowData(tse), so that it is accepted by tse codes
# Extract the `rowData(tse)` as a dataframe
taxonomy_df <- as.data.frame(rowData(tse))

# Reshape the dataframe so that ranks are columns and taxid is the rowname
reshaped_taxonomy <- taxonomy_df %>%
  # Spread the data by rank
  pivot_wider(names_from = rank, values_from = name, values_fn = list(name = first)) %>%
  # Set taxid as rownames
  column_to_rownames("taxid")

# Assign the reshaped dataframe back to rowData(tse)
rowData(tse) <- DataFrame(reshaped_taxonomy)

# Check the result
head(rowData(tse))


## Differential abundance
colData(tse)$Tumor_Group <- ifelse(
  colData(tse)$Tumor_Type %in% c("Ovarian"),
  "Ovarian",
  "Other"
)


#tse <- tse[, rownames(subset(colData(tse), Time == "day_1"))]
tse_rel <- transformAssay(tse, assay.type = "counts", method = "relabundance")
colSums(assay(tse_rel, "relabundance"))

# correct order for your available ranks
correct_order <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus")

# reorder taxonomy columns
rowData(tse_rel) <- rowData(tse_rel)[, correct_order]

tse_genus <- mergeFeaturesByRank(
  tse_rel,
  rank = "genus",
  onRankOnly = TRUE
)


# Function to filter genera present in at least 4 samples
filter_genus <- function(tse_genus) {
  genus_counts <- rowSums(assay(tse_genus) > 0)  
  tse_genus_filtered <- tse_genus[genus_counts >= 4, ]  # Keep genera in ≥4 samples
  return(tse_genus_filtered)
}

# Apply filtering
Other_rel <- filter_genus(tse_genus[, rownames(subset(colData(tse_genus), Tumor_Group == "Other"))])
Ovarian_rel <- filter_genus(tse_genus[, rownames(subset(colData(tse_genus), Tumor_Group == "Ovarian"))])

# Extract unique genus names from each group
genus_other <- unique(rowData(Other_rel)$genus)
genus_ovarian <- unique(rowData(Ovarian_rel)$genus)

# Find overlapping genera
overlapping_genus <- intersect(genus_other, genus_ovarian)

# Create a list for ggVennDiagram
venn_list <- list(Other = genus_other, Ovarian = genus_ovarian)

set.seed(20231214)

# plot Venn diagram
ggVennDiagram(venn_list, label_alpha = 0, label_size = 8) +  
  scale_fill_gradient(low = "#D4A017", high = "#C2A3E0") +  
  theme_void() +  
  theme(
    plot.margin = margin(10, 10, 10, 10), 
    legend.position = "none",  
    text = element_text(size = 16, face = "bold"),  
    plot.title = element_text(hjust = 0.50)
  ) +
  labs(title = "Baseline") +
  coord_flip() 

# Print genera
# Find exclusive genera
only_other_genus <- setdiff(genus_other, overlapping_genus)   
only_ovarian_genus <- setdiff(genus_ovarian, overlapping_genus) 

# Sort the results alphabetically
overlapping_genus <- sort(overlapping_genus)
only_other_genus <- sort(only_other_genus)
only_ovarian_genus <- sort(only_ovarian_genus)

# View results
View(as.data.frame(overlapping_genus))  
View(as.data.frame(only_other_genus))        
View(as.data.frame(only_ovarian_genus))      

#Save in Excel

library(openxlsx)

# Convert to vectors
other <- as.character(only_other_genus)
ovarian <- as.character(only_ovarian_genus)
overlapping <- as.character(overlapping_genus)

# Pad with NA to make equal length
max_len <- max(length(other), length(ovarian), length(overlapping))

pad <- function(x) c(x, rep(NA, max_len - length(x)))

df <- data.frame(
  Other = pad(other),
  Ovarian = pad(ovarian),
  Overlapping = pad(overlapping)
)

# Save to Excel
write.xlsx(df, file = "Z:/PhD CGTG/Experiments/Metagenomics/human_urine/RESULTS/genus_comparison.xlsx", overwrite = TRUE)

# Number of samples used for Venn Diagram
table(colData(tse_genus)$Tumor_Group)



## Taxonomic barplots

# Phylum
taxonomy_df <- as.data.frame(rowData(tse))

reshaped <- taxonomy_df %>%
  pivot_wider(
    names_from = rank,
    values_from = name,
    values_fn = list(name = dplyr::first)
  ) %>%
  column_to_rownames("taxid")

rowData(tse) <- DataFrame(reshaped)

# Correct rank column order
correct_order <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus")
rowData(tse) <- rowData(tse)[, intersect(correct_order, colnames(rowData(tse)))]

#tse <- tse[, rownames(subset(colData(tse), Time == "day_1"))]
tse <- tse[, rownames(subset(colData(tse), Tumor_Type == "Ovarian"))]
#tse <- tse[, rownames(subset(colData(tse), Tumor_Type != "Ovarian"))]
tse_spec <- agglomerateByRank(tse, rank = "phylum", onRankOnly = TRUE)
tse_rel  <- transformAssay(tse_spec, assay.type = "counts", method = "relabundance")

# Define "Other" phyla (<0.1% mean abundance) 
phylum_means <- rowMeans(assay(tse_rel, "relabundance"))
rowData(tse_rel)$phylum <- ifelse(phylum_means >= 0.001,
                                  rowData(tse_rel)$phylum,
                                  "Other")

# Merge small phyla
tse_combined <- mergeFeaturesByRank(
  tse_rel,
  rank = "phylum",
  onRankOnly = TRUE
)

# Drop zero rows
keep <- rowSums(assay(tse_combined, "relabundance")) > 0
tse_filtered <- tse_combined[keep, ]

# Update column names to Patient_day_inj
cd <- colData(tse_combined)
rownames(cd) <- cd$Patient_day_inj
colData(tse_combined) <- cd

# Remove the "Other" row for plotting
tse_filtered <- tse_combined[rowData(tse_combined)$phylum != "Other", ]

# plot
# Build an ordering key and factor it by Tumor_Type -> Patient
cd <- as.data.frame(colData(tse_combined))
stopifnot(all(c("Tumor_Type", "Patient_day") %in% colnames(cd)))

# Key for sorting and for readable labels 
cd$Tumor_Patient <- paste0(cd$Tumor_Type, "_", cd$Patient_day)

# Compute the multi-criteria order
ord <- order(cd$Tumor_Type, cd$Patient_day, na.last = TRUE)

# Make Tumor_Patient an ordered factor according to ord
cd$Tumor_Patient <- factor(cd$Tumor_Patient, levels = cd$Tumor_Patient[ord], ordered = TRUE)

# Write back to colData
colData(tse_combined) <- S4Vectors::DataFrame(cd)

# Change rownames of the colData to the Patient ID
col_data_tse <- colData(tse_combined)
rownames(col_data_tse) <- col_data_tse$Tumor_Patient
colData(tse_combined) <- col_data_tse

# Filter out the "Other"  from rowData
tse_combined <- tse_combined[rowData(tse_combined)$phylum != "Other", ]

# plotAbundance order columns
plotAbundance(
  tse_combined,
  rank = "phylum",
  assay.type = "relabundance",
  order.row.by = "abund",
  order.col.by = "Tumor_Type",  
  add_x_text = TRUE
) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    axis.text.x = element_text(angle = 55, hjust = 1, size = 12),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14)
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Phylum in ovarian cancer", x = "", y = "Relative Abundance (%)") +
  scale_fill_viridis_d(option = "C") +
  labs(fill = "Phylum", colour = NULL) +
  guides(
    colour = "none",
    size = "none",
    shape = "none",
    linetype = "none",
    alpha = "none"
  )



# genus barplot
taxonomy_df <- as.data.frame(rowData(tse))

reshaped <- taxonomy_df %>%
  pivot_wider(
    names_from = rank,
    values_from = name,
    values_fn = list(name = dplyr::first)
  ) %>%
  column_to_rownames("taxid")

rowData(tse) <- DataFrame(reshaped)

# Correct rank column order 
correct_order <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus")
rowData(tse) <- rowData(tse)[, intersect(correct_order, colnames(rowData(tse)))]

#tse <- tse[, rownames(subset(colData(tse), Time == "day_1"))]
tse <- tse[, rownames(subset(colData(tse), Tumor_Type == "Ovarian"))]
#tse <- tse[, rownames(subset(colData(tse), Tumor_Type != "Ovarian"))]
tse_spec <- agglomerateByRank(tse, rank = "genus", onRankOnly = TRUE)
tse_rel  <- transformAssay(tse_spec, assay.type = "counts", method = "relabundance")

# Define "Other" phyla (<1% mean abundance) 
genus_means <- rowMeans(assay(tse_rel, "relabundance"))
rowData(tse_rel)$genus <- ifelse(genus_means >= 0.01,
                                 rowData(tse_rel)$genus,
                                 "Other")

# Merge small phyla
tse_combined <- mergeFeaturesByRank(
  tse_rel,
  rank = "genus",
  onRankOnly = TRUE
)

# Drop zero rows
keep <- rowSums(assay(tse_combined, "relabundance")) > 0
tse_filtered <- tse_combined[keep, ]

# Update column names to Patient_day_inj
cd <- colData(tse_combined)
rownames(cd) <- cd$Patient_day_inj
colData(tse_combined) <- cd

# Remove the "Other" row for plotting
tse_filtered <- tse_combined[rowData(tse_combined)$genus != "Other", ]

# plot
# Build an ordering key and factor it by Tumor_Type -> Patient
cd <- as.data.frame(colData(tse_combined))
stopifnot(all(c("Tumor_Type", "Patient_day") %in% colnames(cd)))

# Key for sorting and for readable labels if you want
cd$Tumor_Patient <- paste0(cd$Tumor_Type, "_", cd$Patient_day)

# Compute the multi-criteria order
ord <- order(cd$Tumor_Type, cd$Patient_day, na.last = TRUE)

# Make Tumor_Patient an ordered factor according to ord
cd$Tumor_Patient <- factor(cd$Tumor_Patient, levels = cd$Tumor_Patient[ord], ordered = TRUE)

# Write back to colData
colData(tse_combined) <- S4Vectors::DataFrame(cd)

# Change rownames of the colData to the Patient ID
col_data_tse <- colData(tse_combined)
rownames(col_data_tse) <- col_data_tse$Tumor_Patient
colData(tse_combined) <- col_data_tse

# Filter out the "Other"  from rowData
tse_combined <- tse_combined[rowData(tse_combined)$genus != "Other", ]

# plotAbundance order columns
plotAbundance(
  tse_combined,
  rank = "genus",
  assay.type = "relabundance",
  order.row.by = "abund",
  order.col.by = "Tumor_Type",  
  add_x_text = TRUE
) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    axis.text.x = element_text(angle = 55, hjust = 1, size = 12),
    axis.title.y = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 14)
  ) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Genus in ovarian cancer", x = "", y = "Relative Abundance (%)") +
  scale_fill_viridis_d(option = "C") +
  labs(fill = "Genus", colour = NULL) +
  guides(
    colour = "none",
    size = "none",
    shape = "none",
    linetype = "none",
    alpha = "none"
  )





## Phylum specific boxplot

# create relabundance
tse_rel <- transformAssay(
  tse,
  assay.type = "counts",
  method = "relabundance"
)

# extract bacteria
target_phylum <- "Bacteroidota"

phylum_row <- which(
  rowData(tse_rel)$rank == "phylum" &
    rowData(tse_rel)$name == target_phylum
)

rowData(tse_rel)[phylum_row, ]

# Filter
#tse_rel <- tse_rel[, rownames(subset(colData(tse_rel), Time == "day_1"))]

# plot
plot_df <- as.data.frame(colData(tse_rel))

plot_df$Abundance <- as.numeric(
  assay(tse_rel, "relabundance")[phylum_row, ]
)

# Collapse all non-ovarian cancers into "Other"
plot_df$Tumor_Type <- ifelse(
  plot_df$Tumor_Type == "Ovarian",
  "Ovarian",
  "Other"
)

plot_df$Tumor_Type <- factor(
  plot_df$Tumor_Type,
  levels = c("Ovarian", "Other")
)

ggplot(
  plot_df,
  aes(x = Tumor_Type,
      y = Abundance,
      fill = Tumor_Type)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    aes(color = Tumor_Type),
    width = 0.15,
    size = 2
  ) +
  #geom_text_repel(
  #  aes(label = Patient_day),
  #  size = 3,
  #  show.legend = FALSE
  #) +
  geom_signif(
    comparisons = list(c("Ovarian", "Other")),
    test = "wilcox.test",
    map_signif_level = function(p)
      paste0("p = ", signif(p, 3))
  ) +
  scale_fill_manual(
    values = c(
      "Ovarian" = "#C2A3E0",
      "Other" = "#D4A017"
    )
  ) +
  scale_color_manual(
    values = c(
      "Ovarian" = "#C2A3E0",
      "Other" = "#D4A017"
    )
  ) +
  scale_y_log10() +
  labs(
    title = paste(target_phylum),
    y = "Relative abundance",
    x = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0.5
    )
  )



## genus specific boxplot

# create relabundance
tse_rel <- transformAssay(
  tse,
  assay.type = "counts",
  method = "relabundance"
)

# extract bacteria
target_genus <- "Anaerococcus"

genus_row <- which(
  rowData(tse_rel)$rank == "genus" &
    rowData(tse_rel)$name == target_genus
)

rowData(tse_rel)[genus_row, ]

# Filter
#tse_rel <- tse_rel[, rownames(subset(colData(tse_rel), Time == "day_1"))]

# plot
plot_df <- as.data.frame(colData(tse_rel))

plot_df$Abundance <- as.numeric(
  assay(tse_rel, "relabundance")[genus_row, ]
)

# Collapse all non-ovarian cancers into "Other"
plot_df$Tumor_Type <- ifelse(
  plot_df$Tumor_Type == "Ovarian",
  "Ovarian",
  "Other"
)

plot_df$Tumor_Type <- factor(
  plot_df$Tumor_Type,
  levels = c("Ovarian", "Other")
)

ggplot(
  plot_df,
  aes(x = Tumor_Type,
      y = Abundance,
      fill = Tumor_Type)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    aes(color = Tumor_Type),
    width = 0.15,
    size = 2
  ) +
  #geom_text_repel(
  #  aes(label = Patient_day),
  #  size = 3,
  #  show.legend = FALSE
  #) +
  geom_signif(
    comparisons = list(c("Ovarian", "Other")),
    test = "wilcox.test",
    map_signif_level = function(p)
      paste0("p = ", signif(p, 3))
  ) +
  scale_fill_manual(
    values = c(
      "Ovarian" = "#C2A3E0",
      "Other" = "#D4A017"
    )
  ) +
  scale_color_manual(
    values = c(
      "Ovarian" = "#C2A3E0",
      "Other" = "#D4A017"
    )
  ) +
  scale_y_log10() +
  labs(
    title = target_genus,
    y = "Relative abundance",
    x = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0.5
    )
  )



########## Ovarian patient cohort
#BiocManager::install("ANCOMBC")

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
library(readr)

## Phyloseq
phyloseq <- readRDS("path/RESULTS_ampliseq/phyloseq/dada2_phyloseq.rds")
head(otu_table(phyloseq))
head(tax_table(phyloseq))

# Read the metadata TSV file
metadata_path <- "path/config/metadata.tsv"
sample_metadata_df <- read_tsv(metadata_path)

# Create the TreeSummarizedExperiment object from the phyloseq object
tse_phylo <- mia::convertFromPhyloseq(phyloseq)

# Add the sample metadata to the colData of the TSE object
colData(tse_phylo) <- DataFrame(sample_metadata_df)

# Check the TSE object
print(tse_phylo)
colData(tse_phylo)

# change rownames of the tse, to our sample_ID
col_data_tse_phylo <- colData(tse_phylo)
rownames(col_data_tse_phylo) <- col_data_tse_phylo$Patient_day
colData(tse_phylo) <- col_data_tse_phylo
col_data_tse_phylo <- as.data.frame(colData(tse_phylo))

colData(tse_phylo)

# change rownames of the tse, to our sample_ID
col_data_phylo <- colData(tse_phylo)
rownames(col_data_phylo) <- col_data_phylo$Patient_day
colData(tse_phylo) <- col_data_phylo
col_data_phylo <- as.data.frame(colData(tse_phylo))



##### Response - Genus
# After creating tse_phylo, replace row names with Genus names
row_data <- rowData(tse_phylo)

# Create unique genus names (add number suffix for duplicates)
genus_names <- as.character(row_data$Genus)
genus_names[is.na(genus_names) | genus_names == ""] <- "Unknown"
genus_names <- make.unique(genus_names, sep = "_")

# Replace row names in the TSE
rownames(tse_phylo) <- genus_names

tse <- tse_phylo


# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Trial == "T563")]
tse <- tse[, which(colData(tse)$Response %in% c("No_DC", "DC"))]

# Transform count assay to relative abundances
tse_rel <- transformAssay(tse,
                          assay.type = "counts",
                          method = "relabundance")

# Convert to phyloseq object
pseq <- makePhyloseqFromTreeSummarizedExperiment(tse_rel)
pseq_Genus <- phyloseq::tax_glom(pseq, taxrank = "Genus")

# Perform ANCOMBC analysis, comparing DC vs No_DC
out_gg2_Genus = ancombc(
  data = pseq_Genus, 
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
res <- out_gg2_Genus$res

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

# Extract the row names and Genus column from rowData(tse_rel)
row_data <- rowData(tse_rel)
row_names_to_Genus <- data.frame(
  row_name = rownames(row_data),
  Genus = row_data$Genus
)

# Create a lookup table for easy matching
lookup_table <- setNames(row_names_to_Genus$Genus, row_names_to_Genus$row_name)

# Replace values in the 'taxon' column of significant_results
significant_results <- significant_results %>%
  mutate(taxon = ifelse(taxon %in% names(lookup_table), 
                        lookup_table[taxon], 
                        taxon))

# Add Response column based on log_fold_change direction
significant_results$Response <- ifelse(significant_results$log_fold_change > 0, "DC", "No_DC")

# Assign colors based on Response only
significant_results$highlight <- ifelse(significant_results$Response == "DC", "blue", "orange")
significant_results$highlight <- factor(significant_results$highlight, 
                                        levels = c("blue", "orange"))

# Create the bar plot with custom colors
bar_plot <- ggplot(significant_results, aes(x = reorder(taxon, log_fold_change), y = log_fold_change, fill = highlight)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  scale_fill_manual(name = NULL,
                    values = c("blue" = "blue", "orange" = "orange"),
                    labels = c("blue" = "DC", 
                               "orange" = "No_DC")) +
  coord_flip() +
  labs(title = "Differential abundant genus 
            (16S rRNAseq)", x = "Taxon", y = "Log Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14, face = "italic"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),  
        plot.title = element_text(size = 16, face = "bold", hjust = 0),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(title = ""))

# Print the plot
print(bar_plot)











##### Response - Species
# After creating tse_phylo, replace row names with Species names
row_data <- rowData(tse_phylo)

# Create unique Species names (add number suffix for duplicates)
Species_names <- as.character(row_data$Species)
Species_names[is.na(Species_names) | Species_names == ""] <- "Unknown"
Species_names <- make.unique(Species_names, sep = "_")

# Replace row names in the TSE
rownames(tse_phylo) <- Species_names

tse <- tse_phylo


# Extract day 1 samples and samples with a known disease outcome
tse <- tse[, which(colData(tse)$Day == "1")]
tse <- tse[, which(colData(tse)$Trial == "T563")]
tse <- tse[, which(colData(tse)$Response %in% c("No_DC", "DC"))]

# Transform count assay to relative abundances
tse_rel <- transformAssay(tse,
                          assay.type = "counts",
                          method = "relabundance")

# Check relative abundance
assay_data <- assay(tse_rel, "relabundance")

col_sums <- colSums(assay_data)

summary(col_sums)

assayNames(tse_rel)

# Convert to phyloseq object
pseq <- makePhyloseqFromTreeSummarizedExperiment(tse_rel)
pseq_Species <- phyloseq::tax_glom(pseq, taxrank = "Species")

# Perform ANCOMBC analysis, comparing DC vs No_DC
out_gg2_Species = ancombc(
  data = pseq_Species, 
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
res <- out_gg2_Species$res

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

# Extract the row names and Species column from rowData(tse_rel)
row_data <- rowData(tse_rel)
row_names_to_Species <- data.frame(
  row_name = rownames(row_data),
  Species = row_data$Species
)

# Create a lookup table for easy matching
lookup_table <- setNames(row_names_to_Species$Species, row_names_to_Species$row_name)

# Replace values in the 'taxon' column of significant_results
significant_results <- significant_results %>%
  mutate(taxon = ifelse(taxon %in% names(lookup_table), 
                        lookup_table[taxon], 
                        taxon))

# Add Response column based on log_fold_change direction
significant_results$Response <- ifelse(significant_results$log_fold_change > 0, "DC", "No_DC")

# Add Response column based on log_fold_change direction
significant_results$Response <- ifelse(significant_results$log_fold_change > 0, "DC", "No_DC")

# Assign colors based on Response only
significant_results$highlight <- ifelse(significant_results$Response == "DC", "lightblue", "#FFC300")

# Ensure 'highlight' is a factor
significant_results$highlight <- factor(significant_results$highlight, 
                                        levels = c("lightblue", "#FFC300"))

# Create the bar plot with custom colors
bar_plot <- ggplot(significant_results, aes(x = reorder(taxon, log_fold_change), y = log_fold_change, fill = highlight)) +
  geom_bar(stat = "identity", show.legend = TRUE) +
  scale_fill_manual(name = NULL,
                    values = c("lightblue" = "lightblue", "#FFC300" = "#FFC300"),
                    labels = c("lightblue" = "DC", 
                               "#FFC300" = "No_DC")) +
  coord_flip() +
  labs(title = "Differential abundant species 
            (16S rRNAseq)",, x = "Taxon", y = "Log Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14, face = "italic"),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),  
        plot.title = element_text(size = 20, face = "bold", hjust = 0),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(title = ""))

# Print the plot
print(bar_plot)






