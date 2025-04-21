## set the working directory 

setwd("C:/Users/ham063/OneDrive - University of Bergen/Desktop/Papers/teeth bacteria/Github")

# upload libraries 

library(phyloseq)
library(microbiome)
library(ggplot2)
library(ggpubr)
library(tidyverse)               
library(readxl)
library(ampvis2)
library(dplyr)
library(table1)
library(openxlsx)
library(writexl)

# making the phyloseq object from otu table, taxa table and metadata table. 

OTU <- read_xlsx("Oral_data_bg_f.xlsx", sheet = "OTU")
TAXA <- read_xlsx("Oral_data_bg_f.xlsx", sheet = "Taxa")

# Convert the OTU tibble to a matrix
otu_matrix <- as.matrix(OTU[, -1])  # Exclude the first column with OTU names

# Convert the matrix to a phyloseq OTU table

otu_physeq <- otu_table(otu_matrix, taxa_are_rows = TRUE)
# Convert taxonomy table to phyloseq object

taxonomy_physeq <- tax_table(as.matrix(TAXA))

# Create a phyloseq object by merging OTU and taxonomy tables

BORALIS <- merge_phyloseq(otu_physeq, taxonomy_physeq)
metadata <-  read.csv("oral_metadata_2.csv", sep =',', header=TRUE, row.names=1)
sample_data(BORALIS) <- metadata


#### Figure 2 ####

## extracting the  alpha diversity indeces from the phyloseq file ### 


# alpha index 

#calculating alpha diversity indeces 
hmp.div <- alpha(BORALIS, index ="all") 
help("alpha")

# get the metadata out as separate object
hmp.meta <- meta(BORALIS)

# Add the rownames as a new colum for easy integration later.
hmp.meta$sam_name <- rownames(hmp.meta)


# Add the rownames to diversity table

hmp.div$sam_name <- rownames(hmp.div)


# merge these two data frames into one

div.df <- merge(hmp.div,hmp.meta, by = "sam_name")

names(div.df)
table(div.df$Rekke)

T0.df <-  div.df %>% filter (Rekke =="T0")
T1.df <-  div.df %>% filter (Rekke =="T1")

my_comparisons = list(c("T0", "T1"))


shannon.rekke = div.df %>% ggplot(aes(x = Rekke, y = diversity_shannon)) + 
  geom_boxplot(width = 0.5, aes(fill = Rekke), alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.7) + 
  labs(x = NULL, y = "Shannon diversity index") + 
  stat_compare_means(comparisons = my_comparisons) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_brewer(name = NULL, palette = "Dark2",label = c("Baseline (T0)","Follow-up (T1)"))



observed.rekke  =  div.df %>% ggplot(aes(x = Rekke, y = observed)) + 
  geom_boxplot(width = 0.5, aes(fill = Rekke), alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.7) + 
  labs(x = NULL, y = "Observed number of taxa") + 
  stat_compare_means(comparisons = my_comparisons) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_brewer(name = NULL, palette = "Dark2",label = c("Baseline (T0)","Follow-up (T1)"))


## PCoA 

library(ampvis2)
mymetadata <- read_excel("oral_metadata_2.xlsx")

# Step 1: Extract the otutable from phyloseq object
otutable <- data.frame(
  OTU = rownames(phyloseq::otu_table(BORALIS, taxa_are_rows = FALSE)@.Data),
  phyloseq::otu_table(BORALIS, taxa_are_rows = FALSE)@.Data,
  phyloseq::tax_table(BORALIS)@.Data,
  check.names = FALSE
)

# Step 2: Assign ASVs as unique identifiers for OTUs
row.names(otutable) <- paste("ASV", 1:nrow(otutable), sep = "")
otutable$OTU <- paste("ASV", 1:nrow(otutable), sep = "")

# Step 3: Remove the 'MGS' column if it exists
if("MGS" %in% colnames(otutable)) {
  otutable <- otutable[, !colnames(otutable) %in% "MGS"]
}

# Step 4: Load metadata
mymetadata <- read_excel("oral_metadata_2.xlsx")

# Step 5: Verify and align sample names
# Extract sample columns (exclude OTU and taxonomy columns)
taxonomy_columns <- c("OTU", "phylum", "class", "order", "family", "genus", "species", "Species")
existing_taxonomy_columns <- taxonomy_columns[taxonomy_columns %in% colnames(otutable)]
otutable_samples <- setdiff(colnames(otutable), existing_taxonomy_columns)

# Find common samples
common_samples <- intersect(otutable_samples, mymetadata$ID)

# Filter otutable and metadata to include only common samples
otutable <- otutable[, c(existing_taxonomy_columns, common_samples)]
mymetadata <- mymetadata[mymetadata$ID %in% common_samples, ]

# Proceed with ampvis
ampvis <- amp_load(otutable = otutable, metadata = mymetadata)


ord.plot = amp_ordinate(ampvis,
                        type = "pcoa",
                        filter_species = 0,
                        transform = "none",
                        sample_color_by = "Rekke",
                        distmeasure = "bray") 


########### beta diversity

library(vegan)
library(pairwiseAdonis)

bray_curtis_dist <- phyloseq::distance(BORALIS, method = "bray")

# Perform pairwise PERMANOVA for FOT Yes vs NO
pairwise_perm_anova <- pairwise.adonis(bray_curtis_dist, as.factor(sample_data(BORALIS)$Rekke))


ord.plot.rekke <- ord.plot + 
  annotate("text", x = 0.01, y = 0.36, 
           label = expression(paste("Pairwise permanova  P = 0.001 / ", R^2, " = 0.163"))) + 
  scale_color_manual(
    name = "Time point",
    values = c("T0" = "#1B9E77", "T1" = "#D95F02"), # Custom green and orange colors
    labels = c("Baseline (T0)", "Follow-up (T1)")
  ) +
  aes(color = factor(Rekke))


# Arrange A and B in a single row
top_row <- ggarrange(shannon.rekke, observed.rekke, 
                     ncol = 2, labels = c("A", "B"))

# Arrange top row and C in a vertical layout
Figure.2 <- ggarrange(top_row, ord.plot.rekke, 
                        ncol = 1, heights = c(1, 1), 
                        labels = c("", "C"))
# Print the final plot
Figure.2


ggsave(filename = "Figure 2.jpeg", plot = Figure.2,
       height = 9, width = 8, units = 'in', dpi = 300)

#### Figure 3 ####

Genus.Rekke  = amp_heatmap(
  ampvis, 
  group_by = "Rekke", 
  tax_aggregate = "Genus", 
  tax_show = 30, 
  plot_colorscale = "sqrt", 
  plot_values = TRUE, 
  min_abundance = 0.1, 
  measure = "mean"
) +
  theme(
    axis.text.x = element_text(size = 11, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_x_discrete(labels = c("T0" = "Baseline (T0)", "T1" = "Follow-up (T1)"))


Species.Rekke  = amp_heatmap(
  ampvis, 
  group_by = "Rekke", 
  tax_aggregate = "Species", 
  tax_show = 30, 
  plot_colorscale = "sqrt", 
  plot_values = TRUE, 
  min_abundance = 0.1, 
  measure = "mean"
) +
  theme(
    axis.text.x = element_text(size = 11, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_x_discrete(labels = c("T0" = "Baseline (T0)", "T1" = "Follow-up (T1)"))

Figure.3  <- ggarrange(Genus.Rekke,Species.Rekke, ncol = 2, nrow = 1, 
                            widths =  c(1, 1), labels = c("A","B"),vjust = 1)

ggsave(filename = "Figure 3.jpeg", plot = Figure.3,
       height = 9, width = 12, units = 'in', dpi = 300)

### Figure 4 ####

# Load required libraries
library(Maaslin2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(readr)
library(RColorBrewer)

# Load metadata and OTU table from your existing Phyloseq object
metadata <- read.csv("oral_metadata_2.csv", sep = ",", header = TRUE, row.names = 1)
otu_table <- as.data.frame(otu_table(BORALIS))
taxonomy <- as.data.frame(tax_table(BORALIS))
taxonomy$feature <- rownames(taxonomy)

# Optional filtering based on abundance
relative_abundance <- rowSums(otu_table) / sum(otu_table) * 100
otu_table_filtered <- otu_table[relative_abundance >= 0.01, ]

# Ensure OTU and metadata alignment
otu_table_filtered <- otu_table_filtered[, colnames(otu_table_filtered) %in% rownames(metadata)]

# Function to run MaAsLin2 and extract significant results only for the main variable
run_maaslin2 <- function(fixed_effect) {
  output_dir <- paste0("maaslin2_", fixed_effect)
  Maaslin2(
    input_data = otu_table_filtered,
    input_metadata = metadata,
    fixed_effects = c(fixed_effect, "Rekke", "sex", "BMI", "age"),
    random_effects = c("ID_1"),
    normalization = "TSS",
    transform = "LOG",
    output = output_dir
  )
  sig <- read_tsv(file.path(output_dir, "significant_results.tsv"))
  sig <- sig %>% filter(metadata == fixed_effect & qval <= 0.05)
  sig$Dataset <- fixed_effect
  return(sig)
}

# Run for R5, R11, R19 and combine
sig_R5 <- run_maaslin2("R5")
sig_R11 <- run_maaslin2("R11")
sig_R19 <- run_maaslin2("R19")

combined_data <- bind_rows(sig_R5, sig_R11, sig_R19)
combined_data <- merge(combined_data, taxonomy, by = "feature")

# Set factor order for plotting
combined_data$Dataset <- factor(combined_data$Dataset, levels = c("R5", "R11", "R19"))

# Genus label placement above bars
genus_label_positions <- combined_data %>%
  group_by(species, genus) %>%
  summarize(y_pos = max(coef) + 0.2, .groups = "drop")

# Plot
bar_dodge <- position_dodge2(width = 0.7, preserve = "single")

Figure.4 <- ggplot() +  
  geom_bar(data = combined_data, aes(x = species, y = coef, fill = Dataset),
           stat = "identity", position = bar_dodge, width = 0.7, color = "black") +  
  geom_text(data = genus_label_positions, aes(x = species, y = y_pos, label = genus),
            color = "black", size = 5.5, angle = 45, hjust = 0.3, vjust = 0) +  
  theme_minimal(base_size = 14) +
  labs(x = "Bacterial Species", y = "Coefficient", fill = NULL) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 12),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey80")) +
  scale_fill_manual(
    values = c("R5" = "#66C2A5", "R11" = "#FC8D62", "R19" = "#8DA0CB"),
    labels = c(expression(R[5]), expression(R[11]), expression(R[19]))
  )

# Save plot
ggsave("Figure 4.jpeg", plot = Figure.4, height = 8, width = 14, dpi = 300)


### Figure 5 ####

# merge these two data frames into one

names(div.df)
table(div.df$baseline_R5)

div.df.above <-  div.df %>% filter (baseline_R5 =="above")
div.df.below <-  div.df %>% filter (baseline_R5 =="equal or below")

table(div.df.above$Rekke)
table(div.df.below$Rekke)

my_comparisons = list(c("T0", "T1"))


shannon.below = div.df.below %>% ggplot(aes(x = Rekke, y = diversity_shannon)) + 
  geom_boxplot(width = 0.5, aes(fill = Rekke), alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.7) + 
  labs(x = NULL, y = "Shannon diversity index") + 
  stat_compare_means(comparisons = my_comparisons) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_brewer(name = NULL, palette = "Dark2",label = c("Baseline (T0)","Follow-up (T1)"))

# 
observed.below  =  div.df.below %>% ggplot(aes(x = Rekke, y = observed)) + 
  geom_boxplot(width = 0.5, aes(fill = Rekke), alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.7) + 
  labs(x = NULL, y = "Observed number of taxa") + 
  stat_compare_means(comparisons = my_comparisons) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_brewer(name = NULL, palette = "Dark2",label = c("Baseline (T0)","Follow-up (T1)"))


shannon.above = div.df.above %>% ggplot(aes(x = Rekke, y = diversity_shannon)) + 
  geom_boxplot(width = 0.5, aes(fill = Rekke), alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.7) + 
  labs(x = NULL, y = "Shannon diversity index") + 
  stat_compare_means(comparisons = my_comparisons) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_brewer(name = NULL, palette = "Dark2",label = c("Baseline (T0)","Follow-up (T1)"))

# 
observed.above =  div.df.above %>% ggplot(aes(x = Rekke, y = observed)) + 
  geom_boxplot(width = 0.5, aes(fill = Rekke), alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.7) + 
  labs(x = NULL, y = "Observed number of taxa") + 
  stat_compare_means(comparisons = my_comparisons) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_brewer(name = NULL, palette = "Dark2",label = c("Baseline (T0)","Follow-up (T1)"))

## calculate the beta diversity 

library(vegan)
library(pairwiseAdonis)


## subset the phyloseq 
BORALIS.below <- subset_samples(BORALIS , baseline_R5 == "equal or below")

BORALIS.above  <- subset_samples(BORALIS , baseline_R5 == "above")

### calcualte the distance matrix 
bray_curtis_dist.below <- phyloseq::distance(BORALIS.below, method = "bray")

bray_curtis_dist.above <- phyloseq::distance(BORALIS.above, method = "bray")

# Perform pairwise PERMANOVA 
pairwise_perm_anova.below <- pairwise.adonis(bray_curtis_dist.below, as.factor(sample_data(BORALIS.below)$Rekke))

pairwise_perm_anova.above <- pairwise.adonis(bray_curtis_dist.above, as.factor(sample_data(BORALIS.above)$Rekke))



#### ordination plot for baseline_R5

ampvis_below <- amp_filter_samples(ampvis, baseline_R5 %in% c("equal or below"))

ampvis_above <- amp_filter_samples(ampvis, baseline_R5 %in% c("above"))


ord.plot.below = amp_ordinate(ampvis_below,
                              type = "pcoa",
                              filter_species = 0,
                              transform = "none",
                              sample_color_by = "Rekke",
                              distmeasure = "bray") 

ord.plot.below <- ord.plot.below + 
  annotate("text", x = 0.01, y = 0.45, 
           label = expression(paste("Pairwise permanova  P = 0.001 / ", R^2, " = 0.187"))) + 
  scale_color_manual(
    name = "Time point",
    values = c("T0" = "#1B9E77", "T1" = "#D95F02"), # Custom green and orange colors
    labels = c("Baseline (T0)", "Follow-up (T1)")
  ) +
  aes(color = factor(Rekke))




ord.plot.above = amp_ordinate(ampvis_above,
                              type = "pcoa",
                              filter_species = 0,
                              transform = "none",
                              sample_color_by = "Rekke",
                              distmeasure = "bray") 

ord.plot.above <- ord.plot.above + 
  annotate("text", x = 0.01, y = 0.37, 
           label = expression(paste("Pairwise permanova  P = 0.001 / ", R^2, " = 0.162"))) +
  scale_color_manual(
    name = "Time point",
    values = c("T0" = "#1B9E77", "T1" = "#D95F02"), # Custom green and orange colors
    labels = c("Baseline (T0)", "Follow-up (T1)")
  ) +
  aes(color = factor(Rekke))

### mergeing the plots 

library(ggpubr)
library(cowplot)

# Create titles using cowplot
titles <- ggdraw() +
  draw_label(expression("Individuals below/equal to the R"[5]*" median"), 
             x = 0.25, y = 0.5, size = 16, fontface = "bold", hjust = 0.5) +
  draw_label(expression("Individuals above the R"[5]*" median"), 
             x = 0.75, y = 0.5, size = 16, fontface = "bold", hjust = 0.5)

# Adjust individual plots with white space padding and add labels
left_column <- plot_grid(
  shannon.below + theme(plot.margin = margin(7, 30, 7, 30)), 
  observed.below + theme(plot.margin = margin(7, 30, 7, 30)),
  ord.plot.below + theme(plot.margin = margin(7, 25, 7, 25)),
  ncol = 1, align = 'v',
  labels = c("A", "C", "E"), label_size = 14, label_fontface = "bold"
)

right_column <- plot_grid(
  shannon.above + theme(plot.margin = margin(7, 30, 7, 30)),
  observed.above + theme(plot.margin = margin(7, 30, 7, 30)),
  ord.plot.above + theme(plot.margin = margin(7, 25, 7, 25)),
  ncol = 1, align = 'v',
  labels = c("B", "D", "F"), label_size = 14, label_fontface = "bold"
)

# Combine left and right columns with spacing
plots_with_space <- plot_grid(
  left_column, NULL, right_column,
  ncol = 3, rel_widths = c(0.45, 0.01, 0.45) # Reduce widths for more space around the red line
)

# Combine titles and plots
combined_plot <- plot_grid(
  titles, plots_with_space,
  ncol = 1, rel_heights = c(0.08, 2) # Reduce title space
)

# Add the red dotted line
Figure.5 <- ggdraw() +
  draw_plot(combined_plot) +
  draw_line(x = c(0.5, 0.5), y = c(0, 1), color = "red", size = 1.5 , linetype = "dotted")

# Display the final plot
print(Figure.5)

# Save the adjusted figure
ggsave(filename = "Figure 5.jpeg", plot = Figure.5,
       height = 17, width = 16, units = 'in', dpi = 300)




### Figure 6 ####


ampvis_below <- amp_filter_samples(ampvis, baseline_R5 %in% c("equal or below"))

ampvis_above <- amp_filter_samples(ampvis, baseline_R5 %in% c("above"))


Genus.below  = amp_heatmap(ampvis_below, group_by = "Rekke",tax_aggregate = "Genus",
                           tax_show = 30, 
                           plot_colorscale = "sqrt", 
                           plot_values = TRUE, 
                           min_abundance = 0.1, 
                           measure = "mean"
) +
  theme(
    axis.text.x = element_text(size = 11, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_x_discrete(labels = c("T0" = "Baseline (T0)", "T1" = "Follow-up (T1)"))



Species.below  = amp_heatmap(ampvis_below, group_by = "Rekke",tax_aggregate = "Species",
                             tax_show = 30, 
                             plot_colorscale = "sqrt", 
                             plot_values = TRUE, 
                             min_abundance = 0.1, 
                             measure = "mean"
) +
  theme(
    axis.text.x = element_text(size = 11, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_x_discrete(labels = c("T0" = "Baseline (T0)", "T1" = "Follow-up (T1)"))



Genus.above  = amp_heatmap(ampvis_above, group_by = "Rekke",tax_aggregate = "Genus",
                           tax_show = 30, 
                           plot_colorscale = "sqrt", 
                           plot_values = TRUE, 
                           min_abundance = 0.1, 
                           measure = "mean"
) +
  theme(
    axis.text.x = element_text(size = 11, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_x_discrete(labels = c("T0" = "Baseline (T0)", "T1" = "Follow-up (T1)"))



Species.above  = amp_heatmap(ampvis_above, group_by = "Rekke",tax_aggregate = "Species",
                             tax_show = 30, 
                             plot_colorscale = "sqrt", 
                             plot_values = TRUE, 
                             min_abundance = 0.1, 
                             measure = "mean"
) +
  theme(
    axis.text.x = element_text(size = 11, hjust = 0.5, angle = 0),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  ) +
  scale_x_discrete(labels = c("T0" = "Baseline (T0)", "T1" = "Follow-up (T1)"))





library(ggpubr)
library(cowplot)

# Create titles using cowplot
titles <- ggdraw() +
  draw_label(expression("Individuals below/equal to the R"[5]*" median"), 
             x = 0.25, y = 0.5, size = 16, fontface = "bold", hjust = 0.5) +
  draw_label(expression("Individuals above the R"[5]*" median"), 
             x = 0.75, y = 0.5, size = 16, fontface = "bold", hjust = 0.5)

# Adjust individual plots with white space padding and add labels
left_column <- plot_grid(
  Genus.below + theme(plot.margin = margin(7, 30, 7, 30)), 
  Species.below + theme(plot.margin = margin(7, 30, 7, 30)),
  ncol = 1, align = 'v',
  labels = c("A", "C"), label_size = 14, label_fontface = "bold"
)

right_column <- plot_grid(
  Genus.above + theme(plot.margin = margin(7, 30, 7, 30)),
  Species.above + theme(plot.margin = margin(7, 30, 7, 30)),
  ncol = 1, align = 'v',
  labels = c("B", "D"), label_size = 14, label_fontface = "bold"
)

# Combine left and right columns with spacing
plots_with_space <- plot_grid(
  left_column, NULL, right_column,
  ncol = 3, rel_widths = c(0.45, 0.01, 0.45) # Reduce widths for more space around the red line
)

# Combine titles and plots
combined_plot <- plot_grid(
  titles, plots_with_space,
  ncol = 1, rel_heights = c(0.08, 2) # Reduce title space
)

# Add the red dotted line
Figure.6 <- ggdraw() +
  draw_plot(combined_plot) +
  draw_line(x = c(0.5, 0.5), y = c(0, 1), color = "red", size = 1.5 , linetype = "dotted")

# Display the final plot
print(Figure.6)

# Save the adjusted figure
ggsave(filename = "Figure 6.jpeg", plot = Figure.6,
       height = 17, width = 15, units = 'in', dpi = 300)

######## ANCOM BC analysis ####

library(ANCOMBC)

output.timepoint = ancombc2(data = BORALIS, fix_formula = "Rekke + sex + BMI + age", rand_formula = NULL,
                                   tax_level = "species",
                                   p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE, # Enable pseudo-count sensitivity test
                                   prv_cut = 0.1, lib_cut = 0, s0_perc = 0.05,
                                   group = "Rekke", struc_zero = TRUE, neg_lb = TRUE,
                                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                                   global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                                   iter_control = list(tol = 1e-2, max_iter = 5, verbose = FALSE), 
                                   em_control = list(tol = 1e-5, max_iter = 5),
                                   lme_control = lme4::lmerControl(),
                                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 5))


names(output.timepoint)

res_prim = output.timepoint$res

write.csv(res_prim,"ANCOM.BC.T1.vs.T0(ref).csv")



### Supplementary figure 1 ####

T0.df <-  div.df %>% filter (Rekke =="T0")

table(T0.df$baseline_R5)

T0.df$baseline_R5 <- factor(T0.df$baseline_R5,levels = c("equal or below","above"), 
                            labels = c("Equal/Below median R5", "Above median R5"))


table(T0.df$baseline_R5)


my_comparisons = list(c("Equal/Below median R5", "Above median R5"))


shannon.R5.T0 = T0.df %>% ggplot(aes(x = baseline_R5, y = diversity_shannon)) + 
  geom_boxplot(width = 0.5, aes(fill = baseline_R5), alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.7) + 
  labs(x = NULL, y = "Shannon diversity index") + 
  stat_compare_means(comparisons = my_comparisons) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_brewer(
    name = NULL, 
    palette = "Dark2",
    labels = c(
      expression("Equal or below R"[5]~"median"), 
      expression("Above R"[5]~"median")
    )
  )



# 
observed.R5.T0  =  T0.df %>% ggplot(aes(x = baseline_R5, y = observed)) + 
  geom_boxplot(width = 0.5, aes(fill = baseline_R5), alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.7) + 
  labs(x = NULL, y = "Observed number of taxa") + 
  stat_compare_means(comparisons = my_comparisons) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")+
  scale_fill_brewer(
    name = NULL, 
    palette = "Dark2",
    labels = c(
      expression("Equal or below R"[5]~"median"), 
      expression("Above R"[5]~"median")
    )
  )


## calculate the beta diversity 


## subset the phyloseq 
BORALIS.T0  <- subset_samples(BORALIS , Rekke == "T0")

### calculate the distance matrix 
bray_curtis_dist.T0 <- phyloseq::distance(BORALIS.T0, method = "bray")

# Perform pairwise PERMANOVA 
pairwise_perm_anova.T0.R5 <- pairwise.adonis(bray_curtis_dist.T0, as.factor(sample_data(BORALIS.T0)$baseline_R5))


## ord plot 

ampvis_T0 <- amp_filter_samples(ampvis, Rekke %in% c("T0"))


ord.plot.T0.R5 = amp_ordinate(ampvis_T0,
                              type = "pcoa",
                              filter_species = 0,
                              transform = "none",
                              sample_color_by = "baseline_R5",
                              distmeasure = "bray") 

# Reorder factor levels in the data
ord.plot.T0.R5$data$baseline_R5 <- factor(
  ord.plot.T0.R5$data$baseline_R5, 
  levels = c("equal or below", "above") # Specify the desired order
)

ord.plot.T0.R5 = ord.plot.T0.R5 + 
  annotate("text", x = 0.01, y = 0.3, 
           label = expression(paste("Pairwise permanova  P = 0.095 / ", R^2, " = 0.026"))) +
  aes(color = factor(baseline_R5)) +  # Ensure correct mapping to baseline_R5
  scale_color_manual(
    name = "",
    values = c("equal or below" = "#1B9E77",   # Green for equal/below R5 median
               "above" = "#D95F02"),          # Orange for above R5 median
    labels = c(
      expression("Equal or below R"[5]~"median"),
      expression("Above R"[5]~"median")
    )
  )


# Arrange A and B in a single row
top_row.T0 <- ggarrange(shannon.R5.T0, observed.R5.T0, 
                        ncol = 2, labels = c("A", "B"))

# Arrange top row and C in a vertical layout
Sup.Figure.1 <- ggarrange(top_row.T0, ord.plot.T0.R5, 
                           ncol = 1, heights = c(1, 1), 
                           labels = c("", "C"))

# Print the final plot
Sup.Figure.1

ggsave(filename = "Supplementary Figure 1.jpeg", plot = Sup.Figure.1,
       height = 9, width = 8, units = 'in', dpi = 300)



### Supplementary figure 3 ####


T1.df <-  div.df %>% filter (Rekke =="T1")

table(T1.df$baseline_R5)

T1.df$baseline_R5 <- factor(T1.df$baseline_R5,levels = c("equal or below","above"), 
                            labels = c("Equal/Below median R5", "Above median R5"))


table(T1.df$baseline_R5)


my_comparisons = list(c("Equal/Below median R5", "Above median R5"))


shannon.R5.T1 = T1.df %>% ggplot(aes(x = baseline_R5, y = diversity_shannon)) + 
  geom_boxplot(width = 0.5, aes(fill = baseline_R5), alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.7) + 
  labs(x = NULL, y = "Shannon diversity index") + 
  stat_compare_means(comparisons = my_comparisons) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") +
  scale_fill_brewer(
    name = NULL, 
    palette = "Dark2",
    labels = c(
      expression("Equal or below R"[5]~"median"), 
      expression("Above R"[5]~"median")
    )
  )



# 
observed.R5.T1  =  T1.df %>% ggplot(aes(x = baseline_R5, y = observed)) + 
  geom_boxplot(width = 0.5, aes(fill = baseline_R5), alpha = 0.7) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 0.7) + 
  labs(x = NULL, y = "Observed number of taxa") + 
  stat_compare_means(comparisons = my_comparisons) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom")+
  scale_fill_brewer(
    name = NULL, 
    palette = "Dark2",
    labels = c(
      expression("Equal or below R"[5]~"median"), 
      expression("Above R"[5]~"median")
    )
  )


## calculate the beta diversity 


## subset the phyloseq 
BORALIS.T1  <- subset_samples(BORALIS , Rekke == "T1")


### calcualte the distance matrix 
bray_curtis_dist.T1 <- phyloseq::distance(BORALIS.T1, method = "bray")

# Perform pairwise PERMANOVA 
pairwise_perm_anova.T1.R5 <- pairwise.adonis(bray_curtis_dist.T1, as.factor(sample_data(BORALIS.T1)$baseline_R5))


## ord plot 

ampvis_T1 <- amp_filter_samples(ampvis, Rekke %in% c("T1"))


ord.plot.T1.R5 = amp_ordinate(ampvis_T1,
                              type = "pcoa",
                              filter_species = 0,
                              transform = "none",
                              sample_color_by = "baseline_R5",
                              distmeasure = "bray") 

# Reorder factor levels in the data
ord.plot.T1.R5$data$baseline_R5 <- factor(
  ord.plot.T1.R5$data$baseline_R5, 
  levels = c("equal or below", "above") # Specify the desired order
)



ord.plot.T1.R5 = ord.plot.T1.R5 + 
  annotate("text", x = 0.01, y = 0.3, 
           label = expression(paste("Pairwise permanova  P = 0.033 / ", R^2, " = 0.033"))) +
  aes(color = factor(baseline_R5)) +  # Ensure correct mapping to baseline_R5
  scale_color_manual(
    name = "",
    values = c("equal or below" = "#1B9E77",   # Green for equal/below R5 median
               "above" = "#D95F02"),          # Orange for above R5 median
    labels = c(
      expression("Equal or below R"[5]~"median"),
      expression("Above R"[5]~"median")
    )
  )




# Arrange A and B in a single row
top_row.T1 <- ggarrange(shannon.R5.T1, observed.R5.T1, 
                        ncol = 2, labels = c("A", "B"))

# Arrange top row and C in a vertical layout
Sup.Fig.3 <- ggarrange(top_row.T1, ord.plot.T1.R5, 
                           ncol = 1, heights = c(1, 1), 
                           labels = c("", "C"))

# Print the final plot
Sup.Fig.3


ggsave(filename = "Supplementary Figure 3.jpeg", plot = Sup.Fig.3,
       height = 9, width = 8, units = 'in', dpi = 300)














