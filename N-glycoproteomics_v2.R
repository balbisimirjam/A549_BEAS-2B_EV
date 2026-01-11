library(tidyverse)
library(dplyr)
library(tidyr)

library(VIM)
library(lawstat)
library(rstatix)
library(ggpattern)
library(gplots)
library(ggplot2)

setwd("D:/EV_GP_241213")

source(file="functions/two_group_test.R")
source(file="functions/glycan_monomer_extractor_from10.R")

# read in GlycReSoft results files ####
setwd("D:/BM/MCP_revision/260110_GP_reeva_214GPsfromprot_top34glycansfromglycomics/glycopeptides_csv")
file_list <- list.files()

#
for (i in 1:length(file_list)) {
  print(i)
  temp_fl <- read.csv(file_list[i])
  temp_fl$Sample <- rep(gsub(".csv", "", file_list[i]), nrow(temp_fl))
  if (i==1) {
    raw_data <- temp_fl
    } else {
    raw_data <- rbind(raw_data, temp_fl)
    }
  }

# ... generate/extract categorical variables ####
# add sample name and grouping columns
raw_data$group <- ifelse(grepl("^A", raw_data$Sample), "A",
                         ifelse(grepl("^B", raw_data$Sample), "B", "M"))

# peptide, glycan, and glycopeptide columns are also generated and added
tmp_peptide_name <- rep(0,nrow(raw_data))
for (i in 1:nrow(raw_data)) {
  print(i)
  tmp_peptide_name[i] <- strsplit(gsub("\\(N-Glycosylation\\)|\\(Oxidation\\)|\\(Carbamidomethyl\\)", "*", raw_data$glycopeptide[i]), "\\{")[[1]][1]
  }
raw_data$Peptide <- tmp_peptide_name

tmp_glycan_name <- rep(0,nrow(raw_data))
for (i in 1:nrow(raw_data)) {
  print(i)
  tmp_1 <- strsplit(gsub("\\(N-Glycosylation\\)|\\(Oxidation\\)|\\(Carbamidomethyl\\)", "*", raw_data$glycopeptide[i]), "\\{")[[1]][2]
  tmp_glycan_name[i] <- gsub("[[:punct:]]","",tmp_1)
  }
raw_data$Glycan <- tmp_glycan_name

raw_data$Glycopeptide <- paste(raw_data$Peptide, raw_data$Glycan)

# extract glycan composition based on glycan name
tmp_glycan_monomers <- matrix(rep(0,4*nrow(raw_data)), ncol=4)
for (i in 1:nrow(raw_data)) {
  print(i)
  tmp_glycan_monomers[i,] <- glycan_monomer_extractor(raw_data$glycopeptide[i])
  }
raw_data$Hex <- tmp_glycan_monomers[,1]
raw_data$HexNAc <- tmp_glycan_monomers[,2]
raw_data$Neu5Ac <- tmp_glycan_monomers[,3]
raw_data$Fuc <- tmp_glycan_monomers[,4]

# calculate glycan type based on monomers
raw_data$glycan_type <- "Misc"
raw_data$glycan_type[raw_data$HexNAc > 3 & raw_data$Hex > 2] <- "Complex"
raw_data$glycan_type[raw_data$HexNAc == 3 & raw_data$Hex > 3] <- "Hybrid"
raw_data$glycan_type[raw_data$HexNAc == 2 & raw_data$Hex > 3] <- "Oligomannose"

#Additional filtering based on ms_score and peptide length
raw_data <- raw_data %>%
  mutate(peptide_length = nchar(gsub("\\*", "", Peptide)))
raw_data%>%
  filter(ms1_score > 3)%>%
  filter(ms2_score > 5)%>%
  filter(peptide_length > 5) -> raw_data_score_filtered

# filter data: throw out NAs
raw_data_score_filtered%>%
  filter(total_signal > 1)%>%
  select(glycopeptide, Peptide, Glycan, glycan_type, protein_name, Sample, group, Hex, HexNAc, Neu5Ac, Fuc, peptide_start, total_signal) -> filtered_data #apex_time-ot kivettem

# Remove duplicates by summing total_signal for reoccurring glycopeptide-Sample pairs
filtered_data_no_d <- filtered_data %>%
  group_by(glycopeptide, Sample) %>%
  summarize(
    total_signal = sum(total_signal, na.rm = TRUE), 
    Peptide = first(Peptide), 
    Glycan = first(Glycan),
    glycan_type = first(glycan_type),
    protein_name = first(protein_name),
    Hex = first(Hex),
    HexNAc = first(HexNAc),
    Neu5Ac = first(Neu5Ac),
    Fuc = first(Fuc),
    peptide_start = first(peptide_start),
    group = first(group), 
    .groups = "drop" 
  )

###
filtered_data_no_d <- filtered_data_no_d

library(stringr)
filtered_data_no_d$glycosylation_site <- filtered_data_no_d$peptide_start + 
  str_locate(filtered_data_no_d$Peptide, "N\\*")[,1] - 
  str_count(substr(filtered_data_no_d$Peptide, 1, str_locate(filtered_data_no_d$Peptide, "N\\*")[,1]), "\\*")

# Remove duplicates by summing TA_normalized_signal for reoccurring Sample, Glycan, protein_name, glycosylation_site pairs
filtered_data_no_d <- filtered_data_no_d %>%
  group_by(Sample, Glycan, protein_name, glycosylation_site) %>%
  summarize(
    total_signal = sum(total_signal, na.rm = TRUE), # Sum total_signal
    Peptide = first(Peptide), # Retain consistent information
    glycopeptide = first(glycopeptide),
    glycan_type = first(glycan_type),
    Hex = first(Hex),
    HexNAc = first(HexNAc),
    Neu5Ac = first(Neu5Ac),
    Fuc = first(Fuc),
    peptide_start = first(peptide_start),
    group = first(group), # Group consistent with Sample
    .groups = "drop" # Ungroup after summarizing
  )


# remove glycopeptides present in at least 2 M samples
group_m <- filtered_data_no_d[filtered_data_no_d$group == "M", ]
glycopeptide_counts <- table(group_m$glycopeptide) 
glycopeptides_to_remove <- names(glycopeptide_counts[glycopeptide_counts >= 2])
filtered_data_no_d <- filtered_data_no_d[!filtered_data_no_d$glycopeptide %in% glycopeptides_to_remove, ]

## remove hits from later statistical results that do not belong to known glycosylation sites
setwd("D:/BM/MCP_revision/260110_GP_reeva_214GPsfromprot_top34glycansfromglycomics/results")
glycosites_to_remove <- read.csv("glycosylation_sites_to_remove.csv", header = TRUE)
filtered_data_no_d <- filtered_data_no_d %>%
  anti_join(
    glycosites_to_remove,
    by = c("glycosylation_site", "protein_name")
  )

# total area normalization of the data ####
filtered_data_no_d%>%
  group_by(Sample)%>%
  summarise(sum(total_signal)) -> total_sample_intensity

for (i in 1:nrow(filtered_data_no_d)) {
  print(i)
  filtered_data_no_d$TA_normalized_signal[i] <- filtered_data_no_d$total_signal[i] / total_sample_intensity$`sum(total_signal)`[total_sample_intensity$Sample %in% filtered_data_no_d$Sample[i]]
}

# imputation of the data ####
# the purpose of this step is ONLY to provide a basis for the statistical comparison of glycopeptides
filtered_data_no_d%>%
  filter(group == "B")%>%
  group_by(glycopeptide)%>%
  summarise(n()) -> test1
filtered_data_no_d%>%
  filter(group == "A")%>%
  group_by(glycopeptide)%>%
  summarise(n()) -> test2
filtered_data_no_d%>%
  filter(group == "M")%>%
  group_by(glycopeptide)%>%
  summarise(n()) -> test3

unique_glycopeptides <- unique(filtered_data_no_d$glycopeptide) 
unique_peptides <- unique(filtered_data_no_d$Peptide) 
unique_proteins <- unique(filtered_data_no_d$protein_name) 

filtered_data_no_d_table <- filtered_data_no_d %>%
  select(glycopeptide, Peptide, Glycan, glycan_type, protein_name, Hex, HexNAc, Neu5Ac, Fuc, glycosylation_site, Sample, TA_normalized_signal) %>%
  pivot_wider(
    names_from = Sample,           
    values_from = TA_normalized_signal, 
    values_fill = NA             
  )
# Combine rows by glycopeptide while retaining all columns in filtered_data_no_d_table
filtered_data_no_d_table2 <- filtered_data_no_d_table %>%
  group_by(glycopeptide) %>%
  summarize(
    across(
      everything(),
      ~ if (is.character(.x)) {
        paste(na.omit(unique(.x)), collapse = "/")  
      } else if (is.numeric(.x)) {
        first(.x[!is.na(.x)])  
      } else {
        first(.x)  
      }
    ),
    .groups = "drop" 
  )
filtered_data_no_d_table2_long <- t(filtered_data_no_d_table2)

# Convert the first row to column names
filtered_data_no_d_table2_long <- as.data.frame(filtered_data_no_d_table2_long)  
colnames(filtered_data_no_d_table2_long) <- filtered_data_no_d_table2_long[1, ] 
filtered_data_no_d_table2_long <- filtered_data_no_d_table2_long[-1, ]          

df_glycopeptides <- filtered_data_no_d_table2_long %>% slice(10:25)
# df_glycopeptides[] <- lapply(df_glycopeptides, function(x) as.numeric(as.character(x)))

group <- c("A549", "A549", "A549", "A549", "A549", "A549", "BEAS-2B","BEAS-2B", "BEAS-2B", "BEAS-2B", "BEAS-2B", "BEAS-2B", "medium", "medium", "medium")

df_glycopeptides$group <- group
df_EV <- df_glycopeptides [1:12,]

# Filtering data ####
grouping_levels <- levels(as.factor(df_EV$group))
df_filtered <- df_EV
for (i in 1:ncol(df_EV)) {
  print(i)
  for (j in 1:length(grouping_levels)) {
    tmp_filtering <- df_filtered[df_filtered$group == grouping_levels[j],i]
    if (sum(is.na(tmp_filtering)) > length(tmp_filtering)/2) {
      df_filtered[df_filtered$group == grouping_levels[j],i] <- NA
    } else {}
  }
}
df_filtered <- df_filtered[,colSums(!is.na(df_filtered)) > 0]

# Imputation of data ####
df_EV_numeric <- df_EV[, 1:(ncol(df_EV)-1)]
df_EV_numeric[] <- lapply(df_EV_numeric, function(x) {
  if (is.factor(x)) {
    as.numeric(as.character(x))  
  } else if (is.character(x)) {
    suppressWarnings(as.numeric(x))  
  } else {
    x  
  }
})

sample_5perc <- apply(df_EV_numeric, 1, function(x) quantile(x, probs = 0.05, na.rm = TRUE))

#... First step of imputation: 5 quantile ####
for (i in 1:length(grouping_levels)) {
  print(paste("grouping", i, sep = "_"))
  tmp_grouping <- df_filtered[df_filtered$group == grouping_levels[i],]
  for (j in 1:(ncol(df_filtered)-1)) {
    print(paste("glycopeptide", j, sep = "_"))
    tmp_glycopeptide <- tmp_grouping[,j]
    tmp_5perc <- sample_5perc[row.names(df_EV) %in% row.names(tmp_grouping)]
    if (sum(is.na(tmp_glycopeptide)) > length(tmp_glycopeptide)/3) {
      tmp_glycopeptide[is.na(tmp_glycopeptide)] <- tmp_5perc[is.na(tmp_glycopeptide)]
      tmp_grouping[,j] <- tmp_glycopeptide
    } else {}}
  if (i == 1) {
    imputed_data_df <- tmp_grouping } else {
      tmp_df <- tmp_grouping
      imputed_data_df <- rbind(imputed_data_df, tmp_df)}
}

imputed_data_df[1:(ncol(imputed_data_df)-1)] <- lapply(imputed_data_df[1:(ncol(imputed_data_df)-1)], function(x) as.numeric(as.character(x)))
#... Second step of imputation: knn ####
full_matrix <- kNN(imputed_data_df, k = 15, trace = T, imp_var = F)
row.names(full_matrix) <- row.names(imputed_data_df)

setwd("D:/BM/MCP_revision/260110_GP_reeva_214GPsfromprot_top34glycansfromglycomics/results")
# write.csv(full_matrix, file = "imputed_full_matrix_251GP.csv", row.names = T)
# full_matrix <- read.csv("D:/ISTA_DIA_NN_241119_12sample6medium/data_eval/imputed_full_matrix_241120_12EV.csv")

# Log2 transformation ####
full_matrix[1:(ncol(full_matrix)-1)] <- lapply(full_matrix[1:(ncol(full_matrix)-1)], function(x) as.numeric(x))

num_cols <- sapply(full_matrix, is.numeric)
numeric_columns <- names(full_matrix)[num_cols]
full_matrix[,numeric_columns] <- log2(full_matrix[,numeric_columns])

# Statistical tests ####
#... A549 vs BEAS tests ####
#...... Normality tests ####
normality_df_group <- as.data.frame(matrix(rep(0, sum(num_cols)*length(levels(as.factor(full_matrix$group)))),
                                           ncol = length(numeric_columns)))
rownames(normality_df_group) <- levels(as.factor(full_matrix$group))
names(normality_df_group) <- numeric_columns

df_long <- gather(full_matrix, Glycopeptide, Value, all_of(numeric_columns))

for (i in 1:sum(num_cols)) {
  print(paste0("Glycopeptide", i))
  for (j in 1:length(levels(as.factor(full_matrix$group)))) {
    print(j)
    unlist(df_long%>%
             filter(Glycopeptide == numeric_columns[i])%>%
             filter(group == levels(as.factor(full_matrix$group))[j])%>%
             dplyr::select(Value)) -> normality_tmp
    shapiro.test(normality_tmp)[[2]] -> normality_df_group[j,i]
  }
}

#...... Equal variance tests ####
#......... Twogroup ####
combinations_group <- as.data.frame(combn(levels(as.factor(full_matrix$group)), 2, collapse=''))
variance_df_group_B <- as.data.frame(matrix(rep(0, sum(num_cols)*ncol(combinations_group)),
                                            ncol = ncol(combinations_group)))
rownames(variance_df_group_B) <- numeric_columns
names(variance_df_group_B) <- apply(combinations_group, 2, paste, collapse = " vs ")
for (i in 1:sum(num_cols)) {
  print(paste0("Comparison", i))
  for (j in 1:ncol(combinations_group)) {
    print(j)
    unlist(df_long%>%
             filter(Glycopeptide == numeric_columns[i])%>%
             filter(group == as.character(combinations_group[1,j]))%>%
             dplyr::select(Value)) -> tmp_variance_a
    unlist(df_long%>%
             filter(Glycopeptide == numeric_columns[i])%>%
             filter(group == as.character(combinations_group[2,j]))%>%
             dplyr::select(Value)) -> tmp_variance_b
    tmp_variance <- rbind(as.matrix(tmp_variance_a), as.matrix(tmp_variance_b))
    tmp_variance_grouping <- as.factor(rbind(as.matrix(rep("A", length(tmp_variance_a))), as.matrix(rep("B", length(tmp_variance_b)))))
    levene.test(tmp_variance, tmp_variance_grouping)[[2]] -> variance_df_group_B[i,j]
  }
}

#...... Two-group comparisons ####
# Student, Welch t-test or Wilcoxon rank sum test for the desired group combinations
twogroup_df_group <- as.data.frame(matrix(rep(0, ncol(combinations_group)*length(numeric_columns)),
                                          ncol = length(numeric_columns)))
names(twogroup_df_group) <- numeric_columns
rownames(twogroup_df_group) <- apply(combinations_group, 2, paste, collapse = " vs ")
twogroup_test_type_group <- twogroup_df_group
twogroup_fc_group <- twogroup_df_group

for (i in 1:sum(num_cols)) {
  print(paste0("Glycopeptide", i))
  for (j in 1:ncol(combinations_group)) {
    print(j)
    unlist(df_long%>%
             filter(Glycopeptide == numeric_columns[i])%>%
             filter(group == as.character(combinations_group[1,j]))%>%
             dplyr::select(Value)) -> tmp_ttest_a
    unlist(df_long%>%
             filter(Glycopeptide == numeric_columns[i])%>%
             filter(group == as.character(combinations_group[2,j]))%>%
             dplyr::select(Value)) -> tmp_ttest_b
    if (sum(normality_df_group[combinations_group[1,j],numeric_columns[i]] > 0.05,
            normality_df_group[combinations_group[2,j],numeric_columns[i]] > 0.05) == 2 &
        sum(variance_df_group_B[numeric_columns[i], j] > 0.05) == 1) {
      t.test(tmp_ttest_a, tmp_ttest_b, var.equal = T)[[3]] -> twogroup_df_group[j,i]
      mean(tmp_ttest_a) - mean(tmp_ttest_b) -> twogroup_fc_group[j,i]
      twogroup_test_type_group[j,i] <- "Student t-test"
    } else if (sum(normality_df_group[combinations_group[1,j],numeric_columns[i]] > 0.05,
                   normality_df_group[combinations_group[2,j],numeric_columns[i]] > 0.05) == 2) {
      t.test(tmp_ttest_a, tmp_ttest_b, var.equal = F)[[3]] -> twogroup_df_group[j,i]
      mean(tmp_ttest_a) - mean(tmp_ttest_b) -> twogroup_fc_group[j,i]
      twogroup_test_type_group[j,i] <- "Welch t-test"
    } else {
      wilcox.test(tmp_ttest_a, tmp_ttest_b, var.equal = F)[[3]] -> twogroup_df_group[j,i]
      mean(tmp_ttest_a) - mean(tmp_ttest_b) -> twogroup_fc_group[j,i]
      twogroup_test_type_group[j,i] <- "Wilcoxon-test"
    }
  }
}

twogroup_df_group_bh <- twogroup_df_group
twogroup_df_group_bh[1,] <- p.adjust(twogroup_df_group[1,], method = "BH")

statistic_results <- rbind("test_type"=twogroup_test_type_group, "p_value"=twogroup_df_group, "p_value_BH"=twogroup_df_group_bh, "FC"=twogroup_fc_group) %>% t()%>% as.data.frame()
significant_changes <- statistic_results[as.numeric(statistic_results$p_value_BH) < 0.05, ]

# write.csv(statistic_results, file = "statistic_results_289GP.csv")

significant_changes <- significant_changes %>%
  rownames_to_column(var = "glycopeptide")

significant_changes <- significant_changes %>%
  left_join(
    filtered_data_no_d_table2 %>% select(glycopeptide, glycan_type, glycosylation_site, protein_name), 
    by = "glycopeptide"
  )

statistic_results <- statistic_results %>%
  rownames_to_column(var = "glycopeptide")

statistic_results <- statistic_results %>%
  left_join(
    filtered_data_no_d_table2 %>% select(glycopeptide, glycan_type, glycosylation_site, protein_name), 
    by = "glycopeptide"
  )

# write.csv(statistic_results, file = "statistic_results_with_protein_names_and_loc.csv")
# write.csv(significant_changes, file = "significant_changes_with_protein_names_and_loc.csv")

#Visualization:
#pca ####
group_factor <- addNA(as.factor(full_matrix$group))

# # full matrix
mcdf <- full_matrix[,1:(ncol(full_matrix)-1)]
mcdf_matrix <- as.matrix(mcdf)

# #significant glycopeptides
# significant_glycopeptides <- significant_changes$glycopeptide
# filtered_full_matrix <- full_matrix[, c(significant_glycopeptides, "group")]
# mcdf <- filtered_full_matrix[,1:(ncol(filtered_full_matrix)-1)]
# mcdf_matrix <- as.matrix(mcdf)

pca <- prcomp(mcdf_matrix, scale. = T)

score.df <- data.frame(pca[[5]])

eigenvalues.vec<-pca[[1]]

pca_data <- data.frame(cbind(score.df$PC1, score.df$PC2))
row.names(pca_data) <- row.names(score.df)

# png(filename="PCA-of-all.png",
#     type="cairo",
#     units="in",
#     width=6,
#     height=6,
#     pointsize=1,
#     res=300)
ggplot(pca_data, aes(x = X1, y = X2)) +
  geom_point(aes(color = full_matrix$group), size = 3.5) +
  scale_color_manual(values = c("#F8766D","#7CAE00")) +
  scale_shape_manual(values = c(22,24,21,23)) +
  # labs(title = "PCA of different EV types") +
  xlab(paste("PC", 1, " (", round(eigenvalues.vec[1]*eigenvalues.vec[1] / sum(eigenvalues.vec*eigenvalues.vec)*100,2), "%)", sep="" )) +
  ylab(paste("PC", 2, " (", round(eigenvalues.vec[2]*eigenvalues.vec[2] / sum(eigenvalues.vec*eigenvalues.vec)*100,2), "%)", sep="" )) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        plot.margin = unit(c(10,15,10,15), "points"),
        plot.subtitle = element_text(colour = "black", size = 16, hjust = 0, family = "sans"),
        plot.tag = element_text(colour = "black", size = 20, family = "sans"),
        axis.title.x = element_text(color = "black", size = 18, family = "sans", vjust = -0.5),
        axis.text = element_text(color = "black", size = 16, family = "sans"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = "black", size = 18, family = "sans", vjust = 2.5),
        panel.background = element_rect(fill = "white"))
# dev.off()

#####################################################################
#heatmap ####
# #full matrix
# mcdf <- full_matrix[]
# heatmap_df <- as.data.frame(t(mcdf[,1:(ncol(mcdf)-1)]))
# heatmap_df <- as.data.frame(t(scale(t(heatmap_df), center = T, scale = T)))

# #significant glycopeptides
significant_glycopeptides <- significant_changes$glycopeptide
filtered_full_matrix <- full_matrix[, c(significant_glycopeptides, "group")]
mcdf <- filtered_full_matrix[]
heatmap_df <- as.data.frame(t(mcdf[,1:(ncol(mcdf)-1)]))
heatmap_df <- as.data.frame(t(scale(t(heatmap_df), center = T, scale = T)))

left.color <- "orange"
right.color <- "darkblue"
group_colors <- as.factor(full_matrix$group)
levels(group_colors) <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
breaks <- quantile(as.matrix(heatmap_df), probs = seq(0,1,0.025), na.rm = T)
ramp <- colorpanel(low  = left.color, high = right.color,n=40)
hm_function <- function(x) hclust(x, method = "ward.D2")

# png(filename="heatmap-of-sign.png",
#     type="cairo",
#     units="in",
#     width=8,
#     height=6,
#     pointsize=1,
#     res=300)
heatmap.2(as.matrix(heatmap_df),  col=ramp, trace = "none", breaks=breaks,
          key = T, density.info = "density", key.title = "", key.ylab = "",
          key.xlab = expression("Normalized intensity (" * italic(Z) * "-scored)"),
          hclustfun = hm_function, ColSideColors = as.character(group_colors),
          labCol = paste(full_matrix$group),
          cexCol = 1.5,  
          # cexRow = 1.5,  
          # cex.lab = 1.5,  
          margins = c(8,1),
          labRow = "", na.color = "darkgrey",
          key.par = list(cex.lab = 1.3))  
# dev.off()

##### Volcano plot ####
# Separate groups
volcano_data <- data.frame(
  Protein = colnames(full_matrix)[1:(ncol(full_matrix)-1)],  
  log2FC = as.numeric(statistic_results$FC),
  p_value = as.numeric(statistic_results$p_value_BH)
)
volcano_data$logP <- -log10(volcano_data$p_value)

volcano_data$expression <- ifelse(volcano_data$log2FC > 0  & volcano_data$p_value < 0.05, "Overexpressed",
                                  ifelse(volcano_data$log2FC < 0 & volcano_data$p_value < 0.05, "Underexpressed", "Not Significant"))

ggplot(volcano_data, aes(x = log2FC, y = logP)) +
  geom_point(aes(color = expression), size = 2.5) +
  scale_color_manual(values = c("Overexpressed" = "red", "Underexpressed" = "blue", "Not Significant" = "grey")) +  
  labs(x = expression(log[2]("Fold-Change")), 
       y = expression(-log[10]("adjusted" ~ italic(p) ~ "-value"))) +  
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed", size = 1) +
  # geom_vline(xintercept = c(-1, 1), col = "grey", linetype = "dashed") +
  annotate(
    "text",
    x = max(volcano_data$log2FC) * 0.9,  
    y = -log10(0.05),  
    label = expression(italic(p) == 0.05),  
    hjust = 0.25, vjust = -0.5,  
    size = 6, color = "black"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16, color = "black") 
  )

#sia-fuc-gal calculation
full_matrix_transposed <- t(full_matrix)

# Extract glycopeptide names (column names of the original matrix)
glycopeptide_names <- colnames(full_matrix)

# Create a matrix to store glycan compositions
tmp_glycan_monomers <- matrix(0, nrow = length(glycopeptide_names), ncol = 4)

# Apply glycan_monomer_extractor to each glycopeptide name
for (i in 1:length(glycopeptide_names)) {
  tmp_glycan_monomers[i, ] <- glycan_monomer_extractor(glycopeptide_names[i])
}

# Create a new data frame for the extracted glycan composition
glycan_composition <- data.frame(
  glycopeptide = glycopeptide_names,
  Hex = tmp_glycan_monomers[, 1],
  HexNAc = tmp_glycan_monomers[, 2],
  Neu5Ac = tmp_glycan_monomers[, 3],
  Fuc = tmp_glycan_monomers[, 4]
)


# # Add sample-specific data to the glycan composition table
sample_data <- as.data.frame(full_matrix_transposed)
sample_data <- sample_data[-nrow(sample_data),]
sample_data[] <- lapply(sample_data, function(x) as.numeric(as.character(x)))
sample_data <- 2^(sample_data)
sample_data <- rbind(sample_data, full_matrix_transposed[nrow(full_matrix_transposed), , drop = FALSE])

glycan_composition <- cbind(glycan_composition, sample_data)

glycan_composition[, 6:17] <- lapply(glycan_composition[, 6:17], function(x) as.numeric(as.character(x)))

glycan_composition[2:5] <- lapply(glycan_composition[2:5], function(x) {
  x[is.na(x)] <- 0
  x
})

glycan_composition <- glycan_composition %>%
  mutate(
    glycan_type = case_when(
      HexNAc == 2 & Hex >= 5                     ~ "oligomannose",
      Hex <= HexNAc + 1 & HexNAc >= 4         ~ "complex",
      !(HexNAc == 2 & Hex >= 5) & !(Hex <= HexNAc + 1 & HexNAc >= 4) & Hex > HexNAc + 1 ~ "hybrid",
      TRUE                                      ~ "misc"
    )
  )

glycan_composition$antennae <- ifelse(
  glycan_composition$glycan_type == "oligomannose",
  0,              # if oligomannose
  glycan_composition$HexNAc-2      # otherwise
)

glycan_composition$galactosylated_antennae <- ifelse(
  glycan_composition$glycan_type == "oligomannose",
  0,              # if oligomannose
  pmin(
    glycan_composition$Hex - 3,
    glycan_composition$HexNAc -2,
    na.rm = TRUE
  )      # otherwise
)

glycan_composition$sialylated_antennae <- ifelse(
  glycan_composition$glycan_type == "oligomannose",
  0,              # if oligomannose
  glycan_composition$Neu5Ac      # otherwise
)

glycan_ratios_df <- glycan_composition %>%
  pivot_longer(cols = 6:17, names_to = "Sample", values_to = "Signal") %>%
  group_by(Sample, glycan_type) %>%
  summarise(total_signal = sum(Signal, na.rm = TRUE), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(ratio = total_signal / sum(total_signal)) %>%
  select(Sample, glycan_type, ratio) %>%
  pivot_wider(names_from = Sample, values_from = ratio)

glycan_ratios_long <- glycan_ratios_df %>%
  pivot_longer(
    cols = -glycan_type,
    names_to = "Sample",
    values_to = "ratio"
  ) %>%
  pivot_wider(
    names_from = glycan_type,
    values_from = ratio
  )

idx <- glycan_composition$antennae > 0

weighted_neu5ac <- 
  colSums(glycan_composition[idx, colnames(sample_data)]*(glycan_composition$sialylated_antennae[idx]/glycan_composition$antennae[idx]),na.rm = TRUE) /
  colSums(glycan_composition[idx, colnames(sample_data)],na.rm = TRUE)

weighted_gal <- 
  colSums(glycan_composition[idx, colnames(sample_data)]*(glycan_composition$galactosylated_antennae[idx]/glycan_composition$antennae[idx]),na.rm = TRUE) /
  colSums(glycan_composition[idx, colnames(sample_data)],na.rm = TRUE)

weighted_fuc <- colSums(glycan_composition[, colnames(sample_data)] * (glycan_composition$Fuc), na.rm=TRUE) / 
  colSums(glycan_composition[, colnames(sample_data)], na.rm = TRUE)


average_glycan_df <- data.frame(
  # weighted_gal,
  weighted_neu5ac,
  weighted_fuc
)
average_glycan_df <- average_glycan_df %>%
  tibble::rownames_to_column("Sample")
average_glycan_df <- average_glycan_df %>%
  left_join(glycan_ratios_long, by = "Sample")
rownames(average_glycan_df) <- average_glycan_df$Sample
average_glycan_df$Sample <- NULL
average_glycan_df$misc <- NULL

print(average_glycan_df)

imputed_matrix_metrics <- data.frame(
  Metric = colnames(average_glycan_df),
  p_value = numeric(ncol(average_glycan_df)),
  test_type = character(ncol(average_glycan_df)),
  FC = numeric(ncol(average_glycan_df)),
  p_value_BH = numeric(ncol(average_glycan_df)),
  stringsAsFactors = FALSE
)

# Perform two-group tests for each metric
for (i in seq_along(colnames(average_glycan_df))) {
  metric <- colnames(average_glycan_df)[i]
  
  # Split data into Group A and Group B
  group_A <- average_glycan_df[1:6, metric]
  group_B <- average_glycan_df[7:12, metric]
  
  # Perform the two-group test
  test_result <- two_group_test(group_A[!is.na(group_A)], group_B[!is.na(group_B)], log_data = FALSE)
  
  # Store the p-value in the results data frame
  imputed_matrix_metrics[i, "p_value"] <- test_result[[1]]
  imputed_matrix_metrics[i, "test_type"] <- test_result[[2]] 
  imputed_matrix_metrics[i, "FC"] <- test_result[[3]] 
}

imputed_matrix_metrics$p_value_BH <- p.adjust(imputed_matrix_metrics$p_value, method = "BH")
# write.csv(imputed_matrix_metrics, file = "imputed_matrix_metrics.csv")

group2 <- c("A549", "A549","A549", "A549","A549", "A549", "BEAS-2B", "BEAS-2B", "BEAS-2B", "BEAS-2B","BEAS-2B", "BEAS-2B")
average_glycan_df$group <- group2

ggplot(average_glycan_df, aes(x = group, y = weighted_neu5ac, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#F8766D","#7CAE00")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  # scale_y_log10() +  # Log10 scale for y-axis
  # facet_wrap(~Component, scales = "free_y") +  # Facet by component
  labs(
    # title = "Fucosylation",
    # x = "Group",
    # y = "Ratio of fucosylated glycans",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # axis.title.x = element_text(size = 18), 
    # axis.title.y = element_text(size = 20, color = "black"),  
    axis.text.x = element_text(size = 18, color = "black"),   
    axis.text.y = element_text(size = 18, color = "black"),  
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 14)    
  )
# write.csv(average_glycan_df,"average_glycan_df.csv")

