library(tidyverse)
library(VIM)
library(lawstat)
library(rstatix)
library(ggpattern)
library(gplots)
library(ggplot2)
library(dplyr)

setwd("D:/ISTA_DIA_NN_241119_12sample6medium/data_eval")
# Reading in data ####
df_original <- read.delim("241120_report.pg_matrix.tsv")
df_raw <- as.data.frame(t(df_original))
names(df_raw) <- df_original$Protein.Group
df_raw <- df_raw [-1,]
df_raw[] <- lapply(df_raw, function(x) as.numeric(as.character(x)))

df_raw[df_raw == 0] <- NA

# Remove proteins with 1 peptide
df_test <- read.delim("D:/ISTA_DIA_NN_241119_12sample6medium/search/report.tsv")
data <- unique(df_test[,c('Protein.Group','Stripped.Sequence')])
protein_peptide_counts <- data %>%
  group_by(Protein.Group) %>%  
  summarise(unique_peptides = n_distinct(Stripped.Sequence), .groups = 'drop')

proteins_with_count_min2 <- protein_peptide_counts %>%
  filter(unique_peptides > 1) 

protein_list <- proteins_with_count_min2$Protein.Group

filtered_df_raw <- df_raw %>%
  select(any_of(protein_list))

#Grouping
group <- c("A549", "A549", "A549", "A549", "A549", "A549", "BEAS-2B", "BEAS-2B", "BEAS-2B", "BEAS-2B", "BEAS-2B", "BEAS-2B", "medium", "medium", "medium", "medium", "medium", "medium")

filtered_df_raw$group <- group
df_EV <- filtered_df_raw [1:12,]

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
sample_5perc <- sapply(as.data.frame(t(df_EV[,1:1739])), quantile, p = 0.05, na.rm = T)

#... First step of imputation: 5 quantile ####
for (i in 1:length(grouping_levels)) {
  print(paste("grouping", i, sep = "_"))
  tmp_grouping <- df_filtered[df_filtered$group == grouping_levels[i],]
  for (j in 1:(ncol(df_filtered)-1)) {
    print(paste("protein", j, sep = "_"))
    tmp_protein <- tmp_grouping[,j]
    tmp_5perc <- sample_5perc[row.names(df_EV) %in% row.names(tmp_grouping)]
    if (sum(is.na(tmp_protein)) > length(tmp_protein)/3) {
      tmp_protein[is.na(tmp_protein)] <- tmp_5perc[is.na(tmp_protein)]
      tmp_grouping[,j] <- tmp_protein
    } else {}}
  if (i == 1) {
    imputed_data_df <- tmp_grouping } else {
    tmp_df <- tmp_grouping
    imputed_data_df <- rbind(imputed_data_df, tmp_df)}
}

#... Second step of imputation: knn ####
full_matrix <- kNN(imputed_data_df, k = 15, trace = T, imp_var = F)
row.names(full_matrix) <- row.names(imputed_data_df)
setwd("D:/ISTA_DIA_NN_241119_12sample6medium/data_eval_2_min2peptides")
# write.csv(full_matrix, file = "imputed_full_matrix_241210_12EV.csv", row.names = T)
# full_matrix <- read.csv("D:/ISTA_DIA_NN_241119_12sample6medium/data_eval/imputed_full_matrix_241120_12EV.csv")

# Log2 transformation ####
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

df_long <- gather(full_matrix, Protein, Value, all_of(numeric_columns))

for (i in 1:sum(num_cols)) {
  print(paste0("Protein", i))
  for (j in 1:length(levels(as.factor(full_matrix$group)))) {
    print(j)
    unlist(df_long%>%
             filter(Protein == numeric_columns[i])%>%
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
             filter(Protein == numeric_columns[i])%>%
             filter(group == as.character(combinations_group[1,j]))%>%
             dplyr::select(Value)) -> tmp_variance_a
    unlist(df_long%>%
             filter(Protein == numeric_columns[i])%>%
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
  print(paste0("Protein", i))
  for (j in 1:ncol(combinations_group)) {
    print(j)
    unlist(df_long%>%
             filter(Protein == numeric_columns[i])%>%
             filter(group == as.character(combinations_group[1,j]))%>%
             dplyr::select(Value)) -> tmp_ttest_a
    unlist(df_long%>%
             filter(Protein == numeric_columns[i])%>%
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

# write.csv(statistic_results, file = "statistic_results.csv")

# Filter cell culture medium components from full_matrix, statistic_results and significant_changes ####
df_medium <- filtered_df_raw [13:18,]
fun_x <- function(x) sum(!is.na(x))
filtering_boolean <- apply(df_medium, MARGIN = 2, FUN = fun_x) > 1
df_medium <- df_medium[, filtering_boolean] 

columns_to_remove <- colnames(full_matrix) %in% colnames(df_medium)
full_matrix_no_medium_comp <- full_matrix[, !columns_to_remove]
full_matrix_no_medium_comp$group <- full_matrix$group
# write.csv(full_matrix_no_medium_comp, file = "imputed_full_matrix_no_medium_945.csv", row.names = T)

rows_to_remove <- rownames(statistic_results) %in% colnames(df_medium)
statistic_results_no_medium_comp <- statistic_results[!rows_to_remove,]
significant_changes_no_medium_comp <- statistic_results_no_medium_comp[as.numeric(statistic_results_no_medium_comp$p_value_BH) < 0.05, ]

# write.csv(statistic_results_no_medium_comp, file = "statistic_results_no_medium_comp.csv")

###plots ####
#pca ####
group_factor <- addNA(as.factor(df_EV$group))

# # full matrix (945 proteins)
mcdf <- full_matrix_no_medium_comp[,1:945]
mcdf_matrix <- as.matrix(mcdf)

# #significant proteins (408 proteins)
# significant_proteins <- rownames(significant_changes_no_medium_comp)
# filtered_full_matrix_no_medium_comp <- full_matrix_no_medium_comp[, c(significant_proteins, "group")]
# mcdf <- filtered_full_matrix_no_medium_comp[,1:408]
# mcdf_matrix <- as.matrix(mcdf)

pca <- prcomp(mcdf_matrix, scale. = T)

score.df <- data.frame(pca[[5]])

eigenvalues.vec<-pca[[1]]

pca_data <- data.frame(cbind(score.df$PC1, score.df$PC2))
row.names(pca_data) <- row.names(score.df)

png(filename="PCA.png",
    type="cairo",
    units="in",
    width=6,
    height=6,
    pointsize=1,
    res=300)
ggplot(pca_data, aes(x = X1, y = X2)) +
  geom_point(aes(color = df_EV$group), size = 3.5) +
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
dev.off()

#####################################################################
#heatmap ####
# #full matrix (945 proteins)
# mcdf <- full_matrix_no_medium_comp[]
# heatmap_df <- as.data.frame(t(mcdf[,1:(ncol(mcdf)-1)]))
# heatmap_df <- as.data.frame(t(scale(t(heatmap_df), center = T, scale = T)))

#significant proteins (408 proteins)
significant_proteins <- rownames(significant_changes_no_medium_comp)
filtered_full_matrix_no_medium_comp <- full_matrix_no_medium_comp[, c(significant_proteins, "group")]
mcdf <- filtered_full_matrix_no_medium_comp[]
heatmap_df <- as.data.frame(t(mcdf[,1:(ncol(mcdf)-1)]))
heatmap_df <- as.data.frame(t(scale(t(heatmap_df), center = T, scale = T)))

left.color <- "orange"
right.color <- "darkblue"
group_colors <- as.factor(df_EV$group)
levels(group_colors) <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
breaks <- quantile(as.matrix(heatmap_df), probs = seq(0,1,0.025), na.rm = T)
ramp <- colorpanel(low  = left.color, high = right.color,n=40)
hm_function <- function(x) hclust(x, method = "ward.D2")

# png(filename="heatmap_all.png",
#     type="cairo",
#     units="in",
#     width=9.15,
#     height=6,
#     pointsize=1,
#     res=300)
heatmap.2(as.matrix(heatmap_df),  col=ramp, trace = "none", breaks=breaks,
          key = T, density.info = "density", key.title = "", key.ylab = "",
          key.xlab = expression("LFQ intensity (" * italic(Z) * "-scored)"),
          hclustfun = hm_function, ColSideColors = as.character(group_colors),
          labCol = paste(df_EV$group),
          cexCol = 1.5,  # Increase column label size
          # cexRow = 1.5,  # Increase column label size
          # cex.lab = 1.5,  # Increase legend text size
          margins = c(8,1),
          labRow = "", na.color = "darkgrey",
          key.par = list(cex.lab = 1.5))  
# dev.off()

##### Volcano plot ####
# Separate groups
volcano_data <- data.frame(
  Protein = colnames(full_matrix_no_medium_comp)[1:945],  
  log2FC = as.numeric(statistic_results_no_medium_comp$FC),
  p_value = as.numeric(statistic_results_no_medium_comp$p_value_BH)
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
    axis.text = element_text(size = 16, color = "black"),
    plot.margin = margin(5, 15, 5, 5)  # Increase right-side margin
  )

#PG core proteins boxplot ####
pg_core_proteins <- read.csv("D:/ISTA_DIA_NN_241119_12sample6medium/data_eval/PG_identified.csv")

# Ensure pg_core_proteins is a data frame and has a PG_protein column
pg_proteins_present <- pg_core_proteins$PG_protein[
  pg_core_proteins$PG_protein %in% colnames(full_matrix_no_medium_comp) 
]

# Subset the matrix
pg_matrix <- full_matrix_no_medium_comp[, pg_proteins_present, drop = FALSE]

# Add grouping information
pg_matrix$group <- full_matrix_no_medium_comp$group

# Reshape the data to long format
library(tidyr)
long_data <- pivot_longer(
  pg_matrix,
  cols = -group,  # Exclude the 'group' column
  names_to = "Protein",
  values_to = "Expression"
)

# Create the boxplots
library(ggplot2)

# Define the output directory
output_dir <- "individual_boxplots" 
dir.create(output_dir, showWarnings = FALSE)  

# Loop through each protein and save a boxplot
for (protein in unique(long_data$Protein)) {
  # Subset data for the specific protein
  protein_data <- subset(long_data, Protein == protein)

  # Create the plot
  p <- ggplot(protein_data, aes(x = group, y = Expression, fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("#F8766D","#7CAE00")) +
    scale_y_log10() +
    labs(
      # title = paste("Boxplot for", protein),
      # x = "Group",
      y = expression(log[2]("LFQ intensity")),
      fill = "Group"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.title.x = element_blank(),
      # axis.title.x = element_text(size = 18),  
      axis.title.y = element_text(size = 20, color = "black"),  
      axis.text.x = element_text(size = 18, color = "black"),   
      axis.text.y = element_text(size = 18, color = "black"),   
      legend.title = element_text(size = 16),  
      legend.text = element_text(size = 14)    
    )

  # Save the plot
  ggsave(
    filename = file.path(output_dir, paste0(protein, "_boxplot.png")),
    plot = p,
    width = 8, height = 6
  )
}
