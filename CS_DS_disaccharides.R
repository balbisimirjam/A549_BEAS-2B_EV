library(tidyverse)
library(ggplot2)
library(ggthemes)
library(gplots)
library(tidyr)
library(scales)
library(dplyr)
library(lawstat)
library(rstatix)

setwd("D:/GAG_EV_12sample_6medium")
filename1<-"AB_EV_CS.csv"
dir<-getwd()
fullpath1<-file.path(dir,filename1)
cs<-read_csv(fullpath1)
cs<-cs[-6,] ##Sample 6 needs to be removed (no disaccharides detected)

#Boxplots
# Generate boxplots for total amount of each components
ggplot(cs, aes(x = group, y = six_four_ratio, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#F8766D","#7CAE00")) +
  labs(
    # title = "Boxplots for CS disaccharides",
    # x = "Group",
    # y = "Total amount [fmol]",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 18, color = "black"),   
    axis.text.y = element_text(size = 18, color = "black"),   
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 14)    
  )

# Total amount of each components ####
# Reshape data from wide to long format
cs_long <- cs %>%
  pivot_longer(
    cols = c(D0a0, D0a4, D0a6, D0a10),  
    names_to = "Component",             
    values_to = "Value"                 
  )

cs_long <- cs_long %>%
  mutate(
    Component = factor(Component, levels = c("D0a0", "D0a4", "D0a6", "D0a10"))
  )

# Generate boxplots for total amount of each components
ggplot(cs_long, aes(x = group, y = Value, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#F8766D","#7CAE00")) +
  # scale_y_log10() +  # Log10 scale for y-axis
  facet_wrap(~Component, scales = "free_y") + 
  labs(
    title = "Boxplots for CS disaccharides",
    x = "Group",
    y = "Total amount [fmol]",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 10)
  )

# Relative amount of each components ####
# Reshape data from wide to long format
cs_long_rel <- cs %>%
  pivot_longer(
    cols = c(D0a0_rel, D0a4_rel, D0a6_rel, D0a10_rel),  
    names_to = "Component",             
    values_to = "Value"                
  )

cs_long_rel <- cs_long_rel %>%
  mutate(
    Component = str_remove(Component, "_rel")
  ) %>% 
  mutate(
    Component = factor(Component, levels = c("D0a0", "D0a4", "D0a6", "D0a10"))
  ) %>%
  mutate(
    Value = Value * 100  
  )

# Generate boxplots for relative amount of each components
ggplot(cs_long_rel, aes(x = group, y = Value, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#F8766D","#7CAE00")) +
  # scale_y_log10() +  # Log10 scale for y-axis
  facet_wrap(~Component, scales = "free_y") +  
  labs(
    # title = "Boxplots for CS disaccharides",
    x = "Group",
    y = "Relative amount (%)",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    # plot.title = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_blank(),
    # axis.title.y = element_blank(),
    # axis.title.x = element_text(size = 18),  
    axis.title.y = element_text(size = 20, color = "black"),  
    axis.text.x = element_text(size = 18, color = "black"),   
    axis.text.y = element_text(size = 18, color = "black"),   
    legend.title = element_text(size = 16),  
    legend.text = element_text(size = 14),    
    strip.text = element_text(size = 18, face = "bold", color = "black")  
  )

#####
#Total amount of CS disaccharides
ggplot(cs, aes(x = group, y = CS_total, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#F8766D","#7CAE00")) +
  # scale_y_log10() +
  labs(
    x = "Group",
    y = "Total amount of CS disaccharides [fmol]",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 10)
  )

#Average sulfation
ggplot(cs, aes(x = group, y = average_CS_sulf, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#F8766D","#7CAE00")) +
  # scale_y_log10() +
  labs(
    x = "Group",
    y = "Average rate of CS sulfation",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 10)
  )

#6S/4S ratio
ggplot(cs, aes(x = group, y = six_four_ratio, fill = group)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#F8766D","#7CAE00")) +
  # scale_y_log10() +
  labs(
    x = "Group",
    y = "6S/4S ratio",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 10)
  )

#barplots ####
cs_summary <- cs %>%
  group_by(group) %>%
  summarise(
    CS_total_mean = mean(CS_total, na.rm = TRUE),
    CS_total_sd = sd(CS_total, na.rm = TRUE),
    average_CS_sulf_mean = mean(average_CS_sulf, na.rm = TRUE),
    average_CS_sulf_sd = sd(average_CS_sulf, na.rm = TRUE),
    six_four_ratio_mean = mean(six_four_ratio, na.rm = TRUE),
    six_four_ratio_sd = sd(six_four_ratio, na.rm = TRUE)
  )

# Total amount of CS disaccharides
ggplot(cs_summary, aes(x = group, y = CS_total_mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7) +  
  geom_errorbar(aes(
    ymin = CS_total_mean - CS_total_sd,
    ymax = CS_total_mean + CS_total_sd
  ), width = 0.2, color = "black") +  
  scale_fill_manual(values = c("#F8766D", "#7CAE00")) +
  labs(
    x = "Group",
    y = "Total amount of CS disaccharides",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 10)
  )

#Average sulfation
ggplot(cs_summary, aes(x = group, y = average_CS_sulf_mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7) +  
  geom_errorbar(aes(
    ymin = average_CS_sulf_mean - average_CS_sulf_sd,
    ymax = average_CS_sulf_mean + average_CS_sulf_sd
  ), width = 0.2, color = "black") +  # Error bars
  scale_fill_manual(values = c("#F8766D", "#7CAE00")) +
  labs(
    x = "Group",
    y = "Average rate of CS sulfation",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 10)
  )

#6S/4S ratio
ggplot(cs_summary, aes(x = group, y = six_four_ratio_mean, fill = group)) +
  geom_bar(stat = "identity", width = 0.7) +  # Use "identity" for actual values
  geom_errorbar(aes(
    ymin = six_four_ratio_mean - six_four_ratio_sd,
    ymax = six_four_ratio_mean + six_four_ratio_sd
  ), width = 0.2, color = "black") +  # Error bars
  scale_fill_manual(values = c("#F8766D", "#7CAE00")) +
  labs(
    x = "Group",
    y = "6S/4S ratio",
    fill = "Group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.text.x = element_text(size = 10)
  )

#######################################################################
#statistics 
numeric_columns <- names(cs)[2:ncol(cs)]

#...... Normality tests ####
normality_df_group <- as.data.frame(matrix(rep(0, length(numeric_columns) * length(unique(cs$group))),
                                           ncol = length(numeric_columns)))
rownames(normality_df_group) <- unique(cs$group)
names(normality_df_group) <- numeric_columns

df_long <- cs %>%
  gather(Key, Value, all_of(numeric_columns))

# Perform Shapiro-Wilk test for each group and numeric column
for (i in numeric_columns) {
  for (j in unique(cs$group)) {
    # Extract values for the given column and group
    normality_tmp <- df_long %>%
      filter(Key == i, group == j) %>%
      pull(Value)
    
    # Run the Shapiro-Wilk test if there are enough values
    if (length(normality_tmp) > 2) { 
      normality_df_group[j, i] <- shapiro.test(normality_tmp)$p.value
    } else {
      normality_df_group[j, i] <- NA 
    }
  }
}

#...... Equal variance tests ####
#......... Twogroup ####
combinations_group <- as.data.frame(combn(levels(as.factor(cs$group)), 2, collapse=''))
variance_df_group_B <- as.data.frame(matrix(rep(0, length(numeric_columns) * ncol(combinations_group)),
                                            ncol = ncol(combinations_group)))
rownames(variance_df_group_B) <- numeric_columns
names(variance_df_group_B) <- apply(combinations_group, 2, paste, collapse = " vs ")

# Perform Levene's test for each numeric column and group pair
for (i in numeric_columns) {
  print(paste0("Testing variance for column: ", i))
  for (j in 1:ncol(combinations_group)) {
    # Extract group pair
    group1 <- as.character(combinations_group[1, j])
    group2 <- as.character(combinations_group[2, j])
    
    # Extract values for the two groups
    tmp_variance_a <- cs %>%
      filter(group == group1) %>%
      pull(i)
    
    tmp_variance_b <- cs %>%
      filter(group == group2) %>%
      pull(i)
    
    # Combine values and group labels
    tmp_variance <- c(tmp_variance_a, tmp_variance_b)
    tmp_variance_grouping <- factor(c(rep("A", length(tmp_variance_a)), rep("B", length(tmp_variance_b))))
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

# Perform comparisons for each numeric column and group pair
for (i in seq_along(numeric_columns)) {
  for (j in 1:ncol(combinations_group)) {
    # Subset data for the two groups being compared
    tmp_ttest_a <- cs %>%
      filter(group == combinations_group[1, j]) %>%
      pull(numeric_columns[i])
    
    tmp_ttest_b <- cs %>%
      filter(group == combinations_group[2, j]) %>%
      pull(numeric_columns[i])
    
    # Determine the test based on normality and variance
    if (sum(normality_df_group[combinations_group[1, j], numeric_columns[i]] > 0.05,
            normality_df_group[combinations_group[2, j], numeric_columns[i]] > 0.05) == 2 &&
        variance_df_group_B[numeric_columns[i], j] > 0.05) {
      # Student's t-test
      t.test(tmp_ttest_a, tmp_ttest_b, var.equal = TRUE)$p.value -> twogroup_df_group[j, i]
      mean(tmp_ttest_a) / mean(tmp_ttest_b) -> twogroup_fc_group[j, i]
      twogroup_test_type_group[j, i] <- "Student t-test"
    } else if (sum(normality_df_group[combinations_group[1, j], numeric_columns[i]] > 0.05,
                   normality_df_group[combinations_group[2, j], numeric_columns[i]] > 0.05) == 2) {
      # Welch's t-test
      t.test(tmp_ttest_a, tmp_ttest_b, var.equal = FALSE)$p.value -> twogroup_df_group[j, i]
      mean(tmp_ttest_a) / mean(tmp_ttest_b) -> twogroup_fc_group[j, i]
      twogroup_test_type_group[j, i] <- "Welch t-test"
    } else {
      # Wilcoxon rank-sum test
      wilcox.test(tmp_ttest_a, tmp_ttest_b)$p.value -> twogroup_df_group[j, i]
      mean(tmp_ttest_a) / mean(tmp_ttest_b) -> twogroup_fc_group[j, i]
      twogroup_test_type_group[j, i] <- "Wilcoxon-test"
    }
  }
}
  
twogroup_df_group_bh <- twogroup_df_group
twogroup_df_group_bh[1,] <- p.adjust(twogroup_df_group[1,], method = "BH")

statistic_results <- rbind("test_type"=twogroup_test_type_group, "p_value"=twogroup_df_group, "p_value_BH"=twogroup_df_group_bh, "FC"=twogroup_fc_group) %>% t()%>% as.data.frame()
significant_changes <- statistic_results[as.numeric(statistic_results$p_value_BH) < 0.05, ]

# write.csv(statistic_results, file = "statistic_results.csv")

#PCA, clustering ####

# gag data heatmap ####

mcdf <- cs[,c(2:5)] #total amounts
# mcdf <- cs[,c(6:9)] %>%   rename_with(~ gsub("_rel$", "", .), everything()) #relative amounts

heatmap_df <- as.data.frame(t(mcdf))
heatmap_df <- as.data.frame(t(scale(t(heatmap_df), center = T, scale = T)))

left.color <- "orange"
right.color <- "darkblue"
group_colors <- as.factor(cs$group)
levels(group_colors) <- c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
breaks <- quantile(as.matrix(heatmap_df), probs = seq(0,1,0.025), na.rm = T)
ramp <- colorpanel(low  = left.color, high = right.color,n=40)
hm_function <- function(x) hclust(x, method = "ward.D2")

heatmap.2(as.matrix(heatmap_df),  col=ramp, trace = "none", breaks=breaks,
          key = T, density.info = "density", key.title = "", key.ylab = "",
          key.xlab = expression("Intensity (" * italic(Z) * "-scored)"),
          hclustfun = hm_function, ColSideColors = as.character(group_colors),
          labCol = paste(cs$group),
          cexCol = 1.5,  
          # cexRow = 1.5,  
          # cex.lab = 1.5,  
          margins = c(8,8),
          labRow = colnames(mcdf), na.color = "darkgrey",
          key.par = list(cex.lab = 1.3))  
#####################################################################

# gag data pca ####

# mcdf <- cs[,c(2:5)] #total amounts
mcdf <- cs[,c(6:9)] #relative amounts

mcdf_matrix <- as.matrix(mcdf)

pca <- prcomp(mcdf_matrix, scale. = T)

score.df <- data.frame(pca[[5]])

eigenvalues.vec<-pca[[1]]

pca_data <- data.frame(cbind(score.df$PC1, score.df$PC2))
row.names(pca_data) <- row.names(score.df)

png(filename="PCA_EV_gag_total.png", 
    type="cairo",
    units="in", 
    width=6, 
    height=6, 
    pointsize=1, 
    res=300)
ggplot(pca_data, aes(x = X1, y = X2, color=cs$group)) +
  geom_point(size = 5) +
  scale_color_manual(values = c("#F8766D","#7CAE00","#00BFC4","#C77CFF")) +  
  scale_shape_manual(values = c(22,24,21,23)) +
  # labs(title = "PCA of different types of lung cancer samples") +
  xlab(paste("PC", 1, " (", round(eigenvalues.vec[1]*eigenvalues.vec[1] / sum(eigenvalues.vec*eigenvalues.vec)*100,2), "%)", sep="" )) +
  ylab(paste("PC", 2, " (", round(eigenvalues.vec[2]*eigenvalues.vec[2] / sum(eigenvalues.vec*eigenvalues.vec)*100,2), "%)", sep="" )) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
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
