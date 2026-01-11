setwd("D:/BM/MCP_revision/260110_GP_reeva_214GPsfromprot_top34glycansfromglycomics/results")
filename1<-"P-GP-corr.csv"
dir<-getwd()
fullpath1<-file.path(dir,filename1)
data <- read.csv(fullpath1,header = TRUE)

# Count the number of proteins associated with each glycan
filtered_data <- data %>%
  group_by(Glycan) %>%
  filter(n() >= 5) %>% # Keep only glycans with 5 or more associated proteins
  # group_by(Protein.name) %>%
  # filter(n() >= 3) %>%    # Keep only sites occurring at least 3x
  group_by(Glycosylation_site) %>%
  filter(n() >= 2) %>%    # Keep only sites occurring at least 3x
  ungroup()

# Add a column to classify positive and negative log values
filtered_data$Sign <- ifelse(filtered_data$log2.GP_norm > 0, "Positive", "Negative")

# Plot the bubble plot
ggplot(filtered_data, aes(x = Glycosylation_site, y = Glycan, size = abs(log2.GP_norm), color=Sign)) +
  geom_point(alpha = 0.7) +
  scale_size_continuous(name = "Log(2) normalized FC") +
  scale_color_manual(name = "Sign", 
                       values = c("Positive" = "red", "Negative" = "blue")) + 
  theme_minimal() +
  labs(
    # title = "Protein-Glycan Overlap",
    x = "Glycosylation site",
    y = "Glycan"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "gray90")
  )




