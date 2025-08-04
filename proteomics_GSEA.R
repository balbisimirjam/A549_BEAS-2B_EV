library(dplyr)
library(MSnID)

df_res <- read.csv("D:/ISTA_DIA_NN_241119_12sample6medium/data_eval_2_min2peptides/statistic_results_for_gsea.csv")
# Add ENTREZID column

conv_tbl <- fetch_conversion_table(organism_name = "Homo sapiens", 
                                  from = "SYMBOL", to = "ENTREZID")
# Merge ENTREZID if required (only if Gene identifiers are not ENTREZID)
df_res <- df_res %>%
  left_join(conv_tbl, by = c("Gene" = "SYMBOL")) 

# Filter out genes without valid ENTREZID or log2FC
geneList <- df_res %>%
  filter(!is.na(ENTREZID), !is.na(log2FC)) %>% 
  mutate(ranking_metric = -log10(adj_p_value) * sign(log2FC)) %>% 
  group_by(ENTREZID) %>% 
  summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>% 
  arrange(-ranking_metric) %>% 
  tibble::deframe() 

# # Check resulting geneList
# head(geneList)
# tail(geneList)


##### clusterProfiler #####
## GO-BP FGSEA with clusterProfiler package
library(clusterProfiler)
set.seed(123)
system.time( # keep track of elapsed time
  cgsea_res <- gseGO(geneList = geneList,
                     ont = "BP", #BP-Biological Process, MF-Molecular Function, CC-cellular Component
                     OrgDb = "org.Hs.eg.db",
                     minGSSize = 15,
                     maxGSSize = 500,
                     pvalueCutoff = 0.1,
                     eps = 0,
                     nPermSimple = 10000,
                     seed = TRUE)
)

# First 8 rows with lowest enrichment p-values
cgsea_res@result %>%
  arrange(NES) %>%
  head(8)

go_results <- data.frame(
  GO = rownames(cgsea_res@result),  
  NES = cgsea_res@result$NES,       
  p.adjust = cgsea_res@result$p.adjust
)
library(GO.db)
go_descriptions <- Term(go_results$GO)
go_results$Description <- go_descriptions

# write.csv(go_results, file = "go_results_cutoff_0.1.csv", row.names = T)

##### VISUALIZATION #####
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html  
# dotplot: https://github.com/junjunlab/GseaVis/wiki/dotplotGsea 

#BiocManager::install("enrichplot")

library(enrichplot)

dotplot(cgsea_res, showCategory=10) + ggtitle("dotplot for GSEA") 

edox <- setReadable(cgsea_res, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, node_label="category",   
         showCategory = 25,    
         color_category='firebrick', 
         color_gene='steelblue',
         cex_label_category = 0.7) 

heatplot(edox, foldChange=geneList, showCategory=15)

edox2 <- pairwise_termsim(edox)
treeplot(edox2, showCategory = 20, label_params = list(size = 1)) #color = "NES", 

cgsea_res_pairwise <- pairwise_termsim(cgsea_res)
emapplot(cgsea_res_pairwise)
emapplot(
  edox2,
  showCategory = 75,  
  color = "NES",      
  layout = "nicely"   
)
