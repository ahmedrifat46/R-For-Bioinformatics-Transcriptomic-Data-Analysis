gc()  # Clear memory to free up space before running analysis

# Check if BiocManager is installed; install if missing
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install CRAN packages
install.packages(c("dplyr", "tibble", "ggplot2"))
#Install clusterProfiler and org.Hs.eg.db from Bioconductor
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("enrichplot")
BiocManager::install("msigdbr")
BiocManager::install("ggupset")

# Load Required Libraries
library(clusterProfiler)  # Main package for functional enrichment analysis
library(org.Hs.eg.db)     # Annotation database for human genes
library(enrichplot)       # Visualization of enrichment results
library(msigdbr)          # Access MSigDB collections such as Hallmark, KEGG, etc.
library(dplyr)            # Data manipulation
library(tibble)
library(tidyr)
library(ggplot2)
library(ggupset)

# Prepare Gene Set Identifiers
# ----------------------------
load("deg_results_GSE79973.RData")

symbols <- as.character(rownames(deg_results))

# hgu133plus2.db

entrezid <- mapIds(x = org.Hs.eg.db,
                   keys = symbols,
                   column = "ENTREZID",
                   keytype = "SYMBOL",
                   multiVals = "first") %>%
  stack() %>%
  dplyr::rename(Entrezid = values, Symbol = ind)

deg_results <- deg_results %>%
  tibble::rownames_to_column("Symbol") %>%
  left_join(entrezid, by = "Symbol") %>%
  relocate(Entrezid, .before = logFC)

# Remove entries without IDs or duplicate mappings
deg_results <- deg_results %>%
  drop_na(Entrezid, adj.P.Val) %>%
  arrange(adj.P.Val) %>%
  distinct(Entrezid, .keep_all = TRUE)

# ----------------------------
# Subset Gene Groups
# ----------------------------

# Separate upregulated and downregulated genes
upregulated <- deg_results %>% filter(threshold == "Upregulated")
downregulated <- deg_results %>% filter(threshold == "Downregulated")
deg_updown <- bind_rows(upregulated, downregulated)


background_genes <- deg_results$Entrezid
sig_genes <- deg_updown$Entrezid

# -------------------------------------
# GO Overrepresentation Analysis (ORA)
# ------------------------------------
# Tests which GO terms (BP, MF, CC) are enriched in significant genes 
# compared to the background set

go_ora <- enrichGO(
  gene = as.character(sig_genes),
  universe = as.character(background_genes),
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL", # for specific ontology = "BP", "MF", or "CC"
  readable = TRUE
)

# View results

GO <- as.data.frame(go_ora)

# Visualize most enriched GO terms
dotplot(go_ora, showCategory = 20, title = "GO Overrepresentation (BP, MF, CC)")
barplot(go_ora, showCategory = 20, title = "GO Overrepresentation")

# Split visualization by ontology (BP, MF, CC)
barplot(go_ora,
        split = "ONTOLOGY",
        showCategory = 5,
        font.size = 10,
        color = "p.adjust",
        label_format = 60) +
  facet_grid(ONTOLOGY ~ ., scale = "free")

# ---------------------------------------
# KEGG Overrepresentation Analysis (ORA)
# ---------------------------------------
# KEGG identifies metabolic or signaling pathways significantly enriched

kegg_ora <- enrichKEGG(
  gene = as.character(sig_genes),
  organism = "hsa", # hsa = Homo sapiens
  keyType = "kegg",
  universe = as.character(background_genes)
)

kegg <- as.data.frame(kegg_ora)

# Visualize KEGG enrichment results
dotplot(kegg_ora, title = "KEGG Pathway Enrichment")

# ------------------------------------
# Gene Set Enrichment Analysis (GSEA)
# ------------------------------------
# GSEA identifies pathways enriched in *ranked* gene lists
# Useful when no strict threshold (like p < 0.05) is used

# Prepare ranked gene list: EntrezID mapped to logFC values
gene_ranks <- deg_results$logFC
names(gene_ranks) <- deg_results$Entrezid
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

# Retrieve Hallmark gene sets from MSigDB (category = H)
hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene)

# Run GSEA
gsea_res <- GSEA(
  geneList = gene_ranks, 
  TERM2GENE = hallmark,
  eps = 0
)

gsea <- as.data.frame(gsea_res)

# Visualization examples for specific pathways
gseaplot(gsea_res, geneSetID = "HALLMARK_BILE_ACID_METABOLISM")

# Compare enrichment profiles of two hallmark gene sets
p1 <- gseaplot(gsea_res, 
               geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
               by = "runningScore",
               title = "Epithelial-Mesenchymal Transition")

p2 <- gseaplot(gsea_res, 
               geneSetID = "HALLMARK_BILE_ACID_METABOLISM",
               by = "runningScore",
               title = "Bile Acid Metabolism")

cowplot::plot_grid(p1, p2, ncol = 1)

# Display top 3 enriched pathways with detailed plots
gseaplot2(gsea_res, geneSetID = 1:3)
gseaplot2(gsea_res, geneSetID = 1:3, pvalue_table = TRUE)

# ------------------------
# Summary and Reporting
# ------------------------

# Convert enrichment results to data frames for filtering and export
GO_df <- as.data.frame(go_ora)
KEGG_df <- as.data.frame(kegg_ora)
GSEA_df <- as.data.frame(gsea_res@result)

# Filter significant results (adjusted p < 0.05)
GO_sig <- GO_df %>% filter(p.adjust < 0.05)

# Separate GO results by ontology category
top_GO_BP <- GO_sig %>% filter(ONTOLOGY == "BP") %>%
  arrange(p.adjust) %>% head(10)
top_GO_MF <- GO_sig %>% filter(ONTOLOGY == "MF") %>%
  arrange(p.adjust) %>% head(10)
top_GO_CC <- GO_sig %>% filter(ONTOLOGY == "CC") %>%
  arrange(p.adjust) %>% head(10)

# Combine top GO terms across categories
top_GO <- bind_rows(
  top_GO_BP %>% mutate(Category = "Biological Process"),
  top_GO_MF %>% mutate(Category = "Molecular Function"),
  top_GO_CC %>% mutate(Category = "Cellular Component")
) %>% relocate(Category, .before = ID)

# KEGG and GSEA significant terms
KEGG_sig <- KEGG_df %>% filter(p.adjust < 0.05)
GSEA_sig <- GSEA_df %>% filter(p.adjust < 0.05)

# Separate GSEA results by direction (up/down)
GSEA_up <- GSEA_sig %>% filter(NES > 0) %>% arrange(desc(NES))
GSEA_down <- GSEA_sig %>% filter(NES < 0) %>% arrange(NES)

# ------------------------
# Create Summary Tables
# ------------------------

# Combine counts of significant pathways across analyses
summary_overall <- data.frame(
  Analysis = c("GO (ORA)", "KEGG (ORA)", "GSEA"),
  Total_Significant = c(nrow(GO_sig), nrow(KEGG_sig), nrow(GSEA_sig)),
  Upregulated = c(NA, NA, nrow(GSEA_up)),
  Downregulated = c(NA, NA, nrow(GSEA_down))
)

summary_overall  # View overview of enrichment findings

# ------------------------
# Export Results
# ------------------------

# Export enrichment tables
write.csv(GO_sig, "Results/GO_significant_terms.csv", row.names = FALSE)
write.csv(KEGG_sig, "Results/KEGG_significant_terms.csv", row.names = FALSE)
write.csv(GSEA_sig, "Results/GSEA_significant_terms.csv", row.names = FALSE)
write.csv(summary_overall, "Results/Summary_Enrichment_Report.csv", row.names = FALSE)

save.image("Module3_class1.RData")

