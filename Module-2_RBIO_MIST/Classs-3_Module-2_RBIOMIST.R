# ==========================================================================
# COURSE: R for Bioinformatics & Transcriptomic Data Analysis (MIST)
# MODULE 2: RNA-Seq Data Analysis
# TOPIC:    From Raw Count Data to Differential Gene Expression
# Instructors: RA Irfan Ahmed Rifat, RA Sanjana Chowdhury Arpa
# ==========================================================================
#   Pipeline: Package Install → GEO Download → Load & Inspect →
#             QC → Normalisation → DGE (DESeq2) → Volcano Plot
# ==========================================================================

# ── Microarray vs RNA-Seq: What changes in the pipeline? ─────────────────
#
#   Step              Microarray             RNA-Seq (this script)
#   ──────────────    ──────────────────     ──────────────────────────────
#   Raw data          .CEL files             Count matrix (integers)
#   Download          GEOquery / .CEL        GEOquery / getGEOSuppFiles()
#   QC                arrayQualityMetrics    DESeq2 size factors + plots
#   Normalisation     RMA (affy/oligo)       DESeq2 VST / rlog
#   DGE engine        limma (eBayes)         DESeq2 (negative binomial)
#   Input to DGE      log2 expression        RAW integer counts (NOT log2)
#   Output            topTable()             results() + lfcShrink()
#
# Key rule: DESeq2 MUST receive raw counts. Never give it log or normalised
#           data — it does its own normalisation internally.
# ─────────────────────────────────────────────────────────────────────────


# ==========================================================================
#   PART 1 │ PACKAGE INSTALLATION
# ==========================================================================

# ── CRAN packages ─────────────────────────────────────────────────────────
install.packages("tidyverse")   # dplyr, ggplot2, tidyr
install.packages("pheatmap")    # heatmaps
install.packages("ggrepel")     # non-overlapping labels on volcano plot

# ── Bioconductor packages ─────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("biomaRt", quietly = TRUE)) 
  BiocManager::install("biomaRt")




BiocManager::install("GEOquery")          # download datasets from NCBI GEO
BiocManager::install("DESeq2")            # RNA-Seq DGE (equivalent of limma)
BiocManager::install("EnhancedVolcano")   # publication volcano plots
BiocManager::install("apeglm")            # log fold change shrinkage for DESeq2

# ── Load all packages ─────────────────────────────────────────────────────
library(tidyverse)
library(GEOquery)
library(DESeq2)
library(pheatmap)
library(ggrepel)
library(EnhancedVolcano)
library(apeglm)
library(org.Hs.eg.db)
library(biomaRt)
library(dplyr)

# Confirm versions
sessionInfo()


# ==========================================================================
#   PART 2 │ DATASET DOWNLOAD FROM GEO
# ==========================================================================

# We use GSE173306 — a real published RNA-Seq study comparing
# colorectal cancer Tumor vs Normal tissue (Homo sapiens).
# It provides a pre-built count matrix as a supplementary file,
# which is the most common format you will find on GEO.

# ── Step 1: Download the GEO series metadata ──────────────────────────────
gse <- getGEO("GSE173306", GSEMatrix = TRUE, getGPL = FALSE)
metadata_raw <- pData(gse[[1]])

# Take a first look at what columns are available
head(metadata_raw)
colnames(metadata_raw)

# ── Step 2: Download the supplementary count matrix ───────────────────────
# GEO stores raw count files as supplementary files attached to the series.
getGEOSuppFiles("GSE173306")

# This creates a folder  GSE173306/  in your working directory.
# List what was downloaded:
list.files("GSE173306/")

# ── Step 3: Read the count matrix ─────────────────────────────────────────
# Adjust the filename to match what was actually downloaded above.
# It is usually a .txt.gz or .csv.gz file.

counts_raw <- read.delim(
  "GSE173306/GSE173306_UCSC_hg19_counts.txt.gz",
  row.names = 1,     # first column = gene IDs → use as row names
  check.names = FALSE
)

head(counts_raw)
dim(counts_raw)      # rows = genes, columns = samples


# ==========================================================================
#   PART 3 │ LOAD & INSPECT
# ==========================================================================

# ── 3A. Inspect the count matrix ──────────────────────────────────────────
dim(counts_raw)               # how many genes × samples?
head(counts_raw)              # first 6 genes
colnames(counts_raw)          # sample names
class(counts_raw[1,1])        # should be integer or numeric

# Check for any negative values (should never exist in raw counts)
any(counts_raw < 0)

# Check for NA values
sum(is.na(counts_raw))

# ── 3B. Build the sample metadata (colData) ───────────────────────────────
# DESeq2 needs a metadata table where rows = samples, columns = variables.
# The row names of colData MUST match the column names of counts_raw exactly.

# 1. Isolate the needed columns from the raw metadata
metadata_cleaned <- metadata_raw[, c("geo_accession", "title", "cell line:ch1", "treatment:ch1")]

# 2. Process, rename, and generate matching matrix keys
metadata <- metadata_cleaned %>%
  dplyr::rename(
    gsm_id      = geo_accession,
    sample_name = title,
    cell_line   = `cell line:ch1`,
    condition   = `treatment:ch1`
  ) %>%
  mutate(
    # Rebuild the exact names found in the count matrix columns
    matrix_key = case_when(
      condition == "untreated" ~ cell_line,
      condition == "crizotinib treated" ~ paste0(cell_line, "960"),
      TRUE ~ cell_line # Fallback baseline
    ),
    
    # Establish your factor levels for DESeq2 (untreated = reference control)
    condition = factor(condition, levels = c("untreated", "crizotinib treated"))
  )

# 3. CRITICAL STEP: Set row names to our newly engineered keys
rownames(metadata) <- metadata$matrix_key

# 4. Clean any whitespace mismatches automatically
colnames(counts_raw) <- trimws(colnames(counts_raw))
rownames(metadata)   <- trimws(rownames(metadata))

# 5. Reorder the counts matrix columns safely to match metadata rows exactly
counts_raw <- counts_raw[, rownames(metadata)]

# 6. Final verification test for DESeq2
all(colnames(counts_raw) == rownames(metadata))
head(metadata)
str(metadata)


# ==========================================================================
#   PART 4 │ QUALITY CONTROL
# ==========================================================================

# ── 4A. Library size (sequencing depth) per sample ────────────────────────
# Total counts per sample should be roughly similar.
# A sample with dramatically fewer counts may have failed during sequencing.

lib_sizes <- colSums(counts_raw)
lib_sizes
summary(lib_sizes)

# Bar plot of library sizes
barplot(lib_sizes / 1e6,
        names.arg = colnames(counts_raw),
        las = 2, cex.names = 0.7,
        col = ifelse(metadata$condition == "crizotinib treated", "#E74C3C", "#2980B9"),
        ylab = "Library Size (millions of reads)",
        main = "Sequencing Depth per Sample")
abline(h = 10, lty = 2, col = "grey50")   # 10M reads is a typical minimum

# ── 4B. Filter out lowly expressed genes ──────────────────────────────────
# Genes with near-zero counts across all samples carry no information and
# inflate the multiple testing correction burden.
#
# Rule of thumb: keep a gene only if it has at least 10 counts
# in at least (number of replicates in smallest group) samples.

n_replicates <- min(table(metadata$condition))   # smallest group size
keep <- rowSums(counts_raw >= 10) >= n_replicates

counts_filtered <- counts_raw[keep, ]

cat("Genes before filtering:", nrow(counts_raw), "\n")
cat("Genes after  filtering:", nrow(counts_filtered), "\n")
cat("Genes removed:         ", nrow(counts_raw) - nrow(counts_filtered), "\n")

# ── 4C. Build the DESeq2 object ───────────────────────────────────────────
# DESeqDataSetFromMatrix() is the entry point — equivalent to reading
# .CEL files into an ExpressionSet in the microarray workflow.

# ── CRITICAL STEP: Round the fractional counts to whole integers ──────────
counts_integer <- round(counts_filtered)

dds <- DESeqDataSetFromMatrix(
  countData = counts_integer,   # raw integer count matrix
  colData   = metadata,          # sample metadata
  design = ~ cell_line + condition       # the variable we want to test
)

dds

# ── 4D. VST — Variance Stabilising Transformation for QC plots ────────────
# Raw counts are very skewed (many zeros, a few very high values).
# VST log-transforms and stabilises variance — used ONLY for QC
# visualisation, NOT for the DGE test (DESeq2 uses raw counts for that).

vst_data <- vst(dds, blind = TRUE)   # blind = TRUE for QC (unbiased)

# PCA plot — do samples cluster by condition?
# If the QC is good: Tumor and Normal samples form two distinct clouds.
plotPCA(vst_data, intgroup = "condition") +
  scale_colour_manual(values = c("untreated" = "#2980B9", "crizotinib treated" = "#E74C3C")) +
  labs(title = "PCA — VST-transformed counts",
       subtitle = "Samples should cluster by condition") +
  theme_bw()

# Sample-to-sample distance heatmap — another QC view
sample_dists <- dist(t(assay(vst_data)))   # Euclidean distance between samples
sample_dist_matrix <- as.matrix(sample_dists)

pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  annotation_col = metadata[ , "condition", drop = FALSE],
  color  = colorRampPalette(c("#2980B9", "white"))(50),
  main   = "Sample-to-Sample Euclidean Distance",
  border_color = NA
)

# ── 4E. Gene count distribution boxplot ───────────────────────────────────
# Are the VST-normalised distributions similar across samples?
# Large differences suggest a normalisation or quality problem.

vst_matrix <- assay(vst_data)

boxplot(vst_matrix,
        las   = 2,
        col   = ifelse(metadata$condition == "crizotinib treated", "#E74C3C", "#2980B9"),
        main  = "VST Expression Distribution per Sample",
        ylab  = "VST Expression",
        cex.axis = 0.7)


# ==========================================================================
#   PART 5 │ NORMALISATION
# ==========================================================================

# DESeq2 normalises internally using the "median of ratios" method.
# This corrects for differences in sequencing depth (library size).
# You do NOT need a separate normalisation step — it happens inside DESeq().

# ── Estimate size factors (sequencing depth checking) ─────────────────────────
dds <- estimateSizeFactors(dds)
sizeFactors(dds)   # values close to 1.0 = samples were similar depth
#Note: Some students thought that we are normalizing here before giving it inside DSEq2.
#But This step simply calculates the library depth correction factors and stores them in a separate slot.




# Visualise: raw vs normalised counts for one gene
# 1. Grab your gene of interest
gene_to_check <- rownames(counts_filtered)[1]

# 2. Split the plotting window into 1 row, 2 columns
par(mfrow = c(1, 2))

# 3. Raw counts (Forced into a vector using unlist)
barplot(unlist(counts_filtered[gene_to_check, ]),
        las = 2, cex.names = 0.6,
        col = ifelse(metadata$condition == "crizotinib treated", "#E74C3C", "#2980B9"),
        main = paste("Raw counts:", gene_to_check),
        ylab = "Raw Count")

# 4. Normalised counts (We use as.numeric just to be completely safe)
norm_counts <- counts(dds, normalized = TRUE)
barplot(as.numeric(norm_counts[gene_to_check, ]),
        names.arg = colnames(norm_counts),  # Restores sample labels to the x-axis
        las = 2, cex.names = 0.6,
        col = ifelse(metadata$condition == "crizotinib treated", "#E74C3C", "#2980B9"),
        main = paste("Normalised counts:", gene_to_check),
        ylab = "Normalised Count")

# 5. Reset layout back to default single screen
par(mfrow = c(1, 1))

# ==========================================================================
#   PART 6 │ DIFFERENTIAL GENE EXPRESSION — DESeq2
# ==========================================================================

# ── 6A. Run the full DESeq2 pipeline ──────────────────────────────────────
# DESeq() performs three steps in one call:
#   1. estimateSizeFactors()   — normalise for sequencing depth
#   2. estimateDispersions()   — estimate gene-wise variance
#   3. nbinomWaldTest()        — statistical test (negative binomial model)
#
# Equivalent to eBayes() + topTable() in limma.

dds <- DESeq(dds)

# Check what comparisons are available
resultsNames(dds)   # "condition_Tumor_vs_Normal" is what we want

# ── 6B. Extract raw results ───────────────────────────────────────────────
res_raw <- results(
  dds,
  contrast = c("condition", "crizotinib treated", "untreated"),  
  alpha    = 0.05   # FDR threshold for the independent filtering step
)

summary(res_raw)
head(res_raw)

# ── 6C. Log Fold Change Shrinkage (apeglm) ────────────────────────────────
# Genes with very low counts have noisy, inflated fold changes.
# lfcShrink() pulls those extreme estimates back towards zero,
# giving a more reliable ranking — especially important for volcano plots.
# This is the DESeq2 equivalent of limma's empirical Bayes smoothing.



res_shrunk <- lfcShrink(
  dds,
  coef = "condition_crizotinib.treated_vs_untreated",  # Updated to match your exact groups!
  type = "apeglm"
)

summary(res_shrunk)

# ── 6D. Convert results to a tidy data frame ──────────────────────────────
deg_results <- as.data.frame(res_shrunk) %>%
  rownames_to_column("Gene") %>%          # move gene IDs from rownames to a column
  arrange(padj) %>%                        # sort by adjusted p-value
  filter(!is.na(padj))                     # remove genes that failed filtering

head(deg_results)
dim(deg_results)

# ── 6E. Classify each gene ────────────────────────────────────────────────
# Standard thresholds:
#   |log2FoldChange| > 1   (2-fold change minimum)
#   padj < 0.05            (5% false discovery rate)

deg_results <- deg_results %>%
  mutate(status = case_when(
    log2FoldChange >  1 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -1 & padj < 0.05 ~ "Downregulated",
    TRUE                               ~ "Not Significant"
  ))

# Summary counts
table(deg_results$status)

# Top 10 most significant genes
deg_results %>%
  filter(status != "Not Significant") %>%
  arrange(padj) %>%
  dplyr::select(Gene, log2FoldChange, padj, status) %>%
  head(10)

# Save full results table
# Create the results directory first
dir.create("results", showWarnings = FALSE)

# Now save your table — this will run perfectly!
write.csv(deg_results, "results/DESeq2_crizotinib_treated_vs_untreated_results.csv", row.names = FALSE)




# Annotation from UCSC
library(biomaRt)

# ── 1. Connect to an alternative stable Ensembl Mirror ─────────────────────
# We use the US East or Asia mirror since the main site is currently down.
mart <- useEnsembl(
  biomart = "genes", 
  dataset = "hsapiens_gene_ensembl", 
  mirror  = "useast" # Options: 'useast', 'uswest', 'asia'
)

# ── 2. Query the mirror server using your clean_id column ──────────────────
annotations <- getBM(
  attributes = c("ucsc", "hgnc_symbol"), 
  filters    = "ucsc", 
  values     = deg_results$clean_id, 
  mart       = mart
)

# ── 3. Clean up the query response and fix data types ──────────────────────
annotations_clean <- annotations %>% 
  filter(hgnc_symbol != "") %>%          
  distinct(ucsc, .keep_all = TRUE) %>%
  mutate(ucsc = as.character(ucsc)) # Keep our previous character fix!

# ── 4. Merge the clean symbols back into your main DEG table ───────────────
deg_results <- deg_results %>%
  left_join(annotations_clean, by = c("clean_id" = "ucsc")) %>%
  mutate(
    # Fallback to the original Gene transcript ID if no symbol was found
    gene_symbol = ifelse(is.na(hgnc_symbol) | hgnc_symbol == "", Gene, hgnc_symbol)
  )

# ── 5. Preview your success! ───────────────────────────────────────────────
head(deg_results[, c("Gene", "gene_symbol", "log2FoldChange", "padj", "status")])







# ==========================================================================
#   PART 7 │ VISUALISATION
# ==========================================================================

# ── PLOT 1: MA Plot — built-in DESeq2 QC plot ────────────────────────────
# X axis = mean expression (A), Y axis = log2FC (M).
# Points shrunk towards the centre by lfcShrink() — this is good.
# Blue dots = significant genes.

plotMA(res_shrunk,
       ylim  = c(-5, 5),
       main  = "MA Plot — Crizotinib Resistant vs Untreated (apeglm shrinkage)",
       alpha = 0.05)

# ── PLOT 2: Volcano Plot (manual ggplot2) ─────────────────────────────────
# X = log2FC (effect size), Y = -log10(padj) (significance).
# Top-right = strongly upregulated in Tumor
# Top-left  = strongly downregulated in Tumor

# Label the top 15 most significant genes
top_genes <- deg_results %>%
  filter(status != "Not Significant") %>%
  arrange(padj) %>%
  head(15)

ggplot(deg_results, aes(x = log2FoldChange,
                        y = -log10(padj),
                        colour = status)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(data  = top_genes,
                  aes(label = Gene),
                  size  = 3,
                  max.overlaps = 20) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c(
    "Upregulated"    = "#E74C3C",
    "Downregulated"  = "#2980B9",
    "Not Significant"= "#BDC3C7"
  )) +
  labs(
    title    = "Volcano Plot — Tumor vs Normal (DESeq2)",
    subtitle = paste0("Up: ", sum(deg_results$status == "Upregulated"),
                      "  |  Down: ", sum(deg_results$status == "Downregulated")),
    x        = "log2 Fold Change",
    y        = "-log10 (adjusted p-value)",
    colour   = "Status"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("results/Volcano_DESeq2_Tumor_vs_Normal.jpeg",
       width = 8, height = 6, dpi = 300)

# ── PLOT 3: EnhancedVolcano (one-line publication version) ────────────────
EnhancedVolcano(
  deg_results,
  lab      = deg_results$Gene,
  x        = "log2FoldChange",
  y        = "padj",
  pCutoff  = 0.05,
  FCcutoff = 1,
  title    = "treated vs crizotinib treated",
  subtitle = "DESeq2 | apeglm shrinkage",
  col      = c("#BDC3C7", "#2980B9", "#F39C12", "#E74C3C")
)

# ── PLOT 4: Heatmap of top 30 significant DEGs ────────────────────────────
top30_genes <- deg_results %>%
  filter(status != "Not Significant") %>%
  arrange(padj) %>%
  head(30) %>%
  pull(Gene)

# Extract VST values for these genes
heatmap_matrix <- assay(vst_data)[top30_genes, ]

# Row-scale so each gene is comparable (z-score)
heatmap_scaled <- t(scale(t(heatmap_matrix)))

# Annotation for columns (samples)
annotation_col <- metadata %>%
  dplyr::select(condition) %>%
  as.data.frame()
rownames(annotation_col) <- rownames(metadata)

ann_colours <- list(
  condition = c(
    "untreated"          = "#2980B9", 
    "crizotinib treated" = "#E74C3C"
  )
)

# ── Run the Heatmap ────────────────────────────────────────────────────────
pheatmap(
  heatmap_scaled,
  annotation_col    = annotation_col,
  annotation_colors = ann_colours,
  color             = colorRampPalette(c("#2980B9", "white", "#E74C3C"))(50),
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_colnames     = FALSE,
  fontsize_row      = 8,
  main              = "Top 30 DEGs — Crizotinib Resistant vs Untreated (z-scored VST)",
  border_color      = NA
)


dev.off()


# ==========================================================================
#   PIPELINE RECAP
# ==========================================================================
#
#   Step   Function / Package    What it does
#   ─────  ────────────────────  ───────────────────────────────────────────
#   1      getGEO()              Download metadata from GEO
#          getGEOSuppFiles()     Download raw count matrix
#   2      read.delim()          Load count matrix into R
#   3      pData()               Extract sample metadata
#   4      rowSums() filter      Remove lowly expressed genes
#          vst()                 VST for QC visualisation only
#          plotPCA()             Check sample clustering
#          pheatmap()            Sample distance heatmap
#   5      DESeq()               Normalise + fit model (all in one)
#          sizeFactors()         Inspect normalisation factors
#   6      results()             Raw DGE results
#          lfcShrink(apeglm)     Stabilise fold changes for low-count genes
#          case_when()           Classify Up / Down / NS
#   7      ggplot2 + ggrepel     Volcano plot
#          EnhancedVolcano       One-line publication volcano
#          pheatmap              Top 30 DEG heatmap
#
#   NEXT:  Functional Enrichment (GO / KEGG) with clusterProfiler
# ==========================================================================
