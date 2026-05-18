# ==========================================================================
# COURSE: R for Bioinformatics & Transcriptomic Data Analysis (MIST)
# MODULE 1: The Foundations of Computational Oncology
# CLASS 3: Packages, Tidyverse, Joins, Reshaping & ggplot2 Visualisation
# ==========================================================================
#   Instructors: RA Irfan Ahmed Rifat | RA Sanjana Chowdhury Arpa
#   Department of Biomedical Engineering, MIST
# ==========================================================================

# ── What we cover today ───────────────────────────────────────────────────
#   PART 1  │  Quick Restart   – one master dataset to use all class
#   PART 2  │  R Packages      – install, load, the Tidyverse family
#   PART 3  │  The Pipe %>%    – select, filter, mutate, arrange
#   PART 4  │  Relational Join – left_join() to merge expression + clinical
#   PART 5  │  Reshaping       – pivot_longer() for the "tidy" format
#   PART 6  │  ggplot2         – histogram, boxplot, violin, scatter, heatmap
#   PART 7  │  gtsummary       – publication Table 1 in one line
#   PART 8  │  Saving Plots    – ggsave(), PDF, high-resolution JPEG
# ─────────────────────────────────────────────────────────────────────────


# ==========================================================================
#   PART 1 │ QUICK RESTART — One Master Dataset for the Whole Class
# ==========================================================================
# We build TWO tables once and reuse them everywhere.
# Table A  →  expression_wide   (genes × patients)
# Table B  →  clinical          (patients × clinical variables)

# ── Table A: Gene expression (wide format — one column per gene) ──────────
expression_wide <- data.frame(
  patient_id  = c("PT_01","PT_02","PT_03","PT_04","PT_05",
                  "PT_06","PT_07","PT_08","PT_09","PT_10"),
  TP53   = c(5.2, 1.1, 6.3, 2.0, 7.1, 1.8, 5.9, 2.4, 6.8, 1.5),
  BRCA1  = c(3.4, 7.8, 2.9, 8.1, 3.1, 7.5, 2.7, 8.4, 3.0, 7.2),
  EGFR   = c(4.1, 4.0, 6.5, 2.2, 5.8, 3.9, 6.0, 2.5, 5.5, 4.3),
  MYC    = c(7.8, 2.3, 8.2, 1.9, 8.5, 2.1, 7.6, 1.7, 8.0, 2.4),
  stringsAsFactors = FALSE
)

# ── Table B: Clinical metadata ────────────────────────────────────────────
clinical <- data.frame(
  patient_id    = c("PT_01","PT_02","PT_03","PT_04","PT_05",
                    "PT_06","PT_07","PT_08","PT_09","PT_10"),
  condition     = c("Tumor","Normal","Tumor","Normal","Tumor",
                    "Normal","Tumor","Normal","Tumor","Normal"),
  disease_stage = c("II","NA","III","NA","IV","NA","II","NA","III","NA"),
  age           = c(45, 60, 52, 38, 67, 55, 49, 42, 63, 57),
  gender        = c("F","M","F","M","F","M","F","M","F","M"),
  smoking       = c("Yes","No","Yes","No","Yes","No","No","No","Yes","No"),
  stringsAsFactors = FALSE
)

# Quick look
head(expression_wide)
str(clinical)


# ==========================================================================
#   PART 2 │ R PACKAGES — Installing, Loading, the Tidyverse Family
# ==========================================================================

# Packages = "apps" or "toolboxes" that extend R's capability.
# Two steps, every session:
#   Step 1 → install.packages()  — download ONCE (like installing an app)
#   Step 2 → library()           — load into session EVERY time

# Install (run once — comment out afterwards)
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("gtsummary")

# Load for this session
library(tidyverse)   # loads dplyr, ggplot2, tidyr, readr, stringr, tibble
library(pheatmap)
library(gtsummary)

# ── Key packages at a glance ─────────────────────────────────────────────
#
#   Package        Purpose
#   ─────────────  ─────────────────────────────────────────────────────
#   dplyr          Data wrangling  (filter, select, mutate …)
#   ggplot2        Publication-quality visualisation
#   tidyr          Reshaping  (pivot_longer / pivot_wider)
#   readr          Fast CSV import
#   stringr        Text / string manipulation
#   tibble         Modern data frames
#   pheatmap       Heatmaps
#   gtsummary      Publication-ready summary tables
#   BiocManager    Gateway to Bioconductor biology packages
#   limma          Microarray differential expression  (Module 2)
#   DESeq2         RNA-Seq differential expression     (Module 2)

#We will Install Bioconductor packages once we start the Module 2 (run once)
# install.packages("BiocManager")
# BiocManager::install("limma")
# BiocManager::install("DESeq2")

sessionInfo()   # confirm tidyverse is attached


# ==========================================================================
#   PART 3 │ THE PIPE  %>%  AND CORE dplyr VERBS
# ==========================================================================

# The pipe  %>%  passes the result on the left into the next function.
# Read it as "and then".

# ── select() — keep or drop columns ──────────────────────────────────────
clinical %>%
  select(patient_id, condition, age)

clinical %>%
  select(-smoking)                     # drop one column with  -

# ── filter() — keep rows that match a condition ───────────────────────────
clinical %>%
  filter(condition == "Tumor")

clinical %>%
  filter(age > 50, smoking == "Yes")   # AND — comma works like &

# ── mutate() — add or transform columns ──────────────────────────────────
clinical %>%
  mutate(age_group = ifelse(age >= 55, "Older", "Younger"))

# ── arrange() — sort rows ─────────────────────────────────────────────────
clinical %>%
  arrange(desc(age))

# ── group_by() + summarise() — group statistics ───────────────────────────
clinical %>%
  group_by(condition) %>%
  summarise(
    n          = n(),
    mean_age   = round(mean(age), 1),
    n_smokers  = sum(smoking == "Yes")
  )

# ── Chaining multiple verbs — reads like a sentence ───────────────────────
# "From clinical data, keep Tumor patients over 50,
#  add an age group label, then sort by age."

clinical %>%
  filter(condition == "Tumor") %>%
  filter(age > 50) %>%
  mutate(age_group = ifelse(age >= 60, "Senior", "Adult")) %>%
  select(patient_id, age, age_group, disease_stage) %>%
  arrange(desc(age))

# ── Without pipe — the same code is almost unreadable ─────────────────────
arrange(
  select(
    mutate(
      filter(filter(clinical, condition == "Tumor"), age > 50),
      age_group = ifelse(age >= 60, "Senior", "Adult")
    ),
    patient_id, age, age_group, disease_stage
  ),
  desc(age)
)


# ==========================================================================
#   PART 4 │ RELATIONAL JOIN — Merging Expression + Clinical Tables
# ==========================================================================

# In real bioinformatics we always have two separate files:
#   1. Expression matrix  (genes in columns, patients in rows)
#   2. Clinical metadata  (patient demographics, stages, outcomes)
#
# left_join() merges them by a shared column (here: patient_id).
# Every row in Table A is kept; matching rows from Table B are added.

combined <- left_join(expression_wide, clinical, by = "patient_id")

combined
str(combined)
dim(combined)   # 10 patients × 10 columns — one unified analysis table


# ==========================================================================
#   PART 5 │ RESHAPING — pivot_longer() — Wide → Tidy (Long) Format
# ==========================================================================

# Wide format  →  one column per gene  →  good for storage
# Long  format  →  one row per gene-patient measurement  →  required by ggplot2

# ── Visualising what "wide" looks like right now ──────────────────────────
head(expression_wide)   # TP53, BRCA1, EGFR, MYC are side-by-side columns

# ── pivot_longer(): rotate gene columns into Gene + Expression rows ────────
long_data <- combined %>%
  pivot_longer(
    cols      = c(TP53, BRCA1, EGFR, MYC),   # which columns to rotate
    names_to  = "Gene",                        # new column: gene name
    values_to = "Expression"                   # new column: expression value
  )

head(long_data, 12)   # notice: 4 rows per patient now  (10 × 4 = 40 rows)
dim(long_data)


# ==========================================================================
#   PART 6 │ DATA VISUALISATION WITH ggplot2
# ==========================================================================

# ── The Grammar of Graphics — three required layers ───────────────────────
#
#   ggplot(data, aes(x, y, colour, fill, …))   ← data + aesthetics
#     + geom_*()                                ← geometry  (the plot type)
#     + labs()                                  ← labels
#     + theme_*()                               ← overall look

# We use long_data and combined for all plots below.

# ─────────────────────────────────────────────────────────────────────────
# PLOT 1 │ Histogram — QC check on raw expression distribution
# ─────────────────────────────────────────────────────────────────────────
# Clinical question: Are my expression values normally distributed?
# If counts are skewed, we need log2 transformation before analysis.

ggplot(long_data, aes(x = Expression)) +
  geom_histogram(bins = 15, fill = "#2980B9", colour = "white") +
  facet_wrap(~ Gene, scales = "free_x") +   # one panel per gene
  labs(
    title    = "Distribution of Gene Expression Values",
    subtitle = "QC check before differential expression analysis",
    x        = "log2 Expression",
    y        = "Count"
  ) +
  theme_bw()

# ─────────────────────────────────────────────────────────────────────────
# PLOT 2 │ Boxplot — Compare expression between Tumor and Normal
# ─────────────────────────────────────────────────────────────────────────
# Clinical question: Which genes show the clearest Tumor vs Normal contrast?

ggplot(long_data, aes(x = condition, y = Expression, fill = condition)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, alpha = 0.8) +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#2980B9", "Tumor" = "#E74C3C")) +
  labs(
    title = "Gene Expression: Tumor vs Normal",
    x     = "Condition",
    y     = "log2 Expression",
    fill  = "Condition"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# ─────────────────────────────────────────────────────────────────────────
# PLOT 3 │ Violin + Jitter — Show distribution shape AND individual points
# ─────────────────────────────────────────────────────────────────────────
# Hybrid layering: violin (distribution shape) + jitter (raw data points)

ggplot(long_data, aes(x = condition, y = Expression, fill = condition)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8, colour = "black") +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#2980B9", "Tumor" = "#E74C3C")) +
  labs(
    title    = "Expression Distribution with Individual Data Points",
    subtitle = "Violin = density shape  |  Dots = individual patients",
    x        = "Condition",
    y        = "log2 Expression",
    fill     = "Condition"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

# ─────────────────────────────────────────────────────────────────────────
# PLOT 4 │ Scatter Plot + Trend Line — Co-expression analysis
# ─────────────────────────────────────────────────────────────────────────
# Clinical question: Is there a co-expression relationship between
# TP53 and MYC? (Anti-correlation is a known oncology finding)

ggplot(combined, aes(x = TP53, y = MYC, colour = condition)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed",
              colour = "grey40") +
  scale_colour_manual(values = c("Normal" = "#2980B9", "Tumor" = "#E74C3C")) +
  labs(
    title    = "TP53 vs MYC Co-expression",
    subtitle = "Negative correlation suggests tumour suppressor / oncogene relationship",
    x        = "TP53 (log2 Expression)",
    y        = "MYC  (log2 Expression)",
    colour   = "Condition"
  ) +
  theme_bw()

# ── Extra: annotate with Pearson R ────────────────────────────────────────
r_val <- round(cor(combined$TP53, combined$MYC), 3)

ggplot(combined, aes(x = TP53, y = MYC, colour = condition)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed",
              colour = "grey40") +
  annotate("text", x = 5.5, y = 8.2,
           label = paste0("Pearson r = ", r_val),
           size = 4.5, fontface = "italic", colour = "grey20") +
  scale_colour_manual(values = c("Normal" = "#2980B9", "Tumor" = "#E74C3C")) +
  labs(
    title  = "TP53 vs MYC Co-expression  (annotated)",
    x      = "TP53 (log2 Expression)",
    y      = "MYC  (log2 Expression)",
    colour = "Condition"
  ) +
  theme_bw()

# ─────────────────────────────────────────────────────────────────────────
# PLOT 5 │ Bar Plot — Mean expression per gene per condition
# ─────────────────────────────────────────────────────────────────────────

mean_expr <- long_data %>%
  group_by(Gene, condition) %>%
  summarise(mean_expr = round(mean(Expression), 3), .groups = "drop")

ggplot(mean_expr, aes(x = Gene, y = mean_expr, fill = condition)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("Normal" = "#2980B9", "Tumor" = "#E74C3C")) +
  labs(
    title = "Mean Gene Expression by Condition",
    x     = "Gene",
    y     = "Mean log2 Expression",
    fill  = "Condition"
  ) +
  theme_bw() +
  theme(legend.position = "top")

# ─────────────────────────────────────────────────────────────────────────
# PLOT 6 │ Heatmap — Visualise expression patterns across all patients
# ─────────────────────────────────────────────────────────────────────────
# pheatmap needs a pure numeric matrix with row names.

# Prepare matrix: patients as rows, genes as columns
expr_matrix <- as.matrix(expression_wide[ , -1])   # drop patient_id column
rownames(expr_matrix) <- expression_wide$patient_id

# Annotation sidebar: colour rows by condition
annotation_rows <- data.frame(
  Condition = clinical$condition,
  row.names = clinical$patient_id
)

# Colour palette for annotation
ann_colours <- list(
  Condition = c(Normal = "#2980B9", Tumor = "#E74C3C")
)

# Basic heatmap — scale = "column" makes each gene comparable
pheatmap(
  expr_matrix,
  scale            = "column",        # z-score per gene
  annotation_row   = annotation_rows,
  annotation_colors = ann_colours,
  color            = colorRampPalette(c("#2980B9", "white", "#E74C3C"))(50),
  main             = "Gene Expression Heatmap (z-scored)",
  fontsize_row     = 9,
  fontsize_col     = 10,
  border_color     = NA,
  cluster_rows     = TRUE,            # hierarchical clustering of patients
  cluster_cols     = TRUE             # clustering of genes
)

# ── Why viridis instead of red/green? ─────────────────────────────────────
# ~8% of men have red-green colour blindness.
# Viridis is perceptually uniform AND readable in greyscale.

pheatmap(
  expr_matrix,
  scale            = "column",
  annotation_row   = annotation_rows,
  annotation_colors = ann_colours,
  color            = viridis::viridis(50),
  main             = "Heatmap with Viridis palette (journal-recommended)",
  fontsize_row     = 9,
  border_color     = NA
)


# ==========================================================================
#   PART 7 │ gtsummary — Publication-Ready Table 1 in One Line
# ==========================================================================

# Table 1 = the patient demographics table found in every clinical paper.
# Manually building it in Excel is error-prone and time-consuming.
# gtsummary does it automatically — and updates if the data changes.

# ── Basic Table 1 ─────────────────────────────────────────────────────────
clinical %>%
  select(age, gender, smoking, disease_stage) %>%
  tbl_summary()

# ── Stratified by condition (Tumor vs Normal) ─────────────────────────────
# This is the standard format for a clinical paper — one column per group.

clinical %>%
  select(condition, age, gender, smoking, disease_stage) %>%
  tbl_summary(by = condition) %>%
  add_p() %>%                      # adds statistical test p-values
  add_overall() %>%                # adds an "Overall" column
  bold_labels()                    # bold row labels for readability

# ── Add mean expression values to the table ───────────────────────────────
clinical_with_expr <- combined %>%
  select(condition, age, gender, smoking, TP53, MYC)

clinical_with_expr %>%
  tbl_summary(
    by = condition,
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",    # mean ± SD for numeric cols
      all_categorical() ~ "{n} ({p}%)"        # n (%) for categorical cols
    ),
    digits = all_continuous() ~ 2
  ) %>%
  add_p() %>%
  bold_labels()


# ==========================================================================
#   PART 8 │ SAVING PLOTS
# ==========================================================================

# ── Save the last printed plot with ggsave() ──────────────────────────────
# ggsave() automatically picks up the most recently drawn ggplot.

# High-resolution JPEG (for IEEE / EMBS submissions — 300 dpi minimum)
ggsave(
  filename = "results/TP53_vs_MYC_scatter.jpeg",
  width    = 7,
  height   = 5,
  dpi      = 300,
  units    = "in"
)

# PDF vector format (infinitely scalable — best for journals)
ggsave(
  filename = "results/TP53_vs_MYC_scatter.pdf",
  width    = 7,
  height   = 5
)

# ── Save a SPECIFIC plot by storing it in a variable first ────────────────
violin_plot <- ggplot(long_data,
                      aes(x = condition, y = Expression, fill = condition)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.8, colour = "black") +
  facet_wrap(~ Gene, scales = "free_y") +
  scale_fill_manual(values = c("Normal" = "#2980B9", "Tumor" = "#E74C3C")) +
  labs(title = "Expression Distribution — Violin + Jitter",
       x = "Condition", y = "log2 Expression", fill = "Condition") +
  theme_bw() +
  theme(legend.position = "bottom")

violin_plot   # display it

ggsave("results/violin_expression.jpeg", plot = violin_plot,
       width = 9, height = 6, dpi = 300)

# ── Save the pheatmap ────────────────────────────────────────────────────
# pheatmap is not a ggplot object — use pdf() / dev.off() instead

pdf("results/expression_heatmap.pdf", width = 7, height = 6)
pheatmap(
  expr_matrix,
  scale            = "column",
  annotation_row   = annotation_rows,
  annotation_colors = ann_colours,
  color            = viridis::viridis(50),
  main             = "Gene Expression Heatmap",
  border_color     = NA
)
dev.off() # IMPORTANT: always close the “image showing” device using dev.off() after an image generation before going to generate a new one — finalises the file
