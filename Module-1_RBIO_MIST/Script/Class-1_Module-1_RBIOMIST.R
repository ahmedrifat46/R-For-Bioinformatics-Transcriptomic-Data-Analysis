# ==========================================================================
# COURSE: R For Bioinformatics & Transcriptomic Data Analysis (MIST)
# MODULE 1: The Grammar of R, Data Wrangling and Publication Grade Visualizations
# CLASS 1: Environment Setup & The Grammar of R
# ==========================================================================
 
 
# 1. ENVIRONMENT SETUP -----------------------------------------------------
# Goal: Checking R version and installing the Bioconductor manager.
 
sessionInfo() # shows R version, OS, and loaded packages
 
 
# 2. R AS A CALCULATOR (OPERATORS) -----------------------------------------
# Goal: Practice Arithmetic and Logical operators.
 
# --- Arithmetic ---
2 + 3 # addition
10 - 4 # subtraction
6 * 7 # multiplication
100 / 4 # division
2 ^ 8 # exponent (2 to the power of 8)
17 %% 5 # remainder (17 divided by 5 leaves 2)
 
# Bioinformatics context: convert a raw count to log2 scale
log2(1000) # log2-transformation is used in every RNA-Seq pipeline
 
# --- Comparison operators → always return TRUE or FALSE ---
10 > 5 # TRUE
3 == 3 # TRUE (double == means "is equal to?")
3 != 4 # TRUE (is NOT equal to?)
7 >= 10 # FALSE
 
# --- Logical operators → combine conditions ---
(5 > 3) & (2 < 4) # AND — both must be TRUE
(5 > 3) | (2 > 4) # OR — at least one must be TRUE
!(5 > 3) # NOT — flips the result
 
 
# 3. VARIABLES AND ASSIGNMENT ----------------------------------------------
# Goal: Storing data in named containers.
 
# Syntax: variable_name <- value
 
patient_age <- 45 # numeric
gene_name <- "BRCA1" # character (text always in quotes)
expression_value <- 6.83 # numeric (log2 expression)
is_tumor <- TRUE # logical
 
# Print them to the Console
patient_age
gene_name
expression_value
is_tumor
 
# Variables can be overwritten at any time
patient_age <- 52
patient_age
 
# R is CASE-SENSITIVE
Gene_Name <- "TP53"
gene_name # still "BRCA1"
Gene_Name # "TP53" — a completely different variable
 
 
# 4. BASIC DATA TYPES ------------------------------------------------------
# Goal: Understanding Numeric, Character, and Logical data.
 
# Use class() to ask R "what type is this?"
 
class(patient_age) # "numeric"
class(gene_name) # "character"
class(is_tumor) # "logical"
class(42L) # "integer" — the L suffix forces integer type
class(TRUE) # "logical"

#Also you may use typeof() function to check the data type
 
# Common trap: numbers inside quotes become TEXT
count_text <- "1050"
class(count_text) # "character" — looks like a number but is NOT
# mean(count_text) # uncomment to see the error this causes
 
count_num <- 1050
class(count_num) # "numeric" — this one works correctly
mean(count_num)
 
# Type conversion functions
as.numeric("3.14") # character → numeric
as.integer(6.99) # numeric → integer (truncates, does NOT round → 6)
as.character(100) # numeric → character
as.logical(0) # 0 = FALSE; anything else = TRUE
 
 
# 5. DATA STRUCTURES: VECTORS ----------------------------------------------
# Goal: Creating lists of genes or expression values.
 
# c() = "combine" — the main way to build a vector
gene_ids <- c("TP53", "BRCA1", "EGFR", "MYC", "PTEN")
expr_values <- c(2.30, 5.61, 4.80, 7.12, 1.95)
is_sig <- c(FALSE, TRUE, TRUE, TRUE, FALSE)
 
gene_ids
expr_values
 
# All values in a vector must be the same type
# If you mix types, R silently converts everything to the most flexible type
mixed <- c("TP53", 5.61, TRUE)
mixed # everything becomes a character — note the quotes
 
# --- Indexing: extract elements with [ ] ---
gene_ids[1] # first gene (R starts counting at 1, not 0)
gene_ids[3] # third gene
expr_values[2:4] # elements 2 through 4 (a range)
gene_ids[c(1, 4)] # elements 1 and 4 (specific positions)
 
# --- Useful vector functions ---
length(gene_ids) # how many elements?
sum(expr_values) # total
mean(expr_values) # average expression
max(expr_values) # highest value
min(expr_values) # lowest value
sort(expr_values) # ascending order
 
 
# 6. DATA STRUCTURES: DATA FRAMES ------------------------------------------
# Goal: The "Bio-Table" — combining clinical and genomic data.
 
# A data frame is a TABLE where each column can be a different type.
# Rows = observations (patients, genes)
# Columns = variables (age, expression, stage …)
 
expression_data <- data.frame(
gene = c("TP53", "BRCA1", "EGFR", "MYC", "PTEN"),
normal = c(2.10, 7.30, 4.50, 1.20, 6.00), # log2 expression — Normal
tumor = c(5.80, 3.20, 4.40, 8.10, 2.10), # log2 expression — Tumor
is_cancer_gene = c(TRUE, TRUE, FALSE, TRUE, TRUE),
stringsAsFactors = FALSE
)
 
expression_data # view the whole table
 
# Quick inspection functions — use these every time you load new data
str(expression_data) # structure: column types and a data preview
head(expression_data) # first 6 rows
dim(expression_data) # dimensions: rows × columns
names(expression_data) # column names only
 
 
# 7. ACCESSING DATA (INDEXING) ---------------------------------------------
# Goal: Mastering the [Row, Column] and $ notation.
 
# --- $ notation: extract a single column as a vector ---
expression_data$gene
expression_data$tumor
 
# --- [ row , column ] notation ---
expression_data[1, ] # entire first row (TP53)
expression_data[ , 2] # entire second column (normal values)
expression_data[2, 3] # row 2, column 3 → BRCA1 tumor value
expression_data[1:3, ] # first three rows
expression_data[c(1,4), c("gene","tumor")] # rows 1 & 4, specific columns
 
# --- Conditional selection: filter rows by a condition ---
expression_data[expression_data$tumor > 5, ] # genes highly expressed in tumor
expression_data[expression_data$is_cancer_gene, ] # known cancer genes only
 
# --- Add a new calculated column ---
expression_data$log2FC <- expression_data$tumor - expression_data$normal
# (subtracting log2 values = dividing the raw counts → this IS log2 Fold Change)
 
expression_data
  
# --- The Help system ---
help(mean) # open full documentation for mean()
?log2 # shorthand — same as help(log2)


# ==========================================================================

#---------------In the next class we will learn-----------------------------

# ==========================================================================
 

#functions
#saving workspace > this is one of the most important thing to know
#Data Wrangling  > Very important to know to handle and analyze big data

#-------------Request------------------
#Everyone don't be late tomorrow
