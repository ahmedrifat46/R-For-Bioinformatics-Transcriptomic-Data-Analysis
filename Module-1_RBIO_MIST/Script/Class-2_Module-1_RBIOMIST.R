# ============================================================
#### R for Bioinformatics & Transcriptomic Data Analysis ####
#### Module 1: Class 2 - Data Wrangling in R ####
#### Dept. of BME, MIST ####
# ============================================================

# creating new directory
dir.create("Results")
dir.create("Scripts")
dir.create("Datasets")

# From last class
gene_ids    <- c("TP53", "BRCA1", "EGFR", "MYC", "PTEN")
expr_normal <- c(2.10,   7.30,   4.50,  1.20,  6.00)
expr_tumor  <- c(5.80,   3.20,   4.40,  8.10,  2.10)

expression_data <- data.frame(
  gene   = gene_ids,
  normal = expr_normal,
  tumor  = expr_tumor,
  stringsAsFactors = FALSE
)

expression_data

# ---------------------------------------------------------
#### Functions in R ####
# ---------------------------------------------------------
# A function is a reusable block of code.

# in-built function
sum(expr_normal)          # total expression
mean(expr_tumor)          # average
max(expr_tumor)           # highest value
min(expr_normal)          # lowest value
log2(1000)                # log2 transformation — critical for RNA-Seq
sqrt(81)                  # square root
round(3.14159, digits = 2)
nchar("BRCA1")            # number of characters in a string


#User-defined Function
# A function in R has 4 key parts:
#   1. Name         -> the name you give to the function
#   2. Arguments    -> the inputs you provide to the function
#   3. Body         -> the set of operations the function performs
#   4. Return Value -> the output the function gives back

# Without Argument
course_info <- function(){
  print("Welcome to R for Bioinformatics")
}

# calling a function
course_info()

#User-defined function (with arguments)

# Example function 1: calculate BMI

calculate_BMI <- function(weight, height) {
  bmi <- weight / (height ^ 2)
  return(bmi)
}

# Call the function
calculate_BMI(60, 1.75)

# Call the function using named arguments
calculate_BMI(weight = 60, height = 1.75)

# Intentional error:
calculate_BMI(60) # This gives an error because height is missing

# Function with a default value
calculate_BMI2 <- function(weight, height = 1.75) {
  bmi <- weight / (height ^ 2)
  return(bmi)
}

calculate_BMI2(60)


calculate_BMI3 <- function(weight, height, age) {
  bmi <- weight / (height ^ 2)
  return(bmi)
}

calculate_BMI3(60, 1.65)

height <- c(1.85, 1.75, 1.5)
weight <- c(50, 75, 45)

calculate_BMI(weight = weight, height = height)
calculate_BMI(height = height, weight = weight)
calculate_BMI(height, weight)

# Example of user-defined function: log2FC
#   Context: in every RNA-Seq / microarray study we compare
#   gene expression in TUMOR vs NORMAL.
#   log2FC > 0  → gene is MORE expressed in tumor  (Upregulated)
#   log2FC < 0  → gene is LESS expressed in tumor  (Downregulated)

calculate_log2FC <- function(normal_expr, tumor_expr) {
  log2fc <- log2(tumor_expr / normal_expr)
  return(log2fc)
}

# Call the function
calculate_log2FC(normal_expr = 2.10, tumor_expr = 5.80)   # TP53

# Apply to entire columns — R is vectorised, so it works on all 5 genes at once
expression_data$log2FC <- calculate_log2FC(expression_data$normal,
                                           expression_data$tumor)
expression_data

# ---------------------------------------------------------
#### if & else: Logical Conditions ####
# ---------------------------------------------------------
gene_expression <- 30

# If statement
if (gene_expression > 30) {
  print("Gene expression is high")
}

if (gene_expression > 50) {
  print("Gene expression is high")
}

# if-else statement
if (gene_expression > 50) {
  print("Gene expression is high")
} else {
  print("Gene expression is low")
}


gene_expression2 <- c(30,50,70,25,15)

if (gene_expression > 50) {
  print("Gene expression is high")
} else {
  print("Gene expression is low")
}

#throws an error. To solve, using for loop

# ---------------------------------------------------------
#### For Loop ####
# ---------------------------------------------------------

for (i in gene_expression2) {
  
  if (i > 30) {
    print("Gene expression if high")
  } else {
    print("Gene expression is low")
  }
}

# example 2
for (i in 1:5){
  print("my roll is: ")
  print(i)
}

for (i in 1:5){
  cat("my roll is:", i, "\n")
}


# example 3

logFC_values <- c(2.5, -1.8, 0.3, 1.2, -0.5)
for(value in logFC_values) {
  if (value > 1) {
    print("Upregulated Gene")
  } else if (value < -1) {
    print("Downregulated Gene")
  } else {
    print("No Major Change")
  }
}

genes <- c("TP53", "BRCA1", "EGFR", "MYC", "PTEN")

for (gene in genes) {
  print(gene)
}


# Loop over row indices (useful when you need to read AND write the same row)
for (i in 1:nrow(expression_data)) {
  if (expression_data$log2FC[i] > 1) {
    expression_data$status[i] <- "Upregulated"
  } else if (expression_data$log2FC[i] < -1) {
    expression_data$status[i] <- "Downregulated"
  } else {
    expression_data$status[i] <- "Not Significant"
  }
}

expression_data[ , c("gene", "log2FC", "status")]


# Save the R script
#   Ctrl + S   (or File → Save As)
#   Saves ONLY your code. When you reopen it, you must re-run everything
#   to rebuild objects in memory. Extension: .R

# Save the entire R workspace
save.image(file = "full_workspace.RData")
save.image(file = "Scripts/full_workspace.RData")

# Save selected objects only
save(expression_data, file = "Scripts/expression_data.RData")

# ---------------------------------------------------------
#### Data Wrangling ####
# ---------------------------------------------------------

data <- data.frame(
  patient_id     = c("P01","P02","P03","P04","P05"),
  age            = c(45, 62, NA, 71, NA),
  gender         = c("F","M","F","M","F"),
  disease_stage  = c("I","III","II","IV","II"),
  smoking        = c("No","Yes","No","Yes","No"),
  expression_MYC = c(7.2, NA, 9.1, NA, 5.5),
  stringsAsFactors = FALSE
)

# Write to CSV to datasets folder
write.csv(data, "Datasets/patient_data.csv", row.names = FALSE)

# Import the file 
data <- read.csv("Datasets/patient_data.csv")

# another technique
data <- read.csv(file.choose())


# To access from a specific file path
data <- read.csv("patient_data.csv") #throws error
raw_file <- "patient_data.csv"
folder <- "Datasets"
file_path <- file.path(folder, raw_file)
data2 <- read.csv(file_path)

# Before analyzing, inspect the dataset first
str(data)
head(data)
tail(data, 2)
dim(data)         # rows × columns
summary(data)     # min, max, mean, NA counts for every column
names(data)       # column names

data$patient_id
data[1:2, c(1, 3)]

# Create a new column
data$new_column <- c(1, 2, 3, 4, 5)

# ---------------------------------------------------------
#### Missing Values ####
# ---------------------------------------------------------

# Identify missing values
is.na(data)
sum(is.na(data))
colSums(is.na(data))
rowSums(is.na(data))

# Remove rows with missing values
clean_data1 <- na.omit(data)

# Remove columns with missing values
clean_data2 <- data[, colSums(is.na(data)) == 0]
View(clean_data2)

# Replace NA with 0
clean_data3 <- data
clean_data3[is.na(clean_data3)] <- 0

# Replace NA with mean
clean_data4 <- data
clean_data4$age[is.na(clean_data4$age)] <- mean(clean_data4$age, na.rm = TRUE)
clean_data4$expression_MYC[is.na(clean_data4$expression_MYC)] <- mean(clean_data4$expression_MYC, na.rm = TRUE)

# Confirm all NAs are gone
colSums(is.na(clean_data4))

# ---------------------------------------------------------
#### Practice Problem ####
# ---------------------------------------------------------
# 1. Define a function to calculate BMI
# Formula: weight (kg) divided by height squared (m^2)

calculate_BMI <- function(weight, height) {
  bmi <- weight / (height ^ 2)
  return(bmi)
}

# 2. Define a function to classify the weight status based on BMI
classify_weight <- function(bmi) {
  if (bmi > 30) {
    return("overweight")
  } else if (bmi < 25) {
    return("underweight")
  } else {
    return("normal weight")
  }
}

# 3. Setup the output folder environment
output_folder <- "Results"
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# 4. Read the raw dataset
data <- read.csv("Datasets/BMI_data_1.csv")
# OR
data <- read.csv(file.choose())

# 5. Clean Data: Replace all missing values (NA) in the table with 1
data[is.na(data)] <- 1

# 6. Initialize new empty columns to store our loop results
data$BMI    <- NA
data$status <- NA

# 7. Loop through the dataset row by row
for (i in 1:nrow(data)) {
  
  # Step A: Calculate BMI using our custom function
  data$BMI[i] <- calculate_BMI(data$weight[i], data$height[i])
  
  # Step B: Classify status using our classification function
  data$status[i] <- classify_weight(data$BMI[i])
  
}

# 8. Save the completed table into the Results folder
full_output_path <- file.path(output_folder, "Clean_BMI_Data_1.csv")
write.csv(data, full_output_path, row.names = FALSE)

# OR
write.csv(data, "Results/Clean_BMI_Data_1.csv")

# save the full workspace
save.image(file = "Scripts/full_workspace_class2.RData")
