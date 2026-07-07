# --- Package Installation ---
install.packages("caret")
install.packages("randomForest")
install.packages("RANN")
install.packages("Boruta")
install.packages("kernlab")
install.packages("pROC")
install.packages("ggplot2")
install.packages("glmnet")   
# install.packages("xgboost")  # Commented out due to small sample size constraints

# Libraries
library(caret)         
library(RANN)          
library(randomForest)  
library(Boruta)        
library(kernlab)       
library(pROC)          
library(ggplot2)       
library(glmnet)        
# library(xgboost)     # Commented out
library(dplyr)

# PHASE 2: Loading the RData file
load("C:/Users/iarif/Downloads/RNA-Seq/Module 4/GSE79973.RData")

dim(processed_data)
str(processed_data)

if(!is.null(dev.list())) dev.off()

groups <- phenotype_data %>%
  mutate(Condition = case_when(
    grepl("adenocarcinoma|tumor", source_name_ch1, ignore.case = TRUE) ~ "Cancer",
    grepl("mucosa|normal", source_name_ch1, ignore.case = TRUE)        ~ "Normal"
  )) %>%
  pull(Condition) %>%
  as.factor()

groups <- factor(groups, levels = c("Normal", "Cancer"))

data <- processed_data

dim(data)
str(data) 
View(data)

dev.off()
boxplot(data,
        outline = FALSE,
        col = "skyblue",
        main = "Before Log Transformation")

par(mfrow = c(1, 2))
hist(data[, 18], main = "Before Log Transform", col = "skyblue")
hist(log10(data + 1)[, 18], main = "After Log Transform", col = "coral")

data_log <- log10(data + 1)  

data_t <- as.data.frame(t(data_log))

# Near-Zero Variance filtering
nzv <- preProcess(data_t, method = "nzv", uniqueCut = 15)
data_t <- predict(nzv, data_t)

# Center and Scale Features
process <- preProcess(data_t, method = c("center", "scale"))
data_t <- predict(process, data_t)

# Looking for missing values and Imputing
anyNA(data_t)
sum(is.na(data_t))

knnImpute <- preProcess(data_t, method="knnImpute")
data_t <- predict(knnImpute, data_t)

sum(is.na(data_t))

class(groups)
levels(groups)
table(groups)

data_t <- cbind(groups, data_t)
# Correcting the typo 'grpups' to 'groups'
# Drop the duplicate and NA columns cleanly using dplyr
data_t <- data_t %>% dplyr::select(-groups.1, -groups.2)
# 1. Force the columns to be dropped using standard base R subsetting
data_t <- data_t[, !(colnames(data_t) %in% c("groups.1", "groups.2"))]

# 2. double-check your column positions to make sure it worked
print(names(data_t)[1:5])
View(data_t)

set.seed(123)
trainIndex <- createDataPartition(data_t$groups, p = 0.7, list = FALSE)
train_data <- data_t[trainIndex, ]
dim(train_data)

test_data  <- data_t[-trainIndex, ]
dim(test_data)

x_train <- train_data[, -1]  
y_train <- train_data$groups 

x_test  <- test_data[, -1]
y_test  <- test_data$groups

table(y_test)

dir.create("Results", showWarnings = FALSE)


# PHASE 1: BIOMARKER EXTRACTION (THREE INDEPENDENT ENGINES)


# Boruta Feature Selection
set.seed(123)
boruta_result <- Boruta(x = x_train, y = y_train, doTrace = 1)
final_boruta <- TentativeRoughFix(boruta_result) 
boruta_features <- getSelectedAttributes(final_boruta)
write.csv(boruta_features, "Results/boruta_selected_genes.csv", row.names = FALSE)


#  RFE Feature Selection 
set.seed(123)
rfe_control <- rfeControl(functions = rfFuncs, method = "cv", number = 5)
rfe_result <- rfe(x = x_train, y = y_train, sizes = seq(5, 100, by = 5), rfeControl = rfe_control)
plot(rfe_result, type = c("g", "o"))
rfe_features <- predictors(rfe_result)
write.csv(rfe_features, "Results/rfe_selected_genes.csv", row.names = FALSE)


# LASSO Feature Selection (Biomarker Extractor) 
set.seed(123)
# Convert data frames to matrix format required natively by glmnet
x_train_matrix <- as.matrix(x_train)
# Convert factor to binary numeric (0 and 1) for binomial LASSO
y_train_binary <- as.numeric(y_train) - 1 

# Run Cross-Validated LASSO to find the optimal penalty value (lambda.min)
cv_lasso <- cv.glmnet(x = x_train_matrix, y = y_train_binary, family = "binomial", alpha = 1)
plot(cv_lasso)

# Extract coefficients at the optimal lambda parameter
lasso_coefs <- coef(cv_lasso, s = "lambda.min")
# Isolate genes whose predictive coefficients were NOT shrunk to absolute zero
lasso_matrix <- as.matrix(lasso_coefs)
lasso_features <- rownames(lasso_matrix)[lasso_matrix[, 1] != 0]
lasso_features <- lasso_features[lasso_features != "(Intercept)"] # Remove intercept baseline entry

write.csv(lasso_features, "Results/lasso_selected_genes.csv", row.names = FALSE)


# --Gene Intersection ---
common_features <- intersect(intersect(boruta_features, rfe_features), lasso_features)
print(paste("Number of core consensus genes matched by all 3 extractors:", length(common_features)))
write.csv(common_features, "Results/three_way_common_genes.csv", row.names = FALSE)



# PHASE 2: ML CLASSIFIER BENCHMARKS (Running Remaining Core Models)

# We will now benchmark the remaining algorithms using the Boruta sets,
# the RFE sets, and your brand-new LASSO selected feature sets.

ctrl <- trainControl(method = "cv", number = 10, verboseIter = TRUE, classProbs = TRUE)
tune_Grid_svm <- expand.grid(C = c(0.25, 0.5, 1, 2, 4, 8, 16), sigma = c(0.01, 0.05, 0.1, 0.5, 1))
tune_Grid_ann <- expand.grid(size = c(1, 3, 5, 10), decay = c(0, 0.001, 0.01, 0.1))
# xgb_grid <- expand.grid(nrounds = c(50, 100, 150), max_depth = c(2, 3, 4), eta = c(0.01, 0.05, 0.1), gamma = 0, colsample_bytree = c(0.6, 0.8), min_child_weight = 1, subsample = c(0.7, 0.9))


# SECTION A: TRAINING WITH BORUTA FEATURES

set.seed(123)
train_boruta <- x_train[, boruta_features]
test_boruta  <- x_test[, boruta_features]

RF_boruta  <- train(x = train_boruta, y = y_train, method = "rf", importance = TRUE, trControl = ctrl)
SVM_boruta <- train(x = train_boruta, y = y_train, method = "svmRadial", trControl = ctrl, tuneGrid = tune_Grid_svm, prob.model = TRUE)
ANN_boruta <- train(x = train_boruta, y = y_train, method = "nnet", trControl = ctrl, tuneGrid = tune_Grid_ann, MaxNWts = 5000, trace = FALSE)
# XGB_boruta <- train(x = train_boruta, y = y_train, method = "xgbTree", trControl = ctrl, tuneGrid = xgb_grid, verbose = FALSE)

# ------------------------------------------------------------------------------
# SECTION B: TRAINING WITH RFE FEATURES
# ------------------------------------------------------------------------------
set.seed(123)
train_rfe <- x_train[, rfe_features]
test_rfe  <- x_test[, rfe_features]

RF_rfe  <- train(x = train_rfe, y = y_train, method = "rf", importance = TRUE, trControl = ctrl, tuneGrid = data.frame(mtry = 1:min(5, ncol(train_rfe))))
SVM_rfe <- train(x = train_rfe, y = y_train, method = "svmRadial", trControl = ctrl, tuneGrid = expand.grid(C = c(0.01, 0.1, 1, 2, 4), sigma = c(0.01, 0.05, 0.1, 0.5)), prob.model = TRUE)
ANN_rfe <- train(x = train_rfe, y = y_train, method = "nnet", trControl = ctrl, tuneGrid = expand.grid(size = c(1,2,3,5), decay = c(0, 0.001, 0.01, 0.1)), MaxNWts = 5000, trace = FALSE)
# XGB_rfe <- train(x = train_rfe, y = y_train, method = "xgbTree", trControl = ctrl, tuneGrid = xgb_grid, verbose = FALSE)


# SECTION C: TRAINING WITH LASSO-EXTRACTED FEATURES

set.seed(123)
train_lasso <- x_train[, lasso_features]
test_lasso  <- x_test[, lasso_features]

RF_lasso  <- train(x = train_lasso, y = y_train, method = "rf", importance = TRUE, trControl = ctrl)
SVM_lasso <- train(x = train_lasso, y = y_train, method = "svmRadial", trControl = ctrl, tuneGrid = tune_Grid_svm, prob.model = TRUE)
ANN_lasso <- train(x = train_lasso, y = y_train, method = "nnet", trControl = ctrl, tuneGrid = tune_Grid_ann, MaxNWts = 5000, trace = FALSE)
# XGB_lasso <- train(x = train_lasso, y = y_train, method = "xgbTree", trControl = ctrl, tuneGrid = xgb_grid, verbose = FALSE)



# PHASE 3: EVALUATION AND BENCHMARK SUMMARY CONTROLLERS

# Run testing cohort predictions across all feature tiers
rf_pred_b   <- predict(RF_boruta, test_boruta);   rf_pred_r   <- predict(RF_rfe, test_rfe);   rf_pred_l   <- predict(RF_lasso, test_lasso)
svm_pred_b  <- predict(SVM_boruta, test_boruta);  svm_pred_r  <- predict(SVM_rfe, test_rfe);  svm_pred_l  <- predict(SVM_lasso, test_lasso)
ann_pred_b  <- predict(ANN_boruta, test_boruta);  ann_pred_r  <- predict(ANN_rfe, test_rfe);  ann_pred_l  <- predict(ANN_lasso, test_lasso)
# xgb_pred_b  <- predict(XGB_boruta, test_boruta);  xgb_pred_r  <- predict(XGB_rfe, test_rfe);  xgb_pred_l  <- predict(XGB_lasso, test_lasso)

# Generate final accuracy benchmarking matrices (XGBoost vectors removed)
final_comparison_table <- data.frame(
  Model = c("Random Forest", "SVM", "ANN"),
  Boruta_Accuracy = c(confusionMatrix(rf_pred_b, y_test)$overall["Accuracy"],
                      confusionMatrix(svm_pred_b, y_test)$overall["Accuracy"],
                      confusionMatrix(ann_pred_b, y_test)$overall["Accuracy"]),
  RFE_Accuracy    = c(confusionMatrix(rf_pred_r, y_test)$overall["Accuracy"],
                      confusionMatrix(svm_pred_r, y_test)$overall["Accuracy"],
                      confusionMatrix(ann_pred_r, y_test)$overall["Accuracy"]),
  LASSO_Accuracy  = c(confusionMatrix(rf_pred_l, y_test)$overall["Accuracy"],
                      confusionMatrix(svm_pred_l, y_test)$overall["Accuracy"],
                      confusionMatrix(ann_pred_l, y_test)$overall["Accuracy"])
)

print("--- Complete Accuracy Benchmarks Across All 3 Biomarker Extraction Strategies ---")
print(final_comparison_table)



# PHASE 4: CONSENSUS ROC CURVE PLOTTING

target_class <- levels(y_test)[2]

# Extract probabilities for the newly introduced LASSO cohort tier
rf_prob_l  <- predict(RF_lasso, newdata = test_lasso, type = "prob")
svm_prob_l <- predict(SVM_lasso, newdata = test_lasso, type = "prob")
ann_prob_l <- predict(ANN_lasso, newdata = test_lasso, type = "prob")
# xgb_prob_l <- predict(XGB_lasso, newdata = test_lasso, type = "prob")

rf_roc_l  <- roc(y_test, rf_prob_l[, target_class], levels=levels(y_test))
svm_roc_l <- roc(y_test, svm_prob_l[, target_class], levels=levels(y_test))
ann_roc_l <- roc(y_test, ann_prob_l[, target_class], levels=levels(y_test))
# xgb_roc_l <- roc(y_test, xgb_prob_l[, target_class], levels=levels(y_test))

# Plot performance curves exclusively for your new LASSO-extracted feature signature
dev.off()
plot(rf_roc_l, col = "black", lwd = 2, main = "ROC Performance: LASSO-Selected Biomarkers")
plot(svm_roc_l, col = "red", lwd = 2, add = TRUE)
plot(ann_roc_l, col = "darkgreen", lwd = 2, add = TRUE)
# plot(xgb_roc_l, col = "darkorange", lwd = 2, add = TRUE)

legend("bottomright", 
       legend=c(paste("RF (AUC =", round(auc(rf_roc_l), 2), ")"),
                paste("SVM (AUC =", round(auc(svm_roc_l), 2), ")"),
                paste("ANN (AUC =", round(auc(ann_roc_l), 2), ")")), 
       col=c("black","red","darkgreen"), lwd=2, bty="n", cex=0.8)