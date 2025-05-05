# Build Rad-score model ####
# Load data ####
# PFV
# train_lassoed <- read.csv("0-2Radiomics_feature_Slected/train_PFV_lassoed.csv", header = TRUE)
# valid_lassoed <- read.csv("0-2Radiomics_feature_Slected/valid_PFV_lassoed.csv", header = TRUE)
# test_lassoed <- read.csv("0-2Radiomics_feature_Slected/test_PFV_lassoed.csv", header = TRUE)
# exvalid1_lassoed <- read.csv("0-2Radiomics_feature_Slected/exvalid1_PFV_lassoed.csv", header = TRUE)
# exvalid2_lassoed <- read.csv("0-2Radiomics_feature_Slected/exvalid2_PFV_lassoed.csv", header = TRUE)
# genevalid_lassoed <- read.csv("0-2Radiomics_feature_Slected/genevalid_PFV_lassoed.csv", header = TRUE)
## Tumor
train_lassoed <- read.csv("0-2Radiomics_feature_Slected/train_T_lassoed.csv", header = TRUE)
valid_lassoed <- read.csv("0-2Radiomics_feature_Slected/valid_T_lassoed.csv", header = TRUE)
test_lassoed <- read.csv("0-2Radiomics_feature_Slected/test_T_lassoed.csv", header = TRUE)
exvalid1_lassoed <- read.csv("0-2Radiomics_feature_Slected/exvalid1_T_lassoed.csv", header = TRUE)
exvalid2_lassoed <- read.csv("0-2Radiomics_feature_Slected/exvalid2_T_lassoed.csv", header = TRUE)
genevalid_lassoed <- read.csv("0-2Radiomics_feature_Slected/genevalid_T_lassoed.csv", header = TRUE)
# Load required packages; install them first if not already installed
library(dplyr)     # For data manipulation
library(caret)     # For tuning parameters and evaluating model metrics
library(pROC)      # For plotting ROC curves
library(ggplot2)   # For plotting
library(ggpubr)    # For enhanced plotting
library(ggprism)   # For Prism-style plots


## Model 1: Logistic Regression (LR) ####
# Load necessary packages; install if not already installed
library(caTools) # For LR model

### Training Set ####
X_train <- train_lassoed[, c(4:ncol(train_lassoed))] # Select all feature columns as independent variables
y_train <- as.factor(train_lassoed$Status) # Select outcome variable as dependent variable
# Build logistic regression model using training set
LR_model <- glm(y_train ~ ., data = X_train, family = binomial)
# Predict using best parameter configuration on training set
train_predictions <- predict(LR_model, newdata = X_train, cluster = "response")

## Export predicted probabilities to CSV
# Convert probabilities to data frame
probabilities_df <- as.data.frame(train_predictions)
# Set column names as class labels
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
# Keep ID and Status columns
train_LRed <- cbind(ID = train_lassoed$ID, Status = train_lassoed$Status, probabilities_df)

## Model Evaluation ##
library(caret) # Powerful package for detailed metrics like SPE, PPV, ACC
library(pROC)
# Compute ROC
roc_obj <- roc(train_lassoed$Status, train_predictions)
auc_ci <- ci.auc(roc_obj)
roc_obj
auc_ci
# Choose optimal threshold based on ROC
best_threshold <- coords(roc_obj, "best", ret = "threshold")
best_threshold # View threshold
event_probability <- -2.186137  # PFV -2.433034; TUMOR -2.186137
# Recompute confusion matrix using ROC-based threshold
confusion_matrix <- confusionMatrix(
  as.factor(ifelse(train_predictions > event_probability, 1, 0)), 
  as.factor(y_train), 
  positive = "1" # Set 1 as the positive class
)
confusion_matrix # Print confusion matrix

## Plot confusion matrix heatmap ###
# install.packages("reshape2")
library(reshape2)
# Convert confusion matrix to data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
# Prepare data for plotting using raw counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
# Plot
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the training cohort") +
  theme_prism(border = T) +
  theme(
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
library(eoffice)  # Load package
# topptx(filename = "4Figure_PPT/LR_PFV_confusion_matrix_train.pptx")
topptx(filename = "4Figure_PPT/LR_T_confusion_matrix_train.pptx")

### Validation Set ####
X_valid <- valid_lassoed[, c(4:ncol(valid_lassoed))] # Select all feature columns as independent variables
y_valid <- as.factor(valid_lassoed$Status) # Select outcome variable
# Predict on validation set using best parameter configuration
valid_predictions <- predict(LR_model, newdata = X_valid, cluster = "response")

## Export predicted probabilities to CSV
# Convert probabilities to data frame
probabilities_df <- as.data.frame(valid_predictions)
# Set column names as class labels
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
# Keep ID and Status columns
valid_LRed <- cbind(ID = valid_lassoed$ID, Status = valid_lassoed$Status, probabilities_df)

## Model Evaluation ##
library(caret) # Powerful package for detailed metrics like SPE, PPV, ACC
# Compute ROC
roc_obj <- roc(valid_lassoed$Status, valid_predictions)
auc_ci <- ci.auc(roc_obj)
roc_obj
auc_ci
# Output confusion matrix
confusion_matrix <- confusionMatrix(
  as.factor(ifelse(valid_predictions > event_probability, 1, 0)), # Stratify predictions using threshold
  y_valid <- as.factor(y_valid), # Convert outcome to factor
  positive = "1" # Set 1 as the positive class
)
confusion_matrix # Print confusion matrix

## Plot confusion matrix heatmap ###
# install.packages("reshape2")
library(reshape2)
# Convert confusion matrix to data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
# Prepare data for plotting using raw counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
# Plot
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the internal validation cohort") +
  theme_prism(border = T) +
  theme(
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
library(eoffice) # Load package
# topptx(filename = "4Figure_PPT/LR_PFV_confusion_matrix_valid.pptx")
topptx(filename = "4Figure_PPT/LR_T_confusion_matrix_valid.pptx")


### Testing Set ####
X_test <- test_lassoed[, c(4:ncol(test_lassoed))] # Select all feature columns as predictors
y_test <- as.factor(test_lassoed$Status) # Select the outcome variable

# Predict using the best parameter configuration on the test set
test_predictions <- predict(LR_model, newdata = X_test, cluster = "response")

## Output probabilities to CSV
# Convert probabilities to data frame
probabilities_df <- as.data.frame(test_predictions)
# Set column names to class names
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
# Retain ID and Status columns
test_LRed <- cbind(ID = test_lassoed$ID, Status = test_lassoed$Status, probabilities_df)

## Model Evaluation ##
library(caret) # caret package provides detailed evaluation metrics like SPE, PPV, ACC, etc.
# Output confusion matrix
confusion_matrix <- confusionMatrix(as.factor(ifelse(test_predictions > event_probability, 1, 0)), # Stratify predictions with threshold = 0.5
                                    y_test <- as.factor(y_test), # Convert outcome variable to factor
                                    positive = "1") # Set class 1 as positive
confusion_matrix # Print confusion matrix

# Calculate ROC
roc_obj <- roc(test_lassoed$Status, test_predictions)
auc_ci <- ci.auc(roc_obj)
roc_obj
auc_ci

## Draw confusion matrix heatmap ###
# install.packages("reshape2")
library(reshape2)
# Convert confusion matrix to data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
# Prepare plotting data using raw counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)

# Plot
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the testing cohort") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

library(eoffice) # Load package
# topptx(filename = "4Figure_PPT/LR_PFV_confusion_matrix_test.pptx")
topptx(filename = "4Figure_PPT/LR_T_confusion_matrix_test.pptx")


### External Testing Set 1 ####
# Predict using the best parameter configuration on the external set
X_exvalid1 <- exvalid1_lassoed[, c(4:ncol(exvalid1_lassoed))] # Select all feature columns as predictors
y_exvalid1 <- as.factor(exvalid1_lassoed$Status) # Select the outcome variable

# Predict using the best parameter configuration
exvalid1_predictions <- predict(LR_model, newdata = X_exvalid1, cluster = "response")

## Output probabilities to CSV
# Convert probabilities to data frame
probabilities_df <- as.data.frame(exvalid1_predictions)
# Set column names to class names
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
# Retain ID and Status columns
exvalid1_LRed <- cbind(ID = exvalid1_lassoed$ID, Status = exvalid1_lassoed$Status, probabilities_df)

## Model Evaluation
library(caret) # caret package provides detailed evaluation metrics like SPE, PPV, ACC, etc.
# Output confusion matrix
confusion_matrix <- confusionMatrix(as.factor(ifelse(exvalid1_predictions > event_probability, 1, 0)), # Stratify predictions with threshold = 0.5
                                    y_exvalid1 <- as.factor(y_exvalid1), # Convert outcome variable to factor
                                    positive = "1") # Set class 1 as positive
confusion_matrix # Print confusion matrix

# Calculate ROC
roc_obj <- roc(exvalid1_lassoed$Status, exvalid1_predictions)
auc(roc_obj)

## Draw confusion matrix heatmap ###
# install.packages("reshape2")
library(reshape2)
# Convert confusion matrix to data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
# Prepare plotting data using raw counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)

# Plot
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the external validation cohort 1") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

library(eoffice) # Load package
# topptx(filename = "4Figure_PPT/LR_PFV_confusion_matrix_exvalid1.pptx")
topptx(filename = "4Figure_PPT/LR_T_confusion_matrix_exvalid1.pptx")

### External test set 2 ####
# Predict on the test set using the best parameter configuration
X_exvalid2<-exvalid2_lassoed[,c(4:ncol(exvalid2_lassoed))] # Select all columns of features as independent variables
y_exvalid2<-as.factor(exvalid2_lassoed$Status) # Select dependent variable
# Predict on the validation set using the best parameter configuration
exvalid2_predictions <- predict(LR_model, newdata =X_exvalid2,cluster="response")

## Output probabilities to CSV
# Convert probabilities to data frame
probabilities_df <- as.data.frame(exvalid2_predictions)
# Set column names to class names
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
# Keep ID and Status columns
exvalid2_LRed <- cbind(ID =exvalid2_lassoed$ID, Status =exvalid2_lassoed$Status,probabilities_df)

## Model evaluation
library(caret) # caret package provides powerful and detailed evaluation metrics like SPE, PPV, ACC, etc.
# Output confusion matrix
confusion_matrix<-confusionMatrix(as.factor(ifelse(exvalid2_predictions>event_probability,1,0)), # Stratify predictions, set threshold to 0.5
                                  y_exvalid2 <- as.factor(y_exvalid2), # Convert outcome variable to factor
                                  positive = "1") # Set 1 as positive
confusion_matrix # Output confusion matrix
# Calculate ROC
roc_obj <- roc(exvalid2_lassoed$Status, exvalid2_predictions)
auc(roc_obj)

## Plot confusion matrix heatmap ###
#install.packages("reshape2")
library(reshape2)
# Convert confusion matrix to data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative","Positive")
rownames(confusion_matrix_df) <- c("Negative","Positive")
# Prepare data for plotting using raw counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
# Plot
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the external validation cohort 2") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
library(eoffice) # Load package
#topptx(filename = "4Figure_PPT/LR_PFV_confusion_matrix_exvalid2.pptx")
topptx(filename = "4Figure_PPT/LR_T_confusion_matrix_exvalid2.pptx")


### Gene validation set ####
X_genevalid<-genevalid_lassoed[,c(4:ncol(genevalid_lassoed))] # Select all columns of features as independent variables
y_genevalid<-as.factor(genevalid_lassoed$Status) # Select dependent variable
# Predict on the validation set using the best parameter configuration
genevalid_predictions <- predict(LR_model, newdata =X_genevalid,cluster="response")

## Output probabilities to CSV
probabilities_df <- as.data.frame(genevalid_predictions) # Convert probabilities to data frame
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df)) # Set column names to class names
# Keep ID and Status columns
genevalid_LRed <- cbind(ID =genevalid_lassoed$ID, Status =genevalid_lassoed$Status,probabilities_df)


## Model evaluation
library(caret) # caret package provides powerful and detailed evaluation metrics like SPE, PPV, ACC, etc.
# Output confusion matrix
confusion_matrix<-confusionMatrix(as.factor(ifelse(genevalid_predictions>event_probability,1,0)), # Stratify predictions, set threshold to 0.5
                                  genevalid <- as.factor(y_genevalid), # Convert outcome variable to factor
                                  positive = "1") # Set 1 as positive
confusion_matrix # Output confusion matrix
# Calculate ROC
roc_obj <- roc(genevalid_lassoed$Status, genevalid_predictions)
auc(roc_obj)

## Plot confusion matrix heatmap ###
#install.packages("reshape2")
library(reshape2)
# Convert confusion matrix to data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative","Positive")
rownames(confusion_matrix_df) <- c("Negative","Positive")
# Prepare data for plotting using raw counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
# Plot
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the gene cohort") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
library(eoffice) # Load package
#topptx(filename = "4Figure_PPT/LR_PFV_confusion_matrix_genevalid.pptx")
topptx(filename = "4Figure_PPT/LR_T_confusion_matrix_genevalid.pptx")


### Save Radiomics_Score to CSV ####
# PFV
write.csv(train_LRed, "0-3Radiomics_Score_Result/train_PFV_LRed.csv", row.names = FALSE)
write.csv(valid_LRed, "0-3Radiomics_Score_Result/valid_PFV_LRed.csv", row.names = FALSE)
write.csv(test_LRed, "0-3Radiomics_Score_Result/test_PFV_LRed.csv", row.names = FALSE)
write.csv(exvalid1_LRed, "0-3Radiomics_Score_Result/exvalid1_PFV_LRed.csv", row.names = FALSE)
write.csv(exvalid2_LRed, "0-3Radiomics_Score_Result/exvalid2_PFV_LRed.csv", row.names = FALSE)
write.csv(genevalid_LRed, "0-3Radiomics_Score_Result/genevalid_PFV_LRed.csv", row.names = FALSE)
# Tumor
#write.csv(train_LRed, "0-3Radiomics_Score_Result/train_T_LRed.csv", row.names = FALSE)
#write.csv(valid_LRed, "0-3Radiomics_Score_Result/valid_T_LRed.csv", row.names = FALSE)
#write.csv(test_LRed, "0-3Radiomics_Score_Result/test_T_LRed.csv", row.names = FALSE)
#write.csv(exvalid1_LRed, "0-3Radiomics_Score_Result/exvalid1_T_LRed.csv", row.names = FALSE)
#write.csv(exvalid2_LRed, "0-3Radiomics_Score_Result/exvalid2_T_LRed.csv", row.names = FALSE)
#write.csv(genevalid_LRed, "0-3Radiomics_Score_Result/genevalid_T_LRed.csv", row.names = FALSE)


##Model 2: Construct SVM Model####
#Support Vector Machine (SVM) is a classic supervised learning algorithm used for classification and regression tasks.
#The principle of SVM is based on the structural risk minimization principle in statistical learning theory, aiming to find an optimal hyperplane (or surface) that separates samples of different classes as much as possible.
#Factor the classification variable and assign labels
#train_lassoed$Status <- factor(train_lassoed$Status, labels = c("N0","N2"))
#valid_lassoed$Status <- factor(valid_lassoed$Status, labels = c("N0","N2"))
#test_lassoed$Status <- factor(test_lassoed$Status, labels = c("N0","N2"))
###Training Set####
#Train the model on the training set
library(e1071) #SVM model usage
#Define feature and target variable for the training set
X_train <- train_lassoed[, c(4:ncol(train_lassoed))] # Select all columns of features as independent variables
y_train <- as.factor(train_lassoed$Status) # Select dependent variable
#Create parameter grid
set.seed(123)
param_grid <- expand.grid(
  sigma = 10^(-5:5), # Sigma (gamma in SVM) is a parameter of the Radial Basis Function (RBF) kernel. It controls the spread of data in feature space. Smaller sigma values lead to potential underfitting, while larger values may cause overfitting.
  #kernel = c('linear', 'polynomial', 'radial', 'sigmoid'), # svmRadial, svmLinear
  C = 10^(-5:5)) # C (cost in SVM) is a penalty parameter that controls the punishment for misclassified samples. Smaller C values allow more misclassifications, making the decision boundary adapt more easily to training data, but potentially leading to poor generalization; larger C values impose stricter penalties for misclassification, resulting in a stricter decision boundary but potentially overfitting the model to the training data.
#Define cross-validation control parameters, here using 5-fold cross-validation
ctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE)
#Perform parameter tuning
tuned_model <- train(x = X_train, y = y_train,
                     method = "svmRadial",
                     probability = TRUE,
                     tuneGrid = param_grid,
                     trControl = ctrl)
#Output the best parameter configuration
print(tuned_model)
##Train the model with the best parameters
svm_final_model <- svm(x = X_train, y = y_train,
                       kernel = 'radial',       #'linear', 'polynomial', 'radial', 'sigmoid'
                       gamma = 0.1,  # PFV: 0.1; Tumor: 0.1
                       cost = 10,   # PFV: 10; Tumor: 10
                       probability = TRUE)
train_predictions <- predict(svm_final_model, newdata = X_train)

##Output decision values and probabilities to CSV for preparation
#Calculate decision values and probabilities:
train_decisionvalues <- attr(predict(svm_final_model, X_train, decision.values = TRUE), "decision.values") #Output SVM predicted decision values, which can be interpreted as confidence. The larger the absolute value, the higher the confidence in predicting the class.
train_pred <- attr(predict(svm_final_model, X_train, probability = TRUE), "probabilities") #Output SVM predicted probabilities, ranging from 0 to 1
#Convert decision values and probabilities to data frames
decision_values_df <- as.data.frame(train_decisionvalues)
probabilities_df <- as.data.frame(train_pred)
#Set the column names of the probability values to class names
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
colnames(decision_values_df) <- paste0("decision_values", colnames(decision_values_df))
#Retain ID and Status columns
train_svmed <- cbind(ID = train_lassoed$ID, Status = train_lassoed$Status, decision_values_df, probabilities_df)

##Model Evaluation##
library(caret) # The caret package provides powerful and detailed output for evaluation metrics like SPE, PPV, ACC, etc.
library(pROC)
#Calculate ROC
roc_obj <- roc(response = y_train, predictor = train_pred[, 2]) #Select the predicted probability for class 1
auc_ci <- ci.auc(roc_obj)
roc_obj
auc_ci
#Select the best threshold based on ROC to adjust predicted values
best_threshold <- coords(roc_obj, "best", ret = "threshold")
best_threshold #View threshold
event_probability <- 0.09179674  # PFV 0.09221858  # T0.09179674
#Construct confusion matrix
confusion_matrix <- confusionMatrix(as.factor(ifelse(train_pred[, 2] > event_probability, 1, 0)), #Stratify the predicted values, set the threshold to 0.5
                                    y_train <- as.factor(y_train), #Convert the outcome variable to a factor
                                    positive = "1") #Set 1 as positive
confusion_matrix #Output confusion matrix


##Plot Confusion Matrix Heatmap###
#install.packages("reshape2")
library(reshape2)
#Convert confusion matrix to data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
#Prepare data for plotting, using original counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
#Plotting
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the training cohort") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
library(eoffice) # Load package
topptx(filename = "4Figure_PPT/SVM_confusion_matrix_train.pptx")


###Validation Set####
#Make predictions on the validation set using the best parameter configuration
X_valid <- valid_lassoed[, c(4:ncol(valid_lassoed))] # Select all columns except the first and second as independent variables
y_valid <- as.factor(valid_lassoed$Status) # Select the second column as the dependent variable
valid_predictions <- predict(svm_final_model, newdata = X_valid) #Use the best parameter configuration from the training set for prediction on the validation set

##Output decision values and probabilities to CSV for preparation
#Calculate decision values and probabilities:
valid_decisionvalues <- attr(predict(svm_final_model, X_valid, decision.values = TRUE), "decision.values") #Output SVM predicted decision values
valid_pred <- attr(predict(svm_final_model, X_valid, probability = TRUE), "probabilities") #Output SVM predicted probabilities
#Convert decision values and probabilities to data frames
decision_values_df <- as.data.frame(valid_decisionvalues)
probabilities_df <- as.data.frame(valid_pred)
#Set the column names of the probability values to class names
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
colnames(decision_values_df) <- paste0("decision_values", colnames(decision_values_df))
#Retain ID and Status columns
valid_svmed <- cbind(ID = valid_lassoed$ID, Status = valid_lassoed$Status, decision_values_df, probabilities_df)

##Model Evaluation
library(caret) # The caret package provides powerful and detailed output for evaluation metrics like SPE, PPV, ACC, etc.
#Calculate ROC
roc_obj <- roc(response = y_valid, predictor = valid_pred[, 2])
auc_ci <- ci.auc(roc_obj)
roc_obj
auc_ci
#Construct confusion matrix
confusion_matrix <- confusionMatrix(as.factor(ifelse(valid_pred[, 2] > event_probability, 1, 0)), #Stratify the predicted values, set the threshold to 0.5
                                    y_valid <- as.factor(y_valid), #Convert the outcome variable to a factor
                                    positive = "1") #Set 1 as positive
confusion_matrix #Output confusion matrix


##Plot Confusion Matrix Heatmap###
#install.packages("reshape2")
library(reshape2)
#Convert confusion matrix to data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
#Prepare data for plotting, using original counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
#Plotting
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the internal validation cohort") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
library(eoffice) # Load package
topptx(filename = "4Figure_PPT/SVM_confusion_matrix_valid.pptx")


###Test Set####
#Make predictions on the test set using the best parameter configuration
X_test <- test_lassoed[, c(4:ncol(test_lassoed))] # Select all columns except the first and second as independent variables
y_test <- as.factor(test_lassoed$Status) # Select the second column as the dependent variable
test_predictions <- predict(svm_final_model, newdata = X_test) #Use the best parameter configuration from the training set for prediction on the validation set

##Output decision values and probabilities to CSV for preparation
#Calculate decision values and probabilities:
test_decisionvalues <- attr(predict(svm_final_model, X_test, decision.values = TRUE), "decision.values") #Output SVM predicted decision values
test_pred <- attr(predict(svm_final_model, X_test, probability = TRUE), "probabilities") #Output SVM predicted probabilities
#Convert decision values and probabilities to data frames
decision_values_df <- as.data.frame(test_decisionvalues)
probabilities_df <- as.data.frame(test_pred)
#Set the column names of the probability values to class names
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
colnames(decision_values_df) <- paste0("decision_values", colnames(decision_values_df))
#Retain ID and Status columns
test_svmed <- cbind(ID = test_lassoed$ID, Status = test_lassoed$Status, decision_values_df, probabilities_df)

##Model Evaluation
library(caret) # The caret package provides powerful and detailed output for evaluation metrics like SPE, PPV, ACC, etc.
#Calculate ROC
roc_obj <- roc(response = y_test, predictor = test_pred[, 2])
auc_ci <- ci.auc(roc_obj)
roc_obj
auc_ci
#Construct confusion matrix
confusion_matrix <- confusionMatrix(as.factor(ifelse(test_pred[, 2] > event_probability, 1, 0)), #Stratify the predicted values, set the threshold to 0.5
                                    y_test <- as.factor(y_test), #Convert the outcome variable to a factor
                                    positive = "1") #Set 1 as positive
confusion_matrix #Output confusion matrix


##Plot Confusion Matrix Heatmap###
#install.packages("reshape2")
library(reshape2)
#Convert confusion matrix to data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
#Prepare data for plotting, using original counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
#Plotting
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the external validation cohort") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
library(eoffice) # Load package
topptx(filename = "4Figure_PPT/SVM_confusion_matrix_test.pptx")


### External Validation Set 1 ####
# Predict on the test set using the best parameter configuration
X_exvalid1 <- exvalid1_lassoed[c(4:ncol(exvalid1_lassoed))] # Select all columns except the first and second as independent variables
y_exvalid1 <- as.factor(exvalid1_lassoed$Status) # Select the second column as the dependent variable
exvalid1_predictions <- predict(svm_final_model, newdata = X_exvalid1)

## Output decision values and probabilities to CSV for preparation
# Calculate decision values and probabilities:
exvalid1_decisionvalues <- attr(predict(svm_final_model, X_exvalid1, decision.values = TRUE), "decision.values") # Output SVM predicted decision values; decision values can also be understood as confidence, with higher absolute values indicating higher confidence for the predicted class
exvalid1_pred <- attr(predict(svm_final_model, X_exvalid1, probability = TRUE), "probabilities") # Output SVM predicted probabilities, between 0 and 1
# Convert decision values and probabilities to data frames
decision_values_df <- as.data.frame(exvalid1_decisionvalues)
probabilities_df <- as.data.frame(exvalid1_pred)
# Set the column names of the probability values to class names
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
colnames(decision_values_df) <- paste0("decision_values", colnames(decision_values_df))
# Retain ID and Status columns
exvalid1_svmed <- cbind(ID = exvalid1_lassoed$ID, Status = exvalid1_lassoed$Status, decision_values_df, probabilities_df)

## Model Evaluation ##
library(caret) # The caret package provides powerful and detailed output for evaluation metrics like SPE, PPV, ACC, etc.
# Calculate ROC
roc_obj <- roc(response = y_exvalid1, predictor = exvalid1_pred[, 2])
auc(roc_obj)
# Construct confusion matrix
confusion_matrix <- confusionMatrix(as.factor(ifelse(exvalid1_pred[, 2] > event_probability, 1, 0)), # Stratify the predicted values, set the threshold to 0.5
                                    y_exvalid1 <- as.factor(y_exvalid1), # Convert the outcome variable to a factor
                                    positive = "1") # Set 1 as positive
confusion_matrix # Output confusion matrix

## Plot Confusion Matrix Heatmap ###
# install.packages("reshape2")
library(reshape2)
# Convert confusion matrix to data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
# Prepare data for plotting, using original counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
# Plotting
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the internal validation cohort") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
library(eoffice) # Load package
topptx(filename = "4Figure_PPT/SVM_confusion_matrix_exvalid1.pptx")


### External Validation Set 2 ####
# Predict on the test set using the best parameter configuration
X_exvalid2 <- exvalid2_lassoed[c(4:ncol(exvalid2_lassoed))] # Select all columns except the first and second as independent variables
y_exvalid2 <- as.factor(exvalid2_lassoed$Status) # Select the second column as the dependent variable
exvalid2_predictions <- predict(svm_final_model, newdata = X_exvalid2)

## Output decision values and probabilities to CSV for preparation
# Calculate decision values and probabilities:
exvalid2_decisionvalues <- attr(predict(svm_final_model, X_exvalid2, decision.values = TRUE), "decision.values") # Output SVM predicted decision values; decision values can also be understood as confidence, with higher absolute values indicating higher confidence for the predicted class
exvalid2_pred <- attr(predict(svm_final_model, X_exvalid2, probability = TRUE), "probabilities") # Output SVM predicted probabilities, between 0 and 1
# Convert decision values and probabilities to data frames
decision_values_df <- as.data.frame(exvalid2_decisionvalues)
probabilities_df <- as.data.frame(exvalid2_pred)
# Set the column names of the probability values to class names
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
colnames(decision_values_df) <- paste0("decision_values", colnames(decision_values_df))
# Retain ID and Status columns
exvalid2_svmed <- cbind(ID = exvalid2_lassoed$ID, Status = exvalid2_lassoed$Status, decision_values_df, probabilities_df)

## Model Evaluation ##
library(caret) # The caret package provides powerful and detailed output for evaluation metrics like SPE, PPV, ACC, etc.
# Calculate ROC
roc_obj <- roc(response = y_exvalid2, predictor = exvalid2_pred[, 2])
auc(roc_obj)
# Construct confusion matrix
confusion_matrix <- confusionMatrix(as.factor(ifelse(exvalid2_pred[, 2] > event_probability, 1, 0)), # Stratify the predicted values, set the threshold to 0.5
                                    y_exvalid2 <- as.factor(y_exvalid2), # Convert the outcome variable to a factor
                                    positive = "1") # Set 1 as positive
confusion_matrix # Output confusion matrix

## Plot Confusion Matrix Heatmap ###
# install.packages("reshape2")
library(reshape2)
# Convert confusion matrix to data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
# Prepare data for plotting, using original counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
# Plotting
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the internal validation cohort") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
library(eoffice) # Load package
topptx(filename = "4Figure_PPT/SVM_confusion_matrix_exvalid2.pptx")


### Gene Validation Set ###
## Validation Set
# Predict using the best parameter configuration on the test set
X_genevalid <- genevalid_lassoed[c(4:ncol(genevalid_lassoed))] # Select all columns except the first and second as independent variables
y_genevalid <- as.factor(genevalid_lassoed$Status) # Select the second column as the dependent variable
genevalid_predictions <- predict(svm_final_model, newdata = X_genevalid)

## Prepare decision values and probabilities for CSV output
# Calculate decision values and probabilities:
genevalid_decisionvalues <- attr(predict(svm_final_model, X_genevalid, decision.values = TRUE), "decision.values") # Output the decision values predicted by SVM, where the absolute value indicates the confidence, with higher absolute values indicating higher confidence in the predicted class
genevalid_pred <- attr(predict(svm_final_model, X_genevalid, probability = TRUE), "probabilities") # Output the predicted probabilities by SVM, between 0 and 1
# Convert decision values and probabilities to data frames
decision_values_df <- as.data.frame(genevalid_decisionvalues)
probabilities_df <- as.data.frame(genevalid_pred)
# Set the column names of the probabilities data frame to class names
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
colnames(decision_values_df) <- paste0("decision_values", colnames(decision_values_df))
# Retain the ID and Status columns
genevalid_svmed <- cbind(ID = genevalid_lassoed$ID, Status = genevalid_lassoed$Status, decision_values_df, probabilities_df)

## Model Evaluation ##
library(caret) # The caret package is powerful and provides detailed evaluation metrics like SPE, PPV, ACC, etc.
# Compute ROC
roc_obj <- roc(response = y_genevalid, predictor = genevalid_pred[, 2])
auc(roc_obj)
# Construct confusion matrix
confusion_matrix <- confusionMatrix(as.factor(ifelse(genevalid_pred[, 2] > event_probability, 1, 0)), # Stratify predicted values, setting threshold to 0.5
                                    y_genevalid <- as.factor(y_genevalid), # Convert outcome variable to factor
                                    positive = "1") # Set 1 as the positive class
confusion_matrix # Output confusion matrix

## Plot Confusion Matrix Heatmap ###
# install.packages("reshape2")
library(reshape2)
# Convert confusion matrix to data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
# Prepare plot data using raw counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
# Plot
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the internal validation cohort") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
library(eoffice) # Load package
topptx(filename = "4Figure_PPT/SVM_confusion_matrix_genevalid.pptx")

### Save Radiomics_Score to CSV ###
# PFV
write.csv(train_svmed, "0-3Radiomics_Score_Result/train_PFV_svmed.csv", row.names = FALSE)
write.csv(valid_svmed, "0-3Radiomics_Score_Result/valid_PFV_svmed.csv", row.names = FALSE)
write.csv(test_svmed, "0-3Radiomics_Score_Result/test_PFV_svmed.csv", row.names = FALSE)
write.csv(exvalid1_svmed, "0-3Radiomics_Score_Result/exvalid1_PFV_svmed.csv", row.names = FALSE)
write.csv(exvalid2_svmed, "0-3Radiomics_Score_Result/exvalid2_PFV_svmed.csv", row.names = FALSE)
write.csv(genevalid_svmed, "0-3Radiomics_Score_Result/genevalid_PFV_svmed.csv", row.names = FALSE)
# Tumor
# write.csv(train_svmed, "0-3Radiomics_Score_Result/train_T_svmed.csv", row.names = FALSE)
# write.csv(valid_svmed, "0-3Radiomics_Score_Result/valid_T_svmed.csv", row.names = FALSE)
# write.csv(test_svmed, "0-3Radiomics_Score_Result/test_T_svmed.csv", row.names = FALSE)
# write.csv(exvalid1_svmed, "0-3Radiomics_Score_Result/exvalid1_T_svmed.csv", row.names = FALSE)
# write.csv(exvalid2_svmed, "0-3Radiomics_Score_Result/exvalid2_T_svmed.csv", row.names = FALSE)
# write.csv(genevalid_svmed, "0-3Radiomics_Score_Result/genevalid_T_svmed.csv", row.names = FALSE)


### Extracting relevant parameters from the decision function ####
# For a linear kernel, the decision function can be written as `f(x) = w^T x + b`. The weight and bias can be extracted using the following code:
# Extract the coefficients (weights) of the support vectors
coefs <- t(svm_final_model$coefs) %*% svm_final_model$SV
# Extract the bias term
bias <- svm_final_model$rho
# Print the weights and bias
print(coefs)
print(bias) #-5.267795
# For a nonlinear kernel (e.g., radial basis kernel `radial`), the decision function is more complex and is typically not easy to express directly in a formula.
# However, you can use the code to make predictions and view the decision function value.

## Model3: XGBoost ####
# XGBoost (eXtreme Gradient Boosting) is a machine learning algorithm based on Gradient Boosting Trees, suitable for classification and regression problems
# install.packages("xgboost")
# install.packages("rBayesianOptimization")
library(xgboost) # Model usage
library(rBayesianOptimization) # Bayesian optimization for parameters
str(train_lassoed)
### Training set ####
# Train the model on the training set
# Define the feature and target variables for the training set
X_train <- train_lassoed[,c(4:ncol(train_lassoed))] # Select all columns as features
y_train <- as.numeric(train_lassoed$Status) # Select the target variable and convert to numeric
# Convert features and target variables to DMatrix format
dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)

# Define the evaluation function (objective function) for Bayesian optimization, which will evaluate the effect of each parameter set
xgb_cv_bayes <- function(max_depth, min_child_weight, subsample, colsample_bytree, gamma, eta, nrounds) {
  params <- list(
    booster = "gbtree",
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = as.integer(max_depth),  # Tree depth
    min_child_weight = as.integer(min_child_weight),  # Minimum weight of leaf nodes
    subsample = subsample,  # Sample ratio
    colsample_bytree = colsample_bytree,  # Feature sampling ratio
    gamma = gamma,  # Minimum loss reduction for tree splitting
    eta = eta  # Learning rate
  )
  
  cv_results <- xgb.cv(
    params = params,
    data = dtrain,
    nrounds = as.integer(nrounds),  # Number of iterations
    nfold = 5,  # 5-fold cross-validation
    verbose = 0,
    early_stopping_rounds = 10,  # Early stopping
    maximize = TRUE  # The goal is to maximize AUC
  )
  # Return the maximum AUC value as the optimization target
  return(list(Score = max(cv_results$evaluation_log$test_auc_mean), Pred = 0))
}
# Set the parameter range for Bayesian optimization
bayes_bound <- list(
  max_depth = c(2L, 10L),
  min_child_weight = c(1L, 10L),
  subsample = c(0.5, 1),
  colsample_bytree = c(0.5, 1),
  gamma = c(0, 5),
  eta = c(0.01, 0.3),
  nrounds = c(50L, 200L)
)
# Run Bayesian optimization
set.seed(123)
bayes_opt <- BayesianOptimization(
  FUN = xgb_cv_bayes,  # Objective function
  bounds = bayes_bound,  # Parameter bounds
  init_points = 10,  # Initial random search iterations
  n_iter = 10,  # Number of optimization iterations
  acq = "ei",  # Acquisition function
  verbose = TRUE
)
# Output the best parameter combination
print(bayes_opt)
# Train the final model using the best parameters found by Bayesian optimization
best_params <- bayes_opt$Best_Par
xgb_final_model <- xgboost(
  data = dtrain,
  max_depth = as.integer(best_params["max_depth"]),
  min_child_weight = as.integer(best_params["min_child_weight"]),
  subsample = best_params["subsample"],
  colsample_bytree = best_params["colsample_bytree"],
  gamma = best_params["gamma"],
  eta = best_params["eta"],
  nrounds = as.integer(best_params["nrounds"]),
  objective = 'binary:logistic',
  eval_metric = 'auc',
  nthread = parallel::detectCores(),
  verbose = 0
)
# Manual parameter tuning
xgb_final_model <- xgboost(
  data = dtrain,
  max_depth = 6,
  min_child_weight = 2,
  # subsample = 0.3，
  # colsample_bytree = best_params["colsample_bytree"],
  gamma = 5,
  eta = 0.1,
  nrounds = 10,
  objective = 'binary:logistic',
  eval_metric = 'auc',
  # nthread = parallel::detectCores(),
  verbose = 0
)

# Make predictions on the training set
train_predictions <- predict(xgb_final_model, newdata = dtrain)

## Output probabilities to CSV ##
# Calculate probabilities:
train_pred <- predict(xgb_final_model, dtrain, type = "response") # Output predicted probabilities between 0-1
# Convert probabilities to a data frame
probabilities_df <- as.data.frame(train_pred)
# Set the column names of the probabilities to match the class name
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
# Retain ID and Status columns
train_XGBoosted <- cbind(ID = train_lassoed$ID, Status = train_lassoed$Status, probabilities_df)


## Model evaluation ##
library(caret) # The caret package provides detailed evaluation metrics like SPE, PPV, ACC, etc.
library(pROC)
# Calculate ROC
roc_obj <- roc(train_lassoed$Status, train_predictions)
auc_ci <- ci.auc(roc_obj)
roc_obj
auc_ci
# Choose the best threshold based on ROC
best_threshold <- coords(roc_obj, "best", ret = "threshold")
best_threshold # View the threshold  
event_probability <- 0.2074544  # PFV0.1264574   T0.2074544
# Recalculate the confusion matrix
# Output confusion matrix based on the best threshold chosen from ROC
confusion_matrix <- confusionMatrix(as.factor(ifelse(train_predictions > event_probability, 1, 0)), 
                                    as.factor(y_train), 
                                    positive = "1") # Set 1 as positive class
confusion_matrix # Output confusion matrix


## Plot confusion matrix heatmap ###
# install.packages("reshape2")
library(reshape2)
# Convert the confusion matrix to a data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
# Prepare data for plotting using raw counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
# Plotting
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the training cohort") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
library(eoffice) # Load package
topptx(filename = "4Figure_PPT/XGB_confusion_matrix_train.pptx")



### Validation set ###
# Predict on the test set using the best parameter configuration
# Define features and target variable for training set
X_valid <- valid_lassoed[, c(4:ncol(valid_lassoed))]  # Select all columns of features as independent variables
y_valid <- as.numeric(valid_lassoed$Status)  # Select the dependent variable and convert to numeric
# Convert features and target variable to DMatrix format
dvalid <- xgb.DMatrix(data = as.matrix(X_valid), label = y_valid)
# Predict on the test set
valid_predictions <- predict(xgb_final_model, newdata = dvalid)

## Output probabilities to CSV ##
# Compute probabilities:
valid_pred <- predict(xgb_final_model, dvalid, type = "response")  # Output predicted probabilities between 0 and 1
# Convert probabilities to a data frame
probabilities_df <- as.data.frame(valid_pred)
# Set column names for probabilities
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
# Retain ID and Status columns
valid_XGBoosted <- cbind(ID = valid_lassoed$ID, Status = valid_lassoed$Status, probabilities_df)

## Model evaluation ##
library(caret)  # caret package is powerful and provides detailed evaluation metrics like SPE, PPV, ACC, etc.
# Compute ROC
roc_obj <- roc(valid_lassoed$Status, valid_predictions)
auc_ci <- ci.auc(roc_obj)
roc_obj
auc_ci
# Choose the best threshold based on ROC and output confusion matrix
confusion_matrix <- confusionMatrix(as.factor(ifelse(valid_predictions > event_probability, 1, 0)), 
                                    as.factor(y_valid), 
                                    positive = "1")  # Set 1 as the positive class
confusion_matrix  # Output confusion matrix

## Plot confusion matrix heatmap ##
# install.packages("reshape2")
library(reshape2)
# Convert confusion matrix to a data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
# Prepare data for plotting using raw counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
# Plotting
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the internal validation cohort") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
library(eoffice)  # Load the package
topptx(filename = "4Figure_PPT/XGB_confusion_matrix_valid.pptx")

### Test set ###
# Predict on the test set using the best parameter configuration
# Define features and target variable for training set
X_test <- test_lassoed[, c(4:ncol(test_lassoed))]  # Select all columns of features as independent variables
y_test <- as.numeric(test_lassoed$Status)  # Select the dependent variable and convert to numeric
# Convert features and target variable to DMatrix format
dtest <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)
# Predict on the test set
test_predictions <- predict(xgb_final_model, newdata = dtest)

## Output probabilities to CSV ##
# Compute probabilities:
test_pred <- predict(xgb_final_model, dtest, type = "response")  # Output predicted probabilities between 0 and 1
# Convert probabilities to a data frame
probabilities_df <- as.data.frame(test_pred)
# Set column names for probabilities
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
# Retain ID and Status columns
test_XGBoosted <- cbind(ID = test_lassoed$ID, Status = test_lassoed$Status, probabilities_df)

## Model evaluation ##
library(caret)  # caret package is powerful and provides detailed evaluation metrics like SPE, PPV, ACC, etc.
# Compute ROC
roc_obj <- roc(test_lassoed$Status, test_predictions)
auc_ci <- ci.auc(roc_obj)
roc_obj
auc_ci
# Choose the best threshold based on ROC and output confusion matrix
confusion_matrix <- confusionMatrix(as.factor(ifelse(test_predictions > event_probability, 1, 0)), 
                                    as.factor(y_test), 
                                    positive = "1")  # Set 1 as the positive class
confusion_matrix  # Output confusion matrix

## Plot confusion matrix heatmap ##
# install.packages("reshape2")
library(reshape2)
# Convert confusion matrix to a data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
# Prepare data for plotting using raw counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
# Plotting
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the test cohort") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
library(eoffice)  # Load the package
topptx(filename = "4Figure_PPT/XGB_confusion_matrix_test.pptx")

### External validation set 1 ###
# Define features and target variable for training set
X_exvalid1 <- exvalid1_lassoed[, c(4:ncol(exvalid1_lassoed))]  # Select all columns of features as independent variables
y_exvalid1 <- as.numeric(exvalid1_lassoed$Status)  # Select the dependent variable and convert to numeric
# Convert features and target variable to DMatrix format
dexvalid1 <- xgb.DMatrix(data = as.matrix(X_exvalid1), label = y_exvalid1)
# Predict on the test set
exvalid1_predictions <- predict(xgb_final_model, newdata = dexvalid1)

## Output probabilities to CSV ##
# Compute probabilities:
exvalid1_pred <- predict(xgb_final_model, dexvalid1, type = "response")  # Output predicted probabilities between 0 and 1
# Convert probabilities to a data frame
probabilities_df <- as.data.frame(exvalid1_pred)
# Set column names for probabilities
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
# Retain ID and Status columns
exvalid1_XGBoosted <- cbind(ID = exvalid1_lassoed$ID, Status = exvalid1_lassoed$Status, probabilities_df)

## Model evaluation ##
library(caret)  # caret package is powerful and provides detailed evaluation metrics like SPE, PPV, ACC, etc.
# Compute ROC
roc_obj <- roc(exvalid1_lassoed$Status, exvalid1_predictions)
auc_ci <- ci.auc(roc_obj)
roc_obj
auc_ci
# Choose the best threshold based on ROC and output confusion matrix
confusion_matrix <- confusionMatrix(as.factor(ifelse(exvalid1_predictions > event_probability, 1, 0)), 
                                    as.factor(y_exvalid1), 
                                    positive = "1")  # Set 1 as the positive class
confusion_matrix  # Output confusion matrix

## Plot confusion matrix heatmap ##
# install.packages("reshape2")
library(reshape2)
# Convert confusion matrix to a data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
# Prepare data for plotting using raw counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
# Plotting
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the external validation cohort 1") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
library(eoffice)  # Load the package
topptx(filename = "4Figure_PPT/XGB_confusion_matrix_exvalid1.pptx")

### External validation set 2 ###
# Define features and target variable for training set
X_exvalid2 <- exvalid2_lassoed[, c(4:ncol(exvalid2_lassoed))]  # Select all columns of features as independent variables
y_exvalid2 <- as.numeric(exvalid2_lassoed$Status)  # Select the dependent variable and convert to numeric
# Convert features and target variable to DMatrix format
dexvalid2 <- xgb.DMatrix(data = as.matrix(X_exvalid2), label = y_exvalid2)
# Predict on the test set
exvalid2_predictions <- predict(xgb_final_model, newdata = dexvalid2)

## Output probabilities to CSV ##
# Compute probabilities:
exvalid2_pred <- predict(xgb_final_model, dexvalid2, type = "response")  # Output predicted probabilities between 0 and 1
# Convert probabilities to a data frame
probabilities_df <- as.data.frame(exvalid2_pred)
# Set column names for probabilities
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))
# Retain ID and Status columns
exvalid2_XGBoosted <- cbind(ID = exvalid2_lassoed$ID, Status = exvalid2_lassoed$Status, probabilities_df)

## Model evaluation ##
library(caret)  # caret package is powerful and provides detailed evaluation metrics like SPE, PPV, ACC, etc.
# Compute ROC
roc_obj <- roc(exvalid2_lassoed$Status, exvalid2_predictions)
auc_ci <- ci.auc(roc_obj)
roc_obj
auc_ci
# Choose the best threshold based on ROC and output confusion matrix
confusion_matrix <- confusionMatrix(as.factor(ifelse(exvalid2_predictions > event_probability, 1, 0)), 
                                    as.factor(y_exvalid2), 
                                    positive = "1")  # Set 1 as the positive class
confusion_matrix  # Output confusion matrix

## Plot confusion matrix heatmap ##
# install.packages("reshape2")
library(reshape2)
# Convert confusion matrix to a data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")
# Prepare data for plotting using raw counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)
# Plotting
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the external validation cohort 2") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")
library(eoffice)  # Load the package
topptx(filename = "4Figure_PPT/XGB_confusion_matrix_exvalid2.pptx")



### Gene validation set ####
# Define training set features and target variable
X_genevalid <- genevalid_lassoed[, c(4:ncol(genevalid_lassoed))] # Select all columns of features as independent variables
y_genevalid <- as.numeric(genevalid_lassoed$Status) # Select the target variable and convert it to numeric

# Convert features and target variable to DMatrix format
dgenevalid <- xgb.DMatrix(data = as.matrix(X_genevalid), label = y_genevalid)

# Make predictions on the validation set
genevalid_predictions <- predict(xgb_final_model, newdata = dgenevalid)

## Output probabilities to CSV ##
# Calculate probabilities:
genevalid_pred <- predict(xgb_final_model, dgenevalid, type = "response") # Output predicted probabilities between 0 and 1

# Convert probabilities to data frame
probabilities_df <- as.data.frame(genevalid_pred)

# Set column names of probabilities to class names
colnames(probabilities_df) <- paste0("Prob_", colnames(probabilities_df))

# Retain ID and Status columns
genevalid_XGBoosted <- cbind(ID = genevalid_lassoed$ID, Status = genevalid_lassoed$Status, probabilities_df)

## Model evaluation ##
library(caret) # The caret package is powerful and provides detailed evaluation metrics like SPE, PPV, ACC, etc.

# Calculate ROC
roc_obj <- roc(genevalid_lassoed$Status, genevalid_predictions)
auc(roc_obj)

# Choose the best threshold based on ROC and output confusion matrix
confusion_matrix <- confusionMatrix(as.factor(ifelse(genevalid_predictions > event_probability, 1, 0)), 
                                    as.factor(y_genevalid), 
                                    positive = "1") # Set 1 as positive class
confusion_matrix # Output confusion matrix


## Plot confusion matrix heatmap ###
# install.packages("reshape2")
library(reshape2)

# Convert confusion matrix to data frame
confusion_matrix_df <- as.data.frame.matrix(confusion_matrix$table)
colnames(confusion_matrix_df) <- c("Negative", "Positive")
rownames(confusion_matrix_df) <- c("Negative", "Positive")

# Prepare data for plotting, using raw counts
draw_data <- confusion_matrix_df
draw_data$real <- rownames(draw_data)
draw_data <- melt(draw_data)

# Plotting
ggplot(draw_data, aes(real, variable, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value)) +  # Display raw counts
  scale_fill_gradient(low = "skyblue", high = "#3575b5") +
  labs(x = "Predictive", y = "Actual", title = "Occult N2 prediction in the genevalid cohort") +
  theme_prism(border = T) +
  theme(panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none")

library(eoffice) # Load the package
topptx(filename = "4Figure_PPT/XGB_confusion_matrix_genevalid.pptx")


### Save Radiomics_Score to CSV ####
# PFV
write.csv(train_XGBoosted, "0-3Radiomics_Score_Result/train_PFV_XGBoosted.csv", row.names = FALSE)
write.csv(valid_XGBoosted, "0-3Radiomics_Score_Result/valid_PFV_XGBoosted.csv", row.names = FALSE)
write.csv(test_XGBoosted, "0-3Radiomics_Score_Result/test_PFV_XGBoosted.csv", row.names = FALSE)
write.csv(exvalid1_XGBoosted, "0-3Radiomics_Score_Result/exvalid1_PFV_XGBoosted.csv", row.names = FALSE)
write.csv(exvalid2_XGBoosted, "0-3Radiomics_Score_Result/exvalid2_PFV_XGBoosted.csv", row.names = FALSE)
write.csv(genevalid_XGBoosted, "0-3Radiomics_Score_Result/genevalid_PFV_XGBoosted.csv", row.names = FALSE)

# Tumor
# write.csv(train_XGBoosted, "0-3Radiomics_Score_Result/train_T_XGBoosted.csv", row.names = FALSE)
# write.csv(valid_XGBoosted, "0-3Radiomics_Score_Result/valid_T_XGBoosted.csv", row.names = FALSE)
# write.csv(test_XGBoosted, "0-3Radiomics_Score_Result/test_T_XGBoosted.csv", row.names = FALSE)
# write.csv(exvalid1_XGBoosted, "0-3Radiomics_Score_Result/exvalid1_T_XGBoosted.csv", row.names = FALSE)
# write.csv(exvalid2_XGBoosted, "0-3Radiomics_Score_Result/exvalid2_T_XGBoosted.csv", row.names = FALSE)
# write.csv(genevalid_XGBoosted, "0-3Radiomics_Score_Result/genevalid_T_XGBoosted.csv", row.names = FALSE)
