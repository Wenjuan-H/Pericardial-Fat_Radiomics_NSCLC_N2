#Radiomics fature selection(Person+LASSO)####
#Load data####
#PFV
#train_normalized<-read.csv("0-2Radiomics_feature_Slected/train_PFV_normalized.csv",header=TRUE)
#valid_normalized<-read.csv("0-2Radiomics_feature_Slected/valid_PFV_normalized.csv",header=TRUE)
#test_normalized<-read.csv("F:0-2Radiomics_feature_Slected/test_PFV_normalized.csv",header=TRUE)
#exvalid1_normalized<-read.csv("0-2Radiomics_feature_Slected/exvalid1_PFV_normalized.csv",header=TRUE)
#exvalid2_normalized<-read.csv("0-2Radiomics_feature_Slected/exvalid2_PFV_normalized.csv",header=TRUE)
#genevalid_normalized<-read.csv("0-2Radiomics_feature_Slected/genevalid_PFV_normalized.csv",header=TRUE)
#OS_normalized<-read.csv("0-2Radiomics_feature_Slected/OS_PFV_normalized.csv",header=TRUE)
#Tumor
train_normalized<-read.csv("0-2Radiomics_feature_Slected/train_T_normalized.csv",header=TRUE)
valid_normalized<-read.csv("0-2Radiomics_feature_Slected/valid_T_normalized.csv",header=TRUE)
test_normalized<-read.csv("F:0-2Radiomics_feature_Slected/test_T_normalized.csv",header=TRUE)
exvalid1_normalized<-read.csv("0-2Radiomics_feature_Slected/exvalid1_T_normalized.csv",header=TRUE)
exvalid2_normalized<-read.csv("0-2Radiomics_feature_Slected/exvalid2_T_normalized.csv",header=TRUE)
genevalid_normalized<-read.csv("0-2Radiomics_feature_Slected/genevalid_T_normalized.csv",header=TRUE)
OS_normalized<-read.csv("0-2Radiomics_feature_Slected/OS_T_normalized.csv",header=TRUE)

#Feature redundancy method 1 (Pearson)####
#Remove redundancy based on Pearson correlation
#Check the dimension (i.e., number of features) of training and validation sets
dim(train_normalized)
dim(valid_normalized)
dim(test_normalized)
##feature selection
#Perform normality test on features from column 4 to the last column
norm_result <- apply(train_normalized[, 4:ncol(train_normalized)], 2, function(x) {
  if (length(na.omit(x)) < 3) {
    return(NA)  # Return NA if valid data points are fewer than 3
  }
  return(shapiro.test(x)$p.value)
})
# Select features that meet normality (p > 0.05)
norm_feature <- train_normalized[, 4:ncol(train_normalized)][, norm_result >= 0.05, drop = FALSE]

# Calculate correlation coefficient: r
# Use Pearson for normally distributed features
cor_nor <- cor(norm_feature, method = "pearson")
# Use Spearman for non-normally distributed features
cor_all <- cor(train_normalized[, 4:ncol(train_normalized)], method = "spearman")
# Replace Spearman values with Pearson values for normal features
num_nor <- dim(cor_nor)[1]
cor_all[1:num_nor, 1:num_nor] <- cor_nor
# Set upper triangle and diagonal to 0
cor_all[upper.tri(cor_all)] <- 0
diag(cor_all) <- 0
# Remove features with correlation coefficient greater than 0.9
reduced_features <- train_normalized[, 4:ncol(train_normalized)][, !apply(cor_all, 2, function(x) any(abs(x) > 0.9))]
# Keep original columns 1–3 and combine with reduced features
data_reduce <- cbind(train_normalized[, 1:3], reduced_features)
# Check dimension after reduction
dim(data_reduce)


#Lasso feature selection####
#install.packages("glmnet") # for lasso visualization
#install.packages("broom")  # for improved plotting
library(glmnet)
library(broom)  
#Matrix conversion (set predictors and response for lasso)
#x <- as.matrix(train_normalized[,c(4:ncol(train_normalized))])  #  predictors
#y <- as.matrix(train_normalized$Status)   # response
x <- as.matrix(data_reduce[,c(4:ncol(data_reduce))])  # predictors
y <- as.matrix(data_reduce$Status)   # response
#x <- as.matrix(train_mRMRed[,c(3:ncol(train_mRMRed))])  # predictors
#y <- as.matrix(train_mRMRed$Status)   # response

#Plot regression coefficient paths
# Perform LASSO using glmnet
lasso_model <- glmnet(x, y,
                      family="binomial",
                      nlambda=100, # 100 points along the regularization path
                      alpha=1)  # alpha=1 for LASSO; alpha=0 for Ridge
#family parameter sets the type of regression model:
#family="gaussian": univariate continuous response
#family="mgaussian": multivariate continuous response
#family="poisson": non-negative count data
#family="binomial": binary response
#family="multinomial": categorical response
print(lasso_model)
#View(print(lasso_model))#Df: degrees of freedom, %Dev: deviance explained, Lambda
#Check number of iterations to find optimal lambda
NROW(lasso_model$lambda)
#Number of iterations performed
#Stops at 92nd lambda, optimal solution found, no further improvement
#Select optimal lambda value
lasso_model$lambda[100]#Optimal lambda PFV100:6.708869e-06 ;Tumor100:9.21998e-06
#Plot coefficient paths
#windows()
plot(lasso_model, xvar = "lambda", label = TRUE)
#X-axis: log of lambdas, Y-axis: coefficients, some shrink to 0 as lambda increases

#Export editable figure to PPT format
library(eoffice) # load package
topptx(filename = "4Figure_PPT/Lasso_PFV_feature_slected_1.pptx")
#topptx(filename = "4Figure_PPT/Lasso_T_feature_slected_1.pptx")


# Plot cross-validation error curve
#Set seed
set.seed(1234) #PFV1234, Tumor1234
#Build cross-validation model
cvmodel <-  cv.glmnet(x, y,family="binomial", 
                      nlambda=100, #PFV100, Tumor
                      alpha=1)#default nfolds = 10            
plot(cvmodel)
print(cvmodel) 
#PFV:1se Lambda0.024110  Index12 Measure0.5495 SE0.04774 Nonzero7
#Tumor:1se Lambda0.04380  Index9 Measure0.5350  SE0.04571 Nonzero4
#Check Nonzero column for number of variables at each λ
#X-axis: log λ, Y-axis: deviance (lower = better fit)
#Top: number of variables left at each λ
#Left dashed line: λ min, best fit
#Right dashed line: λ 1se, 1 SE to the right of min
#λ 1se is usually preferred clinically for simpler model

###Export editable figure to PPT format
library(eoffice) # load package
topptx(filename = "4Figure_PPT/Lasso_PFV_feature_slected_2.pptx")
#topptx(filename = "4Figure_PPT/Lasso_T_feature_slected_2.pptx")


##Extract features with non-zero coefficients and their values
coefPara <-coef(object=cvmodel,s="lambda.1se") #Choose either "lambda.min" or "lambda.1se"
lasso_values <- as.data.frame( which(coefPara !=0, arr.ind = T))
lasso_names <-rownames(lasso_values)[-1]
Lasso_coef <- data.frame(Feature =rownames(lasso_values),
                         Coef = coefPara[which(coefPara !=0,arr.ind = T)])
Lasso_coef

#Organize the reduced feature sets
train_set_lasso<- data.frame(x)[lasso_names]
valid_set_lasso<-valid_normalized[names(train_set_lasso)] 
test_set_lasso<-test_normalized[names(train_set_lasso)] 
exvalid1_set_lasso<-exvalid1_normalized[names(train_set_lasso)] 
exvalid2_set_lasso<-exvalid2_normalized[names(train_set_lasso)] 
genevalid_set_lasso<-genevalid_normalized[names(train_set_lasso)] 
OS_set_lasso<-OS_normalized[names(train_set_lasso)]
# save ID and Status 
train_lassoed <- cbind(ID = train_normalized$ID, Status = train_normalized$Status,train_set_lasso)
valid_lassoed <- cbind(ID = valid_normalized$ID, Status = valid_normalized$Status,valid_set_lasso)
test_lassoed <- cbind(ID = test_normalized$ID, Status = test_normalized$Status,test_set_lasso)
exvalid1_lassoed <- cbind(ID = exvalid1_normalized$ID, Status = exvalid1_normalized$Status,exvalid1_set_lasso)
exvalid2_lassoed <- cbind(ID = exvalid2_normalized$ID, Status = exvalid2_normalized$Status,exvalid2_set_lasso)
genevalid_lassoed <- cbind(ID = genevalid_normalized$ID, Status =genevalid_normalized$Status,genevalid_set_lasso)
OS_lassoed <- cbind(ID = OS_normalized$ID, Status = OS_normalized$Status,OS_set_lasso)
#save sets
#PFV
#write.csv(train_lassoed,"0-2Radiomics_feature_Slected/train_PFV_lassoed.csv")    #训练集文件保存
#write.csv(valid_lassoed,"0-2Radiomics_feature_Slected/valid_PFV_lassoed.csv")    #验证集文件保存
#write.csv(test_lassoed,"0-2Radiomics_feature_Slected/test_PFV_lassoed.csv")    #验证集文件保存
#write.csv(exvalid1_lassoed,"0-2Radiomics_feature_Slected/exvalid1_PFV_lassoed.csv")    #验证集文件保存
#write.csv(exvalid2_lassoed,"0-2Radiomics_feature_Slected/exvalid2_PFV_lassoed.csv")    #验证集文件保存
#write.csv(genevalid_lassoed,"0-2Radiomics_feature_Slected/genevalid_PFV_lassoed.csv")    #验证集文件保存
#write.csv(OS_lassoed,"0-2Radiomics_feature_Slected/OS_PFV_lassoed.csv")    #验证集文件保存
#Tumor
write.csv(train_lassoed,"0-2Radiomics_feature_Slected/train_T_lassoed.csv")    #训练集文件保存
write.csv(valid_lassoed,"0-2Radiomics_feature_Slected/valid_T_lassoed.csv")    #验证集文件保存
write.csv(test_lassoed,"0-2Radiomics_feature_Slected/test_T_lassoed.csv")    #验证集文件保存
write.csv(exvalid1_lassoed,"0-2Radiomics_feature_Slected/exvalid1_T_lassoed.csv")    #验证集文件保存
write.csv(exvalid2_lassoed,"0-2Radiomics_feature_Slected/exvalid2_T_lassoed.csv")    #验证集文件保存
write.csv(genevalid_lassoed,"0-2Radiomics_feature_Slected/genevalid_T_lassoed.csv")    #验证集文件保存
write.csv(OS_lassoed,"0-2Radiomics_feature_Slected/OS_T_lassoed.csv")    #验证集文件保存

