###Splitting sets + Z-score
#Importing datasets####
#PFV
#trainall<-read.csv("0-2Radiomics_feature_Slected/三院PFVradiomics_icced.csv",header=TRUE)
#exvalid1<-read.csv("0-2Radiomics_feature_Slected/二院PFVradiomics_icced.csv",header=TRUE)
#exvalid2<-read.csv("0-2Radiomics_feature_Slected/一院PFVradiomics_icced.csv",header=TRUE)
#genevalid<-read.csv("0-2Radiomics_feature_Slected/genePFVradiomics_icced.csv",header=TRUE)
#OS<-read.csv("0-2Radiomics_feature_Slected/OS三院PFVradiomics_icced.csv",header=TRUE)
#Tumor
trainall<-read.csv("0-2Radiomics_feature_Slected/三院Tradiomics_icced.csv",header=TRUE)
exvalid1<-read.csv("0-2Radiomics_feature_Slected/二院Tradiomics_icced.csv",header=TRUE)
exvalid2<-read.csv("0-2Radiomics_feature_Slected/一院Tradiomics_icced.csv",header=TRUE)
genevalid<-read.csv("0-2Radiomics_feature_Slected/geneTradiomics_icced.csv",header=TRUE)
OS<-read.csv("0-2Radiomics_feature_Slected/OS三院Tradiomics_icced.csv",header=TRUE)

#Splitting training and validation sets####
library(splitTools)
library(ranger)
#Set random seed and splitting ratios
set.seed(3781)  #3781
Datasplit<-partition(trainall$ID,p = c(train = 0.6, valid = 0.2, test = 0.2)) 
str(Datasplit)
#Rename the resulting datasets
train<-trainall[Datasplit$train,]
valid<-trainall[Datasplit$valid,]
test<-trainall[Datasplit$test,]
#save sets
#PFV
#write.csv(train,"0-2Radiomics_feature_Slected/trainPFVradiomics_icced.csv")    #训练集文件保存
#write.csv(valid,"0-2Radiomics_feature_Slected/validPFVradiomics_icced.csv")    #验证集文件保存
#write.csv(test,"0-2Radiomics_feature_Slected/testPFVradiomics_icced.csv")    #验证集文件保存
#Tumor
write.csv(train,"0-2Radiomics_feature_Slected/trainTradiomics_icced.csv")    #训练集文件保存
write.csv(valid,"0-2Radiomics_feature_Slected/validTradiomics_icced.csv")    #验证集文件保存
write.csv(test,"0-2Radiomics_feature_Slected/testTradiomics_icced.csv")    #验证集文件保存


#Z-score####
#install.packages("caret")
#Z-score normalization
#Z-score normalization is a method to scale the data to a standard normal distribution.
#It preserves the original distribution of the data and transforms it to have a mean of 0 and standard deviation of 1.
#The mathematical formula for Z-score normalization is:
#z = (x-μ)/σ
library(caret)
### Extract parameters from the training set. method = c("center", "scale") is for standardization; method = c("range") is for normalization.
# Use the preProcess() function for normalization
#normal_para<- preProcess(x=train[, 4:829], method = c("center","scale"))  #PFV
normal_para<- preProcess(x=train[, 4:792], method = c("center","scale")) #Tumor
# Apply the transformation to the data
train_normal <- predict(object = normal_para, newdata =train[, 4:829])
valid_normal <-predict(object= normal_para, newdata =valid[, 4:829])
test_normal<-predict(object=normal_para, newdata =test[, 4:829])
exvalid1_normal <-predict(object= normal_para, newdata =exvalid1[, 4:829])
exvalid2_normal <-predict(object= normal_para, newdata =exvalid2[, 4:829])
genevalid_normal <-predict(object= normal_para, newdata =genevalid[, 4:829])
OS_normal <- predict(object = normal_para, newdata =OS[, 4:829])
# Keep ID and Status columns
train_normalized <- cbind(ID = train$ID, Status = train$Status,train_normal)
valid_normalized <- cbind(ID = valid$ID, Status = valid$Status,valid_normal)
test_normalized <- cbind(ID = test$ID, Status = test$Status,test_normal)
exvalid1_normalized <- cbind(ID = exvalid1$ID, Status = exvalid1$Status,exvalid1_normal)
exvalid2_normalized <- cbind(ID = exvalid2$ID, Status = exvalid2$Status,exvalid2_normal)
genevalid_normalized <- cbind(ID = genevalid$ID, Status = genevalid$Status,genevalid_normal)
OS_normalized <- cbind(ID = OS$ID, Status = OS$Status,OS_normal)
# Save the data
#PFV
#write.csv(train_normalized,"0-2Radiomics_feature_Slected/train_PFV_normalized.csv")    #训练集文件保存
#write.csv(valid_normalized,"0-2Radiomics_feature_Slected/valid_PFV_normalized.csv")    #验证集文件保存
#write.csv(test_normalized,"0-2Radiomics_feature_Slected/test_PFV_normalized.csv")    #验证集文件保存
#write.csv(exvalid1_normalized,"0-2Radiomics_feature_Slected/exvalid1_PFV_normalized.csv")    #验证集文件保存
#write.csv(exvalid2_normalized,"0-2Radiomics_feature_Slected/exvalid2_PFV_normalized.csv")    #验证集文件保存
#write.csv(genevalid_normalized,"0-2Radiomics_feature_Slected/genevalid_PFV_normalized.csv")    #验证集文件保存
#write.csv(OS_normalized,"0-2Radiomics_feature_Slected/OS_PFV_normalized.csv")    #训练集文件保存
#Tumor
write.csv(train_normalized,"0-2Radiomics_feature_Slected/train_T_normalized.csv")    #训练集文件保存
write.csv(valid_normalized,"0-2Radiomics_feature_Slected/valid_T_normalized.csv")    #验证集文件保存
write.csv(test_normalized,"0-2Radiomics_feature_Slected/test_T_normalized.csv")    #验证集文件保存
write.csv(exvalid1_normalized,"0-2Radiomics_feature_Slected/exvalid1_T_normalized.csv")    #验证集文件保存
write.csv(exvalid2_normalized,"0-2Radiomics_feature_Slected/exvalid2_T_normalized.csv")    #验证集文件保存
write.csv(genevalid_normalized,"0-2Radiomics_feature_Slected/genevalid_T_normalized.csv")    #验证集文件保存
write.csv(OS_normalized,"0-2Radiomics_feature_Slected/OS_T_normalized.csv")    #训练集文件保存
