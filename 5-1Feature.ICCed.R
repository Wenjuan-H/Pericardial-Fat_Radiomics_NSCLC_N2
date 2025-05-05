#ICC####
#install.packages("psych")
#install.packages("irr")
library(irr)
library(psych)


#PFV dataset####
# Importing datasets
PFV_icc1 <- read.csv("0-1Radiomics_feature_Data/ICC-PFV-radiomics1.csv", header=TRUE)
PFV_icc2 <- read.csv("0-1Radiomics_feature_Data/ICC-PFV-radiomics2.csv", header=TRUE)

# Checking data structure
str(PFV_icc1)
str(PFV_icc2)

# Data preparation
n <- dim(PFV_icc1)[1]
df_intra <- rbind(PFV_icc1, PFV_icc2)

# Computing ICC for intra-consistency
icc_all_intra <- apply(df_intra, 2, function(x) ICC(x = data.frame(x[1:n], x[(n+1):(2*n)]), lmer = FALSE)$results[1,2])

# Saving ICC results to CSV
df_icc_all_intra <- data.frame(feature = names(df_intra), icc = icc_all_intra, row.names = 1:ncol(df_intra))

# Filtering datasets based on ICC results
trainall_icc <- trainall[, columns_intra, drop = FALSE]
exvalid1_icc <- exvalid1[, columns_intra, drop = FALSE]
exvalid2_icc <- exvalid2[, columns_intra, drop = FALSE]
genevalid_icc <- genevalid[, columns_intra, drop = FALSE]
trainOS_icc <- trainOS[, columns_intra, drop = FALSE]

# Merging ID and Status columns with selected features
trainall_icced <- cbind(ID = trainall$ID, Status = trainall$StatusN2, trainall_icc)
exvalid1_icced <- cbind(ID = exvalid1$ID, Status = exvalid1$StatusN2, exvalid1_icc)
exvalid2_icced <- cbind(ID = exvalid2$ID, Status = exvalid2$StatusN2, exvalid2_icc)
genevalid_icced <- cbind(ID = genevalid$ID, Status = genevalid$StatusN2, genevalid_icc)
trainOS_icced <- cbind(ID = trainOS$ID, Status = trainOS$StatusN2, trainOS_icc)

# Saving datasets
write.csv(trainall_icced, "0-2Radiomics_feature_Slected/三院PFVradiomics_icced.csv")
write.csv(exvalid1_icced, "0-2Radiomics_feature_Slected/二院PFVradiomics_icced.csv")
write.csv(exvalid2_icced, "0-2Radiomics_feature_Slected/一院PFVradiomics_icced.csv")
write.csv(genevalid_icced, "0-2Radiomics_feature_Slected/genePFVradiomics_icced.csv")
write.csv(trainOS_icced, "0-2Radiomics_feature_Slected/OS三院PFVradiomics_icced.csv")


# Tumor dataset####
# Importing datasets
T_icc1 <- read.csv("1Clinial_Data/ICC-T-radiomics1.csv", header=TRUE)
T_icc2 <- read.csv("1Clinial_Data/ICC-T-radiomics2.csv", header=TRUE)

# Checking data structure
str(T_icc1)
str(T_icc2)

# Data preparation
n <- dim(T_icc1)[1]
df_intra <- rbind(T_icc1, T_icc2)

# Computing ICC for intra-consistency
icc_all_intra <- apply(df_intra, 2, function(x) ICC(x = data.frame(x[1:n], x[(n+1):(2*n)]), lmer = FALSE)$results[1,2])

# Saving ICC results to CSV
df_icc_all_intra <- data.frame(feature = names(df_intra), icc = icc_all_intra, row.names = 1:ncol(df_intra))

# Filtering datasets based on ICC results
trainall_icc <- trainall[, columns_intra, drop = FALSE]
exvalid1_icc <- exvalid1[, columns_intra, drop = FALSE]
exvalid2_icc <- exvalid2[, columns_intra, drop = FALSE]
genevalid_icc <- genevalid[, columns_intra, drop = FALSE]
trainOS_icc <- trainOS[, columns_intra, drop = FALSE]

# Merging ID and Status columns with selected features
trainall_icced <- cbind(ID = trainall$ID, Status = trainall$StatusN2, trainall_icc)
exvalid1_icced <- cbind(ID = exvalid1$ID, Status = exvalid1$StatusN2, exvalid1_icc)
exvalid2_icced <- cbind(ID = exvalid2$ID, Status = exvalid2$StatusN2, exvalid2_icc)
genevalid_icced <- cbind(ID = genevalid$ID, Status = genevalid$StatusN2, genevalid_icc)
trainOS_icced <- cbind(ID = trainOS$ID, Status = trainOS$StatusN2, trainOS_icc)

# Saving datasets
write.csv(trainall_icced, "0-2Radiomics_feature_Slected/三院Tradiomics_icced.csv")
write.csv(exvalid1_icced, "0-2Radiomics_feature_Slected/二院Tradiomics_icced.csv")
write.csv(exvalid2_icced, "0-2Radiomics_feature_Slected/一院Tradiomics_icced.csv")
write.csv(genevalid_icced, "0-2Radiomics_feature_Slected/geneTradiomics_icced.csv")
write.csv(trainOS_icced, "0-2Radiomics_feature_Slected/OS三院Tradiomics_icced.csv")