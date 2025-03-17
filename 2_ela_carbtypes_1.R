#-------------------------------------------------------------------------------
#  October 24, 2024
#  To derive types of Carb's metabolites signature from LVS:
# 1. wgrain_ela_293.RData
# 2. rgrain_ela_293.RData
# 3. addsugar_ela_293.RData
#-------------------------------------------------------------------------------

library(glmnet)
library(ggplot2)  
library(dplyr)    
library(readr)  

load(file="/udd/n2xwa/carb_met/temp_lvs.RData")

#-------------------------------------------------------------------------------

#                      Elastic net for different carbohydrate types

#-------------------------------------------------------------------------------
names(temp_lvs)
lvs_ela <- temp_lvs[, c(1:294,355,379:388)]
names(lvs_ela)

lvs_ela_type <- lvs_ela
names(lvs_ela_type)

#-------------------------------------------------------------------------------

#                      Elastic net for whole grain: whole_grain_r

#-------------------------------------------------------------------------------

set.seed(1234)
n=floor(nrow(lvs_ela_type)*0.7)
train_ind=sample(seq_len(nrow(lvs_ela_type)), size = n)
train=lvs_ela_type[train_ind,]
test=lvs_ela_type[-train_ind,]

#-------------------------------------------------------------------------------

#                      PREDICTION Training (whole grain)
#                          LOO CV APPROACH

#-------------------------------------------------------------------------------
record_carbo = data.frame(accumulateid=NA,newid=NA,NoMetabs=NA,accumulate_cor=NA)
x=0

for (i in 1:dim(train)[1]) {
  
  Training = cv.glmnet(as.matrix(train[-i,c(2:294)]), train[-i,"whole_grain_r"], nfolds=10, alpha=0.5, family="gaussian")
  lambda_min_10F = Training$lambda.min
  
  cmin = coef(Training, s=lambda_min_10F)
  metabsmin = data.frame(cmin[which(cmin[,1]!=0),])
  names(metabsmin) = "Test.min"
  
  #### apply
  
  if(i==1) {
    train[i,dim(train)[2]+1] = predict(Training, as.matrix(train[i,c(2:294)]), type="response", s=lambda_min_10F) - metabsmin$Test.min[1]
  } else if(i>1) {
    
    train[i,dim(train)[2]] = predict(Training, as.matrix(train[i,c(2:294)]), type="response", s=lambda_min_10F) - metabsmin$Test.min[1]
  }
  
  #### test statistics
  
  x = x+1
  
  record_carbo[x,"accumulateid"] = i
  record_carbo[x,"newid"] = train[i,"id"]
  record_carbo[x,"NoMetabs"] = dim(metabsmin)[1]-1
  
  if(i<3) {
    record_carbo[x,"accumulate_cor"] = NA
  } else if(i>=3) {
    
    tmp = train[1:i,]
    record_carbo[x,"accumulate_cor"] = cor(tmp[,"whole_grain_r"],
                                           tmp[,dim(tmp)[2]], use="complete.obs")
  }
  
}

record_carbo


#-------------------------------------------------------------------------------

#                   PREDICTION TESTING (whole grain)

#-------------------------------------------------------------------------------

record_carbo = data.frame(repNo=NA,NoMetabs=NA,carbo_cor=NA, mse=NA)
met_coef_carbo = data.frame(met=c("intercept",colnames(train[,c(2:294)])))

for (inrep in 1:100) { 
  # traning set on the testing sets
  #cross-validation
  Training_CV = cv.glmnet(as.matrix(train[,c(2:294)]), train[,"whole_grain_r"], nfolds=10, alpha=0.5, family="gaussian")
  lambda_min_10F = Training_CV$lambda.min
  
  Training_M = glmnet(as.matrix(train[,c(2:294)]), train[,"whole_grain_r"], family="gaussian", alpha=0.5) #normal distribution
  cmin = coef(Training_M, s=lambda_min_10F)
  metabsmin = data.frame(cmin[which(cmin[,1]!=0),])
  names(metabsmin) = "Test.min"
  met_coef_carbo[,dim(met_coef_carbo)[2]+1]=cmin[,1]
  names(met_coef_carbo)[dim(met_coef_carbo)[2]] = paste("inrep",inrep,sep='')
  
  #### apply
  #############
  
  test[,dim(test)[2]+1] = as.numeric(predict(Training_M, as.matrix(test[,c(2:294)]), type="response", s=lambda_min_10F)[,1])
  
  test[,dim(test)[2]] = test[,dim(test)[2]] - metabsmin$Test.min[1]
  names(test)[dim(test)[2]] = paste("carbo",inrep,sep='')
  
  
  #### test statistics
  ####################
  
  record_carbo[inrep,"repNo"] = paste("inrep",inrep,sep='')
  record_carbo[inrep,"NoMetabs"] = dim(metabsmin)[1]-1
  record_carbo[inrep,"carbo_cor"] = cor(test[,"whole_grain_r"],test[,dim(test)[2]])
  record_carbo[inrep,"mse"] = min(Training_CV$cvm)
  
}

record_carbo
met_coef_carbo


### find the score for carbo
table(record_carbo$NoMetabs) #  highest 51 times for 52 metabolites
# 18 21 24 36 40 47 52 53 58  #
#  1  1  1  1  1 14 51 26  4  #

sorted_data <- record_carbo %>%
  filter(NoMetabs == 52) %>%
  arrange(mse)

met_carbo=as.character(met_coef_carbo[-1, ] %>% filter(inrep74!=0) %>% pull(met)) #NoMetabs is 49 when inrep3
met_carbo
carbo_coef=as.character(met_coef_carbo[-1, ] %>% filter(inrep74!=0) %>% pull(inrep74)) 
carbo_coef

# create signature file
carbo_score_testing=test[,c("id", "whole_grain_r","carbo74")] #carbo3 should be consistent with (inrep3
colnames(carbo_score_testing)[3]="carbo_score"
carbo_score_training <- train[, c("id", "whole_grain_r", grep("^V", names(train), value = TRUE))]
colnames(carbo_score_training)[3]="carbo_score"

identical(colnames(carbo_score_testing), colnames(carbo_score_training))
carbo_score_testing$set="Testing"
carbo_score_training$set="Training"

signature_carbo=rbind(carbo_score_testing, carbo_score_training) %>% arrange(id)
signature_carbo
write.csv(signature_carbo, "/udd/n2xwa/carb_met/signature_wgrain.csv", na="")

# correlation between metabolite profile score and carbohydrate intake #
cor_pearson_train <- cor(carbo_score_training$carbo_score, carbo_score_training$whole_grain_r, method = "pearson")
cor_spearman_train <- cor(carbo_score_training$carbo_score, carbo_score_training$whole_grain_r, method = "spearman")
cor_pearson_test <- cor(carbo_score_testing$carbo_score, carbo_score_testing$whole_grain_r, method = "pearson")
cor_spearman_test <- cor(carbo_score_testing$carbo_score, carbo_score_testing$whole_grain_r, method = "spearman")

print(cor_pearson_train)
print(cor_spearman_train)
print(cor_pearson_test)
print(cor_spearman_test)


# save ela results to "XXX.RData" and later load
carbo_ela <- list(test, train, record_carbo, met_coef_carbo, met_carbo, carbo_coef, Training_M)
save(carbo_ela, file = "/udd/n2xwa/carb_met/cohort/wgrain_ela_293.RData")


########################    refined grain      ##########################

names(lvs_ela_type)

set.seed(1234)
n=floor(nrow(lvs_ela_type)*0.7)
train_ind=sample(seq_len(nrow(lvs_ela_type)), size = n)
train=lvs_ela_type[train_ind,]
test=lvs_ela_type[-train_ind,]

#-------------------------------------------------------------------------------

#                      PREDICTION Training (refined grain)
#                          LOO CV APPROACH

#-------------------------------------------------------------------------------

record_carbo = data.frame(accumulateid=NA,newid=NA,NoMetabs=NA,accumulate_cor=NA)
x=0

for (i in 1:dim(train)[1]) {
  
  Training = cv.glmnet(as.matrix(train[-i,c(2:294)]), train[-i,"refined_grain_r"], nfolds=10, alpha=0.5, family="gaussian")
  lambda_min_10F = Training$lambda.min
  
  cmin = coef(Training, s=lambda_min_10F)
  metabsmin = data.frame(cmin[which(cmin[,1]!=0),])
  names(metabsmin) = "Test.min"
  
  #### apply
  
  if(i==1) {
    train[i,dim(train)[2]+1] = predict(Training, as.matrix(train[i,c(2:294)]), type="response", s=lambda_min_10F) - metabsmin$Test.min[1]
  } else if(i>1) {
    
    train[i,dim(train)[2]] = predict(Training, as.matrix(train[i,c(2:294)]), type="response", s=lambda_min_10F) - metabsmin$Test.min[1]
  }
  
  #### test statistics
  
  x = x+1
  
  record_carbo[x,"accumulateid"] = i
  record_carbo[x,"newid"] = train[i,"id"]
  record_carbo[x,"NoMetabs"] = dim(metabsmin)[1]-1
  
  if(i<3) {
    record_carbo[x,"accumulate_cor"] = NA
  } else if(i>=3) {
    
    tmp = train[1:i,]
    record_carbo[x,"accumulate_cor"] = cor(tmp[,"refined_grain_r"],
                                           tmp[,dim(tmp)[2]], use="complete.obs")
  }
  
}

record_carbo


#-------------------------------------------------------------------------------

#                   PREDICTION TESTING (refined grain)

#-------------------------------------------------------------------------------


record_carbo = data.frame(repNo=NA,NoMetabs=NA,carbo_cor=NA, mse=NA)
met_coef_carbo = data.frame(met=c("intercept",colnames(train[,c(2:294)])))

for (inrep in 1:100) { 
  # traning set on the testing sets
  #cross-validation
  Training_CV = cv.glmnet(as.matrix(train[,c(2:294)]), train[,"refined_grain_r"], nfolds=10, alpha=0.5, family="gaussian")
  lambda_min_10F = Training_CV$lambda.min
  
  Training_M = glmnet(as.matrix(train[,c(2:294)]), train[,"refined_grain_r"], family="gaussian", alpha=0.5) #normal distribution
  cmin = coef(Training_M, s=lambda_min_10F)
  metabsmin = data.frame(cmin[which(cmin[,1]!=0),])
  names(metabsmin) = "Test.min"
  met_coef_carbo[,dim(met_coef_carbo)[2]+1]=cmin[,1]
  names(met_coef_carbo)[dim(met_coef_carbo)[2]] = paste("inrep",inrep,sep='')
  
  #### apply
  #############
  
  test[,dim(test)[2]+1] = as.numeric(predict(Training_M, as.matrix(test[,c(2:294)]), type="response", s=lambda_min_10F)[,1])
  
  test[,dim(test)[2]] = test[,dim(test)[2]] - metabsmin$Test.min[1]
  names(test)[dim(test)[2]] = paste("carbo",inrep,sep='')
  
  
  #### test statistics
  ####################
  
  record_carbo[inrep,"repNo"] = paste("inrep",inrep,sep='')
  record_carbo[inrep,"NoMetabs"] = dim(metabsmin)[1]-1
  record_carbo[inrep,"carbo_cor"] = cor(test[,"refined_grain_r"],test[,dim(test)[2]])
  record_carbo[inrep,"mse"] = min(Training_CV$cvm)
  
}

record_carbo
met_coef_carbo

### find the score for carbo
table(record_carbo$NoMetabs) #  highest 40 times for 61 metabolites
# 16 38 43 47 50 54 61 64  #
#  2  1  3  5  3 25 40 21 #

sorted_data <- record_carbo %>%
  filter(NoMetabs == 61) %>%
  arrange(mse)


met_carbo=as.character(met_coef_carbo[-1, ] %>% filter(inrep64!=0) %>% pull(met)) #NoMetabs is 49 when inrep3
met_carbo
carbo_coef=as.character(met_coef_carbo[-1, ] %>% filter(inrep64!=0) %>% pull(inrep64)) 
carbo_coef

# create signature file
carbo_score_testing=test[,c("id", "refined_grain_r","carbo64")] #carbo3 should be consistent with (inrep3
colnames(carbo_score_testing)[3]="carbo_score"
carbo_score_training <- train[, c("id", "refined_grain_r", grep("^V", names(train), value = TRUE))]
colnames(carbo_score_training)[3]="carbo_score"

identical(colnames(carbo_score_testing), colnames(carbo_score_training))
carbo_score_testing$set="Testing"
carbo_score_training$set="Training"

signature_carbo=rbind(carbo_score_testing, carbo_score_training) %>% arrange(id)
signature_carbo
write.csv(signature_carbo, "/udd/n2xwa/lvs_met/signature_rgrain.csv", na="")

# correlation between metabolite profile score and carbohydrate intake #
cor_pearson_train <- cor(carbo_score_training$carbo_score, carbo_score_training$refined_grain_r, method = "pearson")
cor_spearman_train <- cor(carbo_score_training$carbo_score, carbo_score_training$refined_grain_r, method = "spearman")
cor_pearson_test <- cor(carbo_score_testing$carbo_score, carbo_score_testing$refined_grain_r, method = "pearson")
cor_spearman_test <- cor(carbo_score_testing$carbo_score, carbo_score_testing$refined_grain_r, method = "spearman")

print(cor_pearson_train)
print(cor_spearman_train)
print(cor_pearson_test)
print(cor_spearman_test)


# save ela results to "XXX.RData" and later load
carbo_ela <- list(test, train, record_carbo, met_coef_carbo, met_carbo, carbo_coef, Training_M)
save(carbo_ela, file = "/udd/n2xwa/lvs_met/cohort/rgrain_ela_293.RData")



#-------------------------------------------------------------------------------

#                      Elastic net for added sugar: addsugar_r

#-------------------------------------------------------------------------------

set.seed(1234)
n=floor(nrow(lvs_ela_type)*0.7)
train_ind=sample(seq_len(nrow(lvs_ela_type)), size = n)
train=lvs_ela_type[train_ind,]
test=lvs_ela_type[-train_ind,]

#-------------------------------------------------------------------------------

#                      PREDICTION Training (addsugar_r)
#                          LOO CV APPROACH

#-------------------------------------------------------------------------------
record_carbo = data.frame(accumulateid=NA,newid=NA,NoMetabs=NA,accumulate_cor=NA)
x=0

for (i in 1:dim(train)[1]) {
  
  Training = cv.glmnet(as.matrix(train[-i,c(2:294)]), train[-i,"addsugar_r"], nfolds=10, alpha=0.5, family="gaussian")
  lambda_min_10F = Training$lambda.min
  
  cmin = coef(Training, s=lambda_min_10F)
  metabsmin = data.frame(cmin[which(cmin[,1]!=0),])
  names(metabsmin) = "Test.min"
  
  #### apply
  
  if(i==1) {
    train[i,dim(train)[2]+1] = predict(Training, as.matrix(train[i,c(2:294)]), type="response", s=lambda_min_10F) - metabsmin$Test.min[1]
  } else if(i>1) {
    
    train[i,dim(train)[2]] = predict(Training, as.matrix(train[i,c(2:294)]), type="response", s=lambda_min_10F) - metabsmin$Test.min[1]
  }
  
  #### test statistics
  
  x = x+1
  
  record_carbo[x,"accumulateid"] = i
  record_carbo[x,"newid"] = train[i,"id"]
  record_carbo[x,"NoMetabs"] = dim(metabsmin)[1]-1
  
  if(i<3) {
    record_carbo[x,"accumulate_cor"] = NA
  } else if(i>=3) {
    
    tmp = train[1:i,]
    record_carbo[x,"accumulate_cor"] = cor(tmp[,"addsugar_r"],
                                           tmp[,dim(tmp)[2]], use="complete.obs")
  }
  
}

record_carbo


#-------------------------------------------------------------------------------

#                   PREDICTION TESTING (addsugar_r)

#-------------------------------------------------------------------------------

record_carbo = data.frame(repNo=NA,NoMetabs=NA,carbo_cor=NA, mse=NA)
met_coef_carbo = data.frame(met=c("intercept",colnames(train[,c(2:294)])))

for (inrep in 1:100) { 
  # traning set on the testing sets
  #cross-validation
  Training_CV = cv.glmnet(as.matrix(train[,c(2:294)]), train[,"addsugar_r"], nfolds=10, alpha=0.5, family="gaussian")
  lambda_min_10F = Training_CV$lambda.min
  
  Training_M = glmnet(as.matrix(train[,c(2:294)]), train[,"addsugar_r"], family="gaussian", alpha=0.5) #normal distribution
  cmin = coef(Training_M, s=lambda_min_10F)
  metabsmin = data.frame(cmin[which(cmin[,1]!=0),])
  names(metabsmin) = "Test.min"
  met_coef_carbo[,dim(met_coef_carbo)[2]+1]=cmin[,1]
  names(met_coef_carbo)[dim(met_coef_carbo)[2]] = paste("inrep",inrep,sep='')
  
  #### apply
  #############
  
  test[,dim(test)[2]+1] = as.numeric(predict(Training_M, as.matrix(test[,c(2:294)]), type="response", s=lambda_min_10F)[,1])
  
  test[,dim(test)[2]] = test[,dim(test)[2]] - metabsmin$Test.min[1]
  names(test)[dim(test)[2]] = paste("carbo",inrep,sep='')
  
  
  #### test statistics
  ####################
  
  record_carbo[inrep,"repNo"] = paste("inrep",inrep,sep='')
  record_carbo[inrep,"NoMetabs"] = dim(metabsmin)[1]-1
  record_carbo[inrep,"carbo_cor"] = cor(test[,"addsugar_r"],test[,dim(test)[2]])
  record_carbo[inrep,"mse"] = min(Training_CV$cvm)
  
}

record_carbo
met_coef_carbo


### find the score for carbo
table(record_carbo$NoMetabs) #  highest 42 times for 68 metabolites
# 56 62 65 68 73 82    #
#  1  5 14 42 35  3      #
sorted_data <- record_carbo %>%
  filter(NoMetabs == 68) %>%
  arrange(mse)

met_carbo=as.character(met_coef_carbo[-1, ] %>% filter(inrep57!=0) %>% pull(met)) #NoMetabs is 49 when inrep3
met_carbo
carbo_coef=as.character(met_coef_carbo[-1, ] %>% filter(inrep57!=0) %>% pull(inrep57)) 
carbo_coef

# create signature file
carbo_score_testing=test[,c("id", "addsugar_r","carbo57")] #carbo3 should be consistent with (inrep3
colnames(carbo_score_testing)[3]="carbo_score"
carbo_score_training <- train[, c("id", "addsugar_r", grep("^V", names(train), value = TRUE))]
colnames(carbo_score_training)[3]="carbo_score"

identical(colnames(carbo_score_testing), colnames(carbo_score_training))
carbo_score_testing$set="Testing"
carbo_score_training$set="Training"

signature_carbo=rbind(carbo_score_testing, carbo_score_training) %>% arrange(id)
signature_carbo


# correlation between metabolite profile score and carbohydrate intake #

cor_pearson_train <- cor(carbo_score_training$carbo_score, carbo_score_training$addsugar_r, method = "pearson")
cor_spearman_train <- cor(carbo_score_training$carbo_score, carbo_score_training$addsugar_r, method = "spearman")
cor_pearson_test <- cor(carbo_score_testing$carbo_score, carbo_score_testing$addsugar_r, method = "pearson")
cor_spearman_test <- cor(carbo_score_testing$carbo_score, carbo_score_testing$addsugar_r, method = "spearman")

print(cor_pearson_train)
print(cor_spearman_train)
print(cor_pearson_test)
print(cor_spearman_test)


# save ela results to "XXX.RData" and later load
carbo_ela <- list(test, train, record_carbo, met_coef_carbo, met_carbo, carbo_coef, Training_M)
save(carbo_ela, file = "/udd/n2xwa/carb_met/cohort/addsugar_ela_293.RData")

# Bar figure: Carbo Coef by selected metabolites 

ela_3 <- data.frame(met_carbo, carbo_coef)
ela_3<-ela_3 %>% mutate(met=met_carbo)

re <- read.csv("/udd/n2xwa/carb_met/results_lvs_20241205.csv")
ela_4 <- left_join(ela_3, re[, c("met", "metabolite_name","class_metabolon")], by = "met")
ela_4$carbo_coef = as.numeric(ela_4$carbo_coef)
ela_4$class_metabolon = as.character(ela_4$class_metabolon)
ela_4$metabolite_name = as.character(ela_4$metabolite_name)

ela_4_sorted <- ela_4[order(ela_4$class_metabolon), ] 
ela_4_sorted$metabolite_name <- factor(ela_4_sorted$metabolite_name, levels = unique(ela_4_sorted$metabolite_name))

pp <- ggplot(ela_4_sorted, aes(x = metabolite_name, y = carbo_coef, fill = class_metabolon, group = class_metabolon)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Carbo Coef by selected metabolites", x = "Metabolite Name", y = "Carbo Coef") +
  scale_y_continuous(expand = c(0, 0), limits = c(-1, 1) * max(abs(ela_4_sorted$carbo_coef))) +
  theme_minimal() +
  theme(axis.text.x = element_text(hjust = 1)) +
  theme(axis.text.y = element_text(hjust = 1)) +
  coord_flip() 

ggsave("/udd/n2xwa/carb_met/elabar_addsugar.png", plot = pp, width = 8, height = 8, dpi = 1200)
