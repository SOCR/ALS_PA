---
title: "ALS project"
author: "Ming Tang, Chao Gao, Ivo Dinov"
date: "October 14, 2016"
output: html_document
---

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


### set directory######
```{r}
setwd("/User/project/ALS")
```

### Request more heap memory ###
```{r}
options( java.parameters = "-Xmx12g" )
```

### Check Current Java Version ###
```{r, results='hide'}
Sys.getenv("JAVA_HOME")
```


### load packages###
```{r, message=FALSE, warning=FALSE}
library(dplyr)
library(plyr)
library(mi)
library(reshape)
library(reshape2)
library(ggplot2)
library(magrittr)
library(psych)
#library(bartMachine)
library(SuperLearner)
library(caret)
library(e1071)
```

### Personalized functions ###
```{r}
##Missing plot
ggplot_missing <- function(x){
  
  x %>% 
    is.na %>%
    melt %>%
    ggplot(data = .,aes(x = Var2,y = Var1)) +
    geom_raster(aes(fill = value)) +
    scale_fill_grey(name = "",labels = c("Present","Missing")) +
    theme_minimal() + 
    theme(axis.text.x  = element_text(angle=45, vjust=0.5)) + 
    coord_flip()+
    labs(x = "Variables in Dataset",y = "Rows / observations")
}
##Replace missing with the median of the column
f<-function(x){
  x[is.na(x)] =median(x, na.rm=TRUE) #convert the item with NA to median value from the column
  x #display the column
}
```
###Load datasets###
#### Load raw PRO-ACT training dataset (3103974 obs. of 6 variables) ####
```{r}
df <- read.table("Data/all_forms_PROACT.txt", header=TRUE, sep="|", na.strings=c("NA","","NaN"), fill = TRUE, quote = "", stringsAsFactors = FALSE, comment.char="")
```

#### Load PRO-ACT training ALSFRS_slope dataset (2424 obs.) ####
```{r}
df_slope <- read.table("Data/ALSFRS_slope_PROACT.txt", header=TRUE, sep="|", na.strings=c("NA","","NaN"), fill = TRUE, quote = "", stringsAsFactors = FALSE, comment.char="")
#remove study subjects with NA in ALSFRS_slope
df_slope <- df_slope[!is.na(df_slope$ALSFRS_slope), ]
```

#### Load raw PRO-ACT testing dataset (133724 obs. of 6 variables)
```{r}
df_v <- read.table("Data/all_forms_validate_spike.txt", header=TRUE, sep="|", na.strings=c("NA","","NaN"), fill = TRUE, quote = "", stringsAsFactors = FALSE, comment.char="")
```

#### Load PRO-ACT testing ALSFRS_slope dataset (101 obs.)###
```{r}
df_v_slope <- read.table("Data/ALSFRS_slope_validate_spike.txt",header=TRUE,sep="|", na.strings=c("NA","","NaN"), fill = TRUE, quote = "", stringsAsFactors = FALSE, comment.char="")
```
### Generate tile plot to present missing pattern of raw PRO-ACT training dataset ###
```{r,fig.keep='all', echo=FALSE, eval=FALSE}
#this chuch of code generates a tile plot containing 2424 variables and 103 variables
#load dataset
df0 <- df

#list of features
feature_list <- c("Age", "Gender", "Race", "onset_delta", "onset_site", "Q1_Speech", "Q2_Salivation","Q3_Swallowing", "mouth","Q4_Handwriting", "Q5a_Cutting_without_Gastrostomy", "Q5b_Cutting_with_Gastrostomy", "Q5_Cutting", "hands",
"Q6_Dressing_and_Hygiene", "Q7_Turning_in_Bed","trunk", "Q8_Walking", "Q9_Climbing_Stairs", "leg", "Q10_Respiratory","R1_Dyspnea","respiratory","R2_Orthopnea", "R3_Respiratory_Insufficiency", "respiratory_R","ALSFRS_Total", "ALSFRS_R_Total","fvc_normal","fvc_percent","fvc", "svc", "svc_percent","bp_diastolic","bp_systolic","height","weight", "BMI","temperature","pulse", "respiratory_rate","if_use_Riluzole", "Urine Ph","Urine Specific Gravity", "Urine Protein", "Urine WBCs","Urine Glucose", "Urine blood", "Urine Casts","Urine Ketones","Urine Appearance", "Urine Color", "Urine Clarity", 
"Urine RBCs","Albumin", "Protein", "Sodium", "Potassium","Bicarbonate", "Chloride","Anion Gap", "Magnesium",
"Blood Urea Nitrogen (BUN)","Uric Acid", "Creatinine","Alkaline Phosphatase","ALT(SGPT)", "Gamma-glutamyltransferase", "AST(SGOT)","Bilirubin (Total)", "White Blood Cell (WBC)", "Neutrophils","Band Neutrophils","Lymphocytes", "Monocytes",
"Eosinophils","Basophils","RBC Morphology","Red Blood Cells (RBC)","Hemoglobin", "Hematocrit","Mean Corpuscular Hemoglobin Concentration","Mean Corpuscular Volume","Mean Corpuscular Hemoglobin","Platelets","Transferrin","CK","Triglycerides","Total Cholesterol","Lactate Dehydrogenase","Glucose","HbA1c (Glycated Hemoglobin)", "Calcium", "Phosphorus", "Hepatitis A Antibody","Hepatitis B Antibody", "Hepatitis C Antibody","Hepatitis B Antigen","IMMUNOGLOBULIN A","IMMUNOGLOBULIN G","IMMUNOGLOBULIN M", "GAMMA-GLOBULIN", "TSH", "Free T3", "Free T4", "Beta HCG", "International Normalized Ratio (clotting)", "Amylase", "Salivary Amylase", "Pancreatic Amylase")

#only display information of study subjects with valid ALSFRS_slope
id_lst <- unique(df_slope$SubjectID)
df_counts <- data.frame(id_lst)
df_temp <- df0[df0$SubjectID %in% df_slope$SubjectID, ]
colnames(df_counts) <- "SubjectID"
for (f in feature_list){
    df_feature <- df_temp[df_temp$feature_name == f, ]
    counts <- as.data.frame(table(df_feature$SubjectID))
    colnames(counts) <- c("SubjectID", f)
    df_counts <- merge(df_counts, counts, by = "SubjectID", all = TRUE)
    
}

df_counts[is.na(df_counts)] <- 0 
aaa <- df_counts
aaa$SubjectID <- as.factor(aaa$SubjectID)

breakpoints <- levels(aaa$SubjectID)[seq(1, length(levels(aaa$SubjectID)), 20)]
aaa.m <- melt(aaa)
colnames(aaa.m)[3] <- "Counts"

png("plot2.png", width=9000, height=6000, res = 300)


ppp <- ggplot(aaa.m, aes(variable, SubjectID)) + geom_tile(aes(fill = Counts))
ppp <- ppp + scale_fill_gradientn(colours=c("black",heat.colors(30)), values=c(seq(0,0.99),seq(1,30,length.out=0.5)))
ppp <- ppp + theme(plot.title = element_text(size=45),axis.text.x = element_text(angle = 45, hjust = 1))
ppp <- ppp + scale_y_discrete(breaks = breakpoints)
ppp <- ppp + labs(title= "Observation Counts",x="Features",y="SubjectID(Partially labeled)") 
print(ppp)

dev.off()
```



###Preprocessing raw training data(2424 obs. of 172 variables)###
```{r, warning=FALSE}
##read data
df0 <- df

#remove Adverse Event
df1 <- df0[df0$form_name != "Adverse Event", ] 

# select features
df1 <- df1[df1$feature_name %in% c("Age","Gender", "onset_delta" ,"onset_site","ALSFRS_Total","mouth","hands", "trunk","leg","respiratory", "fvc_percent","bp_diastolic","bp_systolic", "BMI", "pulse", "if_use_Riluzole", "Urine Ph","Albumin", "Protein", "Sodium", "Potassium","Bicarbonate", "Chloride", "Blood Urea Nitrogen (BUN)", "Creatinine","Alkaline Phosphatase","ALT(SGPT)", "Gamma-glutamyltransferase", "AST(SGOT)", "Bilirubin (total)", "White Blood Cell (WBC)", "Neutrophils","Lymphocytes", "Monocytes", "Eosinophils","Basophils","Red Blood Cells (RBC)","Hemoglobin", "Hematocrit","Platelets","CK","Triglycerides", "Total Cholesterol", "Glucose","HbA1c (Glycated Hemoglobin)", "Calcium", "Phosphorus", "Race"),]
### Organized the data
# change "feature delta" as numeric feature
df1$feature_delta<-as.numeric(df1$feature_delta)

#Organized the data format 
colnames(df1)[1]<-"SubjectID"
#Read slope values 
slope <- df_slope
dff<-merge(slope,df1,by="SubjectID",all.x=TRUE)
feature_lst <- unique(dff$feature_name)
list_length<-length(feature_lst) 
lengthh<-c(1:list_length)

IID<-order(dff$feature_delta)
dff<-dff[IID,]

#deal with constant numeric features
constant_id<-c(30,37,41,45)
for (i in constant_id){ 
  aaa<-dff[which(dff$feature_name==feature_lst[i]),]
  colnames(aaa)[1]<-"SubjectID"
  ID_lst <- unique(aaa$ID)
  aaa$feature_value<-as.numeric(as.factor(aaa$feature_value))
  cof<-ddply(aaa,~SubjectID,summarise,mean=mean(feature_value) )
  colnames(cof)[2]<-paste(feature_lst[i],"mean",sep="_") 
  slope<-merge(slope,cof,by="SubjectID",all=TRUE)
} 

# deal with constant categorical feature
for (i in c(32,11)){ 
  aaa<-dff[which(dff$feature_name==feature_lst[i]),]
  colnames(aaa)[1]<-"SubjectID"
  ID_lst <- unique(aaa$ID)
  aaa$feature_value<-as.numeric(aaa$feature_value)
  cof<-ddply(aaa,~SubjectID,summarise,mean=mean(feature_value) )
  colnames(cof)[2]<-paste(feature_lst[i],"mean",sep="_") 
  slope<-merge(slope,cof,by="SubjectID",all=TRUE)
} 


# deal with time-varying features
for (i in lengthh[-c(30,32,37,41,45,11)]){
  aaa<-dff[which(dff$feature_name==feature_lst[i]),]
  colnames(aaa)[1]<-"SubjectID"
  ID_lst <- unique(aaa$ID)
  aaa$feature_value<-as.numeric(aaa$feature_value)
  cof<-ddply(aaa,~SubjectID,summarise,max=max(feature_value),min=min(feature_value),
             median=median(feature_value)
  )
  #range 1
  cof1<-ddply(aaa,~SubjectID,summarise,range=range(feature_delta))
  nn<-nrow(cof1)
  nn_2<-seq(2,nn,by=2)
  nn_1<-nn_2-1
  xx<-cof1[nn_2,2]
  xx_1<-cof1[nn_1,2]
  range_1<-xx-xx_1
  #range 2
  cof1<-ddply(aaa,~SubjectID,summarise,range=range(feature_value))
  nn<-nrow(cof1)
  nn_2<-seq(2,nn,by=2)
  nn_1<-nn_2-1
  xx<-cof1[nn_2,2]
  xx_1<-cof1[nn_1,2]
  range_2<-xx-xx_1
  range<-range_2/range_1
  cof$range<-range
  colnames(cof)[2:5]<-paste(feature_lst[i],c("max","min","median","range"),sep="_")
  slope<-merge(slope,cof,by="SubjectID",all=TRUE)
}
## write dataset into local directory 
write.csv(slope,"Data/slope_raw_training_new.csv", row.names = FALSE)

```

### Preprocessing raw testing dataset (101 obs. of 172 variables)###
```{r, warning= FALSE}
df0_v <- df_v

#remove Adverse Event
df1_v <- df0_v[df0_v$form_name != "Adverse Event", ] 

df1_v <- df1_v[df1_v$feature_name %in% c("Age","Gender", "onset_delta" ,"onset_site","ALSFRS_Total","mouth","hands","trunk","leg","respiratory", "fvc_percent","bp_diastolic","bp_systolic", "BMI",
"pulse", "if_use_Riluzole", "Urine Ph","Albumin", "Protein", "Sodium", "Potassium","Bicarbonate", "Chloride","Blood Urea Nitrogen (BUN)", "Creatinine","Alkaline Phosphatase","ALT(SGPT)", "Gamma-glutamyltransferase", "AST(SGOT)","Bilirubin (total)", "White Blood Cell (WBC)", "Neutrophils","Lymphocytes", "Monocytes","Eosinophils","Basophils","Red Blood Cells (RBC)","Hemoglobin", "Hematocrit","Platelets","CK","Triglycerides","Total Cholesterol", "Glucose","HbA1c (Glycated Hemoglobin)", "Calcium", "Phosphorus", "Race"),]

df1<-df1_v
df1$feature_delta<-as.numeric(df1$feature_delta)

colnames(df1)[1]<-"SubjectID"
slope<- df_v_slope

dff<-merge(slope,df1,by="SubjectID",all.x=TRUE)
feature_lst <- unique(dff$feature_name)
list_length<-length(feature_lst) 
lengthh<-c(1:list_length)
#library(lme4)

IID<-order(dff$feature_delta)
dff<-dff[IID,]

## deal with constant numeric features
constant_id<-c(7,12,13,4)
for (i in constant_id){ 
  aaa<-dff[which(dff$feature_name==feature_lst[i]),]
  colnames(aaa)[1]<-"SubjectID"
  ID_lst <- unique(aaa$ID)
  aaa$feature_value<-as.numeric(as.factor(aaa$feature_value))
  cof<-ddply(aaa,~SubjectID,summarise,mean=mean(feature_value) )
  colnames(cof)[2]<-paste(feature_lst[i],"mean",sep="_") 
  slope<-merge(slope,cof,by="SubjectID",all=TRUE)
} 

#deal with constant categorical features
for (i in c(8,11)){ 
  aaa<-dff[which(dff$feature_name==feature_lst[i]),]
  colnames(aaa)[1]<-"SubjectID"
  ID_lst <- unique(aaa$ID)
  aaa$feature_value<-as.numeric(aaa$feature_value)
  cof<-ddply(aaa,~SubjectID,summarise,mean=mean(feature_value) )
  colnames(cof)[2]<-paste(feature_lst[i],"mean",sep="_") 
  slope<-merge(slope,cof,by="SubjectID",all=TRUE)
} 


#deal with time-varying features
for (i in lengthh[-c(7,12,13,4,8,11)])  {
  aaa<-dff[which(dff$feature_name==feature_lst[i]),]
  colnames(aaa)[1]<-"SubjectID"
  ID_lst <- unique(aaa$ID)
  aaa$feature_value<-as.numeric(aaa$feature_value)
  cof<-ddply(aaa,~SubjectID,summarise,max=max(feature_value),min=min(feature_value),
             median=median(feature_value)
  )
  #range 1
  cof1<-ddply(aaa,~SubjectID,summarise,range=range(feature_delta))
  nn<-nrow(cof1)
  nn_2<-seq(2,nn,by=2)
  nn_1<-nn_2-1
  xx<-cof1[nn_2,2]
  xx_1<-cof1[nn_1,2]
  range_1<-xx-xx_1
  #range 2
  cof1<-ddply(aaa,~SubjectID,summarise,range=range(feature_value))
  nn<-nrow(cof1)
  nn_2<-seq(2,nn,by=2)
  nn_1<-nn_2-1
  xx<-cof1[nn_2,2]
  xx_1<-cof1[nn_1,2]
  range_2<-xx-xx_1
  range<-range_2/range_1
  cof$range<-range
  colnames(cof)[2:5]<-paste(feature_lst[i],c("max","min","median","range"),sep="_")
  slope<-merge(slope,cof,by="SubjectID",all=TRUE)
}
# write file 
write.csv(slope,"Data/slope_raw_validate_new.csv", row.names = FALSE)
```

###Deal with missing data###
```{r}
data1<-read.csv("Data/slope_raw_training_new.csv",header=TRUE)
data2<-read.csv("Data/slope_raw_validate_new.csv",header=TRUE)
```

###Raw analysis###
####Plot raw missing pattern plot#### 
```{r, echo=FALSE,eval=FALSE}
ggplot_missing(data1)

```

```{r, echo=FALSE, eval=FALSE}
ggplot_missing(data2)

```

```{r, eval=FALSE}
##Descriptives 
write.csv(describe(data1), "Data/Descriptive_training_Ncomplete.csv")
write.csv(describe(data2), "Data/Descriptive_testing_Ncomplete.csv")
```


####deal with missing values in PRO-ACT testing data first (78 obs. of 100 variables) #####
```{r}
data2<-data2[,order(colnames(data2))]

ID<-which(apply(is.na(data2),1,sum)<136)
data2_1<-data2[ID,]

col_ID<-which(apply(is.na(data2_1),2,sum)>9)
data2_2<-data2_1[,-col_ID]
for(i in 1:ncol(data2_2)){
  m<-median(data2_2[,i])
  for (j in 1:nrow(data2_2)){
    if(is.na(data2_2[j,i])){data2_2[j,i]<-m}
    if(is.infinite(data2_2[j,i])){data2_2[j,i]<-m}
  }
}

data2_2=data.frame(apply(data2_2,2,f))

#delete useless covariate, and unchanged BMI
data2_2<-data2_2[,-c(133,35,36)]
col_2<-colnames(data2_2)
```

#### Then deal with missing values in PRO-ACT training data (2223 obs. of 100 variables)####
```{r}
data1<-data1[,order(colnames(data1))]
col_chose<-is.element(colnames(data1),col_2)
data1_1<-data1[,col_chose]
data1_2<-data1_1[,which(apply(is.na(data1_1),2,sum)<400)]
data1_3<-data1_2[which(apply(is.na(data1_2),1,sum)<20),]
for(i in 1:ncol(data1_3)){
  m<-median(data1_3[,i])
  for (j in 1:nrow(data1_3)){
    if(is.na(data1_3[j,i])){data1_3[j,i]<-m}
    if(is.infinite(data1_3[j,i])){data1_3[j,i]<-m}
  }
}
data1_3<-data.frame(apply(data1_3,2,f))

data2_2<-data2_2[,is.element(colnames(data2_2),colnames(data1_3))] 
data2_2<-data2_2[,order(colnames(data2_2))]
data1_3<-data1_3[,order(colnames(data1_3))]
#Delete SubjectID and Race
write.csv(data2_2[,-85],"Data/complete_testing.csv", row.names = FALSE)
write.csv(data1_3[,-85],"Data/complete_training.csv", row.names = FALSE)
```


### Predictive Analytics###

#### Linear Regression ####
```{r, warning=FALSE, results=FALSE}
training <- read.csv("Data/complete_training.csv",header=TRUE)
testing <- read.csv("Data/complete_testing.csv",header=TRUE)
obj<-lm(ALSFRS_slope~.,data=training) 
obj2<-step(obj, trace = FALSE) 
y_p<-predict(obj2,testing) 
test.y<-testing$ALSFRS_slope 
write.csv(y_p,"LR_predict.csv", row.names = FALSE)
```

####Random forest ####
```{r, message=FALSE, warning=FALSE}
set.seed(7788888)
#load datasets
df_train <- read.csv("Data/complete_training.csv")
df_test <- read.csv("Data/complete_testing.csv")

#data type transformation
df_train$onset_site_mean <- as.factor(df_train$onset_site_mean)
df_test$onset_site_mean <- as.factor(df_test$onset_site_mean)
levels(df_train$onset_site_mean)[4] <- "4"

df_train$Gender_mean <- as.factor(df_train$Gender_mean)
df_test$Gender_mean <- as.factor(df_test$Gender_mean)

rf.fit <-  train(ALSFRS_slope~., data = df_train[,-93], method = "rf", trControl = 
                   trainControl(method = "cv", number = 5), allowParallel=TRUE)

rf.result <- data.frame(df_test$ALSFRS_slope)
colnames(rf.result)[1] <- "ALSFRS_slope"
rf.result$prediction <- predict(rf.fit, df_test)

write.csv(rf.result,"RF_predict.csv",row.names = FALSE)
```

#### Bayesian Additive Regression Trees (BART) ####
```{r, eval=FALSE, message=FALSE, results=FALSE}
#########BART###############
# set_bart_machine_memory(10000)
set.seed(614)
data1<-read.csv('Data/complete_training.csv',header=TRUE)
train.x<-data1[,-c(6,93)]
train.y<-data1[,-c(6,93)]
data2<-read.csv('Data/complete_testing.csv',header=TRUE)
test.x<-data2[,-c(6,93)]
test.y<-data2[,-c(6,93)]
options(java.parameters="-Xmx68000m")
set_bart_machine_num_cores(8)
bart_machine_cv <- bartMachineCV(train.x, train.y, mem_cache_for_speed = FALSE)##CV
summary(bart_machine_cv)
y_p<-bart_predict_for_test_data(bart_machine_cv,test.x,test.y)$y_hat
write.csv(y_p,"bart_predict.csv")
```

####Super Learner####
```{r, warning=FALSE, message=FALSE, results=FALSE}
set.seed(23432)
data1<-read.csv("Data/complete_training.csv", header =TRUE)
train.x<-data1[,-c(6,93)]
train.y<-data1[,-c(6,93)]
data2<-read.csv("Data/complete_testing.csv",header=TRUE)
test.x<-data2[,-c(6,93)]
test.y<-data2[,-c(6,93)]
test.x<-test.x[,order(colnames(test.x))]
train.x<-train.x[,order(colnames(train.x))]

# generate Library and run Super Learner
SL.library <- c("SL.svm", "SL.ridge", "SL.randomForest","SL.mean","SL.caret", "SL.rpart","SL.stepAIC","SL.step.forward", "SL.step")

test_predict <- SampleSplitSuperLearner(Y = train.y, X = train.x, newX = test.x, SL.library = SL.library, verbose = FALSE, method = "method.NNLS")

y_p_sl<-test_predict$SL.predict
write.csv(y_p_sl,"SL_predict.csv",row.names = FALSE)
```



### Measurement of accuracy###
```{r, fig.keep='all',echo = FALSE, warning=FALSE, message=FALSE}
#aggregate prediction results from 4 different methods
result_lr <- read.csv("LR_predict.csv")
result_rf <- read.csv("RF_predict.csv")
result_rf <- result_rf[,-1]
result_bt <- read.csv("bart_predict.csv")
result_bt <- result_bt[,-1]
result_sl <- read.csv("SL_predict.csv")

results <- data.frame(cbind(testing$ALSFRS_slope, result_lr, result_rf,result_bt, result_sl))
colnames(results) <- c("ALSFRS_slope","LinearRegression","RandomForests","BART","SuperLearner")

#calculate measurements of accuracy
##LinearRegression
ss_t<-sum((results$ALSFRS_slope-mean(results$ALSFRS_slope))^2) 
ss_r<-sum((results$ALSFRS_slope-results$LinearRegression)^2)
lr_r2 <- 1-ss_r/ss_t #0.4185877
lr_rmse <- sqrt(mean((results$ALSFRS_slope - results$LinearRegression)^2)) #0.4920588
lr_cor <- cor(results$ALSFRS_slope,results$LinearRegression) #0.6635726
lr <- c(lr_r2, lr_rmse,lr_cor)
##RandomForests
ss_r<-sum((results$ALSFRS_slope-results$RandomForests)^2)
rf_r2 <- 1-ss_r/ss_t #0.6112357
rf_rmse <- sqrt(mean((results$ALSFRS_slope - results$RandomForests)^2)) #0.4023631
rf_cor <- cor(results$ALSFRS_slope,results$RandomForests) #0.7858035
rf <- c(rf_r2, rf_rmse,rf_cor)
##BART
ss_r<-sum((results$ALSFRS_slope-results$BART)^2)
bt_r2 <- 1-ss_r/ss_t #0.4185877
bt_rmse <- sqrt(mean((results$ALSFRS_slope - results$BART)^2)) #0.4920588
bt_cor <- cor(results$ALSFRS_slope,results$BART) #0.6635726
bt <- c(bt_r2, bt_rmse,bt_cor)
##SuperLearner
ss_r<-sum((results$ALSFRS_slope-results$SuperLearner)^2)
sl_r2 <- 1-ss_r/ss_t #0.5363443
sl_rmse <- sqrt(mean((results$ALSFRS_slope - results$SuperLearner)^2)) #0.4394128
sl_cor <- cor(results$ALSFRS_slope,results$SuperLearner) #0.7423889
sl <- c(sl_r2, sl_rmse,sl_cor)

#make a table
measurement_lst <- c("R-squared", "RMSE", "Correlation")
result.table <- data.frame(cbind(measurement_lst, lr, rf, bt, sl))
colnames(result.table) <- c("measurement", "LinearRegression","RandomForests","BART","SuperLearner")
print(result.table)

#Compare error distribution among 4 methods
err_results <- results[,c(2,3,4,5)]-results$ALSFRS_slope
l <- reshape(err_results, 
             varying = c("LinearRegression", "RandomForests","BART","SuperLearner"), 
             v.names = "value",
             times=c("LinearRegression", "RandomForests","BART","SuperLearner"),
             direction = "long")
colnames(l)[1]<-"Technique"
l$Technique<-factor(l$Technique, levels=c("LinearRegression","RandomForests","BART","SuperLearner"))

p <- ggplot(l, aes(x=Technique, y=value, fill=Technique)) + geom_boxplot()+xlab("Method")+
    ylab("Error(Y_hat-Y)")+ggtitle("Difference between Y and its prediction")
print(p)

```


###Extract useful information from adverse events ###
####Adverse Events in raw training dataset####
```{r, warning=FALSE, message=FALSE, results=FALSE}
##train_AE.csv obtained by excel
df1<-read.csv("train_AE.csv", header=TRUE, fill = TRUE,  stringsAsFactors = FALSE)
ID<-grep("Gastrointestinal",df1$feature_name)

df1$resp<-0

ID1<-intersect(grep("failure",df1$feature_value),grep("Respiratory",df1$feature_name)) 
ID2<-intersect(grep("Choking",df1$feature_value),grep("Respiratory",df1$feature_name)) 
ID3<-intersect(grep("insufficiency",df1$feature_value),grep("Respiratory",df1$feature_name)) 
ID4<-intersect(grep("Orthopnea",df1$feature_value),grep("Respiratory",df1$feature_name)) 
ID5<-intersect(grep("distress",df1$feature_value),grep("Respiratory",df1$feature_name)) 
ID6<-intersect(grep("Hoarse voice",df1$feature_value),grep("Respiratory",df1$feature_name)) 

df1$resp[ID1]<-1
df1$resp[ID2]<-df1$resp[ID2]+1
df1$resp[ID3]<-df1$resp[ID3]+1
df1$resp[ID4]<-df1$resp[ID4]+1
df1$resp[ID5]<-df1$resp[ID5]+1
df1$resp[ID6]<-df1$resp[ID6]+1

df1$infection<-0
ID1<-intersect(grep("upper respiratory",df1$feature_value,ignore.case = TRUE),grep("infection",df1$feature_name,ignore.case = TRUE)) 
df1$infection[ID1]<-1

df1$muscle<-0
ID1<-intersect(grep("weakness",df1$feature_value,ignore.case = TRUE),grep("muscle",df1$feature_name,ignore.case = TRUE)) 
df1$muscle[ID1]<-1

df1$other<-0
ID1<-grep("weight decrease",df1$feature_value,ignore.case = TRUE) 
df1$other[ID1]<-1

df1$injuries<-0
ID1<-intersect(grep("fall",df1$feature_value,ignore.case = TRUE),grep("injuries",df1$feature_name,ignore.case = TRUE)) 
df1$injuries[ID1]<-1

df1$move<-0
ID1<-intersect(grep("tremor",df1$feature_value,ignore.case = TRUE),grep("movement",df1$feature_name,ignore.case = TRUE)) 
df1$move[ID1]<-1

df1$musculoskeletal<-0
ID1<-intersect(grep("chewing",df1$feature_value,ignore.case = TRUE),grep("musculoskeletal",df1$feature_name,ignore.case = TRUE)) 
df1$musculoskeletal[ID1]<-1

df1$neurological<-0
ID1<-intersect(grep("speech",df1$feature_value,ignore.case = TRUE),grep("neurological",df1$feature_name,ignore.case = TRUE)) 
df1$neurological[ID1]<-1

df1$neuromuscular<-0
ID1<-intersect(grep("spasticity",df1$feature_value,ignore.case = TRUE),grep("neuromuscular",df1$feature_name,ignore.case = TRUE)) 
df1$neuromuscular[ID1]<-1

df1$cranial_nerve<-0
ID1<-intersect(grep("paralysis",df1$feature_value,ignore.case = TRUE),grep("cranial nerve",df1$feature_name,ignore.case = TRUE)) 
df1$cranial_nerve[ID1]<-1

df1$neurological<-0
ID1<-intersect(grep("dysarthria",df1$feature_value,ignore.case = TRUE),grep("neurological",df1$feature_name,ignore.case = TRUE)) 
df1$neurological[ID1]<-1


df1$peripheral_neuropath<-0
ID1<-intersect(grep("foot drop",df1$feature_value,ignore.case = TRUE),grep("peripheral neuropath",df1$feature_name,ignore.case = TRUE)) 
df1$peripheral_neuropath[ID1]<-1

df1<-as.data.frame(df1)

df_disease1<-df1[,c(2,6:16)]
aaa1<-rowsum(df_disease1,df_disease1$SubjectID)
aaa1$SubjectID<-as.numeric(rownames(aaa1))
write.csv(aaa1,"AE_mapping_train_wizID.csv", row.names = FALSE)
```
####Adverse Events in raw testing dataset####
```{r, warning=FALSE, message=FALSE, results=FALSE}
##test_AE.csv obtained by excel
df1<-read.csv("test_AE.csv", header=TRUE, fill = TRUE,  stringsAsFactors = FALSE)
ID<-grep("Gastrointestinal",df1$feature_name)

df1$resp<-0

ID1<-intersect(grep("failure",df1$feature_value),grep("Respiratory",df1$feature_name)) 
ID2<-intersect(grep("Choking",df1$feature_value),grep("Respiratory",df1$feature_name)) 
ID3<-intersect(grep("insufficiency",df1$feature_value),grep("Respiratory",df1$feature_name)) 
ID4<-intersect(grep("Orthopnea",df1$feature_value),grep("Respiratory",df1$feature_name)) 
ID5<-intersect(grep("distress",df1$feature_value),grep("Respiratory",df1$feature_name)) 
ID6<-intersect(grep("Hoarse voice",df1$feature_value),grep("Respiratory",df1$feature_name)) 

df1$resp[ID1]<-1
df1$resp[ID2]<-df1$resp[ID2]+1
df1$resp[ID3]<-df1$resp[ID3]+1
df1$resp[ID4]<-df1$resp[ID4]+1
df1$resp[ID5]<-df1$resp[ID5]+1
df1$resp[ID6]<-df1$resp[ID6]+1

df1$infection<-0
ID1<-intersect(grep("upper respiratory",df1$feature_value,ignore.case = TRUE),grep("infection",df1$feature_name,ignore.case = TRUE)) 
df1$infection[ID1]<-1

df1$muscle<-0
ID1<-intersect(grep("weakness",df1$feature_value,ignore.case = TRUE),grep("muscle",df1$feature_name,ignore.case = TRUE)) 
df1$muscle[ID1]<-1

df1$other<-0
ID1<-grep("weight decrease",df1$feature_value,ignore.case = TRUE) 
df1$other[ID1]<-1

df1$injuries<-0
ID1<-intersect(grep("fall",df1$feature_value,ignore.case = TRUE),grep("injuries",df1$feature_name,ignore.case = TRUE)) 
df1$injuries[ID1]<-1

df1$move<-0
ID1<-intersect(grep("tremor",df1$feature_value,ignore.case = TRUE),grep("movement",df1$feature_name,ignore.case = TRUE)) 
df1$move[ID1]<-1

df1$musculoskeletal<-0
ID1<-intersect(grep("chewing",df1$feature_value,ignore.case = TRUE),grep("musculoskeletal",df1$feature_name,ignore.case = TRUE)) 
df1$musculoskeletal[ID1]<-1

df1$neurological<-0
ID1<-intersect(grep("speech",df1$feature_value,ignore.case = TRUE),grep("neurological",df1$feature_name,ignore.case = TRUE)) 
df1$neurological[ID1]<-1

df1$neuromuscular<-0
ID1<-intersect(grep("spasticity",df1$feature_value,ignore.case = TRUE),grep("neuromuscular",df1$feature_name,ignore.case = TRUE)) 
df1$neuromuscular[ID1]<-1

df1$cranial_nerve<-0
ID1<-intersect(grep("paralysis",df1$feature_value,ignore.case = TRUE),grep("cranial nerve",df1$feature_name,ignore.case = TRUE)) 
df1$cranial_nerve[ID1]<-1

df1$neurological<-0
ID1<-intersect(grep("dysarthria",df1$feature_value,ignore.case = TRUE),grep("neurological",df1$feature_name,ignore.case = TRUE)) 
df1$neurological[ID1]<-1

df1$peripheral_neuropath<-0
ID1<-intersect(grep("foot drop",df1$feature_value,ignore.case = TRUE),grep("peripheral neuropath",df1$feature_name,ignore.case = TRUE)) 
df1$peripheral_neuropath[ID1]<-1

df1<-as.data.frame(df1)

df_disease2<-df1[,c(2,8:18)]
aaa2<-rowsum(df_disease2,df_disease2$SubjectID)
aaa2$SubjectID<-as.numeric(rownames(aaa2))
write.csv(aaa2,"AE_mapping_test_wizID.csv", row.names = FALSE)
```

####aggregate adverse events information####
```{r, warning=FALSE, message=FALSE, results=FALSE}
training<-read.csv("Data/complete_training.csv",header=TRUE)
testing<-read.csv("Data/complete_testing.csv",header=TRUE)

training_ae<-read.csv("AE_mapping_train_wizID.csv",header=TRUE)
testing_ae<-read.csv("AE_mapping_test_wizID.csv",header=TRUE)

train<-merge(training,training_ae,by="SubjectID")
test<-merge(testing,testing_ae,by="SubjectID")

score_ae<-apply(train[,101:111],1,sum)
train$score_ae<-score_ae

score_ae<-apply(test[,101:111],1,sum)
test$score_ae<-score_ae

write.csv(train,"training_wizAE.csv", row.names = FALSE)
write.csv(test,"testing_wizAE.csv", row.names = FALSE)
```

###Investigation on variable importance(with/without adverse events)###
####RandomForests####
```{r, warning=FALSE, message=FALSE, results=FALSE}
#load datasets
df_train <- read.csv("training_wizAE.csv")
df_test <- read.csv("testing_wizAE.csv")
#remove unnecessary columns
df_train <- df_train[,-c(1,3)]
df_test <- df_test[,-c(1,3)]
#data type transformationï¼onset_site_meanï¼Gender_meanï¼
df_train$onset_site_mean <- as.factor(df_train$onset_site_mean)
df_test$onset_site_mean <- as.factor(df_test$onset_site_mean)
levels(df_train$onset_site_mean)[4] <- "4"

df_train$Gender_mean <- as.factor(df_train$Gender_mean)
df_test$Gender_mean <- as.factor(df_test$Gender_mean)
################RandomForests#######################
set.seed(49239)
#without AE
rf.fit1 <-  train(ALSFRS_slope~., data = df_train[,-c(1,100:111)], method = "rf", trControl = 
                     trainControl(method = "cv", number = 5), allowParallel=TRUE, importance = TRUE)
#with AE
rf.fit2 <-  train(ALSFRS_slope~., data = df_train[,-1], method = "rf", trControl = 
                      trainControl(method = "cv", number = 5), allowParallel=TRUE, importance = TRUE)

#obtain importance measures of features in the fitted models
rf_noAE.ImpMeasure <- data.frame(varImp(rf.fit1)$importance)
rf_noAE.ImpMeasure$Vars <- row.names(rf_noAE.ImpMeasure)
rf_noAE.ImpMeasure <- rf_noAE.ImpMeasure[order(-rf_noAE.ImpMeasure$Overall),]

rf_AE.ImpMeasure <- data.frame(varImp(rf.fit2)$importance)
rf_AE.ImpMeasure$Vars <- row.names(rf_AE.ImpMeasure)
rf_AE.ImpMeasure <- rf_AE.ImpMeasure[order(-rf_AE.ImpMeasure$Overall),]

#output
write.csv(rf_noAE.ImpMeasure, "RF_ImpMeasure_withoutAE.csv", row.names = FALSE)
write.csv(rf_AE.ImpMeasure, "RF_ImpMeasure_withAE.csv", row.names = FALSE)
```

####BART####
```{r, warning=FALSE, message=FALSE, results=FALSE}
options(java.parameters="-Xmx68000m")
library(bartMachine)
set.seed(614)
#without AE
data1<-read.csv('training_wizAE.csv',header=TRUE)
train.x<-data1[,-c(1,7,101:112)]
train.y<-data1[,7]
data2<-read.csv('testing_wizAE.csv',header=TRUE)
test.x<-data2[,-c(1,7,101:112)]
test.y<-data2[,7]
set_bart_machine_num_cores(8)
bart_machine_cv <- bartMachineCV(train.x, train.y, mem_cache_for_speed = FALSE)##CV
aaa1<-investigate_var_importance(bart_machine_cv,type="splits", plot=FALSE, num_replicates_for_avg=5,num_trees_bottleneck=20)
score1<-aaa1$avg_var_props
write.csv(score1,"BART_ImpMeasure_withoutAE.csv")
#with AE
data1<-read.csv('training_wizAE.csv',header=TRUE)
train.x<-data1[,-c(1,7)]
train.y<-data1[,7]
data2<-read.csv('testing_wizAE.csv',header=TRUE)
test.x<-data2[,-c(1,7)]
test.y<-data2[,7]
set_bart_machine_num_cores(8)
bart_machine_cv <- bartMachineCV(train.x, train.y, mem_cache_for_speed = FALSE)##CV
aaa2<-investigate_var_importance(bart_machine_cv,type="splits", plot=FALSE, num_replicates_for_avg=5,num_trees_bottleneck=20)
score2<-aaa2$avg_var_props
write.csv(score2,"BART_ImpMeasure_withAE.csv")
```
