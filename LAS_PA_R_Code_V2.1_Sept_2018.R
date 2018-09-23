
setwd("~/Google Drive/Vibrant/MingTang_ChaoGao_Projects/ALS Project/Paper outline/Code/V2.0")
#Are there any changes in the rates of disease for those 
#with more AEs based on the ALSFRS-R, FVC%, BMI, etc.?_

library(dplyr)
library(plyr)
library(mi)
library(reshape)
library(reshape2)
library(ggplot2)
library(MASS)
library(magrittr)
library(psych) 
library(bartMachine)
library(SuperLearner)
library(caret)
library(e1071)
library(Amelia)
######################
####Data Manipulation####
df <- read.table("all_forms_PROACT.txt", header=TRUE, sep="|", na.strings=c("NA","","NaN"), fill = TRUE, quote = "", stringsAsFactors = FALSE, comment.char="")
df_slope <- read.table("ALSFRS_slope_PROACT.txt", header=TRUE, sep="|", na.strings=c("NA","","NaN"), fill = TRUE, quote = "", stringsAsFactors = FALSE, comment.char="")
#remove study subjects with NA in ALSFRS_slope
df_slope <- df_slope[!is.na(df_slope$ALSFRS_slope), ]
df_v <- read.table("all_forms_validate_spike.txt", header=TRUE, sep="|", na.strings=c("NA","","NaN"), fill = TRUE, quote = "", stringsAsFactors = FALSE, comment.char="")
df_v_slope <- read.table("ALSFRS_slope_validate_spike.txt",header=TRUE,sep="|", na.strings=c("NA","","NaN"), fill = TRUE, quote = "", stringsAsFactors = FALSE, comment.char="")

df0 <- df
#length(unique(df0$feature_name))

#remove Adverse Event
df1 <- df0[df0$form_name != "Adverse Event", ] 
df1$feature_delta<-as.numeric(df1$feature_delta)
df1 <- df1[df1$feature_delta <= 92,]

dfff<-merge(df1,df_slope,by="SubjectID",all.y=TRUE)

#length(unique(dfff$feature_name))

aa<-table(dfff$feature_name)
aa<-as.data.frame(aa)
#as.data.frame(table(aa$Freq))
#variable_list<-as.character(aa$Var1)

#dfff<-dfff[dfff$feature_name %in% variable_list,]


# select features
#df1 <- df1[df1$feature_name %in% c("Age","Gender", "onset_delta" ,"onset_site","ALSFRS_Total","mouth","hands", "trunk","leg","respiratory", "fvc_percent","bp_diastolic","bp_systolic", "BMI", "pulse", "if_use_Riluzole", "Urine Ph","Albumin", "Protein", "Sodium", "Potassium","Bicarbonate", "Chloride", "Blood Urea Nitrogen (BUN)", "Creatinine","Alkaline Phosphatase","ALT(SGPT)", "Gamma-glutamyltransferase", "AST(SGOT)", "Bilirubin (total)", "White Blood Cell (WBC)", "Neutrophils","Lymphocytes", "Monocytes", "Eosinophils","Basophils","Red Blood Cells (RBC)","Hemoglobin", "Hematocrit","Platelets","CK","Triglycerides", "Total Cholesterol", "Glucose","HbA1c (Glycated Hemoglobin)", "Calcium", "Phosphorus", "Race"),]
### Organized the data
# change "feature delta" as numeric feature
#df1$feature_delta<-as.numeric(df1$feature_delta)
dfff <- dfff[dfff$feature_delta <= 92,]


#Organized the data format 
#colnames(df1)[1]<-"SubjectID"
#Read slope values 
slope <- df_slope
#dff<-merge(slope,df1,by="SubjectID",all.x=TRUE)
feature_lst <- unique(dfff$feature_name)
list_length<-length(feature_lst) 
lengthh<-c(1:list_length)

IID<-order(dfff$feature_delta)
dfff<-dfff[IID,]

#deal with constant categorical features treatment_group (86)
constant_id<-c(14,49,73,86,48,89,75,101) #Gender(14),onset_site(49),Race(73),if_use_Riluzole(89) family_ALS_hist(101)

for (i in constant_id){ 
  aaa<-dfff[which(dfff$feature_name==feature_lst[i]),]
  colnames(aaa)[1]<-"SubjectID"
  ID_lst <- unique(aaa$ID)
  aaa$feature_value<-as.numeric(as.factor(aaa$feature_value))
  cof<-ddply(aaa,~SubjectID,summarise,mean=mean(feature_value))
  colnames(cof)[2]<-feature_lst[i]
  slope<-merge(slope,cof,by="SubjectID",all=TRUE)
} 

# deal with constant numeric features
for (i in c(48,75)){ #onset_delta(48),Age(75)
  aaa<-dfff[which(dfff$feature_name==feature_lst[i]),]
  colnames(aaa)[1]<-"SubjectID"
  ID_lst <- unique(aaa$ID)
  aaa$feature_value<-scale(as.numeric(aaa$feature_value) ) #scale it 
  cof<-ddply(aaa,~SubjectID,summarise,mean=mean(feature_value) )
  colnames(cof)[2]<-feature_lst[i]
  slope<-merge(slope,cof,by="SubjectID",all=TRUE)
} 


# deal with time-varying features #ignore 16
for (i in lengthh[-c(14,49,73,86,48,89,75,101,48,75)]){
  #print(i)
  aaa<-dfff[which(dfff$feature_name==feature_lst[i]),]
  aaa<-aaa[,-5]
  colnames(aaa)[1]<-"SubjectID"
  ID_lst <- unique(aaa$ID)
  aaa$feature_value<-scale(as.numeric(aaa$feature_value))
  cof<-ddply(aaa,~SubjectID,summarise,max=max(feature_value),min=min(feature_value),
             median=median(feature_value)
  )
  aaa<-na.omit(aaa)
  if(nrow(aaa)!=0){
    mods = dlply(aaa, .(SubjectID), lm, formula = feature_value ~ feature_delta)
    coefs = ldply(mods, coef)
    cof<-merge(cof,coefs[,c(1,3)],by="SubjectID",all.x=TRUE)
    colnames(cof)[2:5]<-paste(feature_lst[i],c("max","min","median","slope"),sep="_")
    slope<-merge(slope,cof,by="SubjectID",all=TRUE)
  }
}


#delete the variable less than 1500 observation
slope1<-slope[,-which(describe(slope)$n<1600)]
write.csv(slope1,"before_impute.csv",row.names = FALSE)


###################################
#######multiple imputation####

slope<-read.csv("before_impute.csv",header=TRUE)
slope2<-slope[,-c(11,173,174,175)]

for(i in 1:20){
  seed<-1000*i+234
  set.seed(seed)
  imp<-amelia(slope2,m=1,empri=0.01*nrow(slope2),idvars="SubjectID",
              noms=c("Gender","Race","treatment_group","onset_site","if_use_Riluzole"))
  write.amelia(obj=imp,file.stem=paste("outdata",i,sep=""))
}

setwd("~/Google Drive/Vibrant/MingTang_ChaoGao_Projects/ALS Project/Paper outline/Code/V2.0")
library(ggplot2)
library(reshape2)
library(ggplot2)
library(Amelia)
library(dplyr)
before_impute<-read.csv("before_impute.csv",header=TRUE)

ggplot_missing <- function(x){
  
  x %>% 
    is.na %>%
    melt %>%
    ggplot(data = .,
           aes(x = Var2,
               y = Var1)) +
    geom_raster(aes(fill = value)) +
    scale_fill_grey(name = "",
                    labels = c("Present","Missing")) +
    theme_minimal() + 
    theme(axis.text.x  = element_text(angle=45, vjust=0.5)) + 
    labs(x = "Variables in Dataset",
         y = "Rows / observations")+coord_flip()
}
missmap(before_impute)
ggplot_missing(before_impute)

##Prediction####
###BART####
library(dplyr)
library(plyr)
library(mi)
library(reshape)
library(reshape2)
library(ggplot2)
library(MASS)
library(magrittr)
library(psych) 
library(SuperLearner)
library(caret)
library(e1071)
library(glmnet)
library(arm)
library(crossval)
options(java.parameters="-Xmx28000m")
library(bartMachine)
#brary(bartMachine)
predfun.bt = function(train.x, train.y, test.x, test.y, negative) {
  set_bart_machine_num_cores(3)
  bt.fit = bartMachine::bartMachineCV(train.x,train.y,mem_cache_for_speed = FALSE)
  ynew = bart_predict_for_test_data(bt.fit,test.x,test.y)$y_hat
  # count TP, FP etc.
  #out = mean( (ynew - test.y)^2 )
  ss_t<-sum((test.y-mean(test.y))^2) 
  ss_r<-sum((test.y-ynew)^2)
  r2 <- 1-ss_r/ss_t #0.07949168
  rmse <- sqrt(mean((test.y - ynew)^2)) #0.61913986
  cor <- cor(test.y,ynew) #0.34183028
  out <- c(r2, rmse,cor)
  return(out)
}

cv_results<-matrix(NA,10,4)
for(i in 1:20){
  set.seed(7788888)
  data_train<-read.csv(file=paste("outdata",i,"1.csv",sep=""),header=TRUE)
  X <- data_train[, -c(1,2,3)]
  Y <- data_train[, 3]
  cv.out.bt = crossval::crossval(predfun.bt, X, Y, K = 5, B = 1)
  cv.out.bt
  results<-cv.out.bt$stat
  results[4]<-i
  cv_results[i,]<-results
}
cv_results<-as.data.frame(cv_results)
colnames(cv_results)<-c("R2","RMSE","COR","chain")

write.csv(cv_results,"bt_cv_results_1.csv",row.names = FALSE)


####random forest####
library(dplyr)
library(plyr)
library(mi)
library(reshape)
library(reshape2)
library(ggplot2)
library(MASS)
library(magrittr)
library(psych) 
library(bartMachine)
library(SuperLearner)
library(caret)
library(e1071)
library(glmnet)
library(arm)
library(crossval)

predfun.rf = function(train.x, train.y, test.x, test.y, negative) {
  rf.fit = randomForest::randomForest(train.y ~ ., data = train.x)
  ynew = predict(rf.fit, test.x)
  # count TP, FP etc.
  ss_t<-sum((test.y-mean(test.y))^2) 
  ss_r<-sum((test.y-ynew)^2)
  r2 <- 1-ss_r/ss_t #0.07949168
  rmse <- sqrt(mean((test.y - ynew)^2)) #0.61913986
  cor <- cor(test.y,ynew) #0.34183028
  out <- c(r2, rmse,cor)
  return(out)
}

cv_results<-matrix(NA,10,4)
for(i in 1:20){
  set.seed(7788888)
  data_train<-read.csv(file=paste("outdata",i,"1.csv",sep=""),header=TRUE)
  X <- data_train[, -c(1,2,3)]
  Y <- data_train[, 3]
  cv.out.rf = crossval::crossval(predfun.rf, X, Y, K = 5, B = 1)
  cv.out.rf
  results<-cv.out.rf$stat
  results[4]<-i
  cv_results[i,]<-results
}
cv_results<-as.data.frame(cv_results)
colnames(cv_results)<-c("R2","RMSE","COR","chain")

write.csv(cv_results,"rf_cv_results.csv",row.names = FALSE)

########################################
###Feature selection####

#####knock-off######
library(dplyr)
library(plyr)
library(mi)
library(reshape)
library(reshape2)
library(ggplot2)
library(MASS)
library(magrittr)
library(psych) 
library(SuperLearner)
library(caret)
library(e1071)
library(glmnet)
library(arm)
library(crossval)

library(knockoff)
for(j in 1:20){ 
  print(j)
  ud <- read.csv(file=paste("outdata",j,"1.csv",sep=""),header=TRUE)
  ud <- ud[,-c(1,2,6,11)]
  #cat_num <- c(2,3,4)# ,5,7)
  #ud[,cat_num] <- lapply(ud[,cat_num],as.factor)
  X<-ud[,-1]
  Y<-ud[,1]
  #ud <- read.csv("udall_norm_complete.csv")
  my.bt.ko <- function(X,Y,fdr=0.3,r=0.632,mtry=150,iter=1000,random=TRUE){
    counts <- NULL
    print(dim(X))
    nc <- ncol(X)
    for (i in 1:iter){
      if(i%%100 ==0){print(i)}
      bt.id <- sample(1:nc,mtry,replace = FALSE)
      ko_var <- NULL
      try(ko_var <- knockoff.filter(X[,bt.id],fdr=fdr,Y,randomize = random))
      if (length(ko_var$selected)!=0){counts <- c(counts,bt.id[as.numeric(ko_var$selected)])}
    }
    sort(table(counts),decreasing = TRUE)
  }
  set.seed(2017)
  counts_ko_ud <- my.bt.ko(ud[,-1],ud[,1],mtry=150,fdr=0.35,iter=1000) 
  nm_ud <- colnames(ud)[-c(1)]
  
  counts_ko_ud <- data.frame(Feature=nm_ud[names(counts_ko_ud)%>% as.numeric()],Frequency=counts_ko_ud %>% as.matrix())
  #proportion 
  #40 is mtry(number of column selected by each iteration)
  counts_ko_ud <- counts_ko_ud[1:20,] %>% mutate(Prop=Frequency/1000*(248/150))
  print(counts_ko_ud)
  write.csv(counts_ko_ud,file=paste("ko_top20_",j,".csv",sep=""),row.names=FALSE) 
}  

###random forest####
library(dplyr)
library(plyr)
library(mi)
library(reshape)
library(reshape2)
library(ggplot2)
library(MASS)
library(magrittr)
library(psych) 
library(SuperLearner)
library(caret)
library(e1071)
library(glmnet)
library(arm)
library(crossval)

# feature selection

for (j in 1:20){
  #j=1
  ud <- read.csv(file=paste("outdata",j,"1.csv",sep=""),header=TRUE)
  print(j)
  ud <- ud[,-c(1,2)]
  cat_num <- c(2,3,4)# ,5,7)
  ud[,cat_num] <- lapply(ud[,cat_num],as.factor)
  ud<-ud[,-9]
}
  
library(parallel)
library(doParallel)
#cluster <- makeCluster(detectCores() - 1) # convention to leave 1 core for OS
cluster <- makeCluster(4)
registerDoParallel(cluster)
  
#Random Forest (500 iterations)
N <- 100  ##
control <- trainControl(method="cv", number=5,allowParallel=TRUE)
set.seed(2017)
ud_var_lst <- c()

for (i in 1:N){
  print(i)
  samples <- sample(1:nrow(ud), 0.632*nrow(ud),replace = T)
  rf_mod <- train(ALSFRS_slope~., data=ud[samples,], method="rf", trControl=control,importance=TRUE)
  imp_temp <- as.data.frame(varImp(rf_mod)$importance) #%IncMSE
  imp_temp$var <- row.names(imp_temp)
  names(imp_temp) <- c("imp","var")
  imp_temp$imp <- as.numeric(as.character(imp_temp$imp))
  ud_var_lst <- c(ud_var_lst,imp_temp[order(-imp_temp$imp),]$var[1:20]) #top 20
}  
  
#stopCluster(cluster)
#registerDoSEQ()

ud_counts_rf <- as.data.frame(sort(table(ud_var_lst),decreasing=T)[1:20])
names(ud_counts_rf) <- c("Features","Frequency")
ud_counts_rf <- ud_counts_rf %>% mutate(Prop = Frequency/100)
write.csv(ud_counts_rf,file=paste("rf_top20_",j,".csv",sep=""),row.names=FALSE)
  
  
###summarize data##
#setwd("~/Google Drive/Vibrant/MingTang_ChaoGao_Projects/ALS Project/Paper outline/Code/V2.0/0908/results")

a1<-read.csv("bt_cv_results_1.csv",header=TRUE)
a2<-read.csv("bt_cv_results_2.csv") #,header=TRUE)

aaa<-rbind(a1,a2)
apply(aaa,2,mean)


a1<-read.csv("rf_cv_results_1.csv",header=TRUE)
a2<-read.csv("rf_cv_results_2.csv") #,header=TRUE)

aaa<-rbind(a1,a2)
apply(aaa,2,mean)

# top features 
a1<-read.csv("ko_top20_1.csv",header=TRUE)

for(i in 2:20){
  a2<-read.csv(file=paste("ko_top20_",i,".csv",sep=""))
  a1<-rbind(a1,a2)
}

feature_ko<-aggregate(. ~ Feature, a1, median)
feature_ko<-feature_ko[order(feature_ko$Prop,decreasing=TRUE),]
feature_ko<-feature_ko[c(1:20),]
write.csv(feature_ko,"top20_ko.csv",row.names=FALSE)


# top features 
a1<-read.csv("rf_top20_1.csv",header=TRUE)

for(i in 2:20){
  a2<-read.csv(file=paste("rf_top20_",i,".csv",sep=""))
  a1<-rbind(a1,a2)
}

feature_rf<-aggregate(. ~ Features, a1, median)
feature_rf<-feature_rf[order(feature_rf$Prop,decreasing=TRUE),]
feature_rf<-feature_rf[c(1:20),]
write.csv(feature_rf,"top20_rf.csv",row.names=FALSE)
fvc_org<-read.csv("rf_fcv_5folds.csv",header=TRUE)
fvc_wo<-read.csv("rf_fvc_wofvc.csv",header=TRUE)

####Summary Statistics ####
setwd("/Users/gchao/Desktop/ALS_extra")
df <- read.csv("outdata11.csv")
feature_summary <- c()
for (i in c(3,8,10:253)){
  feature_max <- max(df[,i])
  feature_min <- min(df[,i])
  feature_mean <- mean(df[,i])
  feature_sd <- sd(df[,i])
  summary_temp <- c(feature_max,feature_min,feature_mean,feature_sd)
  feature_summary <- rbind(feature_summary,summary_temp)
}
feature_summary <- as.data.frame(cbind(names(df)[c(3,8,10:253)],feature_summary))
row.names(feature_summary) <- 1:nrow(feature_summary)
names(feature_summary) <- c("variable","max","min","mean","sd")
write.csv(feature_summary,"summary_table_outdata11.csv",row.names = F)

###Clustering Analysis####
df2 <- read.csv("before_impute.csv")
feature_summary2 <- c()
for (i in c(2,7,9:256)){
  feature_max <- max(df2[,i],na.rm = T)
  feature_min <- min(df2[,i],na.rm = T)
  feature_mean <- mean(df2[,i],na.rm = T)
  feature_sd <- sd(df2[,i],na.rm = T)
  summary_temp <- c(feature_max,feature_min,feature_mean,feature_sd)
  feature_summary2 <- rbind(feature_summary2,summary_temp)
}
feature_summary2 <- as.data.frame(cbind(names(df2)[c(2,7,9:256)],feature_summary2))
row.names(feature_summary2) <- 1:nrow(feature_summary2)
names(feature_summary2) <- c("variable","max","min","mean","sd")
write.csv(feature_summary2,"summary_table_before_impute.csv",row.names = F)

library(dplyr)
df3 <- df %>% select(-X,-SubjectID,-ALSFRS_slope)
#Determine number of clusters
wss <- (nrow(df3)-1)*sum(apply(df3,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(df3, centers=i)$withinss)
plot(1:15, wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")

set.seed(2017)
# K-Means Cluster Analysis 1
fit0 <- kmeans(df3, 4) # 4 cluster solution
# append cluster assignment 1
df4 <- data.frame(df3, fit0$cluster)
# K-Means Cluster Analysis 2
fit2 <- kmeans(df3, 4) # 4 cluster solution
# append cluster assignment 2
df4 <- data.frame(df4, fit2$cluster)

#robustness of clustering results k=4
Niter <- 1000
cluster0_1 <- which(df4$fit0.cluster==1) #index from original cluster 1
cluster0_2 <- which(df4$fit0.cluster==2) #index from original cluster 2
cluster0_3 <- which(df4$fit0.cluster==3) #index from original cluster 3
cluster0_4 <- which(df4$fit0.cluster==4) #index from original cluster 4
set.seed(666)
stat <- c()
for (i in 1:Niter){
  fit_temp <- kmeans(df3, 4)$cluster #new fit
  c1_new <- table(fit_temp[cluster0_1]) #clustering result among original c1
  c2_new <- table(fit_temp[cluster0_2]) #clustering result among original c2
  c3_new <- table(fit_temp[cluster0_3]) #clustering result among original c3
  c4_new <- table(fit_temp[cluster0_4]) #clustering result among original c4
  stat <-cbind(stat,c(max(prop.table(c1_new)),max(prop.table(c2_new)),
                      max(prop.table(c3_new)),max(prop.table(c4_new))))
}

robust_result <- as.data.frame(cbind(rowMeans(stat),apply(stat,1,sd)))
names(robust_result) <- c("Mean_Proportions","SD")
write.csv(robust_result,"clustering(robustness).csv",row.names = F)
#plot
comparison_matrix <- c()
for (i in 1:nrow(df4)){
  compare <- ifelse(df4[i,"fit1.cluster"] == df4[,"fit1.cluster"]&df4[i,"fit2.cluster"] == df4[,"fit2.cluster"],1,0)
  comparison_matrix <- rbind(comparison_matrix,compare)
}
set.seed(2017)
sample_id <- sort(sample(1:2424,500))
mat <- comparison_matrix[sample_id,sample_id]
png("heatmap.png",width = 3000,height = 3000)
heatmap(mat)
dev.off()




#comparison_matrix <- as.data.frame(comparison_matrix)
#set.seed(2017)
#sample_id <- sort(sample(1:2424,250))
#library(ggplot2)
#library(reshape2)

#library(ggcorrplot)
#png("aaa.png",width = 1000, height = 1000)
#ggcorrplot(comparison_matrix[1:250,1:250])
#dev.off()


##########################################
#For 565 subjects assigned to cluster 1 in the 1st fit
fit1_lab1 <- df4 %>% filter(fit1.cluster==1)
table(fit1_lab1$fit2.cluster) #all 565 have label1 in the 2nd fit
#For 427 subjects assigned to cluster 2 in the 1st fit
fit1_lab2 <- df4 %>% filter(fit1.cluster==2)
table(fit1_lab2$fit2.cluster) #411 have label4 and 16 have label3 in the 2nd fit
#For 699 subjects assigned to cluster 2 in the 1st fit
fit1_lab3 <- df4 %>% filter(fit1.cluster==3)
table(fit1_lab3$fit2.cluster) #621 have label3 and 78 have label 2in the 2nd fit
#For 733 subjects assigned to cluster 2 in the 1st fit
fit1_lab4 <- df4 %>% filter(fit1.cluster==4)
table(fit1_lab4$fit2.cluster) #705 have label2 and 28 have label 2in the 2nd fit

###Trajectory Plots#######
setwd("C:/Users/Chao Gao/Desktop/Gdrive/ALS_extra")
#load library
library(dplyr)
library(ggplot2)
#load data
als1 <-read.csv("outdata11.csv")
als2 <- read.csv("als_kmeans_label.csv")
als3 <- read.table("all_forms_PROACT.txt", header=TRUE, sep="|", na.strings=c("NA","","NaN"), 
                   fill = TRUE, quote = "", stringsAsFactors = FALSE, comment.char="")

id_2424 <- als1$SubjectID 
c1_id <- als1[als2$fit0.cluster==1,]$SubjectID
c2_id <- als1[als2$fit0.cluster==2,]$SubjectID
c3_id <- als1[als2$fit0.cluster==3,]$SubjectID
c4_id <- als1[als2$fit0.cluster==4,]$SubjectID

############ Cluster and features of interest ############
als_onset_trunk <- als1 %>% select(SubjectID,onset_delta.x,trunk_max)
als_onset_trunk$cluster <- NA
als_onset_trunk[als_onset_trunk$SubjectID %in% c1_id, ]$cluster <- "Cluster1" 
als_onset_trunk[als_onset_trunk$SubjectID %in% c2_id, ]$cluster <- "Cluster2" 
als_onset_trunk[als_onset_trunk$SubjectID %in% c3_id, ]$cluster <- "Cluster3" 
als_onset_trunk[als_onset_trunk$SubjectID %in% c4_id, ]$cluster <- "Cluster4" 

png("onset_delta_cluster1.png", width = 1500, height = 1500, res=300)
fill <- "#6699CC"
ggplot(als_onset_trunk[als_onset_trunk$cluster=="Cluster1",], aes(x = onset_delta.x)) + 
  geom_density(fill = fill, alpha = 0.8) + ylab("Density") + 
  theme(text = element_text(size=15), axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15), axis.title.x=element_blank()) +
  xlim(0,1100) + ylim(0,0.005)
dev.off()

png("trunk_cluster1.png", width = 1500, height = 1500, res=300)
fill <- "#6699CC"
ggplot(als_onset_trunk[als_onset_trunk$cluster=="Cluster1",], aes(x = trunk_max)) + 
  geom_density(fill = fill, alpha = 0.8) + ylab("Density")+ 
  theme(text = element_text(size=15), axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15), axis.title.x=element_blank()) +
  xlim(-3,2) + ylim(0,0.8)
dev.off()

png("onset_delta_cluster2.png", width = 1500, height = 1500, res=300)
fill <- "#009E73"
ggplot(als_onset_trunk[als_onset_trunk$cluster=="Cluster2",], aes(x = onset_delta.x)) + 
  geom_density(fill = fill, alpha = 0.8) + ylab("Density")+ 
  theme(text = element_text(size=15), axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15), axis.title.x=element_blank()) +
  xlim(0,1100) + ylim(0,0.005)
dev.off()

png("trunk_cluster2.png", width = 1500, height = 1500, res=300)
fill <- "#009E73"
ggplot(als_onset_trunk[als_onset_trunk$cluster=="Cluster2",], aes(x = trunk_max)) + 
  geom_density(fill = fill, alpha = 0.8) + ylab("Density")+ 
  theme(text = element_text(size=15), axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15), axis.title.x=element_blank()) +
  xlim(-3,2) + ylim(0,0.8)
dev.off()

png("onset_delta_cluster3.png", width = 1500, height = 1500, res=300)
fill <- "#E69F00"
ggplot(als_onset_trunk[als_onset_trunk$cluster=="Cluster3",], aes(x = onset_delta.x)) + 
  geom_density(fill = fill, alpha = 0.8) + ylab("Density")+ 
  theme(text = element_text(size=15), axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15), axis.title.x=element_blank())+
  xlim(0,1100) + ylim(0,0.005)
dev.off()

png("trunk_cluster3.png", width = 1500, height = 1500, res=300)
fill <- "#E69F00"
ggplot(als_onset_trunk[als_onset_trunk$cluster=="Cluster3",], aes(x = trunk_max)) + 
  geom_density(fill = fill, alpha = 0.8) + ylab("Density")+ 
  theme(text = element_text(size=15), axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15), axis.title.x=element_blank()) +
  xlim(-3,2) + ylim(0,0.8)
dev.off()

png("onset_delta_cluster4.png", width = 1500, height = 1500, res=300)
fill <- "#CC79A7"
ggplot(als_onset_trunk[als_onset_trunk$cluster=="Cluster4",], aes(x = onset_delta.x)) + 
  geom_density(fill = fill, alpha = 0.8) + ylab("Density")+ 
  theme(text = element_text(size=15), axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15), axis.title.x=element_blank()) +
  xlim(0,1100) + ylim(0,0.005)
dev.off()

png("trunk_cluster4.png", width = 1500, height = 1500, res=300)
fill <- "#CC79A7"
ggplot(als_onset_trunk[als_onset_trunk$cluster=="Cluster4",], aes(x = trunk_max)) + 
  geom_density(fill = fill, alpha = 0.8) + ylab("Density")+ 
  theme(text = element_text(size=15), axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15), axis.title.x=element_blank())+
  xlim(-3,2) + ylim(0,0.8)
dev.off()



############ ALSFRS_Total vs Delta (LOESS) ############
als_ALSFRS <- als3[als3$SubjectID %in% id_2424,] %>% filter(feature_name=="ALSFRS_Total") %>% 
  select(SubjectID,feature_value,feature_delta) %>% mutate(SubjectID=as.character(SubjectID),
                                                           feature_value=as.numeric(feature_value), feature_delta=as.numeric(feature_delta))
als_ALSFRS$cluster <- NA
als_ALSFRS[als_ALSFRS$SubjectID %in% c1_id, ]$cluster <- "Cluster1" 
als_ALSFRS[als_ALSFRS$SubjectID %in% c2_id, ]$cluster <- "Cluster2" 
als_ALSFRS[als_ALSFRS$SubjectID %in% c3_id, ]$cluster <- "Cluster3" 
als_ALSFRS[als_ALSFRS$SubjectID %in% c4_id, ]$cluster <- "Cluster4" 

png("ALSFRS_Total_trajectory(LOESS)_4cluster.png", width = 3000, height = 1800, res=300)
ggplot(als_ALSFRS, aes(feature_delta, feature_value, colour = cluster)) + geom_smooth(method = "loess")+
  xlab("Delta")+ylab("ALSFRS_Total")+ggtitle("ALSFRS_Total Trajectory (LOESS) for ALS Patients (4 Clusters)")
dev.off()

############ fvc_percent vs Delta (LOESS) ############
als_fvc_percent <- als3[als3$SubjectID %in% id_2424,] %>% filter(feature_name=="fvc_percent") %>% 
  select(SubjectID,feature_value,feature_delta) %>% mutate(SubjectID=as.character(SubjectID),
                                                           feature_value=as.numeric(feature_value), feature_delta=as.numeric(feature_delta)) %>% filter(feature_delta >= 0)

als_fvc_percent$cluster <- NA
als_fvc_percent[als_fvc_percent$SubjectID %in% c1_id, ]$cluster <- "Cluster1" 
als_fvc_percent[als_fvc_percent$SubjectID %in% c2_id, ]$cluster <- "Cluster2" 
als_fvc_percent[als_fvc_percent$SubjectID %in% c3_id, ]$cluster <- "Cluster3" 
als_fvc_percent[als_fvc_percent$SubjectID %in% c4_id, ]$cluster <- "Cluster4" 

png("fvc_percent_trajectory(LOESS)_4cluster.png", width = 3000, height = 1800, res=300)
ggplot(als_fvc_percent, aes(feature_delta, feature_value, colour = cluster)) + geom_smooth(method = "loess")+
  xlab("Delta")+ylab("fvc_percent")+ggtitle("fvc_percent Trajectory (LOESS) for ALS Patients (4 Clusters)")
dev.off()

#################################
## IGNORE THE CODE BLOCK BELOW ##
#################################

p1 <- ggplot(als_ALSFRS[als_ALSFRS$SubjectID %in% c1_id,], aes(feature_delta, feature_value, group=SubjectID))
p1 <- p1+geom_line()+xlab("Delta")+ylab("ALSFRS_Total")+ggtitle("ALSFRS_Total Trajectory for ALS Patients in Cluster 1")+ xlim(0, 650)

p2 <- ggplot(als_ALSFRS[als_ALSFRS$SubjectID %in% c2_id,], aes(feature_delta, feature_value, group=SubjectID))
p2 <- p2+geom_line()+xlab("Delta")+ylab("ALSFRS_Total")+ggtitle("ALSFRS_Total Trajectory for ALS Patients in Cluster 2")+ xlim(0, 650)

p3 <- ggplot(als_ALSFRS[als_ALSFRS$SubjectID %in% c3_id,], aes(feature_delta, feature_value, group=SubjectID))
p3 <- p3+geom_line()+xlab("Delta")+ylab("ALSFRS_Total")+ggtitle("ALSFRS_Total Trajectory for ALS Patients in Cluster 3")+ xlim(0, 650)

p4 <- ggplot(als_ALSFRS[als_ALSFRS$SubjectID %in% c4_id,], aes(feature_delta, feature_value, group=SubjectID))
p4 <- p4+geom_line()+xlab("Delta")+ylab("ALSFRS_Total")+ggtitle("ALSFRS_Total Trajectory for ALS Patients in Cluster 4")+ xlim(0, 650)

png("ALSFRS_Total_trajectory_cluster1.png", width = 1000, height = 500)
p1
dev.off()
png("ALSFRS_Total_trajectory_cluster2.png", width = 1000, height = 500)
p2
dev.off()
png("ALSFRS_Total_trajectory_cluster3.png", width = 1000, height = 500)
p3
dev.off()
png("ALSFRS_Total_trajectory_cluster4.png", width = 1000, height = 500)
p4
dev.off()


############ ALSFRS_Total vs Delta (LOESS) ############
png("LOESS_ALSFRS_Total_trajectory_cluster1.png", width = 1000, height = 500)
g1 <- ggplot(als_ALSFRS[als_ALSFRS$SubjectID %in% c1_id,], aes(feature_delta, feature_value)) + geom_point() + 
  geom_smooth(method = "loess") + xlab("Delta") + ylab("ALSFRS_Total") + ggtitle("LOESS: ALSFRS_Total Trajectory for ALS Patients in Cluster 1")+ xlim(0, 650)
g1
dev.off()

png("LOESS_ALSFRS_Total_trajectory_cluster2.png", width = 1000, height = 500)
g2 <- ggplot(als_ALSFRS[als_ALSFRS$SubjectID %in% c2_id,], aes(feature_delta, feature_value)) + geom_point() + 
  geom_smooth(method = "loess") + xlab("Delta") + ylab("ALSFRS_Total") + ggtitle("LOESS: ALSFRS_Total Trajectory for ALS Patients in Cluster 2")+ xlim(0, 650)
g2
dev.off()

png("LOESS_ALSFRS_Total_trajectory_cluster3.png", width = 1000, height = 500)
g3 <- ggplot(als_ALSFRS[als_ALSFRS$SubjectID %in% c3_id,], aes(feature_delta, feature_value)) + geom_point() + 
  geom_smooth(method = "loess") + xlab("Delta") + ylab("ALSFRS_Total") + ggtitle("LOESS: ALSFRS_Total Trajectory for ALS Patients in Cluster 3")+ xlim(0, 650)
g3
dev.off()

png("LOESS_ALSFRS_Total_trajectory_cluster4.png", width = 1000, height = 500)
g4 <- ggplot(als_ALSFRS[als_ALSFRS$SubjectID %in% c4_id,], aes(feature_delta, feature_value)) + geom_point() + 
  geom_smooth(method = "loess") + xlab("Delta") + ylab("ALSFRS_Total") + ggtitle("LOESS: ALSFRS_Total Trajectory for ALS Patients in Cluster 4")+ xlim(0, 650)
g4
dev.off()

############ fvc_percent vs Delta ############
als_fvc_percent <- als3[als3$SubjectID %in% id_2424,] %>% filter(feature_name=="fvc_percent") %>% 
  select(SubjectID,feature_value,feature_delta) %>% mutate(SubjectID=as.character(SubjectID),
                                                           feature_value=as.numeric(feature_value), feature_delta=as.numeric(feature_delta)) %>% filter(feature_delta >= 0)

p5 <- ggplot(als_fvc_percent[als_fvc_percent$SubjectID %in% c1_id,], aes(feature_delta, feature_value, group=SubjectID))
p5 <- p5+geom_line()+xlab("Delta")+ylab("fvc_percent")+ggtitle("fvc_percent Trajectory for ALS Patients in Cluster 1")+xlim(0, 650)+ylim(0,175)

p6 <- ggplot(als_fvc_percent[als_fvc_percent$SubjectID %in% c2_id,], aes(feature_delta, feature_value, group=SubjectID))+xlim(0, 650)+ylim(0,175)
p6 <- p6+geom_line()+xlab("Delta")+ylab("fvc_percent")+ggtitle("fvc_percent Trajectory for ALS Patients in Cluster 2")

p7 <- ggplot(als_fvc_percent[als_fvc_percent$SubjectID %in% c3_id,], aes(feature_delta, feature_value, group=SubjectID))+xlim(0, 650)+ylim(0,175)
p7 <- p7+geom_line()+xlab("Delta")+ylab("fvc_percent")+ggtitle("fvc_percent Trajectory for ALS Patients in Cluster 3")

p8 <- ggplot(als_fvc_percent[als_fvc_percent$SubjectID %in% c4_id,], aes(feature_delta, feature_value, group=SubjectID))+xlim(0, 650)+ylim(0,175)
p8 <- p8+geom_line()+xlab("Delta")+ylab("fvc_percent")+ggtitle("fvc_percent Trajectory for ALS Patients in Cluster 4")

png("fvc_percent_trajectory_cluster1.png", width = 1000, height = 500)
p5
dev.off()
png("fvc_percent_trajectory_cluster2.png", width = 1000, height = 500)
p6
dev.off()
png("fvc_percent_trajectory_cluster3.png", width = 1000, height = 500)
p7
dev.off()
png("fvc_percent_trajectory_cluster4.png", width = 1000, height = 500)
p8
dev.off()

############ fvc_percent vs Delta (LOESS) ############
png("LOESS_fvc_percent_trajectory_cluster1.png", width = 1000, height = 500)
g5 <- ggplot(als_fvc_percent[als_fvc_percent$SubjectID %in% c1_id,], aes(feature_delta, feature_value)) + geom_point() +xlim(0, 650)+ylim(0,175)+ 
  geom_smooth(method = "loess") + xlab("Delta") + ylab("fvc_percent") + ggtitle("LOESS: fvc_percent Trajectory for ALS Patients in Cluster 1")
g5
dev.off()

png("LOESS_fvc_percent_trajectory_cluster2.png", width = 1000, height = 500)
g6 <- ggplot(als_fvc_percent[als_fvc_percent$SubjectID %in% c2_id,], aes(feature_delta, feature_value)) + geom_point() +xlim(0, 650)+ylim(0,175)+ 
  geom_smooth(method = "loess") + xlab("Delta") + ylab("fvc_percent") + ggtitle("LOESS: fvc_percent Trajectory for ALS Patients in Cluster 2")
g6
dev.off()

png("LOESS_fvc_percent_trajectory_cluster3.png", width = 1000, height = 500)
g7 <- ggplot(als_fvc_percent[als_fvc_percent$SubjectID %in% c3_id,], aes(feature_delta, feature_value)) + geom_point() +xlim(0, 650)+ylim(0,175)+ 
  geom_smooth(method = "loess") + xlab("Delta") + ylab("fvc_percent") + ggtitle("LOESS: fvc_percent Trajectory for ALS Patients in Cluster 3")
g7
dev.off()

png("LOESS_fvc_percent_trajectory_cluster4.png", width = 1000, height = 500)
g8 <- ggplot(als_fvc_percent[als_fvc_percent$SubjectID %in% c4_id,], aes(feature_delta, feature_value)) + geom_point() +xlim(0, 650)+ylim(0,175)+ 
  geom_smooth(method = "loess") + xlab("Delta") + ylab("fvc_percent") + ggtitle("LOESS: fvc_percent Trajectory for ALS Patients in Cluster 4")
g8
dev.off() 
  
  
  
  
  
  
  

  
  