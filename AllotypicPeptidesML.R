library(dplyr)
library(reshape2)
library(caret)
library(kernlab)
library(RSQLite)
library(rPython)
require(pracma)


#Load sqlite db as dataframe
options(stringsAsFactors = F)
con = dbConnect(SQLite(), dbname='Projects/TumorSequencing/DB/TueDB.sqlite')
myQuery <- dbSendQuery(con, "SELECT * FROM Class1")
df <- dbFetch(myQuery, n = -1)

###Load and clean database
#df=read.csv("Projects/TumorSequencing/DB/TueDB_ClassI.csv", sep = ",")
df = df[c('Sequence','HLA','HLA_A1','HLA_A2','HLA_B1','HLA_B2','HLA_C1','HLA_C2','Dignity','Patient..Donor')]
df['HLA_A1'] <- sapply(df[,'HLA_A1'], function(x){gsub(x=strsplit(x=x, split = ":", fixed = T)[[1]][1], pattern = "*", replacement = "", fixed = T)})
df['HLA_A2'] <- sapply(df[,'HLA_A2'], function(x){gsub(x=strsplit(x=x, split = ":", fixed = T)[[1]][1], pattern = "*", replacement = "", fixed = T)})
df['HLA_B1'] <- sapply(df[,'HLA_B1'], function(x){gsub(x=strsplit(x=x, split = ":", fixed = T)[[1]][1], pattern = "*", replacement = "", fixed = T)})
df['HLA_B2'] <- sapply(df[,'HLA_B2'], function(x){gsub(x=strsplit(x=x, split = ":", fixed = T)[[1]][1], pattern = "*", replacement = "", fixed = T)})
df['HLA_C1'] <- sapply(df[,'HLA_C1'], function(x){gsub(x=strsplit(x=x, split = ":", fixed = T)[[1]][1], pattern = "*", replacement = "", fixed = T)})
df['HLA_C2'] <- sapply(df[,'HLA_C2'], function(x){gsub(x=strsplit(x=x, split = ":", fixed = T)[[1]][1], pattern = "*", replacement = "", fixed = T)})
df[,3:8] <- apply(df[,3:8],1:2, function(allele){
  if(is.na(allele)){
    return(NA)
  }else if(allele == ""){
    return(NA)
  }else if(allele == "nd"){
    return(NA)
  }else{
    return(allele)
  }
})
df[which(df[,'HLA_A2']=="B07"),] <- NA
df[which(df[,'HLA_B1']=="B"|df[,'HLA_B2']=="B"),] <- NA
df[which(df[,'HLA_B1']=="C07"|df[,'HLA_B2']=="C07"),] <- NA
df[which(df[,'HLA_B1']=="A66"|df[,'HLA_B2']=="A66"),] <- NA
df=na.omit(df)
df$Dignity <- paste(df$Patient..Donor,df$Dignity)
df['Sequence']<-sapply(df['Sequence'], toupper)


###Define Allotypes to differenciate
allele_count_threshold<-9

A_Allotypes=unique(c(df$HLA_A1,df$HLA_A2))
HLA_A_allotype_sample_count<-data.frame(table(data.frame(cbind(c(df$Patient..Donor,df$Patient..Donor),c(df$HLA_A1,df$HLA_A2)))))
HLA_A_allotype_count=list()
for (allele in A_Allotypes) {
  HLA_A_allotype_count[[allele]]<-length(HLA_A_allotype_sample_count[which(HLA_A_allotype_sample_count$X2==allele & !HLA_A_allotype_sample_count$Freq==0),]$Freq)
}
A_Allotypes<-names(which(HLA_A_allotype_count>allele_count_threshold))

HLA_B_allotype_sample_count<-data.frame(table(data.frame(cbind(c(df$Patient..Donor,df$Patient..Donor),c(df$HLA_B1,df$HLA_B2)))))
B_Allotypes=unique(c(df$HLA_B1,df$HLA_B2))
HLA_B_allotype_count=list()
for (allele in B_Allotypes) {
  HLA_B_allotype_count[[allele]]<-length(HLA_B_allotype_sample_count[which(HLA_B_allotype_sample_count$X2==allele & !HLA_B_allotype_sample_count$Freq==0),]$Freq)
}
B_Allotypes<-names(which(HLA_B_allotype_count>allele_count_threshold))

HLA_C_allotype_sample_count<-data.frame(table(data.frame(cbind(c(df$Patient..Donor,df$Patient..Donor),c(df$HLA_C1,df$HLA_C2)))))
C_Allotypes=unique(c(df$HLA_C1,df$HLA_C2))
HLA_C_allotype_count=list()
for (allele in C_Allotypes) {
  HLA_C_allotype_count[[allele]]<-length(HLA_C_allotype_sample_count[which(HLA_C_allotype_sample_count$X2==allele & !HLA_C_allotype_sample_count$Freq==0),]$Freq)
}
C_Allotypes<-names(which(HLA_C_allotype_count>allele_count_threshold))

All_Allotypes<-c(A_Allotypes,B_Allotypes,C_Allotypes)
All_Allotypes_count<-c(unlist(HLA_A_allotype_count),unlist(HLA_B_allotype_count),unlist(HLA_C_allotype_count))

###Generate Data_Matrix (Rows SampleID, Columns Peptides, Fill by Counts), Multilabel matrix (Rows SampleID, Columns HLA Types, Fill binary)
peptide_count_threshold <- 9

Sequences_list<-list()
Typing_list<-list()
total_sample_count=table(df['Dignity'])
HLA_molten=melt(df,id='Dignity',measure.vars=c('HLA_A1','HLA_A2','HLA_B1','HLA_B2','HLA_C1','HLA_C2'))
total_peptide_count<-table(df['Sequence'])
total_peptide_count<-total_peptide_count[which(total_peptide_count>peptide_count_threshold)]

###initialize zeros data_matrix: peptides x samples
data_matrix=data.frame(matrix(0, ncol = length(total_peptide_count), nrow= length(total_sample_count)))
row.names(data_matrix) <- names(total_sample_count)
colnames(data_matrix) <- names(total_peptide_count)

###initialize zeros label_matrix: alleles x samples
label_matrix=data.frame(matrix(0, ncol = length(All_Allotypes), nrow= length(total_sample_count)))
row.names(label_matrix) <- names(total_sample_count)
colnames(label_matrix) <- All_Allotypes

###fill data and label matrices iterating over samples
counter <- 0
for (sample in names(total_sample_count)){
  if(counter %% 10 == 0){
    print(paste0(counter, "/", length(total_sample_count)))
  }
  
  sample_peptides <- (df[which(df[,'Dignity']==sample), "Sequence"])
  data_matrix[sample, which(colnames(data_matrix) %in% sample_peptides)] <- 1
  sample_alleles=unique(HLA_molten[which(HLA_molten[,'Dignity']==sample), "value"])
  label_matrix[sample, which(colnames(label_matrix) %in% sample_alleles)] <- 1
  counter <- counter + 1
}

###Run cross validation and generate top_n peptide lists
shuffle<-sample(rep(1:length(label_matrix$A02), each=1))
data_matrix=data_matrix[shuffle,]
label_matrix=label_matrix[shuffle,]
rownames(label_matrix)=rownames(data_matrix)
auc_list <- list()

fitControl <- trainControl(method = "repeatedcv",
                           ## 10-fold CV...
                           number = 10,
                           classProbs = TRUE,
                           ## repeated ten times
                           repeats = 10,
                           summaryFunction = twoClassSummary)


for (allele in All_Allotypes){
  ###caret, kernlab packages Train random forrest on selected peptides
  TrainData <- as.data.frame(data_matrix[,names(total_peptide_count)])
  TrainClasses <- as.character(label_matrix[[allele]])
  TrainClasses <- sapply(TrainClasses,function(x) {gsub(x,pattern = "1", replacement = "ClassI")})
  TrainClasses <- as.factor(as.vector(sapply(TrainClasses,function(x) {gsub(x,pattern = "0", replacement = "Class2")})))
  Training <- train(x=TrainData, y=TrainClasses, 
                 method = "cforest", 
                 preProc = c("center"),
                 trControl = fitControl, metric = "ROC")
  auc_list[[allele]]<-Training
}


