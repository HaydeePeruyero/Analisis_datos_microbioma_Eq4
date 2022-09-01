################################################################################ 
###Chapter 12: Modeling Zero-Inflated Microbiome Data 
###Yinglin Xia: September, 2018                                                  
################################################################################ 

################################################################################ 
###12.2. Zero-inflated Models: ZIP and ZINB 
################################################################################ 
setwd("E:/Home/MicrobiomeStatR/Analysis")

abund_table=read.csv("allTD_long.csv",header=TRUE)
head(abund_table)
tail(abund_table)

library(dplyr)
abund_table_28 <- filter(abund_table, Ind == 28)
head(abund_table_28)

abund_table_28$x<- with(abund_table_28,ifelse(as.factor(DX)%in% "0",0, 1))
names(abund_table_28)

abund_table_28$fx <- factor(abund_table_28$x)
names(abund_table_28)

I=is.na(abund_table_28$Y) | is.na(abund_table_28$fx)|is.na(abund_table_28$nReads)
abund_table_28a <- abund_table_28[!I,]

par(mfrow=c(1,2))
plot(table(abund_table_28a$Y),ylab="Frequencies",main="Lactobacillus.vaginalis",
     xlab="Observed read values")
plot(sort(abund_table_28a$Y),ylab="Frequencies",main="Lactobacillus.vaginalis",
     xlab="Observed read values")

abund_table_28a$Offset <- log(abund_table_28a$nReads);
head(abund_table_28a$Offset)

f28 <- formula(Y ~  fx + offset(Offset)|1)
f28 <- formula(Y ~  fx + offset(Offset)|fx)

library(pscl)
ZIP28 <- zeroinfl(formula = f28, dist = "poisson", link = "logit", data = abund_table_28a)
summary(ZIP28)

ZINB28 <- zeroinfl(formula = f28, dist = "negbin", link = "logit", data = abund_table_28a)
summary(ZINB28)


##12.3.3 Modeling using ZHP and ZHNB 

## 12.3.3.2 Conceptual Adjustment for using Zero-hurdle Models
ZHP28 <- hurdle(formula = f28, dist= "poisson", data = abund_table_28a)
summary(ZHP28)

ZHNB28 <- hurdle(formula = f28, dist= "negbin", data = abund_table_28a)
summary(ZHNB28)

    
##12.3.4 Comparing Zero-inflated and Zero-Hurdle Models
##12.3.4.1 Using Likelihood Ratio Test
library(lmtest)

lrtest(ZIP28,ZINB28)
lrtest(ZHP28,ZHNB28)

##12.3.4.2 Using AIC
#lower is better
t(AIC(ZIP28, ZINB28, ZHP28, ZHNB28))

##12.3.4.3 Using BIC
library(nonnest2)

#lower is better
icci(ZIP28,ZINB28)
icci(ZHP28,ZHNB28)

##12.3.4.4 Using Vuong Test
library(pscl)

# compare ZIP vs. ZHP
vuong(ZIP28, ZHP28)

# compare ZIP vs. ZHNB
vuong(ZIP28, ZHNB28)

# compare ZINB vs. ZHP
vuong(ZINB28, ZHP28)

# compare ZINB vs. ZHNB
vuong(ZINB28, ZHNB28)

# compare ZHNB vs. ZHP
vuong(ZHNB28, ZHP28)

  
##12.3.5 Interpreting Main Effects of Modeling Results 
ZINB28 <- zeroinfl(formula = f28, dist = "negbin", link = "logit", data = abund_table_28)
ZHNB28 <- hurdle(formula = f28, dist= "negbin", link="logit",zero.dist="binomial", data = abund_table_28)
summary(ZINB28)
summary(ZHNB28)

expZINB28Coef <- exp(coef((ZINB28)))
expZINB28Coef <- matrix(expZINB28Coef, ncol = 2)
expZINB28Coef 
colnames(expZINB28Coef) <- c("Count_model","Zero_inflation_model")
expZINB28Coef

expZHNB28Coef <- exp(coef((ZHNB28)))
expZHNB28Coef <- matrix(expZHNB28Coef, ncol = 2)
colnames(expZHNB28Coef) <- c("Count_model","Zero_hurdle_model")
expZHNB28Coef

################################################################################
###12.4. Zero-inflated Beta Regression Model with Random Effects
################################################################################

##12.4.4.2 Step-by-Step Implementing ZIBR 
install.packages("devtools")
devtools::install_github("chvlyl/ZIBR")
setwd("E:/Home/MicrobiomeStatR/Analysis")library(ZIBR)

rawfile<- "https://raw.githubusercontent.com/chvlyl/PLEASE/master/
1_Data/Raw_Data/MetaPhlAn/PLEASE/G_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result.xls"
taxa_table <- read.table(rawfile,sep='\t',header=TRUE,row.names = 1,
                         check.names=FALSE,stringsAsFactors=FALSE)
head(taxa_table)

taxa_table_t <- t(taxa_table)
head(taxa_table_t)                         
                         
cat('samples','taxa',dim(taxa_table_t),'\n')
taxa_table_t[1:3,1:3]

totalreadfile <- 'https://raw.githubusercontent.com/chvlyl/PLEASE/master/1_Data/
Raw_Data/MetaPhlAn/Human_Reads/please_combo_human_reads.xls'
abund_table <- read.table(totalreadfile,sep='\t',header=TRUE,
                         row.names=1,stringsAsFactors=FALSE)
head(abund_table)

samplefile <- 'https://raw.githubusercontent.com/chvlyl/PLEASE/master/1_Data/
Processed_Data/Sample_Information/2015_02_13_Processed_Sample_Information.csv'
meta_table <- read.csv(samplefile,row.names=1)
head(meta_table)
tail(meta_table)

meta_table[,'Treatment.Specific']
low_samples <- subset(abund_table,NonHumanReads<10000)
low_samples[,1:5]

rownames(taxa_table_t)[which(rownames(taxa_table_t) %in% rownames(low_samples))]

taxa_table_t <- taxa_table_t[-which(rownames(taxa_table_t) %in% rownames(low_samples)),]


### Before deletion
dim(taxa_table_t) 
### After deletion
dim(taxa_table_t)

filter_index1 <- apply(taxa_table_t,2,function(X){sum(X>0)>0.4*length(X)})
filter_index2 <- apply(taxa_table_t,2,function(X){quantile(X,0.9)>1})
taxa_filter <- taxa_table_t[,filter_index1 & filter_index2]
taxa_filter

taxa_filter <- 100*sweep(taxa_filter, 1, rowSums(taxa_filter), FUN="/")
head(taxa_filter)

cat('after filter:','samples','taxa',dim(taxa_filter),'\n')
cat(colnames(taxa_filter),'\n')

head(rowSums(taxa_filter))

taxa_data <- taxa_filter
head(taxa_data)

dim(taxa_data)

library(dplyr)
# create covariates:Time, antiTNF+EEN
reg_cov <-
  data.frame(Sample=rownames(taxa_data),stringsAsFactors = FALSE) %>% 
  left_join(add_rownames(meta_table,var = 'Sample'),by='Sample')%>%
  # exclude PEN, just keep antiTNF and EEN
  dplyr::filter(Treatment.Specific!='PEN') %>% 
  # subset meta table
  dplyr::select(Sample,Time,Subject,Response,Treatment.Specific) %>%# 
  group_by(Subject) %>% summarise(count = n()) %>% dplyr::filter(count==4) %>%
  dplyr::select(Subject) %>%
  left_join(add_rownames(meta_table,var = 'Sample'),by='Subject') %>%
  # create treatment variable Treat and code antiTNF as 1, EEN as 0
  mutate(Treat=ifelse(Treatment.Specific=='antiTNF',1,0)) %>%
  dplyr::select(Sample,Subject,Time,Response,Treat) %>%
  dplyr::mutate(Subject=paste('S',Subject,sep='')) %>%
  # recode Time variable
  dplyr::mutate(Time=ifelse(Time=='1',0,ifelse(Time=='2',1,ifelse(Time=='3',4,ifelse(Time=='4',8,NA))))) %>%
  # create Time by Treatment interaction term
  dplyr::mutate(Time.X.Treatment=Time*Treat) %>%
  as.data.frame

# take out first time point
reg_cov_t1   <-  subset(reg_cov,Time==0)
rownames(reg_cov_t1) <- reg_cov_t1$Subject
reg_cov_t234 <-  subset(reg_cov,Time!=0)
reg_cov_t234 <- data.frame(
  baseline_sample=reg_cov_t1[reg_cov_t234$Subject,'Sample'],
  baseline_subject=reg_cov_t1[reg_cov_t234$Subject,'Subject'],
  reg_cov_t234,
  stringsAsFactors = FALSE)

head(reg_cov_t234)

library(Matrix)
library(ZIBR)
# create a matrix to hold taxa
taxa_all <- colnames(taxa_data)
taxa_all
# create a list to hold p-values
p_taxa_list_zibr <- list()
p_taxa_list_zibr
#ZIBR independently fit each taxon. 
for (taxa in taxa_all){
  #for example,taxa = "g__Bacteroides"      
  ###create covariates
  X <- data.frame(
    Baseline=taxa_data[reg_cov_t234$baseline_sample, taxa]/100,
    reg_cov_t234[,c('Time','Treat')]
  )
  rownames(X) <- reg_cov_t234$Sample
  Z <- X
  subject_ind <- reg_cov_t234$Subject
  time_ind   <- reg_cov_t234$Time
  ###create a table to summarize statistics
  cat(taxa,'\n')
  Y <- taxa_data[reg_cov_t234$Sample, taxa]/100
  cat('Zeros/All',sum(Y==0),'/',length(Y),'\n')
  if (sum(Y>0)<10 | sum(Y==0) <10 | max(Y)<0.01){
    print('skip')
    next
  }else{
    est <- zibr(logistic.cov=X,beta.cov=Z,Y=Y,
                subject.ind=subject_ind,
                time.ind=time_ind,
                quad.n=30,verbose=TRUE)
    p_taxa_list_zibr[[taxa]] <- est$joint.p
    
  }
  
}


# unadjusted p values
p_taxa_zibr <- t(as.data.frame(p_taxa_list_zibr))
p_taxa_zibr

# adjusted  p values
library(dplyr)
p_taxa_zibr_adj <-
  add_rownames(as.data.frame(p_taxa_zibr),var = 'Taxa') %>% mutate_each(funs(p.adjust(.,'fdr')),-Taxa)
p_taxa_zibr_adj

# make the table
library(xtable)
table <-xtable(p_taxa_zibr_adj, caption="Table of significant taxa", digits=3,
               label="sig_taxa_table")  
print.xtable(table, type="html", file="IBD_Table.html")  
write.csv(p_taxa_zibr_adj,file=paste('Results_antiTNF_EEN_ZIBR.csv',sep=''))

