################################################################################ 
###Chapter 6: Community Diversity Measures and Calculations                      
###Yinglin Xia: September, 2018                                                  
################################################################################ 

#Set working directory  
setwd("E:/Home/MicrobiomeStatUsingR/Analysis/") 

################################################################################ 
###6.3 Alpha Diversity Measures and Calculations                                                   
################################################################################

##6.3.1 Chao 1 Richness Index and Number of Taxa
options(width=65,digits=4)
abund_table=read.csv("VdrGenusCounts.csv",row.names=1,check.names=FALSE)
library(vegan)

head(abund_table) 
abund_table<-t(abund_table)
head(abund_table)

num_genera <- specnumber(abund_table) 
num_genera

index=estimateR(abund_table) 
index
head(index)

chao1_genus=estimateR(abund_table)[2,] 
chao1_genus

##6.3.2 Shannon-Wiener Diversity Index
shannon_genus <- diversity(abund_table,index="shannon",MARGIN=1)
shannon_genus

shannon_genus <- diversity(abund_table,MARGIN=1)
shannon_genus

shannon_genus <- diversity(abund_table)
shannon_genus

#use decostand to convert data into proportions
abund_table_total<-decostand(abund_table, MARGIN=1, method="total")

#multiply that matrix X a natural log transformed matrix - p*ln(p)
abund_table_p_lnp<-abund_table_total*log(abund_table_total)
#sum values by sample and multiply by -1
rowSums(abund_table_p_lnp,na.rm=TRUE)*-1

##6.3.3 Simpson Diversity Index
simp_genus <- diversity(abund_table, "simpson")
simp_genus

#using plain R functions
#use decostand to convert data into proportions
abund_table_total<-decostand(abund_table, MARGIN=1, method="total")
#square the proportions
abund_table_total_p2<-abund_table_total^2
#get the row sums
1-rowSums(abund_table_total_p2, na.rm=TRUE)

inv_simp <- diversity(abund_table, "invsimpson")
inv_simp


##6.3.4 Pielou's Evenness Index
#using specnumber and diversity functions
S <- specnumber(abund_table)
H<-diversity(abund_table, "shannon") 
J <- H/log(S)
J 

##6.3.5 Make a Dataframe of Diversity Indices
#make a data frame of number of genera
N <- specnumber(abund_table) 
df_N <-data.frame(sample=names(N),value=N,measure=rep("Number",length(N))) 

#make a data frame of Chao1 richness
CH=estimateR(abund_table)[2,] 
df_CH <-data.frame(sample=names(CH),value=CH,measure=rep("Chao1",length(CH))) 

#make a data frame of Shannon evenness
H<-diversity(abund_table, "shannon") 
df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon",length(H)))

#make a data frame of Simpson index
df_simp<-data.frame(sample=names(simp_genus),value=simp_genus,measure=rep("Simpson",length(simp_genus)))

#make a data frame of Pielou index
df_J<-data.frame(sample=names(J),value=J,measure=rep("Pielou",length(J)))

df<-rbind(df_N,df_CH,df_H,df_simp,df_J)
rownames(df)<-NULL
df


################################################################################ 
###6.4 Beta Diversity Measures and Calculation                                   
################################################################################ 

library(vegan)
library(BiodiversityR)
betadiver(help=TRUE)

##6.4.1 Binary Similarity Coefficients: Jaccard and Sørensen Indices
jaccard<-vegdist(abund_table, "jaccard")
jaccard

Sørensen<-vegdist(abund_table,binary=TRUE)
Sørensen

##6.4.2 Distance (dissimilarity) Coefficients: Bray-Curtis Index
bray<-vegdist(abund_table, "bray")
bray
