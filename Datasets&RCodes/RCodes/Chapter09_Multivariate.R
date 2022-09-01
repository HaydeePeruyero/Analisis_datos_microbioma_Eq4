################################################################################ 
###Chapter 9: Multivariate Community Analysis 
###Yinglin Xia: September, 2018                                                  
################################################################################ 

################################################################################ 
##9.1 Hypothesis Testing among Groups using Permutational Multivariate Analysis 
###   of Variance (PERMANOVA)
################################################################################ 

setwd("E:/Home/MicrobiomeStatR/Analysis")

##9.1.2 Implementing PERMANOVA using vegan Package
abund_table=read.csv("VdrGenusCounts.csv",row.names=1,check.names=FALSE)

abund_table<-t(abund_table)

library(vegan)
bray<-vegdist(abund_table, "bray")
jaccard<-vegdist(abund_table, "jaccard")
Sørensen<-vegdist(abund_table,binary=TRUE)

grouping<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
grouping$Location <- with(grouping, ifelse(X3%in%"drySt-28F", "Fecal", "Cecal"))
grouping$Group <- with(grouping,ifelse(as.factor(X2)%in% c(11,12,13,14,15),c("Vdr-/-"), c("WT")))
grouping <- grouping[,c(4,5)]
grouping 

set.seed(123)
#adonis="analysis of dissimilarity"  
adonis(bray ~ Group,data=grouping,permutations = 1000) 
adonis(abund_table ~ Group,data=grouping,permutations = 1000, method = "bray") 

adonis(jaccard ~ Group,data=grouping,permutations = 1000) 
adonis(Sørensen ~ Group,data=grouping,permutations = 1000) 

adonis(bray ~ Location,data=grouping,permutations = 1000) 
adonis(jaccard ~ Location,data=grouping,permutations = 1000)
adonis(Sørensen ~ Location,data=grouping,permutations = 1000)

adonis(bray ~ Group*Location,data=grouping,permutations = 1000) 
adonis(bray ~ Location*Group,data=grouping,permutations = 1000) 

grouping$Group4<- with(grouping, interaction(Location,Group))
grouping

adonis(bray ~ Group4,data=grouping,permutations = 1000) 
adonis(bray ~ Group4,data=grouping,permutations = 1000,contr.unordered = "contr.sum") 
adonis(bray ~ Group4,data=grouping,permutations = 1000,contr.unordered = "contr.sum",contr.ordered = "contr.poly") 

adonis(bray ~ grouping$Group4,permutations = 1000) 
adonis(abund_table ~ Group4,data=grouping,permutations = 1000, method = "bray")

##9.1.3 Implementing Pairwise Permutational MANOVA using RVAideMemoire Package
install.packages("RVAideMemoire")
library(RVAideMemoire)

set.seed(0)
pairwise.perm.manova(bray,grouping$Group4,nperm=1000)

# or
pairwise.perm.manova(vegdist(abund_table,"bray"),grouping$Group4,nperm=1000)
  
pairwise.perm.manova(bray, grouping$Group4, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"), nperm = 1000, 
                     progress = TRUE, p.method = "fdr")

pairwise.perm.manova(bray, grouping$Group4, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"), nperm = 1000, 
                     progress = TRUE, p.method = "none")

pairwise.perm.manova(bray, grouping$Group4, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"), nperm = 1000, 
                     progress = TRUE, p.method = "bonferroni")

pairwise.perm.manova(bray, grouping$Group4, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"), nperm = 1000, 
                     progress = TRUE, p.method = "holm")

pairwise.perm.manova(bray, grouping$Group4, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"), nperm = 1000, 
                     progress = TRUE, p.method = "hochberg")

pairwise.perm.manova(bray, grouping$Group4, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"), nperm = 1000, 
                     progress = TRUE, p.method = "hommel")

pairwise.perm.manova(bray, grouping$Group4, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"), nperm = 1000, 
                     progress = TRUE, p.method = "BH")

pairwise.perm.manova(bray, grouping$Group4, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"), nperm = 1000, 
                     progress = TRUE, p.method = "fdr")

pairwise.perm.manova(bray, grouping$Group4, test = c("Pillai", "Wilks","Hotelling-Lawley", "Roy", "Spherical"), nperm = 1000, 
                     progress = TRUE, p.method = "BY")


## 9.1.4 Test Group Homogeneities using the Function betadisper() 
homo <-with(grouping,betadisper(bray, Group4))
homo

plot(homo)  
boxplot(homo)

anova(homo)
permutest(homo)
TukeyHSD(homo)



################################################################################ 
###9.2 Hypothesis Tests  among Group-Differences using Mantel's Test (MANTEL)
################################################################################ 

##9.2.2 Illustrating Mantel Test using vegan Package
##9.2.2.1 Test the Correlation of two Community Distance Matrices

library(GUniFrac)
data(throat.otu.tab)
otu_table <-throat.otu.tab 
data(throat.meta) 
data(throat.tree)

library(dplyr)
throat_meta <- select(throat.meta, SmokingStatus, Age, Sex, PackYears)

group <-select(throat.meta, SmokingStatus)
group$Status <- with(group, ifelse(SmokingStatus%in%"Smoker", 1, 0))

meta <- select(throat.meta, Age,Sex,PackYears)
meta$Gender <- with(meta,ifelse(Sex%in%"Male", 1, 0))
env <-select(meta, Age,Gender,PackYears)

library(vegan)
bray<-vegdist(otu_table, "bray")
jaccard<-vegdist(otu_table, "jaccard")
Sørensen<-vegdist(otu_table,binary=TRUE)

mantel(bray, jaccard,"pearson",permutations=1000) 
plot(bray, jaccard, main="Scatter plot of Bray-Curtis index vs. Jaccard index")

mantel(bray, Sørensen,"spearman",permutations=1000)
mantel(jaccard,Sørensen,"kendall",permutations=1000) 

##9.2.2.2 Test the Correlation of a Community Distance Matrix and a Design Distance Matrix
library(dplyr)
group <-select(throat.meta, SmokingStatus)
group$Status <- with(group, ifelse(SmokingStatus%in%"Smoker", 1, 0))
group <- group[,-1]
group_dist <- vegdist(scale(group), "euclid")

mantel(group_dist,bray,"pearson",permutations = 1000)

##9.2.2.3 Partial Mantel Test the Correlation of Two Distance Matrices Controlling the Third Matrix
throat_meta <- select(throat.meta, SmokingStatus, Age, Sex, PackYears)

meta <- select(throat.meta, Age,Sex,PackYears)
meta$Gender <- with(meta,ifelse(Sex%in%"Male", 1, 0))
env <-select(meta, Age,Gender,PackYears)

env_dist <- vegdist(scale(env), "euclid")

mantel.partial(group_dist,bray,env_dist, method = "pearson", permutations = 1000)

################################################################################ 
###9.3 Hypothesis Tests  among-Group Differences using Analysis of Similarity 
###    (ANOSIM)
################################################################################

##9.3.2 Illustrating Analysis of Similarity (ANOSIM) using vegan Package
abund_table=read.csv("VdrFecalGenusCounts.csv",row.names=1,check.names=FALSE)

abund_table<-t(abund_table)

grouping<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
grouping$Group <- with(grouping,ifelse(as.factor(X2)%in% c(11,12,13,14,15),c("Vdr-/-"), c("WT")))
grouping<- grouping[,c(4)]

set.seed(0)
library(vegan)
bray<-vegdist(abund_table, "bray")
anosim(bray, grouping,permutations = 1000)

anosim(abund_table, grouping, permutations = 1000, distance = "bray")

fit <- anosim(bray, grouping,permutations = 1000)
summary(fit)

plot(fit)

anosim(abund_table, grouping, permutations = 1000, distance = "jaccard")

fit_S <- anosim(Sørensen, grouping, permutations = 1000)
summary(fit_S) 
################################################################################ 
###9.4 Hypothesis Tests of Multi-Response Permutation Procedures (MRPP)
################################################################################ 

##9.4.2 Illustrating MRPP Using Vegan Package

mrpp(bray, grouping,permutations = 1000)

mrpp(abund_table, grouping, permutations = 1000, distance = "bray")
meandist(bray, grouping,permutations = 1000)

bray_mrpp <- meandist(bray, grouping,permutations = 1000,distance = "bray",weight.type = 1)
summary(bray_mrpp)


################################################################################ 
###9.5 Compare Microbiome Communities using the GUniFrac Package 
################################################################################ 

##9.5. 3 Comparing Microbiome Communities Using GUniFrac Package

setwd("E:/Home/MicrobiomeStatUsingR/Analysis/") 
otu_tab <- read.table("td_OTU_tag_mapped_lineage.txt", header=T, sep="\t", 
row.names=1, comment.char="", check.names=FALSE)
head(otu_tab)

taxonomy <- otu_tab$taxonomy
otu_tab <- otu_tab[-length(colnames(otu_tab))]  
otu_tab <- t(as.matrix(otu_tab))
head(otu_tab)

library(vegan)
otu_tab_rarefy <- rrarefy(otu_tab, min(apply(otu_tab,1,sum)))

library(ape)
otu_tree <- read.tree("fasttree_all_seed_OTUs.tre")
otu_tree 

library(phangorn)
otu_tree <- midpoint(otu_tree)
otu_tree

library(GUniFrac)
# Calculate the UniFracs
unifracs <- GUniFrac(otu_tab_rarefy, otu_tree, alpha=c(0, 0.5, 1))$unifracs

dw <- unifracs[, , "d_1"]		# Weighted UniFrac
du <- unifracs[, , "d_UW"]	# Unweighted UniFrac	
dv <- unifracs[, , "d_VAW"]	# Variance adjusted weighted UniFrac
d0 <- unifracs[, , "d_0"]   # GUniFrac with alpha 0  
d5 <- unifracs[, , "d_0.5"] # GUniFrac with alpha 0.5 

meta_tab<- read.table("metadata.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
head(meta_tab)

otu_meta_matched <- match(rownames(meta_tab),rownames(otu_tab))
otu_meta_matched <- otu_meta_matched[!is.na(otu_meta_matched)]
otu_tab <- otu_tab[otu_meta_matched,]
meta_tab <- meta_tab[match(rownames(otu_tab),rownames(meta_tab)),]

set.seed(123)
adonis(as.dist(d5) ~ meta_tab$Gestation)
adonis(as.dist(d5) ~ meta_tab$Gender)
adonis(as.dist(d5) ~ meta_tab$milktype)

# Combine d(0), d(0.5), d(1) for testing
PermanovaG(unifracs[, , c("d_0", "d_0.5", "d_1")] ~ meta_tab$Gestation)


