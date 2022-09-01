################################################################################ 
###Chapter 10: Compositional Analysis of Microbiome Data 
###Yinglin Xia: September, 2018
################################################################################ 
options(width=78,digits=4)

################################################################################ 
###10.3. Exploratory Compositional Data Analysis                                 
################################################################################ 

setwd("E:/Home/MicrobiomeStatR/Analysis") 

##10.3.1 Compositional Biplot 
abund_table=read.csv("VdrFecalGenusCounts.csv",row.names=1,check.names=FALSE)
abund_table_t<-t(abund_table)
head(abund_table_t)

install.packages("zCompositions")
library (zCompositions)
abund_table_r <- t(cmultRepl((abund_table_t), method="CZM", output="counts"))
head(abund_table_r)

abund_table_prop <- apply(abund_table_r, 2, function(x){x/sum(x)})
head(abund_table_prop)

abund_table_prop_f<- abund_table_r[apply(abund_table_prop, 1, min) > 0.001,]
head(abund_table_prop_f)

names_add <- rownames(abund_table_prop_f)[
  order(apply(abund_table_prop_f, 1, sum), decreasing=T) ]

abund_table_prop_reduced <- abund_table_prop_f[names_add,]
head(abund_table_prop_reduced)

abund_clr <- t(apply(abund_table_prop_reduced, 2, function(x){log(x) - mean(log(x))}))
head(abund_clr)

abund_PCX <- prcomp(abund_clr)
abund_PCX$x 

library(compositions)
#Sum the total variances
sum(abund_PCX$sdev[1:2]^2)/mvar(abund_clr)

samples <- c(rep(1, 5,rownames(abund_PCX$x)),
    rep(2, 3,rownames(abund_PCX$x))) 

palette=palette(c(rgb(1,0,0,0.6), rgb(0,0,1,0.6), rgb(.3,0,.3,0.6)))
palette

library(compositions)
coloredBiplot(abund_PCX, col="black", cex=c(0.6, 0.5),xlabs.col=samples, 
              arrow.len=0.05,
              xlab=paste("PC1 ", round (sum(abund_PCX$sdev[1]^2)/mvar(abund_clr),3), sep=""),
              ylab=paste("PC2 ", round (sum(abund_PCX$sdev[2]^2)/mvar(abund_clr),3), sep=""),
              expand=0.8,var.axes=T, scale=1, main="Biplot")

 
    
##10.3.2 Compositional Scree Plot

layout(matrix(c(1,2),1,2, byrow=T), widths=c(6,4), heights=c(6,4))
par(mgp=c(2,0.5,0))
screeplot(abund_PCX, type = "lines", main="Scree plot")
screeplot(abund_PCX, type = "barplot", main="Scree plot")        


##10.3.3 Compositional Cluster Dendrogram

# generate the distance matrix
dist <- dist(abund_clr, method="euclidian")

# cluster the data
hc <- hclust(dist, method="ward.D2")
hc
# plot the dendrogram
plot(hc, cex=1.0)

##10.3.4 Compositional Barplot
re_order <- abund_table_prop_reduced[,hc$order]
re_order

library(compositions)
re_order_acomp <- acomp(t(re_order))
par(mfrow=c(1,2))
#par(mar=c(3,1,1,1)+0.8)
colours <- rainbow(10)
# plot the barplot below
barplot(re_order_acomp, legend.text=F, col=colours, axisnames=F, border=NA, xpd=T)
# and the legend
plot(1,2, pch = 1, lty = 1, ylim=c(-10,10), type = "n", axes = FALSE, ann = FALSE)
legend(x="center", legend=names_add, col=colours, lwd=5, cex=.6, border=NULL)


##plot the dendrogram and barplot together
layout(matrix(c(1,3,2,3),2,2, byrow=T), widths=c(5,2), height=c(3,4))
par(mar=c(3,1,1,1)+0.8)

# plot the dendrogram
plot(hc, cex=0.6)

# plot the barplot below
barplot(re_order_acomp, legend.text=F, col=colours, axisnames=F, border=NA, xpd=T)

# and the legend
plot(1,2, pch = 1, lty = 1, ylim=c(-10,10), type = "n", axes = FALSE, ann = FALSE)
legend(x="center", legend=names_add, col=colours, lwd=5, cex=.6, border=NULL)


################################################################################
###10.4.Comparison between the Groups Using ALDEx2 Package
################################################################################

##10.4.2 Compositional Data Analysis Using ALDEx2 

abund_table=read.csv("VdrSitesGenusCounts.csv",row.names=1,check.names=FALSE)
abund_table_t<-t(abund_table)

ncol(abund_table_t)  # check the number of genera 
nrow(abund_table_t)  # check the number of samples

meta_table <- data.frame(row.names=rownames(abund_table_t),t(as.data.frame(strsplit(rownames(abund_table_t),"_"))))
meta_table

groups <- with(meta_table,ifelse(as.factor(X3)%in% c("drySt-28F"),c("VdrFecal"), c("VdrCecal")))
groups

abund_table[1:3,1:3]

source("https://bioconductor.org/biocLite.R")
biocLite("ALDEx2")

library(ALDEx2)
#This operation is fast.
vdr <- aldex.clr(abund_table, groups, mc.samples=128, verbose=TRUE)

vdr_t <- aldex.ttest(vdr, groups, paired.test=FALSE)
head(vdr_t)

vdr_glm <- aldex.glm(vdr, groups)
 
vdr_effect <- aldex.effect(vdr, groups, include.sample.summary=FALSE, verbose=FALSE)

vdr_all <- data.frame(vdr_t, vdr_glm, vdr_effect)
head(vdr_all)

sig_by_both <- which(vdr_all$we.ep < 0.05 & vdr_all$wi.ep < 0.05)
sig_by_both

sig_by_both_fdr <- which(vdr_all$we.eBH < 0.05 & vdr_all$wi.eBH < 0.05)
sig_by_both_fdr

library(xtable)
table <-xtable(
  vdr_all[sig_by_both,c(12:15,1,3,2,4)], caption="Table of significant taxa", digits=3,
  label="sig.table", align=c("l",rep("r",8) )
)
print.xtable(table, type="html", file="Vdr_Table.html")


vdr_w <- aldex(abund_table, groups, mc.samples=128, test="t", effect=TRUE,
                include.sample.summary=FALSE, denom="iqlr", verbose=FALSE)
head(vdr_w)



##10.4.3 Difference Plot, Effect Size and Effect Plot 

##Bland-Altman Plot
aldex.plot(vdr_all, type="MA", test="welch", cutoff=0.15, all.cex=0.7, called.cex=1.1,
           rare.col="grey", called.col="red")


##Effect Size and Effect Plot
par(mfrow=c(1,2))
aldex.plot(vdr_all, type="MW", test="welch",cutoff=0.15, all.cex=0.7, called.cex=1.1,
           rare.col="black", called.col="red")
aldex.plot(vdr_all, type="MW",test="wilcox", cutoff=0.15, all.cex=0.7, called.cex=1.1,
           rare.col="black", called.col="red")

par(mfrow=c(1,2))
plot(vdr_all$effect, vdr_all$wi.ep, log="y", pch=19, main="Effect",
     cex=0.5, xlab="Effect size", ylab="Expected P value of Wilcoxon rank test")
abline(h=0.05, lty=2,lwd=2, col ='red')
plot(vdr_all$diff.btw, vdr_all$wi.ep, log="y", pch=19, main="Volcano",
     cex=0.5, xlab="Difference", ylab="Expected P value of Wilcoxon rank test")
abline(h=0.05, lty=2,lwd=2, col='red')


################################################################################
###10.5. Proportionality: Correlation Analysis for Relative Data
################################################################################

##10.5.3 Illustrating Proportionality Analysis

##10.5.3.1 Calculating Proportionality
abund_table=read.csv("VdrFecalGenusCounts.csv",row.names=1,check.names=FALSE)
head(abund_table)  
abund_table_t<-t(abund_table)

library(propr)
phi <- phit(abund_table_t, symmetrize = TRUE)
rho <- perb(abund_table_t, ivar = 0)
phs <- phis(abund_table_t, ivar = 0)
phi

head(phi@counts)
head(phi@logratio)
head(phi@pairs)

head(phi@matrix)
head(rho@matrix)
head(phs@matrix)

##10.5.3.2 Identify Proportionally Abundant Taxa 
keep <- apply(abund_table_t, 2, function(x) sum(x >= 10) >= 5)
rho <- perb(abund_table_t, select = keep)

best <- rho[">", 0.80]
best
best@pairs

pirs_taxa<-row.names(abund_table[c(18,92,108,123),])
pirs_taxa

taxa_best <- colnames(best@logratio)
taxa_best 
head(best@matrix)

##10.5.3.3 Visualizing Proportionality

##Index-Aware Plots
install.packages("rlang")
plot(best)

install.packages('ggdendro')
dendrogram(best)

##Index-Naive Plots----#
best <- simplify(best)

grouping<-data.frame(row.names=rownames(abund_table_t),t(as.data.frame(strsplit(rownames(abund_table_t),"_"))))
grouping$Group <- with(grouping,ifelse(as.factor(X2)%in% c(11,12,13,14,15),c("Vdr-/-"), c("WT")))

pca(best, group = grouping$Group)

clusts <- bucket(best, group = grouping$Group, k = 2)

clusts <- prism(best, k = 2)
clusts <- bokeh(best, k = 2)

snapshot(best)

##Down-Stream Plots
sub <- subset(best, select = (clusts == 2))
pca(sub, group = grouping$Group)

taxa <- colnames(sub@logratio)
taxa 

##10.5.3.4 Conducting Differential Proportionality Analysis

library(propr)
pd <- propd(abund_table_t, grouping$Group, alpha = NA, p = 1000)

theta_d <- setDisjointed(pd)
theta_d
theta_e <- setEmergent(pd)
theta_e

theta_d <- updateCutoffs(theta_d, cutoff = seq(0.05, 0.95,0.3))
theta_d

theta_e <- updateCutoffs(theta_e, cutoff = seq(0.05, 0.95,0.3))
theta_e

pd_nn <- propd(abund_table_t, grouping$Group, weighted = FALSE)
pd_wn <- propd(abund_table_t, grouping$Group, weighted = TRUE)
pd_na <- propd(abund_table_t, grouping$Group, weighted = FALSE,alpha = 0.05)
pd_wa <- propd(abund_table_t, grouping$Group, weighted = TRUE, alpha = 0.05)

pd_nn<-updateF(pd_nn,moderated = FALSE)
options(digits=4)
head(pd_nn@theta$Fstat)

pd_wa<-updateF(pd_wa,moderated = FALSE)
options(digits=4)
head(pd_wa@theta$Fstat)

pd_nn<-updateF(pd_nn,moderated = TRUE, ivar = "clr")
pd_wn<-updateF(pd_wn,moderated = TRUE, ivar = "clr")
pd_na<-updateF(pd_na,moderated = TRUE, ivar = "clr")
pd_wa<-updateF(pd_wa,moderated = TRUE, ivar = "clr")


head(pd_nn@theta$Fstat)
head(pd_nn@theta$theta_mod)
head(pd_wa@theta$Fstat)
head(pd_wa@theta$theta_mod)

