################################################################################
###Chapter 7: Exploratory Analysis and Beyond 
###Yinglin Xia: September, 2018  
################################################################################
options(width=78,digits=4)

##Set working directory                                                                  
setwd("E:/Home/MicrobiomeStatR/Analysis")

################################################################################ 
###7.2 Exploratory Analysis with Graphic Summary                                 
################################################################################ 

##7.2.1 Plot Richness
library(phyloseq)
library(ggplot2)
abund_table=read.csv("VdrFecalGenusCounts.csv",row.names=1,check.names=FALSE)
abund_table<-t(abund_table)

meta_table <- data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
meta_table$Group <- with(meta_table,ifelse(as.factor(X2)%in% c(11,12,13,14,15),c("Vdr-/-"), c("WT")))

library(phyloseq)
#Convert the data to phyloseq format
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
SAM = sample_data(meta_table)
physeq<-merge_phyloseq(phyloseq(OTU),SAM)
physeq

library(ggplot2)
plot_richness(physeq, x = "Group ", color = "Group ")

plot_richness(physeq, measures = c("Chao1", "Shannon"),x = "Group ", color = "Group ")


##7.2.2 Plot Abundance Bar
theme_set(theme_bw())

plot_bar(physeq)

plot_bar(physeq, fill="Group")

TopNGenus <- names(sort(taxa_sums(physeq), TRUE)[1:5])
Top5Genus  <- prune_taxa(TopNGenus, physeq)

plot_bar(Top5Genus, fill="Group", facet_grid=~Group)
plot_bar(Top5Genus, fill="Genus", facet_grid=~Group)


##7.2.3 Plot Heat Map                                                                      
TopNGenus <- names(sort(taxa_sums(physeq), TRUE)[1:5])
Top5Genus  <- prune_taxa(TopNGenus, physeq)

plot_heatmap(Top5Genus) 
(p <- plot_heatmap(Top5Genus, "NMDS", "bray"))

plot_heatmap(Top5Genus, "NMDS", "bray", low="#000033", high="#CCFF66")  
plot_heatmap(Top5Genus, "NMDS", "bray", low="#000033", high="#FF3300")
plot_heatmap(Top5Genus, "NMDS", "bray", low="#000033", high="#66CCFF")
plot_heatmap(Top5Genus, "NMDS", "bray", low="#66CCFF", high="#000033", na.value="white")

plot_heatmap(Top5Genus, "PCoA", "bray")

##7.2.4	Plot Network
set.seed(123)
library(igraph)
ig <- make_network(physeq, max.dist=0.8)
plot_network(ig, physeq, color="Group", shape="Group")

plot_net(physeq, maxdist = 0.5, color = "Group", shape="Group")


##7.2.5	Plot Phylogenetic Tree                                                             
library(GUniFrac)
data(throat.otu.tab)
data(throat.tree)
data(throat.meta)

library(phyloseq); packageVersion("phyloseq") 
packageVersion("ggplot2")
#Convert the data to phyloseq format
OTU = otu_table(as.matrix(throat.otu.tab), taxa_are_rows = FALSE)
SAM = sample_data(throat.meta)
TRE <-throat.tree

physeq<-merge_phyloseq(phyloseq(OTU),SAM,TRE)

ntaxa(physeq)   

physeq = prune_taxa(taxa_names(physeq)[1:50], physeq) 

plot_tree(physeq, ladderize = "left", color = "SmokingStatus") 
plot_tree(physeq, ladderize = "left", color = "SmokingStatus",shape = "SmokingStatus") 

plot_tree(physeq, color = "SmokingStatus", 
          shape = "SmokingStatus", ladderize = "left") + coord_polar(theta = "y")

 
################################################################################ 
###7.3	Clustering                                                                   
################################################################################ 
##7.3.2 ClusteringLoad the package and datasets
library(vegan)  

abund_table_norm <- decostand(abund_table, "normalize")
bc_dist<- vegdist(abund_table_norm , method = "bray")

##7.3.2.1 Single Linkage Agglomerative Clustering
cluster_single <- hclust (bc_dist, method = 'single')
plot(cluster_single)

##7.3.2.2 Complete Linkage Agglomerative Clustering
cluster_complete <- hclust (bc_dist, method = 'complete')
plot(cluster_complete)

##7.3.2.3 Average Linkage Agglomerative Clustering
cluster_average <- hclust (bc_dist, method = 'average')
plot(cluster_average)

##7.3.2.4	Ward's Minimum Variance Clustering
cluster_ward <- hclust (bc_dist, method = 'ward.D2')
plot(cluster_ward)

par (mfrow = c(2,2))
plot(cluster_single)
plot(cluster_complete)
plot(cluster_average)
plot(cluster_ward)

par (mfrow = c(1,1))


################################################################################
###7.4	Ordination                          
################################################################################

##7.4.1	Principal Component Analysis (PCA)
abund_table=read.csv("VdrFecalGenusCounts.csv",row.names=1,check.names=FALSE)
abund_table<-t(abund_table)
head(abund_table)

meta_table <- data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
meta_table$Group <- with(meta_table,ifelse(as.factor(X2)%in% c(11,12,13,14,15),c("Vdr-/-"), c("WT")))
 
stand_abund_table <- decostand(abund_table, method = "total")
PCA <-rda(stand_abund_table)
PCA

sum (apply (stand_abund_table, 2, var))

biplot(PCA, display = 'species')
ordiplot(PCA, display = "sites", type = "text")


"cleanplot.pca" <- function(res.pca, ax1=1, ax2=2, point=FALSE, 
                            ahead=0.07, cex=0.5) 
{
  # A function to draw two biplots (scaling 1 and scaling 2) from an object 
  # of class "rda" (PCA or RDA result from vegan's rda() function)
  #
  # License: GPL-2 
  # Authors: Francois Gillet & Daniel Borcard, 24 August 2012
  
  require("vegan")
  
  par(mfrow=c(1,2))
  p <- length(res.pca$CA$eig)
  
  # Scaling 1: "species" scores scaled to relative eigenvalues
  sit.sc1 <- scores(res.pca, display="wa", scaling=1, choices=c(1:p))
  spe.sc1 <- scores(res.pca, display="sp", scaling=1, choices=c(1:p))
  plot(res.pca, choices=c(ax1, ax2), display=c("wa", "sp"), type="n", 
       main="PCA - scaling 1", scaling=1)
  if (point)
  {
    points(sit.sc1[,ax1], sit.sc1[,ax2], pch=20)
    text(res.pca, display="wa", choices=c(ax1, ax2), cex=cex, pos=1, scaling=1)
  }
  else
  {
    text(res.pca, display="wa", choices=c(ax1, ax2), cex=cex, scaling=1)
  }
  text(res.pca, display="sp", choices=c(ax1, ax2), cex=cex, pos=1, col="blue", scaling=1)
  arrows(0, 0, spe.sc1[,ax1], spe.sc1[,ax2], length=ahead, angle=20, col="red")
  pcacircle(res.pca)
  
  # Scaling 2: site scores scaled to relative eigenvalues
  sit.sc2 <- scores(res.pca, display="wa", choices=c(1:p))
  spe.sc2 <- scores(res.pca, display="sp", choices=c(1:p))
  plot(res.pca, choices=c(ax1,ax2), display=c("wa","sp"), type="n", 
       main="PCA - scaling 2")
  if (point) {
    points(sit.sc2[,ax1], sit.sc2[,ax2], pch=20)
    text(res.pca, display="wa", choices=c(ax1 ,ax2), cex=cex, pos=1)
  }
  else
  {
    text(res.pca, display="wa", choices=c(ax1, ax2), cex=cex)
  }
  text(res.pca, display="sp", choices=c(ax1, ax2), cex=cex, pos=1, col="blue")
  arrows(0, 0, spe.sc2[,ax1], spe.sc2[,ax2], length=ahead, angle=20, col="red")
}


"pcacircle" <- function (pca) 
{
  # Draws a circle of equilibrium contribution on a PCA plot 
  # generated from a vegan analysis.
  # vegan uses special constants for its outputs, hence 
  # the 'const' value below.
  
  eigenv <- pca$CA$eig
  p <- length(eigenv)
  n <- nrow(pca$CA$u)
  tot <- sum(eigenv)
  const <- ((n - 1) * tot)^0.25
  radius <- (2/p)^0.5
  radius <- radius * const
  symbols(0, 0, circles=radius, inches=FALSE, add=TRUE, fg=2)
}

cleanplot.pca(PCA)


##7.4.2	Principal Coordinate Analysis (PCoA) 
 
options(digits=4)
options(width=78,digits=4)
library(vegan)

bc_dist <-vegdist(abund_table, "bray")

PCoA <- cmdscale (bc_dist, eig = TRUE,k = 2)
PCoA

explainedvar1 <- round(PCoA$eig[1] / sum(PCoA$eig), 2) * 100
explainedvar1

explainedvar2 <- round(PCoA$eig[2] / sum(PCoA$eig), 2) * 100
explainedvar2

sum_eig <- sum(explainedvar1, explainedvar2)
sum_eig

# Define Plot Parameters
par(mar = c(5, 5, 1, 2) + 0.1)

# Plot Eigenvalues
plot(PCoA$eig, xlab = "PCoA", ylab = "Eigenvalue",
     las = 1, cex.lab = 1.5, pch = 16)

# Add Expectation based on Kaiser-Guttman criterion and Broken Stick Model
abline(h = mean(PCoA$eig), lty = 2, lwd = 2, col = "blue")
b_stick <- bstick(8, sum(PCoA$eig))
lines(1:8, b_stick, type = "l", lty = 4, lwd = 2, col = "red")
# Add Legend
legend("topright", legend = c("Avg Eigenvalue", "Broken-Stick"),
       lty = c(2, 4), bty = "n", col = c("blue", "red"))

# Define Plot Parameters
par(mar = c(5, 5, 1, 2) + 0.1)

# Initiate Plot
plot(PCoA$points[ ,1], PCoA$points[ ,2], ylim = c(-0.5, 0.5),
     xlab = paste("PCoA 1 (", explainedvar1, "%)", sep = ""),
     ylab = paste("PCoA 2 (", explainedvar2, "%)", sep = ""),
     pch = 5, cex = 1.0, type = "n", cex.lab = 1.0, cex.axis = 1.2, axes = FALSE)

# Add Axes
axis(side = 1, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
axis(side = 2, labels = T, lwd.ticks = 2, cex.axis = 1.2, las = 1)
abline(h = 0, v = 0, lty = 3)
box(lwd = 2)

# Add Points & Labels
points(PCoA$points[ ,1], PCoA$points[ ,2],
       pch = 19, cex = 3, bg = "blue", col = "blue")
text(PCoA$points[ ,1], PCoA$points[ ,2],
     labels = row.names(PCoA$points))


# Calculating Relative Abundance
fecalREL <- abund_table
for(i in 1:nrow(abund_table)){
  fecalREL[i, ] = abund_table[i, ] / sum(abund_table[i, ])
}

install.packages("BiodiversityR")
library("pbkrtest")
library("BiodiversityR")

# Calculate and Add Species Scores
PCoA <- add.spec.scores(PCoA,fecalREL,method = "pcoa.scores",Rscale=TRUE,scaling=1, multi=1)
text(PCoA$cproj[ ,1], PCoA $cproj[ ,2],
     labels = row.names(PCoA$cproj),cex=0.5, col = "blue")

genus_corr <- add.spec.scores(PCoA, fecalREL, method = "cor.scores")$cproj
corrcut <- 0.7 # user defined cutoff  
import_genus <- genus_corr[abs(genus_corr[, 1]) >= corrcut | abs(genus_corr[, 2]) >= corrcut, ]
import_genus[complete.cases(import_genus),]


# Permutation Test for Species Abundances Across Axes
envfit(PCoA, fecalREL, perm = 999)

fit <- envfit(PCoA, fecalREL, perm = 999)
plot(fit, p.max = 0.05, cex=0.5, col = "red")


##7.4.3 Non-metric multidimensional scaling (NMDS) 

bc_nmds <- metaMDS(abund_table, dist = "bray")
bc_nmds

ordiplot (bc_nmds, type = 't')

ordiplot(bc_nmds, display = "sites", type = "text")


par (mfrow = c(1,2))
stressplot (bc_nmds)
plot (bc_nmds, display = 'sites', type = 't', main = 'Goodness of fit')
points (bc_nmds, display = 'sites', cex = goodness (bc_nmds)*300)

                  
##7.4.4 Correspondence Analysis (CA)

fecal_genus_cca=cca(abund_table)
fecal_genus_cca

plot(fecal_genus_cca, display="sites") 

plot(fecal_genus_cca, display="sites", type="p") 

ordiplot (fecal_genus_cca)

evplot <- function(ev)
{# Broken stick model (MacArthur 1957)
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p <- 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op <- par(mfrow=c(2,1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}
      
# Plot eigenvalues and % of variance for each axis
ev <- fecal_genus_cca$CA$eig

windows(title="CA eigenvalues")
evplot(ev)
  

##7.4.5 Redundancy Analysis (RDA)

library(GUniFrac)
data(throat.otu.tab)
data(throat.meta)

library(dplyr)
throat_meta <- select(throat.meta, SmokingStatus, Age, Sex, PackYears)

library (vegan)  # if you haven't use it up to now
abund_hell <- decostand (throat.otu.tab, 'hell') 

rda_hell<- rda(abund_hell ~ ., throat_meta)
summary (rda_hell)

coef(rda_hell)

RsquareAdj(rda_hell)

rda_hell$CA$eig[rda_hell$CA$eig > mean(rda_hell$CA$eig)]

plot(rda_hell, display=c("sp", "lc", "cn"),main="Triplot RDA - scaling 2")
taxa_scores <- scores(rda_hell, choices=c(1,2), display="sp")
arrows(0, 0, taxa_scores[,1], taxa_scores[,2], length=0, lty=1, col="red")

set.seed (123)

anova(rda_hell, step=1000)

anova(rda_hell, by="axis", step=1000)

anova(rda_hell, by="terms", step=1000)

step_forward <- ordistep(rda(abund_hell ~ 1, data=throat_meta), 
                         scope=formula(rda_hell), direction="forward", pstep=1000)

step_forward <- ordiR2step(rda(abund_hell ~ 1, data=throat_meta), 
                           scope=formula(rda_hell), direction="forward", pstep=1000)

rda_final<- rda(abund_hell ~ SmokingStatus + Sex, data=throat_meta)
anova(rda_final, step=1000)

anova(rda_final, by="axis", step=1000)
anova(rda_final, by="terms", step=1000)


##7.4.6 Constrained Correspondence Analysis (CCA)

smoker_cca <- cca(throat.otu.tab ~ ., throat_meta)
smoker_cca

summary(smoker_cca)

plot(smoker_cca, scaling=1, display = c("lc","cn"), main="Biplot CCA - scaling 1") 

set.seed (123)

anova(smoker_cca, step=1000)  
anova (smoker_cca, by = 'axis',step=1000)
anova (smoker_cca, by = 'terms')

ordistep(cca(throat.otu.tab ~ 1, data=throat_meta), scope=formula(smoker_cca), 
         direction="forward", pstep=1000)

(smoker_cca_final <- cca(throat.otu.tab ~ SmokingStatus, data=throat_meta))
anova.cca(smoker_cca_final, step=1000)
anova.cca(smoker_cca_final, step=1000, by="terms")


##7.4.7 Constrained Analysis of Principal Coordinates (CAP)
library(dplyr)

throat_meta <- select(throat.meta, SmokingStatus, Age, Sex, PackYears)
throat_cap <- capscale(throat.otu.tab ~ SmokingStatus + Sex + PackYears + Condition(Age), throat_meta,
                       dist="bray")
throat_cap
anova(throat_cap)

plot(throat_cap)

groups <- throat.meta$SmokingStatus
groups

plot(throat_cap, type="n") 
points(throat_cap, col=as.numeric(as.factor(groups)), 
       pch=as.numeric(as.factor(groups))) 
ordispider(throat_cap, groups, lty=2, col="grey", label=T) 
ordiellipse(throat_cap, groups, lty=2, col="grey", label=F) 

set.seed (123)

anova(throat_cap, step=1000)
anova(throat_cap, by="axis", step=1000)
anova(throat_cap, by="terms", step=1000)

step_forward <- ordistep(capscale(throat.otu.tab ~ 1, data=throat_meta), 
                         scope=formula(throat_cap), direction="forward", pstep=1000)

cap_final<- capscale(throat.otu.tab ~ SmokingStatus, throat_meta, dist="bray")
anova(cap_final, step=1000)

