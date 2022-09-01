################################################################################ 
###Chapter 11: Modeling Over-dispersed Microbiome Data 
###Yinglin Xia: September, 2018                                                                   
################################################################################ 

################################################################################ 
###11.2. NB Model in edgeR                                                      
################################################################################ 

##11.3.2 Step-by-Step Implementing edgeR

setwd("F:/Home/MicrobiomeStatR/Analysis")

source("http://www.Bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
biocLite("GUniFrac")
#install.packages("matrixStats")

library(GUniFrac)
data(throat.otu.tab)
head(throat.otu.tab)

throat<-t(throat.otu.tab)
counts<-throat
head(counts)

data(throat.meta)
group <- throat.meta$SmokingStatus
head(group)

dim(counts)
length(group)

library(edgeR)
y <- DGEList(counts=counts,group=group) 
names(y)

head(y$counts)    # original count matrix
y$samples         # contains a summary of samples
sum(y$all.zeros ) # How many genes have 0 counts across all samples

dim(y) 
y_full <- y # keep the old one in case we mess up
head(y$counts)
apply(y$counts, 2, sum) # total OTU(gene) counts per sample
keep <- rowSums(cpm(y)>100) >= 2  
y <- y[keep,]
dim(y)

y$samples$lib.size <- colSums(y$counts)
y$samples

y <- calcNormFactors(y)
y  
y$samples

# effective library sizes
y$samples$lib.size*y$samples$norm.factors

plotMDS(y, method="bcv", main = "MDS Plot for throat Count Data", 
        col=as.numeric(y$samples$group), cex=0.5, labels = colnames(y$counts)) 
legend("topright", as.character(unique(y$samples$group)), col=1:2, cex=0.8, pch=16)

# Output plot as a pdf
pdf("MDS_plot.pdf", width = 7 , height = 7 ) # in inches
plotMDS(y, method="bcv", main = "MDS Plot for throat Count Data", 
        col=as.numeric(y$samples$group), cex=0.5,labels = colnames(y$counts))
legend("topright", as.character(unique(y$samples$group)), col=1:2, cex=0.8, pch=16)
dev.off() # tells R to turn off device and writing to the pdf.  

#estimate the common dispersion
y1 <- estimateCommonDisp(y, verbose=T)
names(y1)

#estimate the tag-wise dispersion
y1 <- estimateTagwiseDisp(y1)
names(y1)

plotBCV(y1)

#use a generalized linear model to estimate the dispersion
design <- model.matrix(~group)
rownames(design) <- colnames(y)
design

install.packages("statmod")
library(statmod)
y2 <- estimateDisp(y, design, robust=TRUE)
y2$common.dispersion

plotBCV(y2)

fit <- glmQLFit(y2, design, robust=TRUE)
plotQLDisp(fit)

##The exactTest() Approach
et <- exactTest(y1,pair = c( "NonSmoker", "Smoker" ))
topTags(et)

et1 <- exactTest(y1, pair=c(1,2))
topTags(et1)

y3 <- y
y3$samples$group <- relevel(y3$samples$group, ref="Smoker")
levels(y3$samples$group)

et <- exactTest(y1)
et <- exactTest(y1,pair = c( "NonSmoker", "Smoker" ))

  
##GLM Appraoch
design <- model.matrix(~group)
rownames(design) <- colnames(y)
design

fit <- glmQLFit(y1, design)
qlf <- glmQLFTest(fit, contrast=c(-1,1))
topTags(qlf)

FDR <- p.adjust(qlf$table$PValue, method="BH")
sum(FDR < 0.05) 
topTags(qlf,n=15)

qlf_lrt <- glmLRT(fit, contrast=c(-1,1))
topTags(qlf_lrt)


#another design
design1 <- model.matrix(~0+group, data=y$samples)
colnames(design1) <- levels(y1$samples$group)
design1

fit1 <- glmQLFit(y1, design1)
qlf1 <- glmQLFTest(fit1, contrast=c(-1,1))
topTags(qlf1)

FDR1 <- p.adjust(qlf1$table$PValue, method="BH")
sum(FDR1 < 0.05) 
topTags(qlf1,n=15)

qlf1_lrt <- glmLRT(fit1, contrast=c(-1,1))
topTags(qlf1_lrt)

##MA-Plot Using plotSmear()                                                 
da = decideTestsDGE(et1 , p.value = 0.1)
da_OTUs = rownames(y1)[as.logical(da)]
plotSmear(et1, de.tags = da_OTUs, cex = 0.5)
abline(h = c(-2, 2), col = "blue")


##Volcano Plot
tab = data.frame(logFC = et1$table[, 1], negLogPval = -log10(et1$table[, 3]))
head(tab)

par(mar = c(5, 4, 4, 4))
plot(tab, pch = 16, cex = 0.6, xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue))

# Log2 fold change and p-value cutoffs
lfc = 2
pval = 0.1

# Selecting interest OTUs
sig_OTUs = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))

# Identifying the selected OTUs
points(tab[sig_OTUs, ], pch = 16, cex = 0.8, col = "red")
abline(h = -log10(pval), col = "green3", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),
      cex = 0.8, line = 0.5)


################################################################################ 
###11.5. The DESeq and DESeq2 Packages                                              
################################################################################ 

##11.5.2 Step-by-Step Implementing DESeq2
library(GUniFrac)
data(throat.otu.tab)
head(throat.otu.tab) 
otu_tab<-throat.otu.tab
head(otu_tab)
  
countData<-as(otu_tab, "matrix")
head(countData)

countData<-(t(countData))#DESeq2 need taxa(genes=rows) by samples(=columns)format
head(countData)

data(throat.meta)
head(throat.meta)

group<-throat.meta$SmokingStatus
head(group)

metaData<-data.frame(row.names=colnames(countData),group=group)
head(metaData)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = metaData,
                              design = ~ group)

dds <- dds[rowSums(counts(dds)) > 0,]
dds

dds <- estimateSizeFactors(dds) 
sizeFactors(dds)

dds<- estimateDispersions(dds)

dds$group <- relevel(dds$group, "NonSmoker")
#or by
dds$group <- factor(dds$group, levels = c("NonSmoker", "Smoker"))

dds <- DESeq(dds)

res <- results(dds)
res

mcols(res, use.names=TRUE)
mcols(dds,use.names=TRUE)[1:4,1:4]

substr(names(mcols(dds)),1,10)
head(assays(dds)[["mu"]])
head(dispersions(dds))
head(mcols(dds)$dispersion)

sizeFactors(dds)
head(coef(dds))

res <- results(dds, contrast = c("group", "Smoker", "NonSmoker") )
res

sum(res$pvalue < 0.01, na.rm=TRUE )
table(is.na(res$pvalue))

table(res[,"padj"] < 0.1)  
sum(res$padj < 0.1, na.rm=TRUE )

res_Sig <- res[which(res$padj < 0.1 ),]
head(res_Sig[order(res_Sig$log2FoldChange),])
tail(res_Sig[order( res_Sig$log2FoldChange ),])

##Diagnostic Plots Using plotMA
plotMA(res)

##Diagnostic Plots Using plotDispEsts
plotDispEsts(dds, ylim = c(1e-2, 1e3))
        
##Clustering with Heatmap
rld <- rlog(dds)
vst <-varianceStabilizingTransformation(dds)

par(mfrow = c(1, 3))
plot(log2( 1+counts(dds, normalized=TRUE)[,1:2] ), main="Ordinary log2",col="#00000020", pch=20, cex=0.3 )
plot(assay(rld)[,1:2], main="Regularized-logarithm", col="#00000020", pch=20, cex=0.3 )
plot(assay(vst)[,1:2], main="Variance stabilizing",col="#00000020", pch=20, cex=0.3 )

head(assay(rld))[,1-3]

install.packages("gplots")
library("gplots" )
library("RColorBrewer" )
library("genefilter" )
library(SummarizedExperiment)
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 10 )
heatmap.2(assay(rld)[ topVarGenes, ], scale="row",
           trace="none", dendrogram="column",
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))


##Histogram of p-Values
hist(res$pvalue, breaks=20, col="grey",
     main = "Smoker vs. NonSmoker", xlab = "p-values")

  
##Independent Filtering
metadata(res)

metadata(res)$alpha
metadata(res)$filterThreshold

plot(metadata(res)$filterNumRej,
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)


##Re-estimate the p-Values
res   

# remove filtered out OTUs by independent filtering 
# they have NA adj. pvals
res <- res[ !is.na(res$padj),]
# with NA pvals (outliers)
res <- res[ !is.na(res$pvalue),]

res <- res[, -which(names(res) == "padj")]

install.packages("fdrtool")
library(fdrtool)
res_fdr <- fdrtool(res$stat, statistic= "normal", plot = T)

head(res_fdr)

res_fdr$param[1, "sd"]
sd

res[,"padj"] <- p.adjust(res_fdr$pval, method = "BH")
 
hist(res_fdr$pval, col = "gray",
     main = "Smoker vs. NonSkoer, correct null model", xlab = "Corrected p-values")

table(res[,"padj"] < 0.1)

res[1:2,]
write.csv( as.data.frame(res), file="results.csv" )


