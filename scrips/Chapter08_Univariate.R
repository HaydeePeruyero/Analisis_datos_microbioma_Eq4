################################################################################
###Chapter 8: Univariate Community Analysis                                      
###Yinglin Xia: September, 2018                                                 
################################################################################

setwd("C:/Users/ANDRESARREDONDOCRUZ/Equipo4/Chapter8") 

##8.1 Comparisons of Diversities between Two Groups 

##8.1.1 Two-sample Welch's t-test
abund_table=read.csv("C:/Users/ANDRESARREDONDOCRUZ/Equipo4/data/VdrGenusCounts.csv",row.names=1,check.names=FALSE)
abund_table<-t(abund_table)

grouping<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
grouping$Location <- with(grouping, ifelse(X3%in%"drySt-28F", "Fecal", "Cecal"))
grouping$Group <- with(grouping,ifelse(as.factor(X2)%in% c(11,12,13,14,15),c("Vdr-/-"), c("WT")))
grouping <- grouping[,c(4,5)]
grouping 

library(vegan)

H<-diversity(abund_table, "shannon") 

df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon",length(H)))
df_G <-cbind(df_H, grouping)
rownames(df_G)<-NULL
df_G

Fecal_G<- subset(df_G, Location=="Fecal")
Fecal_G

library(ggplot2)

p<-ggplot(Fecal_G, aes(x=value))+
  geom_histogram(color="black", fill="black")+
  facet_grid(Group ~ .)

library(plyr)
mu <- ddply(Fecal_G, "Group", summarise, grp.mean=mean(value))
head(mu)

p+geom_vline(data=mu, aes(xintercept=grp.mean, color="red"),
             linetype="dashed")

fit_t <- t.test(value ~ Group, data=Fecal_G)
fit_t

##8.1.2 Wilcoxon Rank Sum Test
fit_w <- wilcox.test(value ~ Group, data=Fecal_G)
fit_w 

################################################################################
###8.2 Comparisons of a Taxon of Interest between Two Groups                                                                     
################################################################################ 

##8.2.1 Comparison of Relative Abundance using Wilcoxon Rank Sum Test
apply(abund_table,1, sum)

relative_abund_table <- decostand(abund_table, method = "total")
apply(relative_abund_table, 1, sum)

relative_abund_table[1:16,1:8]

(Bacteroides <-relative_abund_table[,8]) 

Bacteroides_G <-cbind(Bacteroides, grouping)
rownames(Bacteroides_G)<-NULL
Fecal_Bacteroides_G <- subset(Bacteroides_G, Location=="Fecal")
Fecal_Bacteroides_G

boxplot(Bacteroides ~ Group,data=Fecal_Bacteroides_G, col=rainbow(2),main="Bacteroides in Vdr WT/KO mice")

ggplot(Fecal_Bacteroides_G, aes(x=Group, y=Bacteroides,col=factor(Group))) + 
  geom_boxplot(notch=FALSE)
ggplot(Fecal_Bacteroides_G, aes(x=Group, y=Bacteroides)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4) 

#+ layer(stat_params = list(binwidth = 2))
#Esta opcion da un error

fit_w_b <- wilcox.test(Bacteroides ~ Group,data=Fecal_Bacteroides_G)
fit_w_b


##8.2.2 Comparison of Present or Absent Taxon using Chi-square Test
abund_table[1:16,1:27]
(Parabacteroides <- abund_table[,27])

Parabacteroides_G <-cbind(Parabacteroides, grouping)
rownames(Parabacteroides_G)<-NULL

Cecal_Parabacteroides_G <- subset(Parabacteroides_G, Location=="Cecal")
Cecal_Parabacteroides_G

Cecal_Parabacteroides_G$Present <- ifelse((Cecal_Parabacteroides_G$Parabacteroides > 0), "Present","Absent")
Cecal_Parabacteroides_G

library(MASS) 
tbl = table(Cecal_Parabacteroides_G$Group, Cecal_Parabacteroides_G$Present) 
tbl                 
chisq.test(tbl) 

fisher.test(tbl)

################################################################################ 
###8.3 Comparisons among More Than Two Groups Using ANOVA
################################################################################

##8.3.1 One-way ANOVA 
CH=estimateR(abund_table)[2,] 
df_CH <-data.frame(sample=names(CH),value=CH,measure=rep("Chao1",length(CH))) 
df_CH_G <-cbind(df_CH, grouping)
rownames(df_G)<-NULL
df_CH_G

df_CH_G$Group4<- with(df_CH_G, interaction(Location,Group))
df_CH_G

boxplot(value~Group4, data=df_CH_G, col=rainbow(4), main="Chao1 index")

library(ggplot2)
p <- ggplot(df_CH_G, aes(x=Group4, y=value),col=rainbow(4), main="Chao1 index") + 
  geom_boxplot()
p + coord_flip()
ggplot(df_CH_G, aes(x=Group4, y=value,col=factor(Group4))) + 
  geom_boxplot(notch=FALSE)

library(dplyr)
#Hay un conflicto con el paquete MASS con dplyr por lo que si corremos este codigo da un error
#df_CH_G4 <- select(df_CH_G, Group4, value)
#df_CH_G4
#Se soluciona de esta manera
df_CH_G4 <- df_CH_G %>%
  dplyr::select( Group4, value)
df_CH_G4
#hay un error tambien aqui
bartlett.test(df_CH_G4, Group4) 
qchisq(0.95, 1)

fligner.test(df_CH_G4, Group4)

fit = lm(formula = value~Group4,data=df_CH_G)
anova (fit)

summary(aov(value~Group4, data=df_CH_G))

aov_fit <- aov(value~Group4,data=df_CH_G) 
summary(aov_fit, intercept=T) 

qf(0.95, 12, 3)

install.packages("mnormt")
#error?
library(broom)

tidy(aov_fit)
augment(aov_fit)
#se desordenan ls columnas
glance(aov_fit)

   
##8.3.2 Pairwise and Tukey Multiple Comparisons

#Pairwise tests of mean differences
pairwise.t.test(df_CH_G$value, df_CH_G$Group4, p.adjust="none", pool.sd = T) 
 
#conservative Bonferroni adjustment
pairwise.t.test(df_CH_G$value, df_CH_G$Group4, p.adjust="bonferroni", pool.sd = T)

#Holm method 
pairwise.t.test(df_CH_G$value, df_CH_G$Group4, p.adjust="holm", pool.sd = T)

#Benjamini & Hochberg(BH)
pairwise.t.test(df_CH_G$value, df_CH_G$Group4, p.adjust="BH", pool.sd = T)

#Benjamini & Yekutieli
pairwise.t.test(df_CH_G$value, df_CH_G$Group4, p.adjust="BY", pool.sd = T)

#Tukey multiple comparisons of means
TukeyHSD(aov_fit, conf.level=.95)  

plot(TukeyHSD(aov(df_CH_G$value~df_CH_G$Group4), conf.level=.95))

################################################################################ 
###8.4 Comparisons among More than Two Groups Using Kruskal-Wallis Test                           
################################################################################ 

##8.4.2 Compare Diversities among Groups

library(dplyr)
Data <- mutate(df_CH_G, Group = factor(df_CH_G$Group4, levels=unique(df_CH_G$Group4)))

library(FSA)
Summarize(value ~ Group4, data = df_CH_G)

#Individual plots in panel of 2 columns and 2 rows
library(lattice)
histogram(~ value|Group4, data=df_CH_G,layout=c(2,2)) 

#kruskal wallis test of Chao 1 richness
kruskal.test(value ~ Group4, data = df_CH_G) 

qchisq(0.950, 3)

library(DescTools)
#Tukey method for adjusting p-values
Test_N = NemenyiTest(x = df_CH_G$value,
                 g = df_CH_G$Group4,
                 dist="tukey")
Test_N

library(FSA)
## "bh" suggests Benjamini and Hochberg  method for adjusting p-values
Test_N = dunnTest(df_CH_G$value ~ df_CH_G$Group4,data=df_CH_G, method="bh")
Test_N
              
##8.4.3 Find Significant Taxa among Groups  

data<-log((abund_table+1)/(rowSums(abund_table)+dim(abund_table)[2]))
df<-as.data.frame(data)

df<-as.data.frame(abund_table/rowSums(abund_table))

KW_table <- data.frame()
for (i in 1:dim(df)[2]) {
  #run KW test for each bacterium
  KW_test <- kruskal.test(df[,i], g=df_CH_G$Group4)
  # Store the result in the data frame
  KW_table <- rbind(KW_table,
                    data.frame(id=names(df)[i],
                    p.value=KW_test$p.value
                                ))
  # Report number of bacteria tested
  cat(paste("Kruskal-Wallis test for ",names(df)[i]," ", i, "/", 
            dim(df)[2], "; p-value=", KW_test$p.value,"\n", sep=""))
}

#Check the data frame table
head(KW_table)

##8.4.4 Multiple Testing and E-value, FWER and FDR

##8.4.4.1 E-value
KW_table$E.value <- KW_table$p.value * dim(KW_table)[1]
KW_table$E.value

#check E-value in result data frame
head(KW_table)

##8.4.4.2 FWER
KW_table$FWER <- pbinom(q=0, p=KW_table$p.value,size=dim(KW_table)[1], lower.tail=FALSE)

#check the data frame table
head(KW_table)

##8.4.4.3 FDR
#order p-values from smallest to largest
KW_table <- KW_table[order(KW_table$p.value, decreasing=FALSE), ]
head(KW_table)

#calculate q-value
KW_table$q.value.factor <- dim(KW_table)[1] / 1:dim(KW_table)[1]
head(KW_table$q.value.factor)

KW_table$q.value <- KW_table$p.value * KW_table$q.value.factor
head(KW_table$q.value)

#check to see if q-value added to the result data frame
head(KW_table)

#set up alpha value
KW_alpha=0.05

#identify the last item of the ranked list with a q-value =< alpha 
last.significant.item <- max(which(KW_table$q.value <= KW_alpha))
last.significant.item

#display the chosen results
selected <- 1:5
#selected <- 1:last.significant.item
print(KW_table[selected,])

diff.taxa.factor <- KW_table$id[selected]
diff.taxa <- as.vector(diff.taxa.factor)
diff.taxa

