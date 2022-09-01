################################################################################
###Chapter 5: Power Analysis for Microbiome Data                               
###Yinglin Xia: September, 2018                                                                  
################################################################################

#Set working directory
setwd("E:/Home/MicrobiomeStatR/Analysis")

################################################################################
###5.2 Power Analysis for Testing Differences in Diversity                      
################################################################################

##5.2.2 Diversity Data for ALS Study
options(width=65,digits=4)
#Load abundance table
abund_table=read.csv("ALSG93AGenus.csv",row.names=1,check.names=FALSE) 
abund_table_t<-t(abund_table)

library(vegan)  
#using the diversity function (vegan package) to calculate Shannon index
shannon_genus <- diversity(abund_table,index="shannon",MARGIN=1)
#make a data frame of Shannon index
H<-diversity(abund_table_t, "shannon") 
df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon",length(H)))

#Obtain grouping information from sample data
df_H$Group <- with(df_H,
ifelse(as.factor(sample)%in% c("A11-28F","A12-28F","A13-28F","A14-28F","A15-28F","A16-28F"),c("G93m1"),
ifelse(as.factor(sample)%in% c("A21-28F","A22-28F","A23-28F","A24-28F","A25-28F","A26-28F"),c("WTm1"),   
ifelse(as.factor(sample)%in% c("C11-28F","C12-28F","C13-28F"),c("G93m4"),   
ifelse(as.factor(sample)%in% c("C21-28F","C22-28F","C23-28F"),c("WTm4"),
ifelse(as.factor(sample)%in% c("B11-28F","B12-28F","B13-28F","B14-28F","B15-28F","D11-28F","D12-28F","D13-28F","D14-28F"),c("BUm3to3.5"),
c("NOBUm3to3.5"))))))) 
df_H

install.packages("dplyr")
library(dplyr) 
df_H_G6 <- select(df_H, Group,value)
df_H_G93BUm3  <- filter(df_H_G6,Group=="BUm3to3.5"|Group=="NOBUm3to3.5")
df_H_G93BUm3 

library(ggplot2) 
#Split the plot into multiple panels
p<-ggplot(df_H_G93BUm3, aes(x=value))+
  geom_histogram(color="black", fill="black")+
  facet_grid(Group ~ .)
#Calculate the mean of each group
#Calculate the average Shannon diversity of each group using the package plyr
library(plyr)
mu <- ddply(df_H_G93BUm3, "Group", summarise, grp.mean=mean(value))
head(mu)

#add mean lines
p+geom_vline(data=mu, aes(xintercept=grp.mean, color="red"),
             linetype="dashed")


##5.2.3 Calculating Power or Sample Size Using R Function power.t.test
aggregate(formula = value ~ Group, 
          data = df_H_G93BUm3,
          FUN = mean)
#Group value
#1   BUm3to3.5 2.504
#2 NOBUm3to3.5 2.205

aggregate(formula = value ~ Group, 
          data = df_H_G93BUm3,
          FUN = var)

#Group   value
#1   BUm3to3.5 0.02892
#2 NOBUm3to3.5 0.04349

n1 <- 9
n2 <-7
s1<-sqrt(0.02892)
s1
#[1] 0.1701
s2<-sqrt(0.04349)
s2
#[1] 0.2085
s=sqrt((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2)
s
#[1] 0.05012

power.t.test(n=2:10,delta=2.504-2.205,sd=0.05012)
df_P <-data.frame(n,power)
df_P

n = c(2, 3, 4, 5, 6, 7, 8, 9, 10)  
power = c(0.8324, 0.9994, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000)

power <- sapply(n, function (x) power.t.test(n=x, delta=2.504-2.205,sd=0.05012)$power)
plot(n, power, xlab  = "Sample Size per group", ylab  = "Power to reject null",
     main="Power curve for\n t-test with delta = 0.05",
     lwd=2, col="red", type="l")

abline(h = 0.90, col="blue")


power.t.test(n=2:10,delta=2.504-2.205,sd=0.05012, type = "one.sample" )


################################################################################
###5.3 Power Analysis for Comparing Diversity Across More than Two Groups Using 
###    ANOVA	
################################################################################

##5.3.2 Calculating Power or Sample Size Using R Function pwr.avova.test	
df_H_G93WTm1N4 <- filter(df_H_G6,Group%in%c("G93m1","WTm1","G93m4","WTm4"))
df_H_G93WTm1N4 

fit = lm(formula = value~Group,data=df_H_G93WTm1N4)
anova (fit)

install.packages("pwr")
library(pwr)
pwr.anova.test(f= 0.23,k=4,n=45:55,sig.level=0.05)

################################################################################
###5.4 Power Analysis for Comparing a Taxon of Interest Across Groups 
################################################################################  

##5.4.2 Power Analysis Using R Function power.prop.test 
abund_table_Spe=read.csv("ALSG93AButyrivibrioSpecies.csv",row.names=1,check.names=FALSE)

abund_table_Spe<-t(abund_table_Spe)

#grouping information
grouping<-data.frame(row.names=rownames(abund_table_Spe),t(as.data.frame(strsplit(rownames(abund_table_Spe),"-"))))

grouping$Group <- with(grouping,ifelse(as.factor(X1)%in%c("B11","B12","B13","B14", "B15","D11","D12","D13","D14"),c("Butyrate"), c("Control")))

Butyrivibrio_G <-cbind(abund_table_Spe, grouping)
rownames(Butyrivibrio_G)<-NULL
Butyrivibrio_G

Butyrivibrio_G$Present <- ifelse((Butyrivibrio_G$Butyrivibrio > 0), "Present","Absent")
Butyrivibrio_G

library(MASS)       # load the MASS package 
tbl = table(Butyrivibrio_G$Group, Butyrivibrio_G$Present) 
tbl                 # the contingency table 

p1=1.0
p2=0.57
r=1
alpha=0.05
beta=0.20
(n2=(p1*(1-p1)/r+p2*(1-p2))*((qnorm(1-alpha/2)+qnorm(1-beta))/(p1-p2))^2)
ceiling(n2)
z=(p1-p2)/sqrt(p1*(1-p1)/n2/r+p2*(1-p2)/n2)
(Power=pnorm(z-qnorm(1-alpha/2))+pnorm(-z-qnorm(1-alpha/2)))


power.prop.test(n=10:20,  p1=1,  p2=.57,  sig.level=0.05, power=NULL,  alternative=c("one.sided"), strict  = FALSE)


################################################################################
###5.4.3. Power Analysis Using Chi-square and Fisher's Exact Test
################################################################################    
##4.3.2 Implementing Power Analysis Using the Function pwr.chisq.test  
install.packages("lsr")
library(lsr)
cramersV(tbl)

library(pwr)
pwr.chisq.test(w = 0.3833, N = 45:60, df = 1, sig.level = 0.05, power = NULL)
    

##5.4.3.3 Implementing Power Analysis Using the Functions power.fisher.test and power.exact.test 
install.packages("statmod")
library(statmod)
power.fisher.test(1.0,0.57,15,15,alpha=0.05, nsim=1000)

install.packages("Exact")
library(Exact)
power.exact.test(1.0, 0.57, 15, 15, method="Fisher")

  
################################################################################
###5.5 Comparing the Frequency of All Taxa Across Groups Using Dirichlet
###    -Multinomial Model                            
################################################################################

##5.5.3 Power and Size Calculations Using HMP Package
##5.5.3.1 Preparing Data Sets for Use of HMP Package
install.packages("HMP",repo="http://cran.r-project.org", dep=TRUE)
library(HMP)

setwd("F:/Home/MicrobiomeStatR/Analysis")

Buty=read.csv("ALSG93A3.5mButyrateGenus.csv",row.names=1,check.names=FALSE)
NOButy=read.csv("ALSG93A3.5mNoButyrateGenus.csv",row.names=1,check.names=FALSE) 

head(Buty)
head(NOButy)

Buty_t <- t(Buty)
NOButy_t<-t(NOButy)

head(Buty_t)
head(NOButy_t)

ncol(Buty_t)  # for the number of taxa
nrow(Buty_t)  # for the number of samples

ncol(NOButy_t)  # for the number of taxa
nrow(NOButy_t)  # for the number of samples


##5.5.3.2 Power and Size Calculations Using Taxa Composition Data Analysis
fit_Buty <- DM.MoM(Buty_t);fit_NOButy <- DM.MoM(NOButy_t)
fit_Buty 

numMC <- 1000

#The first number is the number of reads and the second is the number of subjects
nrsGrp1 <- rep(1000, 10)
nrsGrp2 <- rep(1000, 10) 
group_Nrs <- list(nrsGrp1, nrsGrp2)

alphap <- fit_Buty $gamma
pval1 <- MC.Xdc.statistics(group_Nrs, numMC, alphap, "hnull")
pval1

alphap <- rbind(fit_Buty$gamma, fit_NOButy$gamma)
pval2 <- MC.Xdc.statistics(group_Nrs, numMC, alphap)
pval2


##5.5.3.3 Power and Size Calculations Using Rank Abundance Distributions Data Analysis
filter_Buty<- Data.filter(Buty_t, "sample", 1000, 10)
head(filter_Buty)

filter_NOButy<- Data.filter(NOButy_t, "sample", 1000, 10)
head(filter_NOButy)

fit_Buty <- DM.MoM(filter_Buty);fit_NOButy <- DM.MoM(filter_NOButy);

fit_Buty$pi 
fit_NOButy$pi

fit_Buty$theta  
fit_NOButy$theta

numMC <- 1000 
# The irst number is the number of reads and the second is the number of subjects
nrsGrp1 <- rep(1000, 10);nrsGrp2 <- rep(1000, 10)
group_Nrs <- list(nrsGrp1, nrsGrp2)

pi0 <- fit_Buty$pi
group_theta <- c(0.007523, 0.01615)

pval1 <- MC.Xmc.statistics(group_Nrs, numMC, pi0, group.theta=group_theta, type="hnull")
pval1

group_pi <- rbind(fit_Buty$pi, fit_NOButy$pi)
pval2 <- MC.Xmc.statistics(group_Nrs, numMC, pi0, group_pi, group_theta)
pval2

#Generate the number of reads per sample
#The first number is the number of reads and the second is the number of subjects
nrsGrp1 <- rep(1000, 10) ;nrsGrp2 <- rep(1000, 10);
group_Nrs <- list(nrsGrp1, nrsGrp2)

pi0 <- fit_Buty$pi
group_theta <- c(0.007523, 0.01615)

#Computing size of the test statistics (Type I error)
group_theta <- c(fit_Buty$theta, fit_NOButy$theta)
pval1 <- MC.Xmc.statistics(group_Nrs, numMC, pi0, group.theta=group_theta, type="hnull")
pval1

#Computing Power of the test statistics (Type II error)
group_pi <- rbind(fit_Buty$pi, fit_NOButy$pi)
pval2 <- MC.Xmc.statistics(group_Nrs, numMC, pi0, group.pi=group_pi, group.theta=group_theta)
pval2

##5.5.4 Effect Size Calculation Using HMP Package
#Combine the data sets into a single list
group_data <- list(filter_Buty, filter_NOButy)
effect <- Xmcupo.effectsize(group_data)
effect

