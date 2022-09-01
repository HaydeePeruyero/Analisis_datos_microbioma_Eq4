################################################################################
###Chapter 4: Introduction to R, RStudio and ggplot2
###Yinglin Xia: September, 2018                                                                 
################################################################################
 
################################################################################ 
##4.1.Introduction to R and RStudio
################################################################################

##citation
citation () 
citation ("ALDEx2")   
RStudio.Version() 

##4.1.1 Installing R, RStudio, and R Packages    
install.packages("ALDEx2")  
library(ALDEx2) 
installed.packages()[1:5,]

a<-installed.packages()
packages<-a[,1] 
is.element("ALDEx2", packages) 
 

##4.1.2 Set Working Directory in R 
getwd()  
setwd("E:/Home/MicrobiomeStatR/Analysis") 
getwd() 
setwd("~/")
getwd()


##4.1.3 Data Analysis Through RStudio
setwd("E:/Home/MicrobiomeStatR/Analysis")

#Boxplot of writing score by gender 
boxplot(write ~ female,data=hsb2demo, main="High School Students Data", 
        xlab="Gender", ylab="Writing score by gender")

  
##4.1.4	Data Import and Export
tab <- read.table("genus.csv", header=TRUE,row.names=1, sep=",") 
tab <- read.table("genus.txt", header=TRUE,row.names=1, sep="\t") 
raw<- "https://raw.githubusercontent.com/chvlyl/PLEASE/master/1_Data
/Raw_Data/MetaPhlAn/PLEASE/G_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result.xls"  
tab <- read.table(raw,sep='\t',header=TRUE,row.names = 1,check.names=FALSE,stringsAsFactors=FALSE)

tab <- read.delim("genus.txt", header=T, row.names=1) 

tab <- read.csv('table.csv', head = T, row.names = 1, sep = ',', dec = '.')
tab <- read.csv2 ('table.csv',head = T, row.names = 1,sep = ';', dec = ',')
tab <- read.table (file = 'table.csv', head = T, row.names = 1, sep = ';', dec = ',')

install.packages ("gdata")
library(gdata)
tab <- read.xls("table.xlsx",sheet=1,header=TRUE)
tab <- read.xls("table.xlsx", sheet=1,perl="C:/Perl64/bin/perl.exe")

install.packages ("XLConnect")
library (XLConnect)
tab <- readWorksheetFromFile(file = 'table.xlsx', sheet = 1, header = T, rownames = 1)

write.table(genus, file="genus_out.csv", quote=FALSE, row.names=FALSE,sep="\t") 
write.table(genus, file="genus_out.txt", quote=FALSE, col.names=TRUE,sep=",") 

##4.1.5	Data Manipulation
data()
attach(iris)
head(iris)

?data.frame 

#Create data frame using column indices
df <- iris[,c(1,2,3)]
head(df)

#Create data frame using column indices with sequences
df <- iris[,c(1:2,4:5)]
head(df)

#Create data frame using subset() and column indices
df<- subset(iris, select=c(1,2, 4:5))
head(df)

#Create data frame using subset() and column names
df <- subset(iris, select=c("Sepal.Width", "Petal.Length", "Petal.Width"))
head(df)

#Create data frame by selecting column names
df <- iris[,c("Sepal.Width", "Petal.Length", "Petal.Width")]
head(df)

#Create data frame using data.frame()
df <- data.frame(iris$Sepal.Width, iris$Petal.Length, iris$Petal.Width)
head(df)

#Create data frame using c() and data.frame() manually
Sepal.Width = c(3.5, 3.0, 3.2, 3.1,3.6,3.9) 
Petal.Length = c(1.4,1.4,1.3,1.5,1.4,1.7) 
Petal.Width = c(0.2,0.2,0.2,0.2,0.2,0.4) 
df = data.frame(Sepal.Width,Petal.Length,Petal.Width) 
df

head(iris) 
attributes(iris) 
dim(iris) 
nrow(iris)    
ncol(iris)
length(iris[,"Species"])

#check column or row names 
colnames(iris)
rownames(iris)

print(iris)
Species <- iris[,"Species"]
Species

iris[1,3]  
iris["1", "Petal.Length"] 
head(iris[,-c(4:5)])

tab=read.csv("VdrGenusCounts.csv",row.names=1,check.names=FALSE)
#Check total zeros in the table
sum(tab == 0) 
#Check how many non-zeros in the table
sum(tab != 0)


ng <- layout(matrix(c(1,3,2,3),2,2, byrow=TRUE), widths=c(5,2), height=c(3,4))
layout.show(ng) 

options(width=65,digits=4)

##4.1.6	Simple Summary Statistics
summary(iris) 
iris_1 <- (iris[,-5]) 
head(apply(iris_1, 1, mean))
apply(iris_1, 2, mean)
apply(iris_1, 2, mean,na.rm = TRUE)

tab_perc <- apply(tab, 2, function(x){x/sum(x)})  
tab_perc <- apply(tab[,1:ncol(tab)-1], 2, function(x){x/sum(x)})
tab_p1 <- tab[apply(tab_perc, 1, max)>0.01,]
tab_p2 <- tab[apply(tab_perc, 1, min)>0.01,]
head(tab_p2)

count <- 1
tab_min <- data.frame(tab[which(apply(tab, 1, function(x){mean(x)}) > count),], check.names=F) 
 
cutoff = .5
tab_d5 <- data.frame(tab[which(apply(tab, 1, function(x){length(which(x != 0))/length(x)}) > cutoff),])

count = 500
tab_c500 <- data.frame(tab[which(apply(tab, 1, function(x){sum(x)}) > count),])

##4.1.7	Other useful R functions
#Converting data frames
iris_t <-t(iris) 
iris_t[1:5,1:6]

#Sorting and ordering data frames
iris_2 <- (iris[,-c(3:5)])
sorted <- sort(iris_2$Sepal.Length)
ordered <- order(iris_2$Sepal.Length)
new_iris<- data.frame(iris_2,sorted,ordered)
head(new_iris)

#Sorting and ordering data frames
rev_iris <- rev(sort(iris_2$Sepal.Length))
head(rev_iris)

head(iris[order(Sepal.Length),])
head(iris[order(iris[,'Sepal.Length']),])

group <- ifelse(iris$Petal.Length < 4,1,2) 
group_s <- ifelse(iris$Species %in% "setosa",1,
                  ifelse(iris$Species %in% "versicolor",2,3))

##Load abundance table;
tab=read.csv("VdrGenusCounts.csv",row.names=1,check.names=FALSE)
tab_t<-t(tab)
head(tab_t)[1:3,c("Tannerella", "Lactococcus", "Lactobacillus")]

strsplit<-data.frame(row.names=rownames(tab_t),t(as.data.frame(strsplit(rownames(tab_t),"_"))))
head(strsplit)

re <- gsub(pattern = "\\$", replacement = ".", "metacharacters$uses$in$regular$expressions")
re

df = data.frame(how...to...interpret...metacharacters = c(1, 2), in...regular...expressions = c(1,2))  
df
names(df) <- gsub(pattern = "\\.+", replacement = ".", x = names(df))
names(df)


tab_ts <- read.table("tongue_saliva.txt", header=T, row.names=1, sep="\t")
tab_tts <-t(tab_ts)
rownames(tab_tts)[1:3]
rownames(tab_tts)[201:203]

rownames(tab_tts) <- gsub("sa_.+", "saliva", rownames(tab_tts))
rownames(tab_tts) <- gsub("td_.+", "tongue", rownames(tab_tts))
rownames(tab_tts)[1:3]
rownames(tab_tts)[201:203]

tax <- tab_ts$tax.0
tax[1:3]

tax_1 <- gsub(".+__", "", tax)
tax_1[1:3]


data(iris)
head(iris)

grep("[wd]", c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width","Species"))
grep("[wd]", c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width","Species"), value = TRUE)
grep("Width", c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width","Species"), value = TRUE,fixed =TRUE)

group <- data.frame(c(rep(1,length(grep("drySt", colnames(tab)))),
                      rep(2, length(grep("CeSt", colnames(tab))))))

group_1 <- data.frame(c(rep("fecal",length(grep("drySt", colnames(tab)))),
                      rep("cecal", length(grep("CeSt", colnames(tab))))))

fecal <- grep("drySt", colnames(tab))
cecal <- grep("CeSt", colnames(tab))
fecal
cecal

conds <- c(rep("fecal", 8), rep("CeSt", 8))
conds

fecal<- colnames(tab)[grep("drySt", colnames(tab))] # fecal samples
cecal <- colnames(tab)[grep("CeSt", colnames(tab))] # cecal samples
df_mice <- data.frame(tab[,fecal], tab[,cecal]) # make a data frame
head(df_mice)


################################################################################
##4.2. Introduction to the dplyr Package
################################################################################
setwd("E:/Home/MicrobiomeStatR/Analysis")
tab <- read.csv("hsb2demo.csv")
head(tab)

install.packages("dplyr")
library(dplyr)

head(select(tab, id, write))

tab %>% 
select(id,write) %>% 
head 

#Select columns: id, read, write and math
head(select(tab, id, read, write, math))

#Select all columns between read and socst (inclusive)
head(select(tab, read:socst))

#Select all columns except female
head(select(tab, -female))

#Select all columns except those from female to prog (inclusive)
head(select(tab, -(female:prog )))

#Select all columns that start with the character string "s"
head(select(tab, starts_with("s")))


#Filter the rows for students with reading score greater than or equal 70.
filter(tab, read >= 70)

#Filter the rows for students with both reading and math scores greater than or equal 70
filter(tab, read >= 70, math >= 70)

#Re-order by read and write
head(arrange(tab, id, read, write))

#Use desc() to order a column in descending order
head(arrange(tab, desc(read)))

#To re-order rows by a particular column(female)
tab %>% arrange(female) %>% head

#Select three columns id, gender, read from tab
#Arrange the rows by the gender and then by read
#Then return the head of the final data frame
tab%>%select(id, female, read) %>%
          arrange(female, read) %>% 
          head
          
#Filter the rows for read with score greater or equal to 70
tab %>% select(id, female, read) %>%
  arrange(female, read) %>% 
   filter(read >= 70)

#Arrange the rows for read in a descending order
tab %>% select(id, female, read) %>%
   arrange(female, desc(read)) %>% 
   filter(read >= 70)

#Create new columns using mutate()
#Calculate average read and write scores
head(mutate(tab, avg_read = sum(read)/n()))

#To keep only the new variables, use transmute()
head(transmute(tab,avg_read = sum(read)/n()))
 
#Create new columns using mutate() and pipe operator
tab %>% mutate(avg_read = sum(read/n())) %>%
  head

#To collapses a data frame to a single row.
summarise(tab, avg_read = mean(read, na.rm = TRUE))


#Create summaries of the data frame using summarise() and pipe operator
tab %>% summarise(avg_read = mean(read), 
                   min_read = min(read),
                   max_read = max(read),
                   n = n())

#First group by gender, and then get the summary statistics of reading by gender
by_gender <- group_by(tab, female)
read_by_gender <- summarise(by_gender,
                             n = n(),
                             avg_read = mean(read, na.rm = TRUE),
                             min_read = min(read,na.rm = TRUE),
                             max_read = max(read,na.rm = TRUE))
read_by_gender

#Create summaries of the data frame using summarise() and pipe operator
tab %>% group_by(female) %>%
        summarise(n = n(),
                   avg_read  = mean(read), 
                   min_read = min(read),
                   max_read = max(read))


#Use sample_n() to sample a fixed number
sample_n(tab, 5)

#Use sample_frac() to sample a fixed fraction
sample_frac(tab, 0.02)

#Use replace = TRUE to perform a bootstrap sampling
sample_n(tab, 5,replace = TRUE)

################################################################################
###4.3.	Introduction to ggplot2
################################################################################

##4.3.1	ggplot2 and the Grammar of Graphics
citation("ggplot2")


##4.3.2	Simplify Specifications in Creating a Plot Using ggplot() 
library(ggplot2)
ggplot() +
  layer(
    data = iris, mapping = aes(x = Sepal.Width, y = Sepal.Length),
    geom = "point", stat = "identity", position = "identity"
  ) +
  scale_y_continuous() +
  scale_x_continuous() +
  coord_cartesian()


ggplot() +
  layer(
    data = iris, mapping = aes(x = Sepal.Width, y = Sepal.Length),
    geom = "point"
  )

ggplot(iris, aes(Sepal.Width, Sepal.Length)) +
  geom_point()

ggplot(iris, aes(Sepal.Width, Sepal.Length)) +
  geom_point() +
  stat_smooth(method = lm) +
  scale_x_log10() +
  scale_y_log10()

ggplot() +
  layer(
    data = iris, mapping = aes(x = Sepal.Width, y = Sepal.Length),
    geom = "point", stat = "identity", position = "identity"
  ) +
  layer(
    data = iris, mapping = aes(x = Sepal.Width, y = Sepal.Length),
    geom = "smooth", position = "identity",
    stat = "smooth", method = lm
  ) +
  scale_y_log10() +
  scale_x_log10() +
  coord_cartesian()


qplot(Sepal.Width, Sepal.Length, data = iris,
      geom = c("point", "smooth"),
      method = "lm", log = "xy")
      

##4.3.3	Creating a Plot Using ggplot() 
library(ggplot2)

p <- ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) 
# Sepal.Width and Sepal.Length are columns in iris dataframe
p   
summary(p)

#Add scatterplot geom (layer1)
p1 <- p + geom_point()  
summary(p1)

ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) + geom_point()

#Add smoothing geom (layer2)
p2 <- p1 + geom_smooth(method="lm") 
p2  
summary(p2)

#set se = FALSE to turn off confidence bands
p1 + geom_smooth(method="lm", se = FALSE) 

p3 <- ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) + 
#Add scatterplot geom (layer1)
geom_point(col="blue", size=3) + 
#Add smoothing geom (layer2)
geom_smooth(method="lm",col="red",size=2) 
p3 

p4 <- ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) + 
  #Add scatterplot geom (layer1)
  geom_point(aes(col=Species), size=3) + 
  #Add smoothing geom (layer2)
  geom_smooth(method="lm",col="red",size=2) 
p4

p5 <- p4 + coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) # zooms in
plot(p5)

#Add Title and Labels using labs()
p6 <- p5 + labs(title="Sepal width vs sepal length", subtitle="Using iris dataset",
 y="Length of Sepal", x="Width of Sepal")
print(p6)#Or plot(p6)

#Add Title and Labels using ggtitle(), xlab() and ylab()
p7 <-p5 +  ggtitle("Sepal width vs sepal length", subtitle="Using iris dataset")
 + ylab("Length of Sepal") + xlab("Width of Sepal")
print(p7)

library(ggplot2)
ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
geom_point(aes(col=Species), size=3) + 
geom_smooth(method="lm",col="red",size=2) +
coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
labs(title="Sepal width vs sepal length", subtitle="Using iris dataset", 
y="Length of Sepal", x="Width of Sepal")

 
#Spliting plots by rows
ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
geom_point(aes(col=Species), size=3) + 
geom_smooth(method="lm",col="red",size=2) +
coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
# Add Facet Grid
facet_grid(Species ~.) 


#Spliting plots by columns
ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
geom_point(aes(col=Species), size=3) + 
geom_smooth(method="lm",col="red",size=2) +
coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
# Add Facet Grid
facet_grid(.~ Species) 

 
#Spliting plots by columns
ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
geom_point(aes(col=Species), size=3) + 
geom_smooth(method="lm",col="red",size=2) +
coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
# Add Facet Grid
facet_grid(.~ Species, margin=TRUE) 


#Facet Wrap
#Spliting plots by columns
ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) +
geom_point(aes(col=Species), size=3) + 
geom_smooth(method="lm",col="red",size=2) +
coord_cartesian(xlim=c(2.2,4.2), ylim=c(4, 7)) +
#Add Facet Wrap
facet_wrap(~ Species, nrow=2) 

