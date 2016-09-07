#ps6 KEGG Gene Set Enrichment Analysis and Visualization of RNA-seq data
#Sept. 6, 2016


#setwd()
getwd()
setwd("/Users/jsmith/Documents/bi623/160906_PS6_PCRduplicateRemoval/")

#Question 1 
#Read in the edgeR differential gene expression results for 
#pipefish and format KEGG pathway information

#installthe R package GAGE.
source("https://bioconductor.org/biocLite.R")
biocLite("gage")

#install the R package GageData.
source("https://bioconductor.org/biocLite.R")
biocLite("gageData")

#install R package Pathview. 
source("https://bioconductor.org/biocLite.R")
biocLite("pathview")

#load the Gage library 
library("gage")

#load the library GageData
library("gageData")

#load the library Pathview. 
library("pathview")

#Read in the edgeR results.
kegg_pouch <- read.delim("pouch_RNAseq.tsv", sep = '\t', stringsAsFactors = FALSE)
head(kegg_pouch)


#Question 2
#Set up and perform GSEA for all KEGG pathways using gage

#Load the KEGG orthology data from the gageData database
data(kegg.sets.ko)
data(sigmet.idx.ko)
kegg.sets.ko <- kegg.sets.ko[sigmet.idx.ko]
kegg.sets.ko

#what class of object is kegg.sets.ko, and what kind of info does it contain?
class(kegg.sets.ko)
#List class 

#make a variable to hold only the fold_changes column values
pouch_foldchanges <- kegg_pouch$logFC

#names function will set the names for each 
#fold-change value as the KO_ID 
names(pouch_foldchanges) <- kegg_pouch$ko_ID
head(pouch_foldchanges)

#run gage on vector of log2 fold changes to test for 
#enrichment of genes with extreme values 
pouch_test <- gage(pouch_foldchanges, gsets = kegg.sets.ko, 
                   same.dir = FALSE)
#this "same.dir" performed a “two-tailed” test, meaning a significantly 
#enriched KEGG pathway can be defined by extreme fold changes 
#in either direction (up or down).

#use lapply to print the top enties of the pouch_test
lapply(pouch_test, head)
#use head function to list 30 entries  
head(pouch_test$greater, 30)

#columns of pouch_test$greater provide us our hypothesis test
#results, including raw p-values, corrected p-values, 
#sizes of the gene sets, etc.


#use subset() with a logical operator to return only  
#pathways such that our false discovery rate (FDR) = 0.1
subset(pouch_test$greater, pouch_test$greater[,4]<0.1)



#Question 3
#Visualize pregnancy fold change magnitudes for all genes 
#in the “coagulation and complement” KEGG pathway


#Pull out the "coagulaton and complement" pathway and its
#ko ID from the gage() output by simple indexing and the substr() function.
pouch_pathways <- rownames(pouch_test$greater)[4]

#start and stop are the character positions in the 
#string "KO00000" for KO ids
pouch_ids <- substr(pouch_pathways, start = 1, stop = 7)
pouch_ids


#Draw the pathway with a color scale that reflects log2 
#fold change for each gene. 
pathview(gene.data=pouch_foldchanges,
         pathway.id=pouch_ids, species="ko", new.signature=FALSE,
         trans.fun = list(gene = NULL, cpd = NULL),
         low = list(gene = "green", cpd = "blue"),
         mid = list(gene = "yellow", cpd = "gray"),
         high = list(gene = "red", cpd = "yellow"), na.col = "transparent")
#path is saved to a file with KO_ID.pathview.png


#Question 4
#Visualize multi-genic expression patterns in pregnant and
#non-pregnant pouch tissues using non-metric 
#multidimensional scaling.

#install the R packages vegan and MASS
source("https://bioconductor.org/biocLite.R")
biocLite("vegan")
source("https://bioconductor.org/biocLite.R")
biocLite("MASS")

#load the libraries
library("MASS")
library("vegan")

#Read in the normalized pipefish pouch expression values,
#a comma-separated file.
pipe_TMMvals <- read.csv("pouch_TMM_values.csv", header = TRUE, row.names = 1)
dim(pipe_TMMvals)

#transpose the dataset 
pipe_TMMvals <- t(pipe_TMMvals)
head(pipe_TMMvals)

#compute a dissimilarity matrix (distances between pipefish
#samples in “transcript space”) using the default for 
#vegdist(), which is a metric called Bray-Curtis 
#dissimilarity.
pipe.dis <- vegdist(pipe_TMMvals)

#The fit of the new, 2-D distances to the original, 
#can be visualized in a stress plot which plots
#the 2-D versus 15,252-D distance ranks.
pipe.mds0 <- isoMDS(pipe.dis, k=2)
#converged 

#create a pdf
pdf("StressPlot2D.pdf")
#make a stress plot
stressplot(pipe.mds0, pipe.dis, main = "Shepard Diagram for nMDS of Gene Expression \n in Pipefish ")
#close the file
dev.off()

#construct a data frame that links the pouch sample IDs
#with pregnancy status
Targets <- as.data.frame(rownames(pipe_TMMvals))
Targets$PregStat <- factor(c(rep("preg",6),rep("nonpreg",6)))
colnames(Targets)<- c("ID", "PregStat")

#create a PDF
pdf("ordinationPlot.pdf")

#define the ordination plotting parameters. 
par(mgp = c(2.5,1,0))
preg=as.character(Targets$PregStat)

#The ordiplot() function will produce an empty space.
fig <- ordiplot(pipe.mds0, main="Brood Pouches in Transcript Space", 
                ylab="nMDS Dimension 2", xlab="nMDS Dimension 1", 
                font.lab=2, font.axis=2, cex.axis=.7, type="none", 
                cex.main=1, xlim=(c(-.2,0.2)))

#add “confidence ellipses” and the individual samples as 
#points using the ordiellipse() and points() functions.
ordiellipse(pipe.mds0,groups=preg,label=FALSE, lwd=2, 
            show.groups=preg[1:6], col="darkorange4", 
            draw="lines")
ordiellipse(pipe.mds0,groups=preg,label=FALSE, lwd=2, 
            show.groups=preg[7:12], col="dodgerblue4", 
            draw="lines")

#the blue squares are non-pregnant males and the orange circles 
#are pregnant males.
points(fig, "sites", pch=c(rep(19,6),rep(15,6)),
       col=c(rep("darkorange3",6),rep("cornflowerblue",6)),
       cex=1.5)
dev.off()

#Use the ordination to generate 3 nMDS dimensions
#Produce the same type of plot as above, but plot Dimension 2 vs. Dimension 3.

#k=3 for three dimmensions 
pipe.mds3 <- isoMDS(pipe.dis, k=3)

#make a variable to hold the values for dimmensions 2 and 3
dims2and3 <- pipe.mds3$points[,2:3]

#pdf
pdf("stressplot3D.pdf")
#create a stress plot for nMDS with three dimmensions
stressplot(pipe.mds3, pipe.dis, main = "Shepard Diagram for 3D nMDS of Gene Expression \n in Pipefish")
dev.off()

#link pouch_sample IDs with pregnanancy results 
Targets3 <- as.data.frame(rownames(pipe_TMMvals))

#make the factors "preg" and "nonpreg"
Targets3$PregStat <- factor(c(rep("preg",6),rep("nonpreg",6)))
colnames(Targets3) <- c("ID","PregStat")

#pdf 
pdf("OrdinationPlot3D.pdf")
#set up parameters for the ordination plot
par(mgp=c(2.5,1,0))
preg=as.character(Targets3$PregStat)

#create an empty plot with dimmensions 2 and 3 
fig3Dim <- ordiplot(dims2and3, main="Brood Pouches in Transcript Space",
                    ylab="nMDS Dimension 3", xlab="nMDS Dimension 2", 
                    font.lab=2, font.axis=2, cex.axis=.7, type="none", 
                    cex.main=1, xlim=(c(-.2,0.2)))

#create ellipses for dimmension 2 and 3 
ordiellipse(dims2and3,groups=preg,label=FALSE, lwd=2, show.groups=preg[1:6], col="purple", draw="lines")
ordiellipse(dims2and3,groups=preg,label=FALSE, lwd=2, show.groups=preg[7:12], col="navyblue", draw="lines")
#plot the points
points(fig3Dim, "sites", pch=c(rep(19,6),rep(15,6)), col=c(rep("purple",6),rep("navyblue",6)), cex=1.5)
dev.off()



#Question 5
#permutational Multivariate Analysis of Variance (perMANOVA)
#to test for multivariate transcriptional differences 
#between pregnant and non- pregnant males

#The adonis() function  runs the perMANOVA,
#adonis() requires otu.env as a parameter
otu.env <- Targets
adonis(pipe.dis ~ PregStat, otu.env, perm=999)




#question 6
#Constructing a heatmap with clustering dendrograms for 
#Coagulation and Complement Cascade KEGG pathway genes.


#install the packages gplots, RColorBrewer, and dendextend
biocLite("gplots")
biocLite("RColorBrewer")
biocLite("dendextend")

#load the libraries
library("gplots")
library("RColorBrewer")
library("dendextend")

#read in the file with CPM expression data for one pipefish gene
#that was maped to KEGG pathway "Coagulation and Complement Cascade"
pouch_compcoag <- read.delim("CompCoag_pouch_multivar.tsv", sep = '\t', 
                             row.names = 1, header = FALSE)

#add column names to the dataset
colnames(pouch_compcoag) <- c("KO","name","P9","P8","P7","P6","P5",
                              "P3", "NP11","NP10","NP4","NP3","NP2","NP1")

#check that names are added and characterize dataset
head(pouch_compcoag)
dim(pouch_compcoag)

#Define a vector of gene names
names_compcoag <- pouch_compcoag$name
names_compcoag

#The addition of 0.01 is to avoid taking the log of zero
#in the next step. 
pouch_compcoag <- pouch_compcoag[,3:14] + 0.01
head(pouch_compcoag)

#define a new data frame consisting of the CPM columns and log-transform (base 2) each CPM value + 0.01.
pouch_compcoag <- as.data.frame(log2(pouch_compcoag))
class(pouch_compcoag)
dim(pouch_compcoag)

#transpose for row-wise scaling
pouch_compcoag.n <- scale(t(pouch_compcoag))

#put back in original orientation
pouch_compcoag.tn <- t(pouch_compcoag.n)
head(pouch_compcoag.tn)
class(pouch_compcoag.tn)
#matrix class 

#calculate multivariate dissimilarity for all sample pairs,
#using Euclidean Distance
compcoag.d1 <- dist(pouch_compcoag.n, method = "euclidean", diag = FALSE,
                    upper = FALSE)

#view the euclidean distances calculated
round(compcoag.d1,3)

#calculate multivariate dissimilarity for all gene pairs.
compcoag.d2 <- dist(pouch_compcoag.tn,method = "euclidean",
                    diag = FALSE, upper = TRUE)

#cluster samples using Ward linkage clustering
compcoag.c1 <- hclust(compcoag.d1, method = "ward.D2", members = NULL)

#cluster genes using Ward linkage clustering
compcoag.c2 <- hclust(compcoag.d2, method = "ward.D2", members = NULL)


#create a pdf 
pdf("DendrogramsScaled.pdf")
#parameters Make 2 rows, 1 col plot frame and shrink labels
par(mfrow=c(2,1),cex=0.5)  
#plot dendrograms based on the clustering
plot(compcoag.c1, xlab = "Sample Name",
     main = "Ward Clustering of Pregnant and Non-Pregnant Pipefish Samples")
plot(compcoag.c2, xlab = "Ensembl Gene ID",
     main = "Ward Clustering of Gene Expression Levels from Pregnant and Non-Pregnant Pipefish" )
dev.off()

#print the order of samples (left to right) in the tree,
#in case we want to rotate the dendrogram about nodes and
#re-order them later.
compcoag.c1$order

#set the color scale for the heatmap to 299 increments. 
compcoag_pal <- colorRampPalette(c("green","yellow","red")) (299)


#set parameters for the heatmap 
par(cex.main=0.5, cex.axis=0.5, font=2, font.axis=2) 
#plot the heatmap with heatmap.2
heatmap.2(pouch_compcoag.tn, 
          Colv=rotate(as.dendrogram(compcoag.c1), 
          order=c(12,11,10,9,8,7,6,5,4,3,2,1)),
          Rowv=as.dendrogram(compcoag.c2), 
          labRow=names_compcoag, 
          density.info="none",
          trace="none",
          scale="none",
          col = compcoag_pal, cexRow=0.5,cexCol=0.75, margins=c(3,13), lwid=c(.8,3), lhei=c(.8,3), srtCol=45, adjCol=c(1,1),
          keysize=1.3)

#change the values to adjust the appearance of the image. 
par(cex.main=0.5, cex.axis=0.5, font=2, font.axis=2) 
heatmap.2(pouch_compcoag.tn, 
          Colv=rotate(as.dendrogram(compcoag.c1), 
          order=c(12,11,10,9,8,7,6,5,4,3,2,1)),
          Rowv=as.dendrogram(compcoag.c2), 
          labRow=names_compcoag, 
          density.info="none",
          trace="none",
          scale="none",
          col = compcoag_pal, cexRow=0.5,cexCol=0.75, 
          margins=c(3.5,13), lwid=c(.8,3), lhei=c(.8,3), 
          srtCol=45, adjCol=c(1,1),
          keysize=1.0)

#Repeat the heatmap starting without scaling
#transpose the dataset
pouch_compcoag.NoScalen <- t(pouch_compcoag)

#put back in original orientation
pouch_compcoag.NoScaleTn <- t(pouch_compcoag.NoScalen) 

#NS = no scale
#calc dissimilarity of all sample pairs
compcoag.d1.NS  <- dist(pouch_compcoag.NoScalen, method = "euclidean", 
                        diag = FALSE, upper = FALSE) 
#look at distances
round(compcoag.d1.NS,3)

#calc dissimilarity of gene pairs
compcoag.d2.NS <- dist(pouch_compcoag.NoScaleTn, method = "euclidean", 
                       diag = FALSE, upper = TRUE)

#cluster sample, then genes
compcoag.c1.NS <- hclust(compcoag.d1.NS, method = "ward.D2", members = NULL)
compcoag.c2.NS <- hclust(compcoag.d2.NS, method = "ward.D2", members = NULL) 

#pdf
pdf("dendrogramsNotScaled.pdf")
par(mfrow=c(2,1),cex=0.5)
plot(compcoag.c1.NS, xlab = "Sample Name",
     main = "Ward Clustering of Pregnant and Non-Pregnant Pipefish Samples without Scaling")
plot(compcoag.c2.NS, xlab = "Ensembl Gene ID",
     main = "Ward Clustering of Gene Expression Levels from Pregnant and Non-Pregnant Pipefish without Scaling")
dev.off()

#order of the clustering
compcoag.c1.NS$order

#create a heatmap 
par(cex.main=0.5, cex.axis=0.5, font=2, font.axis=2) 
heatmap.2(pouch_compcoag.NoScaleTn, 
          Colv=rotate(as.dendrogram(compcoag.c1.NS), 
          order=c(12,11,10,9,8,7,6,5,4,3,2,1)), 
          Rowv=as.dendrogram(compcoag.c2.NS), 
          labRow=names_compcoag, 
          density.info="none",
          trace="none",
          scale="none",
          col = compcoag_pal, cexRow=0.5,cexCol=0.75, 
          margins=c(3.5,13), lwid=c(.8,3), lhei=c(.8,3), 
          srtCol=45, adjCol=c(1,1),
          keysize=1.0)


#Quesheption 7 
#Constructing heatmaps with clustering dendrograms for 
#stickleback gene expression data.

#For the first heatmap, subset only genes that have the
#Genome_Loc value “groupXIX” AND have a 
#Gene_Start_Pos value between 6000000 and 12000000.

#read in the file.
stickleback <- read.delim("stickleback_CPM.tsv", sep = '\t', 
                               row.names = 1, header = TRUE)

#subset to only chromosome 19 
stickleback_chr19 <- subset(stickleback, stickleback$Genome_Loc == "groupXIX")

#subset to bp 6,000,000 and 12,000,000.
stickleback_chr19.bp <- subset(stickleback_chr19, 
                               stickleback_chr19$Gene_Start_Pos >= 6000000 & stickleback_chr19$Gene_Start_Pos <= 12000000 )

#add 0.01 to each cpm value 
stickleback_chr19.bp <- stickleback_chr19.bp[,3:86] + 0.01

#log2 transform the cpm values
stickleback_chr19.bp <- as.data.frame(log2(stickleback_chr19.bp)) 

#create a vector of sample names for the heat map
names_stickleback <- colnames(stickleback_chr19.bp)

#mean-center and range the data, as well as transpose 
stickleback_chr19.bp.n <- scale(t(stickleback_chr19.bp))

#put back in original orientation
stickleback_chr19.bp.tn <- t(stickleback_chr19.bp.n)

#calc the multivariate dissimilarity for sample pairs, using euclidean distance 
stickleback.d1 <- dist(stickleback_chr19.bp.n, method = "euclidean", diag = FALSE, upper = TRUE)

#calc multivariate dissimilarity for all gene pairs
stickleback.d2 <- dist(stickleback_chr19.bp.tn, method = "euclidean", diag = FALSE, upper = TRUE)

#cluster samples, then genes using ward linkage clustering
stickleback.c1 <- hclust(stickleback.d1, method = "ward.D2",  members = NULL)
stickleback.c2 <- hclust(stickleback.d2, method = "ward.D2", members = NULL)

#Check out the dendrograms
par(mfrow=c(2,1), cex=0.5)
plot(stickleback.c1)
plot(stickleback.c2, cex=0.1)

#print the order of the samples of M and F
stickleback.c1$order

#set the color palette
col_palette <- colorRampPalette(c("green", "yellow", "red")) (n=299)

#Create the heatmap
par(cex.main=0.5, cex.axis=0.5, font=2, font.axis=2)
heatmap.2(stickleback_chr19.bp.tn,
          Colv=as.dendrogram(stickleback.c1),
          Rowv=as.dendrogram(stickleback.c2),
          labRow= FALSE, 
          labCol = names_stickleback, 
          density.info="none",
          trace="none",
          scale="none",
          col = col_palette, cexRow=0.5,cexCol=0.75,
          margins=c(4,1), lwid=c(.8,3), lhei=c(.8,3), 
          srtCol=45, adjCol=c(1,1),
          keysize=1.3)

#For the second heatmap, take a random sample of genes the same size
#as the number of genes in chr19 subset.
#nrow is the number of rows, thus genes, in the counts file
randomSample <- stickleback[sample(nrow(stickleback),314),]

#add 0.01 to each CPM
randSampleLog2 <- randomSample[,3:86] + 0.01
#log2 transform the CPM values
randSampleLog2 <- as.data.frame(log2(randSampleLog2))

#sclae and transpose the data
randomSample.n <- scale(t(randSampleLog2))

#put back in the same order 
randomSample.tn <- t(randomSample.n)

#calculate the multivaraite distances, euclidean, for samples
randomSample.d1 <- dist(randomSample.n, method = "euclidean", diag = FALSE, upper = FALSE)

#calculate the multivaraite distances, euclidean, for genes
randomSample.d2 <- dist(randomSample.tn, method = "euclidean", diag = FALSE, upper = TRUE)

#cluster the samples, then genes 
randomSample.c1 <- hclust(randomSample.d1, method = "ward.D2", members = NULL)
randomSample.c2 <- hclust(randomSample.d2, method = "ward.D2", members = NULL)

#view the cluster dendrograms 
par(mfrow=c(2,1),cex=0.5) 
plot(randomSample.c1)
plot(randomSample.c2)

#create a vector of names for columns
rand_names <- colnames(randSampleLog2)

#create heatmap of randomized sample
par(cex.main=0.5, cex.axis=0.5, font=2, font.axis=2) 
heatmap.2(randomSample.tn,
          Colv=as.dendrogram(randomSample.c1),
          Rowv=as.dendrogram(randomSample.c2),
          labRow= FALSE, 
          labCol = rand_names, 
          density.info="none",
          trace="none",
          scale="none",
          col = col_palette, cexRow=0.5,cexCol=0.75,
          margins=c(4,1), lwid=c(.8,3), lhei=c(.8,3), 
          srtCol=45, adjCol=c(1,1),
          keysize=1.3)



