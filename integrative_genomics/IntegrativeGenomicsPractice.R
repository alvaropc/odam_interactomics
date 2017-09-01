########
#Alvaro Ponce Cabrera
########

setwd('~/Practice')

###
# 1 From the TCGA project choose a specific 
#  cancer type and report how many cases for each technology were processed

#Colon adenocarcinoma	
#Cases per technology
# BCR 1098
#clinical 1097
#CN 1089
#LowP 19
#Methylation 1097
#mRNA 526
#mRNAseq 1093
#miR 0
#miRseq 1078 cases
#RPPA 887
#MAF 977
#raw MFA  0

###
## 2 Using package RTCGAToolbox download from GDAC Firehose three different sorts of data
## (for instance, RNASeq, CNV and miRNA). 
####

library(RTCGAToolbox)

# Valid aliases
getFirehoseDatasets()

# Valid stddata runs
getFirehoseRunningDates()

# Valid analysis running dates (will return 3 recent date)
getFirehoseAnalyzeDates(last=3)


# Get data


# READ some specific data
# readData = getFirehoseData (dataset="COADREAD", runDate="20150821",forceDownload = TRUE, Clinic=TRUE,
#                             RNAseq2_Gene_Norm=TRUE,Methylation  = TRUE,
#                             mRNA_Array = TRUE,miRNASeq_Gene=TRUE)
# save(readData,file='readData.Rdata')
load('readData.Rdata')

#####
# 3  Apply MFA to the downloaded data sets using all cases or 
#just a subset of them in case there are too many. Selection of cases, 
#if applied, must be based on specific characteristics related to clinical 
#variables (for instance pathology state, gender, etc.). Keep in mind that you will 
#require several steps for preparing data by matching selected cases in the databases, 
#selecting necessary variables, transposing, etc. 
#All these steps need to be justified in the script and/or final report.
####

read_clinic<-getData(readData,"Clinical")
read_rna<-getData(readData,'RNASeq2GeneNorm')
read_met<-getData(readData,'Methylation',platform = 1)
read_mrna<-getData(readData,'mRNAArray',platform=1)


library(FactoMineR)
dim(read_rna)
# Main technologies/data downloaded: RNAseq2_Gene_Norm, Methylation, mRNA
rna_samples <- strsplit(colnames(read_rna), "-")
met_samples <- strsplit(colnames(read_met), "-")
mrna_samples <- strsplit(colnames(read_mrna), "-")
clinic_samples <- strsplit(rownames(read_clinic), "\\.")
rna_samples_2 <- unlist(lapply(rna_samples, function(x) x[3]))
met_samples_2 <- unlist(lapply(met_samples, function(x) x[3]))
mrna_samples_2 <- unlist(lapply(mrna_samples, function(x) x[3]))
clinic_samples_2 <- unlist(lapply(clinic_samples, function(x) toupper(x[3])))
# Get samples with information in all the datasets:
common1<-rna_samples_2[rna_samples_2 %in% met_samples_2]
common2<-common1[common1%in% mrna_samples_2]
common_samples <- common2[common2 %in% clinic_samples_2]
common_samples 

# Organize data and check NAs:
read_rna_3<- read_rna[,match(common_samples, rna_samples_2)]
colnames(read_rna_3) <- common_samples
any(is.na(read_rna_3))
read_met_3 <- read_met[,match(common_samples, met_samples_2)]
colnames(read_met_3) <- common_samples
any(is.na(read_met_3))
read_met_3 <- read_met_3[!apply(read_met_3,1,function(x) any(is.na(x))),]
read_mrna_3 <- read_mrna[,match(common_samples, mrna_samples_2)]
colnames(read_mrna_3) <- common_samples
any(is.na(read_mrna_3))
read_mrna_3 <- read_mrna_3[!apply(read_mrna_3,1,function(x) any(is.na(x))),]
read_clinic_2 <- read_clinic[match(common_samples, clinic_samples_2),]
rownames(read_clinic_2) <- common_samples
any(is.na(read_clinic_2))
read_clinic_2 <- read_clinic_2[ , !apply(read_clinic_2,2,function(x) any(is.na(x)))]

# Methylation:
# source("https://bioconductor.org/biocLite.R")
# biocLite("IlluminaHumanMethylation27k.db")
library(IlluminaHumanMethylation27k.db)
CpG_annotation <- as.list(IlluminaHumanMethylation27kCHR[mappedkeys(IlluminaHumanMethylation27kCHR)])

met_subs <- names(CpG_annotation[CpG_annotation %in% c("13","14","17","18")])
read_met_3 <- data.frame(read_met_3[match(met_subs,rownames(read_met_3)),]) # Subsetting samples
read_met_3 <- read_met_3[!apply(read_met_3,1,function(x) any(is.na(x))),]
any(is.na(read_met_3))

dim(read_met_3)
# RNA_Seq:
library(biomaRt)

hs_mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org") # Get the database.
grep("sapiens",listDatasets(hs_mart)[,1], value = T) # Get the dataset.
hs_mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl', host="www.ensembl.org")
head(listAttributes(hs_mart)) # Get the attributes
# Get the parameters of interest:

RNA_Seq_annotation <- getBM(c("external_gene_name", "chromosome_name"),
                            filters = "hgnc_symbol",
                            values = rownames(read_rna_3), mart = hs_mart)

rna_subs <- RNA_Seq_annotation$external_gene_name [RNA_Seq_annotation$chromosome_name %in% c("13","14","17","18")]
read_rna_3 <- data.frame(read_rna_3[match(rna_subs,rownames(read_rna_3)),]) # Subsetting samples

# mRNA:
mRNA_annotation <- getBM(c("external_gene_name", "chromosome_name"),
                         filters = "hgnc_symbol",
                         values = rownames(read_mrna_3)[grep("^[A-Za-z0-9]+$", # Remove some characters in order to solve query error
                                                        rownames(read_mrna_3),
                                                        perl=TRUE)],
                         mart = hs_mart)

mrna_subs <- mRNA_annotation$external_gene_name[mRNA_annotation$chromosome_name %in% c("13","14","17","18")]
length(mrna_subs)
read_mrna_3 <- data.frame(read_mrna_3[match(mrna_subs,rownames(read_mrna_3)),]) # Subsetting samples

# Merge data:
dim(read_rna_3)
dim(read_mrna_3)
dim(read_met_3)

#Eliminating those data with mean==0 due to problems with MFA function
read_mrna_3=t(read_mrna_3)
any(which(apply(read_mrna_3,2,mean) == 0))
read_rna_3= t(read_rna_3)
any(which(apply(read_rna_3,2,mean) == 0))

read_rna_3<-read_rna_3[,-which(apply(read_rna_3,2,mean) == 0)]
read_met_3= t(read_met_3)
read_met_3 <- data.frame(apply(read_met_3,2, as.numeric))
any(which(apply(read_met_3,2,mean) == 0))

DataBinding <- data.frame(read_rna_3, read_mrna_3, read_met_3)
head(DataBinding)
dim(DataBinding)
class(DataBinding[1,5000])
# Apply MFA:
pdf()
MFA_data <- MFA(DataBinding, group=c(2041,1779,2877), type = rep("s",3),
                ncp=5, name.group=c("RNA","mRNA","MET"),
                num.group.sup=NULL) # supplementary group not analyzed in the PCA in order not to have influence.

dev.off()
MFA_data


###
##4  Commentate on the MFA results
#####
# In the last graphic, groups representation, our three kind of data seem to be really close. 
# So they do not look like a random distribution. 
# The first dimension of the global PCA explains almost 22% of the initial variability
# and the second dimension explains the 11%. So in total we have 33% of variability explained by the
# two first dimensions. 
# In the first graphic, partial axes, we see five of the dimensions of the 
# partial PCAs. In general the dimentions hf the data are dispersed inside the pca plot,
# but most of them correlate directly or inversely with with the dimension 1 and 2 of the global pca.
# The graphic of the correlation circle is completely illegible so it makes no sense to
# commentate due to presence of many variables and many dimensions.
# And finally, in the remaining two plots, representing individual factor map,
# the data are really mixed, which is very good for the analysis, because the data is not dispersed between the dimensions,
# being most of them explained by dimension 1.

######
## 5 Generate a CIRCOS plot showing the results of the analysis as well as the
# downloaded data, justifying data filtering and election criteria for each track of the plot.
######

library(OmicCircos)

# Downloaded RNA_Seq:
grep("chr",listAttributes(hs_mart)[,1], value = T)[1]
grep("gene",listAttributes(hs_mart)[,1], value = T)
grep("pos",listAttributes(hs_mart)[,1], value = T)
grep("hgnc",listFilters(hs_mart)[,1], value = T)


RNA_Seq_annotation_2 <- getBM(c("external_gene_name", "chromosome_name", "start_position"),
                              filters = "hgnc_symbol",
                              values = rownames(read_rna), mart = hs_mart)
RNA_Seq_annotation_2 <- RNA_Seq_annotation_2[RNA_Seq_annotation_2$chromosome_name<3,]
RNA_Seq_annotation_2 <- RNA_Seq_annotation_2[!duplicated(RNA_Seq_annotation_2[,"external_gene_name"]),]
read_rna_4 <- data.frame(read_rna[match(RNA_Seq_annotation_2$external_gene_name,rownames(read_rna)),])
read_rna_4 <- cbind(chr=RNA_Seq_annotation_2$chromosome_name, po=RNA_Seq_annotation_2$start_position, NAME=RNA_Seq_annotation_2$external_gene_name,read_rna_4)
head(read_rna_4[,1:4])

# Downloaded  mRNA:
mRNA_annotation_2 <- getBM(c("external_gene_name", "chromosome_name", "start_position"),
                           filters = "hgnc_symbol",
                           values = rownames(read_mrna)[grep("^[A-Za-z0-9]+$", # Remove some characters
                                                            rownames(read_mrna),
                                                            perl=TRUE)],
                           mart = hs_mart)
head(mRNA_annotation_2)
mRNA_annotation_2 <- mRNA_annotation_2[mRNA_annotation_2$chromosome_name<3,]
mRNA_annotation_3 <- mRNA_annotation_2[!duplicated(mRNA_annotation_2$external_gene_name),]
read_mrna_4 <- merge(mRNA_annotation_3,read_mrna,by.x = 'external_gene_name',by.y='row.names')
rownames(read_mrna_4)<-read_mrna_4$external_gene_namecolnames(read_met_4)[1] <- "NAME"
colnames(read_mrna_4)[1] <- "NAME"
colnames(read_mrna_4)[2] <- "chr"
colnames(read_mrna_4)[3] <- "po"
read_mrna_4 <- read_mrna_4[,c(2,3,1,4:dim(read_mrna_4)[2])]
head(read_mrna_4[,1:4])

class(read_mrna_4)
# Downloaded Met:
any(is.na(read_met))
read_met_4 <- read_met[!apply(read_met,1,function(x) any(is.na(x))),]
colnames(read_met_4)[1] <- "NAME"
colnames(read_met_4)[2] <- "chr"
colnames(read_met_4)[3] <- "po"
read_met_4 <- read_met_4[,c(2,3,1,4:dim(read_met_4)[2])]
head(read_met_4[,1:4])


#MFA RESULTS

# Result dim 1, explaining more variance:
dim1 <- MFA_data$global.pca$var$cor[,"Dim.1"]
identical(names(dim1),names(DataBinding))
read_rna_5 <- dim1[1:length(colnames(read_rna_3))]
names(read_rna_5)<- gsub('\\.','-',names(read_rna_5))
identical(names(read_rna_5),colnames(read_rna_3)) 
RNA_Seq_annotation_3 <- getBM(c("external_gene_name", "chromosome_name", "start_position"),
                              filters = "hgnc_symbol",
                              values = names(read_rna_5), mart = hs_mart)
head(RNA_Seq_annotation_3)
RNA_Seq_annotation_3 <- RNA_Seq_annotation_3[RNA_Seq_annotation_3$chromosome_name<3,]
RNA_Seq_annotation_3 <- RNA_Seq_annotation_3[!duplicated(RNA_Seq_annotation_3[,"external_gene_name"]),]
read_rna_5 <- data.frame(read_rna_5[match(RNA_Seq_annotation_3$external_gene_name,names(read_rna_5))])
colnames(read_rna_5) <- "RNA_MFA"
read_rna_5 <- cbind(chr=RNA_Seq_annotation_3$chromosome_name, po=RNA_Seq_annotation_3$start_position, NAME=RNA_Seq_annotation_3$external_gene_name,read_rna_5)
head(read_rna_5[,1:4])


#mRNA

read_mrna_5 <- dim1[as.integer(length(colnames(read_rna_3))+1):as.integer(length(colnames(read_rna_3))+length(colnames(read_mrna_3)))]


names(read_mrna_5) <- gsub('\\..','',names(read_mrna_5))
identical(names(read_mrna_5),colnames(read_mrna_3))
mRNA_annotation_3 <- getBM(c("external_gene_name", "chromosome_name", "start_position"),
                           filters = "hgnc_symbol",
                           values = names(read_mrna_5),
                           mart = hs_mart)
head(mRNA_annotation_3)
mRNA_annotation_3 <- mRNA_annotation_3[mRNA_annotation_3$chromosome_name<3,]
mRNA_annotation_3 <- mRNA_annotation_3[!duplicated(mRNA_annotation_3[,"external_gene_name"]),]
read_mrna_5 <- data.frame(read_mrna_5[match(mRNA_annotation_3$external_gene_name,names(read_mrna_5))]); colnames(read_mrna_5) <- "mRNA_MFA"
read_mrna_5 <- cbind(chr=mRNA_annotation_3$chromosome_name, po=mRNA_annotation_3$start_position, NAME=mRNA_annotation_3$external_gene_name,read_mrna_5)
head(read_mrna_5[,1:4])


##MET

read_met_5 <-  dim1[as.integer(length(colnames(read_rna_3))+length(colnames(read_mrna_3))+1):length(dim1)]
identical(names(read_met_5),colnames(read_met_3))
read_met_idx <- read_met[match(names(read_met_5), rownames(read_met)),c(2,3,1)]
read_met_ann <- read_met_idx[complete.cases(read_met_idx),]
read_met_5 <- read_met_5[complete.cases(read_met_idx)]
read_met_5 <- cbind(read_met_ann,read_met_5); colnames(read_met_5) <- c("chr","po","NAME","MET_MFA")
head(read_met_5[,1:4])


# Plot Circos:
#Initial plot
pdf('circos.pdf')

colors <- rainbow(10, alpha=0.5);
par(mar=c(2, 2, 2, 2));
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="Circos Plot");

#No MFA data mapping
circos(R=325, cir="hg19", W=80, mapping=read_rna_4, col.v=4, type="s", B=F, lwd=1, col=colors[1])
legend(-20,855, c("#1: RNA","#2: mRNA", "#3: Meth", "#4: RNA_res","#5: mRNA_res","#6: Meth_res"), cex=.7,bty = "n")
circos(R=250, cir="hg19", W=80, mapping=read_mrna_4,  col.v=4,  type="heatmap2", 
       cluster=TRUE, col.bar=TRUE, lwd=1, col="blue",col.bar.po="topright")
legend(715,850,"Heatmap #2",cex=.8,text.col="black", box.col="white",bg="white") 
legend(-50,135,"Dendogram #2",cex=.8,text.col="black", box.col="white",bg="white") 
circos (R=170, cir="hg19" , W=80, mapping=read_met_4,col.v=4, type="lh" , B=F, col=colors[ 7 ] , lwd=2, scale=TRUE) ;
circos(R=400, type="chr", cir="hg18", col=colors, print.chr.lab=TRUE, W=4, scale=TRUE); # First circle is the chromosom

#MFA data mapping
circos(R=125, cir="hg19", W=40, mapping=read_rna_5, col.v=4, type="l" , B=TRUE, lwd=1, col=colors[1])
circos(R=85, cir="hg19", W=40, mapping=read_rna_5, col.v=4, type="l" , B=TRUE, lwd=1, col=colors[1])
circos(R=45, cir="hg19", W=40, mapping=read_met_5, col.v=4, type="l" , B=TRUE, lwd=1, col=colors[1])
dev.off()

