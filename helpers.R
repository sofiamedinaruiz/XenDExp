##################################################################################
# LOADING DATA
##################################################################################
library(shiny)
library(RColorBrewer,quietly=T)
library(png,quietly=T)
library(random,quietly=T)
library(reshape2,quietly=T)
library(ggplot2,quietly=T)
library(limma,quietly=T)
library(edgeR,quietly=T)
library(fastcluster, quietly=T)
library(Hmisc, quietly = T)
library(data.table, quietly = T)
library(networkD3, quietly = T)

#To install missing libraries run:
# source("https://bioconductor.org/biocLite.R")
# biocLite('<NAME OF LIBRARY>')
# bceoLite(c("RColorBrewer","png","random","reshape2","ggplot2","limma","edgeR","fastcluster","Hmisc","data.table","networkD3"))

###>>>>>  setwd('/Users/sofiamr/Dropbox/EctoMAP\ figures/DOTPLOTS\ AHMB/EctoMAP\ copy')
##save(data,file = "ALLdata.rda')
##load(file = "/Users/sofiamr/Dropbox/EctoMAP\ figures/DOTPLOTS\ AHMB/EctoMAP\ copy/ALLdata.rda")  ##added new


############ LOADING FILES ############
load(file = "load/ALLdata.rda")

#### #### #### #### NMF #### #### #### #### #### 

#read NMF matrix
nmf.ematrix.st12<-data.frame(fread("data/5.3_NMF_allgenes_st12_expression.txt",header=T,sep="\t") ,row.names=1)
nmf.ematrix.st14<-data.frame(fread("data/5.3_NMF_allgenes_st14_expression.txt",header=T,sep="\t") ,row.names=1)
#gene Names
genenames.st12<-rownames(nmf.ematrix.st12)
genenames.st14<-rownames(nmf.ematrix.st14)
#load drawings & rotate them 90deg clockwise
drawing.st12<-readPNG("data/NMF_drawing_st12_lr.png")
drawing.st12<-t(drawing.st12[nrow(drawing.st12):1,])
drawing.st14<-readPNG("data/NMF_drawing_st14_lr.png")
drawing.st14<-t(drawing.st14[nrow(drawing.st14):1,])
#NMF tissue color order
drawing.st12.order<-c("AE","ANP","PNP","NB","NNE")
drawing.st14.order<-c("AE","NP","ALNB","PLNB","NNE")


#### #### #### #### Gene Descriptions ----> Xenbase.org #### #### #### ####  

#geneData from Xenbase
gene_ManualCuration<-as.matrix(read.table("data/GenePageGeneralInfo_ManuallyCurated.txt",header=F,sep="\t",fill = TRUE))
#gene_ManualCuration<-fread("data/GenePageGeneralInfo_ManuallyCurated.txt", header=F, sep="\t", na.strings=NULL)  #sep="\t",fill = TRUE))

#In situ databse from Xenbase: ftp://ftp.xenbase.org/pub/GenePageReports/XenbaseGeneHumanOrthologMapping.txt
InSitu<-fread("data/InSituHybridization_laevis.txt", header=F,sep="\t")
colnames(InSitu)<-c("GeneID","Gene","NA","Stage","Annot_stage","ImgID","Notes")

### WGCNA ####
#Results from WGCNA
Adj.values <-as.data.frame(fread("WGCNA/Top50_order_pwr22value.txt", header=F, sep="auto", stringsAsFactors=FALSE))
Adj.genes <- as.data.frame(fread("WGCNA/Top50_order_pwr22names.txt", header=F, sep="auto", stringsAsFactors=FALSE))
WGCNAgene_module <-as.data.frame(fread("WGCNA/WGCNAgroups_names_dynModhyb_SP22-bicor-ds4-signed.csv", header=T, sep="\t", stringsAsFactors=FALSE))


### Other info ###
#Groups of genes
NC_genes<-c("snai2","sox10","tfap2b","foxd3","twist") 
NB_genes<-c("pax3","zic1","zic3","tfap2c","msx2")
Ventral_genes<-c("bmp4","bmp7","ventx","gata2","dlx3")
Anterior_genes<-c("otx2","hesx1","six3","rax","pitx1")
NP_genes<-c("foxa4","foxa2","sox2")
Epidermis<-c("foxi1","hmx3")
Posterior<-c("cdx","hoxa1")