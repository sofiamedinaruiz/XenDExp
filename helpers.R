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

###I problems installing libraries type:
#source("https://bioconductor.org/biocLite.R")
##biocLite('<NAME OF LIBRARY>')

#Load raw data and data tables with about the genes/transcripts.
load.rnaseq.data<-function(count.file,sample.file,gene.length.file,gene.hugo.file,manual.hugo.file,evalue.min=100)
{
  #data
  count.table<-read.table(count.file,header=T,sep="\t")
  sample.table<-read.table(sample.file,header=T,stringsAsFactors=F,sep="\t")
  gene.length.table<-read.table(gene.length.file,header=F,stringsAsFactors=F,sep="\t")
  gene.hugo.table<-read.table(gene.hugo.file,header=T,stringsAsFactors=F,sep="\t")
  
  #remove non-gene lines in count table
  rawcounts<-count.table[!(rownames(count.table) %in% c("no_feature","ambiguous","not_aligned","alignment_not_unique","too_low_aqual")),]
  cat("Count table has",nrow(rawcounts),"genes and",ncol(rawcounts),"samples.\n")
  #cat(paste(colnames(rawcounts),sep=" "),"\n")
  
  #extract sample annotation and compute design (use sample order from count table)
  rownames(sample.table)<-sample.table[,1]
  sample.names<-colnames(rawcounts)
  cat("Sample annotation table has",nrow(sample.table),"samples.")
  cat(sum(sample.names==sample.table[sample.names,1]),"names match count table samples.\n")
  stage<-sample.table[sample.names,"Stage"]
  tissue<-sample.table[sample.names,"Tissue"]
  genotype<-sample.table[sample.names,"Genotype"]
  id<-sample.table[sample.names,"ID"]
  sampletype<-paste(stage,tissue,genotype,sep="_")
  design<-model.matrix(~sampletype)
  
  #extract annotations
  #read manual annotation table if it exists and update genehugotable
  #take 4 columns in genehugotable table (name,blast,hugo,blast homology score)
  annotation.table<-gene.hugo.table[,c(1,2,3,7)]
  rownames(annotation.table)<-annotation.table[,1]
  cat("Annotation table has",nrow(annotation.table),"entries.\n")
  #replace with manual annotation if it exists
  if (!is.null(manual.hugo.file))
  {
    manual.hugo<-read.table(manual.hugo.file,header=T,stringsAsFactors=F,sep="\t")[,c(1,2,3,7)]
    rownames(manual.hugo)<-manual.hugo[,1]
    cat("There are",nrow(manual.hugo),"manual annotations.\n")
    selected.manual.names<-rownames(manual.hugo)[(rownames(manual.hugo) %in% rownames(annotation.table))]
    cat(length(selected.manual.names),"annotations to be replaced.\n")
    annotation.table[selected.manual.names,]<-manual.hugo[selected.manual.names,]
  }
  #replace none in blast and hugo with ""
  no.blast.annotation<-(annotation.table[,2]=="none")
  annotation.table[no.blast.annotation,2]<-""
  annotation.table[(annotation.table[,3]=="none"),3]<-""
  #add new column with name and blast annotation in bracket if available (and above e-value threshold:no, just if available)
  name.and.blast<-annotation.table[,1]
  #blast.above.threshold<-(annotation.table[,4]>evalue.min)
  selected.genes.blast<-!no.blast.annotation #(blast.above.threshold&(!no.blast.annotation))
  name.and.blast[selected.genes.blast]<-paste(annotation.table[selected.genes.blast,1],"[",
                                              annotation.table[selected.genes.blast,2],"]",sep="")
  annotation.table<-cbind(annotation.table,name.and.blast)
  
  #add rownames to tables, sort them according to same gene order, and add gene length to annotation
  rownames(annotation.table)<-annotation.table[,1]
  rownames(gene.length.table)<-gene.length.table[,1]
  #gene names in annotation and not in counts: they seem to represent genes that have been merged with other genes during the pipeline merge/reduction step
  annotation.not.in.counts<-annotation.table[!(rownames(annotation.table) %in% rownames(rawcounts)),]
  counts.not.in.annotation<-rownames(rawcounts)[!(rownames(rawcounts) %in% rownames(annotation.table))]
  cat(nrow(annotation.not.in.counts),"genes in annotation not in counts.",
      length(counts.not.in.annotation),"genes in counts not in annotation.\n")
  #restrict and order annotation and gene length according to rawcounts
  annotation.table<-annotation.table[rownames(rawcounts),]
  gene.length.table<-gene.length.table[rownames(rawcounts),]
  annotation.table<-cbind(annotation.table,gene.length.table[,2])
  colnames(annotation.table)[1:6]<-c("gene.name","blast.annotation","hugo.name","blast.score","name.and.blast.annotation","gene.length")
  cat(sum(annotation.table[,2]!="")," genes have blast annotation.\n")
  cat(sum(annotation.table[,3]!="")," genes have hugo annotation.\n")
  
  #change sample names to make them more explicit
  colnames(rawcounts)<-paste(sample.table[colnames(rawcounts),1],
                             sample.table[colnames(rawcounts),"Stage"],
                             sample.table[colnames(rawcounts),"Tissue"],
                             sample.table[colnames(rawcounts),"Genotype"],
                             sep="-")
  
  #order rawcounts and annotation table the same way, replace gene name in count table and annotation table with name+blastannot
  #sum(rownames(annotation.table)!=rownames(rawcounts))
  #rename annotation&rawcounts according to name+blast annot
  names.and.annot<-annotation.table[,5]
  rownames(annotation.table)<-names.and.annot
  rownames(rawcounts)<-names.and.annot
  
  #return as list
  return(list(rawcounts=rawcounts,sampletype=sampletype,design=design,
              stage=stage,tissue=tissue,genotype=genotype,id=id,
              annotations=annotation.table))
}


Adj.genes <- as.data.frame(fread("WGCNA/Top50_order_pwr23names.txt", header=F, sep="auto", stringsAsFactors=FALSE))
Adj.values <-as.data.frame(fread("WGCNA/Top50_order_pwr23value.txt", header=F, sep="auto", stringsAsFactors=FALSE))
WGCNAclusters <-as.data.frame(fread("WGCNA/WGCNA_net_color_Active-SP23-bicor-ds4-signed-signed.csv", header=T, sep="auto", stringsAsFactors=FALSE))

data<-load.rnaseq.data("load/B100B108_htseq_counts_xl91cuff.txt",
                       "load/HTSeq_samples_2.txt",
                       "load/B100B108_cuffmerged_gene_length.txt",
                       "load/B100B108_cuffmerged_gene2hugo_annotation.txt",
                       "load/B100B108_cuffmerged_gene2hugo_annotation_manual.txt")


############ LOADING FILES ############


#read NMF matrix
nmf.ematrix.st12<-data.frame(fread("data/5.3_NMF_allgenes_st12_expression.txt",header=T,sep="\t") ,row.names=1)
nmf.ematrix.st14<-data.frame(fread("data/5.3_NMF_allgenes_st14_expression.txt",header=T,sep="\t") ,row.names=1)

genenames.st12<-rownames(nmf.ematrix.st12)
genenames.st14<-rownames(nmf.ematrix.st14)

#geneData
gene_ManualCuration<-as.matrix(read.table("data/GenePageGeneralInfo_ManuallyCurated.txt",header=F,sep="\t",fill = TRUE))
#gene_ManualCuration<-fread("data/GenePageGeneralInfo_ManuallyCurated.txt", header=F, sep="\t", na.strings=NULL)  #sep="\t",fill = TRUE))

InSitu<-fread("data/InSituHybridization_laevis.txt", header=F,sep="\t")
colnames(InSitu)<-c("GeneID","Gene","NA","Stage","Annot_stage","ImgID","Notes")
### ftp://ftp.xenbase.org/pub/GenePageReports/XenbaseGeneHumanOrthologMapping.txt


#load drawings
drawing.st12<-readPNG("data/NMF_drawing_st12_lr.png")[,,1]
drawing.st14<-readPNG("data/NMF_drawing_st14_lr.png")[,,1]
#rotate them 90deg clockwise
drawing.st12<-t(drawing.st12[nrow(drawing.st12):1,])
drawing.st14<-t(drawing.st14[nrow(drawing.st14):1,])

colnames(WGCNAclusters)<-c('index','colors')
WGCNAclusters$index<-NULL




#NMF tissue color order
drawing.st12.order<-c("AE","ANP","PNP","NB","NNE")
drawing.st14.order<-c("AE","NP","ALNB","PLNB","NNE")

 


###  GENE GROUPS ###
NC_genes<-c("snai2","sox10","tfap2e","tfap2b","foxd3","twist") 
NB_genes<-c("pax3","zic1","zic3","tfap2a","tfap2c")
Ventral_genes<-c("bmp4","bmp7","ventx","gata2","dlx3")
Anterior_genes<-c("otx2","hesx1")
NP_genes<-c("foxa4","foxa2","sox2")
Epidermis<-c("gbx2","hmx3")




