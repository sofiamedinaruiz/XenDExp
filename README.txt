#######################    EctoMAP - v1.2   -    README.txt   #######################
This folder had the required documents to run the Shiny/R app to run locally.
The application is intended obtain the ectodermal expression pattern prediction 
for any gene present in Xenopus laevis during early, mid and late neurulation.

#################################    Updates  #################################
Last updated:  April 26th, 2017, by SMR
Correspondence: anne-helene.monsoro-burq@curie.fr; harland@berkeley.edu
Code available: https://github.com/sofiamedinaruiz/EctoMAP/

#########################        Requirements      #################################
Local installation of Rstudio v1.0.143 and R v.3.3.3 
R packages required to be installed: shiny, RColorBrewer, png, random, reshape2, ggplot2, limma, edgeR, fastcluster, Hmisc, data.table, networkD3

#################  Instructions to run EctoMAP in a Desktop computer  ################
1. Download Rstudio
Rstudio versionv 1.0.143 (https://www.rstudio.com/products/rstudio/download/)

2. Open ‘ui.R' with RStudio

3. To install all the required packages type in the R console:
source("https://bioconductor.org/biocLite.R")
biocLite('<NAME OF LIBRARY>')
###  If running using the App for first time install all required packages:
biocLite("shiny", " RColorBrewer", " png", " random", " reshape2", " ggplot2", " limma", " edgeR", " fastcluster", " Hmisc", " data.table", " networkD3")

4. Three ways to run the app in a web browser:
a) The easiest way to run the app: Press the "Run App" button on the top right corner of the prompt.
b) Go to the working directory by typing the address to the Shiny app:
> setwd ('~/path/to/shiny/EctoMAP')
> library(shiny)
> runApp()
c) Type
> runApp('Directory/Shinysofia2')

5. A web browser will show the app, wait 3 min until all the data has been loaded.




########################  Files in folder   ##############################
README.md 				- this file
server.R				- Shiny file
ui.R					- Shiny file
helpers.R				- Functions and file calling for server.R and ui.R
ExampleOutputShiny.pdf  - Example of an output.

Folders:
data/	<- For the NMF and embryo visualizations
	5.3_NMF_allgenes_st12_expression.txt
	5.3_NMF_allgenes_st14_expression.txt
	GenePageGeneralInfo_ManuallyCurated.txt
	InSituHybridization_laevis.txt
	NMF_drawing_st12_lr.png
	NMF_drawing_st14_lr.png
	
load/	<- RNA-seq expression data
	B100B108_cuffmerged_gene2hugo_annotation.txt
	B100B108_cuffmerged_gene2hugo_annotation_manual.txt
	B100B108_cuffmerged_gene_length.txt
	EctoMap.samples.readcounts.txt
	EctoMap.samples.txt
	data_annotation_hugo_name.txt
	genenames12and14.txt
	

WGCNA/  <- WGCNA results. For each gene we stored the top 50 gene names with the highest adjacency scores.
	Top50_order_pwr22value.txt
	Top50_order_pwr22value.txt
	WGCNAgroups_names_dynModhyb_SP22-bicor-ds4-signed.csv
	
	
##############################  Contributions ##############################
# Tissue collection and sample preparations - JLP, AHMB, SMR
# RNA-seq data analysis and NMF predictions - JLP
# Shiny code implementation and dynamic plot generation - SMR and JP
# Network analysis - SMR
# 

##############################     Funding    ##############################
# NIH
# France
# HHMI
# Conacyt - UCMEXUS


##############################   Publication  ##############################
#  Plouhinec and Medina-Ruiz et al. “A Molecular Atlas of the Developing Ectoderm Defines Neural, Neural Crest, Placode and Non-Neural Progenitor Identity in Vertebrates”  Under Revision, 2017.


