##################    XenDExp - v1   -    README.txt  #######################
This folder had the required documents to run the Shiny/R app to run locally.

The application is intended obtain the ectodermal expression pattern prediction 
for any gene present in Xenopus laevis during early, mid and late neurulation.

###########################    Updates  #####################################
Began: Feb 18th, 2016.
Last updated:  April 14th, 2016
By: Sofia Medina Ruiz - sofiamr@gmail.com
Correspondence, info regarding raw files: anne-helene.monsoro-burq@curie.fr; harland@berkeley.edu
Code available: https://github.com/sofiamedinaruiz/XenDExp

#########################   Requirements    #################################
Local installation of Rstudio v0.99.878 and R v.3.2.3.
R packages: shiny, RColorBrewer, png, random, reshape2, ggplot2, limma, 
edgeR, fastcluster, Hmisc, data.table, networkD3

#################  Instructions to run XenDExp  ##############################
1. Download Rstudio
Rstudio version 0.99.878 (https://www.rstudio.com/products/rstudio/download/)

2. Open 'server.R' or 'ui.R' in RStudio

3. Install all the required packages by typing:
> install.packages ('<name of package>')
OR
>source("https://bioconductor.org/biocLite.R")
>biocLite('<NAME OF LIBRARY>')

4. Three ways to run the app in a web browser:
a) Go to the working directory by typing the address to the Shiny app:
> setwd ('~/path/to/shiny/XenDExp')
> library(shiny)
> runApp()
b) Type
> runApp('Directory/Shinysofia2')
c) The easiest way to run the app: Press the "Run App" bottom on the top right corner of the prompt.

5. A web browser will show the app.


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
	B100B108_htseq_counts_xl91cuff.txt
	HTSeq_samples_2.txt
	data_annotation_hugo_name.txt
	genenames12and14.txt

WGCNA/  <- WGCNA results. For each gene we stored the top 50 gene names with the highest adjacency scores.
	Top50_order_pwr23names.txt
	Top50_order_pwr23value.txt
	WGCNA_net_color_Active-SP23-bicor-ds4-signed-signed.csv

##########################  Contributions ###################################
# Tissue collection and sample preparations - JLP, AHMB, SMR
# RNA-seq data analysis and NMF predictions - JLP
# Shiny code implementation and dynamic plot generation - SMR and JP
# Network analysis - SMR
# Data analysis - AHMB, RH, ME, JL, JP, SMR

##########################     Funding    ###################################
# Conacyt - UCMEXUS
# NIH?
# France?
# HHMI?
