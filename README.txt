##################    XenDExp - v1   -    README.txt  #######################
This folder had the required documents to run the Shiny/R app to run locally.

The application is intended obtain the ectodermal expression pattern prediction for any gene present in Xenopus laevis during early, mid and late neurulation.

###########################    Updates  #####################################
Began: Feb 18th, 2016.
Last updated:  April 5th, 2016
By: Sofia Medina Ruiz - sofiamr@gmail.com
Correspondence: AHMB

#########################   Requirements    #################################
Local installation of Rstudio v0.99.878 and R v.3.2.3.
R packages: shiny, RColorBrewer, png, random, reshape2, ggplot2, limma, edgeR, fastcluster, Hmisc. data.table.

#################  Instructions to run XenDExp  ##############################
1. Download Rstudio
Rstudio version 0.99.878 (https://www.rstudio.com/products/rstudio/download/)

2. Open 'server.R' or 'ui.R' in RStudio

3. Install all the required packages by typing:
> install.packages ('<name of package>')

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
	5.2_NMF_st12_expression.txt
	5.2_NMF_st14_expression.txt
	B100_B108_merged_gene_length.txt
	GenePageGeneralInfo_ManuallyCurated.txt
	InSituHybridization_laevis.txt
	NMF_drawing_st12_lr.png
	NMF_drawing_st14_lr.png
	from the new file not great/   - Ignore, the NMF is not correct.
	rawcounts.txt
load/	<- RNA-seq expression data
	B100B108_cuffmerged_gene2hugo_annotation.txt
	B100B108_cuffmerged_gene2hugo_annotation_manual.txt
	B100B108_cuffmerged_gene_length.txt
	B100B108_htseq_counts_xl91cuff.txt
	HTSeq_samples_2.txt

##########################  Contributions ###################################
# Tissue collection and sample preparations - JLP, AHMB, SMR
# RNA-seq data analysis and NMF predictions - JLP
# Shiny code implementation and dynamic plot generation - SMR and JP
# Data analysis - AHMB, RH, JL, JP, SMR

##########################     Funding    ###################################
# Conacyt - 
# NIH?
# France?

