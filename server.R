library(shiny)
source("helpers.R")
##One day I would like to represent the network this way: http://bl.ocks.org/d3noob/8043434

###### Documentation Ectomap v.1.0
# README.txt for info on how to run the EctoMAP app. 
# Respository: https://github.com/sofiamedinaruiz/EctoMAP/
# To run the app:
# You require to install the libraries found in helper.R
# Updated: Dec 15th, 2016
# To run:
#      runApp('The-Directory-for/EctoMAP')
#setwd('/Users/sofiamedina/Documents/Harland\ Lab/GeneExpression2016/EctoMAP\ copy')

load(file = "load/ALLdata.rda")  ##added new
WGCNAgene_module <-as.data.frame(fread("WGCNA/WGCNAgroups_names_dynModhyb_SP22-bicor-ds4-signed.csv", header=T, sep="\t", stringsAsFactors=FALSE))

shinyServer(function(input, output) {
  
  output$xenPlot1.ui <- renderUI({ 
    if (input$Show_expression_pattern_coexp ==TRUE){
      plotOutput("xenPlot1", width= "100%", height="300px")
    }
  })
  
  output$xenPlot2.ui <- renderUI({ 
    if (input$Show_expression_pattern_coexp ==TRUE){
      plotOutput("xenPlot2", width= "100%", height="300px")
    }
  })
  
  ##It will provide a link to Xenbase for the gene of interest
  output$Gene_link1 <- renderUI({ 
    if (input$Show_expression_pattern_coexp ==TRUE){
      short_gene <- findhugo(input$gene)
      website <- get_website(short_gene)
      GENEID_description<-get_gene_description(short_gene)
      if(length(website)==0){  
        tags$div(class = "header", checked =NA,
               #  br(),
                 tags$p("The diagrams below display the relative ectodermal expression of a given transcript."),br(),
                 tags$h5(strong("1st transcript: "),em(input$gene)), 
                 tags$p(strong("Gene name:   "), em(short_gene)))
      }else{
        tags$div(#class = "header", checked = NA,
                # br(),
                 tags$p("The diagrams below display the relative ectodermal expression of a given transcript."),br(),
                 tags$h5(strong("1st transcript: "),em(input$gene)), 
                 tags$p(strong(em(short_gene)), "- ",GENEID_description),
                 tags$a(href = paste0(website), "Link to Xenbase", target="_blank") )
      }
    }
  })
  
  output$Gene_link2 <- renderUI({ 
    if (input$Show_expression_pattern_coexp ==TRUE){
      short_gene <- findhugo(input$gene2)
      website <- get_website(short_gene)
      GENEID_description<-get_gene_description(short_gene)
      if(length(website)==0){  
        tags$div(class = "header", checked =NA,
                tags$h5(strong("2nd transcript: "),em(input$gene2)),
                tags$p(strong("Gene name:   "), short_gene))
      }else{
        tags$div(class = "header", checked = NA,
                tags$h5(strong("2nd transcript: "),em(input$gene2)),
                tags$p(strong(em(short_gene)), "- ",GENEID_description),
                tags$a(href = paste0(website), "Link to Xenbase", target="_blank") )
      }
    }
  })
  
  output$Text_output_coexp.ui <- renderUI({
    if (input$Show_coexp == TRUE){
      tags$div(class = "header", checked =NA,
               tags$p("The following table displays the 50 transcripts that show close expression profile similarities with the first trascript from Section 2:", strong(em(input$gene)),". This list was obtained by retrieving the top-50 most adjacent transcripts from the unsigned adjacency matrix (2nd column). The Pearson's correlation and associated p-value (3rd and 4th column) result from the pairwise comparisons between two genes across all 79 spatio/temporal datasets. The last column displays the assigned module number, which can range from 0-140"),
               
               #tags$p("Input transcript (from section 2.a): ", em(input$gene)),
               br()
      )
    }
    if (input$Show_coexp == FALSE){
      tags$div(class = "header", checked =NA,
               tags$p("If desired, one can obtain a table that displays the 50 transcripts that show close expression profile similarities with the first trascript input on Section 2."),
               
               #tags$p("Input transcript (from section 2.a): ", em(input$gene)),
               br()
      )
    }
  })
  
  output$Adjacent_matrix1 <- renderDataTable({ 
    if (input$Show_coexp ==TRUE){
      gene1 <- input$gene
      return(obtain_adjacent_and_cor_values(gene1))
    }
  }, options = list(lengthMenu = c(5, 10, 20, 30, 40, 50), pageLength = 10, searching = FALSE))
  
###If you want want to show the co-expression list for Transcript #2 deselect this section and adjust parameters as needed.
##  output$Text_output_coexp2.ui <- renderUI({
##    if (input$Show_coexp == TRUE){
##      tags$div(class = "header", checked =NA,
##               tags$h4(strong("1.b."), "Genes that share spatiotemporal expression with respect to ", input$gene2),
##               tags$p("Adjacency matrix was obtained using the WGCNA package (soft power = 23):")
##      )
##    }
##  })
##  
##  output$Adjacent_matrix2 <- renderDataTable({ 
##    if (input$Show_coexp == TRUE){
##      gene2 <- input$gene2
##      return(obtain_adjacent_values(gene2))
##    }
##  }, options = list(lengthMenu = c(5, 10, 20, 30, 40, 50), pageLength = 5, searching = FALSE))
  
  output$timeSeries.ui <-renderUI({
    #Show_time_series=TRUE
    #if (input$Show_time_series == TRUE){
      selected_gene<-c(input$genes_for_timeserie)
      plotOutput("timeSeries", width= "90%", height= "450px")
    #}
  })
  
  output$tree.ui <-renderUI({
    if(input$Draw_tree_option == TRUE){
      plotOutput("tree", width="95%", height= input$px4/2)
    }
  })
  
  output$Draw_Network.ui <-renderUI({
    if(input$Draw_Network == TRUE){
      forceNetworkOutput("Network", height="1500", width="95%")
    }
  })
 
  output$InSitus12 <- renderUI({ 
    if ((input$In_situ_option == TRUE) & (input$Show_expression_pattern_coexp ==TRUE)){
      short_gene <- findhugo(input$gene)
      insitu_website<-get_insitu_links(short_gene,12)
      if(length(insitu_website)==0){
        tags$div(class = "header", checked =NA,
                 tags$p("      ** No in situ image available from St. 11-13"),
                 tags$p(""))
      }else{  
        tags$div(
          class = "header", checked = NA,
          tags$a(href = paste0(insitu_website), "Mid neurula (St. 11-13) from Xenbase ",target="_blank"),
          tags$img(src=insitu_website,  width ="600px") )
      }
    }
  })
  
  output$InSitus14 <- renderUI({ 
    if ((input$In_situ_option == TRUE) & (input$Show_expression_pattern_coexp ==TRUE)){
      short_gene <- findhugo(input$gene)
      insitu_website<-get_insitu_links(short_gene,14)
      if(length(insitu_website)==0){
        tags$div(class = "header", checked =NA,
                 tags$p("      ** No in situ image available from St. 14-17"),
                 tags$p(""))
      }else{  
        tags$div(
          class = "header", checked = NA,
          tags$a(href = paste0(insitu_website), "Mid neurula (St. 14-17) from Xenbase ",target="_blank"),
          tags$img(src=insitu_website,  width ="600px") )
      }
    }
  })
  
  # Plot the embryo colored by the expression level of a gene
  ###We want to be able to plot all different variants of the genes
  
  #First expression pattern
  output$xenPlot1 <- renderPlot({
    if (input$Show_expression_pattern_coexp ==TRUE){
      
      
      par(mfrow=c(1,2))
      
      if(!is.na(match(input$gene, genenames.st12)))
      { nmf.plot.gene(input$gene,drawing.st12,drawing.st12.order,nmf.ematrix.st12,title.name="Stage 12.5")
      } else { image(matrix(0,nrow=1,ncol=1),col="white",axes=F) }
      
      if(!is.na(match(input$gene,genenames.st14)))
      { nmf.plot.gene(input$gene,drawing.st14,drawing.st14.order,nmf.ematrix.st14,title.name="Stage 14")
      } else { image(matrix(0,nrow=1,ncol=1),col="white",axes=F) }
    }
  }) 
  
  #Second expression pattern
  output$xenPlot2 <- renderPlot({
    if (input$Show_expression_pattern_coexp ==TRUE){
      
      par(mfrow=c(1,2))
      if(!is.na(match(input$gene2, genenames.st12)))
      { nmf.plot.gene(input$gene2,drawing.st12,drawing.st12.order,nmf.ematrix.st12,title.name="Stage 12.5")
      } else { image(matrix(0,nrow=1,ncol=1),col="white",axes=F) }
      
      if(!is.na(match(input$gene2,genenames.st14)))
      { nmf.plot.gene(input$gene2,drawing.st14,drawing.st14.order,nmf.ematrix.st14,title.name="Stage 14")
      } else { image(matrix(0,nrow=1,ncol=1),col="white",axes=F) }
    }
  })

  #1 plot call
  output$timeSeries <-renderPlot({
    
    time.series.plot(c(input$genes_for_timeserie),input$Nomalization_gene,findhugo(input$Nomalization_gene))
  })
  
  output$Absolute2.ui<- renderUI({
    if (input$Expression_levels_option == TRUE) {
      if ((input$Use_prev2 == TRUE) | (length(input$gene_list_for_expression) > 1) | (input$Use_prev3 == TRUE)){
        plotOutput("Absolute2plot", width= "90%", height= input$px)
      }
    }
   })

  output$Absolute2plot <-renderPlot({
    selected_stage<- input$Select_stage11 
    selected_genotype<- "WT"
    selected_gene<-extend.gene.list_samples_from_hugo(c(input$gene_list_for_expression))
    if (input$Use_prev2 == TRUE){
      selected_gene<-append(selected_gene,c(input$genes_for_timeserie))
    }        
    if (input$Use_prev3 == TRUE){
      short_gene <- input$gene
      short_gene2 <- input$gene2
      selected_gene<-append(selected_gene,c(short_gene,short_gene2))
    }
    selected_gene<-unique(sort(c(selected_gene)))
    if  (length(selected_gene)>0){
      absolute2.plot(selected_gene,selected_stage,selected_genotype)
    }else{
      paste('Need to provide a gene name.')
    }
  })
  
  output$heatplot.ui <-renderUI({
    if(input$Make_correlation_option == TRUE){
      selected_gene<-c(input$gene_list_for_expression)
      if((length(input$gene_list_for_expression) > 0) | (input$Use_prev2 == TRUE ) | (input$Use_prev3 == TRUE )) {
        height= length(input$gene_list_for_expression)*50
        plotOutput("heatplot", width="100%", height= input$px3)
      }
    }else{
      paste ('not TRUE - something is not printing')
    }
  })
  
  output$heatplot <- renderPlot({
    if( input$Make_correlation_option == TRUE){
      selected_gene<- c()
      selected_gene<- c(input$gene_list_for_expression)
      Entry<-c()
        if (input$Landmark_genesNC){Entry<-append(Entry,NC_genes ) }
        if (input$Landmark_genesNB){ Entry<-append(Entry,NB_genes) }
        if (input$Landmark_genesVentral){ Entry<-append(Entry,Ventral_genes) }
        if (input$Landmark_genesAnterior){Entry<-append(Entry,Anterior_genes) }
        if (input$Landmark_genesNP){ Entry<-append(Entry,NP_genes) }
        if (input$Landmark_genesEpi){Entry<-append(Entry,Epidermis) }
        if (input$Landmark_genesPos){Entry<-append(Entry,Posterior) }
      if (input$Use_prev2 == TRUE){
        Entry<-append(Entry,reduce.gene.list_samples(c(input$genes_for_timeserie)))
      }        
      if (input$Use_prev3 == TRUE){
        short_gene <- findhugo(input$gene)
        short_gene2 <- findhugo(input$gene2)
        Entry<-append(Entry,c(short_gene,short_gene2))
      }
        selected_gene<- unique(append(Entry,c(input$gene_list_for_expression,reduce.gene.list_samples(c(input$genes_for_timeserie)))))
        if (length(selected_gene)>1){
            heat.plot(unique(selected_gene), input$corr_type)
        }
    }else{}
  })
  
  output$tree <- renderPlot({
    if(input$Draw_tree_option == TRUE){
      selected_gene<- c()
      selected_gene<- c(input$gene_list_for_expression)
      Entry<-c()
      if (input$Landmark_genesNC){Entry<-append(Entry,NC_genes ) }
      if (input$Landmark_genesNB){ Entry<-append(Entry,NB_genes) }
      if (input$Landmark_genesVentral){ Entry<-append(Entry,Ventral_genes) }
      if (input$Landmark_genesAnterior){Entry<-append(Entry,Anterior_genes) }
      if (input$Landmark_genesNP){ Entry<-append(Entry,NP_genes) }
      if (input$Landmark_genesEpi){Entry<-append(Entry,Epidermis) }
      if (input$Landmark_genesPos){Entry<-append(Entry,Posterior) }
      if (input$Use_prev2 == TRUE){
        Entry<-append(Entry,reduce.gene.list_samples(c(input$genes_for_timeserie)))
      }        
      if (input$Use_prev3 == TRUE){
        short_gene <- findhugo(input$gene)
        short_gene2 <- findhugo(input$gene2)
        Entry<-append(Entry,c(short_gene,short_gene2))
      }
      selected_gene<- unique(append(Entry,c(input$gene_list_for_expression,reduce.gene.list_samples(c(input$genes_for_timeserie)) )))
      if (length(selected_gene)>1){
          tree_plot(unique(selected_gene), input$corr_type)
        }
    }
  })

  output$Network <- renderForceNetwork({
    if(input$Draw_Network == TRUE){
      if((length(input$gene_list_for_expression) > 0 ) | (input$Use_prev2 == TRUE) | (input$Use_prev3 == TRUE)){
        selected_gene<-c()
        Entry<-c()
        if (input$Landmark_genesNC){Entry<-append(Entry,NC_genes ) }
        if (input$Landmark_genesNB){ Entry<-append(Entry,NB_genes) }
        if (input$Landmark_genesVentral){ Entry<-append(Entry,Ventral_genes) }
        if (input$Landmark_genesAnterior){Entry<-append(Entry,Anterior_genes) }
        if (input$Landmark_genesNP){ Entry<-append(Entry,NP_genes) }
        if (input$Landmark_genesEpi){Entry<-append(Entry,Epidermis) }
        if (input$Landmark_genesPos){Entry<-append(Entry,Posterior) }
        if (input$Use_prev2 == TRUE){
          Entry<-append(Entry,reduce.gene.list_samples(c(input$genes_for_timeserie)))
        }        
        if (input$Use_prev3 == TRUE){
          short_gene <- findhugo(input$gene)
          short_gene2 <- findhugo(input$gene2)
          Entry<-append(Entry,c(short_gene,short_gene2))
        }
        #selected_gene<- unique(append(Entry,c(input$gene_list_for_expression,input$genes_for_timeserie) ))
        selected_gene<- unique(append(Entry,c(input$gene_list_for_expression,reduce.gene.list_samples(c(input$genes_for_timeserie)))))
        Net<-make_Network(selected_gene, input$p_value_treshold)
        return(Net)
      }
    }  
  })

})


############################################################################## 
########################## Functions for plotting ############################ 
############################################################################## 


##### MAKE SURE TO READ ALL SAMPLESALL
#### DO not filter and normalize again! ... or sould I?






#Dotplot
absolute2.plot <-function(gene_name_str, stage, genotype)
{
  subgroup  <-(data$genotype==genotype)&(data$stage==stage)
  samples   <- filter.n.normalize1(data,subgroup)  #samples are in RPKM values
  gene_name <- gene_name_str
  data.meansONLY <- expr.mean.comp(samples)
  
  cpm_enrichment_values <- obtain.CMP_and_enrichment.values(data.meansONLY, gene_name)
  cpm_enrichment_values1<-cpm_enrichment_values[cpm_enrichment_values$Tissue != 'NB', ]
  cpm_enrichment_values<-cpm_enrichment_values1[cpm_enrichment_values1$Tissue != 'NP', ]
  
  #  qplot(x = log(RPKM,2), y = Tissue, data = cpm_enrichment_values, geom = "point", color=Ectodemal_enrichment, size=RPKM,  ylab = "Ectodermal tissue type", xlab=paste("Log2(Mean RPKM) - St.",stage,sep=" "),  main=toString(gene_name_str), xlim=c(min(log(cpm_enrichment_values$RPKM,2),0),max(log(cpm_enrichment_values$RPKM,2))*1.1)) + 
  
  qplot(x = RPKM, y = Tissue, data = cpm_enrichment_values, geom = "point", color=Ectodemal_enrichment, size=100, show.legend = NA, ylab = "Ectodermal tissue type", xlab=paste("Mean RPKM - St.",stage,sep=" "),  main=toString(reduce.gene.list_samples(unique(cpm_enrichment_values$Gene))), xlim=c(min(cpm_enrichment_values$RPKM,0),max(cpm_enrichment_values$RPKM)*1.1)) + 
    coord_flip()+
    theme_bw() +
    theme( plot.title=element_text(face="bold.italic", size =16, angle=0, hjust = 0.5 ),
           
           axis.title.y=element_text(face="bold", size=16, angle=90 , hjust = 0.5, vjust = 0.5),
           axis.title.x=element_text(face="bold", size=16,angle=0, hjust = 0.5, vjust = 0.5), 
           axis.text.x=element_text(size=12), 
           axis.text.y=element_text(size=12), 
           
           strip.text.x= element_text(angle = 0, size = 11, hjust = 0.5, vjust = 0.5, face="bold.italic"),
           strip.background = element_rect(colour = 'steelblue', fill = 'white', size = 0.2),
           
           legend.title=element_text(size=10), 
           legend.position="top", 
           legend.text =element_text(size=10),
           legend.key = element_rect(colour = NA),
           legend.key.width = unit(1,"cm"),
           legend.key.height = unit(.5,"cm")
           ####10/03 ,       aspect.ratio = 1.1
    ) + 
    geom_vline(xintercept = 1, color="gray", linetype = "longdash") +
    scale_color_gradient( low="darkorange", high="darkgreen", space="Lab", guide = "colorbar", limits = c(0, 100))  + 
    facet_wrap(~ Gene, ncol = 4) 
}


#plot function
nmf.plot.gene<-function(gene.name,embryo.image,draw.order,nmf.ematrix,insitu.color="blue",title.name="",cex.name=1)
{
  #create palette
  draw.palette<-colorRampPalette(c("white",insitu.color))(n=100)
  #colors from NMF values
  #draw.nmf.colors<-draw.palette[(1+floor(99*nmf.ematrix[gene.name,draw.order]/max(nmf.ematrix[gene.name,draw.order])))]
  draw.nmf.colors<-draw.palette[as.matrix(1+floor(nmf.ematrix[gene.name,draw.order]/max(nmf.ematrix[gene.name,draw.order])*99))[,1:5]]
  
  draw.values<-sort(unique(as.vector(embryo.image)))
  draw.values.l<-length(draw.values)
  #breakpoints from image
  draw.breaks<-c(draw.values[1],((draw.values[1:(draw.values.l-1)]+draw.values[2:draw.values.l])/2),draw.values[draw.values.l])
  #plot
  image(embryo.image,col=c("black",draw.nmf.colors,"white"),axes=F, asp=1,breaks=draw.breaks)
  title(main=paste(title.name," - ",gene.name),outer=F,cex.main=cex.name)
}
### The end ###


time.series.plot<-function(gene_name_str,Nomalization_gene,Nomalization_gene_hugo){
  subgroup_time_serieA <-(data$genotype=='WT')&(data$tissue=="WE")&(data$stage!="14")
  samples<-filter.n.normalize1(data,subgroup_time_serieA)
  e.matrix<-2**samples$voom$E
  
#  genes_to_plot<- extend.gene.list_samples(gene_name_str, samples)
  genes_to_plot<-gene_name_str    ### make sure it is present in the dataset ####TODO
 ####HERE#### genes_to_plot <- gene_name_str
  
  if(Nomalization_gene != "no normalization"){
    norm_gene <- Nomalization_gene
    #test1<- t(data.frame(apply(e.matrix[genes_to_plot,],1,function(x) ((x/e.matrix[norm_gene,])))))
    test1<- t(data.frame(apply(e.matrix[which(match(rownames(e.matrix),genes_to_plot)!= "NA"),],1,function(x) ((x/e.matrix[norm_gene,])))))
    
  }else{
    ##test1<- data.frame(e.matrix[genes_to_plot,])
    test1<- data.frame(e.matrix[which(match(rownames(e.matrix),genes_to_plot)!= "NA"),])
  }
  
  colnames(test1)<-samples$stage
  test<- melt(t(as.matrix(test1)))
  test1<-test
  colnames(test1)<-c("Stage","Gene","value")
  test1$Stage <- factor(test1$Stage,levels=c("11-11.5","12","12.5","13","14-15","16-17","18-19"))

  if(Nomalization_gene != "no normalization"){
    qplot(Stage, value, fill=factor(Gene), data=test1, geom=c("point"),  color=Gene, xlab = "NF Stage", ylab = paste("Expression relative to",Nomalization_gene_hugo),  main=toString(reduce.gene.list_samples(test1$Gene)), show.legend=FALSE) +  
      stat_summary(fun.y=mean, geom="path", aes(color=Gene, group=Gene),  alpha= 1) +
      theme_bw() +
      theme( plot.title=element_text(face="bold.italic", size=16, angle=0 ,  hjust = 0.5 ),
             axis.title.y=element_text(face="bold", size=16, angle=90 ),
             axis.title.x=element_text(size=16,angle=0), 
             axis.text.x=element_text(size=14), 
             axis.text.y=element_text(size=14), 
             legend.title=element_blank(), 
             legend.position="top",
             legend.text =element_text(face="italic", size=12),
             legend.key = element_rect(colour = NA) 
      )
  }else{
    qplot(Stage, value, fill=factor(Gene), data=test1, geom=c("point"),  color=Gene, xlab = "NF Stage", ylab = "RPKM",   main=toString(reduce.gene.list_samples(test1$Gene)), show.legend=FALSE) +  
      geom_hline(yintercept = 1, color = "gray", linetype = "longdash") +
      stat_summary(fun.y=mean, geom="path", aes(color=Gene, group=Gene),  alpha= 1) +
      theme_bw() +
      theme( plot.title=element_text(face="bold.italic", size=16, angle=0,  hjust = 0.5  ),
             axis.title.y=element_text(face="bold", size=16, angle=90 ),
             axis.title.x=element_text(size=16,angle=0), 
             axis.text.x=element_text(size=14), 
             axis.text.y=element_text(size=14), 
             legend.title=element_blank(), 
             legend.position="top",
             legend.text =element_text(face="italic", size=12),
             legend.key = element_rect(colour = NA) 
      )
  }
  
}

heat.plot<-function(gene_name_str, corr_type){
  subgroup_time_serie <-(data$genotype=='WT')
  samplesALL<-filter.n.normalize1(data,subgroup_time_serie)
  gene_name<-extend.gene.list_samples_from_hugo(sort(unique(gene_name_str)))
  e.matrix<-2**samplesALL$voom$E[gene_name,]
  gene.expression.matrix<-t(e.matrix[gene_name,])
  pvalue <- melt(rcorr(t(e.matrix[gene_name,]), type=corr_type)$P)
  correlation_data.m <- melt(rcorr(t(e.matrix[gene_name,]), type=corr_type)$r)
  
  correlation_data.m$pvalue<-pvalue$value
  correlation_data.m$stars[pvalue$value>=.1] <- ' '
  correlation_data.m$stars[pvalue$value<.1 & correlation_data.m$pvalue>.01] <- '.'
  correlation_data.m$stars[pvalue$value<.01 & correlation_data.m$pvalue>.001] <- '*'
  correlation_data.m$stars[pvalue$value<.001 & correlation_data.m$pvalue>.0001] <- '* *'
  correlation_data.m$stars[pvalue$value<.0001 & correlation_data.m$pvalue>.00001] <- '* * *'
  correlation_data.m$stars[pvalue$value<.00001] <- '* * * *'
  correlation_data.m$stars[pvalue$value==0] <- '0'
  #correlation_data.m
  
  qplot(x=Var1, y=Var2, data=correlation_data.m, fill=value, geom="tile") + 
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "gene.expression.matrix", 
                         name=paste(corr_type,'correlation', sep="\n")) +
    theme(  axis.text.y=element_text(size=14, face="italic"), 
            axis.text.x=element_text(size=10, face="italic", angle=30, hjust=1,vjust = 1),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()) +
    geom_text(aes(Var1,Var2, label = stars), color = 'black', size=3) + 
    coord_equal(ratio = 1.3)

}

tree_plot<-function(gene_name_str, corr_type){
  subgroup_time_serie <-(data$genotype=='WT')
  samplesALL<-filter.n.normalize1(data,subgroup_time_serie)
  gene_name<-extend.gene.list_samples_from_hugo(sort(unique(gene_name_str)))
  mydata<-2**samplesALL$voom$E[gene_name,]
  e.thresh<-dim(mydata)[2]/2
  filter<-apply(mydata,1,function(x) ((sum(x>e.thresh))>1))
  
  S.correl <-rcorr(t(mydata[filter,]), type=corr_type)$r
  hc<-hclust(dist(S.correl), "mcquitty")
  
  num_colors<-round(length(rownames(mydata[filter,]))/2 + .5)
  if (num_colors>11){
    num_colors<-11
  }

  labelColors =c("mediumorchid4","red","limegreen","royalblue4","hotpink3","darkseagreen4","blue","darkslategray4","darkorange", "seagreen","darkgreen","darkblue")
  
  hcd = as.dendrogram(hc, horiz=T)
  clusMember = cutree(hc, num_colors)
  
  # function to get color labels
  colLab <- function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
      attr(n, "edgePar") <- list(col = labCol)
    }
    n
  }
  
  # using dendrapply
  clusDendro = dendrapply(hcd, colLab)
  
  # make plot
  par(mar=c(2,0,0,33), font=3) 
  plot(clusDendro, type = "rectangle", horiz=T) 
  
  
}



###USEFUL FUNCTIONS
findgn<-function(name,edata){ grep(name,rownames(edata$rawcounts),value=T) }
findhugo<-function(long_name){ data$annotations$hugo.name [data$annotations$name.and.blast.annotation==long_name]}
find_blast_annot <- function(long_name) { data$annotations['name.and.blast.annotation'][data$annotations['hugo.name']==long_name]}


############ Required functions ############ 
### What I need to do is do a file that has the hugo names without having to parse every time. 
### Data has a lot of useless information that I don't need. If I can figure out a was to make it faster that would be ideal.



get_website<-function(gene_name_input){
  #1. remove suffix
  short_gene <- gene_name_input
  #2. Obtaining Gene Description and substracting geneID for website
  GENEID_lookup <- gene_ManualCuration[substr(gene_ManualCuration[,2],1,nchar(short_gene))==short_gene][1]
  PAGE_ID <- gsub("XB-GENEPAGE-","",GENEID_lookup)
  if(is.na(GENEID_lookup)){}else{
    return(paste("http://www.xenbase.org/gene/expression.do?method=displayGenePageExpression&geneId=",PAGE_ID,"&tabId=1",sep=""))
  }
}

get_gene_description<-function(gene_name_input){
  #1. remove suffix
  short_gene <- gene_name_input
  #2. Obtaining Gene Description and substracting geneID for website
  GENEID_lookup <- gene_ManualCuration[substr(gene_ManualCuration[,2],1,nchar(short_gene))==short_gene][1]
  GENEID_description <- gene_ManualCuration[substr(gene_ManualCuration[,2],1,nchar(short_gene))==short_gene][3]
  if(is.na(GENEID_lookup)){}else{
    return(GENEID_description)
  }
}

get_insitu_links<-function(gene_name_input,stage){
  short_gene <- gene_name_input
  if (stage == '12'){
    INSITU_LINK12<- InSitu[substr(InSitu$Gene,1,nchar(short_gene))==short_gene &  (InSitu$Annot_stage == c("NF stage 11","NF stage 12","NF stage 12.5","NF stage 13"))]$ImgID
    
    #INSITU_LINK12<- InSitu[substr(InSitu[,2],1,nchar(short_gene))==short_gene & (InSitu[,5]=="NF stage 11" | InSitu[,5]=="NF stage 12" | InSitu[,5]=="NF stage 12.5" | InSitu[,5]=="NF stage 13") ,6]
    if (length(INSITU_LINK12)>0){
      #randomNum12 <- randomNumbers(n=1,min=1,max=length(INSITU_LINK12),col=1)
      randomNum12<-1
      insitu_website<- paste("http://www.xenbase.org/common/image/standard/",INSITU_LINK12[randomNum12],".jpg",sep="")
      return(insitu_website)
    }else{
      return()
    }
  }
  if (stage == '14'){
    INSITU_LINK14<- InSitu[substr(InSitu$Gene,1,nchar(short_gene))==short_gene &  (InSitu$Annot_stage == c("NF stage 14","NF stage 15","NF stage 16","NF stage 17"))]$ImgID
    
    #    INSITU_LINK14<- InSitu[substr(InSitu[,2],1,nchar(short_gene))==short_gene & (InSitu[,5]=="NF stage 17" | InSitu[,5]=="NF stage 14" |InSitu[,5]=="NF stage 15"|InSitu[,5]=="NF stage 16") ,6]
    if (length(INSITU_LINK14)>0 ){
      randomNum14<-1
      #randomNum14 <- randomNumbers(n=1,min=1,max=length(INSITU_LINK14),col=1)
      insitu_website<- paste("http://www.xenbase.org/common/image/standard/",INSITU_LINK14[randomNum14],".jpg",sep="")
      return(insitu_website)
    }else{
      return()
    }
  }
}

extend.gene.list<-function(gene_list){
  list_of_geness<-c()
  for (gene in gene_list){
    if (length(list_of_geness) == 0) {
      list_of_geness <- findgn(gene,data)
    }else{
      list_of_geness<-append(list_of_geness,findgn(gene,data))
    }
  }
  return(list_of_geness)
}

reduce.gene.list_samples<-function(gene_list){
  list_of_geness<-c()
  for (gene in gene_list){
    if (length(list_of_geness) == 0) {
      list_of_geness <- findhugo(gene)
    }else{
      list_of_geness<-append(list_of_geness,findhugo(gene))
    }
  }
  return(sort(unique(list_of_geness)))
}


extend.gene.list_samples<-function(gene_list,samples){
  list_of_geness<-c()
  for (gene in gene_list){
    if (length(list_of_geness) == 0) {
      list_of_geness <- findgn(gene,samples)
    }else{
      list_of_geness<-append(list_of_geness,findgn(gene,samples))
    }
  }
  return(list_of_geness)
}


extend.gene.list_samples_from_hugo<-function(gene_list){
  list_of_geness<-c()
  for (gene in gene_list){
    if (length(list_of_geness) == 0) {
      list_of_geness <- find_blast_annot(gene)
    }else{
      list_of_geness<-append(list_of_geness,find_blast_annot(gene))
    }
  }
  return(unique(c(list_of_geness)))
}


#Computes average of expression in each sample type (linear, not log2) so we can draw a pattern
expr.mean.comp<-function(edata){
  e.matrix<-2**edata$voom$E
  e.tissues<-sort(unique(edata$tissue))
  e.mean<-c()
  for (t in e.tissues)
  {
    t.col<-(edata$tissue==t)
    t.mean<-apply(e.matrix[,t.col],1,mean)
    e.mean<-cbind(e.mean,t.mean)
  }
  e.result<-cbind(e.mean)
  rownames(e.result)<-rownames(edata$voom$E)
  colnames(e.result)<-c(e.tissues)
  return(e.result)
}


obtain.CMP_and_enrichment.values<-function(data.meansONLY, gene_name){
  num_tissues <- dim(data.meansONLY)[2]-1  #obtain the number of tissues minus 1
  gene_name_1<-rownames(data.meansONLY)[match(gene_name,rownames(data.meansONLY))]
  
#  gene_name_1<-rownames(data.meansONLY)[match(extend.gene.list_samples_from_hugo(gene_name),rownames(data.meansONLY))]
  gene_name_1<-gene_name_1[is.na(gene_name_1)==FALSE]
  ###gene_name_1<-
  
  specificity_in_ectoderm <- data.meansONLY[rownames(data.meansONLY),]/data.meansONLY[rownames(data.meansONLY),][,'WE']
  specificity_matrix1 <- NULL
  specificity_matrix1 <- data.frame(index=rownames(specificity_in_ectoderm)) #NULL#data.frame(colnames=colnames(a))
  for (tissue in colnames(specificity_in_ectoderm)){   #for each different tissue
    specificity_matrix1[tissue]<-rank(specificity_in_ectoderm[,tissue])
  }
  rownames(specificity_matrix1)<-rownames(specificity_in_ectoderm)
  specificity_matrix1$index <-NULL
  
  # specificity_tissueRANK: Will provide the ratio between the percentile
  specificity_tissueRANK <- specificity_matrix1[gene_name_1,colnames(specificity_in_ectoderm)]/dim(specificity_matrix1)[1]*100
  
  to_include<-NULL
  integreated_data_to_plot<-NULL
  
  to_include<- melt(t(as.matrix(specificity_tissueRANK[gene_name_1,1:num_tissues])))   #1:9 prevents WE output
  colnames(to_include)<-c("Stage","Gene","Specificity")
  
  CPM_means <-data.meansONLY[gene_name_1,1:num_tissues] #1:9 prevents WE output
  integreated_data_to_plot<-melt(t(CPM_means))
  colnames(integreated_data_to_plot)<-c("Tissue","Gene","RPKM")
  integreated_data_to_plot['Ectodemal_enrichment']<-to_include$Specificity
  integreated_data_to_plot
  
  return(integreated_data_to_plot)
}


#Widely used function that selects genes with expression values above a treshold of 2^1=2cpm.
filter.n.normalize1<-function(edata,samples,e.thresh=1){
  #crop expression matrix and filter genes with counts under e.thresh in all but 2 samples
  filter<-apply(edata$rawcounts[,samples],1,function(x) ((sum(x>e.thresh))>1))
  rawcounts.n<-edata$rawcounts[filter,samples] /edata$annotations$gene.length[filter] *1000

  #normalize with edgeR/TMM+voom
  vdge<-NULL
  dge<-DGEList(counts=rawcounts.n)
  dge<-calcNormFactors(dge,method="TMM")
  design<-NULL
  design<-model.matrix(~edata$sampletype[samples])
  vdge<-voom(dge,design,plot=FALSE)
  
  #return a new list with all the information
  return(list(rawcounts=rawcounts.n,sampletype=edata$sampletype[samples],design=design,
              stage=edata$stage[samples],tissue=edata$tissue[samples],
              genotype=edata$genotype[samples],id=edata$id[samples],
              annotations=edata$annotations[filter,],
              edger=dge,voom=vdge))
}


obtain_adjacent_and_cor_values<- function(gene_name){
  tmpA<-Adj.genes[1,]==gene_name
  genes_list<- Adj.genes[,tmpA]
  values_list<- Adj.values[,tmpA]
  new_df<- data.frame(row.names = genes_list)
  new_df[[gene_name]] <- genes_list
  new_df[['Adjacency']] <- values_list
  gene_names<- rownames(new_df)
  
  subgroup_time_serie <-(data$genotype=='WT')
  samplesALL<-filter.n.normalize1(data,subgroup_time_serie)  
  e.matrix<-2**samplesALL$voom$E[gene_names,]
  
  new_df[['Pearson']]<-rcorr(t(e.matrix[gene_names,]), type='pearson')$r[,gene_name]
  new_df[['p-val Pear']]<-rcorr(t(e.matrix[gene_names,]), type='pearson')$P[,gene_name]
  new_df<-cbind(new_df, WGCNAgene_module[match(gene_names,WGCNAgene_module$Gene),][2]) 
  return(new_df)
  
}


make_Network <- function (selected_gene,p_value_treshold){
  cat("Retrieving expression data for each input gene")
  subgroup_time_serie <-(data$genotype=='WT')
  samplesALL<-filter.n.normalize1(data,subgroup_time_serie)
  mydata<-2**samplesALL$voom$E
  
  #Select genes with cpm that are expressed above two tresholds (the sum of all conditions > 40cpm or the max of cpm counts is at least 50)
  e.thresh<-round(dim(mydata)[2]/2)
  filter2<-apply(mydata,1,function(x) ( ((max(x)>50)  | (sum(x>e.thresh))>1)))   #dim: 79x10822   new table! (if >50cpr then 79x10073)
  DATA1<- log(t(mydata[filter2,]))   
  dim(DATA1)
  colDATA1<-colnames(DATA1)

  genes_you_want<-WGCNAgene_module$Gene[match(extend.gene.list_samples_from_hugo(selected_gene),WGCNAgene_module$Gene)]
  genes_you_want<-genes_you_want[is.na(genes_you_want)==FALSE]
  
  cat("Drawing network\n\n")
  if (length(genes_you_want)>1){
    long_list <- obtain_long_list(c(genes_you_want), p_value_treshold, mydata, colDATA1)
    MisLinks<-make_link_matrix(unique(long_list))
    MisNodes<-make_node_matrix(unique(long_list), colDATA1)
    MyClickScript <- 'd3.select(this).select("circle").transition()
    .duration(200000)
    .attr("r", 20);
    .style("color","black")
    d3.select(this).select("text").transition()
    .duration(100000)
    .attr("x", 20)
    .style("stroke-width", "30px")
    .style("font", "italic 160px arial")
    .style("color","black")
    .style("opacity", 1)' 
    
    cat("Obtained Links and Node inforamtion. Ready to plot!\n")
    cat(paste("  Links:",dim(MisLinks)[1],"\n"))
    cat(paste("  Nodes:",dim(MisNodes)[1],"\n"))
    
    LinkColors<-MisLinks
    LinkColors[LinkColors$value<=0,"colors"]<-"red"
    LinkColors[LinkColors$value>0,"colors"]<-"gray"
    
    return(forceNetwork(Links = MisLinks, Nodes = MisNodes, Source = "source",
                        Target = "target", Value = "value", NodeID = "name",
                        Group = "group", 
                        opacity = 0.9,
                        Nodesize = "size",  
                        fontFamily="arial", 
                        fontSize =24,         
                        #  linkDistance= JS("function(d){return d.value + 5}"),
                        linkWidth  = JS('function(l){return l.value < 0 ? "2" : ".5" }'),
                        charge =-50, 
                        linkColour = LinkColors$colors,
                        bounded=F, 
                        opacityNoHover = 0, 
                        zoom = T,  
                        clickAction = MyClickScript,
                        radiusCalculation = JS("Math.sqrt(d.nodesize*15)"))
    )
  }
}


obtain_long_list <- function(genes_you_want,p_value_treshold, mydata, colDATA1){
  first<-1
  for (gene in genes_you_want){
    keep<-obtain_adjacent_values2(gene, p_value_treshold, mydata, colDATA1)
    if (first == 1){
      full_list<-keep      
      first<-0
    }else{
      full_list <-rbind(full_list,keep)
    }
  }
  return(full_list)
}


obtain_adjacent_values2<- function(gene_name2,pvalue_treshold, mydata, colDATA1){
  corr_type<-'pearson'
  tmpB<-Adj.genes[1,]==gene_name2
  genes_list<- Adj.genes[,tmpB]
  values_list<- Adj.values[,tmpB]
  new_df<- data.frame(row.names = c(1:length(genes_list)))
  new_df$source <- genes_list
  new_df$target[1:length(genes_list)]<-gene_name2
  new_df[['Adjacency']] <- values_list
  
  e.matrix<-mydata[genes_list,]
  new_df$Pearson<-rcorr(t(e.matrix[genes_list,]), type=corr_type)$r[,gene_name2]
  new_df$p_val<-rcorr(t(e.matrix[genes_list,]), type=corr_type)$P[,gene_name2]
  new_df$color<-WGCNAgene_module$Module[match(genes_list,WGCNAgene_module$Gene)]
  new_df$p_val[is.na(new_df$p_val)] <- 0
  new_df$size[new_df$source==gene_name2] <- 10
  new_df$size[new_df$source!=gene_name2] <- 3
  
  new_df_sel<-new_df[which(new_df$p_val<=pvalue_treshold),]
  new_df_sel$Distance<-(1-abs(new_df_sel$Pearson))
  return(new_df_sel)
  
}

make_link_matrix<-function(long_list){
  new_long_list<-data.frame(row.names = c(1:dim(long_list)[1]))
  new_long_list$source<-match(long_list[,1],unique(long_list[,c(1)]))-1
  new_long_list$target<-match(long_list[,2],unique(long_list[,c(1)]))-1
  new_long_list$value<-long_list[,'Pearson']
  new_long_list$Adjacency<-NULL
  MisLinks<-new_long_list
  return(MisLinks)
}

make_node_matrix<-function(long_list, colDATA1){
  new_long_list2<-data.frame(row.names = c(1:dim(unique(long_list[1]))[1]))
  new_long_list2$name<-as.matrix(unique(long_list[1]))
  new_long_list2$group<-WGCNAgene_module$Module[match(new_long_list2$name,WGCNAgene_module$Gene)]
  new_long_list2$size[1:length(new_long_list2$name)]<-5
  new_long_list2$size[which(match(new_long_list2$name, long_list$target)>0)]<-20
  colnames(new_long_list2)<-c('name','group','size')
  return(new_long_list2)
}
