library(shiny)
library(data.table, quietly=T)

###### Documentation EctoMap v.1.2.0 - 04/2017-1
# README.txt for info on how to run the EctoMap app. 
# Respository: https://github.com/sofiamedinaruiz/EctoMAP/
# To run the app: 
# You require to install the libraries found in helper.R
# To run:
#      runApp('The-Directory-for/EctoMAP')
# e.g. runApp('~/Dropbox/EctoMAP figures/DOTPLOTS AHMB/EctoMAP v.1.2.1')
#setwd('/Users/sofiamedina/Documents/Harland\ Lab/GeneExpression2016/EctoMAP')
#To update the app to the Shyserver:  
#library(rsconnect)
#rsconnect::deployApp('/Users/sofiamr/Dropbox/EctoMAP\ figures/DOTPLOTS\ AHMB/EctoMAP\ copy')
###### 


#Load gene manes for for Input choices.
#genenames<-unique(rownames(data$annotations['name.and.blast.annotation']))
#write(genenames, 'load/genenames.txt')

load(file = "load/ALLdata.rda")

genenames<-fread('load/genenames.txt', header=F)$V1
names<-fread('load/data_annotation_hugo_name.txt', header=F)$V1


shinyUI(fluidPage(
  theme = "bootstrap.css",
  titlePanel("",windowTitle ="EctoMap tool"),
  h1(strong("EctoMap v1.2"), align="center" , style="color:darkblue" ), 
  h1("- Explore the spatiotemporal expression program of", em("Xenopus"), "embryos -",  align="center" , style="color:darkblue"),br(),br(),
  p("This tool provides the expression profile of genes expressed in", em("Xenopus laevis"),"during embryonic development."),
  p("This app has two main panels, one for input queries and the other for data plots. The results for the genes queried on the left panel, will be displayed inside the tabs on the right that matches the type of desired query (eg. time series, expression pattern, etc...)."),
  p("The results of each query are provided in four separate tabs. The first section, provides the developmental progression of gene expression throughout neurula stages. The second, section displays the expression pattern prediction of the genes input on the left panel."),
  br(),
  h5("Dissected ectodermal tissues:"),
 # tags$img(src="Tissue_dissection.jpg", type = "jpeg", width="50%", align="right", alt="Tissue dissections"),#, height = 200, width = 200")),

  sidebarLayout(
    sidebarPanel(
      fixedRow(
        h3(strong('Gene input panel')),
        p('This long tab will serve as gene query inputs. Please read the descriptions of each of the sections on the right tabs.'),
        p("Note: Please press [enter] after typing a gene name"),
        h6(strong("Please wait while data is being loaded.")
        ),br(),
        h4(strong('1.'), 'Time series plot'),
        column(12, wellPanel(
          #checkboxInput("Show_time_series", "Show expression levels during neurula stages from Xenopus laevis embryos (St. 11-19)", value =TRUE, width = NULL),
   ########       selectizeInput('genes_for_timeserie', 'Input one or more gene names:', choices=names, multiple=TRUE, selected='pax3'),
          selectizeInput("genes_for_timeserie", label = "Input one or more transcript names:", choices = genenames, multiple=TRUE, selected = c("bmp4.l[bmp4(0.0)]","bmp4.s[bmp4(0.0)]")),
        #  p("e.g. ",em("pax3, myod1")),
          selectInput("Nomalization_gene", label = h5("Normalize by:"),  choices=c("no normalization","ddx3x.l[ddx3x(0.0)]","ddx3x.s[ddx3x(0.0)]","odc1.l[odc1(0.0)]","odc1.s[odc1(0.0)]","xelaev18004064m.g[eef1a1o(0.0)]","xetrov90028866m.1[eef1a1(0.0)]","xelaev18012630m.g[eef1a1o(0.0)]","xelaev18034420m.g[gapdh(0.0)]","xelaev18036725m.g[gapdh(0.0)]"), selected="no normalization"),
          submitButton("Draw predictions")
        )#wellPanel
        )#column
      ),#fixed row 
      br(),
      
      fixedRow(
        h4(strong('2.'), 'Expression pattern predictions'),
        column(12, wellPanel(
          p("Expression pattern predictions for Stage 12.5 and St. 14. "),br(),
           # Select gene1
          selectInput("gene", label = p("Transcript #1:"),
                      choices = genenames, selected = "bmp4.l[bmp4(0.0)]"),
          checkboxInput("Show_coexp", "Show genes that co-express with Transcript 1 (2.b)", value =FALSE, width = NULL),
          checkboxInput("In_situ_option", "If available, show an example of a published in situ image deposited in www.xenbase.org (2.c)",  value =FALSE, width = NULL), 
          
          # Select gene2
          selectInput("gene2", label = p("Transcript #2:"),
                      choices = genenames, selected = "bmp4.s[bmp4(0.0)]"),br(),
          checkboxInput("Show_expression_pattern_coexp", "Show the expression pattern prediction", value =TRUE, width = NULL),
          submitButton("Draw predictions")
          )
        )#column
      ),#fixedRow
      br(),
      
      fixedRow(
        h4(strong('3.'), 'Mean expression levels and correlations'),
        column(12, wellPanel(
    #      br(), 
          radioButtons("Select_stage11", label = h5("Select the stage to obtain the mean expression values:"),  choices=sort(c("12.5","14")), selected = ('14')),
          h5(strong('Please select the genes you would like to include:')),
          checkboxInput("Use_prev2", "Include genes plotted in Time Series (section #1)", value =FALSE, width = NULL),
          checkboxInput("Use_prev3", "Include NMF prediction transcripts (section #2)", value =FALSE, width = NULL),
          p('You can also select genes or transcripts from previous sections.'),
          checkboxInput("Expression_levels_option", "Show the average RPKM values from the ectodermal tissue dissections (Section 2)", value = FALSE, width = NULL),
          selectizeInput('gene_list_for_expression', h5('Genes names (optional):'), choices=names, multiple=TRUE ),
           br(),
  #        h5(strong('3.b.'), 'Correlation plots'),
  #        p('This plots are created with the genes and/or transcripts chosen above'),
  #        checkboxInput("Make_correlation_option", "3.b.1. Show correlation matrix", value = FALSE, width = NULL),
  #        checkboxInput("Draw_tree_option", "3.b.2. Show distance tree ", value = FALSE, width = NULL),
  #OK        checkboxInput("Draw_Network", "3.d. Show Co-Expression Network (can take as long as 5 min, you can download the data instead)", value = FALSE, width = NULL),
  #        h6("For correlation-based plots include landmark genes:"),
  #        checkboxInput("Landmark_genesNC", "Neural crest", value = FALSE, width = NULL ),
  #        checkboxInput("Landmark_genesNB", "Neural border", value = FALSE, width = NULL ),
  #        checkboxInput("Landmark_genesVentral", "Ventral", value = FALSE, width = NULL ),
  #        checkboxInput("Landmark_genesAnterior", "Anterior", value = FALSE, width = NULL ),
  #        checkboxInput("Landmark_genesNP", "Neural Plate", value = FALSE, width = NULL ),
  #        checkboxInput("Landmark_genesEpi", "Epidermis", value = FALSE, width = NULL ),
  #        checkboxInput("Landmark_genesPos", "Posterior", value = FALSE, width = NULL ),
  #        selectInput('corr_type', 'Type of correlation', choices=c("spearman","pearson"), selected='pearson'),
          submitButton("Submit")
        )#wellPanel
       )#column
      ),#fixedRow
      br(),
  
      fixedRow(
        h4(strong('4.'), 'Co-Expression Network'),
        column(12, wellPanel(
          
          p('All the genes input from sections 1 through 4 will become input nodes (including landmark genes). The number of associated connected to the input nodes will be determined by the p-value cutoff. Please give some time for the Network to be made'),
          checkboxInput('Draw_Network', 'Display co-expression Network',  value = FALSE, width = NULL),
          h6("Landmark genes to include in the network:"),
          checkboxInput("Landmark_genesNC", "Neural crest", value = FALSE, width = NULL ),
          checkboxInput("Landmark_genesNB", "Neural border", value = FALSE, width = NULL ),
          checkboxInput("Landmark_genesVentral", "Ventral", value = FALSE, width = NULL ),
          checkboxInput("Landmark_genesAnterior", "Anterior", value = FALSE, width = NULL ),
          checkboxInput("Landmark_genesNP", "Neural Plate", value = FALSE, width = NULL ),
          checkboxInput("Landmark_genesEpi", "Epidermis", value = FALSE, width = NULL ),
          checkboxInput("Landmark_genesPos", "Posterior", value = FALSE, width = NULL ),
          selectInput('corr_type', 'Type of correlation', choices=c("spearman","pearson"), selected='pearson'),
          
          selectInput('p_value_treshold', 'Type of correlation', choices=c("1e-3","1e-4","1e-5","1e-6","1e-7","1e-8","1e-9","1e-9","1e-10","1e-11","1e-12","1e-13","1e-14"), selected='1e-10'),
          submitButton("Submit")
        )#Well panel
        )#column
      ),#fixedRow 
      br(),
      
      fixedRow(
        column(12, wellPanel(
          p(strong("Sample codes of dissected tissues:")),
          p(strong("Eca:"), "Ectoderm - anterior"),
          p(strong("Ecl:"), "Ectoderm - lateral"),
          p(strong("Ecv:"), "Ectoderm - ventral"),
          p(strong("NPa:"), "Neural Plate anterior"),
          p(strong("NPp:"), "Neural Plate posterior"),
          p(strong("ANB:"), "Anterior Neural Border/Preplacodal ectoderm"),
          p(strong("NBa:"), "Neural Border anterior"),
          p(strong("NBp:"), "Neural Border posterior "),
          p(strong("WE:"), "Whole Embryo")
        )#well panel
        )#column
      ),#fixed row
      br(),br(),
      
      p("This project was carried out by the", a("Monsoro Burq (Curie Institute)", href="http://umr3347.curie.fr/fr/equipes-de-recherche/ah-monsoro-burq",target="_blank"), ",", a("Harland  (UC Berkeley)",href="https://mcb.berkeley.edu/labs/harland/",target="_blank"),', and',a("Eisen  (UC Berkeley/HHMI)",href="http://eisenlab.org/",target="_blank"),'labs.'),
      h6("Please cite: Plouhinec and Medina-Ruiz et al., 2017 - In review"),
      br(),br(),br(),
      h5("Send your comments or feedback to:"),
      h5(p("Prof. Monsoro-Burq - anne-helene.monsoro-burq@curie.fr")),

      br(),
      h6("Most recent update:"), 
      h5("   04/2017 - ",a("SMR",href="https://www.linkedin.com/in/sofia-medina-ruiz-3189904a",target="_blank")),
    
    width = 3),#sidebarPanel //sidebarPanel
    
  #######################   Main Panel  #######################    

    mainPanel(
      
      tabsetPanel(
        tabPanel("1. Time series plot", br(),br(),
##1          h3(strong("1."), "Transcript levels during neurula stages"),
          p("Each data point represents the expression level (RFPKM units) of a particular transcript, collected from a single whole embryo at a determined Nieuwkoop Faber (NF) developmental stage."),
          p("The human gene(s) input from the left panel will be shown as title. However, if genes happens to be duplicated in the X. laevis genome, each of the transcripts  will be differentiated by color in the figure legend."),br(),
          uiOutput("timeSeries.ui"),
          p(strong("RPKM:"),"Reads Per Kilobase of transcript per Million"),br(),
          HTML('<center><img src="Temporal_expression.jpg" width="90%"></center>'), 
          br()
        ),
          
      tabPanel("2. Expression pattern  predictions",
        tabsetPanel(
          tabPanel("2.a. Ectodermal expression patterns prediction",
###2            h3(strong("2.a."), "Ectodermal expression patterns prediction"), 
            br(),p(strong("About this section:")," Seven ectodermal domain tissue samples were obtained from ",em("Xenopua leavis")," embryos at Stages 12.5 and 14. The colored diagram show the different dissected tissues viewed from the dorsal and lateral sides."),br(),
            HTML('<center><img src="Disection_DiagramAH.png" width="80%"></center>'), 
            br(),

            p("Non-Negative Matrix Factorization was used to predict the expression patterns by reducing data dimensionality [Plouhinec and Medina-Ruiz et al., 2017]"),br(),
            
            h5(htmlOutput("Gene_link1")), uiOutput("xenPlot1.ui"),
            h5(htmlOutput("Gene_link2")), uiOutput("xenPlot2.ui"),
            h1(textOutput("test1102")),br()),
          tabPanel("2.b. List of co-expressed transcripts", 
            br(),
            p(strong("About this section:"),"Weighted correlation network analysis (WGCNA) [Langfelder and Horvath, 2008. BMC], was used to cluster genes into modules based correlations. A signed and unsigned adjacency matrices were generated using the expression values of all transcripts from both spatial and temporal datasets. The co-expression modules were obtained by introducing the signed adjacency matrix and the following parameters to the WGCNA analysis package: soft power = 22, method = hybrid, deepSplit = 4 [Plouhinec and Medina-Ruiz et al., 2017] ."),br(),
            #p("The following table will display the 50 transcripts that show expression profile similarities with transcipt #1 (Section 2). This list was obtained by retrieving the top-50 most adjacent transcripts from the unsigned adjacency matrix (2nd column). The Pearson's correlation and associated p-value (3rd and 4th column) result from the pairwise comparisons between two genes across all 79 spatio/temporal datasets. The last column displays the assigned module number, which can range from 0-140"),
            uiOutput("Text_output_coexp.ui"), 
            dataTableOutput("Adjacent_matrix1"),
            uiOutput("Text_output_coexp2.ui"),
            dataTableOutput("Adjacent_matrix2"),br(),br()),
          tabPanel("2.c. In situ from Xenbase",
            h3(strong("2.c."), "Examples of ",em("in situ"),"expression patterns from", em("Xenopus laevis"),"deposited in Xenbase.org:"),br(),
            htmlOutput("InSitus12"),br(),
            htmlOutput("InSitus14"),br(),br())
        )
      ),
          
  
      tabPanel("3. Mean expression per tissue",
        #h3(strong("3."), "Expression patterns"), 
#####        tabsetPanel(
#tP          tabPanel("3.a. Mean expression values per tissue",
            #h3(strong("3.a"), "Transcript levels across different ectodermal tissues", style="color:black"),
              br(),p(strong("About this section:"), "This section is complementary to the NMF pattern predictions, however it provides additional information about transcript enrichment in the ectoderm. The graphs obtained here show the mean transcript levels (RFPKM) across different ectodermal tissues described below."),
              p("Ectodermal enrichment values are ranked from 0-100, and are represented by a color gradient (orange through green). 
                       Dark-green is equivalent to an ectodermal enrichment value of 100, the transcript is likely to be enriched in ectodermal tissues. 
                       Orange means the transcript could be enriched elsewhere (e.g. endoderm or mesoderm). 
                       The ectodermal enrichment value was obtained by calculating the ratio between the mean expression of a given tissue over 
                       the mean expression of the whole embryo (same stage), and further normalized by the highest ectodermal enrichment 
                       value and multiplied by 100%."),
              #p("Samples below 1 RPKM may not be detected by", em("in situ "),"(hybridization gray dotted line)."),
              p(strong("Note:"), "Dinamically adjust height of the dot-plot at the bottom of the page"),
            
              uiOutput("Absolute2.ui"), 
              br(),
              HTML('<center><img src="Disection_DiagramAH.png" width="90%"></center>'), 

              wellPanel(
                h5(strong("Abreviation of tissue samples, Stage 12.5:")),
                p(strong("Non-Neural Ectoderm --> "), strong("Eca: "), "anterior ectoderm , ",strong(" Ecl:"), "lateral ectoderm , ",strong(" Ecv:"), "ventral ectoderm "), 
                p(strong("Neural Border --------> "), strong("ANB: "), "Anterior Neural Border/Preplacodal ectoderm, ", strong(" NBa:"), "Neural Border anterior, ",strong(" NBp:"),"Neural Border posterior "),
                p(strong("Neural Plate ---------> "), strong("NPa: "), "Neural Plate anterior, ",strong(" NPp:"), "Neural Plate posterior")
              ),br(),
              wellPanel(
                h5(strong("Abreviation of tissue samples, Stage 14:")),
                p(strong("Epidermis (E) ----> "), strong("Ea: "), "anterior epidermis , ",strong(" Ep:"), "posterior epidermis"), 
                p(strong("Neural Border (NB) ----> "), strong("AANB: "), "Anterior Neural Border/Preplacodal ectoderm, ", strong(" NBa:"), "Neural Border anterior, ",strong(" NBp:"),"Neural Border posterior "),
                p(strong("Neural Plate (NP) ----> "), strong("NPa: "), "Neural Plate anterior, ",strong(" NPp:"), "Neural Plate posterior")
              ),
              
              
              sliderInput("px", h5("Adjust pixel height for Dotplot"), value = 1000, min = 300,max = 2500),
              submitButton("OK"),
              br(),br()
      ), #tabPanel #3     
  #       tabPanel("3.b. Correlation",
  #           h3(strong("3.b."), "Correlation matrix", style="color:black"),
  #           p("Correlation of expression values across all WT data obtained from tissue dissections and whole embryos."),
  #           p("Output will only be generated for 3 or more transcripts."),
  #           p("p-value codes: [. = 0.1-0.01; * = 0.01-0.001; ** = 0.001-0.0001; *** = 0.0001-0.00001; **** = 0.00001-0]"),
  #           uiOutput("heatplot.ui"),br(),br(),
  #           sliderInput("px3", h5("Adjust pixel height for Heatplot"), value = 700, min = 300,max = 2500),
  #           submitButton("OK"),
  #           #   plotOutput("heatplot", width="90%"),#, height= 600),
  #           br(),br()),
  #       tabPanel("3.c. Dendrogram",
  #           h3(strong("3.c."), "Distance of expression based on correlations.", style="color:black"),
  #           p("Dendrogram of similarity in gene expression profiles across 79 wildtype samples (WE and tissues)."),
  #           uiOutput("tree.ui"),br(),br(),
  #           sliderInput("px4", h5("Adjust pixel height for Dengrogram"), value = 1000, min = 300,max = 2500),
  #           submitButton("OK"),br())
 #tP       )#tabsePanel
 #     ), #tabPanel #3     
          
      tabPanel("4. Co-expression Network",
        h3(strong("4."), "Gene Network", style="color:black"),  br(),br(),
        p("Note: The visualization of the network visualization may take up to 5 min.",  style = "color:gray"),
        uiOutput("Draw_Network.ui"),br(),br() 
      ) #tabPanel #4
  
    ), #tabsetPanel
        
    br(), 
    width = 9
    )#mainPanel
  )#sidebarLayout
) #fluidPage
) #shinyUI





