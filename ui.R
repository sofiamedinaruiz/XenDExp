library(shiny)
library(data.table, quietly=T)

###### Documentation XenDExp v.1.2
# README.txt for info on how to run the XenDExp app. 
# Respository: https://github.com/sofiamedinaruiz/XenDExp/
# To run the app: 
# You require to install the libraries found in helper.R
# To run:
#      runApp('The-Directory-for/XenDExp')
#setwd('/Users/sofiamedina/Documents/Harland\ Lab/GeneExpression2016/XenDExp')
#To update the app to the Shyserver:  
#library(rsconnect)
#rsconnect::deployApp('../XenDExp/')


#Load gene manes for for Input choices.
genenames<-as.list(fread("load/genenames12and14.txt", header=F))$V1
names<-as.list(fread("load/data_annotation_hugo_name.txt", header=F))$V1
#It is faster if it only reads it once: 
#write(unique(data$annotations$hugo.name), "data_annotation_hugo_name.txt")
#write(sort(unique(c(genenames.st12,genenames.st14))), "genenames12and14.txt")


shinyUI(fluidPage(
  theme = "bootstrap.css",
  titlePanel("",windowTitle ="XenDExp tool"),
  h1(strong("XenDExp"), align="center" , style="color:darkgreen" ), 
  h1("- Explore the spatiotemporal expression program of", em("Xenopus"), "embryos -",  align="center" , style="color:darkgreen"),br(),br(),
  p("The XenDExp tool provides the expression pattern prediction of most genes present in", em("Xenopus laevis"),", including pseudo-alleles coming from the short (.s) or long(.l) chromosomes."),
  p("We provide a lists of co-expressed genes in neurulation based on a single gene. Ideal for exploratory/screening experiments. "),

  sidebarLayout(
    sidebarPanel(
      
      h4(strong('1.'), 'Obtain expression pattern predictions'),
      
      fixedRow(
        column(12, wellPanel(
          checkboxInput("Show_expression_pattern_coexp", "1.a. Show expression pattern and co-expression data.", value =FALSE, width = NULL),
          checkboxInput("Show_coexp", "1.b. Show list of co-expressed genes (Max. 50 genes)", value =FALSE, width = NULL),
          
          
          
          # Select gene
          selectInput("gene", label = h5("Reference Gene "),
                      choices = genenames, selected = "bmp4.l[bmp4(0.0)]"),
          
          # Select gene
          selectInput("gene2", label = h5("Compare with..."),
                      choices = genenames, selected = "bmp4.s[bmp4(0.0)]"),br(),
          h5("Optional:"),
          checkboxInput("In_situ_option", "Show examples of published in situ data (retrieved from Xenbase.org).  Images at the bottom of the page.",  value =FALSE, width = NULL), 
          p("To prevent website from crashing while website is loading, the user needs to Check this box after loading is complete."),

          submitButton("Draw predictions")
          )
          
        )), br(),br(),
      
      fixedRow(
        h4(strong('2.'), 'Input for time series'),
        p("e.g. snai2, sox10, pax3, ventx, foxd3, odc1, myod1"),
        p("Press enter after typing gene name"),
        column(12, wellPanel(
          checkboxInput("Use_prev", "Include the two genes shown in the expression pattern prediction", value =FALSE, width = NULL),
          selectizeInput('genes_for_timeserie', 'Write down genes to search:', choices=names, multiple=TRUE),
          selectInput("Nomalization_gene", label = h5("Normalize againts a specific gene:"),  choices=c("no normalization","ddx3x.l[ddx3x(0.0)]","ddx3x.s[ddx3x(0.0)]","odc1.l[odc1(0.0)]","odc1.s[odc1(0.0)]","xelaev18004064m.g[eef1a1o(0.0)]","xetrov90028866m.1[eef1a1(0.0)]","xelaev18012630m.g[eef1a1o(0.0)]","xelaev18034420m.g[gapdh(0.0)]","xelaev18036725m.g[gapdh(0.0)]"), selected="no normalization"),
          submitButton("Draw predictions")
        ))
      ), br(),br(),

      fixedRow(
        h4(strong('3.'), 'Show expression levels, heatplot and co-expresion relationships'),
        p("e.g. snai2, sox10, pax3, ventx, foxd3, odc1, myod1"),
        column(12, wellPanel(
          radioButtons("Select_stage11", label = h5("Select stage of interest"),  choices=sort(c("12.5","14")), selected = ('14')),
          checkboxInput("Use_prev2", "Include all previously selected genes", value =FALSE, width = NULL),
          selectizeInput('gene_list_for_expression', h5('Input gene names:'), choices=names, multiple=TRUE),
          br(),
          h4(strong('3a.'), 'Show expression levels, heatplot and co-expresion relationships'),
          
          checkboxInput("Expression_levels_option", "3.a. Show expression levels", value = FALSE, width = NULL),
          sliderInput("px", h5("Adjust height (pixels)"), value = 1000, min = 300,max = 2000),

         br(),
        h4(strong('3.b.c.'), 'Update information for heatplot and dendrograms'),
         checkboxInput("Make_correlation_option", "3.b. Show correlation matrix", value = FALSE, width = NULL),
         selectInput('corr_type', 'Type of correlation', choices=c("spearman","pearson"), selected='pearson'),
         checkboxInput("Draw_tree_option", "3.c. Show distance tree ", value = FALSE, width = NULL),
         checkboxInput("Draw_Network", "3.d. Show Co-Expression Network (can take as long as 5 min, you can download the data instead)", value = FALSE, width = NULL),
         h5(" Would you like to also include landmark genes?"),
         checkboxInput("Landmark_genesNC", "Neural crest", value = FALSE, width = NULL ),
         checkboxInput("Landmark_genesNB", "Neural border", value = FALSE, width = NULL ),
         checkboxInput("Landmark_genesVentral", "Ventral", value = FALSE, width = NULL ),
         checkboxInput("Landmark_genesAnterior", "Anterior", value = FALSE, width = NULL ),
         checkboxInput("Landmark_genesNP", "Neural Plate", value = FALSE, width = NULL ),
         checkboxInput("Landmark_genesEpi", "Epidermis", value = FALSE, width = NULL ),
         sliderInput("px2", h5("Adjust height (pixels)"), value = 700, min = 300,max = 2000),
          submitButton("Submit")
        ))), br(),br(),
      
##      fixedRow(
##        column(12, wellPanel(
##          h4(strong('4.'), 'Download Network information'),
##          checkboxInput('Download_Network', 'Download Network',  value = FALSE, width = NULL),
##          submitButton("Submit")
##        )
##      )), br(),br(),
      
      fixedRow(
        
        column(12, wellPanel(
          # Select gene
  
          p(strong("Sample codes:")),
          p(strong("Eca:"), "Ectoderm - anterior"),
          p(strong("Ecl:"), "Ectoderm - lateral"),
          p(strong("Ecv:"), "Ectoderm - ventral"),
          p(strong("NPa:"), "Neural Plate anterior"),
          p(strong("NPp:"), "Neural Plate posterior"),
          p(strong("ANB:"), "Anterior Neural Border/Preplacodal ectoderm"),
          p(strong("NBa:"), "Neural Border anterior"),
          p(strong("NBp:"), "Neural Border posterior "),
          p(strong("NB:"), "Neural Border (ant+post)"), 
          p(strong("WE:"), "Whole Embryo"), 
          
          br(),
          br(),
          p(strong("Other acronyms:")),
          p(strong("CPM:"),"Counts per Million mapped reads"),
          p(strong("NF:"),"Nieuwkoop Faber Stages")
          
        )
        
      )), br(),br(),
       
      h5('About this project:'),
      p("This project was carried out by the", a("Monsoro Burq lab (Curie Institute)", href="http://umr3347.curie.fr/fr/equipes-de-recherche/ah-monsoro-burq",target="_blank"), "and the", a("Harland lab (UC Berkeley)",href="https://mcb.berkeley.edu/labs/harland/",target="_blank")),
      h6("REF: Manuscript of preparation. "),
      h6("* For more information about in situ expression patterns displayed in this page, go to the links from Xenbase.org"),
      br(),br(),br(),
      h3(strong("Help us improve this app!")),
      h5("Please send your feedback to:"),
      h5("sofiamr@berkeley.edu"),
      br(),
      h6("Most recent update:"), 
      h5("   4/14/16 - ",a("Sofia Medina",href="https://www.linkedin.com/in/sofia-medina-ruiz-3189904a",target="_blank")),
    width = 3),

    
  #######################   Main Panel  #######################    

    mainPanel(
      tabsetPanel(
        tabPanel("1. Expression pattern  predictions",
                 tabsetPanel(
                   
                             tabPanel("1a. NMF prediction",
                                      h3(strong("1."), "Compare expression patterns between two genes"), 
                                      h5(htmlOutput("Gene_link1")), uiOutput("xenPlot1.ui"),
                                      h5(htmlOutput("Gene_link2")), uiOutput("xenPlot2.ui"),br(),br()),
                             tabPanel("1.b Adjacency table",
                                      uiOutput("Text_output_coexp.ui"), 
                                      dataTableOutput("Adjacent_matrix1"),
                                      uiOutput("Text_output_coexp2.ui"),
                                      dataTableOutput("Adjacent_matrix2"),br(),br()),
                             tabPanel("1.c In situ from Xenbase",
                                      h3(strong("1.c."), "Examples of ",em("in situ"),"expression patterns from", em("Xenopus laevis"),"deposited in Xenbase.org:"),br(),
                                      htmlOutput("InSitus12"),br(),
                                      htmlOutput("InSitus14"),br(),br())
                 )
        ),
        
                 
        
        tabPanel("2. Time series", br(),br(),
                 h3(strong("2."), "Transcript levels during neurula stages"),
                 p(strong("Note:"), "Samples under 5 CPM (red dotted line) may not be detected by", em("in situ"),"hybridization."), 
                 uiOutput("timeSeries.ui"),br(),br()),
        
        tabPanel("3. Tissue expression",
                 #h3(strong("3."), "Expression patterns"), 
                 tabsetPanel(
                             tabPanel("3.a. Normalized CMP in all tissues",
                                h3(strong("3.a"), "Transcript levels across different ectodermal tissues", style="color:black"),
                                p(strong("Note:"), "Samples below 5 CPM may not be detected by", em("in situ "),"(hybridization dotted red line)."),
                                uiOutput("Absolute2.ui"), br(),
                                br(),br()),
                             tabPanel("3.b. Correlation",
                                 h3(strong("3.b."), "Correlation matrix", style="color:black"),
                                 p("Correlation of expression values across all WT data obtained from tissue dissections and whole embryos."),
                                 p("Output will only be generated for 3 or more transcripts."),
                                 p("p-value codes: [. = 0.1-0.01; * = 0.01-0.001; ** = 0.001-0.0001; *** = 0.0001-0.00001; **** = 0.00001-0]"),
                                 uiOutput("heatplot.ui"),br(),br(),
                                 #plotOutput("heatplot", width="90%"),#, height= 600),
                                 br(),br()),
                            tabPanel("3.c. Tree",
                                     h3(strong("3.c."), "Tree of expression distance based on correlations.", style="color:black"),
                                     p("Dendrogram of similarity in gene expression profiles across 79 wildtype samples (WE and tissues)."),
                                     uiOutput("tree.ui"),br(),br())
                 )
        ),
        
        tabPanel("3.d*. Co-expression Network",
                 h3(strong("3.d."), "Gene Network", style="color:black"),
                 p("The calculation and network construction may take up to 10 min",  style = "color:red"),
                 #,
                 uiOutput("Draw_Network.ui"),br(),br() 
                 )
      ),
      
      br(), 
    
    width = 9
    )
  )
))


###Need to fix: Warning: Removed 543 rows containing missing values (geom_text).
###Warning: Removed 543 rows containing missing values (geom_text). 
### From 3.b
####Another warning!:
####Found more than one class "connection" in cache; using the first, from namespace 'BiocGenerics'
#### ???
