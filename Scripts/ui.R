library(shiny)

# Define UI for random distribution application 
shinyUI(

### TO see the title in a seperate block, uncomment the below line and the bracket at the end of the code ##
fluidPage(#theme = "shiny.css",
       
       list(tags$head(HTML('<link rel="icon", href="", 
                                    type="image/jpg" />'))),
      div(style="padding: 0px 0px; width: '100%'",
          titlePanel(
            title="", windowTitle="My Window Title"
          )
     ),
      navbarPage(#div(img(src="WIN.jpg")),
        title=div(tags$a(HTML('<img src="WINLOGO_backwhite.png" height="50" width="80" 
                              <font color = "white"> <strong> SIMS: Simplified Interventional Mapping System </strong></font>')                         
                        )),

    sidebarLayout(
    sidebarPanel(
      fileInput("mytable","Upload a csv file containing mRNA fold change"),
      #conditionalPanel(condition="output.checkUpload!=null",
       #                uiOutput("experiments")			
      #),
      tags$p(""),
      fileInput("mytable_miRNA","Upload a csv file containing miRNA fold change"),
      #conditionalPanel(condition="output.checkUpload!=null",
      #                uiOutput("experiments")  		
      #),
      tags$p("")
    ),
    
    # Show a tabset that includes a plot, summary, and table view
    # of the generated distribution
    mainPanel(
      tabsetPanel(#type = "tabs",                   
                  tabPanel('Home', 
                           p('Non-small cell lung cancer (NSCLC) is a leading cause of death worldwide. Targeted monotherapies produce high regression rates, albeit in only small oncogene-driven subsets of patients with metastatic NSCLC. Further, responders inevitably develop resistance and succumb. Combination therapies may overcome resistance, but a scientific methodology for determining individualized combinations that can be used to prosecute lung cancer remains an urgent unmet need. We present a novel strategy for identifying customized combinations of triple regimen targeted agents, deploying a simplified interventional mapping system (SIMS) that merges knowledge about existent drugs and their impact on the hallmarks of cancer. Based on interrogation of matched lung tumor and normal tissue (183 genes in 121 patients with NSCLC) and an integrative assessment of targeted genomic sequencing, copy number variation, transcriptomics and miRNA expression, the activation status of 24 interventional nodes (target/gene/pathway(s)) were elucidated. An algorithm was developed to create a 1-10 scoring system, enabling ranking of the activated interventional nodes in each patient. Based on the trends of co-activation at interventional points, tailored combinations of drug triplets were defined, which will inform a prospective trial facilitated by cooperation between academia and industry in the Worldwide Innovative Network (WIN) consortium. This unique approach will enable selection of personalized combination therapies in next generation clinical trials, whose objective should be to significantly impact survival in metastatic NSCLC and other malignancies.',
                             align = "justify",
                             style = "font-family: 'Source Sans Pro';")
                           ),
                  navbarMenu("mRNA",
                             tabPanel("Contents", tableOutput("mRNAcontents")),
                             tabPanel("Summary", verbatimTextOutput('mRNAsummary')),
                             tabPanel("Plot", plotOutput("mRNAplot"))
                             ),
                  navbarMenu("miRNA",
                             tabPanel("Contents", tableOutput("miRNAcontents")),
                             tabPanel("Summary", verbatimTextOutput('miRNAsummary')),
                             tabPanel("Plot", plotOutput("plot"))
                  ) ,
                  navbarMenu("CNA",
                             tabPanel("Contents", tableOutput("CNAcontents")),
                             tabPanel("Summary", verbatimTextOutput('CNAsummary')),
                             tabPanel("Plot", plotOutput("CNAplot"))
                  ) ,
                  navbarMenu("Mutation",
                             tabPanel("Contents", tableOutput("Mutation_contents")),
                             tabPanel("Summary", verbatimTextOutput('Mutation_summary')),
                             tabPanel("Plot", plotOutput("Mutation_plot"))
                  )
                  )
                  
      )
    )
  )
)
)
