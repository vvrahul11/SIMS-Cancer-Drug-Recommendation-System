library(shiny)
#setwd("/srv/shiny-server/sims/")
#library(googleVis)

# Define UI for random distribution application 
shinyUI(
  
  ### TO see the title in a seperate block, uncomment the below line and the bracket at the end of the code ##
  fluidPage(theme = "shiny.css",    
            list(tags$head(HTML('<link rel="icon", href="", 
                                type="image/jpg" />')
            )
            ),
            div(style="padding: 0px 0px; width: '100%'",
                titlePanel(
                  title="", windowTitle="SIMS - Simplified Interventional Mapping System"
                )
            ),
            navbarPage(title=div(tags$div(HTML('<img src="WINLOGO_backwhite.png" height="60" width="80" <br>
                                               <font size = "40"> <strong> SIMS: Simplified Interventional Mapping System </strong></font>'), collapsible = TRUE, fluid = FALSE                         
            )),   
            
            sidebarLayout(position = "left",
                          sidebarPanel(textInput("userid", "User ID:"),
                                       passwordInput("password", "Password:"),      
                                       actionButton("login", "Log in"),
                                       actionButton("signup", "Sign up"),
                                       p(""),  
                                       tags$hr(),
                                       tags$div(HTML('<font size = "4"><strong> SIMS score </strong></font>')),
                                       tags$div(HTML('')),
                                       fileInput("mytable_mrna","Upload a csv file containing mRNA fold change", accept = c('text/csv',
                                                                                                                            'text/comma-seperated-values',
                                                                                                                            'text/tab-seperated-values',
                                                                                                                            'text/plain',
                                                                                                                            '.csv')),
                                       conditionalPanel(condition="output.checkUpload!=null",
                                                        uiOutput("mRNA")    	
                                       ),                    
                                       fileInput("mytable_mirna","Upload a csv file containing miRNA fold change", accept = c('text/csv',
                                                                                                                              'text/comma-separated-values',
                                                                                                                              'text/tab-separated-values',
                                                                                                                              'text/plain',
                                                                                                                              '.csv')), # accept defines the file type accepted by the program
                                       conditionalPanel(condition="output.checkUpload!=null",
                                                        uiOutput("miRNA")  		
                                       ),
                                       fileInput("mytable_CNA","Upload a csv file containing CNA", accept = c('text/csv',
                                                                                                              'text/comma-separated-values',
                                                                                                              'text/tab-separated-values',
                                                                                                              'text/plain',
                                                                                                              '.csv')), # accept defines the file type accepted by the program
                                       conditionalPanel(condition="output.checkUpload!=null",
                                                        uiOutput("CNA")    	
                                       ),                    
                                       fileInput("mytable_mutation","Upload a csv file containing mutation", accept = c('text/csv',
                                                                                                                        'text/comma-separated-values',
                                                                                                                        'text/tab-separated-values',
                                                                                                                        'text/plain',
                                                                                                                        '.csv')), # accept defines the file type accepted by the program
                                       conditionalPanel(condition="output.checkUpload!=null",
                                                        uiOutput("Mutation")  		
                                       ),
                                       checkboxInput('header', 'Header', TRUE),
                                       radioButtons('sep', 'Separator',
                                                    c(Comma=',',
                                                      Semicolon=';',
                                                      Tab='\t'),
                                                    ','),      
                                       #numericInput("n", "N:", min = 0, max = 100, value = 50),
                                       #actionButton("goButton", "Go!"),
                                       actionButton("go", "Calculate score"),
                                       tags$hr(),
                                       p('Developed by', a(href = 'http://www.winconsortium.org/', 'Win Consortium')
                                       )           
                          ),
                          
                          # Show a tabset that includes a plot, summary, and table view
                          # of the generated distribution
                          mainPanel(#verbatimTextOutput("nText"),
                            tabsetPanel(#type = "tabs",              
                              tabPanel('Home', 
                                       p('                                                                                 '),
                                       p('                                                                                 '),
                                       tags$div(HTML('<p style = "text-align:justify" style = "font-family: Source Sans Pro;"><font size = "5"> <strong>Welcome to SIMS</strong></font></p>')),
                                       tags$div(HTML('
                                                     <p style = "text-align:justify" style = "font-family: Source Sans Pro;"><font size = "5"> <color = "red">Non-small cell lung cancer (NSCLC) is a leading cause of death worldwide. Targeted monotherapies produce high regression rates, albeit in 
                                                     only small oncogene-driven subsets of patients with metastatic NSCLC. Further, responders inevitably develop resistance and succumb. Combination
                                                     therapies may overcome resistance, but a scientific methodology for determining individualized combinations that can be used to prosecute lung
                                                     cancer remains an urgent unmet need. We present a novel strategy for identifying customized combinations of triple regimen targeted agents, deploying
                                                     a simplified interventional mapping system (SIMS) that merges knowledge about existent drugs and their impact on the hallmarks of cancer.
                                                     Based on interrogation of matched lung tumor and normal tissue (183 genes in 121 patients with NSCLC) and an integrative assessment of targeted
                                                     genomic sequencing, copy number variation, transcriptomics and miRNA expression, the activation status of 24 interventional nodes 
                                                     (target/gene/pathway(s)) were elucidated. An algorithm was developed to create a 1-10 scoring system, enabling ranking of the activated interventional
                                                     nodes in each patient. Based on the trends of co-activation at interventional points, tailored combinations of drug triplets were defined, which will
                                                     inform a prospective trial facilitated by cooperation between academia and industry in the Worldwide Innovative Network (WIN) consortium. This unique
                                                     approach will enable selection of personalized combination therapies in next generation clinical trials, whose objective should be to significantly
                                                     impact survival in metastatic NSCLC and other malignancies.</color></font></p>')                      
                                       ),
                                       tags$div(HTML('<p style = "text-align:justify" style = "font-family: Source Sans Pro;"><font size = "3"> Please cite: Lazar, V., et al. "A simplified interventional mapping system (SIMS) for the selection of combinations of targeted treatments in non-small cell lung cancer." Oncotarget (2015).</font></p>'))
                                       ),
                              navbarMenu("Score",
                                         tabPanel("Calculate Score", h4("SIMS Score"), tableOutput("calculateScore")),
                                         tabPanel("Analysis", plotOutput("scoreAnalysis"))
                                         #tabPanel("go", verbatimTextOutput("nText"))
                              ),
                              navbarMenu("mRNA",
                                         tabPanel("Contents", tableOutput("mRNAcontents")),
                                         tabPanel("Summary", verbatimTextOutput('mRNAsummary')),
                                         tabPanel("Plot", plotOutput("mRNAplot"))
                              ),
                              navbarMenu("miRNA",
                                         tabPanel("Contents", tableOutput("miRNAcontents")),
                                         tabPanel("Summary", verbatimTextOutput('miRNAsummary')),
                                         tabPanel("Plot", plotOutput("miRNAplot"))
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
