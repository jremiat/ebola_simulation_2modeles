#Appli de simulation avec les 2 modèles ensembles

library(reactable)
library(shiny)
library(deSolve)
library(tidyverse)
library(GGally)
library(fda)
library(dplyr)
library(DT)
library(pbapply)
library(shinythemes)
library(openxlsx)
library(shinydashboard)



##############
#     UI     #
##############
titre <-span(img(src = "inserm.jpg", height = 30), "Modèles de réponse immunitaire Ebola",img(src='Logo-Inria.png', height="60", width="100"))


ui <- tagList(
  dashboardPage(
  dashboardHeader(titleWidth = "800px",title = titre),
  
  dashboardSidebar(
    disable = TRUE),
  dashboardBody(
    tags$head(tags$style(HTML('
         .skin-blue .left-side, .skin-blue .wrapper {
                        background-color: #ecf0f5;
                        }
         '))),
    tabsetPanel(
      tabPanel(h4("Modèle avec Mdelay"),
               tabsetPanel(
                 tabPanel("Présentation du modèle",
                          
                          h1("Les variables du modèle",align="center"),
                          fluidRow(
                            column(3,HTML('<p><img src="Modèle.png"width="500" height="300"/></p>'),),
                            column(4,offset=4,
                                   HTML('<p>Les variables du modèles sont :</p>'),
                                   tags$li('A : la population des antigènes induits par la vaccination'),
                                   tags$li('M : la population de cellules B spécifiques (avec un temps de latence représenté par la variable Mdelay).'),
                                   tags$li("S : cellules sécrétrices d'anticorps à courte vie ."),
                                   tags$li("L : cellules sécrétrices d'anticorps à longue vie"),
                                   tags$li('Ab : les anticorps sécrétés par les populations S et L'),
                            ),),
                          tags$hr(),
                          HTML("<center>Tableau descriptif des différents paramètres</center>"),
                          HTML('<center><img src="Tableau Mdelay.png"width="600" height="400"/></center>'),
                          tags$hr(),
                          fluidRow(withMathJax("Dont les évolutions de concentration entre les injections i et i+1 sont modélisées par le modèle mécaniste suivant:
                      $$ \\begin{cases}\\frac{dM_{delay}}{dt}=\\rho_ie^{-\\delta_{A,i}(t-t_i)}-\\gamma_i M_{delay}
                                  \\\\ \\frac{dM}{dt}=\\gamma_i M_{delay}-(\\mu_{S,i}+\\mu_{L,i})e^{-\\delta_{A,i}(t-t_i)}M-\\delta_{M,i} M
                                  \\\\ \\frac{dS}{dt}=\\mu_{S,i}e^{-\\delta_{A,i}(t-t_i)}-\\delta_{S,i} S
                                  \\\\ \\frac{dL}{dt}=\\mu_{L,i}e^{-\\delta_{A,i}(t-t_i)}-\\delta_{L,i} L
                                  \\\\ \\frac{dAb}{dt}=\\theta_{S,i} S+\\theta_{L,i} L-\\delta_{Ab,i} Ab \\end{cases} $$")),
                          tags$hr(),
                          fluidRow(withMathJax(" Et dont les paramètres demandés sont transformés pour un individu j de la manière suivante :
                      $$\\begin{cases}\\log(\\rho^{(j)}=log(\\rho)+\\epsilon_{\\rho}^{(j)}
                      \\\\log(\\delta_A^{(j)})=log(\\delta_A)+\\epsilon_{\\delta_A}^{(j)}
                      \\\\log(\\delta_{Ab}^{(j)})=log(\\delta_{Ab})+\\epsilon_{\\delta_{Ab}}^{(j)}
                      \\\\log(\\delta_M^{(j)})=log(\\delta_M)+\\epsilon_{\\delta_M}^{(j)}
                      \\\\log(\\delta_S^{(j)})=log(\\delta_S)+\\epsilon_{\\delta_S}^{(j)}
                      \\\\log(\\delta_L^{(j)})=log(\\delta_L)+\\epsilon_{\\delta_L}^{(j)}
                      \\\\log(\\mu_S^{(j)})=log(\\mu_S)+\\epsilon_{\\mu_S}^{(j)}
                      \\\\log(\\mu_L^{(j)})=log(\\mu_L)+\\epsilon_{\\mu_L}^{(j)}
                      \\\\log(\\gamma^{(j)})=log(\\gamma)+\\epsilon_{\\gamma}^{(j)}
                      \\\\log(\\theta_S^{(j)})=log(\\theta_S)+\\epsilon_{\\theta_S}^{(j)}
                      \\\\log(\\theta_L^{(j)})=log(\\theta_S)+\\epsilon_{\\theta_L}^{(j)} 
                       \\end{cases} $$")),
                     
                 ),
                 tabPanel("Entrée des paramètres",withMathJax(),
                          tags$style(type="text/css", "#inline3 label{ display: table-cell; text-align: center; vertical-align: middle; } 
                 #inline3 .form-group { display: table-row;}"),
                          box(title=h3("Entrée rapide du schéma vaccinal", align = "center"),width = 12,collapsible = FALSE,solidHeader= TRUE,status='primary',
                              helpText("Vous pouvez entrer ici le schéma vaccinal avec les valeurs par défaut pour les 3 premières doses",align="center"),
                              fluidRow(column(width=12, offset=4,div(style="height: 40px; width: 200px",tags$div(id = "inline3",numericInput("nbinj", label = "nombre d'injections", value = 3,min =1,step=1,width="10%"))))), 
                              tags$hr(),                           
                              fluidRow(column(width=5,offset=1,uiOutput("tinj")),
                                       column(width=5,offset=0,uiOutput("aleatinj")),
                              ),
                              
                              helpText("Le type de vaccin fait varier \\(\\delta_A\\) et \\(\\delta_S\\), attention la modalité 'autre' ne prédéfinit pas de valeur pour ces paramètres",align="center"),
                              fluidRow(column(width=12, align="center", uiOutput("vaccin"))),
                              
                              fluidRow(column(width=12, align="center", sliderInput("times", label = "Temps simulation", min = 0, 
                                                                                    max = 3000, value =  1000))),
                              fluidRow(column(width=12, align="center", actionButton("action", label = "Go !",style='padding:4px; font-size:80%')),)),
                          
                          box(title=h3("Entrée détaillée des paramètres", align = "center"),width = 12,collapsible = TRUE,solidHeader= TRUE,
                            helpText("Toutes les variabilités sont représentées par des paramètres \\(\\epsilon_{i}\\sim N(0,\\sigma_{i}^{2})\\) et l'utilisateur doit entrer la valeur de \\(\\sigma_i\\). Si l'écart-type entré est 0, pas de variabilité. Cette variable aléatoire est ajoutée au logarithme des paramètres .",align="center"),
                            
                            helpText("\\(\\rho\\) : Taux de génération des lymphocytes B",align="center"),
                            fluidRow(column(3,offset=3,uiOutput("rho")),
                                     column(3,uiOutput("alealrho")),),
                            
                            helpText("\\(\\gamma\\) : Taux d'activation latent des lymphocytes B",align="center"),
                            fluidRow(column(3,offset=3,uiOutput("gamma")),
                                     column(3,uiOutput("alealgamma"))),
                            helpText(" \\(\\mu\\) : Taux de différenciation des lymphocytes B en plasmocytes ",align="center"),
                            
                            fluidRow(column(3,offset=3,uiOutput("muS")),
                                     column(3,uiOutput("alealmuS")),),
                            fluidRow(column(3,offset=3,uiOutput("muL")),
                                     column(3,uiOutput("alealmuL"),),
                            ),
                            helpText("\\(\\delta\\) : taux de dégradation des antigènes et anticorps et de mortalité des cellules immunitaires",align="center"),
                            fluidRow(column(3,offset=3,uiOutput("deltaA")),
                                     column(3,uiOutput("aleadeltaA"),),
                            ),
                            fluidRow(column(3,offset=3,uiOutput("deltaM")),
                                     column(3,uiOutput("aleadeltaM"),),
                            ),
                            fluidRow(column(3,offset=3,uiOutput("deltaS")),
                                     column(3,uiOutput("aleadeltaS"),),
                            ),
                            fluidRow(column(3,offset=3,uiOutput("deltaL")),
                                     column(3,uiOutput("aleadeltaL"),),
                            ), 
                            
                            fluidRow(column(3,offset=3,uiOutput("deltaAb")),
                                     column(3,uiOutput("aleadeltaAb"),),
                            ), 
                            
                            helpText(" \\(\\theta\\) : taux de production des anticorps par les différents plasmocytes",align="center"),
                            fluidRow(column(3,offset=3,uiOutput("thetaS")),
                                     column(3, uiOutput("alealthetaS"),),
                            ),
                            fluidRow(column(3,offset=3,uiOutput("thetaL")),
                                     column(3,uiOutput("alealthetaL"),),
                            ),
                            
                            
                            
                          ),
                          
                 ),
                 tabPanel("Graphes",
                          # textOutput("vec"),
                          waiter::use_waiter(),
                          
                          
                          textInput("nommodel", label="nom de cette simulation", value ="", width = "20%", placeholder = NULL),
                          actionButton("keep", label = "Garder pour comparer",width="20%",style='padding:4px; font-size:80%'),
                          fluidRow(
                            box(status='info',width=6,title="Graphe de M", plotOutput("grapheM")),  
                            box(status='info',width=6,title="Graphe de S",plotOutput("grapheS"),)
                          ),
                          fluidRow(
                            box(status='info',width=6,title="Graphe de L", plotOutput("grapheL")),
                            box(status='info',width=6,title="Graphe de Ab",  plotOutput("grapheAb"),),
                            
                          ),
                          downloadButton("allgraphs", "Download",width="50%",style='padding:4px; font-size:80%')
                          
                 ),
                 tabPanel("Comparaison aux autres simulations",withMathJax(),
                          box(status='info',width=12,title="Tableau des paramètres des simulation enregistrées",
                          reactableOutput("table_param"),
                          actionButton("clear", label = "Clear",style='padding:4px; font-size:80%')),
                        
                          fluidRow(
                            box(status='info',width=6,title="Graphe de M",
                                plotOutput("grapheMcompare"),
                                checkboxInput(inputId="logM", label="log-transform", value = FALSE)),
                            box(status='info',width=6,title="Graphe de S",
                                plotOutput("grapheScompare"),
                                checkboxInput(inputId="racS", label=paste0("\\(\\sqrt[4]{S}\\)"), value = FALSE)),),
                          fluidRow(
                            box(status='info',width=6,title="Graphe de L",
                                plotOutput("grapheLcompare"),
                                checkboxInput(inputId="racL", label=paste0("\\(\\sqrt[4]{L}\\)"), value = FALSE)),
                            box(status='info',width=6,title="Graphe de Ab",
                                plotOutput("grapheAbcompare"),
                                checkboxInput(inputId="logAb", label="log-transform", value = FALSE),)),
                          downloadButton("dl","Export in Excel",width="50%",style='padding:4px; font-size:80%')
                 ),
                 
               )),
      tabPanel(h4("Modèle avec régulation des populations"),
               tabsetPanel(
                 tabPanel("Présentation du modèle",
                          
                          h1("Les variables du modèle",align="center"),
                          fluidRow(
                            column(3,HTML('<p><img src="Modèle1.png"width="400" height="200"/></p>'),),
                            column(4,offset=4,
                                   HTML('<p>Les variables du modèles sont :</p>'),
                                   tags$li('A : la population des antigènes induits par la vaccination'),
                                   tags$li('M : la population de cellules B spécifiques.'),
                                   tags$li("S : cellules sécrétrices d'anticorps à courte vie ."),
                                   tags$li("L : cellules sécrétrices d'anticorps à longue vie"),
                                   tags$li('Ab : les anticorps sécrétés par les populations S et L'),
                            ),),
                          tags$hr(),
                          HTML("<center>Tableau descriptif des différents paramètres</center>"),
                          HTML('<center><img src="Tableau homeostasie.png"width="600" height="400"/></center>'),
                          tags$hr(),
                          
                          fluidRow(withMathJax("Dont les évolutions de concentration entre les injections i et i+1 sont modélisées par le modèle mécaniste suivant :
                      
                      $$\\begin{cases} \\frac{dM}{dt}=\\rho_i(M)e^{-\\delta_{A,i}(t-t_i)}-(\\mu_{S,i}(S)+\\mu_{L,i}(L))e^{-\\delta_{A,i}(t-t_i)}M-\\delta_{M,i} M
                      \\\\ \\frac{dS}{dt}=\\mu_{S,i}(S)e^{-\\delta_{A,i}(t-t_i)}M-\\delta_{S,i} S
                      \\\\ \\frac{dL}{dt}=\\mu_{L,i}(L)e^{-\\delta_{A,i}(t-t_i)}M-\\delta_{L,i} L
                      \\\\ \\frac{dAb}{dt}=\\theta_{S,i} S+\\theta_{L,i} L-\\delta_{Ab,i} Ab \\end{cases} $$")),
                          
                          tags$hr(),
                          fluidRow(withMathJax("Et dont les paramètres utilisés sont transformés par les fonctions suivantes :
                      $$\\begin{cases} \\rho(M)=\\rho_i\\frac{\\beta_i M+1}{\\alpha_{M,i} M^2+1}
                      \\\\ \\mu_S(S)=\\frac{\\mu_{S,i}}{\\alpha_{S,i} S^2+1}
                      \\\\ \\mu_S(S)=\\frac{\\mu_{L,i}}{\\alpha_{L,i} L^2+1}\\end{cases} $$")),
                          fluidRow(withMathJax("Pour pouvoir prendre en compte :")),
                          withMathJax(tags$li("le phénomène d'auto-régulation des populations avec les paramètres \\(\\alpha\\)")),
                          withMathJax(tags$li("le phénomène de prolifération de lymphocytes en db=ébut de réaction immunitaire avec le paramètre \\(\\beta\\)")),
                          
                 ),
                 tabPanel("Entrée des paramètres",withMathJax(),
                          tags$style(HTML("
  input[type=\"number\"] {
    height: 35px;}")),
                          tags$style(type="text/css", "#inline3 label{ display: table-cell; text-align: center; vertical-align: middle; } 
                 #inline3 .form-group { display: table-row;}"),
                          box(title=h3("Entrée rapide du schéma vaccinal", align = "center"),width = 12,collapsible = FALSE,solidHeader= TRUE,status='primary',
                              helpText("Vous pouvez entrer ici le schéma vaccinal avec les valeurs par défaut pour les 3 premières doses",align="center"),
                              fluidRow(column(width=12, offset=4,div(style="height: 40px; width: 200px",tags$div(id = "inline3",numericInput("hnbinj", label = "nombre d'injections", value = 3,min =1,step=1,width="10%"))))), 
                              tags$hr(),                           
                              fluidRow(column(width=5,offset=1,uiOutput("htinj")),
                                       column(width=5,offset=0,uiOutput("haleatinj")),
                              ),
                              
                              helpText("Le type de vaccin fait varier \\(\\delta_A\\) et \\(\\delta_S\\), attention la modalité 'autre' ne prédéfinit pas de valeur pour ces paramètres",align="center"),
                              fluidRow(column(width=12, align="center", uiOutput("hvaccin"))),
                              
                              fluidRow(column(width=12, align="center", sliderInput("htimes", label = "Temps simulation", min = 0, 
                                                                                    max = 3000, value =  1000))),
                              fluidRow(column(width=12, align="center", actionButton("haction", label = "Go !",style='padding:4px; font-size:80%')),)),
                          
                          box(title=h3("Entrée détaillée des paramètres", align = "center"),width = 12,collapsible = TRUE,solidHeader= TRUE,
                            helpText("Toutes les variabilités sont représentées par des paramètres \\(\\epsilon_{i}\\sim N(0,\\sigma_{i}^{2})\\) et l'utilisateur doit entrer la valeur de \\(\\sigma_i\\). Si l'écart-type entré est 0, pas de variabilité. Cette variable aléatoire est ajoutée au logarithme des paramètres ."),
                            
                            helpText("\\(\\rho\\) : Taux de génération des lymphocytes B",align="center"),
                            fluidRow(column(3,offset=3,uiOutput("hrho")),
                                     column(3,uiOutput("halealrho")),),
                            fluidRow(column(3,offset=3,uiOutput("hbeta"),),
                                     column(3,uiOutput("haleabeta")),),
                            fluidRow(column(3,offset=3,uiOutput("halphaM"),),
                                     column(3,uiOutput("haleaalphaM")),),
                            
                            
                            
                            helpText(" \\(\\mu\\) : Taux de différenciation des lymphocytes B en plasmocytes ",align="center"),
                            
                            fluidRow(column(3,offset=3,uiOutput("hmuS")),
                                     column(3,uiOutput("halealmuS")),),
                            fluidRow(column(3,offset=3,uiOutput("halphaS")),
                                     column(3,uiOutput("haleaalphaS")),),
                            fluidRow(column(3,offset=3,uiOutput("hmuL")),
                                     column(3,uiOutput("halealmuL"),),),
                            fluidRow(column(3,offset=3,uiOutput("halphaL")),
                                     column(3,uiOutput("haleaalphaL")),),
                            helpText("\\(\\delta\\) : taux de dégradation des antigènes et anticorps et de mortalité des cellules immunitaires",align="center"),
                            fluidRow(column(3,offset=3,uiOutput("hdeltaA")),
                                     column(3,uiOutput("haleadeltaA"),),
                            ),
                            fluidRow(column(3,offset=3,uiOutput("hdeltaM")),
                                     column(3,uiOutput("haleadeltaM"),),
                            ),
                            fluidRow(column(3,offset=3,uiOutput("hdeltaS")),
                                     column(3,uiOutput("haleadeltaS"),),
                            ),
                            fluidRow(column(3,offset=3,uiOutput("hdeltaL")),
                                     column(3,uiOutput("haleadeltaL"),),
                            ), 
                            
                            fluidRow(column(3,offset=3,uiOutput("hdeltaAb")),
                                     column(3,uiOutput("haleadeltaAb"),),
                            ), 
                            
                            helpText(" \\(\\theta\\) : taux de production des anticorps par les différents plasmocytes",align="center"),
                            fluidRow(column(3,offset=3,uiOutput("hthetaS")),
                                     column(3, uiOutput("halealthetaS"),),
                            ),
                            fluidRow(column(3,offset=3,uiOutput("hthetaL")),
                                     column(3,uiOutput("halealthetaL"),),
                            ),
                            
                            
                            
                          ),
                          ),
                 
                 tabPanel("Graphes",
                          waiter::use_waiter(),
                          # uiOutput("printMyDynamicInputs"),
                          # textOutput("text"),
                          #textOutput("hvar"),
                          #textOutput("hvec"),
                          #dataTableOutput("table"),
                          textInput("hnommodel", label="nom de cette modélisation", value ="", width = "20%", placeholder = NULL),
                          actionButton("hkeep", label = "Garder pour comparer",width="20%",style='padding:4px; font-size:80%'),
                          fluidRow(
                            box(status='info',width=6,title="Graphe de M", plotOutput("hgrapheM")),
                            box(status='info',width=6,title="Graphe de S",  plotOutput("hgrapheS"),)
                          ),
                          fluidRow(
                            box(status='info',width=6,title="Graphe de L", plotOutput("hgrapheL")),
                            box(status='info',width=6,title="Graphe de Ab",  plotOutput("hgrapheAb"),),
                            
                          ),
                          downloadButton("hallgraphs", "Download",width="50%",style='padding:4px; font-size:80%')
                          
                 ),
                 tabPanel("Comparaison aux autres simulations",withMathJax(),

                          reactableOutput("htable_param"),
                          actionButton("hclear", label = "Clear",style='padding:4px; font-size:80%'),
                          
                          fluidRow(
                            box(status='info',width=6,title="Graphe de M",
                                plotOutput("hgrapheMcompare"),
                                checkboxInput(inputId="hlogM", label="log-transform", value = FALSE)),
                            box(status='info',width=6,title="Graphe de S",
                                plotOutput("hgrapheScompare"),
                                checkboxInput(inputId="hracS", label=paste0("\\(\\sqrt[4]{S}\\)"), value = FALSE)),),
                          fluidRow(
                            box(status='info',width=6,title="Graphe de L",
                                plotOutput("hgrapheLcompare"),
                                checkboxInput(inputId="hracL", label=paste0("\\(\\sqrt[4]{L}\\)"), value = FALSE)),
                            box(status='info',width=6,title="Graphe de Ab",
                                plotOutput("hgrapheAbcompare"),
                                checkboxInput(inputId="hlogAb", label="log-transform", value = FALSE),)),
                          downloadButton("hdl","Export in Excel",width="50%",style='padding:4px; font-size:80%')
                 ),
               )),
      tabPanel(h4("Comparaison des modèles"),
               tabsetPanel(
                 tabPanel("Tableaux des paramètres",
                          tags$li('Modélisations enregistrées pour le modèle avec Mdelay'),
                          reactableOutput("table_param2"),
                          tags$li('Modélisations enregistrées pour le modèle avec régulation des populations cellulaires'),
                          reactableOutput("htable_param2"),),
                 tabPanel("Graphes des simulations sélectionnées.",
                          #textOutput("noms"),
                          #textOutput("hnoms"),
                          #reactableOutput("tente"),
                          fluidRow(
                            box(status='info',width=6,title="Graphe de M",
                                plotOutput("grapheMtotal"),
                                checkboxInput(inputId="logMtotal", label="log-transform", value = FALSE)),
                            box(status='info',width=6,title="Graphe de S",plotOutput("grapheStotal"),
                                checkboxInput(inputId="racStotal", label=paste0("\\(\\sqrt[4]{S}\\)"), value = FALSE)),),
                          fluidRow(
                            box(status='info',width=6,title="Graphe de L",
                                plotOutput("grapheLtotal"),
                                checkboxInput(inputId="racLtotal", label=paste0("\\(\\sqrt[4]{L}\\)"), value = FALSE)),
                            box(status='info',width=6,title="Graphe de Ab",
                                plotOutput("grapheAbtotal"),
                                checkboxInput(inputId="logAbtotal", label="log-transform", value = FALSE),)),
                          downloadButton("compare", "Download",width="50%",style='padding:4px; font-size:80%')
                 ),
               )),
    )
  )
),
tags$footer(HTML("Modèles et valeurs des paramètres basés sur l'article :<br/>Balelli I, Pasin C, Prague M, Crauste F, Effelterre TV, Bockstal V, et al. A model for establishment, maintenance and reactivation of the immune response after vaccination against Ebola virus. J Theor Biol. 2020;495:110254. pmid:32205143")),
            #tags$a(href="www.rstudio.com", "Click here!",target="_blank")),
)
