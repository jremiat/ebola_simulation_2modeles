#Appli de simulation avec les 2 modèles ensembles


#couleur bleue claire : #8AA9DB


#########################################
#                                       #
#           Code avant serveur          #
#                                       #
#########################################

library(shiny)
library(reactable)
library(deSolve)
library(tidyverse)
library(GGally)
library(fda)
library(dplyr)
library(DT)
library(data.table)
library(waiter)
library(openxlsx)
library(viridis)
#library(pbapply)
###################################################################
#                                                                 #
#   création des objets et fonctions pour le modèle avec Mdelay   #
#                                                                 #
###################################################################
#fonction qui prend les paramètres, établit une valeur selon laa variable aléatoire et fait le champ de vecteur des équas diff du modèle :
model.MdMSLAb <-function (t,x,parms){
  Mdelay <- x[1]
  M <-x[2]
  S <- x[3]
  L <- x[4]
  Ab <- x[5]
  nbinj<-as.numeric(parms[1])
  if (length(parms)!=12*nbinj){
    stop("erreur dans le vecteur des paramètres")}
  if (NA %in% parms){
    stop("veuillez renseigner toutes les valeurs demandées")
  }
  else {
    
    
    dMdelay <- parms[2]*exp(-parms[3]*t)-parms[4]*Mdelay
    dM <- parms[4]*Mdelay-(parms[5]+parms[6])*exp(-parms[3]*t)*M-parms[7]*M
    dS <- parms[5]*exp(-parms[3]*t)*M-parms[8]*S
    dL <- parms[6]*exp(-parms[3]*t)*M-parms[9]*L
    dAb <- parms[10]*S+parms[11]*L-parms[12]*Ab
    if (nbinj>=2){
      for (k in 2:nbinj){
        if (t>=as.numeric(parms[13+(k-2)*12])){
          dMdelay <- parms[14+(k-2)*12]*exp(-parms[15+(k-2)*12]*(t-parms[13+(k-2)*12]))-parms[16+(k-2)*12]*Mdelay
          dM <- parms[16+(k-2)*12]*Mdelay-(parms[17+(k-2)*12]+parms[18+(k-2)*12])*exp(-parms[15+(k-2)*12]*(t-parms[13+(k-2)*12]))*M-parms[19+(k-2)*12]*M
          dS <- parms[17+(k-2)*12]*exp(-parms[15+(k-2)*12]*(t-parms[13+(k-2)*12]))*M-parms[20+(k-2)*12]*S
          dL <- parms[18+(k-2)*12]*exp(-parms[15+(k-2)*12]*(t-parms[13+(k-2)*12]))*M-parms[21+(k-2)*12]*L
          dAb <- parms[22+(k-2)*12]*S+parms[23+(k-2)*12]*L-parms[24+(k-2)*12]*Ab
        }
      }
    }
    dX <- c(dMdelay,dM,dS,dL,dAb)
    list(dX)
  }
}

#création de la fonction qui ajoute les effets aléatoires et résoud l'ode :
resolve.ode<-function(i,param,times,xstart){
  nbinj<-as.numeric(param[1])
  parms<-vector(length=12*nbinj)
  parms[1]<-nbinj
  parms[2]<-exp(log(as.numeric(param[2]))+rnorm(1,0,as.numeric(param[3])))
  parms[3]<-exp(log(as.numeric(param[4]))+rnorm(1,0,as.numeric(param[5])))
  parms[4]<-exp(log(as.numeric(param[6]))+rnorm(1,0,as.numeric(param[7])))
  parms[5]<-exp(log(as.numeric(param[8]))+rnorm(1,0,as.numeric(param[9])))
  parms[6]<-exp(log(as.numeric(param[10]))+rnorm(1,0,as.numeric(param[11])))
  parms[7]<-exp(log(as.numeric(param[12]))+rnorm(1,0,as.numeric(param[13])))
  parms[8]<-exp(log(as.numeric(param[14]))+rnorm(1,0,as.numeric(param[15])))
  parms[9]<-exp(log(as.numeric(param[16]))+rnorm(1,0,as.numeric(param[17])))
  parms[10]<-exp(log(as.numeric(param[18]))+rnorm(1,0,as.numeric(param[19])))
  parms[11]<-exp(log(as.numeric(param[20]))+rnorm(1,0,as.numeric(param[21])))
  parms[12]<-exp(log(as.numeric(param[22]))+rnorm(1,0,as.numeric(param[23])))
  if (nbinj>=2){
    for (i in 2:nbinj){
      parms[13+12*(i-2)]<-as.numeric(param[24+24*(i-2)])+rnorm(1,0,as.numeric(param[25+24*(i-2)]))
      parms[14+12*(i-2)]<-exp(log(as.numeric(param[26+24*(i-2)]))+rnorm(1,0,as.numeric(param[27+24*(i-2)])))
      parms[15+12*(i-2)]<-exp(log(param[28+24*(i-2)])+rnorm(1,0,param[29+24*(i-2)]))
      parms[16+12*(i-2)]<-exp(log(as.numeric(param[30+24*(i-2)]))+rnorm(1,0,as.numeric(param[31+24*(i-2)])))
      parms[17+12*(i-2)]<-exp(log(as.numeric(param[32+24*(i-2)]))+rnorm(1,0,as.numeric(param[33+24*(i-2)])))
      parms[18+12*(i-2)]<-exp(log(as.numeric(param[34+24*(i-2)]))+rnorm(1,0,as.numeric(param[35+24*(i-2)])))
      parms[19+12*(i-2)]<-exp(log(as.numeric(param[36+24*(i-2)]))+rnorm(1,0,as.numeric(param[37+24*(i-2)])))
      parms[20+12*(i-2)]<-exp(log(as.numeric(param[38+24*(i-2)]))+rnorm(1,0,as.numeric(param[39+24*(i-2)])))
      parms[21+12*(i-2)]<-exp(log(as.numeric(param[40+24*(i-2)]))+rnorm(1,0,as.numeric(param[41+24*(i-2)])))
      parms[22+12*(i-2)]<-exp(log(as.numeric(param[42+24*(i-2)]))+rnorm(1,0,as.numeric(param[43+24*(i-2)])))
      parms[23+12*(i-2)]<-exp(log(as.numeric(param[44+24*(i-2)]))+rnorm(1,0,as.numeric(param[45+24*(i-2)])))
      parms[24+12*(i-2)]<-exp(log(as.numeric(param[46+24*(i-2)]))+rnorm(1,0,as.numeric(param[47+24*(i-2)])))
    }
  }
  simul<- ode(
    func=model.MdMSLAb,
    y=xstart,
    times=times,
    parms=parms
  ) %>%
    as.data.frame()
}
#création de la fonction qui teste si présence de variabilité et si oui, resoud 100 fois et si non, résoud 1 fois 
test.resoud <- function(var,i,param,times,xstart){
  
}

#point de départ de la simulation
xstart <- c(Mdelay=0,M=0, S=0, L=0, Ab=0)

#création des tables dans lesquelles sont enregistrées les simulation et paramètres associés
parametres <- data.frame(matrix(nrow=0,ncol=1))
colnames(parametres)<-c("nom")
param<-c()
simul_brute <- data.frame()
simulation_MC <- data.frame(matrix(nrow=0,ncol=18))
simulation_brute <-data.frame()
ordre<-c()

colnames(simulation_MC)<-c("time","Mdelay","M","S","L","Ab","Mdelayinf","Minf","Sinf","Linf","Abinf","Mdelaysup","Msup","Ssup","Lsup","Absup","nom","modele")

#############################################################################################
#
#fonctions et objets crées pour le modèle homéostatique
#
###############################################################################################
#fonction et des équas diff du modèle homéostatique:
model.homeostatic <- function(t,x,parms){
  M <-x[1]
  S <- x[2]
  L <- x[3]
  Ab <- x[4]
  nbinj<-as.numeric(parms[1])
  if (length(parms)!=15*nbinj){
    stop("erreur dans le vecteur des paramètres")}
  if (NA %in% parms){
    stop("veuillez renseigner toutes les valeurs demandées")
  }
  else {
    
    dM<-parms[2]*(parms[3]*M+1)/(parms[4]*M^2+1)*exp(-parms[5]*t)-(parms[6]/(parms[7]*S^2+1)+parms[8]/(parms[9]*L^2+1))*exp(-parms[5]*t)*M-parms[10]*M
    dS <- parms[6]/(parms[7]*S^2+1)*exp(-parms[5]*t)*M-parms[11]*S
    dL <- parms[8]/(parms[9]*L^2+1)*exp(-parms[5]*t)*M-parms[12]*L
    dAb <- parms[13]*S+parms[14]*L-parms[15]*Ab
    if (nbinj>=2){
      for (k in 2:nbinj){
        if (t>=as.numeric(parms[16+(k-2)*15])){
          dM<-parms[17+(k-2)*15]*(parms[18+(k-2)*15]*M+1)/(parms[19+(k-2)*15]*M^2+1)*exp(-parms[20+(k-2)*15]*(t-parms[16+(k-2)*15]))-(parms[21+(k-2)*15]/(parms[22+(k-2)*15]*S^2+1)+parms[23+(k-2)*15]/(parms[24+(k-2)*15]*L^2+1))*exp(-parms[20+(k-2)*15]*(t-parms[16+(k-2)*15]))*M-parms[25+(k-2)*15]*M
          dS <- parms[21+(k-2)*15]/(parms[22+(k-2)*15]*S^2+1)*exp(-parms[20+(k-2)*15]*(t-parms[16+(k-2)*15]))*M-parms[26+(k-2)*15]*S
          dL <- parms[23+(k-2)*15]/(parms[24+(k-2)*15]*L^2+1)*exp(-parms[20+(k-2)*15]*(t-parms[16+(k-2)*15]))*M-parms[27+(k-2)*15]*L
          dAb <- parms[28+(k-2)*15]*S+parms[29+(k-2)*15]*L-parms[30+(k-2)*15]*Ab
        }
      }
    }
    dX <- c(dM,dS,dL,dAb)
    list(dX)
  }
}


#création de la fonction qui ajoute les effets aléatoires et résoud l'ode :
hresolve.ode<-function(i,param,times,hxstart){
  nbinj<-as.numeric(param[1])
  parms<-vector(length=15*nbinj)
  parms[1]<-nbinj
  parms[2]<-exp(log(as.numeric(param[2]))+rnorm(1,0,as.numeric(param[3])))
  parms[3]<-exp(log(as.numeric(param[4]))+rnorm(1,0,as.numeric(param[5])))
  parms[4]<-exp(log(as.numeric(param[6]))+rnorm(1,0,as.numeric(param[7])))
  parms[5]<-exp(log(as.numeric(param[8]))+rnorm(1,0,as.numeric(param[9])))
  parms[6]<-exp(log(as.numeric(param[10]))+rnorm(1,0,as.numeric(param[11])))
  parms[7]<-exp(log(as.numeric(param[12]))+rnorm(1,0,as.numeric(param[13])))
  parms[8]<-exp(log(as.numeric(param[14]))+rnorm(1,0,as.numeric(param[15])))
  parms[9]<-exp(log(as.numeric(param[16]))+rnorm(1,0,as.numeric(param[17])))
  parms[10]<-exp(log(as.numeric(param[18]))+rnorm(1,0,as.numeric(param[19])))
  parms[11]<-exp(log(as.numeric(param[20]))+rnorm(1,0,as.numeric(param[21])))
  parms[12]<-exp(log(as.numeric(param[22]))+rnorm(1,0,as.numeric(param[23])))
  parms[13]<-exp(log(as.numeric(param[24]))+rnorm(1,0,as.numeric(param[25])))
  parms[14]<-exp(log(as.numeric(param[26]))+rnorm(1,0,as.numeric(param[27])))
  parms[15]<-exp(log(as.numeric(param[28]))+rnorm(1,0,as.numeric(param[29])))
  if (nbinj>=2){
    for (i in 2:nbinj){
      parms[16+15*(i-2)]<-as.numeric(param[30+30*(i-2)])+rnorm(1,0,as.numeric(param[31+30*(i-2)]))
      parms[17+15*(i-2)]<-exp(log(as.numeric(param[32+30*(i-2)]))+rnorm(1,0,as.numeric(param[33+30*(i-2)])))
      parms[18+15*(i-2)]<-exp(log(param[34+30*(i-2)])+rnorm(1,0,param[35+30*(i-2)]))
      parms[19+15*(i-2)]<-exp(log(as.numeric(param[36+30*(i-2)]))+rnorm(1,0,as.numeric(param[37+30*(i-2)])))
      parms[20+15*(i-2)]<-exp(log(as.numeric(param[38+30*(i-2)]))+rnorm(1,0,as.numeric(param[39+30*(i-2)])))
      parms[21+15*(i-2)]<-exp(log(as.numeric(param[40+30*(i-2)]))+rnorm(1,0,as.numeric(param[41+30*(i-2)])))
      parms[22+15*(i-2)]<-exp(log(as.numeric(param[42+30*(i-2)]))+rnorm(1,0,as.numeric(param[43+30*(i-2)])))
      parms[23+15*(i-2)]<-exp(log(as.numeric(param[44+30*(i-2)]))+rnorm(1,0,as.numeric(param[45+30*(i-2)])))
      parms[24+15*(i-2)]<-exp(log(as.numeric(param[46+30*(i-2)]))+rnorm(1,0,as.numeric(param[47+30*(i-2)])))
      parms[25+15*(i-2)]<-exp(log(as.numeric(param[48+30*(i-2)]))+rnorm(1,0,as.numeric(param[49+30*(i-2)])))
      parms[26+15*(i-2)]<-exp(log(as.numeric(param[50+30*(i-2)]))+rnorm(1,0,as.numeric(param[51+30*(i-2)])))
      parms[27+15*(i-2)]<-exp(log(as.numeric(param[52+30*(i-2)]))+rnorm(1,0,as.numeric(param[53+30*(i-2)])))
      parms[28+15*(i-2)]<-exp(log(as.numeric(param[54+30*(i-2)]))+rnorm(1,0,as.numeric(param[55+30*(i-2)])))
      parms[29+15*(i-2)]<-exp(log(as.numeric(param[56+30*(i-2)]))+rnorm(1,0,as.numeric(param[57+30*(i-2)])))
      parms[30+15*(i-2)]<-exp(log(as.numeric(param[58+30*(i-2)]))+rnorm(1,0,as.numeric(param[59+30*(i-2)])))
    }
  }
  simul<- ode(
    func=model.homeostatic,
    y=hxstart,
    times=times,
    parms=parms
  ) %>%
    as.data.frame()
}





#création des valeurs de départs, tableaux vides etc...


hxstart <- c(M=0, S=0, L=0, Ab=0)
hparametres <- data.frame(matrix(nrow=0,ncol=1))
colnames(hparametres)<-c("nom")
hsimul_brute <- data.frame()
hsimulation_MC <- data.frame(matrix(nrow=0,ncol=15))
hsimulation_brute <- data.frame()
hparam<-c()
colnames(hsimulation_MC)<-c("time","M","S","L","Ab","Minf","Sinf","Linf","Abinf","Msup","Ssup","Lsup","Absup","nom","modele")
hordre<-c()

#création de la fonction qui regarde si une valeur est dans une colonne de data frame 
#utilisé dans les 2 modèles
identique<- function(i,rang,tableau){
  return(which(tableau[,i+1] ==rang[i]))
}



####creation tableau réunissant les 2 modèles
comparaison <- data.frame(matrix(nrow=0,ncol=15))

#création de la fonction identique: 
identique<- function(i,rang,tableau){
  return(which(tableau[,i+1] ==rang[i]))
}


######################################
#                                    #
#               Serveur              #
#                                    #
######################################
shinyServer(function(input, output) {
  
  
  #liste de tous les renderUI à créer pour le modèle Mdelay:
  #uioutput("vaccin),
  output$vaccin <- renderUI({ lapply(1:input$nbinj, function(i) {
    column(3,
           withMathJax(radioButtons(paste0("vaccin", i), label = paste0("type de vaccin", i), c("Ad26"=1,"MVA"=2,"Autre"=3))
           ))})
    
  })
  # uiOutput("tinj"),
  output$tinj <- renderUI({
    nbinj <- as.integer(input$nbinj)
    if (nbinj>=2){
      lapply(2:nbinj, function(i) {
        
        numericInput(paste0("tinj", i), label = paste0("jour de l'injection", i), value = switch(i,0,56,365),min=0,step=1)
      })
    }
  })
  # uiOutput("aleatinj"),
  output$aleatinj <- renderUI({
    nbinj <- as.integer(input$nbinj)
    if (nbinj>=2){
      lapply(2:nbinj, function(i) {
        
        numericInput(paste0("aleatinj", i), label = paste0("\\(\\sigma\\)(jour de l'injection)", i), value = 0,min=0)
      })
    }
  })
  # uiOutput("rho"),
  output$rho <- renderUI({ lapply(1:input$nbinj, function(i) {
    
    withMathJax(numericInput(paste0("rho", i), label = paste0(" \\(\\rho\\)", i), value = switch(i,4.4,36.6,167.3),step=0.5,min=0)
    )})
    
  })
  # uiOutput("alealrho"),
  output$alealrho <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("alealrho", i), label = paste0("\\(\\sigma_{l\\rho}\\)", i), value = 0,step=0.1,min=0)
      )})
    
  })
  # uiOutput("deltaA"),
  output$deltaA <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("deltaA", i), label = paste0(" \\(\\delta_{A}\\)", i), value =ifelse(input[[paste0("vaccin",i)]]==1,0.064,ifelse(input[[paste0("vaccin",i)]]==2,0.21,0)),step=0.05,min=0)
      )})
    
  })
  # uiOutput("aleadeltaA"),
  output$aleadeltaA <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("aleadeltaA", i), label = paste0("\\(\\sigma_{\\delta_{A}}\\)", i), value = 0,min=0,step=0.05)
      )})
  })
  # uiOutput("gamma"),
  output$gamma <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("gamma", i), label = paste0(" \\(\\gamma\\)", i), value = 1,min=0)
      )})
    
  })
  # uiOutput("alealgamma"),
  output$alealgamma <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("alealgamma", i), label = paste0("\\(\\sigma_{\\gamma}\\)", i), value = 0,min=0)
      )})
  })
  # uiOutput("muS"),
  output$muS <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("muS", i), label = paste0(" \\(\\mu_{S}\\)", i), value =switch(i,0.13,1.28,0.25),step=0.1,min=0)
      )})
    
  })
  # uiOutput("alealmuS"),
  output$alealmuS <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("alealmuS", i), label = paste0("\\(\\sigma_{l\\mu_{S}}\\)", i), value =0,min=0,step=0.1)
      )})
  })
  # uiOutput("muL"),
  output$muL <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("muL", i), label = paste0(" \\(\\mu_{L}\\)", i), value =14*10^(-4),step=0.0001,min=0)
      )})
  })
  # uiOutput("alealmuL"),
  output$alealmuL <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("alealmuL", i), label = paste0("\\(\\sigma_{l\\mu_{L}}\\)", i), value =0,step=1,min=0)
      )})
  })
  # uiOutput("deltaM"),
  output$deltaM <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("deltaM", i), label = paste0(" \\(\\delta_{M}\\)", i), value =1/(60*365.25),step=0.00001,min=0)
      )})
  })
  # uiOutput("aleadeltaM"),
  output$aleadeltaM <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("aleadeltaM", i), label = withMathJax(paste0("\\(\\sigma_{\\delta_{M}}\\)"), i), value =0,min=0,step=0.00001)
      )})
  })
  # uiOutput("deltaS"),
  output$deltaS <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("deltaS", i), label = paste0(" \\(\\delta_{S}\\)", i), value =ifelse(input[[paste0("vaccin",i)]]==1,0.34,ifelse(input[[paste0("vaccin",i)]]==2,0.23,0)),step=0.1,min=0)
      )})
  })
  # uiOutput("aleadeltaS"),
  output$aleadeltaS <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("aleadeltaS", i), label = withMathJax(paste0("\\(\\sigma_{\\delta_{S}}\\)"), i), value =0,min=0,step=0.1)
      )})
  })
  # uiOutput("deltaL"),
  output$deltaL <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("deltaL", i), label = paste0(" \\(\\delta_{L}\\)", i), value =1/(8.5*365.25),step=0.0001,min=0)
      )})
  })
  # uiOutput("aleadeltaL"),
  output$aleadeltaL <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("aleadeltaL", i), label = withMathJax(paste0("\\(\\sigma_{\\delta_{L}}\\)"), i), value =0,min=0,step=0.0001)
      )})
  })
  # uiOutput("thetaS"),
  output$thetaS <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("thetaS", i), label = paste0(" \\(\\theta_{S}\\)", i), value =13.5,min=0)
      )})
  })
  # uiOutput("alealthetaS"),
  output$alealthetaS <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("alealthetaS", i), label = withMathJax(paste0("\\(\\sigma_{l\\theta_{S}}\\)"), i), value =0,min=0)
      )})
  })
  # uiOutput("thetaL"),
  output$thetaL <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("thetaL", i), label = paste0(" \\(\\theta_{L}\\)", i), value =13.5,min=0)
      )})
  })
  # uiOutput("alealthetaL"),
  output$alealthetaL <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("alealthetaL", i), label = withMathJax(paste0("\\(\\sigma_{l\\theta_{L}}\\)"), i), value =0,min=0)
      )})
  })
  # uiOutput("deltaAb"),
  output$deltaAb <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("deltaAb", i), label = paste0(" \\(\\delta_{Ab}\\)", i), value =0.029,step=0.01,min=0)
      )})
  })
  # uiOutput("aleadeltaAb"),
  output$aleadeltaAb <- renderUI({
    lapply(1:input$nbinj, function(i) {
      
      withMathJax(numericInput(paste0("aleadeltaAb", i), label = withMathJax(paste0("\\(\\sigma_{\\delta_{Ab}}\\)"), i), value =0,min=0,step=0.01)
      )})
  })
  
  
  
  #liste de tous les renderUI à créer pour le modèle homeostatique:
  #uioutput("hvaccin),
  output$hvaccin <- renderUI({ lapply(1:input$hnbinj, function(i) {
    column(3,
           withMathJax(radioButtons(paste0("hvaccin", i), label = paste0("type de vaccin", i), c("Ad26"=1,"MVA"=2,"Autre"=3))
           ))})
    
  })
  # uiOutput("htinj"),
  output$htinj <- renderUI({
    nbinj <- as.integer(input$hnbinj)
    if (nbinj>=2){
      lapply(2:nbinj, function(i) {
        
        numericInput(paste0("htinj", i), label = paste0("jour de l'injection", i), value = switch(i,0,56,365),min=0)
      })
    }
  })
  # uiOutput("haleatinj"),
  output$haleatinj <- renderUI({
    nbinj <- as.integer(input$hnbinj)
    if (nbinj>=2){
      lapply(2:nbinj, function(i) {
        
        numericInput(paste0("haleatinj", i), label = paste0("\\(\\sigma\\) (jour injection)", i), value = 0,min=0)
      })
    }
  })
  # uiOutput("hrho"),
  output$hrho <- renderUI({ lapply(1:input$hnbinj, function(i) {
    
    withMathJax(numericInput(paste0("hrho", i), label = paste0(" \\(\\rho\\)", i), value = switch(i,4.4,36.6,167.3),min=0)
    )})
    
  })
  # uiOutput("halealrho"),
  output$halealrho <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("halealrho", i), label = paste0("\\(\\sigma_{l\\rho}\\)", i), value = 0,min=0)
      )})
    
  })
  # uiOutput("hbeta"),
  output$hbeta <- renderUI({ lapply(1:input$hnbinj, function(i) {
    
    withMathJax(numericInput(paste0("hbeta", i), label = paste0(" \\(\\beta\\)", i), value =switch(i,30,40,700),step=10^(-4))
    )})
    
  })
  # uiOutput("haleabeta"),
  output$haleabeta <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("haleabeta", i), label = paste0("\\(\\sigma_{\\beta}\\)", i), value = 0,min=0,step=.1)
      )})
    
  })
  
  # uiOutput("halphaM"),
  output$halphaM <- renderUI({ lapply(1:input$hnbinj, function(i) {
    
    withMathJax(numericInput(paste0("halphaM", i), label = paste0(" \\(\\alpha_M\\)", i), value = 5,min=0)
    )})
    
  })
  # uiOutput("haleaalphaM"),
  output$haleaalphaM <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("haleaalphaM", i), label = paste0("\\(\\sigma_{\\alpha_M}\\)", i), value = 0,min=0,step=0.1)
      )})
    
  })
  
  # uiOutput("hdeltaA"),
  output$hdeltaA <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("hdeltaA", i), label = paste0(" \\(\\delta_{A}\\)", i), value = switch(i,0.064,0.21,0.064),min=0,step=0.05)
      )})
    
  })
  # uiOutput("haleadeltaA"),
  output$haleadeltaA <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("haleadeltaA", i), label = paste0("\\(\\sigma_{\\delta_{A}}\\)", i), value = 0,min=0,step=0.1)
      )})
  })
  # uiOutput("hmuS"),
  output$hmuS <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("hmuS", i), label = paste0(" \\(\\mu_{S}\\)", i), value =switch(i,0.13,1.28,0.25),min=0,step=0.1)
      )})
    
  })
  # uiOutput("halealmuS"),
  output$halealmuS <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("halealmuS", i), label = paste0("\\(\\sigma_{l\\mu_{S}}\\)", i), value =0,min=0,step=0.1)
      )})
  })
  
  # uiOutput("halphaS"),
  output$halphaS <- renderUI({ lapply(1:input$hnbinj, function(i) {
    
    withMathJax(numericInput(paste0("halphaS", i), label = paste0(" \\(\\alpha_S\\)", i), value = 10^(-4),min=0,step=10^(-4))
    )})
    
  })
  # uiOutput("haleaalphaS"),
  output$haleaalphaS <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("haleaalphaS", i), label = paste0("\\(\\sigma_{\\alpha_S}\\)", i), value = 0,min=0,step=10^(-4))
      )})
    
  })
  
  # uiOutput("hmuL"),
  output$hmuL <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("hmuL", i), label = paste0(" \\(\\mu_{L}\\)", i), value =14*10^(-4),min=0,step=10^(-4))
      )})
  })
  # uiOutput("halealmuL"),
  output$halealmuL <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("halealmuL", i), label = paste0("\\(\\sigma_{l\\mu_{L}}\\)", i), value =0,min=0,step=10^(-4))
      )})
  })
  
  # uiOutput("halphaL"),
  output$halphaL <- renderUI({ lapply(1:input$hnbinj, function(i) {
    
    withMathJax(numericInput(paste0("halphaL", i), label = paste0(" \\(\\alpha_L\\)", i), value = 10^(-4),min=0,step=10^(-4))
    )})
    
  })
  # uiOutput("haleaalphaL"),
  output$haleaalphaL <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("haleaalphaL", i), label = paste0("\\(\\sigma_{\\alpha_L}\\)", i), value = 0,min=0,step=10^(-4))
      )})
    
  })
  # uiOutput("hdeltaM"),
  output$hdeltaM <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("hdeltaM", i), label = paste0("\\(\\delta_{M}\\)", i), value =1/(60*365.25),min=0,step=10^(-5))
      )})
  })
  # uiOutput("haleadeltaM"),
  output$haleadeltaM <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("haleadeltaM", i), label = withMathJax(paste0("\\(\\sigma_{\\delta_{M}}\\)"), i), value =0,min=0,step=10^(-5))
      )})
  })
  # uiOutput("hdeltaS"),
  output$hdeltaS <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("hdeltaS", i), label = paste0(" \\(\\delta_{S}\\)", i), value =switch(i,0.34,0.23,0.34),min=0,step=.1)
      )})
  })
  # uiOutput("haleadeltaS"),
  output$haleadeltaS <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("haleadeltaS", i), label = withMathJax(paste0("\\(\\sigma_{\\delta_{S}}\\)"), i), value =0,min=0,step=0.1)
      )})
  })
  # uiOutput("hdeltaL"),
  output$hdeltaL <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("hdeltaL", i), label = paste0(" \\(\\delta_{L}\\)", i), value =1/(8.5*365.25),min=0,step=10^(-4))
      )})
  })
  # uiOutput("haleadeltaL"),
  output$haleadeltaL <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("haleadeltaL", i), label = withMathJax(paste0("\\(\\sigma_{\\delta_{L}}\\)"), i), value =0,min=0,step=10^(-4))
      )})
  })
  # uiOutput("hthetaS"),
  output$hthetaS <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("hthetaS", i), label = paste0(" \\(\\theta_{S}\\)", i), value =13.5,min=0)
      )})
  })
  # uiOutput("halealthetaS"),
  output$halealthetaS <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("halealthetaS", i), label = withMathJax(paste0("\\(\\sigma_{l\\theta_{S}}\\)"), i), value =0,min=0)
      )})
  })
  # uiOutput("hthetaL"),
  output$hthetaL <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("hthetaL", i), label = paste0(" \\(\\theta_{L}\\)", i), value =13.5,min=0)
      )})
  })
  # uiOutput("halealthetaL"),
  output$halealthetaL <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("halealthetaL", i), label = withMathJax(paste0("\\(\\sigma_{l\\theta_{L}}\\)"), i), value =0,min=0)
      )})
  })
  # uiOutput("hdeltaAb"),
  output$hdeltaAb <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("hdeltaAb", i), label = paste0("\\(\\delta_{Ab}\\)", i), value =0.029,min=0,step=0.1)
      )})
  })
  # uiOutput("haleadeltaAb"),
  output$haleadeltaAb <- renderUI({
    lapply(1:input$hnbinj, function(i) {
      
      withMathJax(numericInput(paste0("haleadeltaAb", i), label = withMathJax(paste0("\\(\\sigma_{\\delta_{Ab}}\\)"), i), value =0,min=0,step=0.1)
      )})
  })
  #####################################################################################################"
  # #coin des tests pour le modèle avec Mdelay
  #  output$text <- renderText({param()[1]+param()[12]})
  #  
  # output$repa <- renderText({repar(param())
  # 
  # })
  # 
  # 
  # #test 
  #   output$var <-renderText({var <- as.numeric(input$alealgamma)+as.numeric(input$alealrho)+as.numeric(input$aleadeltaA1)+as.numeric(input$alealmuS1)+as.numeric(input$alealmuL1)+as.numeric(input$aleadeltaM1)+as.numeric(input$aleadeltaS1)+as.numeric(input$aleadeltaL1)+as.numeric(input$alealthetaS1)+as.numeric(input$alealthetaL1)+as.numeric(input$aleadeltaAb1)
  # if (input$nbinj>=2){
  #   for (i in 2:input$nbinj){
  #     #checker si tout le monde existe en dessous et penser à copier coller
  #     var <- var+as.numeric(input[[paste0("aleatinj",i)]]) + as.numeric(input[[paste0("aleafrho",i)]])+as.numeric(input[[paste0("aleadeltaA",i)]])+as.numeric(input[[paste0("aleafmuS",i)]])+as.numeric(input[[paste0("alealmuL",i)]])+as.numeric(input[[paste0("aleadeltaM",i)]])+as.numeric(input[[paste0("aleadeltaS",i)]])+as.numeric(input[[paste0("aleadeltaL",i)]])+as.numeric(input[[paste0("alealthetaS",i)]])+as.numeric(input[[paste0("alealthetaL",i)]])+as.numeric(input[[paste0("aleadeltaAb",i)]])
  #   }
  # }
  #  var})
  # #test de mon vecteur param: 
  # output$vec<- renderText({nbinj<-input$nbinj
  #   par<- c(as.numeric(nbinj),as.numeric(input[[paste0("rho",1)]]),as.numeric(input[[paste0("alealrho",1)]]),as.numeric(input[[paste0("deltaA",1)]]),as.numeric(input[[paste0("aleadeltaA",1)]]),as.numeric(input[[paste0("gamma",1)]]),as.numeric(input[[paste0("alealgamma",1)]]),as.numeric(input[[paste0("muS",1)]]),as.numeric(input[[paste0("alealmuS",1)]]),as.numeric(input[[paste0("muL",1)]]),as.numeric(input[[paste0("alealmuL",1)]]),as.numeric(input[[paste0("deltaM",1)]]),as.numeric(input[[paste0("aleadeltaM",1)]]),as.numeric(input[[paste0("deltaS",1)]]),as.numeric(input[[paste0("aleadeltaS",1)]]),as.numeric(input[[paste0("deltaL",1)]]),as.numeric(input[[paste0("aleadeltaL",1)]]),as.numeric(input[[paste0("thetaS",1)]]),as.numeric(input[[paste0("alealthetaS",1)]]),as.numeric(input[[paste0("thetaL",1)]]),as.numeric(input[[paste0("alealthetaL",1)]]),as.numeric(input[[paste0("deltaAb",1)]]),as.numeric(input[[paste0("aleadeltaAb",1)]]))
  # if (nbinj>=2){
  #   for (i in 2:nbinj){
  #     inji<-c(as.numeric(input[[paste0("tinj",i)]]),as.numeric(input[[paste0("aleatinj",i)]]),as.numeric(input[[paste0("rho",i)]]),as.numeric(input[[paste0("alealrho",i)]]),as.numeric(input[[paste0("deltaA",i)]]),as.numeric(input[[paste0("aleadeltaA",i)]]),as.numeric(input[[paste0("gamma",i)]]),as.numeric(input[[paste0("alealgamma",i)]]),as.numeric(input[[paste0("muS",i)]]),as.numeric(input[[paste0("alealmuS",i)]]),as.numeric(input[[paste0("muL",1)]]),as.numeric(input[[paste0("alealmuL",1)]]),as.numeric(input[[paste0("deltaM",1)]]),as.numeric(input[[paste0("aleadeltaM",1)]]),as.numeric(input[[paste0("deltaS",i)]]),as.numeric(input[[paste0("aleadeltaS",i)]]),as.numeric(input[[paste0("deltaL",1)]]),as.numeric(input[[paste0("aleadeltaL",1)]]),as.numeric(input[[paste0("thetaS",1)]]),as.numeric(input[[paste0("alealthetaS",1)]]),as.numeric(input[[paste0("thetaL",1)]]),as.numeric(input[[paste0("alealthetaL",1)]]),as.numeric(input[[paste0("deltaAb",1)]]),as.numeric(input[[paste0("aleadeltaAb",1)]]))
  #     par<-c(par,inji)
  #   }}
  # par
  # })

  #####################################################################################################################
  # Coin des tests pour le modèle avec homéostasie
  # output$htext <- renderText({as.numeric(input$halealrho1)})
  # 
  # 
  # #test pour voir si var est bien calculé 
  # output$hvar <-renderText({
  #   nbinj<-as.numeric(input$hnbinj)
  #   var <- as.numeric(input[[paste0("halealrho",1)]])+as.numeric(input[[paste0("haleabeta",1)]])+as.numeric(input[[paste0("haleaalphaM",1)]])+as.numeric(input[[paste0("haleadeltaA",1)]])+as.numeric(input[[paste0("halealmuS",1)]])+as.numeric(input[[paste0("haleaalphaS",1)]])+as.numeric(input[[paste0("halealmuL",1)]])+as.numeric(input[[paste0("haleaalphaL",1)]])+as.numeric(input[[paste0("haleadeltaM",1)]])+as.numeric(input[[paste0("haleadeltaS",1)]])+as.numeric(input[[paste0("haleadeltaL",1)]])+as.numeric(input[[paste0("halealthetaS",1)]])+as.numeric(input[[paste0("halealthetaL",1)]])+as.numeric(input[[paste0("haleadeltaAb",1)]])
  #   if (nbinj>=2){
  #     for (i in 2:nbinj){
  #       var <- var+as.numeric(input[[paste0("aleatinj",i)]])+ as.numeric(input[[paste0("alealrho",i)]])+as.numeric(input[[paste0("aleabeta",i)]])+as.numeric(input[[paste0("aleaalphaM",i)]])+as.numeric(input[[paste0("aleadeltaA",i)]])+as.numeric(input[[paste0("alealmuS",i)]])+as.numeric(input[[paste0("aleaalphaS",i)]])+as.numeric(input[[paste0("alealmuL",i)]])+as.numeric(input[[paste0("aleaalphaL",i)]])+as.numeric(input[[paste0("aleadeltaM",i)]])+as.numeric(input[[paste0("aleadeltaS",i)]])+as.numeric(input[[paste0("aleadeltaL",i)]])+as.numeric(input[[paste0("alealthetaS",i)]])+as.numeric(input[[paste0("alealthetaL",i)]])+as.numeric(input[[paste0("aleadeltaAb",i)]])
  #     }
  #   }
  #   var})
  # # #test de mon vecteur 
  # output$hvec<- renderText({vecteur<- c(as.numeric(input$hnbinj),as.numeric(input[[paste0("hrho",1)]]),as.numeric(input[[paste0("halealrho",1)]]),as.numeric(input[[paste0("hbeta",1)]]),as.numeric(input[[paste0("haleabeta",1)]]),as.numeric(input[[paste0("halphaM",1)]]),as.numeric(input[[paste0("haleaalphaM",1)]]),as.numeric(input[[paste0("hdeltaA",1)]]),as.numeric(input[[paste0("haleadeltaA",1)]]),as.numeric(input[[paste0("hmuS",1)]]),as.numeric(input[[paste0("halealmuS",1)]]),as.numeric(input[[paste0("halphaS",1)]]),as.numeric(input[[paste0("haleaalphaS",1)]]),as.numeric(input[[paste0("hmuL",1)]]),as.numeric(input[[paste0("halealmuL",1)]]),as.numeric(input[[paste0("halphaL",1)]]),as.numeric(input[[paste0("haleaalphaL",1)]]),as.numeric(input[[paste0("hdeltaM",1)]]),as.numeric(input[[paste0("haleadeltaM",1)]]),as.numeric(input[[paste0("hdeltaS",1)]]),as.numeric(input[[paste0("haleadeltaS",1)]]),as.numeric(input[[paste0("hdeltaL",1)]]),as.numeric(input[[paste0("haleadeltaL",1)]]),as.numeric(input[[paste0("hthetaS",1)]]),as.numeric(input[[paste0("halealthetaS",1)]]),as.numeric(input[[paste0("hthetaL",1)]]),as.numeric(input[[paste0("halealthetaL",1)]]),as.numeric(input[[paste0("hdeltaAb",1)]]),as.numeric(input[[paste0("haleadeltaAb",1)]]))
  # if (input$nbinj>=2){
  #   for (i in 2:input$nbinj){
  #     inji<-c(as.numeric(input[[paste0("htinj",i)]]),as.numeric(input[[paste0("haleatinj",i)]]),as.numeric(input[[paste0("hrho",i)]]),as.numeric(input[[paste0("halealrho",i)]]),as.numeric(input[[paste0("hbeta",i)]]),as.numeric(input[[paste0("haleabeta",i)]]),as.numeric(input[[paste0("halphaM",i)]]),as.numeric(input[[paste0("haleaalphaM",i)]]),as.numeric(input[[paste0("hdeltaA",i)]]),as.numeric(input[[paste0("haleadeltaA",i)]]),as.numeric(input[[paste0("hmuS",i)]]),as.numeric(input[[paste0("halealmuS",i)]]),as.numeric(input[[paste0("halphaS",i)]]),as.numeric(input[[paste0("haleaalphaS",i)]]),as.numeric(input[[paste0("hmuL",i)]]),as.numeric(input[[paste0("halealmuL",i)]]),as.numeric(input[[paste0("halphaL",i)]]),as.numeric(input[[paste0("haleaalphaL",i)]]),as.numeric(input[[paste0("hdeltaM",i)]]),as.numeric(input[[paste0("haleadeltaM",i)]]),as.numeric(input[[paste0("hdeltaS",i)]]),as.numeric(input[[paste0("haleadeltaS",i)]]),as.numeric(input[[paste0("hdeltaL",i)]]),as.numeric(input[[paste0("haleadeltaL",i)]]),as.numeric(input[[paste0("hthetaS",i)]]),as.numeric(input[[paste0("halealthetaS",i)]]),as.numeric(input[[paste0("hthetaL",i)]]),as.numeric(input[[paste0("halealthetaL",i)]]),as.numeric(input[[paste0("hdeltaAb",i)]]),as.numeric(input[[paste0("haleadeltaAb",i)]]))
  #     vecteur<-c(vecteur,inji)
  #   }}
  # 
  # vecteur})
  
  
  #################################################################
  #chez Mdelay
  
  #création du tableau gardant en mémoire les valeurs des paramètres pour les simulations enregistrées.
  #création du data frame gardant les modélisations enregistrées.
  param <- reactiveVal(param)
  ordre <- reactiveVal(ordre)
  parametres <- reactiveVal(parametres)
  simulation_MC<-reactiveVal(simulation_MC)
  simulation_brute <-reactiveVal(simulation_brute)
  simul_brute <- reactiveVal(simul_brute)
  
  #########################################################################
  #chez homéostasie
  #création du tableau gardant en mémoire les valeurs des paramètres pour les simulations enregistrées.
  #création du data frame gardant les simulations enregistrées.
  hparam<- reactiveVal(hparam)
  hparametres<- reactiveVal(hparametres)
  hordre <- reactiveVal(hordre)
  hsimulation_MC <- reactiveVal(hsimulation_MC)
  hsimulation_brute <-reactiveVal(hsimulation_brute)
  hsimul_brute <- reactiveVal(hsimul_brute)
  
  
  ################################################################################
  #Chez MDelay:
  
  #Code déclenché par le bouton go : 
  #checke et corrige le calendrier d'injection
  #-stocke les valeurs de paremtres dans param, regarde si simulation déjà faite et si non :
  #-teste si les ecarts-type sont à 0
  #-selon si variabilité, fait 1 ou 100 ode et Monte carlo
  #crée le data.frame simul_MC contenant les valeurs moyennes et bornes de l'IC95% de chaque compartiment
  #et le data fram simul_brut contenant toutes les simulations faites avec ces valeurs de paramètres (1 si pas de varia et 100 si varia)

  simul_MC<-eventReactive(input$action, 
  { nbinj<-as.numeric(input$nbinj)
    
    #création du vecteur des paramètres à enregistrer et de la chronologie de la vaccination
    #création du vecteur des paramètres à enregistrer et de la chronologie des injections
    par<- c(as.numeric(nbinj),as.numeric(input[[paste0("rho",1)]]),as.numeric(input[[paste0("alealrho",1)]]),as.numeric(input[[paste0("deltaA",1)]]),as.numeric(input[[paste0("aleadeltaA",1)]]),as.numeric(input[[paste0("gamma",1)]]),as.numeric(input[[paste0("alealgamma",1)]]),as.numeric(input[[paste0("muS",1)]]),as.numeric(input[[paste0("alealmuS",1)]]),as.numeric(input[[paste0("muL",1)]]),as.numeric(input[[paste0("alealmuL",1)]]),as.numeric(input[[paste0("deltaM",1)]]),as.numeric(input[[paste0("aleadeltaM",1)]]),as.numeric(input[[paste0("deltaS",1)]]),as.numeric(input[[paste0("aleadeltaS",1)]]),as.numeric(input[[paste0("deltaL",1)]]),as.numeric(input[[paste0("aleadeltaL",1)]]),as.numeric(input[[paste0("thetaS",1)]]),as.numeric(input[[paste0("alealthetaS",1)]]),as.numeric(input[[paste0("thetaL",1)]]),as.numeric(input[[paste0("alealthetaL",1)]]),as.numeric(input[[paste0("deltaAb",1)]]),as.numeric(input[[paste0("aleadeltaAb",1)]]))
    ord<-c()
    if (nbinj>=2){
      for (i in 2:nbinj){
        ord<-c(ord,input[[paste0("tinj",i)]])
        inji<-c(as.numeric(input[[paste0("tinj",i)]]),as.numeric(input[[paste0("aleatinj",i)]]),as.numeric(input[[paste0("rho",i)]]),as.numeric(input[[paste0("alealrho",i)]]),as.numeric(input[[paste0("deltaA",i)]]),as.numeric(input[[paste0("aleadeltaA",i)]]),as.numeric(input[[paste0("gamma",i)]]),as.numeric(input[[paste0("alealgamma",i)]]),as.numeric(input[[paste0("muS",i)]]),as.numeric(input[[paste0("alealmuS",i)]]),as.numeric(input[[paste0("muL",1)]]),as.numeric(input[[paste0("alealmuL",1)]]),as.numeric(input[[paste0("deltaM",1)]]),as.numeric(input[[paste0("aleadeltaM",1)]]),as.numeric(input[[paste0("deltaS",i)]]),as.numeric(input[[paste0("aleadeltaS",i)]]),as.numeric(input[[paste0("deltaL",1)]]),as.numeric(input[[paste0("aleadeltaL",1)]]),as.numeric(input[[paste0("thetaS",1)]]),as.numeric(input[[paste0("alealthetaS",1)]]),as.numeric(input[[paste0("thetaL",1)]]),as.numeric(input[[paste0("alealthetaL",1)]]),as.numeric(input[[paste0("deltaAb",1)]]),as.numeric(input[[paste0("aleadeltaAb",1)]]))
        par<-c(par,inji)
      }}
    param(par)
    ordre(ord)
    #on teste si les doses sont dans un ordre croissant strict
    if(any(diff(ordre()) <= 0)){
      #si le calendrier d'injections n'est pas chronologique, on le rend chronologique
      ord<-ordre()
      ord <- ord[order(ord)]
      ordre(ord)
      #la on regarde si 2 injections sont prévues le même jour, on fait un while
      while (anyDuplicated(ordre())>0){
        #on récupère les valeurs en double
        double<-ordre()[duplicated(ordre())]
        #on récupère les rangs ou apparaissent la première valeur dupliquée
        rangs<-which(ordre() %in% double[1])
        #et on rajoute +1 à sa deuxième apparition
        nouveau<-ordre()
        nouveau[rangs[2]]<-nouveau[rangs[2]]+1
        ordre(nouveau)
        #et on recommence si besoin
      }
      showNotification("Le calendrier vaccinal a été modifié car la chronologie présentait une anomalie")
      #on modifie param() pour prendre les valeurs de ordre comme jours d'injections.
      nouveauparam<-param()
      for (i in 1:length(ordre())){
        nouveauparam[i*24]<-ordre()[i]
      }
      param(nouveauparam)
    }
    
      #test pour voir si param() déjà enregistré dans parametres()
      #premier test : on voit s'il y a des paramètres d'enregistrés et si le nombre de paramètres est au moins égal à celui de la nouvelle simulation
      if (nrow(parametres())!=0 & length(param())<=length(parametres())-2){
        #on ne compare que pour le même nombre d'injections
        comparaison <-subset(parametres(),nbinj==nbinj)
        #on garde le nombre de colonnes qui correspond au nombre d'inj ( enleve les colonnes vides)
        comparaison <- comparaison[,1:(2+23+24*(nbinj-1))]
        double<-c(input$times,param())
        #on trouve les numéros des rangs identiques
        rang<-Reduce(intersect,lapply(1:min(ncol(comparaison)-2,ncol(double)),identique,double,comparaison))
        if(length(rang)!=0){
          showNotification("Cette simulation est déjà enregistrée")
          simul_MC <- subset(simulation_MC(),nom==comparaison$nom[min(rang)])
          brut <- subset(simulation_brute(),nom==comparaison$nom[min(rang)])
          simul_brute(brut)
        }
        else {
          #on fait marcher le pgm
          var <- as.numeric(input[[paste0("alealrho",1)]])+as.numeric(input[[paste0("aleadeltaA",1)]])+as.numeric(input[[paste0("alealgamma",1)]])+as.numeric(input[[paste0("alealmuS",1)]])+as.numeric(input[[paste0("alealmuL",1)]])+as.numeric(input[[paste0("aleadeltaM",1)]])+as.numeric(input[[paste0("aleadeltaS",1)]])+as.numeric(input[[paste0("aleadeltaL",1)]])+as.numeric(input[[paste0("alealthetaS",1)]])+as.numeric(input[[paste0("alealthetaL",1)]])+as.numeric(input[[paste0("aleadeltaAb",1)]])
          if (nbinj>=2){
            for (i in 2:nbinj){
              var <- var+as.numeric(input[[paste0("aleatinj",i)]]) + as.numeric(input[[paste0("alealrho",i)]])+as.numeric(input[[paste0("aleadeltaA",i)]])+as.numeric(input[[paste0("alealgamma",i)]])+as.numeric(input[[paste0("alealmuS",i)]])+as.numeric(input[[paste0("alealmuL",i)]])+as.numeric(input[[paste0("aleadeltaM",i)]])+as.numeric(input[[paste0("aleadeltaS",i)]])+as.numeric(input[[paste0("aleadeltaL",i)]])+as.numeric(input[[paste0("alealthetaS",i)]])+as.numeric(input[[paste0("alealthetaL",i)]])+as.numeric(input[[paste0("aleadeltaAb",i)]])
            }
          }
          #si absence de variabilité, une seule ode
          if (var==0){
            
            simul_MC<-resolve.ode(1,param=param(),times=seq(from=0,to=input$times),xstart=xstart)
            simul_brute(simul_MC)
            simul_MC$Mdelaysup<-simul_MC$Mdelay
            simul_MC$Msup<-simul_MC$M
            simul_MC$Ssup<-simul_MC$S
            simul_MC$Lsup<-simul_MC$L
            simul_MC$Absup<-simul_MC$Ab
            simul_MC$Mdelayinf<-simul_MC$Mdelay
            simul_MC$Minf<-simul_MC$M
            simul_MC$Sinf<-simul_MC$S
            simul_MC$Linf<-simul_MC$L
            simul_MC$Abinf<-simul_MC$Ab
            
          }
          else {
            #si variabilité, cela va prendre plus de temps donc on fait marcher le petit truc de chargement
            waiter <- waiter::Waiter$new()
            waiter$show()
            on.exit(waiter$hide())
            
            
            #on fait les 100 modélisations
            recaps<-lapply(seq_len(100), resolve.ode,param(),seq(from=0,to=input$times),xstart)
            #on les met ensemble
            brut<-data.table::rbindlist(recaps)
            
            #calcul des moyennes et quantiles
            brut<-as.data.table(brut)
            moy<-brut[,lapply(.SD,mean),by=time]
            inf <-brut[,lapply(.SD,quantile,probs=.025),by=time]
            sup<-brut[,lapply(.SD,quantile,probs=.975),by=time]
            colnames(inf)<-c("time","Mdelayinf","Minf","Sinf","Linf","Abinf")
            colnames (sup)<-c("time","Mdelaysup","Msup","Ssup","Lsup","Absup")
            simul_MC<-data.table::merge.data.table(moy,sup,by="time")
            simul_MC<-data.table::merge.data.table(simul_MC,inf,by="time")
            simul_brute(brut)
          }
          p<-c(input$times,param())
          param(p)
          simul_MC<-as.data.frame(simul_MC)
        }
      }
      else{
        #sinon
        #somme de tous les écarts-type :
        var <- as.numeric(input[[paste0("alealrho",1)]])+as.numeric(input[[paste0("aleadeltaA",1)]])+as.numeric(input[[paste0("alealgamma",1)]])+as.numeric(input[[paste0("alealmuS",1)]])+as.numeric(input[[paste0("alealmuL",1)]])+as.numeric(input[[paste0("aleadeltaM",1)]])+as.numeric(input[[paste0("aleadeltaS",1)]])+as.numeric(input[[paste0("aleadeltaL",1)]])+as.numeric(input[[paste0("alealthetaS",1)]])+as.numeric(input[[paste0("alealthetaL",1)]])+as.numeric(input[[paste0("aleadeltaAb",1)]])
        if (nbinj>=2){
          for (i in 2:nbinj){
            var <- var+as.numeric(input[[paste0("aleatinj",i)]]) + as.numeric(input[[paste0("alealrho",i)]])+as.numeric(input[[paste0("aleadeltaA",i)]])+as.numeric(input[[paste0("alealgamma",i)]])+as.numeric(input[[paste0("alealmuS",i)]])+as.numeric(input[[paste0("alealmuL",i)]])+as.numeric(input[[paste0("aleadeltaM",i)]])+as.numeric(input[[paste0("aleadeltaS",i)]])+as.numeric(input[[paste0("aleadeltaL",i)]])+as.numeric(input[[paste0("alealthetaS",i)]])+as.numeric(input[[paste0("alealthetaL",i)]])+as.numeric(input[[paste0("aleadeltaAb",i)]])
          }
        }
        #si absence de variabilité, une seule ode
        if (var==0){
          
          simul_MC<-resolve.ode(1,param=param(),times=seq(from=0,to=input$times),xstart=xstart)
          simul_brute(simul_MC)
          simul_MC$Mdelaysup<-simul_MC$Mdelay
          simul_MC$Msup<-simul_MC$M
          simul_MC$Ssup<-simul_MC$S
          simul_MC$Lsup<-simul_MC$L
          simul_MC$Absup<-simul_MC$Ab
          simul_MC$Mdelayinf<-simul_MC$Mdelay
          simul_MC$Minf<-simul_MC$M
          simul_MC$Sinf<-simul_MC$S
          simul_MC$Linf<-simul_MC$L
          simul_MC$Abinf<-simul_MC$Ab
          
        }
        else {
          #si variabilité, cela va prendre plus de temps donc on fait marcher le petit truc de chargement
          waiter <- waiter::Waiter$new()
          waiter$show()
          on.exit(waiter$hide())
          
          
          #on fait les 100 modélisations
          recaps<-lapply(seq_len(100), resolve.ode,param(),seq(from=0,to=input$times),xstart)
          #on les met ensemble
          brut<-data.table::rbindlist(recaps)
          
          #calcul des moyennes et quantiles
          brut<-as.data.table(brut)
          moy<-brut[,lapply(.SD,mean),by=time]
          inf <-brut[,lapply(.SD,quantile,probs=.025),by=time]
          sup<-brut[,lapply(.SD,quantile,probs=.975),by=time]
          colnames(inf)<-c("time","Mdelayinf","Minf","Sinf","Linf","Abinf")
          colnames (sup)<-c("time","Mdelaysup","Msup","Ssup","Lsup","Absup")
          simul_MC<-data.table::merge.data.table(moy,sup,by="time")
          simul_MC<-data.table::merge.data.table(simul_MC,inf,by="time")
          simul_brute(brut)
        }
        p<-c(input$times,param())
        param(p)
        simul_MC<-as.data.frame(simul_MC)
      }
    })
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #construction des ggplot
  #on fait en 2 morceaux pour pouvoir télécharger :
  grapheM<-reactive({ggplot(simul_MC())+
      geom_line(aes(x=time,y=M),col='blue',size=0.4)+
      geom_ribbon(aes(ymin=Minf,ymax=Msup, x= time), fill="skyblue3",alpha=0.3)+
      theme_classic()+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      #scale_y_continuous(trans="log10")+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      labs(x='time post vaccination(days)',y="Memory B-cells concentration (ASCS/millions)")+
      theme(legend.position = "none") 
      #ggtitle("graphe de M")
      })
  
  output$grapheM <- renderPlot({grapheM()
  })
  
  grapheS <- reactive({ggplot(simul_MC())+
      geom_line(aes(x=time,y=S),col='blue',size=0.4)+
      geom_ribbon(aes(ymin=Sinf,ymax=Ssup, x= time), fill="skyblue3",alpha=0.3)+
      #scale_y_continuous(trans="log10")+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      
      theme_classic()+
      labs(x='time post vaccination(days)',y="ASCS concentration (ASCS/millions)")+
      theme(legend.position = "none") 
      #ggtitle("graphe de S")
      })
  
  output$grapheS <- renderPlot ({grapheS()
  })
  
  grapheL<-reactive({ggplot(simul_MC())+
      geom_line(aes(x=time,y=L),col='blue',size=0.4)+
      geom_ribbon(aes(ymin=Linf,ymax=Lsup, x= time), fill="skyblue3",alpha=0.3)+
      theme_classic()+
      #scale_y_continuous(trans="log10")+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      labs(x='time post vaccination(days)',y="ASCS concentration (ASCS/millions)")+
      theme(legend.position = "none") 
      #ggtitle("graphe de L")
      })
  output$grapheL <- renderPlot ({
    grapheL()})
  
  
  grapheAb<-reactive({ggplot(simul_MC())+
      geom_line(aes(x=time,y=Ab),col='blue',size=0.4)+
      geom_ribbon(aes(ymin=Abinf,ymax=Absup, x= time), fill="skyblue3",alpha=0.3)+
      theme_classic()+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      #scale_y_continuous(trans="log10")+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      labs(x='time post vaccination(days)',y="Ab concentration (ELISA units)")+
      theme(legend.position = "none") 
      #ggtitle("graphe de Ab")
      })
  
  output$grapheAb <- renderPlot({grapheAb()
  })
  
  output$allgraphs = downloadHandler(
    filename = function() {
      'all_images.zip'
    }, 
    content = function(fname) {
      fs <- replicate(4, tempfile(fileext = ".png"))
      ggsave(fs[1], grapheM())
      ggsave(fs[2], grapheS())
      ggsave(fs[3], grapheL())
      ggsave(fs[4], grapheAb())
      zip::zipr(zipfile=fname, files=fs)
    },
    contentType = "application/zip")
  
  
  
  #ajout d'un enregistrement
  
  observeEvent(input$keep,{
    nom<-ifelse(input$nommodel=="",paste0("simulation ",input$keep+input$hkeep),input$nommodel)
    nbinj<-param()[1]
    if (nrow(parametres())==0){
      p <- c (nom,param())
      param(p)
      t <- rbind (parametres(),param())
      parametres(t)
      #On ajoute les valeurs de chaque compartiment de la simulation dans le tableau simulation
      simul<-as.data.frame(simul_MC())
      simul$nom<-nom
      simul$modele<- paste("modèle avec Mdelay",nom,sep=" ")
      j <-rbind(simulation_MC(),simul)
      simulation_MC(j)
      brut <- as.data.table(simul_brute())
      brut$nom<-nom
      brut$modele<- paste("modèle avec Mdelay",nom,sep=" ")
      r <- rbind (simulation_brute(),brut)
      simulation_brute(r)
    }
    
    else {
      #test pour voir si le nom existe déjà :
      if(nom %in% parametres()[,1]){
        showNotification("Vous avez déja utlisé ce nom")
      }
      else{
        #ajoute des NA au nouvel enregistrement si nb injection inférieur pour avoir les même tailles de lignes
        if(length(param())<length(parametres())-1){
          n <- (length(parametres())-1) - length(param())
          a <- rep(NA,n)
          par<-c(param(),a)
          param(par)
        }
        
        if (length(param())==length(parametres())-1){
          #nom<-ifelse(input$nommodel=="",paste0("simulation ",input$keep+input$hkeep),input$nommodel)
          p<-c(nom,param())
          t<-rbind(parametres(),p)
          parametres(t)
        }
        if (length(param())>length(parametres())){
          n <- length(param())-(length(parametres())-1)
          t=parametres()
          ajout<-data.frame(matrix(nrow=1,ncol=n))
          t <- cbind (t,ajout)
          nom<-ifelse(input$nommodel=="",paste0("simulation ",input$keep+input$hkeep),input$nommodel)
          p<-c(nom,param())
          t <- rbind (t,p)
          parametres(t)
        }
    
    #On ajoute les valeurs de chaque compartiment de la simulation dans le tableau simulation
        simul<-as.data.frame(simul_MC())
        simul$nom<-nom
        simul$modele<- paste("modèle avec Mdelay",nom,sep=" ")
    j <-rbind(simulation_MC(),simul)
    simulation_MC(j)
    brut <- as.data.table(simul_brute())
    brut$nom<-nom
    brut$modele<- paste("modèle avec Mdelay",nom,sep=" ")
    r <- rbind (simulation_brute(),brut)
    simulation_brute(r)
    }
    }
    #on gère les colnames de tout le monde :
    
    n<-length(parametres())
    
    k<-((n-25)/24)+1
    
    noms_param<-c("nom","temps","nbinj","ρ","δ(lρ)","δA1","σ(δA1)","γ1","σ(lγ1)","μS1","σ(lμS1)","μL1","σ(lμL1)", "δM1","σ(δM1)","δS1","σ(δS1)","δL1","σ(δL1)","ϴS1","σ(lϴS1)","ϴL1","σ(lϴL1)","δAb1","σ(δAb1)")
    if (k>=2){
      for (i in 2:(k)){
        ajout<-c(paste0("tinj",i),paste0("σ(tinj)",i),paste0("ρ",i),paste0("σ(lρ)",i),paste0("δA",i),paste0("σ(δA)",i),paste0("γ",i),paste0("σ(γ)",i),paste0("μS",i),paste0("σ(lμS)",i),paste0("μL",i),paste0("σ(lμL)",i),paste0("δM",i),paste0("σ(δM)",i),paste0("δS",i),paste0("σ(δS)",i),paste0("δL",i),paste0("σ(δL)",i),paste0("ϴS",i),paste0("σ(lϴS)",i),paste0("ϴL",i),paste0("σ(lϴL)",i),paste0("δAb",i),paste0("σδAb",i))
        noms_param<-c(noms_param,ajout)
      }}
    a=parametres()
    b=simulation_MC()
    colnames(a)<-noms_param
    colnames(b) <- c("time","Mdelay","M","S","L","Ab","Mdelayinf","Minf","Sinf","Linf","Abinf","Mdelaysup","Msup","Ssup","Lsup","Absup","nom","modele")
    parametres(a)
    simulation_MC(b)
    
  })
  
  
  
  #tableau avec cases à cocher
  output[["table_param"]] <-renderReactable({reactable(parametres(),
                                                     selection = "multiple", 
                                                     onClick = "select")}
  )
  #suppression des enregistrements :
  observeEvent(input$clear,{
    p<-parametres()[0,]
    parametres(p)
    t<-simulation_MC()[0,]
    simulation_MC(t)
    r<-simulation_brute()[0,]
    simulation_brute(r)
  })
  #récuperer les rangs sélectionnés : getReactableState()
  selected <- reactive(getReactableState("table_param", "selected"))
  #étape d'après: récupérer les noms des modélisation pour subset dans simulation() et plot que ceux là
  noms_simul<-reactive({lapply(selected(),function(i){parametres()[i,1]})})
  

  
  
  #Faire les moyennes et IC des données si les variables sont transformées
  comparaison_M<-reactive({
    if (input$logM==TRUE){
      transfo<-as.data.table(simulation_brute()[,c("time","M","nom")])
      transfo$M<-log10(transfo$M)
      moy<-transfo[,lapply(.SD,mean),by=.(time,nom)]
      inf <-transfo[,lapply(.SD,quantile,probs=.025),by=.(time,nom)]
      colnames(inf)<-c("time","nom","Minf")
      sup<-transfo[,lapply(.SD,quantile,probs=.975),by=.(time,nom)]
      colnames(sup)<-c("time","nom","Msup")
      comparaison_M <- data.table::merge.data.table(moy,sup,by=c("time","nom"))
      comparaison_M<-data.table::merge.data.table(comparaison_M,inf,by=c("time","nom"))
    }
    else{

      comparaison_M <- simulation_MC()
    }
    comparaison_M
  })

  comparaison_S<-reactive({if (input$racS==TRUE){
    transfo<-as.data.table(simulation_brute()[,c("time","S","nom")])
    transfo$S<-(transfo$S)^(1/4)
    moy<-transfo[,lapply(.SD,mean),by=.(time,nom)]
    inf <-transfo[,lapply(.SD,quantile,probs=.025,na.rm=TRUE),by=.(time,nom)]
    colnames(inf)<-c("time","nom","Sinf")
    sup<-transfo[,lapply(.SD,quantile,probs=.975,na.rm=TRUE),by=.(time,nom)]
    colnames(sup)<-c("time","nom","Ssup")
    comparaison_S <- data.table::merge.data.table(moy,sup,by=c("time","nom"))
    comparaison_S<-data.table::merge.data.table(comparaison_S,inf,by=c("time","nom"))
  }
    else{

      comparaison_S <- simulation_MC()
    }
    comparaison_S})

  comparaison_L<-reactive({if (input$racL==TRUE){
    transfo<-as.data.table(simulation_brute()[,c("time","L","nom")])
    transfo$L<-(transfo$L)^(1/4)
    moy<-transfo[,lapply(.SD,mean),by=.(time,nom)]
    inf <-transfo[,lapply(.SD,quantile,probs=.025,na.rm=TRUE),by=.(time,nom)]
    colnames(inf)<-c("time","nom","Linf")
    sup<-transfo[,lapply(.SD,quantile,probs=.975,na.rm=TRUE),by=.(time,nom)]
    colnames(sup)<-c("time","nom","Lsup")
    comparaison_L <- data.table::merge.data.table(moy,sup,by=c("time","nom"))
    comparaison_L<-data.table::merge.data.table(comparaison_L,inf,by=c("time","nom"))
  }
     else{
  
       comparaison_L <- simulation_MC()
     }
     comparaison_L})

   comparaison_Ab<-reactive({if (input$logAb==TRUE){
     transfo<-as.data.table(simulation_brute()[,c("time","Ab","nom")])
     transfo$Ab<-log10(transfo$Ab)
     #comparaison_M<-transfo
     moy<-transfo[,lapply(.SD,mean),by=.(time,nom)]
     inf <-transfo[,lapply(.SD,quantile,probs=.025,na.rm=TRUE),by=.(time,nom)]
     colnames(inf)<-c("time","nom","Abinf")
     sup<-transfo[,lapply(.SD,quantile,probs=.975,na.rm=TRUE),by=.(time,nom)]
     colnames(sup)<-c("time","nom","Absup")
     comparaison_Ab <- data.table::merge.data.table(moy,sup,by=c("time","nom"))
     comparaison_Ab<-data.table::merge.data.table(comparaison_Ab,inf,by=c("time","nom"))
   }
     else{
  
       comparaison_Ab <- simulation_MC()
     }
     comparaison_Ab})
  


  
  
  
  
  
  output$grapheMcompare <- renderPlot({
    partie<-subset(comparaison_M(),nom %in% noms_simul())
    ggplot(partie)+
      geom_line(aes(x=time,y=M,color=nom),size=0.4)+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      #scale_y_continuous(trans="log10")+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      geom_ribbon(aes(ymin=Minf,ymax=Msup, x= time, fill=nom),alpha=0.3)+
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE)+
      theme_classic()+
      labs(x='time post vaccination(days)',y="Memory B-cells concentration (ASCS/millions)")+
      theme(legend.position="top")
      #theme(legend.position = "none") +
      #ggtitle("graphe de M")
      }) 
  
  
  output$grapheScompare <- renderPlot({
    partie<-subset(comparaison_S(),nom %in% noms_simul())
    ggplot(partie)+
      geom_line(aes(x=time,y=S,color=nom),size=0.4)+
      geom_ribbon(aes(ymin=Sinf,ymax=Ssup, x= time, fill=nom),alpha=0.3)+
      theme_classic()+
      #scale_y_continuous(trans="log10")+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      labs(x='time post vaccination(days)',y="ASCS concentration (ASCS/millions)")+
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE)+
      theme(legend.position="top")
      #theme(legend.position = "none") +
      #ggtitle("graphe de S")
      }) 
  
  output$grapheLcompare <- renderPlot({
    partie<-subset(comparaison_L(),nom %in% noms_simul())
    ggplot(partie)+
      geom_line(aes(x=time,y=L,color=nom),size=0.4)+
      geom_ribbon(aes(ymin=Linf,ymax=Lsup, x= time, fill=nom),alpha=0.3)+
      theme_classic()+
      #scale_y_continuous(trans="log10")+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      labs(x='time post vaccination(days)',y="ASCS concentration (ASCS/millions)")+
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE)+
      theme(legend.position="top")
      #theme(legend.position = "none") +
      #ggtitle("graphe de L")
      })
  
  output$grapheAbcompare <- renderPlot({
    partie<-subset(comparaison_Ab(),nom %in% noms_simul())
    ggplot(partie)+
      geom_line(aes(x=time,y=Ab,color=nom),size=0.4)+
      geom_ribbon(aes(ymin=Abinf,ymax=Absup, x= time, fill=nom),alpha=0.3)+
      theme_classic()+
      #scale_y_continuous(trans="log10") + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      scale_x_continuous(expand = c(0, 0), limits = c(0, NA))+
      labs(x='time post vaccination(days)',y="Ab concentration (ELISA units)")+
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE)+
      theme(legend.position="top")
      #theme(legend.position = "none") +
      #ggtitle("graphe de Ab")
      })
  
  
  output$dl <- downloadHandler(
    
    filename = function() {
      "Mdelay.xlsx"
    },
    content = function(filename){
      
      df_list <- list(parametres(), simulation_MC(), simulation_brute())
      write.xlsx(x = df_list , file = filename, rowNames = FALSE)
    }
  ) 
  ######################################################################################
  #Code chez homéostasie
  #Construction des modélisations 
  hsimul_MC<-eventReactive(input$haction,
                         { nbinj<-as.numeric(input$hnbinj)
                         
                         #création du vecteur des paramètres à enregistrer et du calendrier vaccinal 
                         par<- c(as.numeric(nbinj),as.numeric(input[[paste0("hrho",1)]]),as.numeric(input[[paste0("halealrho",1)]]),as.numeric(input[[paste0("hbeta",1)]]),as.numeric(input[[paste0("haleabeta",1)]]),as.numeric(input[[paste0("halphaM",1)]]),as.numeric(input[[paste0("haleaalphaM",1)]]),as.numeric(input[[paste0("hdeltaA",1)]]),as.numeric(input[[paste0("haleadeltaA",1)]]),as.numeric(input[[paste0("hmuS",1)]]),as.numeric(input[[paste0("halealmuS",1)]]),as.numeric(input[[paste0("halphaS",1)]]),as.numeric(input[[paste0("haleaalphaS",1)]]),as.numeric(input[[paste0("hmuL",1)]]),as.numeric(input[[paste0("halealmuL",1)]]),as.numeric(input[[paste0("halphaL",1)]]),as.numeric(input[[paste0("haleaalphaL",1)]]),as.numeric(input[[paste0("hdeltaM",1)]]),as.numeric(input[[paste0("haleadeltaM",1)]]),as.numeric(input[[paste0("hdeltaS",1)]]),as.numeric(input[[paste0("haleadeltaS",1)]]),as.numeric(input[[paste0("hdeltaL",1)]]),as.numeric(input[[paste0("haleadeltaL",1)]]),as.numeric(input[[paste0("hthetaS",1)]]),as.numeric(input[[paste0("halealthetaS",1)]]),as.numeric(input[[paste0("hthetaL",1)]]),as.numeric(input[[paste0("halealthetaL",1)]]),as.numeric(input[[paste0("hdeltaAb",1)]]),as.numeric(input[[paste0("haleadeltaAb",1)]]))
                         hord<-c()
                         if (nbinj>=2){
                           for (i in 2:nbinj){
                             hord<-c(hord,as.numeric(input[[paste0("htinj",i)]]))
                             inji<-c(as.numeric(input[[paste0("htinj",i)]]),as.numeric(input[[paste0("haleatinj",i)]]),as.numeric(input[[paste0("hrho",i)]]),as.numeric(input[[paste0("halealrho",i)]]),as.numeric(input[[paste0("hbeta",i)]]),as.numeric(input[[paste0("haleabeta",i)]]),as.numeric(input[[paste0("halphaM",i)]]),as.numeric(input[[paste0("haleaalphaM",i)]]),as.numeric(input[[paste0("hdeltaA",i)]]),as.numeric(input[[paste0("haleadeltaA",i)]]),as.numeric(input[[paste0("hmuS",i)]]),as.numeric(input[[paste0("halealmuS",i)]]),as.numeric(input[[paste0("halphaS",i)]]),as.numeric(input[[paste0("haleaalphaS",i)]]),as.numeric(input[[paste0("hmuL",i)]]),as.numeric(input[[paste0("halealmuL",i)]]),as.numeric(input[[paste0("halphaL",i)]]),as.numeric(input[[paste0("haleaalphaL",i)]]),as.numeric(input[[paste0("hdeltaM",i)]]),as.numeric(input[[paste0("haleadeltaM",i)]]),as.numeric(input[[paste0("hdeltaS",i)]]),as.numeric(input[[paste0("haleadeltaS",i)]]),as.numeric(input[[paste0("hdeltaL",i)]]),as.numeric(input[[paste0("haleadeltaL",i)]]),as.numeric(input[[paste0("hthetaS",i)]]),as.numeric(input[[paste0("halealthetaS",i)]]),as.numeric(input[[paste0("hthetaL",i)]]),as.numeric(input[[paste0("halealthetaL",i)]]),as.numeric(input[[paste0("hdeltaAb",i)]]),as.numeric(input[[paste0("haleadeltaAb",i)]]))
                             par<-c(par,inji)
                           }}
                         hparam(par)
                         hordre(hord)
                         #on teste si les doses sont dans un ordre croissant strict
                         if(any(diff(hordre()) <= 0)){
                           hord<-hordre()
                           hord <- hord[order(hord)]
                           hordre(hord)
                           #la on regarde si 2 injections sont prévues le même jour, on fait un while
                           while (anyDuplicated(hordre())>0){
                             #on récupère les valeurs en double
                             hdouble<-hordre()[duplicated(hordre())]
                             #on récupère les rangs ou apparaissent la première valeur dupliquée
                             hrangs<-which(hordre() %in% hdouble[1])
                             #et on rajoute +1 à sa deuxième apparition
                             hnouveau<-hordre()
                             hnouveau[rangs[2]]<-hnouveau[hrangs[2]]+1
                             hordre(hnouveau)
                             #et on recommence si besoin
                           }
                           showNotification("Le calendrier vaccinal a été modifié car la chronologie présentait une anomalie")
                           #on modifie param() pour prendre les valeurs de ordre comme jours d'injections.
                           hnouveauparam<-hparam()
                           for (i in 1:length(hordre())){
                             hnouveauparam[30*i]<-hordre()[i]
                           }
                           hparam(hnouveauparam)
                         }
                         #test pour voir si param() déjà enregistré dans parametres()
                         #premier test : on voit s'il y a des paramètres d'enregistrés et si le vecteur de paramètres n'est pas plus grand que le tableau 
                           if (nrow(hparametres())!=0 & length(hparam())<=length(hparametres())-2){
                             #on ne compare que pour le même nombre d'injections
                             comparaison <-subset(hparametres(),nbinj==nbinj)
                             #on garde le nombre de colonnes qui correspond au nombre d'inj ( enleve les colonnes vides)
                             comparaison <- comparaison[,1:(2+29+30*(nbinj-1))]
                             double<-c(input$htimes,hparam())
                             #on trouve les numéros des rangs identiques
                             rang<-Reduce(intersect,lapply(1:min(ncol(comparaison)-2,ncol(double)),identique,double,comparaison))
                             if(length(rang)!=0){
                               showNotification("Cette simulation est déja enregistrée")
                               hsimul_MC <- subset(hsimulation_MC(),nom==comparaison$nom[min(rang)])
                               hbrut <-subset(hsimulation_brute(),nom==comparaison$nom[min(rang)])
                               hsimul_brute(hbrut)
                             }
                         
                           else{
                             #si pas d'enregistrement identique
                             #somme de tous les écarts-type :
                             var <- as.numeric(input[[paste0("halealrho",1)]])+as.numeric(input[[paste0("haleabeta",1)]])+as.numeric(input[[paste0("haleaalphaM",1)]])+as.numeric(input[[paste0("haleadeltaA",1)]])+as.numeric(input[[paste0("halealmuS",1)]])+as.numeric(input[[paste0("haleaalphaS",1)]])+as.numeric(input[[paste0("halealmuL",1)]])+as.numeric(input[[paste0("haleaalphaL",1)]])+as.numeric(input[[paste0("haleadeltaM",1)]])+as.numeric(input[[paste0("haleadeltaS",1)]])+as.numeric(input[[paste0("haleadeltaL",1)]])+as.numeric(input[[paste0("halealthetaS",1)]])+as.numeric(input[[paste0("halealthetaL",1)]])+as.numeric(input[[paste0("haleadeltaAb",1)]])
                             if (nbinj>=2){
                               for (i in 2:nbinj){
                                 var <- var+as.numeric(input[[paste0("haleatinj",i)]])+ as.numeric(input[[paste0("halealrho",i)]])+as.numeric(input[[paste0("haleabeta",i)]])+as.numeric(input[[paste0("haleaalphaM",i)]])+as.numeric(input[[paste0("haleadeltaA",i)]])+as.numeric(input[[paste0("halealmuS",i)]])+as.numeric(input[[paste0("haleaalphaS",i)]])+as.numeric(input[[paste0("halealmuL",i)]])+as.numeric(input[[paste0("haleaalphaL",i)]])+as.numeric(input[[paste0("haleadeltaM",i)]])+as.numeric(input[[paste0("haleadeltaS",i)]])+as.numeric(input[[paste0("haleadeltaL",i)]])+as.numeric(input[[paste0("halealthetaS",i)]])+as.numeric(input[[paste0("halealthetaL",i)]])+as.numeric(input[[paste0("haleadeltaAb",i)]])
                               }
                             }
                             #si absence de variabilité, une seule ode
                             if (var==0){
                               
                               hsimul_MC<-hresolve.ode(1,param=hparam(),times=seq(from=0,to=input$htimes),hxstart=hxstart)
                               hsimul_brute(hsimul_MC)
                               
                               hsimul_MC$Msup<-hsimul_MC$M
                               hsimul_MC$Ssup<-hsimul_MC$S
                               hsimul_MC$Lsup<-hsimul_MC$L
                               hsimul_MC$Absup<-hsimul_MC$Ab
                               hsimul_MC$Minf<-hsimul_MC$M
                               hsimul_MC$Sinf<-hsimul_MC$S
                               hsimul_MC$Linf<-hsimul_MC$L
                               hsimul_MC$Abinf<-hsimul_MC$Ab
                               
                             }
                             else {
                               #si variabilité, cela va prendre plus de temps donc on fait marcher le petit truc de chargement
                               waiter <- waiter::Waiter$new()
                               waiter$show()
                               on.exit(waiter$hide())
                               
                               
                               #on fait les 100 modélisations
                               hrecaps<-lapply(seq_len(100), hresolve.ode,hparam(),seq(from=0,to=input$htimes),hxstart)
                               #on les met ensemble
                               recap<-data.table::rbindlist(hrecaps)
                               
                               #calcul des moyennes et quantiles
                               truc<-as.data.table(recap)
                               moy<-truc[,lapply(.SD,mean),by=time]
                               inf <-truc[,lapply(.SD,quantile,probs=.025),by=time]
                               sup<-truc[,lapply(.SD,quantile,probs=.975),by=time]
                               colnames(inf)<-c("time","Minf","Sinf","Linf","Abinf")
                               colnames (sup)<-c("time","Msup","Ssup","Lsup","Absup")
                               hsimul_MC<-data.table::merge.data.table(moy,sup,by="time")
                               hsimul_MC<-data.table::merge.data.table(hsimul_MC,inf,by="time")
                             }
                           }
                           
                         }
                         #else de si il n'y avait pas d'enregistrement 
                         else {
                           #somme de tous les écarts-type :
                           var <- as.numeric(input[[paste0("halealrho",1)]])+as.numeric(input[[paste0("haleabeta",1)]])+as.numeric(input[[paste0("haleaalphaM",1)]])+as.numeric(input[[paste0("haleadeltaA",1)]])+as.numeric(input[[paste0("halealmuS",1)]])+as.numeric(input[[paste0("haleaalphaS",1)]])+as.numeric(input[[paste0("halealmuL",1)]])+as.numeric(input[[paste0("haleaalphaL",1)]])+as.numeric(input[[paste0("haleadeltaM",1)]])+as.numeric(input[[paste0("haleadeltaS",1)]])+as.numeric(input[[paste0("haleadeltaL",1)]])+as.numeric(input[[paste0("halealthetaS",1)]])+as.numeric(input[[paste0("halealthetaL",1)]])+as.numeric(input[[paste0("haleadeltaAb",1)]])
                           if (nbinj>=2){
                             for (i in 2:nbinj){
                               var <- var+as.numeric(input[[paste0("haleatinj",i)]])+ as.numeric(input[[paste0("halealrho",i)]])+as.numeric(input[[paste0("haleabeta",i)]])+as.numeric(input[[paste0("haleaalphaM",i)]])+as.numeric(input[[paste0("haleadeltaA",i)]])+as.numeric(input[[paste0("halealmuS",i)]])+as.numeric(input[[paste0("haleaalphaS",i)]])+as.numeric(input[[paste0("halealmuL",i)]])+as.numeric(input[[paste0("haleaalphaL",i)]])+as.numeric(input[[paste0("haleadeltaM",i)]])+as.numeric(input[[paste0("haleadeltaS",i)]])+as.numeric(input[[paste0("haleadeltaL",i)]])+as.numeric(input[[paste0("halealthetaS",i)]])+as.numeric(input[[paste0("halealthetaL",i)]])+as.numeric(input[[paste0("haleadeltaAb",i)]])
                             }
                           }
                           #si absence de variabilité, une seule ode
                           if (var==0){
                             
                             hsimul_MC<-hresolve.ode(1,param=hparam(),times=seq(from=0,to=input$htimes),hxstart=hxstart)
                             hsimul_brute(hsimul_MC)
                             
                             hsimul_MC$Msup<-hsimul_MC$M
                             hsimul_MC$Ssup<-hsimul_MC$S
                             hsimul_MC$Lsup<-hsimul_MC$L
                             hsimul_MC$Absup<-hsimul_MC$Ab
                             hsimul_MC$Minf<-hsimul_MC$M
                             hsimul_MC$Sinf<-hsimul_MC$S
                             hsimul_MC$Linf<-hsimul_MC$L
                             hsimul_MC$Abinf<-hsimul_MC$Ab
                             
                           }
                           else {
                             #si variabilité, cela va prendre plus de temps donc on fait marcher le petit truc de chargement
                             waiter <- waiter::Waiter$new()
                             waiter$show()
                             on.exit(waiter$hide())
                             
                             
                             #on fait les 100 modélisations
                             hrecaps<-lapply(seq_len(100), hresolve.ode,hparam(),seq(from=0,to=input$htimes),hxstart)
                             #on les met ensemble
                             recap<-data.table::rbindlist(hrecaps)
                             
                             #calcul des moyennes et quantiles
                             brut<-as.data.table(recap)
                             moy<-brut[,lapply(.SD,mean),by=time]
                             inf <-brut[,lapply(.SD,quantile,probs=.025),by=time]
                             sup<-brut[,lapply(.SD,quantile,probs=.975),by=time]
                             colnames(inf)<-c("time","Minf","Sinf","Linf","Abinf")
                             colnames (sup)<-c("time","Msup","Ssup","Lsup","Absup")
                             hsimul_MC<-data.table::merge.data.table(moy,sup,by="time")
                             hsimul_MC<-data.table::merge.data.table(hsimul_MC,inf,by="time")
                             hsimul_brute(brut)
                           }
                         }
                         p<-c(input$htimes,hparam())
                         hparam(p)
                         hsimul_MC<-as.data.frame(hsimul_MC)

                         })
  
  
  

  
  #construction des ggplot
  #on fait en 2 morceaux pour pouvoir télécharger :
  hgrapheM<-reactive({ggplot(hsimul_MC())+
      geom_line(aes(x=time,y=M),col="blue",size=0.4)+
      geom_ribbon(aes(ymin=Minf,ymax=Msup, x= time), fill="skyblue3",alpha=0.3)+
      theme_classic()+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      #scale_y_continuous(trans="log10")+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      labs(x='time post vaccination(days)',y="Memory B-cells concentration (ASCS/millions)")+
      theme(legend.position = "none") 
      #ggtitle("graphe de M")
      })
  
  output$hgrapheM <- renderPlot({hgrapheM()
  })
  
  hgrapheS<-reactive({ggplot(hsimul_MC())+
      geom_line(aes(x=time,y=S),col="blue",size=0.4)+
      #scale_y_continuous(trans="log10")+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      geom_ribbon(aes(ymin=Sinf,ymax=Ssup, x= time), fill="skyblue3",alpha=0.3)+
      theme_classic()+
      labs(x='time post vaccination(days)',y="ASCS concentration (ASCS/millions)")+
      theme(legend.position = "none") 
      #ggtitle("graphe de S")
      })
  
  output$hgrapheS <- renderPlot ({hgrapheS()
  })
  
  hgrapheL<-reactive({ggplot(hsimul_MC())+
      geom_line(aes(x=time,y=L),col="blue",size=0.4)+
      geom_ribbon(aes(ymin=Linf,ymax=Lsup, x= time), fill="skyblue3",alpha=0.3)+
      theme_classic()+
      #scale_y_continuous(trans="log10")+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      labs(x='time post vaccination(days)',y="ASCS concentration (ASCS/millions)")+
      theme(legend.position = "none") 
      #ggtitle("graphe de L")
      })
  output$hgrapheL <- renderPlot ({
    hgrapheL()})
  
  
  hgrapheAb<-reactive({ggplot(hsimul_MC())+
      geom_line(aes(x=time,y=Ab),col="blue",size=0.4)+
      geom_ribbon(aes(ymin=Abinf,ymax=Absup, x= time), fill="skyblue3",alpha=0.3)+
      theme_classic()+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      #scale_y_continuous(trans="log10")+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      labs(x='time post vaccination(days)',y="Ab concentration (ELISA units)")+
      theme(legend.position = "none") 
      #ggtitle("graphe de Ab")
      })
  
  output$hgrapheAb <- renderPlot({hgrapheAb()
  })
  
  output$hallgraphs = downloadHandler(
    filename = function() {
      'all_images.zip'
    }, 
    content = function(fname) {
      fs <- replicate(4, tempfile(fileext = ".png"))
      ggsave(fs[1], hgrapheM())
      ggsave(fs[2], hgrapheS())
      ggsave(fs[3], hgrapheL())
      ggsave(fs[4], hgrapheAb())
      zip::zipr(zipfile=fname, files=fs)
    },
    contentType = "application/zip")
  
  
  
  #ajout d'un enregistrement
  #un enregistrement=param()
  #table des enregistrements : parametres()
 
  observeEvent(input$hkeep,{
    nom<-ifelse(input$hnommodel=="",paste0("simulation ",input$hkeep+input$keep),input$hnommodel)
    nbinj<-hparam()[1]
    if (nrow(hparametres())==0){
      p <- c (nom,hparam())
      hparam(p)
      t <- rbind (hparametres(),hparam())
      hparametres(t)
      simul<-as.data.frame(hsimul_MC())
      simul$nom<-nom
      simul$modele<- paste("modèle avec régulation",nom,sep=" ")
      j <-rbind(hsimulation_MC(),simul)
      hsimulation_MC(j)
      hbrut <- as.data.frame(hsimul_brute())
      hbrut$nom<-nom
      hbrut$modele<- paste("modèle avec régulation",nom,sep=" ")
      r <- rbind (hsimulation_brute(),hbrut)
      hsimulation_brute(r)
    }
    
    else {
      #test pour voir si le nom existe déjà :
      if(nom %in% hparametres()[,1]){
        showNotification("Vous avez déja utlisé ce nom")
      }
      else{
        #ajoute des NA au nouvel enregistrement si nb injection inférieur pour avoir les même tailles de lignes
        if(length(hparam())<length(hparametres())-1){
          n <- (length(hparametres())-1) - length(hparam())
          a <- rep(NA,n)
          par<-c(hparam(),a)
          hparam(par)
        }
        
        if (length(hparam())==length(hparametres())-1){
          nom<-ifelse(input$hnommodel=="",paste0("simulation ",input$keep+input$hkeep),input$hnommodel)
          p<-c(nom,hparam())
          t<-rbind(hparametres(),p)
          hparametres(t)
        }
        if (length(hparam())>length(hparametres())){
          n <- length(hparam())-(length(hparametres())-1)
          t=hparametres()
          ajout<-data.frame(matrix(nrow=1,ncol=n))
          t <- cbind (t,ajout)
          nom<-ifelse(input$hnommodel=="",paste0("simulation ",input$keep+input$hkeep),input$hnommodel)
          p<-c(nom,hparam())
          t <- rbind (t,p)
          hparametres(t)
        }
        
      
    
    
    
        simul<-as.data.frame(hsimul_MC())
        simul$nom<-nom
        simul$modele<- paste("modèle avec régulation",nom,sep=" ")
    j <-rbind(hsimulation_MC(),simul)
    hsimulation_MC(j)
    hbrut <- as.data.frame(hsimul_brute())
    hbrut$nom<-nom
    hbrut$modele<- paste("modèle avec régulation",nom,sep=" ")
    r <- rbind (hsimulation_brute(),hbrut)
    hsimulation_brute(r)
      }
    }
    #on gère les colnames de tout le monde :
    
    n<-length(hparametres())
    
    k<-((n-31)/30)+1
    
    noms_param<-c("nom","temps","nbinj","ρ","δ(lρ)","β1","σ(β1)","αM1","σ(αM1)","δA1","σ(δA1)","μS1","σ(lμS1)","αS1","σ(αS1)","μL1","σ(lμL1)","αL1","σ(αL1)", "δM1","σ(δM1)","δS1","σ(δS1)","δL1","σ(δL1)","ϴS1","σ(lϴS1)","ϴL1","σ(lϴL1)","δAb1","σ(δAb1)")
    if (k>=2){
      for (i in 2:(k)){
        ajout<-c(paste0("tinj",i),paste0("σ(tinj)",i),paste0("ρ",i),paste0("σ(lρ)",i),paste0("β",i),paste0("σ(β)",i),paste0("αM",i),paste0("σ(αM)",i),paste0("δA",i),paste0("σ(δA)",i),paste0("μS",i),paste0("σ(lμS)",i),paste0("αS",i),paste0("σ(αS)",i),paste0("μL",i),paste0("σ(lμL)",i),paste0("αL",i),paste0("σ(αL)",i),paste0("δM",i),paste0("σ(δM)",i),paste0("δS",i),paste0("σ(δS)",i),paste0("δL",i),paste0("σ(δL)",i),paste0("ϴS",i),paste0("σ(lϴS)",i),paste0("ϴL",i),paste0("σ(lϴL)",i),paste0("δAb",i),paste0("σδAb",i))
        noms_param<-c(noms_param,ajout)
      }}
    a=hparametres()
    b=hsimulation_MC()
    colnames(a)<-noms_param
    colnames(b) <- c("time","M","S","L","Ab","Minf","Sinf","Linf","Abinf","Msup","Ssup","Lsup","Absup","nom","modele")
    hparametres(a)
    
  })
  
  
  
  #tableau avec cases à cocher
  output[["htable_param"]] <-renderReactable({reactable(hparametres(),
                                                      selection = "multiple", 
                                                      onClick = "select"
                                                      
  )
  }
  )
  
  
  #suppression des enregistrements :
  observeEvent(input$hclear,{
    p<-hparametres()[0,]
    hparametres(p)
    t<-hsimulation_MC()[0,]
    hsimulation_MC(t)
    r <- hsimul_brute()[0,]
    hsimul_brute(r)
  })
  
  
  
  #récuperer les rangs sélectionnés : getReactableState()
  hselected <- reactive(getReactableState("htable_param", "selected"))
  #étape d'après: récupérer les noms des modélisation pour subset dans simulation() et plot que ceux là
  
  hnoms_simul<-reactive({lapply(hselected(),function(i){hparametres()[i,1]})})
  
  #Faire les moyennes et IC des données si les variables sont transformées
  hcomparaison_M<-reactive({
    if (input$hlogM==TRUE){
      transfo<-as.data.table(hsimulation_brute()[,c("time","M","nom")])
      transfo$M<-log10(transfo$M)
      moy<-transfo[,lapply(.SD,mean),by=.(time,nom)]
      inf <-transfo[,lapply(.SD,quantile,probs=.025),by=.(time,nom)]
      colnames(inf)<-c("time","nom","Minf")
      sup<-transfo[,lapply(.SD,quantile,probs=.975),by=.(time,nom)]
      colnames(sup)<-c("time","nom","Msup")
      hcomparaison_M <- data.table::merge.data.table(moy,sup,by=c("time","nom"))
      hcomparaison_M<-data.table::merge.data.table(hcomparaison_M,inf,by=c("time","nom"))
    }
    else{

      hcomparaison_M <- hsimulation_MC()
    }
    hcomparaison_M
  })

  hcomparaison_S<-reactive({if (input$hracS==TRUE){
    transfo<-as.data.table(hsimulation_brute()[,c("time","S","nom")])
    transfo$S<-(transfo$S)^(1/4)
    moy<-transfo[,lapply(.SD,mean),by=.(time,nom)]
    inf <-transfo[,lapply(.SD,quantile,probs=.025,na.rm=TRUE),by=.(time,nom)]
    colnames(inf)<-c("time","nom","Sinf")
    sup<-transfo[,lapply(.SD,quantile,probs=.975,na.rm=TRUE),by=.(time,nom)]
    colnames(sup)<-c("time","nom","Ssup")
    hcomparaison_S <- data.table::merge.data.table(moy,sup,by=c("time","nom"))
    hcomparaison_S<-data.table::merge.data.table(hcomparaison_S,inf,by=c("time","nom"))
  }
    else{

      hcomparaison_S <- hsimulation_MC()
    }
    hcomparaison_S})

  hcomparaison_L<-reactive({if (input$hracL==TRUE){
    transfo<-as.data.table(simulation_brute()[,c("time","L","nom")])
    transfo$L<-(transfo$L)^(1/4)
    moy<-transfo[,lapply(.SD,mean),by=.(time,nom)]
    inf <-transfo[,lapply(.SD,quantile,probs=.025,na.rm=TRUE),by=.(time,nom)]
    colnames(inf)<-c("time","nom","Linf")
    sup<-transfo[,lapply(.SD,quantile,probs=.975,na.rm=TRUE),by=.(time,nom)]
    colnames(sup)<-c("time","nom","Lsup")
    hcomparaison_L <- data.table::merge.data.table(moy,sup,by=c("time","nom"))
    hcomparaison_L<-data.table::merge.data.table(hcomparaison_L,inf,by=c("time","nom"))
  }
    else{

      hcomparaison_L <- hsimulation_MC()
    }
    hcomparaison_L})

  hcomparaison_Ab<-reactive({if (input$hlogAb==TRUE){
    transfo<-as.data.table(hsimulation_brute()[,c("time","Ab","nom")])
    transfo$Ab<-log10(transfo$Ab)
    #comparaison_M<-transfo
    moy<-transfo[,lapply(.SD,mean),by=.(time,nom)]
    inf <-transfo[,lapply(.SD,quantile,probs=.025,na.rm=TRUE),by=.(time,nom)]
    colnames(inf)<-c("time","nom","Abinf")
    sup<-transfo[,lapply(.SD,quantile,probs=.975,na.rm=TRUE),by=.(time,nom)]
    colnames(sup)<-c("time","nom","Absup")
    hcomparaison_Ab <- data.table::merge.data.table(moy,sup,by=c("time","nom"))
    hcomparaison_Ab<-data.table::merge.data.table(hcomparaison_Ab,inf,by=c("time","nom"))
  }
    else{

      hcomparaison_Ab <- hsimulation_MC()
    }
    hcomparaison_Ab})
  
  
  output$hgrapheMcompare <- renderPlot({
    partie<-subset(hcomparaison_M(),nom %in% hnoms_simul())
    ggplot(partie)+
      geom_line(aes(x=time,y=M,color=nom),size=0.4)+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      #scale_y_continuous(trans="log10")+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      geom_ribbon(aes(ymin=Minf,ymax=Msup, x= time, fill=nom),alpha=0.3)+
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE)+
      theme_classic()+
      labs(x='time post vaccination(days)',y="Memory B-cells concentration (ASCS/millions)")+
      theme(legend.position="top")
      #theme(legend.position = "none") +
      #ggtitle("graphe de M")
      }) 
  
  
  output$hgrapheScompare <- renderPlot({
    partie<-subset(hcomparaison_S(),nom %in% hnoms_simul())
    ggplot(partie)+
      geom_line(aes(x=time,y=S,color=nom),size=0.4)+
      geom_ribbon(aes(ymin=Sinf,ymax=Ssup, x= time, fill=nom),alpha=0.3)+
      theme_classic()+
      #scale_y_continuous(trans="log10")+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      labs(x='time post vaccination(days)',y="ASCS concentration (ASCS/millions)")+
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE)+
      theme(legend.position="top")
      #theme(legend.position = "none") +
      #ggtitle("graphe de S")
      }) 
  
  output$hgrapheLcompare <- renderPlot({
    partie<-subset(hcomparaison_L(),nom %in% hnoms_simul())
    ggplot(partie)+
      geom_line(aes(x=time,y=L,color=nom),size=0.4)+
      geom_ribbon(aes(ymin=Linf,ymax=Lsup, x= time, fill=nom),alpha=0.3)+
      theme_classic()+
      #scale_y_continuous(trans="log10")+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      labs(x='time post vaccination(days)',y="ASCS concentration (ASCS/millions)")+
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE)+
      theme(legend.position="top")
      #theme(legend.position = "none") +
      #ggtitle("graphe de L")
      })
  
  output$hgrapheAbcompare <- renderPlot({
    partie<-subset(hcomparaison_Ab(),nom %in% hnoms_simul())
    ggplot(partie)+
      geom_line(aes(x=time,y=Ab,color=nom),size=0.4)+
      geom_ribbon(aes(ymin=Abinf,ymax=Absup, x= time, fill=nom),alpha=0.3)+
      theme_classic()+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) + 
      #scale_y_continuous(trans="log10")+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      labs(x='time post vaccination(days)',y="Ab concentration (ELISA units)")+
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE)+
      theme(legend.position="top")
      #theme(legend.position = "none") +
      #ggtitle("graphe de Ab")
      })
  
  
  output$hdl <- downloadHandler(
    
    filename = function() {
      "homeo.xlsx"
    },
    content = function(filename){
      
      df_list <- list(hparametres(), hsimulation_MC(),hsimulation_brute())
      write.xlsx(x = df_list , file = filename, rowNames = FALSE)
    }
  ) 
  
  
  
  
  
  ################# mise en commun sur le dernier onglet
  
  output$table_param2 <-renderReactable({reactable(parametres(),
                                                 selection = "multiple",
                                                 onClick = "select")}
  )
  selected2 <- reactive(getReactableState("table_param2", "selected"))
  noms_simul2<-reactive({lapply(selected2(),function(i){parametres()[i,1]})})

  
  
  output$htable_param2 <-renderReactable({reactable(hparametres(),
                                                  selection = "multiple",
                                                  onClick = "select")}
  )
  hselected2 <- reactive(getReactableState("htable_param2", "selected"))
  hnoms_simul2<-reactive({lapply(hselected2(),function(i){hparametres()[i,1]})})

  comparaison_brute<-reactiveVal(0)
  
  comparaison_MC<-reactive({
    partie1<-data.table()
    partie2<-data.table()
    partie1_brute<-data.table()
    partie2_brute<-data.table()
    if (length(noms_simul2())!=0){
      partie1<-subset(simulation_MC(),nom %in% noms_simul2())
      partie1_brute <-subset(simulation_brute(),nom %in% noms_simul2())}
    if (length(hnoms_simul2())!=0){
      partie2<-subset(hsimulation_MC(),nom %in% hnoms_simul2())
      partie2_brute<-subset(hsimulation_brute(),nom %in% hnoms_simul2())
    }
    # partie1<-subset(simulation_MC(),nom %in% noms_simul2())
    # partie1_brute <-subset(simulation_brute(),nom %in% noms_simul2())
    # partie2<-subset(hsimulation_MC(),nom %in% hnoms_simul2())
    # partie2_brute<-subset(hsimulation_brute(),nom %in% hnoms_simul2())
    # partie2$nom<-as.character(partie2$nom)
    # #partie2_brute$nom<-as.character(partie2_brute$nom)
    # partie1$nom<-as.character(partie1$nom)
    #partie1_brute$nom<-as.character(partie1_brute$nom)
    if (nrow(partie1)!=0|nrow(partie2)!=0){
      ensemble_MC<-dplyr::bind_rows(partie1,partie2)
      ensemble_brute<-dplyr::bind_rows(partie1_brute,partie2_brute)
      comparaison_brute(ensemble_brute)
    }
    else {
      ensemble_MC <- data.frame(matrix(nrow=0,ncol=15))
      ensemble_brute <- data.frame(matrix(nrow=0,ncol=7))
      colnames(ensemble_brute)<-c("time","M","S","L","Ab","nom","modele")
      colnames(ensemble_MC)<-c("time","M","S","L","Ab","Minf","Sinf","Linf","Abinf","Msup","Ssup","Lsup","Absup","nom","modele")
      comparaison_brute(ensemble_brute)
    }
    
    ensemble_MC
  })
  
  #Faire les moyennes et IC des données si les variables sont transformées
  comparaison_M_total<-reactive({
    if (input$logMtotal==TRUE){
      transfo<-as.data.table(comparaison_brute()[,c("time","M","nom","modele")])
      transfo$M<-log10(transfo$M)
      moy<-transfo[,lapply(.SD,mean),by=.(time,nom,modele)]
      inf <-transfo[,lapply(.SD,quantile,probs=.025),by=.(time,nom,modele)]
      colnames(inf)<-c("time","nom","modele","Minf")
      sup<-transfo[,lapply(.SD,quantile,probs=.975),by=.(time,nom,modele)]
      colnames(sup)<-c("time","nom","modele","Msup")
      comparaison_M_total <- data.table::merge.data.table(moy,sup,by=c("time","nom","modele"))
      comparaison_M_total<-data.table::merge.data.table(comparaison_M_total,inf,by=c("time","nom","modele"))
    }
    else{
      
      comparaison_M_total <- comparaison_MC()
    }
    comparaison_M_total
  })
  
  comparaison_S_total<-reactive({
    if (input$racStotal==TRUE){
    transfo<-as.data.table(comparaison_brute()[,c("time","S","nom","modele")])
    transfo$S<-(transfo$S)^(1/4)
    moy<-transfo[,lapply(.SD,mean),by=.(time,nom,modele)]
    inf <-transfo[,lapply(.SD,quantile,probs=.025,na.rm=TRUE),by=.(time,nom,modele)]
    colnames(inf)<-c("time","nom","modele","Sinf")
    sup<-transfo[,lapply(.SD,quantile,probs=.975,na.rm=TRUE),by=.(time,nom,modele)]
    colnames(sup)<-c("time","nom","modele","Ssup")
    comparaison_S_total <- data.table::merge.data.table(moy,sup,by=c("time","nom","modele"))
    comparaison_S_total<-data.table::merge.data.table(comparaison_S_total,inf,by=c("time","nom","modele"))
  }
    else{
      
      comparaison_S_total <- comparaison_MC()
    }
    comparaison_S_total})
  
  comparaison_L_total<-reactive({if (input$racLtotal==TRUE){
    transfo<-as.data.table(comparaison_brute()[,c("time","L","nom","modele")])
    transfo$L<-(transfo$L)^(1/4)
    moy<-transfo[,lapply(.SD,mean),by=.(time,nom,modele)]
    inf <-transfo[,lapply(.SD,quantile,probs=.025,na.rm=TRUE),by=.(time,nom,modele)]
    colnames(inf)<-c("time","nom","modele","Linf")
    sup<-transfo[,lapply(.SD,quantile,probs=.975,na.rm=TRUE),by=.(time,nom,modele)]
    colnames(sup)<-c("time","nom","modele","Lsup")
    comparaison_L_total <- data.table::merge.data.table(moy,sup,by=c("time","nom","modele"))
    comparaison_L_total<-data.table::merge.data.table(comparaison_L_total,inf,by=c("time","nom","modele"))
  }
    else{
      
      comparaison_L_total <- comparaison_MC()
    }
    comparaison_L_total})
  
  comparaison_Ab_total<-reactive({if (input$logAbtotal==TRUE){
    transfo<-as.data.table(comparaison_brute()[,c("time","Ab","nom","modele")])
    transfo$Ab<-log10(transfo$Ab)
    #comparaison_M<-transfo
    moy<-transfo[,lapply(.SD,mean),by=.(time,nom,modele)]
    inf <-transfo[,lapply(.SD,quantile,probs=.025,na.rm=TRUE),by=.(time,nom,modele)]
    colnames(inf)<-c("time","nom","modele","Abinf")
    sup<-transfo[,lapply(.SD,quantile,probs=.975,na.rm=TRUE),by=.(time,nom,modele)]
    colnames(sup)<-c("time","nom","modele","Absup")
    comparaison_Ab_total <- data.table::merge.data.table(moy,sup,by=c("time","nom","modele"))
    comparaison_Ab_total<-data.table::merge.data.table(comparaison_Ab_total,inf,by=c("time","nom","modele"))
  }
    else{
      
      comparaison_Ab_total <- comparaison_MC()
    }
    comparaison_Ab_total})
  
  
  
  grapheMtotal<-reactive({
    ggplot(comparaison_M_total())+
      geom_line(aes(x=time,y=M,color=modele),size=0.4)+
      
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) +
      #scale_y_continuous(trans="log10")+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      geom_ribbon(aes(ymin=Minf,ymax=Msup, x= time, fill=modele),alpha=0.3)+
      theme_classic()+
      labs(x='time post vaccination(days)',y="Memory B-cells concentration (ASCS/millions)")+
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE)+
      theme(legend.position="top")
      #theme(legend.position = "none") +
      #ggtitle("graphe de M")
      })
  
  output$grapheMtotal <- renderPlot({grapheMtotal()})
  
  grapheStotal<-reactive({
    ggplot(comparaison_S_total())+
      geom_line(aes(x=time,y=S,color=modele),size=0.4)+
      geom_ribbon(aes(ymin=Sinf,ymax=Ssup, x= time, fill=modele),alpha=0.3)+
      theme_classic()+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) +
      #scale_y_continuous(trans="log10")+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE)+
      labs(x='time post vaccination(days)',y="ASCS concentration (ASCS/millions)")+
      theme(legend.position="top")
      #theme(legend.position = "none") +
      #ggtitle("graphe de S")
      })
  
  output$grapheStotal <- renderPlot({grapheStotal()})
  
  grapheLtotal<-reactive({
    ggplot(comparaison_L_total())+
      geom_line(aes(x=time,y=L,color=modele),size=0.4)+
      geom_ribbon(aes(ymin=Linf,ymax=Lsup, x= time, fill=modele),alpha=0.3)+
      theme_classic()+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) +
      #scale_y_continuous(trans="log10")+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE)+
      labs(x='time post vaccination(days)',y="ASCS concentration (ASCS/millions)")+
      theme(legend.position="top")
      #theme(legend.position = "none") +
      #ggtitle("graphe de L")
      })
  
  output$grapheLtotal <- renderPlot({grapheLtotal()})
  
  grapheAbtotal<-reactive({
    ggplot(comparaison_Ab_total())+
      geom_line(aes(x=time,y=Ab,color=modele),size=0.4)+
      geom_ribbon(aes(ymin=Abinf,ymax=Absup, x= time, fill=modele),alpha=0.3)+
      theme_classic()+
      scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) +
      #scale_y_continuous(trans="log10")+
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+
      scale_color_viridis(discrete = TRUE, option = "D")+
      scale_fill_viridis(discrete = TRUE)+
      labs(x='time post vaccination(days)',y="Ab concentration (ELISA units)")+
      theme(legend.position="top")
      #theme(legend.position = "none") +
      #ggtitle("graphe de Ab")
      })
  
  output$grapheAbtotal <- renderPlot({grapheAbtotal()})
  
  #enregistrementt de ces petits plots
  output$compare = downloadHandler(
    filename = function() {
      'graphes_compares.zip'
    }, 
    content = function(fname) {
      fs <- replicate(4, tempfile(fileext = ".png"))
      ggsave(fs[1], grapheMtotal())
      ggsave(fs[2], grapheStotal())
      ggsave(fs[3], grapheLtotal())
      ggsave(fs[4], grapheAbtotal())
      zip::zipr(zipfile=fname, files=fs)
    },
    contentType = "application/zip")
})