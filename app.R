##### PERIFIX INTERACTIVE WEB APP #####
##### BY: SAMANTHA GLEICH #####

# Load packages
library(shiny)
library(tidyverse)
library(ggplot2)
library(plotly)
library(reshape2)
library(decontam)
library(randomcoloR)
library(scales)

# Load data
finalASVs <- read.csv("data.csv",header=TRUE)

# Set up ShinyApp fluid page
ui <- fluidPage(
  titlePanel("PERIFIX 18S AMPLICON DATA"),
  sidebarLayout(
    sidebarPanel(
      h4("PERIFIX 18S AMPLICON DATA"),
      selectInput("Group","Select Taxonomic Group:",c("Chlorophyte","Ciliate","Diatom","Dinoflagellate","Haptophyte","Metazoa","All Eukaryotes"), selected="All Eukaryotes"),
      tags$hr(style="border-color: black;"),),
    mainPanel(
      plotlyOutput("plot"))))

# Set up server
server <- function (input, output, session){
  output$plot <- renderPlotly({
    # Chlorophyte
    if (input$Group=="Chlorophyte"){
      # Wrangle data - calculate averages
      new <- subset(finalASVs,Group=="Chlorophyte")
      new$Group <- NULL
      new$Feature.ID <- NULL
      cols <- colsplit(new$Taxon,";",c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
      new$Taxon <- cols$Class
      new$Taxon <- ifelse(new$Taxon=="","Unknown Chlorophyte",new$Taxon)
      newMelt <- melt(new,id.vars="Taxon")
      newMelt <- newMelt %>% group_by(variable,Taxon) %>% summarize(value=sum(value))
      samps <- colsplit(newMelt$variable,"\\.",c("Treatment","Rep","Month","Day"))
      newMelt$Treatment <- samps$Treatment
      newMelt$Date <- paste(samps$Month,samps$Day,"2021",sep="-")
      newMelt$variable <- NULL
      newMeltAvg <- newMelt %>% group_by(Treatment,Taxon,Date) %>% summarize(m=mean(value))
      newMeltAvg$Date <- as.Date(newMeltAvg$Date,format='%m-%d-%Y')
      
      colrs <- distinctColorPalette(length(unique(newMeltAvg$Taxon)))
      
      # Rename treatments
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="C","Control",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="F","+Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="N","+N",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NF","+N +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NP","+N +P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NPF","+N +P +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="P","+P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="PF","+P +Fe",newMeltAvg$Treatment)
      
      # Make the tote the T0 sample for all treatments
      toteC <- subset(newMeltAvg,Treatment=="Tote")
      toteC$Treatment <- "Control"
      
      toteF <- subset(newMeltAvg,Treatment=="Tote")
      toteF$Treatment <- "+Fe"
      
      toteN <- subset(newMeltAvg,Treatment=="Tote")
      toteN$Treatment <- "+N"
      
      toteNF <- subset(newMeltAvg,Treatment=="Tote")
      toteNF$Treatment <- "+N +Fe"
      
      toteNP <- subset(newMeltAvg,Treatment=="Tote")
      toteNP$Treatment <- "+N +P"
      
      toteNPF <- subset(newMeltAvg,Treatment=="Tote")
      toteNPF$Treatment <- "+N +P +Fe"
      
      toteP <- subset(newMeltAvg,Treatment=="Tote")
      toteP$Treatment <- "+P"
      
      totePF <- subset(newMeltAvg,Treatment=="Tote")
      totePF$Treatment <- "+P +Fe"
      
      tote <- rbind(toteC,toteF,toteN,toteNF,toteNP,toteNPF,toteP,totePF)
      newMeltAvg <- subset(newMeltAvg,Treatment!="Tote")
      newMeltAvg <- rbind(newMeltAvg,tote)
    
      
      # Ggplot 
      p<- newMeltAvg %>% ggplot(aes(x=Date,y=m,fill=Taxon))+ geom_bar(stat="identity",position="fill")+theme_classic()+facet_wrap(~Treatment)+scale_x_date(labels = date_format("%m-%d-%Y"))+theme(axis.text.x = element_text(angle = 90, hjust=1))+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+ylab("Relative Abundance")+xlab("")}
    
    # Ciliate
    if (input$Group=="Ciliate"){
      # Wrangle data - calculate averages
      new <- subset(finalASVs,Group=="Ciliate")
      new$Group <- NULL
      new$Feature.ID <- NULL
      cols <- colsplit(new$Taxon,";",c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
      new$Taxon <- cols$Class
      new$Taxon <- ifelse(new$Taxon=="","Unknown Ciliate",new$Taxon)
      newMelt <- melt(new,id.vars="Taxon")
      newMelt <- newMelt %>% group_by(variable,Taxon) %>% summarize(value=sum(value))
      samps <- colsplit(newMelt$variable,"\\.",c("Treatment","Rep","Month","Day"))
      newMelt$Treatment <- samps$Treatment
      newMelt$Date <- paste(samps$Month,samps$Day,"2021",sep="-")
      newMelt$variable <- NULL
      newMeltAvg <- newMelt %>% group_by(Treatment,Taxon,Date) %>% summarize(m=mean(value))
      newMeltAvg$Date <- as.Date(newMeltAvg$Date,format='%m-%d-%Y')
      
      colrs <- distinctColorPalette(length(unique(newMeltAvg$Taxon)))
      
      # Rename treatments
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="C","Control",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="F","+Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="N","+N",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NF","+N +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NP","+N +P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NPF","+N +P +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="P","+P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="PF","+P +Fe",newMeltAvg$Treatment)
      
      # Make the tote the T0 sample for all treatments
      toteC <- subset(newMeltAvg,Treatment=="Tote")
      toteC$Treatment <- "Control"
      
      toteF <- subset(newMeltAvg,Treatment=="Tote")
      toteF$Treatment <- "+Fe"
      
      toteN <- subset(newMeltAvg,Treatment=="Tote")
      toteN$Treatment <- "+N"
      
      toteNF <- subset(newMeltAvg,Treatment=="Tote")
      toteNF$Treatment <- "+N +Fe"
      
      toteNP <- subset(newMeltAvg,Treatment=="Tote")
      toteNP$Treatment <- "+N +P"
      
      toteNPF <- subset(newMeltAvg,Treatment=="Tote")
      toteNPF$Treatment <- "+N +P +Fe"
      
      toteP <- subset(newMeltAvg,Treatment=="Tote")
      toteP$Treatment <- "+P"
      
      totePF <- subset(newMeltAvg,Treatment=="Tote")
      totePF$Treatment <- "+P +Fe"
      
      tote <- rbind(toteC,toteF,toteN,toteNF,toteNP,toteNPF,toteP,totePF)
      newMeltAvg <- subset(newMeltAvg,Treatment!="Tote")
      newMeltAvg <- rbind(newMeltAvg,tote)
      
      # Ggplot 
      p<- newMeltAvg %>% ggplot(aes(x=Date,y=m,fill=Taxon))+ geom_bar(stat="identity",position="fill")+theme_classic()+facet_wrap(~Treatment)+scale_x_date(labels = date_format("%m-%d-%Y"))+theme(axis.text.x = element_text(angle = 90, hjust=1))+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+ylab("Relative Abundance")+xlab("")}
    
    # Diatom
    if (input$Group=="Diatom"){
      # Wrangle data - calculate averages
      new <- subset(finalASVs,Group=="Diatom")
      new$Group <- NULL
      new$Feature.ID <- NULL
      cols <- colsplit(new$Taxon,";",c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
      new$Taxon <- cols$Genus
      new$Taxon <- ifelse(new$Taxon=="","Unknown Diatom",new$Taxon)
      newMelt <- melt(new,id.vars="Taxon")
      newMelt <- newMelt %>% group_by(variable,Taxon) %>% summarize(value=sum(value))
      samps <- colsplit(newMelt$variable,"\\.",c("Treatment","Rep","Month","Day"))
      newMelt$Treatment <- samps$Treatment
      newMelt$Date <- paste(samps$Month,samps$Day,"2021",sep="-")
      newMelt$variable <- NULL
      newMeltAvg <- newMelt %>% group_by(Treatment,Taxon,Date) %>% summarize(m=mean(value))
      newMeltAvg$Date <- as.Date(newMeltAvg$Date,format='%m-%d-%Y')
      
      colrs <- distinctColorPalette(length(unique(newMeltAvg$Taxon)))
      
      # Rename treatments
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="C","Control",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="F","+Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="N","+N",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NF","+N +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NP","+N +P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NPF","+N +P +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="P","+P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="PF","+P +Fe",newMeltAvg$Treatment)
      
      # Make the tote the T0 sample for all treatments
      toteC <- subset(newMeltAvg,Treatment=="Tote")
      toteC$Treatment <- "Control"
      
      toteF <- subset(newMeltAvg,Treatment=="Tote")
      toteF$Treatment <- "+Fe"
      
      toteN <- subset(newMeltAvg,Treatment=="Tote")
      toteN$Treatment <- "+N"
      
      toteNF <- subset(newMeltAvg,Treatment=="Tote")
      toteNF$Treatment <- "+N +Fe"
      
      toteNP <- subset(newMeltAvg,Treatment=="Tote")
      toteNP$Treatment <- "+N +P"
      
      toteNPF <- subset(newMeltAvg,Treatment=="Tote")
      toteNPF$Treatment <- "+N +P +Fe"
      
      toteP <- subset(newMeltAvg,Treatment=="Tote")
      toteP$Treatment <- "+P"
      
      totePF <- subset(newMeltAvg,Treatment=="Tote")
      totePF$Treatment <- "+P +Fe"
      
      tote <- rbind(toteC,toteF,toteN,toteNF,toteNP,toteNPF,toteP,totePF)
      newMeltAvg <- subset(newMeltAvg,Treatment!="Tote")
      newMeltAvg <- rbind(newMeltAvg,tote)
      
      # Ggplot 
      p<- newMeltAvg %>% ggplot(aes(x=Date,y=m,fill=Taxon))+ geom_bar(stat="identity",position="fill")+theme_classic()+facet_wrap(~Treatment)+scale_x_date(labels = date_format("%m-%d-%Y"))+theme(axis.text.x = element_text(angle = 90, hjust=1))+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+ylab("Relative Abundance")+xlab("")}
    
    # Dinoflagellate
    if (input$Group=="Dinoflagellate"){
      # Wrangle data - calculate averages
      new <- subset(finalASVs,Group=="Dinoflagellate")
      new$Group <- NULL
      new$Feature.ID <- NULL
      cols <- colsplit(new$Taxon,";",c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
      new$Taxon <- cols$Family
      new$Taxon <- ifelse(new$Taxon=="","Unknown Dinoflagellate",new$Taxon)
      newMelt <- melt(new,id.vars="Taxon")
      newMelt <- newMelt %>% group_by(variable,Taxon) %>% summarize(value=sum(value))
      samps <- colsplit(newMelt$variable,"\\.",c("Treatment","Rep","Month","Day"))
      newMelt$Treatment <- samps$Treatment
      newMelt$Date <- paste(samps$Month,samps$Day,"2021",sep="-")
      newMelt$variable <- NULL
      newMeltAvg <- newMelt %>% group_by(Treatment,Taxon,Date) %>% summarize(m=mean(value))
      newMeltAvg$Date <- as.Date(newMeltAvg$Date,format='%m-%d-%Y')
      
      colrs <- distinctColorPalette(length(unique(newMeltAvg$Taxon)))
      
      # Rename treatments
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="C","Control",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="F","+Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="N","+N",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NF","+N +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NP","+N +P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NPF","+N +P +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="P","+P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="PF","+P +Fe",newMeltAvg$Treatment)
      
      # Make the tote the T0 sample for all treatments
      toteC <- subset(newMeltAvg,Treatment=="Tote")
      toteC$Treatment <- "Control"
      
      toteF <- subset(newMeltAvg,Treatment=="Tote")
      toteF$Treatment <- "+Fe"
      
      toteN <- subset(newMeltAvg,Treatment=="Tote")
      toteN$Treatment <- "+N"
      
      toteNF <- subset(newMeltAvg,Treatment=="Tote")
      toteNF$Treatment <- "+N +Fe"
      
      toteNP <- subset(newMeltAvg,Treatment=="Tote")
      toteNP$Treatment <- "+N +P"
      
      toteNPF <- subset(newMeltAvg,Treatment=="Tote")
      toteNPF$Treatment <- "+N +P +Fe"
      
      toteP <- subset(newMeltAvg,Treatment=="Tote")
      toteP$Treatment <- "+P"
      
      totePF <- subset(newMeltAvg,Treatment=="Tote")
      totePF$Treatment <- "+P +Fe"
      
      tote <- rbind(toteC,toteF,toteN,toteNF,toteNP,toteNPF,toteP,totePF)
      newMeltAvg <- subset(newMeltAvg,Treatment!="Tote")
      newMeltAvg <- rbind(newMeltAvg,tote)
      
      # Ggplot 
      p<- newMeltAvg %>% ggplot(aes(x=Date,y=m,fill=Taxon))+ geom_bar(stat="identity",position="fill")+theme_classic()+facet_wrap(~Treatment)+scale_x_date(labels = date_format("%m-%d-%Y"))+theme(axis.text.x = element_text(angle = 90, hjust=1))+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+ylab("Relative Abundance")+xlab("")}
    
    # Metazoa
    if (input$Group=="Metazoa"){
      # Wrangle data - calculate averages
      new <- subset(finalASVs,Group=="Metazoa")
      new$Group <- NULL
      new$Feature.ID <- NULL
      cols <- colsplit(new$Taxon,";",c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
      new$Taxon <- cols$Family
      new$Taxon <- ifelse(new$Taxon=="","Unknown Metazoan",new$Taxon)
      newMelt <- melt(new,id.vars="Taxon")
      newMelt <- newMelt %>% group_by(variable,Taxon) %>% summarize(value=sum(value))
      samps <- colsplit(newMelt$variable,"\\.",c("Treatment","Rep","Month","Day"))
      newMelt$Treatment <- samps$Treatment
      newMelt$Date <- paste(samps$Month,samps$Day,"2021",sep="-")
      newMelt$variable <- NULL
      newMeltAvg <- newMelt %>% group_by(Treatment,Taxon,Date) %>% summarize(m=mean(value))
      newMeltAvg$Date <- as.Date(newMeltAvg$Date,format='%m-%d-%Y')
      
      colrs <- distinctColorPalette(length(unique(newMeltAvg$Taxon)))
      
      # Rename treatments
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="C","Control",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="F","+Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="N","+N",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NF","+N +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NP","+N +P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NPF","+N +P +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="P","+P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="PF","+P +Fe",newMeltAvg$Treatment)
      
      # Make the tote the T0 sample for all treatments
      toteC <- subset(newMeltAvg,Treatment=="Tote")
      toteC$Treatment <- "Control"
      
      toteF <- subset(newMeltAvg,Treatment=="Tote")
      toteF$Treatment <- "+Fe"
      
      toteN <- subset(newMeltAvg,Treatment=="Tote")
      toteN$Treatment <- "+N"
      
      toteNF <- subset(newMeltAvg,Treatment=="Tote")
      toteNF$Treatment <- "+N +Fe"
      
      toteNP <- subset(newMeltAvg,Treatment=="Tote")
      toteNP$Treatment <- "+N +P"
      
      toteNPF <- subset(newMeltAvg,Treatment=="Tote")
      toteNPF$Treatment <- "+N +P +Fe"
      
      toteP <- subset(newMeltAvg,Treatment=="Tote")
      toteP$Treatment <- "+P"
      
      totePF <- subset(newMeltAvg,Treatment=="Tote")
      totePF$Treatment <- "+P +Fe"
      
      tote <- rbind(toteC,toteF,toteN,toteNF,toteNP,toteNPF,toteP,totePF)
      newMeltAvg <- subset(newMeltAvg,Treatment!="Tote")
      newMeltAvg <- rbind(newMeltAvg,tote)
      
      # Ggplot 
      p<- newMeltAvg %>% ggplot(aes(x=Date,y=m,fill=Taxon))+ geom_bar(stat="identity",position="fill")+theme_classic()+facet_wrap(~Treatment)+scale_x_date(labels = date_format("%m-%d-%Y"))+theme(axis.text.x = element_text(angle = 90, hjust=1))+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+ylab("Relative Abundance")+xlab("")}
    
    
    # Haptophyte
    if (input$Group=="Haptophyte"){
      # Wrangle data - calculate averages
      new <- subset(finalASVs,Group=="Haptophyte")
      new$Group <- NULL
      new$Feature.ID <- NULL
      cols <- colsplit(new$Taxon,";",c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
      new$Taxon <- cols$Order
      new$Taxon <- ifelse(new$Taxon=="","Unknown Haptophyte",new$Taxon)
      newMelt <- melt(new,id.vars="Taxon")
      newMelt <- newMelt %>% group_by(variable,Taxon) %>% summarize(value=sum(value))
      samps <- colsplit(newMelt$variable,"\\.",c("Treatment","Rep","Month","Day"))
      newMelt$Treatment <- samps$Treatment
      newMelt$Date <- paste(samps$Month,samps$Day,"2021",sep="-")
      newMelt$variable <- NULL
      newMeltAvg <- newMelt %>% group_by(Treatment,Taxon,Date) %>% summarize(m=mean(value))
      newMeltAvg$Date <- as.Date(newMeltAvg$Date,format='%m-%d-%Y')
      
      colrs <- distinctColorPalette(length(unique(newMeltAvg$Taxon)))
      
      # Rename treatments
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="C","Control",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="F","+Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="N","+N",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NF","+N +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NP","+N +P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NPF","+N +P +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="P","+P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="PF","+P +Fe",newMeltAvg$Treatment)
      
      # Make the tote the T0 sample for all treatments
      toteC <- subset(newMeltAvg,Treatment=="Tote")
      toteC$Treatment <- "Control"
      
      toteF <- subset(newMeltAvg,Treatment=="Tote")
      toteF$Treatment <- "+Fe"
      
      toteN <- subset(newMeltAvg,Treatment=="Tote")
      toteN$Treatment <- "+N"
      
      toteNF <- subset(newMeltAvg,Treatment=="Tote")
      toteNF$Treatment <- "+N +Fe"
      
      toteNP <- subset(newMeltAvg,Treatment=="Tote")
      toteNP$Treatment <- "+N +P"
      
      toteNPF <- subset(newMeltAvg,Treatment=="Tote")
      toteNPF$Treatment <- "+N +P +Fe"
      
      toteP <- subset(newMeltAvg,Treatment=="Tote")
      toteP$Treatment <- "+P"
      
      totePF <- subset(newMeltAvg,Treatment=="Tote")
      totePF$Treatment <- "+P +Fe"
      
      tote <- rbind(toteC,toteF,toteN,toteNF,toteNP,toteNPF,toteP,totePF)
      newMeltAvg <- subset(newMeltAvg,Treatment!="Tote")
      newMeltAvg <- rbind(newMeltAvg,tote)
      
      # Ggplot 
      p<- newMeltAvg %>% ggplot(aes(x=Date,y=m,fill=Taxon))+ geom_bar(stat="identity",position="fill")+theme_classic()+facet_wrap(~Treatment)+scale_x_date(labels = date_format("%m-%d-%Y"))+theme(axis.text.x = element_text(angle = 90, hjust=1))+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+ylab("Relative Abundance")+xlab("")}
    
    
    # All
    if (input$Group=="All Eukaryotes"){
      # Wrangle data - calculate averages
      new <- finalASVs
      new$Taxon <- NULL
      new$Feature.ID <- NULL
      new$Taxon <- new$Group
      new$Group <- NULL
      newMelt <- melt(new,id.vars="Taxon")
      newMelt <- newMelt %>% group_by(variable,Taxon) %>% summarize(value=sum(value))
      samps <- colsplit(newMelt$variable,"\\.",c("Treatment","Rep","Month","Day"))
      newMelt$Treatment <- samps$Treatment
      newMelt$Date <- paste(samps$Month,samps$Day,"2021",sep="-")
      newMelt$variable <- NULL
      newMeltAvg <- newMelt %>% group_by(Treatment,Taxon,Date) %>% summarize(m=mean(value))
      newMeltAvg$Date <- as.Date(newMeltAvg$Date,format='%m-%d-%Y')
      
      colrs <- distinctColorPalette(length(unique(newMeltAvg$Taxon)))
      
      # Rename treatments
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="C","Control",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="F","+Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="N","+N",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NF","+N +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NP","+N +P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="NPF","+N +P +Fe",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="P","+P",newMeltAvg$Treatment)
      newMeltAvg$Treatment <- ifelse(newMeltAvg$Treatment=="PF","+P +Fe",newMeltAvg$Treatment)
      
      # Make the tote the T0 sample for all treatments
      toteC <- subset(newMeltAvg,Treatment=="Tote")
      toteC$Treatment <- "Control"
      
      toteF <- subset(newMeltAvg,Treatment=="Tote")
      toteF$Treatment <- "+Fe"
      
      toteN <- subset(newMeltAvg,Treatment=="Tote")
      toteN$Treatment <- "+N"
      
      toteNF <- subset(newMeltAvg,Treatment=="Tote")
      toteNF$Treatment <- "+N +Fe"
      
      toteNP <- subset(newMeltAvg,Treatment=="Tote")
      toteNP$Treatment <- "+N +P"
      
      toteNPF <- subset(newMeltAvg,Treatment=="Tote")
      toteNPF$Treatment <- "+N +P +Fe"
      
      toteP <- subset(newMeltAvg,Treatment=="Tote")
      toteP$Treatment <- "+P"
      
      totePF <- subset(newMeltAvg,Treatment=="Tote")
      totePF$Treatment <- "+P +Fe"
      
      tote <- rbind(toteC,toteF,toteN,toteNF,toteNP,toteNPF,toteP,totePF)
      newMeltAvg <- subset(newMeltAvg,Treatment!="Tote")
      newMeltAvg <- rbind(newMeltAvg,tote)
      
      # Ggplot 
      p<- newMeltAvg %>% ggplot(aes(x=Date,y=m,fill=Taxon))+ geom_bar(stat="identity",position="fill")+theme_classic()+facet_wrap(~Treatment)+scale_x_date(labels = date_format("%m-%d-%Y"))+theme(axis.text.x = element_text(angle = 90, hjust=1))+scale_fill_manual(name="Taxonomic Group",values=c(colrs))+ylab("Relative Abundance")+xlab("")}
    
    
    return(p)})}


# Run app
shinyApp(ui = ui, server = server)
# runApp()
