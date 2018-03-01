#
# This is a Shiny web application to estimated total mortality rate from length frequency observations.
# The app is based on the Ehrhardt Ault mortality estimation model. 
# Ehrhardt, N. M., and Ault, J. S. 1992. Analysis of two length-based mortality models applied to 
# bounded catch length frequencies. Transactions of the American Fisheries Society, 121: 115-122.
#

library(shiny)
library(plotrix)
# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel(
     fluidRow(column(12,h1("Length based mortality rate estimator"),
              h4("This tool is designed to estimate total mortality rate
                 from the average length of observed individuals within a 
                 fixed size/age window. It utilizes the mortality estimation 
                 model of Ehrhardt & Ault 1992."),
              h6("Ehrhardt, N. M. and J. S. Ault. 1992. Analysis of two 
                 length-based mortality models applied to bounded catch 
                 length frequencies. Transactions of the American Fisheries 
                 Society. 121(1): 115-122.")))
     ),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
        
        fluidRow(column(12,
                         radioButtons("lengSource","Source of length frequence",c("Data","Simulation"),selected="Simulation"))),
        
        fluidRow( 
          column(12,
                 textInput("SimSampSize","Simulated annual sample sizes",value="100")),
          column(12,
                 textInput("SimYears","Simulated years",value="2014")),
          column(12,
                 textInput("SimZs","Simulated mortality rates",value="0.3"))
        ),
        
        fluidRow(
          column(12,downloadButton("saveSim","Download simulation data"))
        ),
        
        fileInput("lengFreq", "Choose length frequency CSV file",
                  multiple = TRUE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
        
        
        fluidRow(
          column(12,
                 strong("Life history parameter inputs"))
        ),
        
        fluidRow( 
          column(12,
                 textInput("SpeciesName","Species name",value="Example fish 1"))
          ),
        
        fluidRow( column(4,
                 textInput("t0Val","t0",value="0")),
          column(4,
                 textInput("KVal","K",value="0.2")),
          column(4,
                 textInput("LinfVal","Linf",value="70")),
          column(4,
                 textInput("MaxAgeVal","Max Age",value="18")),
          column(4,
                 textInput("LcVal","Min Size",value="16"))
        ),
       
        fluidRow(
          column(12,downloadButton("saveLifeHist","Download life history values"))
        ),
        
        fileInput("lifeHist", "Choose life history CSV file",
                  multiple = TRUE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        
         sliderInput("bins",
                     "Number of bins for histogram display:",
                     min = 1,
                     max = 100,
                     value = 30),
        
        fluidRow(
          column(12,downloadButton("saveResults","Download results"))
        )
        
        
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("lengthHist"),
         plotOutput("growth"),
         plotOutput("meanLength"),
         plotOutput("mortality"),
         tableOutput('testTable')
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
   
      lengFreqValReact<-reactiveVal(value=matrix(NA,nrow=1,ncol=3))
      
      output$saveSim<-downloadHandler(
        filename = function() {
          paste(input$SpeciesName,"_",input$SimSampSize,"_",input$SimZs,"_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          write.table(round(lengFreqVals,2), file,row.names = FALSE,col.names=c("length","count","year"),sep=",",append=FALSE)
        },
        contentType = "text/csv"
      )
      
      getLifeHist<-function(){
        lifeHistValsOut<-matrix(NA,ncol=6,nrow=1)
        lifeHistValsOut[1,]<-c(input$SpeciesName,input$KVal,input$LinfVal,input$t0Val,input$MaxAgeVal,input$LcVal)
        return(lifeHistValsOut)
      }
      
      output$saveLifeHist<-downloadHandler(
        filename = function() {
          paste(input$SpeciesName,"_LifeHist_", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
          write.table(getLifeHist(), file,row.names = FALSE,col.names=c("Species","K","Linf","t0","MaxAge","MinSize"),sep=",",append=FALSE)
        },
        contentType = "text/csv"
      )
      
      LengFreqHist<-function(){}
      growFunc<-function(){}
      meanLengBox<-function(){}
      zBox<-function(){}
      
      output$saveResults<-downloadHandler(
        filename = function() {
          paste(input$SpeciesName,"_Results_", Sys.Date(), ".pdf", sep="")
        },
        content = function(file) {
          pdf(file)
          print( LengFreqHist() )
          print( growFunc() )
          print( meanLengBox() )  
          print( zBox() )
          dev.off()
        }
      )
  
      CalcLbar<-function(Z,Lc,Lmax,Linf,K,Lbar)
      {
        ratio<-(Linf-Lmax)/(Linf-Lc)
        EstimLbar<-Linf-(Linf-Lc)*(Z/(Z+K))*((1-ratio^((Z+K)/K))/(1-ratio^(Z/K)))
        diff<-abs(Lbar-EstimLbar)
      }
   
      observeEvent(input$lifeHist,{
        if(!is.null(input$lifeHist))
        {
           lifeHistVals <<- read.csv(input$lifeHist$datapath,header=TRUE)
           updateTextInput(session=session, inputId = "SpeciesName", value = lifeHistVals$Species)
           updateTextInput(session=session, inputId = "t0Val", value = lifeHistVals$t0)
           updateTextInput(session=session, inputId = "KVal", value = lifeHistVals$K)
           updateTextInput(session=session, inputId = "LinfVal", value = lifeHistVals$Linf)
           updateTextInput(session=session, inputId = "MaxAgeVal", value = lifeHistVals$MaxAge)
           updateTextInput(session=session, inputId = "LcVal", value = lifeHistVals$MinSize)
        }
      })
      
      
      output$lengthHist <- renderPlot({
     #req(input$lengFreq)
        test<-input$SimSampSize
        test<-input$SimYears
        test<-input$SimZs
        test<-input$KVal
        test<-input$t0Val
        test<-input$LinfVal
        test<-input$LcVal
        test<-input$MaxAgeVal
        test<-input$SpeciesName
        LengFreqHist<<-function(){
      if(input$lengSource=="Simulation")
      {
        temp<-strsplit(input$SimZs,"[, ]+")[[1]]
        mortHist<-as.numeric(strsplit(input$SimZs,"[, ]+")[[1]])
        mortHist<-c(rep(mortHist[1],as.numeric(input$MaxAgeVal)),mortHist)
        mortHist<-rep(mortHist/100,each=100)
        cumMortHist<-cumsum(rev(mortHist))
        lengAtAge<-as.numeric(input$LinfVal)*(1-exp(-as.numeric(input$KVal)*((c(1:(100*(as.numeric(input$MaxAgeVal))))/100)-as.numeric(input$t0Val))))
        abunAtAge<-exp(-cumMortHist) 
        abunAtAge<-abunAtAge[1:length(lengAtAge)]
        abunAtAge[lengAtAge<as.numeric(input$LcVal)]<-0
        lengFreqVals<-matrix(NA,nrow=as.numeric(input$SimSampSize),ncol=3)
        lengFreqVals[,1]<-sample(x=lengAtAge,size=as.numeric(input$SimSampSize),replace=TRUE,prob=abunAtAge[1:length(lengAtAge)])
        lengFreqVals[,2]<-rep(1,length(lengFreqVals[,2]))
        lengFreqVals[,3]<-rep(as.numeric(input$SimYears),length(lengFreqVals[,3]))
        xGlobal    <- lengFreqVals[,1]
        bins <- seq(min(xGlobal), max(xGlobal), length.out = input$bins + 1)
        HistValsGlobal<-hist(xGlobal, breaks = bins, plot=FALSE)
        barplot(HistValsGlobal$counts,ylim=c(0,(max(HistValsGlobal$counts)+1)), col=gray(0),border="red",names.arg=round((bins[1:(length(bins)-1)]+bins[2:(length(bins))])/2),ylab="count",xlab="length",main=paste0("Simulated length frequency for ",input$SpeciesName))
        box()
      }
      if(input$lengSource=="Data")
      {
        lengFreqVals <- read.csv(input$lengFreq$datapath,header=TRUE)
        lengFreqSubset<-lengFreqVals
        Years<-unique(lengFreqVals[,3])
        numYears<-length(Years)
        xGlobal    <- rep(lengFreqVals[,1],lengFreqVals[,2])
        bins <- seq(min(xGlobal), max(xGlobal), length.out = input$bins + 1)
        HistValsGlobal<-hist(xGlobal, breaks = bins, plot=FALSE)
        barplot(HistValsGlobal$counts,ylim=c(0,(max(HistValsGlobal$counts)+1)), col=gray(0),border="red",names.arg=round(x=((bins[1:(length(bins)-1)]+bins[2:(length(bins))])/2),digits=2),ylab="count",xlab="length",main=paste0("Observed length frequency for ",input$SpeciesName))
        if(numYears>1)
        {
          for(i in 1:(numYears-1))
          {
            par(new=TRUE)
            lengFreqSubset<-lengFreqSubset[lengFreqSubset[,3]!=Years[i],]
            x    <- rep(lengFreqSubset[,1],lengFreqSubset[,2])
            HistVals<-hist(x, breaks = bins, plot=FALSE)
            barplot(HistVals$counts,ylim=c(0,(max(HistValsGlobal$counts)+1)),col=gray((0+(i)*(1/(numYears-1)))),border="red",names.arg=round(x=((bins[1:(length(bins)-1)]+bins[2:(length(bins))])/2),digits=2),ylab="count",xlab="length",main=paste0("Observed length frequency for ",input$SpeciesName))
            box()
          }
        }
      }
          lengFreqValReact(lengFreqVals)
      }
        LengFreqHist()
   })
   
   output$growth <- renderPlot({
     #req(input$lengFreq)
     test<-lengFreqVals
     test<-input$SimSampSize
     test<-input$SimYears
     test<-input$SimZs
     test<-input$KVal
     test<-input$t0Val
     test<-input$LinfVal
     test<-input$LcVal
     test<-input$MaxAgeVal
     test<-input$SpeciesName
     
     growFunc<<-function(){
     ages<-seq(as.numeric(input$t0Val),as.numeric(input$MaxAgeVal),0.01)
     ageRange<-c((as.numeric(input$t0Val)-abs(as.numeric(input$t0Val)*0.05)),as.numeric(input$MaxAgeVal)*1.1)
     plot(x=ages,y=as.numeric(input$LinfVal)*(1-exp(-as.numeric(input$KVal)*(ages-as.numeric(input$t0Val)))),type="l",xlab="Age", ylab="Length",xlim=ageRange,main=paste(input$SpeciesName," Growth Function",sep=""))
     }
     growFunc()
   })
   
   output$meanLength <- renderPlot({
     #req(input$lengFreq)
     test<-lengFreqVals
     test<-input$SimSampSize
     test<-input$SimYears
     test<-input$SimZs
     test<-input$KVal
     test<-input$t0Val
     test<-input$LinfVal
     test<-input$LcVal
     test<-input$MaxAgeVal
     test<-input$SpeciesName
     meanLengBox<<-function(){
     lengFreqVals<-lengFreqValReact()
     Years<-unique(lengFreqVals[,3])
     numYears<-length(Years)
     meanBox<-list()
     meanBox$stats<-matrix(NA,nrow=5,ncol=numYears)
     meanBox$n<-rep(0,numYears)
     meanBox$conf<-matrix(NA,nrow=2,ncol=numYears)
     meanBox$out<-0
     meanBox$group<-0
     meanBox$names<-as.character(Years)
     
     for(i in 1:numYears)
     {
       lengFreqSubset<-lengFreqVals[lengFreqVals[,3]==Years[i],]
       lengths<-rep(lengFreqSubset[,1],lengFreqSubset[,2])
       lengths<-lengths[lengths>=input$LcVal]
       meanlength<-mean(lengths)
       SEmeanlength<-sqrt(var(lengths))/sqrt(length(lengths))
       meanBox$stats[1,i]<-meanlength-1.96*SEmeanlength
       meanBox$stats[2,i]<-meanlength-SEmeanlength
       meanBox$stats[3,i]<-meanlength
       meanBox$stats[4,i]<-meanlength+SEmeanlength
       meanBox$stats[5,i]<-meanlength+1.96*SEmeanlength
       meanBox$n[i]<-length(lengths)
     }
    bxp(meanBox,ylim=c(0.95*min(meanBox$stats),1.05*max(meanBox$stats)),ylab="Length",xlab="Year",main=paste0("Annual mean length of ",input$SpeciesName),show.names = TRUE)
     }
     meanLengBox()
  })
   
   output$mortality <- renderPlot({
     #req(input$lengFreq)
     test<-input$SimSampSize
     test<-input$SimYears
     test<-input$SimZs
     test<-input$KVal
     test<-input$t0Val
     test<-input$LinfVal
     test<-input$LcVal
     test<-input$MaxAgeVal
     test<-input$SpeciesName
      zBox<<-function(){
     lengFreqVals<-lengFreqValReact()
     Years<-unique(lengFreqVals[,3])
     numYears<-length(Years)
     ZBox<-list()
     ZBox$stats<-matrix(NA,nrow=5,ncol=numYears)
     ZBox$n<-rep(0,numYears)
     ZBox$conf<-matrix(NA,nrow=2,ncol=numYears)
     ZBox$out<-0
     ZBox$group<-0
     ZBox$names<-as.character(Years)
     ZVect<-matrix(NA,nrow=numYears,ncol=7)
     for(i in 1:numYears)
     {
       lengFreqSubset<-lengFreqVals[lengFreqVals[,3]==Years[i],]
       lengths<-rep(lengFreqSubset[,1],lengFreqSubset[,2])
       lengths<-lengths[lengths>=input$LcVal]
       meanlength<-mean(lengths)
       SEmeanlength<-sqrt(var(lengths))/sqrt(length(lengths))
       ZVect[i,1]<-meanlength
       ZVect[i,2]<-SEmeanlength
       BHZ<-as.numeric(input$KVal)*(as.numeric(input$LinfVal)-meanlength)/(meanlength-as.numeric(input$LcVal))
       Lmax<-as.numeric(input$LinfVal)*(1-exp(-as.numeric(input$KVal)*(as.numeric(input$MaxAgeVal)-as.numeric(input$t0Val))))
       meanZ<-optim(par=BHZ,CalcLbar,method="Brent",lower=0.001,upper=2*BHZ,Lc=as.numeric(input$LcVal),Lmax=Lmax,Linf=as.numeric(input$LinfVal),K=as.numeric(input$KVal),Lbar=meanlength)
       ZVect[i,5]<-meanZ$par
       BHZlow<-as.numeric(input$KVal)*(as.numeric(input$LinfVal)-(meanlength+1*SEmeanlength))/((meanlength+1*SEmeanlength)-as.numeric(input$LcVal))
       lowZ<-optim(par=BHZlow,CalcLbar,method="Brent",lower=0.001,upper=BHZlow,Lc=as.numeric(input$LcVal),Lmax=Lmax,Linf=as.numeric(input$LinfVal),K=as.numeric(input$KVal),Lbar=(meanlength+1*SEmeanlength))
       ZVect[i,4]<-lowZ$par
       BHZhigh<-as.numeric(input$KVal)*(as.numeric(input$LinfVal)-(meanlength-1*SEmeanlength))/((meanlength-1*SEmeanlength)-as.numeric(input$LcVal))
       highZ<-optim(par=BHZhigh,CalcLbar,method="Brent",lower=0.001,upper=BHZhigh,Lc=as.numeric(input$LcVal),Lmax=Lmax,Linf=as.numeric(input$LinfVal),K=as.numeric(input$KVal),Lbar=(meanlength-1*SEmeanlength))
       ZVect[i,6]<-highZ$par
       BHZlow95<-as.numeric(input$KVal)*(as.numeric(input$LinfVal)-(meanlength+1.96*SEmeanlength))/((meanlength+1.96*SEmeanlength)-as.numeric(input$LcVal))
       lowZ95<-optim(par=BHZlow95,CalcLbar,method="Brent",lower=0.001,upper=BHZlow95,Lc=as.numeric(input$LcVal),Lmax=Lmax,Linf=as.numeric(input$LinfVal),K=as.numeric(input$KVal),Lbar=(meanlength+1.96*SEmeanlength))
       ZVect[i,3]<-lowZ95$par
       BHZhigh95<-as.numeric(input$KVal)*(as.numeric(input$LinfVal)-(meanlength-1.96*SEmeanlength))/((meanlength-1.96*SEmeanlength)-as.numeric(input$LcVal))
       highZ95<-optim(par=BHZhigh95,CalcLbar,method="Brent",lower=0.001,upper=BHZhigh95,Lc=as.numeric(input$LcVal),Lmax=Lmax,Linf=as.numeric(input$LinfVal),K=as.numeric(input$KVal),Lbar=(meanlength-1.96*SEmeanlength))
       ZVect[i,7]<-highZ95$par
       ZBox$stats[,i]<-ZVect[i,3:7]
     }
     bxp(ZBox,ylim=c(0.95*min(ZBox$stats),1.05*max(ZBox$stats)),ylab="Mortality rate (Z)",xlab="Year",main=paste0("Annual mortality rate of ",input$SpeciesName),show.names = TRUE)
      }
      zBox()
   })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

