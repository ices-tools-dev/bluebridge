library(DLMtool)
library(FLCore)
setup()

#Read in the DLMtoolkit data object you want to extract data from.
myDataDLMtool<-new("Data","C:/Test/Noodle_Fish.csv")

#Create FLR stock object and fill with Catch, natural  mortality, and maturity data.
myDataFLStock<-FLStock()
myDataFLStock@name<-myDataDLMtool@Name
myDataFLStock@m<-FLQuant(myDataDLMtool@Mort)
myDataFLStock@mat<-FLQuant(c(0,0.5,0.95,1),dimnames=list(length=c(0,myDataDLMtool@L50,myDataDLMtool@L95,myDataDLMtool@vbLinf),year="all"))
myDataFLStock@catch<-FLQuant(as.vector(myDataDLMtool@Cat[1,]),units=myDataDLMtool@Units,dimnames=list(catch="total",year=myDataDLMtool@Year))

#Create FLR index object and fill with index of abundance data
myDataFLIndex<-FLIndex()
myDataFLIndex@name<-myDataDLMtool@Name
myDataFLIndex@index<-FLQuant(as.vector(myDataDLMtool@Ind[1,]),dimnames=list(Index="Abundance",year=myDataDLMtool@Year))
myDataFLIndex@index.var<-FLQuant(myDataDLMtool@CV_Ind[1],dimnames=list(Index="CV",year="all"))

#Create a variance covariance matrix for growth data and then create and fill an a4a growth object.
#Make an empty cor matrix
cm <- diag(c(1,1,1))
#Assume k and linf are negatively correlated while t0 is independent
cm[1,2] <- cm[2,1] <- -0.5
#Scale cor to var using cv values from DLMtoolkit data
cv <- c(myDataDLMtool@CV_vbLinf,myDataDLMtool@CV_vbK,myDataDLMtool@CV_vbt0)
p <- c(linf=myDataDLMtool@vbLinf, k=myDataDLMtool@vbK, t0=myDataDLMtool@vbt0)
vc <- matrix(1, ncol=3, nrow=3)
l <- vc
l[1,] <- l[,1] <- p[1]*cv[1]
k <- vc
k[,2] <- k[2,] <- p[2]*cv[2]
t <- vc
t[3,] <- t[,3] <- p[3]*cv[3]
mm <- t*k*l
diag(mm) <- diag(mm)^2
mm <- mm*cm
#Check that the intended correlation matrix was created
all.equal(cm, cov2cor(mm))

#Create the a4a growth object and fill with data for von Bertallanfy growth model
vbObj <- a4aGr(
  grMod=~linf*(1-exp(-k*(t-t0))),      
  grInvMod=~t0-1/k*log(1-len/linf),      
  params=FLPar(linf=myDataDLMtool@vbLinf, k=myDataDLMtool@vbK, t0=myDataDLMtool@vbt0, units=c('cm','year-1','year')),vcov=mm)

