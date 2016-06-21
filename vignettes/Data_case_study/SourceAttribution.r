setwd('D:/SourceAttribution')
winbugs<-'C:/Program Files/WinBUGS1.4.3'# The location of WinBugs
library(R2WinBUGS)#You may need to download this R package
# Set up prevalence MCMC parameters
pburnin<-2000
piters<-10000
pchains<-5
pthin<-10
# Import data
TSources<-7
TTypes<-73
alldata<-matrix(0,TSources,TTypes)
# If you put your data in like this (ie without the 1) then you can easily spot mistakes and change the prior/data if you want to!
alldata[1,]<-c(0, 1, 0, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 2, 0,  1, 0, 0, 1, 0, 0, 0, 0, 0, 4, 0, 0, 1, 2, 0, 0, 0, 1, 1, 1, 0, 3, 0, 0, 1, 1, 7, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 2, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1)#Water
alldata[2,]<-c(0, 2, 0, 4,33, 1, 5, 5,10, 0, 0, 7,10, 1, 0, 0, 0, 23, 5, 8, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)#Tegel
alldata[3,]<-c(1, 0, 0, 0,15,26,17, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0,  2, 0, 0, 1, 6, 5, 0,11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,15, 0, 0, 1, 0, 0, 0, 0)#Inghams
alldata[4,]<-c(0, 2, 0, 0,23, 7,13, 1, 4, 0, 0, 0, 0, 0, 0, 0,20,  2, 2, 4, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)#Other Poultry
alldata[5,]<-c(0, 0, 0, 0,12, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,  0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)#Wild Bird
alldata[6,]<-c(3, 0, 1,11, 2, 0, 7, 0,18, 4, 0, 4, 0, 0, 2, 2, 0,  3, 5, 0, 0, 0, 0, 4, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)#Cow
alldata[7,]<-c(1, 0, 0,15, 2, 0,18, 0, 1, 9, 0, 4, 0, 0,16, 3, 0,  2, 2, 1, 1, 7, 0,15, 0, 0, 4, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0)#Sheep
human <-     c(4, 1,11,12,35,39,19,15,17, 7, 1,19,11,19, 2, 3, 7,118, 5, 7, 3, 4, 1, 7, 2, 1, 3, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)#Human
#              1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17  18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 
allsources<-c("Water","Tegel","Inghams","Other Poultry","Wild Bird","Bovine origin","Origin")
alltypes<-c("21", "25","38", "42", "45", "48", "50", "52", "53", "61", "137", "190", "257", "354", "422", "436", "451", "474", "520", "583", "677", "1517", "1581", "2026", "2345", "3072", "NEW", "u21d", "u354", "u42", "177", "227", "233", "393", "526", "694", "1030", "1191", "1225", "1911", "2309", "2347", "2354", "2381", "2392", "2397", "2535", "2584", "2619", "3093", "3230", "3232", "3301", "u1275", "u1275b", "u2", "u21/48/206", "u22", "u2389", "u257", "u2620", "u3", "u4", "u403", "u45", "u48", "u5", "u61", "u692", "u692b", "u692c", "u7", "u8")
totSamples<-c(332,239,196,127,192,595,552)# Number of samples  tested  for c. jejuni, for each source.
posSamples<-c(62,181,113,109,24,97,165)   # Number of samples positive for c. jejuni, for each source.
# This function puts together tighter datasets from the complete data.
# sources should be a list of sources to use eg sources=list(1,c(2,3)) uses only 2 sources: water and Tegel and Inghams combined
# sourcelist summarises your choices!
# types should be a vector of single types to be included in the model eg 1:30 includes them all.
# combineTypes should be a list of *additional* types to combine
# EG createData(list(1,2,3,4,5,6,7),1:30,list(31:40,41:73))
# uses 32 types.
createData <- function(sources,types,combineTypes=list()) {
  NSources<<-length(sources)
  NTypes<<-length(types)
  Pdata<<-matrix(0,NSources,TTypes)
  Pos<<-rep(0,NSources)
  Tot<-rep(0,NSources)
  sourcelist<<-rep("",NSources)
  for (i in 1:NSources) {
    for (j in sources[[i]]) {
      Pdata[i,]<<-Pdata[i,]+alldata[j,]
      Pos[i]<<-Pos[i]+posSamples[j]
      Tot[i]<- Tot[i]+totSamples[j]
      sourcelist[i]<<-paste(sourcelist[i],allsources[j],sep=" ")
    }
  }
  Neg<<-Tot-Pos
  # now deal with types
  Pdat<<-Pdata[,types]
  typelist<<-alltypes[types]
  o<<-human[types]
  if (length(combineTypes)>0) {
    for (i in 1:length(combineTypes)) {
      NTypes<<-NTypes+1
      Pdat<<-cbind(Pdat,apply(Pdata[,combineTypes[[i]]],1,sum))
      typelist<<-c(typelist,paste(alltypes[combineTypes[[i]]]))
      o<<-c(o,sum(human[combineTypes[[i]]]))
    }
  }
}
calculatePrevalence <- function(dirprior=1,betaprior=c(1,1),filename="PrevalenceLog.odc") {
  # Add priors
  Pdat<<-dirprior+Pdat
  Pos<<-betaprior[1]+Pos
  Neg<<-betaprior[2]+Neg
  # Prepare for WinBugs for prevalence model
  dat<-list("NSources","NTypes","Pdat","Pos","Neg")
  parameters<-list("Pi","p")
  # Call WinBugs debug=FALSE means that it closes Winbugs and moves on without waiting for you to see the output.
  output<-bugs(dat,inits=NULL,parameters,"PrevalenceModel.odc",n.thin=pthin,n.burnin=pburnin,n.chains=pchains,n.iter=piters,bugs.directory=winbugs,debug=T,DIC=FALSE)
  file.rename("log.odc",filename)
  #
  # Do the method of moments on the output
  #
  pa<<-matrix(0,NSources,NTypes)
  pb<<-matrix(0,NSources,NTypes)
  for (i in 1:NSources) {
    for (j in 1:NTypes) {
      m<-output$mean$p[i,j]
      v<-output$sd$p[i,j]^2
      aaddb<-m*(1-m)/v-1
      if (aaddb*m<1) {
        pa[i,j]<<-1
        pb[i,j]<<-1/m-1
      } else {
        pa[i,j]<<-aaddb*m
        pb[i,j]<<-aaddb*(1-m)
      }
    }
  }
  return(output)
}
calculateAttribution <- function(filename="AttributionLog.odc") {
  # Prepare for WinBugs for attribution model
  dat<-list("NSources","NTypes","pa","pb","o","atau","btau","aPrior")
  parameters<-list("source","proportion","lambda","a","q","tau","chicken")
  i1<-list(a=rep(10,NSources),tau=0.12)
  i2<-list(a=rep(1,NSources),tau=0.2)
  i3<-list(a=rep(50,NSources),tau=0.05)
  i4<-list(a=rep(2,NSources),tau=0.1)
  i5<-list(a=rep(25,NSources),tau=0.2)
  initials<-list(i1,i2,i3,i4,i5)
  # call winbugs and return the results
  final<-bugs(dat,inits=initials,parameters,"AttributionModel.odc",n.thin=thin,n.burnin=burnin,n.chains=5,n.iter=iters,bugs.directory=winbugs,debug=T)
  file.rename("log.odc",filename)
  return(final)
}
# Priors for attribution model
atau<-0.01
btau<-0.01
aPrior<-0.002
# MCMC parameters for the attribution model
burnin<-200
iters<-1200
thin<-10
#
#
# THIS BIT DOES THE WORK!!!!!!
#
#
#createData(list(1,2,3,4,5,6,7),1:30)# environment split up in water and WB
#calculatePrevalence()
#HumanTypesAnalysis<-calculateAttribution("HumanTypesAnalysis.odc")
#
#createData(list(c(1,5),c(2,3,4),6,7),1:30)# combine water with wild bird and all chicken sources together
#calculatePrevalence()
#AllTog<-calculateAttribution("AllTog.odc")
#
#createData(list(c(1,5),2,3,4,6,7),1:30)# combine water with wild bird 
#calculatePrevalence()
#EnvTog<-calculateAttribution("EnvTog.odc")
#
createData(list(1,2,3,4,5,6,7),1:73)# combine water with wild bird and add non-human ST's individually
calculatePrevalence()
ManyExtraST<-calculateAttribution("TestST.odc")
#
#createData(list(1,2,3,4,5,6,7),1:30,list(31:73))# add non-human ST's as an 'other' type.
#calculatePrevalence()
#SingleExtraST<-calculateAttribution("SingleExtraST.odc")
#
#aPrior<-0.01
#
#createData(list(1,2,3,4,5,6,7),1:30)# environment split up in water and WB
#calculatePrevalence()
#HumanTypesAnalysis2<-calculateAttribution("HumanTypesAnalysis2.odc")
#
#createData(list(c(1,5),c(2,3,4),6,7),1:30)# combine water with wild bird and all chicken sources together
#calculatePrevalence()
#AllTog2<-calculateAttribution("AllTog2.odc")
#
#createData(list(c(1,5),2,3,4,6,7),1:30)# combine water with wild bird 
#calculatePrevalence()
#EnvTog2<-calculateAttribution("EnvTog2.odc")
#
#createData(list(1,2,3,4,5,6,7),1:73)# combine water with wild bird and add non-human ST's individually
#calculatePrevalence()
#ManyExtraST2<-calculateAttribution("ManyExtraST2.odc")
#
#createData(list(1,2,3,4,5,6,7),1:30,list(31:73))# add non-human ST's as an 'other' type.
#calculatePrevalence()
#SingleExtraST2<-calculateAttribution("SingleExtraST2.odc")
#
# Now type 'ManyExtraST$summary' or more simply 'ManyExtraST$median$proportion' into R to see the results.