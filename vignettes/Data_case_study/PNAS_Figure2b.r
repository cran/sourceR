try(dev.off(),T)
#setwd("C:/Documents and Settings/jcmarsha/My Documents/Massey/SVN/Other/Petra/attributiongraph")
#PNAS_Figure2 <- read.table("PNAS_Figure2.csv", header=TRUE, sep = ",")
PNAS_Figure2 <- read.table (file.choose (), header=TRUE, sep = ",")
# pdf ("PNAS_Figure2b.pdf", 7, 5)
X11()
d<-rbind(PNAS_Figure2$Dutch,PNAS_Figure2$Danish, PNAS_Figure2$Island)
cis<-rbind(PNAS_Figure2$Lower,PNAS_Figure2$Upper,PNAS_Figure2$DLower,PNAS_Figure2$DUpper, PNAS_Figure2$DULower,PNAS_Figure2$DUUpper)
#cols<-c("yellow","red","blue","green")
cols<-c("#FFFF00FF", "#FFaF00FF", "#FF7F00FF", "#FF0000FF", "#0000FFFF", "#00FF00FF")
mp<-barplot (d, ylab="Proportion of cases attributed", ylim=c(0, 1), names=c("Poultry\nSupplier A", "Poultry\nSupplier B", "Poultry\nSupplier C", "Bovine\n", "Ovine\n", "Environment\n"), beside=TRUE, cex.lab=1, cex.names=0.9, cex.axis=0.85,col=rbind(cols,cols,cols), mgp=c(3,2,0))
# mp is the midpoints (on the x axis) of the bars.
segments(mp[1,],cis[1,],mp[1,],cis[2,])# Confidence intervals for the left bars
segments(mp[2,],cis[3,],mp[2,],cis[4,])# Confidence intervals for the middle bars
segments(mp[3,],cis[5,],mp[3,],cis[6,])# Confidence intervals for the right bars
for(k in 1:6)
{
  mtext("I",1,line=0,at=mp[1,k],cex=0.8)
  mtext("II",1,line=0,at=mp[2,k],cex=0.8)
  mtext("III",1,line=0,at=mp[3,k],cex=0.8)
}

dev.off()

#