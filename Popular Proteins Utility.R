##### histogram of publication counts
setwd("~/Documents/Graphing Temp")
Data <-read.table("pubcount.txt",header=TRUE)
pdf(file="popbarplot.pdf")
par(mfrow=c(2,3))
attach(Data)
barplot(Brain, ylim=c(0,500),main="Brain")
barplot(Heart, ylim=c(0,500),main="Heart")
barplot(Liver, ylim=c(0,500),main="Liver")
barplot(Kidney, ylim=c(0,500),main="Kidney")
barplot(Lung, ylim=c(0,500),main="Lung")
barplot(Gut, ylim=c(0,500),main="Gut")
dev.off()

#### By how much do publication counts drop?
a <- sum(Data$Brain[51:100])/sum(Data$Brain[1:50])
b <- sum(Data$Heart[51:100])/sum(Data$Heart[1:50])
c <- sum(Data$Kidney[51:100])/sum(Data$Kidney[1:50])
d <- sum(Data$Gut[51:100])/sum(Data$Gut[1:50])
e <- sum(Data$Lung[51:100])/sum(Data$Lung[1:50])
f <- sum(Data$Liver[51:100])/sum(Data$Liver[1:50])
(a+b+c+d+e+f)/6


### similarity vs frequency

setwd("~/Documents/Graphing Temp")
Data <-read.table("simvspub.txt",header=TRUE)
attach(Data)
pdf(file="semanticsimilarity.pdf")
par(mfrow=c(2,3))
plot(Data$brain.sim,log10(Data$brain.pub),main="brain",pch=1,cex=.5,ylim=c(0.8,3.2),xlim=c(0.5,1.5))
text(Data$brain.sim[1:5],log10(Data$brain.pub[1:5]),label=Data$brain.gn[1:5],pos=4,col="red",cex=0.5)
plot(Data$heart.sim,log10(Data$heart.pub),main="heart",pch=1,cex=.5,ylim=c(0.8,3.2),xlim=c(0.5,1.5))
text(Data$heart.sim[1:5],log10(Data$heart.pub[1:5]),label=Data$heart.gn[1:5],pos=4,col="red",cex=0.5)
plot(Data$kidney.sim,log10(Data$kidney.pub),main="kidney",pch=1,cex=.5,ylim=c(0.8,3.2),xlim=c(0.5,1.5))
text(Data$kidney.sim[1:5],log10(Data$kidney.pub[1:5]),label=Data$kidney.gn[1:5],pos=4,col="red",cex=0.5)
plot(Data$gut.sim,log10(Data$gut.pub),main="gut",pch=1,cex=.5,ylim=c(0.8,3.2),xlim=c(0.5,1.5))
text(Data$gut.sim[1:5],log10(Data$gut.pub[1:5]),label=Data$gut.gn[1:5],pos=4,col="red",cex=0.5)
plot(Data$liver.sim,log10(Data$liver.pub),main="liver",pch=1,cex=.5,ylim=c(0.8,3.2),xlim=c(0.5,1.5))
text(Data$liver.sim[1:5],log10(Data$liver.pub[1:5]),label=Data$liver.gn[1:5],pos=4,col="red",cex=0.5)
plot(Data$lung.sim,log10(Data$lung.pub),main="lung",pch=1,cex=.5,ylim=c(0.8,3.2),xlim=c(0.5,1.5))
text(Data$lung.sim[1:5],log10(Data$lung.pub[1:5]),label=Data$lung.gn[1:5],pos=4,col="red",cex=0.5)


dev.off()


### human vs. mouse ranks
Data <-read.table("hvsm.txt",header=TRUE)
attach(Data)
cor.test(brain.h.sim,brain.m.sim,method=c("spearman"))
cor.test(heart.h.sim,heart.m.sim,method=c("spearman"))
cor.test(liver.h.sim,liver.m.sim,method=c("spearman"))
cor.test(kidney.h.sim,kidney.k.sim,method=c("spearman"))
cor.test(gut.h.sim,gut.m.sim,method=c("spearman"))
cor.test(lung.h.sim,lung.m.sim,method=c("spearman"))
