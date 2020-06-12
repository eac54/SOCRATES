valuestoplot<-read.csv("TableforR_BlissvsRanking.txt",sep="\t",header=TRUE,stringsAsFactors = FALSE)
library(magrittr)
library(dplyr)
valuestoplot=valuestoplot %>%
  mutate(RankbyQuartile = case_when(Rank <= 8 ~ 'Q1',
                                    Rank >= 22 ~ 'Q4',
                                    TRUE ~ 'Q2 or Q3'))

Q1dist=c(valuestoplot[valuestoplot$Bliss<=0.06,])

Q4dist=c(valuestoplot[valuestoplot$Bliss>=0.32,])

ks.test(Q1dist$Rank,Q4dist$Rank,alternative="less")

pvalues<- vector()

for(i in 1:10000){
  dist1=c(valuestoplot[sample(nrow(valuestoplot),25),]['Rank'])
  dist2=c(valuestoplot[sample(nrow(valuestoplot),23),]['Rank'])
  temp=ks.test(dist1$Rank,dist2$Rank,alternative="less")$p.value
  pvalues[i] <- temp
}

hist(pvalues)


breakpoints=pretty(1:23,n=15)

q4hist=hist(Q4dist$Rank,plot=FALSE,breaks=breakpoints)
q1hist=hist(Q1dist$Rank,plot=FALSE,breaks=breakpoints)
c1<-"lightblue1"

#tiff("Plot_Bliss_vs_Prediction_Distof10000randompermutations.tiff", width = 8, height = 4, units = 'in', res = 300, compression = 'none')
pvalue_hist=hist(pvalues,main='Distribution of P values (10000 random permutations)',col=c1)
abline(v=0.0243,col='blue', lwd=3)
text(0.005, 1000, "p value of our predictions",col='blue',srt=90,cex=1.8)
#dev.off()