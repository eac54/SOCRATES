values<-read.csv('SupplementaryTable4.csv')

library(dplyr)
set.seed(13)
values$Quartile=ntile(values$Bliss.score,4)

q4dist=values$Rank[values$Quartile=="4"]
q1dist=values$Rank[values$Quartile=="1"]

wilcox.test(q1dist,q4dist) # Gives the p value of Mann-whitney U test for predictions

pvalues=c()

for(i in 1:10000){
  dist1=c(values[sample(nrow(values),25),]['Rank'])
  dist2=c(values[sample(nrow(values),23),]['Rank'])
  temp=wilcox.test(dist1$Rank,dist2$Rank)$p.value
  pvalues[i] <- temp
}

breakpoints=pretty(1:23,n=15)

q4hist=hist(Q4dist$Rank,plot=FALSE,breaks=breakpoints)
q1hist=hist(Q1dist$Rank,plot=FALSE,breaks=breakpoints)

c1<-"mediumpurple"
pvalue_hist=hist(pvalues,main='Distribution of p values (10000 random permutations)',col=c1)
abline(v=0.0038,col='brown1', lwd=3)
text(0.005, 500, "p value of \nour predictions",col='brown1',srt=90,cex=1.8)
