library(ape)
library(RColorBrewer)
library(dendextend)

cont_sens_Allphosphochanges=read.csv(‘AllData_includingNonsignificant_ContinuousSens_ReadyforAnalysis.txt’,sep=“\t”)
only_changes=cont_sens_Allphosphochanges_data[,5:]

plot(as.phylo(hclust(dist(scale(t(only_changes)), method = "euclidean"), method = "ward.D2")), type = "fan", tip.color=brewer.pal(8, "Dark2")[cutree(hc, h=heights_per_k.dendrogram(dend)["8"])],cex=0.65,font=4)