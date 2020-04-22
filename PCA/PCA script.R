library(pca3d)
data_read=read.table("Merged_Fold_change_OnlySignificant_File.txt",sep="\t",header=TRUE)
data_read_df=data_read[,4:ncol(data_read)]
data_read_df.pca=prcomp(data_read_df,center=T,scale.=T)
type_group=factor(data_read_df$SampleType)
pca3d(data_read_df.pca, group= type_group)

snapshotPCA3d(file="PCA3d.png")
