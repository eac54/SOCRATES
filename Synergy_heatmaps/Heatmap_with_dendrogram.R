library(superheat)
data<-read.csv("Bliss_values.csv")
rownames(data)<-data[,1]
data<-data[, -1]

superheat(data,
          pretty.order.rows = TRUE,
          row.dendrogram = TRUE,
          pretty.order.cols = TRUE,
          col.dendrogram = TRUE,
          left.label.text.size = 2.9,
          bottom.label.text.size = 3.5,
          grid.hline.col = "white",
          grid.vline.col = "white",
          heat.pal = c("white", "#542788"),
          heat.pal.values = c(0, 1),
          row.title = "Combination",
          column.title = "Cell line",
          title = "Bliss scores",
          scale=FALSE,
          legend.height = 0.25,
          legend.width = 2,
          legend.text.size = 10)