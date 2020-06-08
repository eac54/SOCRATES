data<-read.table("TableforR_BlissvsRanking.txt", sep="\t", header = TRUE)

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adjusted R^2 = ",signif(summary(fit)$adj.r.squared, 5)), xlab= "Bliss score", ylab= "EPS rank")
}


ggplotRegression(lm(Rank ~ Bliss, data = data))
