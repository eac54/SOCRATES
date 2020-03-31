library(ggplot2)
data<-read.csv("Lisa_data_bliss_option1 copy.csv")
ggplot(data, aes(Combinations, Cell.line)) +
  ggtitle('Bliss scores') +
  theme_bw() +
  xlab('Combination') +
  ylab('Cell line') +
  geom_tile(aes(fill = Bliss.score), color='white') +
  scale_fill_gradient(low = 'white', high = 'darkblue', space = 'Lab') +
  theme(axis.text.x=element_text(angle=90),
        axis.ticks=element_blank(),
        axis.line=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_line(color='#eeeeee'))

