#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(gridExtra)
library(viridis)
require(cowplot)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
  }

plot_saliency_vs_distance <- function(df, graphname){
  textsize <- 7
  colorscale <- c(brewer.pal(12,"Set1"))
  p <- ggplot(data=df,aes(x=distance,y=saliency)) +
          #geom_hex(bins = 70) +
          #scale_fill_continuous(type = "viridis") +
          geom_point(pch=16, size=0.3, color='black', alpha=0.3) +
          #geom_smooth(method="loess", se=FALSE, fullrange=FALSE, level=0.95, color=colorscale[2], linewidth=0.3) +
          theme_cowplot(12) +
          theme(plot.background = element_rect(fill = "white"),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=textsize,face="bold"),
                legend.title=element_blank(),
                legend.key.height = unit(0.15, 'in'),
                legend.key.width = unit(0.05, 'in'),
                legend.text=element_text(size=textsize,face="bold"),
                legend.position='right') +
          ylim(0,0.9) +
          xlab("Distance (Å)") +
          ylab("Saliency score") +
          scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8))
  ggsave(graphname, p, width=1.8, height=1.8, dpi=600)
  }

df <- read_csv('data/distance_saliency_correlation.csv')
plot_saliency_vs_distance(df, "graph/saliency_vs_dist.png")
print (paste('Pearson correlation between distance and saliency score:', cor(df$distance, df$saliency, method='pearson')))
print (paste('Spearman correlation between distance and saliency score:', cor(df$distance, df$saliency, method='spearman')))
