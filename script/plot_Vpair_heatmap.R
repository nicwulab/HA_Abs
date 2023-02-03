# Title     : to plot basic statistics of Ab table
# Created by: yiquan
# Created on: 10/22/21

library(ggplot2)
library(readxl)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(dplyr)
library(ggrepel)
library(gridExtra)
library(data.table)
library(GGally)
library(e1071)
library(ggforce)
library(ggbeeswarm)
library(ggExtra)
require(cowplot)
require(ggseqlogo)

plot_point_heatmap <- function (df,path,x_lab,y_lab,title){
  p  <- ggplot(df, aes(x=LV,y=HV,fill=epitope,size=abs(freq))) +
          geom_point(alpha = 0.5,color='black',pch=21) +
          theme_classic() + 
          scale_fill_brewer(palette = "Set2") +
          scale_size_continuous(name = "Freq (%)",
                                range = c(0.2, 3.8)) +
          theme(plot.title=element_text(size=7,face="bold",hjust = 0.5),
                text = element_text(size=7,face="bold"),
                legend.position = 'right',
                legend.key.size = unit(1, 'lines'),
                legend.background=element_blank(),
                legend.text=element_text(size=7,face="bold",color='black'),
                legend.spacing.y = unit(0.02, 'lines'),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                axis.text.y=element_text(size=5,face="bold",color='black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=5,face="bold",color='black')) +
          ggtitle(title) +
          xlab(x_lab) +
          ylab(y_lab) +
          labs(fill='Epitope',size='Freq')
  ggsave(path,p,height = 4,width = 6,bg='white')
  }

df <- read_tsv('result/IGV_pair_freq.tsv') %>%
        mutate(HV = mapply(function(s){return (str_split(s, '__')[[1]][1])}, gene)) %>%
        mutate(LV = mapply(function(s){return (str_split(s, '__')[[1]][2])}, gene)) 
plot_point_heatmap(df,'graph/Vpair_heatmap.png','IGK(L)V','IGHV','')
