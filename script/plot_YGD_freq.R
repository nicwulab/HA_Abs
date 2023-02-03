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
library(ggforce)
library(ggbeeswarm)
library(sinaplot)
library(readxl)
require(cowplot)

plot_motif_usage <- function(data,w,h,graphname){
  textsize <- 7
  width_dodge <- 0.8
  p <-  ggplot() +
          geom_bar(data=data, aes(x=Ab_type, y=freq*100), 
                   stat='identity', width=width_dodge, position=position_dodge(width_dodge), fill='gray50') +
          geom_errorbar(data=data,aes(x=Ab_type, ymin=(freq-SE)*100, ymax=(freq+SE)*100),
                        size=0.25, width=0,position=position_dodge(width_dodge), color='gray10') +
          theme_cowplot(12) +
          theme(plot.title=element_text(size=textsize,face="bold",colour = 'black',hjust = 0.5),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=45,hjust=1,vjust=1,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_blank(),
                axis.title.y=element_text(size=textsize,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                legend.position = "right",
                legend.title    = element_blank(),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.size = unit(0.3,"line")) +
          ylab("frequency (%)") +
          ggtitle('CDR H3 with YGD motif')
  ggsave(graphname,p,width=w, height=h, bg='white', dpi=1200)
}

exclude <- function(list,names){
     ## return the elements of the list not belonging to names
     member..names <- names(list)
     index <- which(!(member..names %in% names))
     list[index]    
     }  

set.seed(5)
df <- read_tsv('result/YGD_motif_freq.tsv')
plot_motif_usage(df, 1.5, 2, 'graph/YGD_motif_freq.png')
