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

plot_IG_usage <- function(data,w,h,graphname){
  palette  <- c(brewer.pal(3,"Set2"))
  textsize <- 7
  width_dodge <- 0.75
  p <-  ggplot() +
          geom_bar(data=data, aes(x=gene, y=freq*100, group=epitope, fill=epitope), 
                   stat='identity', width=width_dodge, position=position_dodge(width_dodge)) +
          scale_fill_manual(values=palette,drop=FALSE) +
          geom_errorbar(data=data,aes(x=gene, group=epitope, ymin=(freq-SE)*100, ymax=(freq+SE)*100),
                        size=0.25, width=0,position=position_dodge(width_dodge), color='gray30') +
          theme_cowplot(12) +
          theme(axis.text=element_text(size=textsize-1,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_blank(),
                axis.title.y=element_text(size=textsize,face="bold"),
                axis.line = element_line(colour = 'black', size = 0),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                legend.position = "top",
                legend.title    = element_blank(),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.size = unit(0.3,"line")) +
          ylab("frequency (%)")
  ggsave(graphname,p,width=w, height=h, bg='white', dpi=1200)
}

exclude <- function(list,names){
     ## return the elements of the list not belonging to names
     member..names <- names(list)
     index <- which(!(member..names %in% names))
     list[index]    
     }  

set.seed(5)
IGHV_data <- read_tsv('result/IGHV_freq.tsv')
plot_IG_usage(IGHV_data, 6.5, 1.3, 'graph/IGHV_usage.png')

IGLV_data <- read_tsv('result/IGLV_freq.tsv') 
IGLV_level <- unique(IGLV_data$gene)
IGLV_level <- c(IGLV_level[-unlist(lapply(c('IGLV10-54'), function(x) grep(x, IGLV_level)))], 'IGLV10-54')
IGLV_data <- mutate(IGLV_data, gene=factor(gene,levels=IGLV_level))
plot_IG_usage(IGLV_data, 6.5, 1.3, 'graph/IGLV_usage.png')

IGHD_data <- read_tsv('result/IGHD_freq.tsv')
plot_IG_usage(IGHD_data, 4.2, 1.3, 'graph/IGHD_usage.png')

IGHJ_data <- read_tsv('result/IGHJ_freq.tsv')
plot_IG_usage(IGHJ_data, 2, 1.16, 'graph/IGHJ_usage.png')

IGHJ_data <- read_tsv('result/IGLJ_freq.tsv')
plot_IG_usage(IGHJ_data, 2.7, 1.16, 'graph/IGLJ_usage.png')
