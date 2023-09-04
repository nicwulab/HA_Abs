#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(data.table)
library(reshape)
library(stringr)
library(dplyr)
library(qualpalr)
library(readxl)
require(cowplot)

plot_gene_dist <- function(df_plot, c, gene){
  df_plot_V <- filter(df, saliency_cluster==c) %>%
                 separate(Heavy_V_gene, into=c('IGHV_gene'), sep='\\*', extra='drop')
  df_plot_D <- filter(df, saliency_cluster==c) %>%
                 filter(!grepl(',', Heavy_D_gene)) %>%
                 separate(Heavy_D_gene, into=c('IGHD_gene'), sep='\\*', extra='drop')
  if (gene == 'IGHV'){df_plot <- data.table(table(df_plot_V$IGHV_gene))}
  if (gene == 'IGHD'){df_plot <- data.table(table(df_plot_D$IGHD_gene))}
  textsize <- 7
  p <-  ggplot() +
          geom_bar(data=df_plot, aes(x=V1, y=N),
                   stat='identity', width=0.8, position=position_dodge(0.8), fill='gray50') +
          theme_cowplot(12) +
          theme(plot.title=element_text(size=textsize+1,face="bold",colour='black',hjust=0.5),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title.x=element_blank(),
                axis.title.y=element_text(size=textsize,face="bold"),
                legend.position = "right",
                legend.title    = element_blank(),
                legend.text=element_text(size=textsize,face="bold"),
                legend.justification='center',
                legend.key.size = unit(0.3,"line")) +
          ylab("count") +
          ggtitle(paste('Cluster ', c, sep=''))
  ggsave(paste('graph/', gene, 'count_cluster_', c, '.png', sep=''), p, width=length(df_plot$V1)/12+0.8, height=1.3, bg='white', dpi=600)
  }

df <- read_tsv('data/HA_Stem_saliency_clustered.tsv')
print ('cluster distribution:')
print (table(df$saliency_cluster))
for (c in unique(df$saliency_cluster)){
  gene <- 'IGHV'
  plot_gene_dist(df_plot, c, gene)
  gene <- 'IGHD'
  plot_gene_dist(df_plot, c, gene)
  }
