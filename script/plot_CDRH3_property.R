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

plot_CDRH3_property <- function(data_table, graphname, ylab){
  p_value_1 <- t.test(filter(data_table,epi=='Stem')$parameter, filter(data_table,epi=='Head')$parameter)$p.value
  p_value_2 <- t.test(filter(data_table,epi=='Stem')$parameter, filter(data_table,epi=='GenBank')$parameter)$p.value
  p_value_3 <- t.test(filter(data_table,epi=='Head')$parameter, filter(data_table,epi=='GenBank')$parameter)$p.value
  print (ylab)
  print (paste('Stem vs Head:', p_value_1))
  print (paste('Stem vs GenBank:', p_value_2))
  print (paste('Head vs GenBank:', p_value_3))
  baseline   <- mean(filter(data_table, epi=='GenBank')$parameter)
  data_table <- filter(data_table, epi!='GenBank')
  textsize <- 7
  p <- ggplot() +
          geom_boxplot(data=data_table,aes(x=epi,y=parameter),width=0.7, outlier.shape=NA, size=0.3) +
          geom_sina(data=data_table,aes(x=epi,y=parameter),
                    pch=16, size=0.1,method="counts", bin_limit=0.4, scale="width", maxwidth=0.5, color='gray40', alpha=0.4) +
          geom_hline(yintercept=baseline, lty="11") +
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
          guides(colour = guide_legend(override.aes = list(size=0.5))) +
          xlab("") +
          ylab(ylab)
  ggsave(graphname, p, width=1, height=2, dpi=600)
  }

df <- read_tsv('result/CDRH3_property.tsv') %>%
        filter(epi!='unknown')

df_select <- mutate(df, parameter=len_CDRH3)
plot_CDRH3_property(df_select, 'graph/CDRH3_length.png', 'CDR H3 length')

df_select <- mutate(df, parameter=CDRH3_H_score)
plot_CDRH3_property(df_select, 'graph/CDRH3_hydrophobicity.png', 'CDR H3 hydrophobic score')

#df_select <- mutate(df, parameter=tip_C_score)
#plot_CDRH3_property(df_select, 'graph/CDRH3_charge.png', 'CDR H3 charge score')

#df_select <- mutate(df, parameter=tip_A_score)
#plot_CDRH3_property(df_select, 'graph/CDRH3_aromatic.png', 'CDR H3 aromatic score')

df_select <- mutate(df, parameter=tip_H_score)
plot_CDRH3_property(df_select, 'graph/tip_hydrophobicity.png', 'CDR H3 tip hydrophobic score')

#df_select <- mutate(df, parameter=tip_C_score)
#plot_CDRH3_property(df_select, 'graph/tip_charge.png', 'CDR H3 tip charge score')

#df_select <- mutate(df, parameter=tip_A_score)
#plot_CDRH3_property(df_select, 'graph/tip_aromatic.png', 'CDR H3 tip aromatic score')
