#R code
library(ggplot2)
library(scales)
library(RColorBrewer)
library(readr)
library(tidyr)
library(reshape)
library(stringr)
library(gridExtra)
library(ggforce)
library(ggbeeswarm)
library(sinaplot)
library(readxl)
library(dplyr)
require(cowplot)

plot_model_perform <- function(data_table, graphname, xlab, ylab, xmin, xmax, ymin, ymax){
  textsize <- 8
  palette  <- c(brewer.pal(7,"Dark2"))
  p <- ggplot(data=data_table,aes(x=param_x,y=param_y,color=model)) +
          geom_line(size=0.3) +
          geom_point(size=0.5, pch=16) +
          scale_color_manual(values=palette,drop=FALSE) +
          theme_cowplot(12) +
          theme(plot.background = element_rect(fill = "white"),
                axis.text=element_text(size=textsize,face="bold",colour = 'black'),
                axis.text.x=element_text(angle=0,hjust=0.5,vjust=0.5,colour = 'black'),
                #axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,colour = 'black'),
                axis.text.y=element_text(hjust=0.5,vjust=0.5,colour = 'black'),
                axis.title=element_text(size=textsize,face="bold"),
                legend.title=element_blank(),
                legend.key.height = unit(0.15, 'in'),
                legend.key.width = unit(0.05, 'in'),
                legend.text=element_text(size=textsize,face="bold"),
                legend.position='right') +
          guides(colour = guide_legend(override.aes = list(size=0.5))) +
          xlim(xmin, xmax) +
          ylim(ymin, ymax) +
          ylab(ylab) +
          xlab(xlab)
  ggsave(graphname, p, width=4, height=1.5, dpi=600)
  }

plot_model_result <- function(data_table, graphname, ylab){
  textsize <- 7
  data_table <- filter(data_table, param_x==1)
  p <- ggplot(data=data_table,aes(x=model,y=param_y)) +
          geom_bar(stat="identity",position=position_dodge(), fill='grey30') +
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
          ylim(0, 1) +
          ylab(ylab) +
          xlab("")
  ggsave(graphname, p, width=1, height=1.5, dpi=600)
  }

exp1 <- read_excel('model/results/experiment1/repeated_with_Yiquan/results.xlsx') %>%
          mutate(model=str_replace(model,'original','Transform-6 (CDRs)')) %>%
          mutate(model=str_replace(model,'single','Transform-1 (CDRs)')) %>%
          mutate(model=str_replace(model,'dense','MLP (CDRs)')) %>%
          mutate(model=str_replace(model,'tree','Tree (CDRs)')) %>%
          rename(param_x=ratio)
exp2 <- read_excel('model/results/experiment2/repeated/results.xlsx') %>%
          mutate(model=str_replace(model,'single','Transform-1 (VH/VL)')) %>%
          mutate(model=str_replace(model,'dense','MLP (VH/VL)')) %>%
          mutate(model=str_replace(model,'tree','Tree (VH/VL)')) %>%
          rename(param_x=ratio)
exp_HAvsRBD <- rbind(exp1,exp2)
plot_model_perform(mutate(exp_HAvsRBD, param_y=prc), 'graph/model_HAvsRBD_PRAUC.png', 'sample size', 'PR AUC', 0, 1500, 0.4, 1)
plot_model_perform(mutate(exp_HAvsRBD, param_y=auc), 'graph/model_HAvsRBD_ROCAUC.png', 'sample size', 'ROC AUC', 0, 1500, 0.4, 1)

exp3 <- read_excel('model/results/experiment3/repeated/results.xlsx') %>%
          mutate(model=str_replace(model,'single','Transform-1')) %>%
          mutate(model=str_replace(model,'dense','MLP')) %>%
          mutate(model=str_replace(model,'tree','Tree')) %>%
          mutate(model=factor(model, levels=c('Transform-1','MLP','Tree'))) %>%
          rename(param_x=ratio)
exp4 <- read_excel('model/results/experiment4/repeated/results.xlsx') %>%
          mutate(model=str_replace(model,'single','Transform-1')) %>%
          mutate(model=str_replace(model,'dense','MLP')) %>%
          mutate(model=str_replace(model,'tree','Tree')) %>%
          mutate(model=factor(model, levels=c('Transform-1','MLP','Tree'))) %>%
          rename(param_x=ratio)
plot_model_perform(mutate(exp3, param_y=prc), 'graph/model_HAvsOthers_compressed_PRAUC.png', 'Others:HA ratio', 'PR AUC', 0.5, 7.5, 0.3, 0.9)
plot_model_perform(mutate(exp3, param_y=auc), 'graph/model_HAvsOthers_compressed_ROCAUC.png', 'Others:HA ratio', 'ROC AUC', 0.5, 7.5, 0.3, 0.9)
plot_model_perform(mutate(exp4, param_y=prc), 'graph/model_HAvsOthers_all_PRAUC.png', 'Others:HA ratio', 'PR AUC', 0.2, 3.8, 0.3, 0.9)
plot_model_perform(mutate(exp4, param_y=auc), 'graph/model_HAvsOthers_all_ROCAUC.png', 'Others:HA ratio', 'ROC AUC', 0.2, 3.8, 0.3, 0.9)

plot_model_result(mutate(exp4,param_y=prc), 'graph/model_HAvsOthers_compressed_PRAUC_NR1.png', 'PR AUC')
plot_model_result(mutate(exp4,param_y=auc), 'graph/model_HAvsOthers_compressed_ROCAUC_NR1.png', 'ROC AUC')
