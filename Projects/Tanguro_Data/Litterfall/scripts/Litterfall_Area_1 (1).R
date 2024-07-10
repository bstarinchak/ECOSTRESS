###########################
#Author: Maracahipes-Santos, Leonardo; Silverio, Divino Vicente; #Instituto de Pesquisa Ambiental da Amazonia (IPAM)
#Junho 09, 2017
#?ltima atualiza??o Fevereiro 14, 2018
#E-mail: leonardo.maracahipes@ipam.org.br 
#        maracahipesbio@gmail.com

### Script para master de dados de liteira da area 1
### Pacotes
library(plyr)
library(data.table) #1.9.5+
library(tidyverse)
library(ggplot2)
# install.packages("plyr", dependencies = TRUE)

## diretorio Maracahipes-Santos - pasta com dados


path_rd = "/ECOSTRESS/Projects/Tanguro_Data/Litterfall/raw_data"
path_pd = "/ECOSTRESS/Projects/Tanguro_Data/Literfall/proc_data"
setwd("C:/Users/leonardo.maracahipes/Dropbox (DadosTanguro)/trabalho/projetos/Savanizacao/Liteira_AR1/dados_atuais")
dir()

#####
#####

all_1 = read.csv("1_Master_Liteira_Area_1_Abr2024_For_Bela.csv",h=T)

#####
#####

plot1=all_1 %>%
  #  filter(var == variables[NN]) %>%
  #  mutate(COLOR = ifelse(FR_CR_Dif > 0, "red", "blue")) %>%
  na.omit() %>%
  ggplot(aes(x = weight_1, 
             fill = plot, color = plot)) +
  geom_density(alpha = .3) +
  scale_fill_manual(values = c("darkgreen", "orange","brown")) +
  scale_color_manual(values = c("darkgreen", "orange","brown")) +
  facet_wrap(~years, scales = "free", ncol = 4) +
  theme_light();plot1

#####
#####
library(ggplot2)
library(gridExtra)
library(extrafont)

#####
#####

pdf("Additional_documents_Figure_litter_AR1.pdf",width=9,height=18)

grid.arrange(plot1, ncol=1)

dev.off()

###
rm(list=ls())

### The end!