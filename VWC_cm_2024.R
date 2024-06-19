
library(tidyverse)
library(Hmisc)
library(akima)
library(Hmisc)
library(pracma)

setwd("C:/Users/antonio.silva/DadosTanguro Dropbox/Antônio Carlos Silveiro da Silva/projetos/Savanizacao/TDR_pocos")


tdr.cm = read.csv("master_TDR_vwc_cm_2023.csv")


tdr.cm2 = tdr.cm%>%
  filter(site %in%c('G26','P25','H19','C02','M08','P02', 'K04'))%>%
  mutate(data=as.Date(date),
         ano = year(data),
         mes = month(data),
         dia = day(data),
         hour = substr(date,11,19),
         treat = ifelse(site%in%c('G26','P25'),'B1yr', 
                        ifelse(site%in%c('C02','M08','P02', 'K04'),'Control','B3Yr')),
         hydro_yr = ifelse(mes > 9, ano + 1, ano),
         Drought = ifelse(hydro_yr %in% 2016, 'D' , 'N')) %>%ungroup()


tdr.cm2$vwc000cm = as.numeric(as.character(tdr.cm2$vwc000cm))
names(tdr.cm2)

vwc = tdr.cm2 %>% 
  dplyr::  select(date:vwc900cm, total:Drought, ano, mes) %>%  
  mutate(#ano = year(date_wk),
    #mes = month(date_wk),
    data2 = ymd(paste(ano, mes, 15, sep = "_")),
    VWC_8m = rowSums(.[3:13]),
    VWC_6_8m = rowSums(.[12:13]),
    VWC_4_6m = rowSums(.[10:11]),
    VWC_2_4m = rowSums(.[8:9]),
    VWC_7m = rowSums(.[3:12]),
    VWC_6m = rowSums(.[3:11]),
    VWC_5m = rowSums(.[3:10]), 
    VWC_4m = rowSums(.[3:9]),
    VWC_3m = rowSums(.[3:8]),
    VWC_0_2m = rowSums(.[3:7]),
    VWC_2m = rowSums(.[3:7])) %>%
  
  #mutate(data2 = ymd(paste(ano, mes, 15, sep = "_"))) %>% 
  group_by(data2, treat) %>% 
  #summarise_at(vars("VWC_8m", "VWC_6_8m", "VWC_4_6m", "VWC_2_4m", "VWC_4m", "VWC_2m"), mean,na.rm = TRUE) %>% 
  summarise_at(vars("VWC_8m", "VWC_4m", "VWC_2m"), mean,na.rm = TRUE) %>% 
  gather(var, vwc, -data2, -treat) 


vwc = vwc %>% 
  mutate(var2 = ifelse(var%in%c("VWC_8m"), "8 m",
                       ifelse(var %in%c("VWC_4m"), "4 m", '2 m')))
         


vwc %>% 
  filter(!treat == 'Soy') %>% 
  #filter(var %in%c('VWC_8m')) %>% 
  group_by(data2, treat, var2) %>% 
  summarise(vwc_m = smean.cl.boot(vwc, na.rm = T)[1],
            vwc_l = smean.cl.boot(vwc, na.rm = T)[2],
            vwc_h = smean.cl.boot(vwc, na.rm = T)[3]) %>% 
  ggplot(aes(x = data2, y = vwc_m, color =treat))+
  geom_point(size = 4)+
  geom_line()+
  geom_errorbar(aes(ymin = vwc_l,
                    ymax = vwc_h), width = 0.2)+
  theme_light(base_size = 20)+
  scale_color_manual( values = c("darkred","darkorange",   "darkgreen"))+
  facet_wrap(~var2, ncol = 1, strip.position = "right", scales = 'free')+
  annotate("rect", xmin = ymd('2015-10-01'), xmax = ymd('2016-09-30'), 
           ymin = -Inf, ymax = Inf, fill = "#fdae6b", alpha = 0.2) +
  guides(colour = guide_legend(override.aes = list(size=6)))+ # aumenta a legenda
  labs(x = 'Ano', y = 'VWC (cm)')+
  facet_grid(~var2)+
  theme(legend.position= c(.06,.33),
        legend.title =    element_blank(),
        axis.title.x = element_blank(),
        
        #text = element_text(family = "Times New Roman", colour = "black", size = 26),
        axis.text = element_text(family = "Times New Roman",  colour = "black", size = 26),
        text = element_text(family = "Times New Roman", size = 26),
        strip.text = element_text( size=22),
        # axis.text.x=element_blank(),
        # axis.ticks.x=element_blank(),
        legend.background = element_blank())
