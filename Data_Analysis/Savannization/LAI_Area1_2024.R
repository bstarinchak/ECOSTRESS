


# This script create some simple charts of LAI dataset. 
# ACS 2024


library(tidyverse)
library(lubridate)
library(anytime)



setwd("C:/Users/antonio.silva/DadosTanguro Dropbox/Antônio Carlos Silveiro da Silva/projetos/Savanizacao/LAI/Area1/LAI_area1_ABC")

lai = read.csv("MASTER_LAI_Area1_ABC_may2024.csv")


##

lai3 = lai %>% 
  #filter(ano %in%c(2023, 2024)) %>% 
  mutate(date = paste(ano, mes, 15, sep="_"),
         date2 = anydate(date),
         treat = ifelse(linhas %in%c(1:10), "Control",
                        ifelse(linhas %in%c(11:20), "B3yr", "B1yr"))) 

lai3 = lai3 %>% filter(!transecto == "BORDA") # remove
lai3 = lai3 %>% filter(!transecto == "") 

unique(lai3$linhas)
unique(lai3$transecto)
names(lai3)

plot(lai3$LAI)


lai3.AA = lai3 %>% filter(transecto == "D")


lai3$transecto2 = factor(lai3$transecto, levels = c("Edge", "A", "AA", "AB", "B", "C", "D", "E", "F","G", "H","I", "J", "K", "L", "M","N",
                                                    "O","P","Q", "R", "S","T", "U" ))



lai3 %>% 
  filter(LAI > 0) %>% 
  #group_by()
  ggplot(aes(linhas, LAI, colour = transecto2, fill = transecto2))+
  geom_smooth()+
  geom_point()+
  facet_wrap(~transecto2)+
  theme_light(base_size = 20)


lai3 %>% 
  filter(LAI > 0) %>% 
  filter(transecto2%in% c("Edge", "A",  "AB", "B", "C", "F",   "K", "P","U" )) %>% 
  #filter(transecto2 == "Edge") %>% 
  #group_by()
  ggplot(aes(date2, LAI, colour = transecto2, fill = transecto2))+
  geom_smooth()+
  geom_point()+
  labs(x = "Date", y = "LAI (m².m²)")+
  #scale_colour_manual(values = c('#A50026', '#D73027', '#F46D43', 'darkorange', '#FDAE61',   '#A6D96A', '#66BD63', '#1A9850', '#006837'))+
  #scale_fill_manual(values = c('#A50026', '#D73027', '#F46D43', 'darkorange', '#FDAE61',   '#A6D96A', '#66BD63', '#1A9850', '#006837'))+
  facet_wrap(~transecto2, scales = "free")+
  theme_light(base_size = 20)


lai3 %>% 
  filter(LAI > 0) %>% 
  filter(transecto2%in% c("Edge", "A",  "AB", "B", "C", "F",   "K", "P","U" )) %>% 
  
  ggplot(aes(ano, LAI, colour = transecto2, fill = transecto2))+
  #geom_line()+
  stat_summary(fun.data = "mean_se")
  geom_point()+
  labs(x = "Date", y = "LAI (m².m²)")+
  #scale_colour_manual(values = c('#A50026', '#D73027', '#F46D43', 'darkorange', '#FDAE61',   '#A6D96A', '#66BD63', '#1A9850', '#006837'))+
  #scale_fill_manual(values = c('#A50026', '#D73027', '#F46D43', 'darkorange', '#FDAE61',   '#A6D96A', '#66BD63', '#1A9850', '#006837'))+
  facet_wrap(~transecto2, scales = "free")+
  theme_light(base_size = 20)



lai3 %>% 
  filter(LAI > 0) %>% 
  filter(transecto2%in% c("Edge", "A",  "AB", "B", "C", "F",   "K", "P","U" )) %>% 
  #filter(transecto2 == "Edge") %>% 
  group_by(date2, transecto2) %>% 
  summarise(lai.mean = mean(LAI, na.rm = T),
            lai.inf = quantile(LAI, na.rm = T)[2],
            lai.sup = quantile(LAI, na.rm = T)[4]) %>% 
  ggplot(aes(date2, lai.mean, colour = transecto2, fill = transecto2))+
  #geom_smooth()+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin = lai.inf,
                    ymax = lai.sup), width = 0.2)+
  #geom_ribbon(aes(ymin = lai.inf, ymax = lai.sup), alpha = 0.5) +
  labs(x = "Date", y = "LAI (m².m²)")+
  #scale_colour_manual(values = c('#A50026', '#D73027', '#F46D43', 'darkorange', '#FDAE61',   '#A6D96A', '#66BD63', '#1A9850', '#006837'))+
  #scale_fill_manual(values = c('#A50026', '#D73027', '#F46D43', 'darkorange', '#FDAE61',   '#A6D96A', '#66BD63', '#1A9850', '#006837'))+
  facet_wrap(~transecto2, scales = "free")+
  theme_light(base_size = 20)



### Borda e Interior


lai3 %>% 
  filter(LAI > 0) %>% 
  filter(transecto2%in% c("Edge", "A",  "AB", "B", "C", "F",   "K", "P","U" )) %>% 
  mutate(edge.interior = ifelse(transecto2%in% c("Edge", "A",  "AB", "B", "C"), "Edge", "Interior")) %>% 
  #filter(transecto2 == "Edge") %>% 
  group_by(ano, edge.interior, treat) %>% 
  summarise(lai.mean = mean(LAI, na.rm = T),
            lai.inf = quantile(LAI, na.rm = T)[2],
            lai.sup = quantile(LAI, na.rm = T)[4]) %>% 
  ggplot(aes(ano, lai.mean, colour = edge.interior, fill = edge.interior))+
  geom_smooth()+
  geom_point()+
  #geom_line()+
  #geom_errorbar(aes(ymin = lai.inf, ymax = lai.sup), width = 0.2)+
  
  labs(x = "Date", y = "LAI (m².m²)")+
  facet_wrap(~treat, scales = "free")+
  theme_light(base_size = 20)



lai3 %>% 
  filter(LAI > 0) %>% 
  #filter(!ano == 2024) %>% 
  filter(transecto2%in% c("Edge", "A",  "AB", "B", "C", "F",   "K", "P","U" )) %>% 
  mutate(edge.interior = ifelse(transecto2%in% c("Edge", "A",  "AB", "B", "C"), "Edge", "Interior")) %>% 
  ggplot(aes(ano, LAI, colour = treat, group = treat))+
  facet_wrap(~edge.interior, scales = "free")+
  stat_summary(fun.data = "mean_se", geom = "line")+
  stat_summary(fun.data = "mean_se", size = 1)+
  labs(x = "Date", y = "LAI (m².m²)")+
  scale_colour_manual(values = c("darkred", "darkorange", "darkgreen"))+
  theme_minimal(base_size = 22)+
  theme(legend.background = element_blank(),
        legend.position = c(0.3,0.9),
        legend.title = element_blank(),
        axis.text = element_text( colour = "black"))


###
###

# Recent data
lai2 = lai %>% 
  filter(ano %in%c(2023, 2024)) %>% 
  mutate(date = paste(ano, mes, 15, sep="_"),
         date2 = anydate(date),
         treat = ifelse(linhas %in%c(1:10), "Control",
                        ifelse(linhas %in%c(11:20), "B3yr", "B1yr")))



lai2$transecto2 = factor(lai2$transecto, levels = c("Edge", "A",  "AB", "B", "C", "F",   "K", "P","U" ))


lai2 %>% 
  filter(LAI > 0) %>% 
  #group_by()
  ggplot(aes(linhas, LAI, colour = transecto2, fill = transecto2))+
  geom_smooth()+
  geom_point()+
  labs(x = "Linhas", y = "LAI (m².m²)")+
  #scale_x_continuous(breaks = c(1:30))+
  geom_vline(xintercept = c(11,21), linetype = 2)+
  scale_colour_manual(values = c('#A50026', '#D73027', '#F46D43', 'darkorange', '#FDAE61',   '#A6D96A', '#66BD63', '#1A9850', '#006837'))+
  scale_fill_manual(values = c('#A50026', '#D73027', '#F46D43', 'darkorange', '#FDAE61',   '#A6D96A', '#66BD63', '#1A9850', '#006837'))+
  facet_wrap(~transecto2)+
  theme_light(base_size = 20)



lai2 %>% 
  filter(LAI > 0) %>% 
  #filter(transecto2 == "Edge") %>% 
  #group_by()
  ggplot(aes(date2, LAI, colour = transecto2, fill = transecto2))+
  geom_smooth()+
  geom_point()+
  labs(x = "Date", y = "LAI (m².m²)")+
  scale_colour_manual(values = c('#A50026', '#D73027', '#F46D43', 'darkorange', '#FDAE61',   '#A6D96A', '#66BD63', '#1A9850', '#006837'))+
  scale_fill_manual(values = c('#A50026', '#D73027', '#F46D43', 'darkorange', '#FDAE61',   '#A6D96A', '#66BD63', '#1A9850', '#006837'))+
  facet_wrap(~transecto2)+
  theme_light(base_size = 20)




lai2 = lai2 %>% 
  mutate(dd = anydate(paste(ano, mes, sep = "_")))

lai2 %>% 
  filter(LAI > 0) %>% 
   
  dplyr::group_by( transecto2, dd) %>% 
  dplyr::summarise(lai.mean = mean(LAI, na.rm=T),
                   lai.lower = quantile(LAI, na.rm = T)[2],
                   lai.upper = quantile(LAI, na.rm = T)[4]) %>% 
  ggplot(aes(dd, lai.mean, colour = transecto2, group = transecto2))+
  geom_point(size = 4)+
  geom_line()+
  labs(x = "Date", y = "LAI (m².m²)")+
  facet_wrap(~transecto2)+
  scale_colour_manual(values = c('#A50026', '#D73027', '#F46D43', 'darkorange', '#FDAE61',   '#A6D96A', '#66BD63', '#1A9850', '#006837'))+
  scale_fill_manual(values = c('#A50026', '#D73027', '#F46D43', 'darkorange', '#FDAE61',   '#A6D96A', '#66BD63', '#1A9850', '#006837'))+
  geom_ribbon(aes(ymin = lai.lower, ymax =  lai.upper , fill =transecto2), alpha = 0.3) +
  theme_light(base_size = 20)
  

  
  lai2 %>% 
    filter(LAI > 0) %>% 
    #filter(transecto2 == "Edge") %>% 
    #group_by()
    ggplot(aes( LAI, colour = transecto2))+
    geom_density()+
    labs( x = "LAI (m².m²)")+
    scale_colour_manual(values = c('#A50026', '#D73027', '#F46D43', 'darkorange', '#FDAE61',   '#A6D96A', '#66BD63', '#1A9850', '#006837'))+
    #facet_wrap(~transecto2)+
    theme_light(base_size = 20)
  
  

rm(list = ls())
