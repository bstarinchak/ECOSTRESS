#library(dplyr)
library(tidyverse)
#library(ggplot2)
library(dplyrUtil, quietly = TRUE)
#library(sapproc)
library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
#library(tidyverse)
library(purrr, quietly = TRUE)
#library(dplyrUtil, quietly = TRUE)
library(ggplot2, quietly = TRUE)

###analyse of the 14C in the sugars and starch form the Tanguro data

setwd("/Users/_dherrera/Documents/Yale_2024/Manuscripts/Lateral_mobilization _NSC_stemwood/analyzes/")

Starch_may_2019=read.csv("data/14C_starch_may_2019_nostandards.csv", header=T)
Starch_may_2019=Starch_may_2019%>%separate(Probe, c("ID_tree", "wood_core", "depth"), sep="-")
# starch_02=read.csv2("Starch_14C_02_jul_18.csv", header = T)
# sugars_02=read.csv2("sugars_14C_02_jul18.csv")
# sugars_24=read.csv2("sugars_14C_24_jul18.csv")

names(Starch_may_2019)
Starch_may_2019$wood_core_mean_age=NA

# sugars_02[!(sugars_02$ID_tree%in%sugars_24$ID_tree),]
# sugars_24[!(sugars_24$ID_tree%in%sugars_02$ID_tree),]
# setdiff(sugars_02$ID_tree, sugars_24$ID_tree)
# 
# sugars_14C=merge(sugars_02,sugars_24, by=names(sugars_24), all=T)
# 
# starch_14C=Starch_02_merge[,c("P.Nr", "species", "ID_tree", "depth", "sugar", 'F14C', "err", "correction_14C_starch",
#                               "errorD14C")]
# #colnames(starch_14C)=names(sugars_24)
# 
# #sugars_14C_m=rbind(sugars_14C, starch_14C)
# names(sugars_14C)
# names(starch_14C)
# sugars_14C_merge=sugars_14C[,c(1,2,3,4,5,6,7,8, 9)]
# names(starch_14C)[c(1,8,9)]=c("P.Nr.","D14C", "err.1")
# names(sugars_14C_merge)
# names(starch_14C)
# sugars_starch_14C=rbind(sugars_14C_merge, starch_14C)
# sugars_starch_14C=filter(sugars_starch_14C, D14C>-100)
# Starch_may_2019$species[which(Starch_may_2019$species=="dac")]="Dac"
# sugars_starch_14C$species[which(sugars_starch_14C$species=="L.Prat")]="L.prat"
# sugars_starch_14C$species[which(sugars_starch_14C$species=="L.Part")]="L.prat"
# sugars_starch_14C$species[which(sugars_starch_14C$species=="sac")]="Sac"
# sugars_starch_14C$species[which(sugars_starch_14C$species=="sbnha")]="Sbnha"
# sugars_starch_14C$species[which(sugars_starch_14C$species=="sbra")]="Sbra"
# sugars_starch_14C$species[which(sugars_starch_14C$species=="sch")]="Sch"
# sugars_starch_14C$species[which(sugars_starch_14C$species=="Tab")]="Tap"


ggplot(Starch_may_2019, aes(x=Species, y=D14C))+geom_point()+
  facet_wrap(~depth)
# 
# sugars_14C$species[which(sugars_14C$species=="dac")]="Dac"
# sugars_14C$species[which(sugars_14C$species=="L.Prat")]="L.prat"
# sugars_14C$species[which(sugars_14C$species=="L.Part")]="L.prat"
# sugars_14C$species[which(sugars_14C$species=="sac")]="Sac"
# sugars_14C$species[which(sugars_14C$species=="sbnha")]="Sbnha"
# sugars_14C$species[which(sugars_14C$species=="sbra")]="Sbra"
# sugars_14C$species[which(sugars_14C$species=="sch")]="Sch"
# sugars_14C$species[which(sugars_14C$species=="Tab")]="Tap"



ggplot(Starch_may_2019, aes(x=depth, y=D14C, col=ID_tree, label=ID_tree))+geom_point()+geom_text()+
  facet_wrap(~Species)
  

#sugars_14C[sugars_14C$D14C<0 & !is.na(sugars_14C$D14C),c(6:9)]=NA


str(Starch_may_2019)
Starch_may_2019$depth2=as.numeric(as.factor(Starch_may_2019$depth))
Starch_may_2019$ID_sp=Starch_may_2019$species

sp_names=c("S. guianensis",
           "O. guianansis",
           "D. microcarpa"
            )
Starch_may_2019$ID_sp=Starch_may_2019$Species

library(plyr)
Starch_may_2019$ID_sp=mapvalues(Starch_may_2019$ID_sp,
                           unique(Starch_may_2019$Species), sp_names)
detach("package:plyr", unload=TRUE)

Starch_may_2019$ID_sp=factor(Starch_may_2019$ID_sp, 
                        levels = c("D. microcarpa",
                                   "S. guianensis",
                                   "O. guianansis"
                                   ))

#### standaridizoing by the age of the wood

path2="/Users/_dherrera/SOIL-R/Theses/David_Herrera/sugar extraction analysis /second_field_campaign/cellulose_14C/"
cellulose_14C_tanguro_trees=read.csv2(file.path(path2, "14C_calibrated.csv"), header = T, sep = ";")
cellulose_14C_tanguro_trees=cellulose_14C_tanguro_trees%>%separate(Probe,c("Species", "ID_tree","Sample_type", "ring_count", "sampling_date"), sep="_")
cellulose_14C_tanguro_trees=cellulose_14C_tanguro_trees%>%
  mutate(year_estimated=2018-cellulose_14C_tanguro_trees$ring_number)

cellulose_14C_tanguro_trees=cellulose_14C_tanguro_trees%>%separate(ID_tree,c("ID_tree", "core_No"))
unique(cellulose_14C_tanguro_trees$Species)

cellulose_sac=cellulose_14C_tanguro_trees%>%filter(Species %in% c("Sac", "ab"))
cellulose_sac$depth=cellulose_sac$ring_count
cellulose_sac=cellulose_sac[-c(1,2),]

str(cellulose_sac)
str(Starch_may_2019) 

unique(Starch_may_2019$depth)

Starch_may_2019$depth=ifelse(Starch_may_2019$depth=="2","0-2", "2-4")
Starch_may_2019$depth=factor(Starch_may_2019$depth, levels = c("0-2", "2-4"))
     
unique(Starch_may_2019$depth)
unique(cellulose_sac$ring_count)

unique(Starch_may_2019$ID_tree)
unique(cellulose_sac$ID_tree)

Starch_may_2019=Starch_may_2019[order(Starch_may_2019$depth2),]
rownames(Starch_may_2019)=NULL

Starch_may_2019$wood_core_tree_rings_start=NA
for(i in 1:nrow(Starch_may_2019)){
  if(Starch_may_2019$depth[i]=="0-2"){
    Starch_may_2019$wood_core_tree_rings_start[i]=0
  }
  else{
    for(j in 1:nrow(Starch_may_2019[Starch_may_2019$depth=="0-2",])){
      if(Starch_may_2019$ID_tree[i]==Starch_may_2019$ID_tree[j]){
        Starch_may_2019$wood_core_tree_rings_start[i]=Starch_may_2019$wood_core_tree_rings[j]
      }
    }
  }
}

Starch_may_2019$wood_core_mean_age=with(Starch_may_2019, (wood_core_tree_rings+wood_core_tree_rings_start)/2)



for(i in 1:nrow(Starch_may_2019)){
  if(is.na(Starch_may_2019$wood_core_mean_age[i])){
    for(j in 1:nrow(cellulose_sac)){
      if(Starch_may_2019$ID_tree[i]==cellulose_sac$ID_tree[j] &&
         Starch_may_2019$depth[i]== cellulose_sac$ring_count[j]){
        Starch_may_2019$wood_core_mean_age[i]=cellulose_sac$Years_since_collected[j]
      }
    }
  }
}



library("SoilR")

calib_curv=filter(Hua2021$SHZone3, Year>1970)
calib_curv2=Hua2013$NHZone3

C14_Trpcs <- read.table("C14Graven_Tropics_IntCal20")
C14_trpcs_graven = read.table("C14Graven_Tropics")
C14_NH <- read.table("C14Graven_PNAS_GMD_Int20_2100_-53050")


max(calib_curv[,4])
min(calib_curv[,4])

plot(calib_curv$Year.AD, calib_curv$mean.F14C)

#Sugars_cal_age=Starch_may_2019[,c(3,6,7, 8,9)]

calibrating_F14C_to_calendar_dates=function(df, calcurve){##df must me a data frame with the F14 and errors data in teh second and third column respectively
  df$calibrated_dates=NA
  for(i in 1:nrow(df)){
    possible_Dates=filter(calcurve, mean.F14C>df[i,"F14C"]-df[i,"err"] & mean.F14C<df[i,"F14C"]+df[i,"err"])
    if(nrow(possible_Dates)==0){
      df$calibrated_dates[i]=2253.3-230.8*df$F14C[i]
    }else{
      df$calibrated_dates[i]=mean(possible_Dates[,1])
    }
  }
  return(df)
  
}

calibrating_D14C_to_calendar_dates=function(df, calcurve){##df must me a data frame with the D14 and errors data in teh second and third column respectively
  df$calibrated_dates=NA
  for(i in 1:nrow(df)){
    possible_Dates=filter(calcurve, mean.Delta14C>df[i,"D14C"]-df[i,"err.1"] & mean.Delta14C<df[i,"D14C"]+df[i,"err.1"])
    if(nrow(possible_Dates)==0){
      df$calibrated_dates[i]=2020.52-0.2257*df$D14C[i]
    }else{
      df$calibrated_dates[i]=mean(possible_Dates[,1])
    }
  }
  return(df)
  
}

#df=sugars_14C[,c(3,6,7, 8,9)]
# i=1
# filter(C14_NH_filtered, mean.Delta14C>df[i,4]-df[i,5]*4 & mean.Delta14C<df[i,4]+df[i,5]*4)
# filter(calib_curv, mean.F14C>df[i,2]-df[i,3] & mean.F14C<df[i,2]+df[i,3])
# 2020.52-0.2257*df$D14C[i]

Sugars_cal_age_calib_curve=calibrating_F14C_to_calendar_dates(Starch_may_2019, calib_curv)

Sugars_cal_age_calib_curve_D14C=calibrating_D14C_to_calendar_dates(Starch_may_2019, calib_curv)
# colnames(C14_NH)=c("Year.AD", "mean.Delta14C")
# C14_NH_filtered=C14_NH%>%filter(Year.AD>1970 & Year.AD<2050)
# Sugars_cal_age_C14_NH=calibrating_D14C_to_calendar_dates(Sugars_cal_age, C14_NH_filtered)
# Sugars_cal_age_calib_curve$C14_NH_Cal=Sugars_cal_age_C14_NH$calibrated_dates


# ggplot(calib_curv)+geom_point(aes(x=Year.AD, y=mean.Delta14C))+
#   geom_point(data=C14_NH_filtered, aes(x=Year.AD, y=mean.Delta14C), col="red")+
#   geom_point(data=calib_curv, aes(x=Year.AD, y=mean.Delta14C), col="blue")+
#   geom_line(data=calib_curv, aes(x=Year.AD, y=lm_modeled), col="green")

#ggplot(C14_NH_filtered, aes(x=Year.AD, y=mean.Delta14C))+geom_point()


#sugars_14C$sugar_cal_age_calib_curv=Sugars_cal_age_calib_curve$calibrated_dates
Sugars_cal_age_calib_curve$sugars_mean_age_calib_curv=2019-Sugars_cal_age_calib_curve$calibrated_dates
Sugars_cal_age_calib_curve$Calibrated_D14C=Sugars_cal_age_calib_curve_D14C$calibrated_dates
Sugars_cal_age_calib_curve$sugars_mean_age_calib_curv_D14C=2019-Sugars_cal_age_calib_curve$Calibrated_D14C


png("sugar_mean_age_per_wood_depth_may2019.png", width = 800, height = 600)
ggplot(Sugars_cal_age_calib_curve, aes(x=depth, y=sugars_mean_age_calib_curv_D14C))+
  geom_point()+
  geom_text(aes(label=ID_tree))+
  geom_line(aes(x=depth2, y=sugars_mean_age_calib_curv_D14C, col=ID_tree))+
  facet_wrap(~ID_sp)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(2)),
        axis.text.y=element_text(size = rel(2)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        plot.title = element_text(size = rel(1.5), face = "italic"))
dev.off()


# sugars_14C$std_sugar_age=sugars_14C$sugars_mean_age_calib_curv/sugars_14C$wood_core_mean_age
# 
# ggplot(sugars_14C, aes(x=depth, y=std_sugar_age))+
#   geom_point()+
#   geom_text(aes(label=ID_tree))+
#   geom_line(aes(x=depth2, y=std_sugar_age, col=ID_tree))+
#   facet_wrap(~ID_sp)+
#   theme_bw()+
#   theme(axis.text.x=element_text(size = rel(2)),
#         axis.text.y=element_text(size = rel(2)),
#         axis.title.x = element_text(size = rel(1.5)),
#         axis.title.y = element_text(size = rel(1.5)),
#         strip.text = element_text(size = rel(1.5)),
#         legend.position="none",
#         plot.title = element_text(size = rel(1.5), face = "italic"))

ggplot(Sugars_cal_age_calib_curve, aes(x=wood_core_mean_age, y=sugars_mean_age_calib_curv_D14C, col=ID_sp))+ geom_point()+
  geom_smooth(method="lm", se=F)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(2)),
        axis.text.y=element_text(size = rel(2)),
        axis.title.x = element_text(size = rel(1.5)),
        axis.title.y = element_text(size = rel(1.5)),
        strip.text = element_text(size = rel(1.5)),
        #legend.position="none",
        plot.title = element_text(size = rel(1.5), face = "italic"))


#png("/Users/_dherrera/SOIL-R/Theses/David_Herrera/NSC_dynamics_TT_age_and_seasonality/calculations/means_Age_wood_vs_mean_age_sugars_fiber.png", width = 800, height = 600)
# ggplot(filter(sugars14C_lm_coef_all, ID_sp!="S. morototoni"), aes(x=wood_core_tree_rings, y=sugars_mean_age_calib_curv, col=ID_sp_slope))+ geom_point()+
#   geom_smooth(method="lm", se=F)+
#   scale_color_discrete("Species = slope")+
#   geom_abline(intercept = 0, slope = 1, linetype="dashed", color="black")+
#   xlab("Wood core mean age (yr)")+ylab("sugars mean age (yr)")+
#   theme_bw()+
#   theme(axis.text.x=element_text(size = rel(3)),
#         axis.text.y=element_text(size = rel(3)),
#         axis.title.x = element_text(size = rel(2.5)),
#         axis.title.y = element_text(size = rel(2.5)),
#         strip.text = element_text(size = rel(2.5)),
#         legend.text = element_text(size=rel(2)),
#         legend.title = element_text(size = rel(3)),
#         legend.key = element_rect(size = rel(2)),
#         #legend.position="none",
#         plot.title = element_text(size = rel(2.5), face = "italic"))
# #dev.off()


for (i in 1:nrow(Sugars_cal_age_calib_curve)){
  if(is.na(Sugars_cal_age_calib_curve$wood_core_tree_rings[i])){
    if(Sugars_cal_age_calib_curve$depth[i]=="0-2"){
      Sugars_cal_age_calib_curve$wood_core_tree_rings[i]=2*Sugars_cal_age_calib_curve$wood_core_mean_age[i]
    }
    else{
      for(j in 1:nrow(Sugars_cal_age_calib_curve[Sugars_cal_age_calib_curve$depth=="0-2",])){
        if(Sugars_cal_age_calib_curve$ID_tree[i]==Sugars_cal_age_calib_curve$ID_tree[j]){
          Sugars_cal_age_calib_curve$wood_core_tree_rings[i]=2*Sugars_cal_age_calib_curve$wood_core_mean_age[i]-Sugars_cal_age_calib_curve$wood_core_tree_rings[j]
        }
      }
    }
  }
}

#Sugars_cal_age_calib_curve$std_by_treering_sugar_age=Sugars_cal_age_calib_curve$sugars_mean_age_calib_curv/Sugars_cal_age_calib_curve$wood_core_tree_rings

Sugars_cal_age_calib_curve=Sugars_cal_age_calib_curve%>%
  mutate(storage_strategy=ifelse(Species%in%c("Dac"), "fiber_storage", "parenchyma_storage"),
         leaf_habit=ifelse(Species%in%c("Ab"), "evergreen", "semi-deciduous"))

library(broom)
Starch_may_14C_lm=Sugars_cal_age_calib_curve%>%group_by(ID_sp)%>%nest()%>%
  mutate(lm=map(data, ~lm(sugars_mean_age_calib_curv_D14C~wood_core_mean_age, data=.x)),
         tidy_lm=map(lm, tidy))

Starch_may_14C_lm_coef=Starch_may_14C_lm%>%unnest(tidy_lm)
Starch_may_14C_lm_coef=Starch_may_14C_lm_coef%>%filter(term=="wood_core_mean_age")

Starch_may_14C_lm_coef_all=merge(Sugars_cal_age_calib_curve, Starch_may_14C_lm_coef[,c(1,5, 6, 7, 8)], by="ID_sp")


Starch_may_14C_lm_coef_all$ID_sp_slope=paste(Starch_may_14C_lm_coef_all$ID_sp,  round(Starch_may_14C_lm_coef_all$estimate, 2), sep = "= ")
Starch_may_14C_lm_coef_all$ID_sp_slope=paste(Starch_may_14C_lm_coef_all$ID_sp_slope,"(",  round(Starch_may_14C_lm_coef_all$p.value, 2), ")", sep = "")
label_storage_strategy=c("Fiber storage", "Parenchyma storage")
names(label_storage_strategy)=unique(Starch_may_14C_lm_coef_all$storage_strategy)

png("means_Age_wood_vs_mean_age_sugars_may19_storage_streatgy.png", width = 1000, height = 800)
ggplot(Starch_may_14C_lm_coef_all, aes(x=wood_core_mean_age, y=sugars_mean_age_calib_curv_D14C, col=ID_sp_slope))+ geom_point(size=3)+
  geom_smooth(method="lm", se=F, size=2)+
  scale_color_discrete("Species = slope (p.value)")+
  geom_abline(intercept = 0, slope = 1, linetype="dashed", color="black")+
  xlab("Wood core mean tree rings age (yr)")+ylab(expression(paste("sugars", " C"^14, " age (yr)")))+
  theme_bw()+
  facet_wrap(~storage_strategy, labeller = labeller(storage_strategy=label_storage_strategy))+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(2.5)),
        axis.title.y = element_text(size = rel(2.5)),
        strip.text = element_text(size = rel(2.5)),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size=rel(2)),
        legend.title = element_text(size = rel(3)),
        legend.key = element_rect(size = rel(2)),
        #legend.position="none",
        plot.title = element_text(size = rel(2.5), face = "italic"))
dev.off()

png("means_Age_wood_vs_mean_age_sugars_leaf_habit_may2019.png", width = 1000, height = 800)
ggplot(Starch_may_14C_lm_coef_all, aes(x=wood_core_mean_age, y=sugars_mean_age_calib_curv_D14C, col=ID_sp_slope))+ geom_point(size=3)+
  geom_smooth(method="lm", se=F, size=2)+
  scale_color_discrete("Species = slope(p.value)")+
  geom_abline(intercept = 0, slope = 1, linetype="dashed", color="black")+
  xlab("Wood core mean tree rings age (yr)")+ylab(expression(paste("Sugars", ""^14,"C", " age (yr)")))+
  theme_bw()+
  facet_wrap(leaf_habit~storage_strategy, labeller = labeller(storage_strategy=label_storage_strategy))+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(2.5)),
        axis.title.y = element_text(size = rel(2.5)),
        strip.text = element_text(size = rel(2.5)),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size=rel(1.5)),
        legend.title = element_text(size = rel(2.5)),
        legend.key = element_rect(size = rel(1.5)),
        #legend.position="none",
        plot.title = element_text(size = rel(2.5), face = "italic"))
dev.off()
