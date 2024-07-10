library(tidyverse)
library(dplyr)
library(dplyrUtil, quietly = TRUE)
library(ggplot2)
library(lme4)
library(nlme)
library(broom)
library(ggpubr)
library(rstatix)
library(reshape2)
###reading the data of parenchyma and starch out of the ImageJ files
source("/Users/_dherrera/SOIL-R/Theses//David_Herrera/R_function_stat_smooth.R")
source("/Users/_dherrera/SOIL-R/Theses/David_Herrera/NSC_dynamics_TT_age_and_seasonality/calculations/functions.R")
setwd("/Users/_dherrera/Documents/balzan_project/Manuscript_submissions/NSC_seasonality_and_growth/Journal_of_ecology/revision1/data_processing/")
load("growth_data_tanguro_2018_2020.RData")

###### ab


load("ab_profile_data")

names(ab_profile_data_Jan18)
str(ab_profile_data_Jan18)
names(ab_profile_data_Feb20)

### selecting the varables that are gonna be used in this analysis

ab_May19=ab_profile_data_May19%>%select("ID", "depth_mm", "measurement_type", "X.Area", "sampling_month", "month")
ab_Aug19=ab_profile_data_Aug19%>%select("ID", "depth_mm", "measurement_type", "X.Area", "sampling_month", "month")
ab_Nov19=ab_profile_data_Nov19%>%select("ID", "depth_mm", "measurement_type", "X.Area", "sampling_month", "month")
ab_feb19=ab_profile_data_Feb20%>%select("ID", "depth_mm", "measurement_type", "X.Area", "sampling_month", "month")

### unify the data frame for the species and choosing the measurement type that will be analyzed (woodstarch in this case)

Ab=rbind(ab_May19, ab_Aug19, ab_Nov19, ab_feb19)
Ab$ID[which(Ab$ID=="ab34283 1")]="ab34288 1"

Ab=Ab%>%separate(ID, c("ID", "radii_No"), sep=" ")


Ab_woodStarch=filter(Ab, measurement_type=="woodStarch")
Ab_woodStarch$depth_mm=as.numeric(Ab_woodStarch$depth_mm)
str(Ab_woodStarch)

### generate the cero values for the inner layesr of wood that were not measured because the absence of starch 

Ab_woodStarch=Ab_woodStarch %>% group_by(ID, month)%>% 
  mapGroups(function(groups_ID){
    max_depth_spp=max(Ab_woodStarch$depth_mm)
    max_depth_grou=max(groups_ID$depth_mm)
    while(max_depth_grou<max_depth_spp){
      i=length(groups_ID$depth_mm)
      groups_ID[c(i+1, i+2, i+3), "depth_mm"]=max(groups_ID$depth_mm)+5
      groups_ID$ID[c(i+1, i+2, i+3)]=groups_ID$ID[1]
      groups_ID$measurement_type="woodStarch"
      groups_ID$month=groups_ID$month[1]
      groups_ID$sampling_month[c(i+1, i+2, i+3)]=groups_ID$sampling_month[1]
      max_depth_grou=max(groups_ID$depth_mm)
    }
    groups_ID$X.Area[is.na(groups_ID$X.Area)]=0
    return(groups_ID)
  }
  )%>%ungroup() 


ab_check_filling=Ab_woodStarch%>%filter(depth_mm=="60")

save(Ab_woodStarch, file = "Ab_woodStarch.Rda")
save(Ab_woodStarch, file = "Ab_woodStarch.csv")

############## Calculating mean statistics per individual, month and sampling depth ####################

## mean statistics by tree 

Ab_woodStarch_means_by_tree=Ab_woodStarch%>%group_by(ID, depth_mm, month, sampling_month)%>%
  summarize(mean_starch_area=mean(X.Area, na.rm=T),
            median_starch_are=median(X.Area, na.rm=T),
            Int_quart_range=IQR(X.Area, na.rm=T),
            sd_starch_area=sd(X.Area, na.rm=T),
            sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area)))


### mean statistics by species 

Ab_woodStarch_means_by_specie=Ab_woodStarch%>%group_by(depth_mm, month, sampling_month)%>%
  summarize(mean_starch_area=mean(X.Area, na.rm=T),
            median_starch_are=median(X.Area, na.rm=T),
            Int_quart_range=IQR(X.Area, na.rm=T),
            sd_starch_area=sd(X.Area, na.rm=T),
            sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area)))

Ab_woodStarch_means_by_tree$month=factor(Ab_woodStarch_means_by_tree$month, 
                                         levels = c("May19", "Aug19", "Nov19", "Feb20"))
Ab_woodStarch_means_by_specie$month=factor(Ab_woodStarch_means_by_specie$month, 
                                           levels = c("May19", "Aug19", "Nov19", "Feb20"))

Ab_woodStarch_means_by_specie$month_season=Ab_woodStarch_means_by_specie$month

library(plyr)

months=unique(Ab_woodStarch_means_by_specie$month)
months_season=c("May 19 (transition)", "Aug 19 (Dry season)", "Nov 19 (transition)", "Feb 20 (wet season)")

Ab_woodStarch_means_by_specie$month_season=mapvalues(Ab_woodStarch_means_by_specie$month_season,
                                                     months,
                                                     months_season)
detach("package:plyr", unload=TRUE)

Ab_woodStarch_means_by_specie$month_season=factor(Ab_woodStarch_means_by_specie$month_season,
                                                  labels = c("May 19 (transition)", 
                                                             "Aug 19 (Dry season)", 
                                                             "Nov 19 (transition)", 
                                                             "Feb 20 (wet season)"))

png("radial_starch_content_O-leu_scale_seasons_withlengend.png", width = 800, height = 600)
ggplot()+
  geom_line(data=Ab_woodStarch_means_by_specie , aes(x=depth_mm, y=mean_starch_area, col=month_season), size=1)+
  geom_ribbon(data=Ab_woodStarch_means_by_specie, aes(x=depth_mm, 
                                                      ymin=mean_starch_area-sd_starch_area,
                                                      ymax=mean_starch_area+sd_starch_area,
                                                      fill= month_season),
              #fill="gray80",
              alpha=0.1)+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  xlab("depth (mm)")+ylab("starch %")+
  labs(color="Seasons", fill="Seasons")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(2.5)),
        axis.title.y = element_text(size = rel(2.5)),
        strip.text = element_text(size = rel(3)),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(3)),
        legend.title = element_text(size = rel(3)),
        #legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

png("radial_starch_content_O-leu_scale_seasons.png", width = 800, height = 600)
ggplot()+
  geom_line(data=Ab_woodStarch_means_by_specie , aes(x=depth_mm, y=mean_starch_area, col=month_season), size=1)+
  geom_ribbon(data=Ab_woodStarch_means_by_specie, aes(x=depth_mm, 
                                                      ymin=mean_starch_area-sd_starch_area,
                                                      ymax=mean_starch_area+sd_starch_area,
                                                      fill= month_season),
              #fill="gray80",
              alpha=0.1)+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  xlab("depth (mm)")+ylab("starch %")+
  labs(color="Seasons", fill="Seasons")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(3)),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = rel(3)),
        legend.title = element_text(size = rel(3)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


## statistical comparisons at different depth ranges (ANOVA)

Ab_pairwise_comparison=Ab_woodStarch_means_by_tree
Ab_pairwise_comparison=Ab_pairwise_comparison[-c(1:24), c(1,2,3,5)]
Ab_pairwise_comparison$time=as.factor(Ab_pairwise_comparison$month)
Ab_pairwise_comparison=Ab_pairwise_comparison[order(Ab_pairwise_comparison$month),]
rownames(Ab_pairwise_comparison)=NULL

Ab_pairwise_comparison$ID=as.factor(Ab_pairwise_comparison$ID)
AB_Anva2=Ab_pairwise_comparison%>%group_by(depth_mm)%>%nest()%>%
  mutate(ANOVA=map(data, ~ aov(mean_starch_area~month, data=.x)),
         tukey_comparison=map(ANOVA, TukeyHSD),
         tidied=map(ANOVA, tidy),
         tidied2=map(tukey_comparison, tidy), 
  )


names(AB_Anva2$tidied)=AB_Anva2$depth_mm
names(AB_Anva2$tidied2)=AB_Anva2$depth_mm
AB_Anva2$tidied
AB_Anva2$tidied2
estimates_ANOVA_by_depth=AB_Anva2%>%unnest(tidied)

Ab_pairwise_comparison2=filter(Ab_pairwise_comparison, ID!="ab34288")
Ab_pairwise_comparison2$ID=factor(Ab_pairwise_comparison2$ID)
Ab_pairwise_comparison2$depth_mm2=as.factor(Ab_pairwise_comparison2$depth_mm)
table(Ab_pairwise_comparison2$ID)

str(Ab_pairwise_comparison2)
AB_Anva3=Ab_pairwise_comparison2%>%group_by(depth_mm2)%>%
  wilcox_test(mean_starch_area~month, data=., paired = T)
  

ggplot(Ab_woodStarch_means_by_tree, aes(x=month, y=mean_starch_area))+
  geom_boxplot(data=Ab_woodStarch_means_by_tree, aes(x=month, y=mean_starch_area))+
  geom_point(aes(col=ID))+
  geom_point(data=Ab_woodStarch_means_by_specie, aes(x=month, y=mean_starch_area), col="black")+
  facet_wrap(~depth_mm)


### estimate the maximum and relative seasonal range in starch mass during the year of measurements in each wood depth

Ab_woodStarch_radial_activity=Ab_woodStarch_means_by_tree%>%group_by(ID, depth_mm)%>%
  summarise(activ_mean=max(mean_starch_area)-min(mean_starch_area),
            rel_activ_mean=(max(mean_starch_area)-min(mean_starch_area))/max(mean_starch_area)*100,
            activ_median=max(median_starch_are)-min(median_starch_are),
            activ_IQR=max(Int_quart_range)-min(Int_quart_range),
            max_starch=max(mean_starch_area),
            max_month=.data$month[max(which(.data$mean_starch_area==max(.data$mean_starch_area)))],
            min_starch=min(mean_starch_area),
            min_month=.data$month[min(which(.data$mean_starch_area==min(.data$mean_starch_area)))]
  )

Ab_woodStarch_radial_activity=na.omit(Ab_woodStarch_radial_activity)
Ab_woodStarch_radial_activity=Ab_woodStarch_radial_activity%>%filter(ID!="ab26820")

Ab_woodStarch_radial_Activity_by_species=Ab_woodStarch_radial_activity%>%group_by(depth_mm)%>%
  summarise(mean_active_mean=mean(activ_mean),
            rel_mean_active_mean=mean(rel_activ_mean, na.rm=T),
            SD_rel_mean_active_mean=sd(rel_activ_mean, na.rm=T),
            mean_active_median=mean(activ_median),
            mean_active_IQR=mean(activ_IQR),
            SD_active_mean=sd(activ_mean),
            SD_active_median=sd(activ_median),
            Sd_active_IQR=sd(activ_IQR), 
            max=max(max_starch),
            mean_max=mean(max_starch),
            max_month_mode=names(which.max(table(max_month))),
            min=min(min_starch),
            mean_min=mean(min_starch),
            min_month_mode=names(which.max(table(min_month)))
  )

png("metabolic_activity_O.leu.png", width = 800, height=600)
ggplot(data=Ab_woodStarch_radial_Activity_by_species, aes(x=depth_mm, y=mean_active_mean))+
  geom_line(data=Ab_woodStarch_radial_activity, aes(x=depth_mm, y=activ_mean, factor=ID), col="gray95")+
  geom_line( col="turquoise3")+
  geom_ribbon( aes(ymin=mean_active_mean-SD_active_mean, ymax=mean_active_mean+SD_active_mean,
                   fill= "var"),
               fill="turquoise3",
               alpha=0.3)+
  ggtitle("O. leucoxylon\n(evergree/parenchyma-storing-species)")+
  xlab("depth (mm)")+ylab("seasonal amplitud starch %")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(4)),
        axis.text.y=element_text(size = rel(4)),
        axis.title.x = element_text(size = rel(4)),
        axis.title.y = element_text(size = rel(4)),
        strip.text = element_text(size = rel(4)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

png("relative_sorage_change_O.leu.png", width = 800, height = 600)
ggplot(data=Ab_woodStarch_radial_Activity_by_species, aes(x=depth_mm, y=rel_mean_active_mean))+
  geom_line(data=Ab_woodStarch_radial_activity, aes(x=depth_mm, y=rel_activ_mean, factor=ID), col="gray95")+
  geom_line( col="turquoise3")+
  geom_ribbon( aes(ymin=rel_mean_active_mean-SD_rel_mean_active_mean, ymax=rel_mean_active_mean+SD_rel_mean_active_mean,
                   fill= "var"),
               fill="turquoise3",
               alpha=0.3)+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  xlab("depth (mm)")+ylab("relative seasonal amplitud")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(4)),
        axis.text.y=element_text(size = rel(4)),
        axis.title.x = element_text(size = rel(4)),
        axis.title.y = element_text(size = rel(4)),
        strip.text = element_text(size = rel(4)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

## mean content of starch in the entire core 
## https://www.fao.org/3/w4095e/w4095e0c.htm- source of the wood density Ocotea leucoxylon 0.45 g/cm3

Ab_woodStarch_means_by_tree$segment_volume_cm3=pi*(0.5/2)^2*0.5 ### this in in cm3

Ab_woodStarch_means_by_tree$grams_wood=Ab_woodStarch_means_by_tree$segment_volume_cm3*0.45
Ab_woodStarch_means_by_tree$grams_starch=Ab_woodStarch_means_by_tree$grams_wood*Ab_woodStarch_means_by_tree$mean_starch_area/100


Ab_woodStarch_means_entire_tree=Ab_woodStarch_means_by_tree%>%group_by(month, ID)%>%
  summarise(mean_starch_percentage=mean(mean_starch_area),
            median_starch_percentage=median(median_starch_are),
            mean_IQR=mean(Int_quart_range),
            starch_woodcore_content_by_mean=sum(mean_starch_area),
            starch_woodcore_content_by_median=sum(median_starch_are),
            starch_woodcore_content_grams=sum(grams_starch))### estimate the variansby avegarng the 

Ab_woodStarch_means_entire_tree$starch_woodcore_content_mg=Ab_woodStarch_means_entire_tree$starch_woodcore_content_grams*1000



Ab_woodStarch_means_entire_tree_by_specie=Ab_woodStarch_means_entire_tree%>%group_by(month)%>%
  summarise(mean_starch_percentage=mean(mean_starch_percentage),
            sd_mean_starch_percentage=sd(mean_starch_percentage),
            median_starch_percentage=mean(median_starch_percentage),
            sd_median_starch_percentage=sd(median_starch_percentage),
            mean_IQR=mean(mean_IQR), 
            mean_starch_wood_content_by_mean=mean(starch_woodcore_content_by_mean),
            sd_starch_wood_content_by_mean=sd(starch_woodcore_content_by_mean),
            mean_starch_wood_content_by_median=mean(starch_woodcore_content_by_median),
            sd_starch_wood_content_by_medial=sd(starch_woodcore_content_by_median),
            mean_starch_wood_content_grams=mean(starch_woodcore_content_grams),
            sd_mean_starch_wood_content_grams=sd(starch_woodcore_content_grams),
            mean_starch_wood_content_mg=mean(starch_woodcore_content_mg)
  )

png("Ab_starch_wood_content_comp_month.png", width = 800, height = 600)
ggplot(Ab_woodStarch_means_entire_tree, aes(x=month, y=starch_woodcore_content_mg))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_point(data=Ab_woodStarch_means_entire_tree_by_specie, aes(x=month, y=mean_starch_wood_content_mg),  col="black", pch=17, size=3)+
  xlab("")+
  ylab("Starch woodcore content (mg)")+
  ggtitle("O. leucoxylon\n(evergree/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


### pairwised comparison with paired samples for the content of sugar in the entire tree. The differences between months are lost mainly because the influence of the lack 
## of seasonality deeper in the stemwood or the contraty seasonality, and then all the averaging and the integration. 

Ab_pairwise_comparison_entire_tree=Ab_woodStarch_means_entire_tree
Ab_pairwise_comparison_entire_tree$month2=as.factor(as.numeric(Ab_pairwise_comparison_entire_tree$month))
Ab_pairwise_comparison_entire_tree=Ab_pairwise_comparison_entire_tree[-c(which(Ab_pairwise_comparison_entire_tree$ID=="ab26820")),]
names(Ab_pairwise_comparison_entire_tree)
Ab_pairwise_comparison_entire_tree=Ab_pairwise_comparison_entire_tree[ ,c(2,10,9)]
Ab_pairwise_comparison_entire_tree$ID=as.factor(Ab_pairwise_comparison_entire_tree$ID)
str(Ab_pairwise_comparison_entire_tree)

ab_pwc_starch_mg <- Ab_pairwise_comparison_entire_tree %>%
  pairwise_t_test(starch_woodcore_content_mg ~ month2, paired = T,
                  p.adjust.method = "BH"
  )
ab_pwc_starch_mg


wilcox_test(starch_woodcore_content_mg~month2, data=Ab_pairwise_comparison_entire_tree, paired = T)


#### add growth data to AB_woodStarch_content_woodcore and compare with total starch, mean starch, starch in the first two cm... 

monthly_growth2_summarized=monthly_growth2%>%filter(date3%in%c(8:20))%>%group_by(ID)%>%
  summarise(annual_growth_2019=sum(rate_growth),
            mean_monthly_growth_2019=mean(rate_growth),
            max_monthly_growth_2019=max(rate_growth),
            min_monthly_growth_2019=min(rate_growth)
  )


Ab_woodStarch_entire_tree=merge(Ab_woodStarch_means_entire_tree, monthly_growth2, by=c("ID", "month"))
Ab_woodStarch_entire_tree=merge(Ab_woodStarch_entire_tree, monthly_growth2_summarized, by=c("ID"))

Ab_woodStarch_entire_tree$significant_rel=as.character(Ab_woodStarch_entire_tree$month)
key_months_1=unique(as.character(Ab_woodStarch_entire_tree$month))
ID_key_signif2=c("Not_significant",
                "Not_significant",
                "Not_significant",
                "Not_significant")

library(plyr)

Ab_woodStarch_entire_tree$significant_rel=mapvalues(Ab_woodStarch_entire_tree$significant_rel,
                                                     key_months_1,
                                                     ID_key_signif2)
detach("package:plyr", unload=TRUE)


png("annual_growth_vs_starch_content_O_leu.png",
    width = 800, height = 600)
ggplot(Ab_woodStarch_entire_tree, aes(x=starch_woodcore_content_mg, y=annual_growth_2019))+
  geom_point(color="turquoise2")+
  geom_smooth(aes(linetype=significant_rel), method = "lm", se=F, color="turquoise2", linetype = "dashed")+
  stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE, size=8, xpos = 1,
                             ypos = 1.1, xpos2 = 1,
                             ypos2= 1)+
  xlab("starch content (mg)")+ ylab("annual growth 2019 (cm)")+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  facet_wrap(~month)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(2)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


# changes in starch content between the months. 

Ab_woodStarch_entire_tree=Ab_woodStarch_entire_tree[order(Ab_woodStarch_entire_tree$ID, Ab_woodStarch_entire_tree$date3),]
Ab_woodStarch_entire_tree$starch_change=NA
Ab_woodStarch_entire_tree$starch_change_percent=NA
Ab_woodStarch_entire_tree$rel_starch_change=NA

Ab_woodStarch_entire_tree2=Ab_woodStarch_entire_tree%>%group_by(ID)%>%
  mapGroups(function(x){
    for(i in 2:nrow(x)){
      x$starch_change[i]=x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1]
      x$starch_change_percent[i]=x$mean_starch_percentage[i]-x$mean_starch_percentage[i-1]
      x$rel_starch_change[i]=(x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1])/x$starch_woodcore_content_mg[i]
    }
    return(x)
  }
  )%>%ungroup()


## pairwise comparison between changes in the starch content between seasons 

png("Ab_starch_mass_change_comp_month.png", width = 800, height = 600)
ggplot(Ab_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=starch_change))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (mg)")+
  ggtitle("O. leucoxylon\n(evergree/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


png("Ab_starch_mass_change_rel_comp_month.png", width = 800, height = 600)
ggplot(Ab_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=rel_starch_change))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (%)")+
  ylim(-4, 2)+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

png("Ab_starch_mass_change_percent_month.png", width = 800, height = 600)
ggplot(Ab_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=starch_change_percent))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (%)")+
  ylim(-4, 2)+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

Ab_woodStarch_entire_tree2_wilcox_pair_comp=filter(Ab_woodStarch_entire_tree2, month!="May19")
Ab_woodStarch_entire_tree2_wilcox_pair_comp$month=factor(Ab_woodStarch_entire_tree2_wilcox_pair_comp$month,
                                                         levels = c("Aug19", "Nov19", "Feb20"))
Ab_woodStarch_entire_tree2_wilcox_pair_comp$month2=as.factor(as.numeric(Ab_woodStarch_entire_tree2_wilcox_pair_comp$month))


Ab_woodStarch_entire_tree2_wilcox_pair_comp=Ab_woodStarch_entire_tree2_wilcox_pair_comp%>%
  filter(ID!="ab26820")
names(Ab_woodStarch_entire_tree2_wilcox_pair_comp)
Ab_woodStarch_entire_tree2_wilcox_pair_comp=Ab_woodStarch_entire_tree2_wilcox_pair_comp[,c("ID", "month2", "starch_change")]
Ab_woodStarch_entire_tree2_wilcox_pair_comp$ID=as.factor(Ab_woodStarch_entire_tree2_wilcox_pair_comp$ID)
str(Ab_woodStarch_entire_tree2_wilcox_pair_comp)

wilcox_test(starch_change~month2, data=Ab_woodStarch_entire_tree2_wilcox_pair_comp, paired = T)

### test for the relationship growth vs starch changes

Ab_woodStarch_entire_tree2=filter(Ab_woodStarch_entire_tree2,month!="May19")
Ab_woodStarch_entire_tree2$month_changes=as.character(Ab_woodStarch_entire_tree2$month)
month_chages_key=unique(as.character(Ab_woodStarch_entire_tree2$month))
ID_key=c("starch (May19-Aug19)\ngrowth (Aug19-Nov19)", 
         "starch (Aug19-Nov19)\ngrowth (Nov19-Feb20)", 
         "starch (Nov19-Feb20)\n growth (Feb20-May20)")

Ab_woodStarch_entire_tree2$significant_rel=as.character(Ab_woodStarch_entire_tree2$month)
ID_key_signif=c("Not_significant",
                "significant",
                "Not_significant")

library(plyr)
Ab_woodStarch_entire_tree2$month_changes=mapvalues(Ab_woodStarch_entire_tree2$month_changes,
                                                   month_chages_key,
                                                   ID_key)
Ab_woodStarch_entire_tree2$significant_rel=mapvalues(Ab_woodStarch_entire_tree2$significant_rel,
                                                   month_chages_key,
                                                   ID_key_signif)
detach("package:plyr", unload=TRUE)

Ab_woodStarch_entire_tree2$month_changes=factor(Ab_woodStarch_entire_tree2$month_changes,
                                                levels = c("starch (May19-Aug19)\ngrowth (Aug19-Nov19)", 
                                                           "starch (Aug19-Nov19)\ngrowth (Nov19-Feb20)", 
                                                           "starch (Nov19-Feb20)\n growth (Feb20-May20)"))

png("starch_change_and_growth_O_leu.png",
    width = 800, height = 600)
ggplot(Ab_woodStarch_entire_tree2, aes(x=starch_change, y=rate_growth_summ_three_months_after))+
  geom_point(color="turquoise2", size=4, alpha=0.5)+
  geom_smooth(aes(linetype=significant_rel),method = "lm", se=F, color="turquoise2")+
  scale_linetype_manual(values=c("dashed", "solid"))+ 
  stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE, size=7, xpos = -15,
                             ypos = 0.7, xpos2 = -15,
                             ypos2= 0.65)+
  #geom_hline(yintercept = 0, linetype="dashed")+
  geom_vline(xintercept = 0, linetype="dashed")+
  xlab("starch change (mg)")+ ylab("accumulated growth (cm)")+
  ggtitle("O. leucoxylon\n(evergreen/parenchyma-storing-species)")+
  facet_wrap(~month_changes)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(2)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),
        plot.title = element_text(size = rel(2.5), face = "italic"),
        strip.text.x=element_text(size = 20, margin = margin(10), vjust = 1))
dev.off()



##### Sacoglotis 

load("Sac_profile_data") 


names(Sac_profile_data_Jan18)
str(Sac_profile_data_Jan18)
names(Sac_profile_data_Aug_2019)

### selecting the varSacles that are gonna be used in this analysis

Sac_May19=Sac_profile_data_May_2019%>%select("ID", "depth_mm", "measurement_type", "X.Area", "sampling_month", "month")
Sac_Aug19=Sac_profile_data_Aug_2019%>%select("ID", "depth_mm", "measurement_type", "X.Area", "sampling_month", "month")
Sac_Aug19$month="Aug19"
Sac_Nov19=Sac_profile_data_Nov_2019%>%select("ID", "depth_mm", "measurement_type", "X.Area", "sampling_month", "month")
Sac_feb19=Sac_profile_data_Feb_2020%>%select("ID", "depth_mm", "measurement_type", "X.Area", "sampling_month", "month")

### unify the data frame for the species and choosing the measurement type that will be analyzed (woodstarch in this case)

Sac=rbind(Sac_May19, Sac_Aug19, Sac_Nov19, Sac_feb19)
unique(Sac$ID)
Sac=Sac%>%separate(ID, c("ID", "radii_No"))
Sac$IDspecies=gsub("[[:digit:]]","",Sac$ID)
Sac$IDspecies="Sac"
Sac$IDnumber=gsub("[^[:digit:]]","", Sac$ID)
Sac$ID=paste(Sac$IDspecies, Sac$IDnumber, sep="")
Sac$ID[which(Sac$ID=="Sac34283")]="Sac34288"
Sac$ID[which(Sac$ID=="Sac27376")]="Sac25376"

Sac_woodStarch=filter(Sac, measurement_type=="woodStarch")
Sac_woodStarch$depth_mm=as.numeric(Sac_woodStarch$depth_mm)

### generate the cero values for the inner layesr of wood that were not measured because the Sacsence of starch 

Sac_woodStarch=Sac_woodStarch %>% group_by(ID, month) %>% 
  mapGroups(function(groups_ID){
    max_depth_spp=max(Sac_woodStarch$depth_mm)
    max_depth_grou=max(groups_ID$depth_mm)
    while(max_depth_grou<max_depth_spp){
      i=length(groups_ID$depth_mm)
      groups_ID[c(i+1, i+2, i+3), "depth_mm"]=max(groups_ID$depth_mm)+5
      groups_ID$ID[c(i+1, i+2, i+3)]=groups_ID$ID[1]
      groups_ID$measurement_type="woodStarch"
      groups_ID$month=groups_ID$month[1]
      groups_ID$sampling_month[c(i+1, i+2, i+3)]=groups_ID$sampling_month[1]
      max_depth_grou=max(groups_ID$depth_mm)
    }
    groups_ID$X.Area[is.na(groups_ID$X.Area)]=0
    return(groups_ID)
  }
  ) %>% ungroup() 

save(Sac_woodStarch, file = "Sac_woodStarch.Rda")
save(Sac_woodStarch, file = "Sac_woodStarch.csv")
############## Calculating mean statistics per individual, month and sampling depth ####################

## mean statistics by tree 

Sac_woodStarch_means_by_tree=Sac_woodStarch%>%group_by(ID, depth_mm, month, sampling_month)%>%
  summarize(mean_starch_area=mean(X.Area, na.rm=T),
            median_starch_are=median(X.Area, na.rm=T),
            Int_quart_range=IQR(X.Area, na.rm=T),
            sd_starch_area=sd(X.Area, na.rm=T),
            sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area)))

table(Sac_woodStarch_means_by_tree$ID)
### mean statistics by species 

Sac_woodStarch_means_by_specie=Sac_woodStarch%>%group_by(depth_mm, month, sampling_month)%>%
  summarize(mean_starch_area=mean(X.Area, na.rm=T),
            median_starch_are=median(X.Area, na.rm=T),
            Int_quart_range=IQR(X.Area, na.rm=T),
            sd_starch_area=sd(X.Area, na.rm=T),
            sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area)))

Sac_woodStarch_means_by_tree$month=factor(Sac_woodStarch_means_by_tree$month, 
                                          levels = c("May19", "Aug19", "Nov19", "Feb20"))
Sac_woodStarch_means_by_specie$month=factor(Sac_woodStarch_means_by_specie$month, 
                                            levels = c("May19", "Aug19", "Nov19", "Feb20"))

png("radial_starch_content_S-gui.png", width = 900, height = 600)
ggplot()+
  geom_line(data=filter(Sac_woodStarch_means_by_specie,month!="Jan18" && month !="Jul18") , aes(x=depth_mm, y=mean_starch_area, col=month), size=1)+
  #facet_wrap(~month)+
  #geom_line(data=filter(Sac_woodStarch_means_by_tree, month!="Jan18" && month !="Jul18"), aes(x=depth_mm, y=mean_starch_area, col=ID), col="gray95")+
  geom_ribbon(data=filter(Sac_woodStarch_means_by_specie,month!="Jan18" && month !="Jul18"), aes(x=depth_mm, 
                                                                                                 ymin=abs(mean_starch_area-sd_starch_area),
                                                                                                 ymax=mean_starch_area+sd_starch_area,
                                                                                                 fill= month),
              #fill="gray80",
              alpha=0.1)+
  
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  xlab("depth (mm)")+ylab("starch %")+
  #xlim(0,110)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(3)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

## statistical comparisons at different depth ranges (ANOVA)

Sac_pairwise_comparison=Sac_woodStarch_means_by_tree
Sac_pairwise_comparison=Sac_pairwise_comparison[, c(1,2,3,5)]
Sac_pairwise_comparison$time=as.factor(Sac_pairwise_comparison$month)
Sac_pairwise_comparison=Sac_pairwise_comparison[order(Sac_pairwise_comparison$month),]
rownames(Sac_pairwise_comparison)=NULL

Sac_pairwise_comparison$ID=as.factor(Sac_pairwise_comparison$ID)
Sac_Anva2=Sac_pairwise_comparison%>%group_by(depth_mm)%>%nest()%>%
  mutate(ANOVA=map(data, ~ aov(mean_starch_area~month, data=.x)),
         tukey_comparison=map(ANOVA, TukeyHSD),
         tidied=map(ANOVA, tidy),
         tidied2=map(tukey_comparison, tidy), 
  )

names(Sac_Anva2$tidied)=Sac_Anva2$depth_mm
names(Sac_Anva2$tidied2)=Sac_Anva2$depth_mm
Sac_Anva2$tidied
Sac_Anva2$tidied2
estimates_ANOVA_by_depth=Sac_Anva2%>%unnest(tidied)

ggplot(Sac_woodStarch_means_by_tree, aes(x=month, y=mean_starch_area))+
  geom_boxplot(data=Sac_woodStarch_means_by_tree, aes(x=month, y=mean_starch_area))+
  geom_point(aes(col=ID))+
  geom_point(data=Sac_woodStarch_means_by_specie, aes(x=month, y=mean_starch_area), col="black")+
  facet_wrap(~depth_mm)


### estimate the maximum and relative seasonal range in starch mass during the year of measurements in each wood depth

Sac_woodStarch_radial_activity=Sac_woodStarch_means_by_tree%>%group_by(ID, depth_mm)%>%
  summarise(activ_mean=max(mean_starch_area)-min(mean_starch_area),
            rel_activ_mean=(max(mean_starch_area)-min(mean_starch_area))/max(mean_starch_area)*100,
            activ_median=max(median_starch_are)-min(median_starch_are),
            activ_IQR=max(Int_quart_range)-min(Int_quart_range),
            max_starch=max(mean_starch_area),
            max_month=.data$month[max(which(.data$mean_starch_area==max(.data$mean_starch_area)))],
            min_starch=min(mean_starch_area),
            min_month=.data$month[min(which(.data$mean_starch_area==min(.data$mean_starch_area)))]
  )

Sac_woodStarch_radial_activity=na.omit(Sac_woodStarch_radial_activity)

Sac_woodStarch_radial_Activity_by_species=Sac_woodStarch_radial_activity%>%group_by(depth_mm)%>%
  summarise(mean_active_mean=mean(activ_mean, na.rm=T),
            rel_mean_active_mean=mean(rel_activ_mean, na.rm=T),
            SD_rel_mean_active_mean=sd(rel_activ_mean, na.rm=T),
            mean_active_median=mean(activ_median),
            mean_active_IQR=mean(activ_IQR),
            SD_active_mean=sd(activ_mean),
            SD_active_median=sd(activ_median),
            Sd_active_IQR=sd(activ_IQR), 
            max=max(max_starch),
            mean_max=mean(max_starch),
            max_month_mode=names(which.max(table(max_month))),
            min=min(min_starch),
            mean_min=mean(min_starch),
            min_month_mode=names(which.max(table(min_month)))
  )

png("metabolic_activity_Sac.png", width = 800, height=600)
ggplot(data=Sac_woodStarch_radial_Activity_by_species, aes(x=depth_mm, y=mean_active_mean))+
  geom_line(data=Sac_woodStarch_radial_activity, aes(x=depth_mm, y=activ_mean, factor=ID), col="gray95")+
  geom_line( col="turquoise3")+
  geom_ribbon( aes(ymin=mean_active_mean-SD_active_mean, ymax=mean_active_mean+SD_active_mean,
                   fill= "var"),
               fill="turquoise3",
               alpha=0.3)+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  xlab("depth (mm)")+ylab("seasonal amplitud starch %")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(4)),
        axis.text.y=element_text(size = rel(4)),
        axis.title.x = element_text(size = rel(4)),
        axis.title.y = element_text(size = rel(4)),
        strip.text = element_text(size = rel(4)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

png("relative_sorage_change_Sac.png", width = 900, height = 600)
ggplot(data=Sac_woodStarch_radial_Activity_by_species, aes(x=depth_mm, y=rel_mean_active_mean))+
  geom_line(data=Sac_woodStarch_radial_activity, aes(x=depth_mm, y=rel_activ_mean, factor=ID), col="gray95")+
  geom_line( col="turquoise3")+
  geom_ribbon( aes(ymin=rel_mean_active_mean-SD_rel_mean_active_mean, ymax=rel_mean_active_mean+SD_rel_mean_active_mean,
                   fill= "var"),
               fill="turquoise3",
               alpha=0.3)+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  xlab("depth (mm)")+ylab("relative seasonal amplitud")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(4)),
        axis.text.y=element_text(size = rel(4)),
        axis.title.x = element_text(size = rel(4)),
        axis.title.y = element_text(size = rel(4)),
        strip.text = element_text(size = rel(4)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

## mean content of starch in the entire core 
## wood density 0.81 g/cm3 was taken from the recors from tanguro

Sac_woodStarch_means_by_tree$segment_volume_cm3=pi*(0.5/2)^2*0.5 ### this in in cm3

Sac_woodStarch_means_by_tree$grams_wood=Sac_woodStarch_means_by_tree$segment_volume_cm3*0.8
Sac_woodStarch_means_by_tree$grams_starch=Sac_woodStarch_means_by_tree$grams_wood*Sac_woodStarch_means_by_tree$mean_starch_area/100
Sac_woodStarch_means_by_tree$mg_starch=Sac_woodStarch_means_by_tree$grams_starch*1000

Sac_woodStarch_means_entire_tree=Sac_woodStarch_means_by_tree%>%group_by(month, ID)%>%
  summarise(mean_starch_percentage=mean(mean_starch_area),
            median_starch_percentage=mean(median_starch_are),
            mean_IQR=mean(Int_quart_range),
            starch_woodcore_content_by_mean=sum(mean_starch_area),
            starch_woodcore_content_by_median=sum(median_starch_are),
            starch_woodcore_content_grams=sum(grams_starch),
            starch_woodcore_content_mg=sum(mg_starch))


Sac_woodStarch_means_entire_tree=Sac_woodStarch_means_by_tree%>%group_by(month, ID)%>%
  summarise(mean_starch_percentage=mean(mean_starch_area),
            median_starch_percentage=median(median_starch_are),
            mean_IQR=mean(Int_quart_range),
            starch_woodcore_content_by_mean=sum(mean_starch_area),
            starch_woodcore_content_by_median=sum(median_starch_are),
            starch_woodcore_content_grams=sum(grams_starch))

Sac_woodStarch_means_entire_tree$starch_woodcore_content_mg=Sac_woodStarch_means_entire_tree$starch_woodcore_content_grams*1000



Sac_woodStarch_means_entire_tree_by_specie=Sac_woodStarch_means_entire_tree%>%group_by(month)%>%
  summarise(mean_starch_percentage=mean(mean_starch_percentage),
            sd_mean_starch_percentage=sd(mean_starch_percentage),
            median_starch_percentage=mean(median_starch_percentage),
            sd_median_starch_percentage=sd(median_starch_percentage),
            mean_IQR=mean(mean_IQR), 
            mean_starch_wood_content_by_mean=mean(starch_woodcore_content_by_mean),
            sd_starch_wood_content_by_mean=sd(starch_woodcore_content_by_mean),
            mean_starch_wood_content_by_median=mean(starch_woodcore_content_by_median),
            sd_starch_wood_content_by_medial=sd(starch_woodcore_content_by_median),
            mean_starch_wood_content_grams=mean(starch_woodcore_content_grams),
            sd_mean_starch_wood_content_grams=sd(starch_woodcore_content_grams),
            mean_starch_wood_content_mg=mean(starch_woodcore_content_mg)
  )
Sac_woodStarch_means_entire_tree$starch_woodcore_content_mg=Sac_woodStarch_means_entire_tree$starch_woodcore_content_grams*1000

png("Sac_starch_wood_content_comp_month.png", width = 800, height = 600)
ggplot(Sac_woodStarch_means_entire_tree, aes(x=month, y=starch_woodcore_content_mg))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_point(data=Sac_woodStarch_means_entire_tree_by_specie, aes(x=month, y=mean_starch_wood_content_mg),  col="black", pch=17, size=3)+
  xlab("")+
  ylab("Starch woodcore content (mg)")+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


### pairwised comparison with paired samples for the content of sugar in the entire tree. The differences between months are lost mainly because the influence of the lack 
## of seasonality deeper in the stemwood or the contraty seasonality, and then all the averaging and the integration. 

Sac_pairwise_comparison_entire_tree=Sac_woodStarch_means_entire_tree
Sac_pairwise_comparison_entire_tree$month2=as.factor(as.numeric(Sac_pairwise_comparison_entire_tree$month))
#Sac_pairwise_comparison_entire_tree=Sac_pairwise_comparison_entire_tree[-c(which(Sac_pairwise_comparison_entire_tree$ID=="Sac26820")),]
Sac_pairwise_comparison_entire_tree=Sac_pairwise_comparison_entire_tree[-c(which(Sac_pairwise_comparison_entire_tree$ID=="Sac1054")),]
names(Sac_pairwise_comparison_entire_tree)
Sac_pairwise_comparison_entire_tree=Sac_pairwise_comparison_entire_tree[ ,c(2,10,9)]
Sac_pairwise_comparison_entire_tree$ID=as.factor(Sac_pairwise_comparison_entire_tree$ID)
str(Sac_pairwise_comparison_entire_tree)

Sac_pwc_starch_mg <- Sac_pairwise_comparison_entire_tree %>%
  pairwise_t_test(starch_woodcore_content_mg ~ month2, paired = T,
                  p.adjust.method = "BH"
  )
Sac_pwc_starch_mg


wilcox_test(starch_woodcore_content_mg~month2, data=Sac_pairwise_comparison_entire_tree, paired = T)


#### add growth data to Sac_woodStarch_content_woodcore and compare with total starch, mean starch, starch in the first two cm... 

monthly_growth2_summarized=monthly_growth2%>%filter(date3%in%c(8:20))%>%group_by(ID)%>%
  summarise(annual_growth_2019=sum(rate_growth),
            mean_monthly_growth_2019=mean(rate_growth),
            max_monthly_growth_2019=max(rate_growth),
            min_monthly_growth_2019=min(rate_growth)
  )


Sac_woodStarch_entire_tree=merge(Sac_woodStarch_means_entire_tree, monthly_growth2, by=c("ID", "month"))
Sac_woodStarch_entire_tree=merge(Sac_woodStarch_entire_tree, monthly_growth2_summarized, by=c("ID"))

Sac_woodStarch_entire_tree$significant_rel=as.character(Sac_woodStarch_entire_tree$month)
key_months_1=unique(as.character(Sac_woodStarch_entire_tree$month))
ID_key_signif2=c("Not_significant",
                 "significant",
                 "Not_significant",
                 "significant")

library(plyr)

Sac_woodStarch_entire_tree$significant_rel=mapvalues(Sac_woodStarch_entire_tree$significant_rel,
                                                    key_months_1,
                                                    ID_key_signif2)
detach("package:plyr", unload=TRUE)

png("annual_growth_vs_starch_content_Sac.png",
    width = 800, height = 600)
ggplot(Sac_woodStarch_entire_tree, aes(x=starch_woodcore_content_mg, y=annual_growth_2019))+
  geom_point(color="turquoise2")+
  geom_smooth(aes(linetype=significant_rel), method = "lm", se=F, color="turquoise2")+
  stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE, size=8, xpos = 1,
                             ypos = 2.4, xpos2 = 1,
                             ypos2= 2.1)+
  xlab("starch content (mg)")+ ylab("annual growth 2019 (cm)")+
  ylim(0, 2.5)+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  facet_wrap(~month)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(2)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


# changes in starch content between the months. 

Sac_woodStarch_entire_tree=Sac_woodStarch_entire_tree[order(Sac_woodStarch_entire_tree$ID, Sac_woodStarch_entire_tree$date3),]
Sac_woodStarch_entire_tree$starch_change=NA
Sac_woodStarch_entire_tree$starch_change_percent=NA
Sac_woodStarch_entire_tree$rel_starch_change=NA

Sac_woodStarch_entire_tree2=Sac_woodStarch_entire_tree%>%group_by(ID)%>%
  mapGroups(function(x){
    for(i in 2:nrow(x)){
      x$starch_change[i]=x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1]
      x$starch_change_percent[i]=x$mean_starch_percentage[i]-x$mean_starch_percentage[i-1]
      x$rel_starch_change[i]=(x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1])/x$starch_woodcore_content_mg[i]
    }
    return(x)
  }
  )%>%ungroup()


## pairwise comparison between changes in the starch content between seasons 

png("Sac_starch_mass_change_comp_month.png", width = 800, height = 600)
ggplot(Sac_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=starch_change))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (mg)")+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


png("Sac_starch_mass_change_rel_comp_month.png", width = 800, height = 600)
ggplot(Sac_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=rel_starch_change))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (%)")+
  #ylim(-4, 2)+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

png("Sac_starch_mass_change_percentage_month.png", width = 800, height = 600)
ggplot(Sac_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=starch_change_percent))+geom_boxplot(fill="turquoise2")+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (%)")+
  #ylim(-4, 2)+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

Sac_woodStarch_entire_tree2_wilcox_pair_comp=filter(Sac_woodStarch_entire_tree2, month!="May19")
Sac_woodStarch_entire_tree2_wilcox_pair_comp$month=factor(Sac_woodStarch_entire_tree2_wilcox_pair_comp$month,
                                                          levels = c("Aug19", "Nov19", "Feb20"))
Sac_woodStarch_entire_tree2_wilcox_pair_comp$month2=as.factor(as.numeric(Sac_woodStarch_entire_tree2_wilcox_pair_comp$month))
names(Sac_woodStarch_entire_tree2_wilcox_pair_comp)
Sac_woodStarch_entire_tree3_wilcox_pair_comp=Sac_woodStarch_entire_tree2_wilcox_pair_comp[,c("ID", "month2", "starch_change")]
Sac_woodStarch_entire_tree4_wilcox_pair_comp=Sac_woodStarch_entire_tree2_wilcox_pair_comp[,c("ID", "month2", "rel_starch_change")]
Sac_woodStarch_entire_tree3_wilcox_pair_comp$ID=as.factor(Sac_woodStarch_entire_tree3_wilcox_pair_comp$ID)
Sac_woodStarch_entire_tree4_wilcox_pair_comp$ID=as.factor(Sac_woodStarch_entire_tree4_wilcox_pair_comp$ID)
str(Sac_woodStarch_entire_tree3_wilcox_pair_comp)
str(Sac_woodStarch_entire_tree4_wilcox_pair_comp)

wilcox_test(starch_change~month2, data=Sac_woodStarch_entire_tree3_wilcox_pair_comp, paired = T)

wilcox_test(rel_starch_change~month2, data=Sac_woodStarch_entire_tree4_wilcox_pair_comp, paired = T)


### test for the relationship growth vs starch changes

Sac_woodStarch_entire_tree2=filter(Sac_woodStarch_entire_tree2,month!="May19")
Sac_woodStarch_entire_tree2$month_changes=as.character(Sac_woodStarch_entire_tree2$month)
month_chages_key=unique(as.character(Sac_woodStarch_entire_tree2$month))
ID_key=c("starch (May19-Aug19)\ngrowth (Aug19-Nov19)", 
         "starch (Aug19-Nov19)\ngrowth (Nov19-Feb20)", 
         "starch (Nov19-Feb20)\n growth (Feb20-May20)")


library(plyr)
Sac_woodStarch_entire_tree2$month_changes=mapvalues(Sac_woodStarch_entire_tree2$month_changes,
                                                   month_chages_key,
                                                   ID_key)
detach("package:plyr", unload=TRUE)

Sac_woodStarch_entire_tree2$month_changes=factor(Sac_woodStarch_entire_tree2$month_changes,
                                                levels = c("starch (May19-Aug19)\ngrowth (Aug19-Nov19)", 
                                                           "starch (Aug19-Nov19)\ngrowth (Nov19-Feb20)", 
                                                           "starch (Nov19-Feb20)\n growth (Feb20-May20)"))

png("three_month_growth_after_vs_starch_change_S_gui.png",
    width = 800, height = 600)
ggplot(Sac_woodStarch_entire_tree2, aes(x=starch_change, y=rate_growth_summ_three_months_after))+
  geom_point(color="turquoise2", size=4, alpha=0.5)+
  #geom_text(aes(label=ID))+
  geom_smooth(linetype=c("dashed"), method = "lm", se=F, color="turquoise")+
  stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE, size=7, xpos = -75,
                             ypos = 0.66, xpos2 = -75,
                             ypos2= 0.63)+
  xlab("starch change (mg)")+ylab("accumulated growth (cm)")+
  ggtitle("S. guianensis\n(semi-deciduous/parenchyma-storing-species)")+
  ylim(c(0.0, 0.7))+
  geom_vline(xintercept = 0, linetype="dashed")+
  facet_wrap(~month_changes)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(3)),
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),
        legend.position="none",
        plot.title = element_text(size = rel(2.5), face = "italic"),
        strip.text.x = element_text(size= 20, margin = margin(10), vjust = 1),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()



############################################################################################################################################################


##### Dacryodes

load("Dac_profile_data")

names(Dac_profile_data_Jan18)
str(Dac_profile_data_Jan18)
names(Dac_profile_data_Aug_2019)

### selecting the varDacles that are gonna be used in this analysis

Dac_May19=Dac_profile_data_May_2019%>%select("ID", "depth_mm", "measurement_type", "X.Area", "sampling_month", "month")
Dac_Aug19=Dac_profile_data_Aug_2019%>%select("ID", "depth_mm", "measurement_type", "X.Area", "sampling_month", "month")
Dac_Aug19$month="Aug19"
Dac_Nov19=Dac_profile_data_Nov_2019%>%select("ID", "depth_mm", "measurement_type", "X.Area", "sampling_month", "month")
Dac_feb19=Dac_profile_data_Feb_2020%>%select("ID", "depth_mm", "measurement_type", "X.Area", "sampling_month", "month")

### unify the data frame for the species and choosing the measurement type that will be analyzed (woodstarch in this case)

Dac=rbind(Dac_May19, Dac_Aug19, Dac_Nov19, Dac_feb19)

unique(Dac$ID)
Dac=Dac%>%separate(ID, c("ID", "radii_No"))
Dac$IDspecies=gsub("[[:digit:]]","",Dac$ID)
Dac$IDspecies="Dac"
Dac$IDnumber=gsub("[^[:digit:]]","", Dac$ID)
Dac$ID=paste(Dac$IDspecies, Dac$IDnumber, sep="")

Dac$ID[which(Dac$ID=="Dac22924")]="Dac27924"
Dac$ID[which(Dac$ID=="Dac25629")]="Dac25624"

Dac_woodStarch=filter(Dac, measurement_type=="woodStarch")
Dac_woodStarch$depth_mm2=as.numeric(Dac_woodStarch$depth_mm)
Dac_woodStarch[which(is.na(Dac_woodStarch$depth_mm2)), c("depth_mm","depth_mm2")]=35
Dac_woodStarch=Dac_woodStarch[,-c(10)]
Dac_woodStarch$depth_mm=as.numeric(Dac_woodStarch$depth_mm)


### generate the cero values for the inner layesr of wood that were not measured because the Dacsence of starch 

Dac_woodStarch=Dac_woodStarch %>% group_by(ID, month) %>% 
  mapGroups(function(groups_ID){
    max_depth_spp=max(Dac_woodStarch$depth_mm)
    max_depth_grou=max(groups_ID$depth_mm)
    while(max_depth_grou<max_depth_spp){
      i=length(groups_ID$depth_mm)
      groups_ID[c(i+1, i+2, i+3), "depth_mm"]=max(groups_ID$depth_mm)+5
      groups_ID$ID[c(i+1, i+2, i+3)]=groups_ID$ID[1]
      groups_ID$measurement_type="woodStarch"
      groups_ID$month=groups_ID$month[1]
      groups_ID$sampling_month[c(i+1, i+2, i+3)]=groups_ID$sampling_month[1]
      max_depth_grou=max(groups_ID$depth_mm)
    }
    groups_ID$X.Area[is.na(groups_ID$X.Area)]=0
    return(groups_ID)
  }
  ) %>% ungroup() 


save(Dac_woodStarch, file="Dac_woodStarch.Rda")
save(Dac_woodStarch, file="Dac_woodStarch.csv")

############## Calculating mean statistics per individual, month and sampling depth ####################

## mean statistics by tree 

Dac_woodStarch_means_by_tree=Dac_woodStarch%>%group_by(ID, depth_mm, month, sampling_month)%>%
  summarize(mean_starch_area=mean(X.Area, na.rm=T),
            median_starch_are=median(X.Area, na.rm=T),
            Int_quart_range=IQR(X.Area, na.rm=T),
            sd_starch_area=sd(X.Area, na.rm=T),
            sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area)))

### mean statistics by species 

Dac_woodStarch_means_by_specie=Dac_woodStarch%>%group_by(depth_mm, month, sampling_month)%>%
  summarize(mean_starch_area=mean(X.Area, na.rm=T),
            median_starch_are=median(X.Area, na.rm=T),
            Int_quart_range=IQR(X.Area, na.rm=T),
            sd_starch_area=sd(X.Area, na.rm=T),
            sterr=var(X.Area, na.rm=T)/sqrt(length(X.Area)))

Dac_woodStarch_means_by_tree$month=factor(Dac_woodStarch_means_by_tree$month, 
                                          levels = c("May19", "Aug19", "Nov19", "Feb20"))
Dac_woodStarch_means_by_specie$month=factor(Dac_woodStarch_means_by_specie$month, 
                                            levels = c("May19", "Aug19", "Nov19", "Feb20"))

png("radial_starch_content_Dmic.png", width = 900, height = 600)
ggplot()+
  geom_line(data=filter(Dac_woodStarch_means_by_specie,month!="Jan18" && month !="Jul18") , aes(x=depth_mm, y=mean_starch_area, col=month), size=1)+
  #facet_wrap(~month)+
  #geom_line(data=filter(Dac_woodStarch_means_by_tree, month!="Jan18" && month !="Jul18"), aes(x=depth_mm, y=mean_starch_area, col=ID), col="gray95")+
  geom_ribbon(data=filter(Dac_woodStarch_means_by_specie,month!="Jan18" && month !="Jul18"), aes(x=depth_mm, 
                                                                                                 ymin=mean_starch_area-sd_starch_area,
                                                                                                 ymax=mean_starch_area+sd_starch_area,
                                                                                                 fill= month),
              #fill="gray80",
              alpha=0.1)+
  
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  xlab("depth (mm)")+ylab("starch %")+
  #xlim(0,110)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(3)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

## statistical comparisons at different depth ranges (ANOVA)

Dac_pairwise_comparison=Dac_woodStarch_means_by_tree
Dac_pairwise_comparison=Dac_pairwise_comparison[, c(1,2,3,5)]
Dac_pairwise_comparison$time=as.factor(Dac_pairwise_comparison$month)
Dac_pairwise_comparison=Dac_pairwise_comparison[order(Dac_pairwise_comparison$month),]
rownames(Dac_pairwise_comparison)=NULL

Dac_pairwise_comparison$ID=as.factor(Dac_pairwise_comparison$ID)
Dac_Anva2=Dac_pairwise_comparison%>%group_by(depth_mm)%>%nest()%>%
  mutate(ANOVA=map(data, ~ aov(mean_starch_area~month, data=.x)),
         tukey_comparison=map(ANOVA, TukeyHSD),
         tidied=map(ANOVA, tidy),
         tidied2=map(tukey_comparison, tidy), 
  )

names(Dac_Anva2$tidied)=Dac_Anva2$depth_mm
names(Dac_Anva2$tidied2)=Dac_Anva2$depth_mm
Dac_Anva2$tidied
Dac_Anva2$tidied2
estimates_ANOVA_by_depth=Dac_Anva2%>%unnest(tidied)

ggplot(Dac_woodStarch_means_by_tree, aes(x=month, y=mean_starch_area))+
  geom_boxplot(data=Dac_woodStarch_means_by_tree, aes(x=month, y=mean_starch_area))+
  geom_point(aes(col=ID))+
  geom_point(data=Dac_woodStarch_means_by_specie, aes(x=month, y=mean_starch_area), col="black")+
  facet_wrap(~depth_mm)


### estimate the maximum and relative seasonal range in starch mass during the year of measurements in each wood depth

Dac_woodStarch_radial_activity=Dac_woodStarch_means_by_tree%>%group_by(ID, depth_mm)%>%
  summarise(activ_mean=max(mean_starch_area)-min(mean_starch_area),
            rel_activ_mean=(max(mean_starch_area)-min(mean_starch_area))/max(mean_starch_area)*100,
            activ_median=max(median_starch_are)-min(median_starch_are),
            activ_IQR=max(Int_quart_range)-min(Int_quart_range),
            max_starch=max(mean_starch_area),
            max_month=.data$month[max(which(.data$mean_starch_area==max(.data$mean_starch_area)))],
            min_starch=min(mean_starch_area),
            min_month=.data$month[min(which(.data$mean_starch_area==min(.data$mean_starch_area)))]
  )

Dac_woodStarch_radial_activity=na.omit(Dac_woodStarch_radial_activity)
Dac_woodStarch_radial_activity=Dac_woodStarch_radial_activity%>%filter(ID!="Dac26820")

Dac_woodStarch_radial_Activity_by_species=Dac_woodStarch_radial_activity%>%group_by(depth_mm)%>%
  summarise(mean_active_mean=mean(activ_mean, na.rm=T),
            rel_mean_active_mean=mean(rel_activ_mean, na.rm=T),
            SD_rel_mean_active_mean=sd(rel_activ_mean, na.rm=T),
            mean_active_median=mean(activ_median),
            mean_active_IQR=mean(activ_IQR),
            SD_active_mean=sd(activ_mean),
            SD_active_median=sd(activ_median),
            Sd_active_IQR=sd(activ_IQR), 
            max=max(max_starch),
            mean_max=mean(max_starch),
            max_month_mode=names(which.max(table(max_month))),
            min=min(min_starch),
            mean_min=mean(min_starch),
            min_month_mode=names(which.max(table(min_month)))
  )

png("metabolic_activity_Dac.png", width = 800, height=600)
ggplot(data=Dac_woodStarch_radial_Activity_by_species, aes(x=depth_mm, y=mean_active_mean))+
  geom_line(data=Dac_woodStarch_radial_activity, aes(x=depth_mm, y=activ_mean, factor=ID), col="gray95")+
  geom_line( col="orangered2")+
  geom_ribbon( aes(ymin=mean_active_mean-SD_active_mean, ymax=mean_active_mean+SD_active_mean,
                   fill= "var"),
               fill="orangered2",
               alpha=0.3)+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  xlab("depth (mm)")+ylab("seasonal amplitud starch %")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(4)),
        axis.text.y=element_text(size = rel(4)),
        axis.title.x = element_text(size = rel(4)),
        axis.title.y = element_text(size = rel(4)),
        strip.text = element_text(size = rel(4)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

png("relative_sorage_change_Dac.png", width = 900, height = 600)
ggplot(data=Dac_woodStarch_radial_Activity_by_species, aes(x=depth_mm, y=rel_mean_active_mean))+
  geom_line(data=Dac_woodStarch_radial_activity, aes(x=depth_mm, y=rel_activ_mean, factor=ID), col="gray95")+
  geom_line( col="orangered2")+
  geom_ribbon( aes(ymin=rel_mean_active_mean-SD_rel_mean_active_mean, ymax=rel_mean_active_mean+SD_rel_mean_active_mean,
                   fill= "var"),
               fill="orangered2",
               alpha=0.3)+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  xlab("depth (mm)")+ylab("relative seasonal amplitud (%)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(4)),
        axis.text.y=element_text(size = rel(4)),
        axis.title.x = element_text(size = rel(4)),
        axis.title.y = element_text(size = rel(4)),
        strip.text = element_text(size = rel(4)),
        legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

## mean content of starch in the entire core 
## wood density 0.81 g/cm3 was taken from the recors from tanguro

Dac_woodStarch_means_by_tree$segment_volume_cm3=pi*(0.5/2)^2*0.5 ### this in in cm3

Dac_woodStarch_means_by_tree$grams_wood=Dac_woodStarch_means_by_tree$segment_volume_cm3*0.51
Dac_woodStarch_means_by_tree$grams_starch=Dac_woodStarch_means_by_tree$grams_wood*Dac_woodStarch_means_by_tree$mean_starch_area/100
Dac_woodStarch_means_by_tree$mg_starch=Dac_woodStarch_means_by_tree$grams_starch*1000

Dac_woodStarch_means_entire_tree=Dac_woodStarch_means_by_tree%>%group_by(month, ID)%>%
  summarise(mean_starch_percentage=mean(mean_starch_area),
            median_starch_percentage=mean(median_starch_are),
            mean_IQR=mean(Int_quart_range),
            starch_woodcore_content_by_mean=sum(mean_starch_area),
            starch_woodcore_content_by_median=sum(median_starch_are),
            starch_woodcore_content_grams=sum(grams_starch),
            starch_woodcore_content_mg=sum(mg_starch))


Dac_woodStarch_means_entire_tree=Dac_woodStarch_means_by_tree%>%group_by(month, ID)%>%
  summarise(mean_starch_percentage=mean(mean_starch_area),
            median_starch_percentage=median(median_starch_are),
            mean_IQR=mean(Int_quart_range),
            starch_woodcore_content_by_mean=sum(mean_starch_area),
            starch_woodcore_content_by_median=sum(median_starch_are),
            starch_woodcore_content_grams=sum(grams_starch))

Dac_woodStarch_means_entire_tree$starch_woodcore_content_mg=Dac_woodStarch_means_entire_tree$starch_woodcore_content_grams*1000



Dac_woodStarch_means_entire_tree_by_specie=Dac_woodStarch_means_entire_tree%>%group_by(month)%>%
  summarise(mean_starch_percentage=mean(mean_starch_percentage),
            sd_mean_starch_percentage=sd(mean_starch_percentage),
            median_starch_percentage=mean(median_starch_percentage),
            sd_median_starch_percentage=sd(median_starch_percentage),
            mean_IQR=mean(mean_IQR), 
            mean_starch_wood_content_by_mean=mean(starch_woodcore_content_by_mean),
            sd_starch_wood_content_by_mean=sd(starch_woodcore_content_by_mean),
            mean_starch_wood_content_by_median=mean(starch_woodcore_content_by_median),
            sd_starch_wood_content_by_medial=sd(starch_woodcore_content_by_median),
            mean_starch_wood_content_grams=mean(starch_woodcore_content_grams),
            sd_mean_starch_wood_content_grams=sd(starch_woodcore_content_grams),
            mean_starch_wood_content_mg=mean(starch_woodcore_content_mg)
  )

png("Dac_starch_wood_content_comp_month.png", width = 800, height = 600)
ggplot(Dac_woodStarch_means_entire_tree, aes(x=month, y=starch_woodcore_content_mg))+geom_boxplot(fill="orangered2",  alpha=0.3)+
  geom_point(aes(col=ID))+
  geom_point(data=Dac_woodStarch_means_entire_tree_by_specie, aes(x=month, y=mean_starch_wood_content_mg),  col="black", pch=17, size=3)+
  xlab("")+
  ylab("Starch woodcore content (mg)")+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


### pairwised comparison with paired samples for the content of sugar in the entire tree. The differences between months are lost mainly because the influence of the lack 
## of seasonality deeper in the stemwood or the contraty seasonality, and then all the averaging and the integration. 

Dac_pairwise_comparison_entire_tree=Dac_woodStarch_means_entire_tree
Dac_pairwise_comparison_entire_tree$month2=as.factor(as.numeric(Dac_pairwise_comparison_entire_tree$month))
Dac_pairwise_comparison_entire_tree=Dac_pairwise_comparison_entire_tree[-c(which(Dac_pairwise_comparison_entire_tree$ID=="Dac32001")),]
names(Dac_pairwise_comparison_entire_tree)
Dac_pairwise_comparison_entire_tree=Dac_pairwise_comparison_entire_tree[ ,c(2,10,9)]
Dac_pairwise_comparison_entire_tree$ID=as.factor(Dac_pairwise_comparison_entire_tree$ID)
str(Dac_pairwise_comparison_entire_tree)


Dac_pwc_starch_mg <- Dac_pairwise_comparison_entire_tree %>%
  pairwise_t_test(starch_woodcore_content_mg ~ month2, paired = T,
                  p.adjust.method = "BH"
  )
Dac_pwc_starch_mg


wilcox_test(starch_woodcore_content_mg~month2, data=Dac_pairwise_comparison_entire_tree, paired = T)


#### add growth data to Dac_woodStarch_content_woodcore and compare with total starch, mean starch, starch in the first two cm... 

monthly_growth2_summarized=monthly_growth2%>%filter(date3%in%c(8:20))%>%group_by(ID)%>%
  summarise(annual_growth_2019=sum(rate_growth),
            mean_monthly_growth_2019=mean(rate_growth),
            max_monthly_growth_2019=max(rate_growth),
            min_monthly_growth_2019=min(rate_growth)
  )


Dac_woodStarch_entire_tree=merge(Dac_woodStarch_means_entire_tree, monthly_growth2, by=c("ID", "month"))
Dac_woodStarch_entire_tree=merge(Dac_woodStarch_entire_tree, monthly_growth2_summarized, by=c("ID"))


Dac_woodStarch_entire_tree$significant_rel=as.character(Dac_woodStarch_entire_tree$month)
key_months_1=unique(as.character(Dac_woodStarch_entire_tree$month))
ID_key_signif2=c("Not_significant",
                 "Not_significant",
                 "Not_significant",
                 "Not_significant")

library(plyr)

Dac_woodStarch_entire_tree$significant_rel=mapvalues(Dac_woodStarch_entire_tree$significant_rel,
                                                     key_months_1,
                                                     ID_key_signif2)
detach("package:plyr", unload=TRUE)

png("annual_growth_vs_starch_content_Dac.png",
    width = 800, height = 600)
ggplot(Dac_woodStarch_entire_tree, aes(x=starch_woodcore_content_mg, y=annual_growth_2019))+
  geom_point(color="orangered2")+
  geom_smooth(method = "lm", se=F, color="orangered2", linetype="dashed")+
  stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE, size=8, xpos = 1,
                             ypos = 2.4, xpos2 = 1,
                             ypos2= 2.1)+
  xlab("starch content (mg)")+ ylab("annual growth 2019 (cm)")+
  ylim(0, 2.5)+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  facet_wrap(~month)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(2)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


# changes in starch content between the months. 

Dac_woodStarch_entire_tree=Dac_woodStarch_entire_tree[order(Dac_woodStarch_entire_tree$ID, Dac_woodStarch_entire_tree$date3),]
Dac_woodStarch_entire_tree$starch_change=NA
Dac_woodStarch_entire_tree$starch_change_percent=NA
Dac_woodStarch_entire_tree$rel_starch_change=NA

Dac_woodStarch_entire_tree2=Dac_woodStarch_entire_tree%>%group_by(ID)%>%
  mapGroups(function(x){
    for(i in 2:nrow(x)){
      x$starch_change[i]=x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1]
      x$starch_change_percent[i]=x$mean_starch_percentage[i]-x$mean_starch_percentage[i-1]
      x$rel_starch_change[i]=(x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1])/x$starch_woodcore_content_mg[i]
    }
    return(x)
  }
  )%>%ungroup()


## pairwise comparison between changes in the starch content between seasons 

png("Dac_starch_mass_change_comp_month.png", width = 800, height = 600)
ggplot(Dac_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=starch_change))+geom_boxplot(fill="orangered2",  alpha=0.3)+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (mg)")+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


png("Dac_starch_mass_change_rel_comp_month.png", width = 800, height = 600)
ggplot(Dac_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=rel_starch_change))+geom_boxplot(fill="orangered2",  alpha=0.3)+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (%)")+
  #ylim(-4, 2)+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()

png("Dac_starch_mass_change_percentage_month.png", width = 800, height = 600)
ggplot(Dac_woodStarch_entire_tree2%>%filter(month!="May19"), aes(x=month, y=starch_change_percent))+geom_boxplot(fill="orangered2",  alpha=0.3)+
  geom_point(aes(col=ID))+
  geom_hline(yintercept=0, linetype="dashed", color="red")+
  xlab("")+
  ylab("Starch mass change (%)")+
  #ylim(-4, 2)+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(1.5)),
        legend.position="none",
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),#legend.position="none",
        plot.title = element_text(size = rel(3), face = "italic"))
dev.off()


Dac_woodStarch_entire_tree2_wilcox_pair_comp=filter(Dac_woodStarch_entire_tree2, month!="May19")
Dac_woodStarch_entire_tree2_wilcox_pair_comp$month=factor(Dac_woodStarch_entire_tree2_wilcox_pair_comp$month,
                                                          levels = c("Aug19", "Nov19", "Feb20"))
Dac_woodStarch_entire_tree2_wilcox_pair_comp$month2=as.factor(as.numeric(Dac_woodStarch_entire_tree2_wilcox_pair_comp$month))
names(Dac_woodStarch_entire_tree2_wilcox_pair_comp)
Dac_woodStarch_entire_tree3_wilcox_pair_comp=Dac_woodStarch_entire_tree2_wilcox_pair_comp[,c("ID", "month2", "starch_change")]
Dac_woodStarch_entire_tree4_wilcox_pair_comp=Dac_woodStarch_entire_tree2_wilcox_pair_comp[,c("ID", "month2", "rel_starch_change")]
Dac_woodStarch_entire_tree3_wilcox_pair_comp$ID=as.factor(Dac_woodStarch_entire_tree3_wilcox_pair_comp$ID)
Dac_woodStarch_entire_tree4_wilcox_pair_comp$ID=as.factor(Dac_woodStarch_entire_tree4_wilcox_pair_comp$ID)
str(Dac_woodStarch_entire_tree3_wilcox_pair_comp)
str(Dac_woodStarch_entire_tree4_wilcox_pair_comp)

wilcox_test(starch_change~month2, data=Dac_woodStarch_entire_tree3_wilcox_pair_comp, paired = T)

wilcox_test(rel_starch_change~month2, data=Dac_woodStarch_entire_tree4_wilcox_pair_comp, paired = T)


### test for the relationship growth vs starch changes

Dac_woodStarch_entire_tree2=filter(Dac_woodStarch_entire_tree2,month!="May19")
Dac_woodStarch_entire_tree2$month_changes=as.character(Dac_woodStarch_entire_tree2$month)
month_chages_key=unique(as.character(Dac_woodStarch_entire_tree2$month))
ID_key=c("starch (May19-Aug19)\ngrowth (May19-Aug19)",
         "starch (Aug19-Nov19)\ngrowth (Aug19-Nov20)", 
         "starch (Nov19-Feb20)\n growth (Nov19-Feb20)"
         )

Dac_woodStarch_entire_tree2$significant_rel=as.character(Dac_woodStarch_entire_tree2$month)
ID_key_dac_signif=c("Not_significant",
                    "Not_significant",
                    "significant")

library(plyr)
Dac_woodStarch_entire_tree2$month_changes=mapvalues(Dac_woodStarch_entire_tree2$month_changes,
                                                    month_chages_key,
                                                    ID_key)

Dac_woodStarch_entire_tree2$significant_rel=mapvalues(Dac_woodStarch_entire_tree2$significant_rel,
                                                    month_chages_key,
                                                    ID_key_dac_signif)
detach("package:plyr", unload=TRUE)

Dac_woodStarch_entire_tree2$month_changes=factor(Dac_woodStarch_entire_tree2$month_changes,
                                                 levels = c("starch (May19-Aug19)\ngrowth (May19-Aug19)",
                                                            "starch (Aug19-Nov19)\ngrowth (Aug19-Nov20)", 
                                                            "starch (Nov19-Feb20)\n growth (Nov19-Feb20)"
                                                            ))

png("three_month_growth_after_vs_starch_change_Dac.png",
    width = 800, height = 600)
ggplot(Dac_woodStarch_entire_tree2, aes(x=starch_change, y=rate_growth_summ_three_months_prior))+
  geom_point(color="orangered2", size=4, alpha=0.5)+
  geom_smooth(aes(linetype=significant_rel), method = "lm", se=F, color="orangered")+
  scale_linetype_manual(values=c("dashed", "solid")) +
  stat_smooth_func_with_pval(geom="text",method="lm",hjust=0,parse=TRUE, size=7, xpos = -25,
                             ypos = 0.6, xpos2 = -25,
                             ypos2= 0.57)+
  ylim(-0.1,0.60)+
  xlab("starch change (mg)")+ylab("accumulated growth (cm)")+
  ggtitle("D. microcarpa\n(semi-deciduous/fiber-storing-species)")+
  geom_vline(xintercept = 0, linetype="dashed")+
  facet_wrap(~month_changes)+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(3)),
        legend.key = element_rect(size=rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size=rel(1)),
        legend.position="none",
        plot.title = element_text(size = rel(2.5), face = "italic"),
        strip.text.x = element_text(size= 20, margin = margin(10), vjust = 1),
        plot.margin = margin(10, 30, 10, 10,))
dev.off()

str(Dac_woodStarch_entire_tree2)
Dac_lm_starch_change_growth=Dac_woodStarch_entire_tree2%>%group_by(month_changes)%>%nest()%>%
mutate(lm_dac_starch_growth=map(data, ~lm(starch_change~rate_growth_summ_three_months_prior, data = .x)))
Dac_lm_starch_change_growth$lm_dac_starch_growth[[1]]
summary(Dac_lm_starch_change_growth$lm_dac_starch_growth[[3]])
plot(Dac_lm_starch_change_growth$lm_dac_starch_growth[[3]])

Dac_lm_starch_change_growth2=lm(starch_change~rate_growth_summ_three_months_prior*month_changes, data = Dac_woodStarch_entire_tree2)
summary(Dac_lm_starch_change_growth2)  
plot(Dac_lm_starch_change_growth2)
(Dac_lm_starch_change_growth2)
#########################################################################################################################################################


### joining the species together 


woodstarch_all_species=rbind(Ab_woodStarch_entire_tree,
                             Sac_woodStarch_entire_tree,
                             Dac_woodStarch_entire_tree
                            )
woodstarch_all_species$date=woodstarch_all_species$month

date_vector=c("2019-05-01",
              "2019-08-01",
              "2019-11-01",
              "2020-02-01")
month_vector=unique(woodstarch_all_species$month)

library(plyr)
woodstarch_all_species$date=mapvalues(woodstarch_all_species$date,
                                      month_vector,
                                      date_vector)

woodstarch_all_species$date=as.POSIXct(strptime(as.character(woodstarch_all_species$date),
                                                format = "%Y-%m-%d"))

woodstarch_all_species$species=gsub("[[:digit:]]","", woodstarch_all_species$ID)

species_tag=unique(woodstarch_all_species$species)
species_name=c("O.leucoxylon",
               "S. guianensis",
               "D. microcarpa"
               )

woodstarch_all_species$species=mapvalues(woodstarch_all_species$species,
                                         species_tag,
                                         species_name)

woodstarch_all_species$species_trait=woodstarch_all_species$species
species_tag=unique(woodstarch_all_species$species)
species_name2=c("O. leucoxylon\n(evergreen/parenchyma-storing-species)",
               "S. guianensis\n(semi-deciduous/parenchyma-storing-species)",
               "D. microcarpa\n(semi-deciduous/fiber-storing-species)"
)

woodstarch_all_species$species_trait=mapvalues(woodstarch_all_species$species_trait,
                                         species_tag,
                                         species_name2)


detach("package:plyr", unload=TRUE)

woodstarch_all_species$storage_strategy=ifelse(woodstarch_all_species$species%in%c("D. microcarpa", "T. burserifolia", "T. glaziovii", "T. guianensis"), "Fiber storage", "Parenchyma storage")
woodstarch_all_species$storage_lipids=ifelse(woodstarch_all_species$species%in%c("D. microcarpa", "T. burserifolia", "T. glaziovii", "V. vismiifolia"), "Lipid storage", "No lipid storage")
woodstarch_all_species$leaf_habit=ifelse(woodstarch_all_species$species%in%c("D. microcarpa", "V. vismiifolia", "S. guianensis", "T. guianensis"), "semi-deciduous", "evergreen")


load("climatic_Data_2018_2020.RData")

woodstarch_all_species_mean_species=woodstarch_all_species%>%group_by(month, species)%>%
  summarise(mean_starch_percentage2=mean(mean_starch_percentage),
            sd_mean_starch=sd(mean_starch_percentage),
            sterr_mean_starch =var(mean_starch_percentage)/sqrt(length(mean_starch_percentage)),
            median_starch_percentage=mean(median_starch_percentage),
            mean_IQR=mean(mean_IQR),
            starch_woodcore_content_by_mean2=mean(starch_woodcore_content_by_mean),
            sd_starch_woodcore_content_by_mean=sd(starch_woodcore_content_by_mean),
            starch_woodcore_content_by_median=mean(starch_woodcore_content_by_median),
            mean_starch_content_mg=mean(starch_woodcore_content_mg),
            sd_starch_content_mg=sd(starch_woodcore_content_mg),
            date=first(.data$date),
            storage_strategy=first(.data$storage_strategy,),
            storage_lipids=first(.data$storage_lipids),
           
  )

woodstarch_all_species_mean_species$cdl=c("a", "a","a",
                                          "b", "a", "a",
                                          "a", "a", "ab",
                                          "ab", "a", "b")

str(woodstarch_all_species_mean_species)
names(woodstarch_all_species)


sp_trait_key=unique(woodstarch_all_species$species_trait)
woodstarch_all_species_mean_species$species_trait=woodstarch_all_species_mean_species$species
species_key=c("O.leucoxylon", "S. guianensis", "D. microcarpa")

library(plyr)
woodstarch_all_species_mean_species$species_trait=mapvalues(woodstarch_all_species_mean_species$species_trait,
                                               species_key,
                                               sp_trait_key)  
detach("package:plyr", unload=TRUE)

png("total_starch_seasonality_D_S_O_separated_mg_boxplot.png",width = 2000, height = 800)
ggplot(woodstarch_all_species)+
  geom_area(data=climatic_data_months_2018_2020%>%filter(Year=="2019"| Year=="2020"),
            aes(x=date, y=precip_month/3), fill="grey80")+
  geom_boxplot(aes(x=date, y=starch_woodcore_content_mg, fill=storage_strategy, factor=as.factor(date)), )+
  geom_point(aes(x=date, y=starch_woodcore_content_mg))+
  geom_text(data=woodstarch_all_species_mean_species, aes(label=cdl, x=date, y= mean_starch_content_mg+(sd_starch_content_mg+20)), size=8)+
  #geom_errorbar(aes(x=date, ymin=mean_starch_content_mg-sd_starch_content_mg,
  #                 ymax=mean_starch_content_mg+sd_starch_content_mg, factor=species),
  #             position=position_dodge(), width=3000000)+
  scale_x_time(name="",
               breaks = c(unique(c(climatic_data_months_2018_2020%>%filter(Year=="2019"| Year=="2020"))$date))[c(3,4,9,10,15,16, 21, 22)],
               labels = c("Feb19", "Feb20", "May19", "May20", "Aug19", "Aug20", "Nov19", "Nov20"),
               limits = as.POSIXct(strptime(c("2019-01-01", "2020-06-01"), format = "%Y-%m-%d")))+
  #ylab("Starch content (g)")+
  scale_y_continuous(name="Starch mass (mg)",
                     sec.axis=sec_axis(~.*3, name="Rainfall"), limits = c(0,150))+
  #scale_fill_manual("Species")+ #values=c("Brown","darkolivegreen3", "aquamarine4"))+
  #ylim(0,250)+
  facet_wrap(~species_trait, scales = "free")+
  theme_classic()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(2.4), face = "italic"),
        legend.position = "none",
        legend.title = element_text(size=rel(3)),
        legend.text = element_text(size = rel(3))
  )
dev.off()


aov_storage_strategy=aov(starch_woodcore_content_mg~storage_strategy, data=
                    woodstarch_all_species)
summary(aov_storage_strategy)
TukeyHSD(aov_storage_strategy)
ggplot(woodstarch_all_species)+geom_boxplot(aes(x=storage_strategy, y=starch_woodcore_content_mg))

aov_storage_strategy=aov(mean_starch_percentage~storage_strategy, data=
                           woodstarch_all_species)
summary(aov_storage_strategy)
TukeyHSD(aov_storage_strategy)
ggplot(woodstarch_all_species)+geom_boxplot(aes(x=storage_strategy, y=mean_starch_percentage))

aov_leaf_habit=aov(starch_woodcore_content_mg~leaf_habit, data=
                           woodstarch_all_species)
summary(aov_leaf_habit)
TukeyHSD(aov_leaf_habit)
ggplot(woodstarch_all_species)+geom_boxplot(aes(x=leaf_habit, y=starch_woodcore_content_mg))


lm_allsp=aov(starch_woodcore_content_by_mean~storage_strategy*leaf_habit, data=
               woodstarch_all_species)
summary(lm_allsp)
TukeyHSD(lm_allsp)

###use this test to argue that storage strategy in combination with leaf habit has a effect on the amount of starch accumulated. 
lm_allsp=aov(starch_woodcore_content_mg~storage_strategy*leaf_habit, data=
               woodstarch_all_species)
summary(lm_allsp)
TukeyHSD(lm_allsp)
ggplot(woodstarch_all_species)+geom_boxplot(aes(x=storage_strategy, y=starch_woodcore_content_mg))+
  facet_wrap(~leaf_habit)


##### comparison between the changes in starch content between all the months examined. 
woodstarch_all_species2=woodstarch_all_species
woodstarch_all_species2=woodstarch_all_species2%>%filter(ID!="ab26820" & ID!="Sac1054" & ID!="Dac32001")

woodstarch_all_species2$starch_change=NA
woodstarch_all_species2$rel_starch_change=NA
woodstarch_all_species2$changed_months=as.character("NA")

woodstarch_all_species2=woodstarch_all_species2%>%group_by(ID)%>%
  mapGroups(function(x){
    for(i in 2:nrow(x)){
      x$starch_change[i]=x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1]
      x$rel_starch_change[i]=(x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-1])/x$starch_woodcore_content_mg[i]
      x$changed_months[i]=paste(x$month[i-1],"-",x$month[i], sep = "")
       }
    return(x)
  }
  )%>%ungroup()


woodstarch_all_species2$starch_change2=NA
woodstarch_all_species2$rel_starch_change2=NA
woodstarch_all_species2$changed_months2=as.character("NA")

woodstarch_all_species2=woodstarch_all_species2%>%group_by(ID)%>%
  mapGroups(function(x){
    for(i in 3:nrow(x)){
      x$starch_change2[i]=x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-2]
      x$rel_starch_change2[i]=(x$starch_woodcore_content_mg[i]-x$starch_woodcore_content_mg[i-2])/x$starch_woodcore_content_mg[i]
      x$changed_months2[i]=paste(x$month[i-2],"-", x$month[i], sep = "") 
       }
    return(x)
  }
  )%>%ungroup()


woodstarch_all_species2$starch_change3=NA
woodstarch_all_species2$rel_starch_change3=NA
woodstarch_all_species2$changed_months3=as.character("NA")

woodstarch_all_species2=woodstarch_all_species2%>%group_by(ID)%>%
  mapGroups(function(x){
      x$starch_change3[4]=x$starch_woodcore_content_mg[4]-x$starch_woodcore_content_mg[4-3]
      x$rel_starch_change3[4]=(x$starch_woodcore_content_mg[4]-x$starch_woodcore_content_mg[4-3])/x$starch_woodcore_content_mg[4]
      x$changed_months3[4]=paste(x$month[4-3],"-", x$month[4], sep = "")     
      return(x)
  }
  )%>%ungroup()

#### melt the columns together

names(woodstarch_all_species2)
id_vars=c("ID", "month", "starch_woodcore_content_mg", "Scientific_name", "species_trait", "storage_strategy", "storage_lipids", "leaf_habit")

woodstarch_all_species2_melt=woodstarch_all_species2%>%melt(id.vars = id_vars, measure.vars =c("starch_change","starch_change2", "starch_change3"),
                                                              value.name = "starch_change" 
                                                              )
  woodstarch_all_species2_melt2=woodstarch_all_species2%>%melt(id.vars = id_vars, measure.vars =c("rel_starch_change","rel_starch_change2", "rel_starch_change3"), 
                                                                    value.name = "rel_starch_change")
  woodstarch_all_species2_melt3=woodstarch_all_species2%>%melt(id.vars = id_vars, measure.vars =c("changed_months","changed_months2", "changed_months3"), 
                                                                    value.name = "changed_months")

  woodstarch_all_species2_melt=cbind(woodstarch_all_species2_melt,rel_starch_change=woodstarch_all_species2_melt2$rel_starch_change, changed_months=woodstarch_all_species2_melt3$changed_months)
  
  woodstarch_all_species2_melt=woodstarch_all_species2_melt[!is.na(woodstarch_all_species2_melt$starch_change),]
  
  unique(woodstarch_all_species2_melt$changed_months)
  woodstarch_all_species2_melt$changed_months=factor( woodstarch_all_species2_melt$changed_months, 
                                                      levels = c("May19-Aug19",
                                                                 "Aug19-Nov19",
                                                                 "Nov19-Feb20",
                                                                 "May19-Nov19",
                                                                 "Aug19-Feb20",
                                                                 "May19-Feb20"))
  
  aov_woodstarch_allsp_melt=aov(formula = rel_starch_change~changed_months*Scientific_name, data = woodstarch_all_species2_melt)
  summary(aov_woodstarch_allsp_melt)
  TukeyHSD(aov_woodstarch_allsp_melt)
  
  str(woodstarch_all_species2_melt)
  
  table(woodstarch_all_species2_melt$ID)
  
  woodstarch_all_species2_melt$wilcox_factor=paste(woodstarch_all_species2_melt$Scientific_name, woodstarch_all_species2_melt$changed_months, sep= "_")
  woodstarch_all_species2_melt$wilcox_factor=as.factor(woodstarch_all_species2_melt$wilcox_factor)
  
  wilcox_starch_change=woodstarch_all_species2_melt%>%group_by(Scientific_name)%>%
  wilcox_test(rel_starch_change~changed_months, data=.)
  
  wilcox_1=wilcox_test(rel_starch_change~wilcox_factor, data=woodstarch_all_species2_melt)
 
  library(boot) 
  
  # function to obtain the mean
  Bmean <- function(data, indices) {
    d <- data[indices] # allows boot to select sample 
    return(mean(d))
  } 
  
  woodstarch_all_species2_melt_boot_ci=woodstarch_all_species2_melt%>%filter(rel_starch_change>-2)%>%group_by(changed_months, Scientific_name)%>%nest()%>%
    mutate(boot_results=map(data, ~ boot(data=.x$rel_starch_change, statistic=Bmean, R=1000)),
           boostrap_var=map(boot_results, ~var(.x$t)),
           vector_bootstrap=unlist(boostrap_var),
           ci_boot=map(boot_results, boot.ci, var.t0=vector_bootstrap))

   names(woodstarch_all_species2_melt_boot_ci$ci_boot)=c(paste(woodstarch_all_species2_melt_boot_ci$Scientific_name,woodstarch_all_species2_melt_boot_ci$changed_months))
  
  mean_groupe=woodstarch_all_species2_melt%>%group_by(changed_months, Scientific_name)%>%
    summarise(mean_groups=mean(rel_starch_change))
  
  boot_results_unnest=woodstarch_all_species2_melt_boot_ci$boot_results
  boot_results_unnest_df=lapply(seq_along(boot_results_unnest),
                                function(i){
                                  unlist(boot_results_unnest[[i]][1])
                                })
  boot_results_unnest_df=as.data.frame(do.call(rbind, boot_results_unnest_df))

  ci_boot_unnest=woodstarch_all_species2_melt_boot_ci$ci_boot
  ci_boot_unnest_df= lapply(seq_along(ci_boot_unnest), function(i){
    unlist(ci_boot_unnest[[i]][7])
  })
  ci_boot_unnest_df=do.call("rbind", ci_boot_unnest_df)
  ci_boot_unnest_df=as.data.frame(ci_boot_unnest_df)
  ci_boot_unnest_df=cbind(ci_boot_unnest_df, sp=woodstarch_all_species2_melt_boot_ci$Scientific_name, month_changes=woodstarch_all_species2_melt_boot_ci$changed_months, 
                          mean=boot_results_unnest_df)

  sp_names=c("D. microcarpa\n(semi-deciduous/fiber-storing-species)",
             "S. guianensis\n(semi-deciduous/parenchyma-storing-species)",
             "O. leucoxylon\n(evergreen/parenchyma-storing-species)"
            )
  names(sp_names)=unique(ci_boot_unnest_df$sp)
  
  ci_boot_unnest_df$storage_strategy=ifelse(ci_boot_unnest_df$sp=="Dacryodes microcarpa", "Fiber_storage", "parenchyma_storage")
  
png("CI_rel_starch_changes_between_months_D_S_O_separated_boxplot.png",width = 2600, height = 800)
ggplot(ci_boot_unnest_df, aes(x=month_changes, y=t0, col=storage_strategy))+
  geom_point(size=6)+
  geom_errorbar(aes(ymin=bca4, ymax=bca5), size=3)+
  #ylim(-3, 2)+
  xlab("")+ylab("Mean relative starch change (%)")+
  geom_hline(yintercept = 0, linetype="dashed", color="red")+
  facet_wrap(~sp, labeller = labeller(sp=sp_names))+
  theme_classic()+
  theme(axis.text.x=element_text(size = rel(2.4)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(3)),
        axis.title.y = element_text(size = rel(3)),
        strip.text = element_text(size = rel(3), face = "italic"),
        legend.position = "none",
        legend.title = element_text(size=rel(3)),
        legend.text = element_text(size = rel(3))
  )
dev.off()
  
 
png("rel_starch_changes_between_months_D_S_O_separated_boxplot.png",width = 1500, height = 2000)
ggplot(woodstarch_all_species2_melt)+
  # geom_area(data=climatic_data_months_2018_2020%>%filter(Year=="2019"| Year=="2020"),
  #           aes(x=date, y=precip_month/3), fill="grey80")+
  geom_boxplot(aes(x=changed_months, y=rel_starch_change, fill=storage_strategy, alpha=0.3))+
  geom_hline(yintercept = 0, linetype="dashed", color="red")+
  xlab("")+ylab("Relative starch change (%)")+
  #geom_point(aes(x=month, y=starch_change))+
  #geom_errorbar(aes(x=date, ymin=mean_starch_content_mg-sd_starch_content_mg,
  #                 ymax=mean_starch_content_mg+sd_starch_content_mg, factor=species),
  #             position=position_dodge(), width=3000000)+
  # scale_x_time(name="",
  #              breaks = c(unique(c(climatic_data_months_2018_2020%>%filter(Year=="2019"| Year=="2020"))$date))[c(3,4,9,10,15,16, 21, 22)],
  #              labels = c("Feb19", "Feb20", "May19", "May20", "Aug19", "Aug20", "Nov19", "Nov20"),
  #              limits = as.POSIXct(strptime(c("2019-01-01", "2020-06-01"), format = "%Y-%m-%d")))+
  # #ylab("Starch content (g)")+
  # scale_y_continuous(name="Starch Content (mg)",
  #                    sec.axis=sec_axis(~.*3, name="Rainfall"), limits = c(0,150))+
  # #scale_fill_manual("Species")+ #values=c("Brown","darkolivegreen3", "aquamarine4"))+
  ylim(-4,1)+
  facet_wrap(~species_trait, scales = "free", ncol=1)+
  theme_classic()+
  theme(axis.text.x=element_text(size = rel(3.5)),
        axis.text.y=element_text(size = rel(6)),
        axis.title.x = element_text(size = rel(6)),
        axis.title.y = element_text(size = rel(6)),
        strip.text = element_text(size = rel(5), face = "italic"),
        legend.position = "none",
        legend.title = element_text(size=rel(3)),
        legend.text = element_text(size = rel(3))
  )
dev.off()


###### seasonality in starch content using the seasonal amplitud (maximum amount of starch vs the minimum no matter what is the moth)
### this method shows the maximun change in starch content per individual along the year, but it do not consider the sincronicy between trees as for some trees the max and min starch content occurred at different months. 

woodstarch_all_species3=woodstarch_all_species
woodstarch_all_species3=woodstarch_all_species3%>%filter(ID!="ab26820" & ID!="Sac1054" & ID!="Dac32001")


seasonal_amplitud_woodstarch=woodstarch_all_species3%>%
  group_by(ID, species_trait)%>%
  summarise(seasonal_amplitud_concentrations=max(mean_starch_percentage)-min(mean_starch_percentage),
            rel_seasonal_amplitud_concentrations= (max(mean_starch_percentage)-min(mean_starch_percentage))/max(mean_starch_percentage),
            seasonal_amplitud_content= max(starch_woodcore_content_mg)-min(starch_woodcore_content_mg),
            rel_seasonal_amplitud_content= (max(starch_woodcore_content_mg)-min(starch_woodcore_content_mg))/max(starch_woodcore_content_mg),
            max_month=.data$month[max(which(.data$starch_woodcore_content_mg==max(.data$starch_woodcore_content_mg)))],
            min_month=.data$month[min(which(.data$starch_woodcore_content_mg==min(.data$starch_woodcore_content_mg)))],
            sp=first(.data$species),
            storage_strategy=first(.data$storage_strategy)
  )

ggplot(seasonal_amplitud_woodstarch, aes(x=sp, y=rel_seasonal_amplitud_content))+
  geom_boxplot(aes(fill=storage_strategy))+
  geom_point()+
  geom_text(aes(label=max_month))
  

###### seasonal amplitud at each depth point. ( this may be more related to the mobility of starch in the stem wood than wiht the seasonality of the starch content)

Ab_woodStarch_radial_activity$species="O.leucoxylon"
Ab_woodStarch_radial_activity$storage_strategy="parenchyma_storage"
Ab_woodStarch_radial_activity$leaf_habit="evergreen"
Ab_woodStarch_radial_activity$lipid_storage="No_lipids"

Sac_woodStarch_radial_activity$species="S.guianenesis"
Sac_woodStarch_radial_activity$storage_strategy="parenchyma_storage"
Sac_woodStarch_radial_activity$leaf_habit="semi_deciduous"
Sac_woodStarch_radial_activity$lipid_storage="No_lipids"

Dac_woodStarch_radial_activity$species="D.microcarpa"
Dac_woodStarch_radial_activity$storage_strategy="fiber_storage"
Dac_woodStarch_radial_activity$leaf_habit="semi_deciduous"
Dac_woodStarch_radial_activity$lipid_storage="lipids"


radial_activity_2019=rbind(Ab_woodStarch_radial_activity, 
                           Sac_woodStarch_radial_activity,
                           Dac_woodStarch_radial_activity)  


seasonal_woodstarch_amplitud_indepth=radial_activity_2019%>%
  group_by(ID, species)%>%
  summarise(mean_activ_mean=mean(activ_mean),
            sd_active_mean=sd(activ_mean),
            mean_rel_activ_mean=mean(rel_activ_mean),
            sd_mean_rel_active_mean=sd(rel_activ_mean),
            mean_activ_median=mean(activ_median)
  )

ggplot(seasonal_woodstarch_amplitud_indepth, aes(x=as.factor(species), y=mean_activ_mean))+geom_boxplot()
ggplot(seasonal_woodstarch_amplitud_indepth, aes(x=as.factor(species), y=mean_rel_activ_mean))+geom_boxplot()
ggplot(seasonal_woodstarch_amplitud_indepth, aes(x=as.factor(species), y=mean_activ_median))+geom_boxplot()


ggplot(seasonal_woodstarch_amplitud_indepth, aes(x=as.factor(species), y=mean_rel_activ_mean))+geom_boxplot()+
  geom_point()+
  geom_text(aes(label=ID))

aoc_seasonal_range2=aov(mean_rel_activ_mean~as.factor(species), data=seasonal_woodstarch_amplitud_indepth)
summary(aoc_seasonal_range2)
TukeyHSD(aoc_seasonal_range2)



##########################################################################################################################################
##########################################################################################################################################
##################################################################### #####################################################################
##########################################################################################################################################
##########################################################################################################################################

### climatic analyzsis 

climatic_data_months_2018_2020

climatic_data_2019=climatic_data_months_2018_2020[climatic_data_months_2018_2020$Year!=2018, ]
climatic_data_2019=climatic_data_2019[order(climatic_data_2019$date),]

climatic_data_2019$PP_three_month_average=NA
climatic_data_2019$PP_three_month_accumulated=NA
climatic_data_2019$RH_three_month_average=NA
for(i in 2:22){
  climatic_data_2019$PP_three_month_average[i]=mean(climatic_data_2019$precip_month[c(i-1,i,i+1)])
  climatic_data_2019$PP_three_month_accumulated[i]=sum(climatic_data_2019$precip_month[c(i-1,i,i+1)])
  climatic_data_2019$RH_three_month_average[i]=mean(climatic_data_2019$mean_RH_month[c(i-1,i,i+1)])
  
}

names(climatic_data_2019)

climatic_data_2019$PP_three_month_prior_average=NA
climatic_data_2019$PP_three_month_prior_accumulated=NA
climatic_data_2019$RH_three_month_prior_average=NA
for(i in 3:23){
  climatic_data_2019$PP_three_month_prior_average[i]=mean(climatic_data_2019$precip_month[c(i-2,i-1,i)])
  climatic_data_2019$PP_three_month_prior_accumulated[i]=sum(climatic_data_2019$precip_month[c(i-2,i-1,i)])
  climatic_data_2019$RH_three_month_prior_average[i]=mean(climatic_data_2019$mean_RH_month[c(i-2,i-1,i)])
  
}

climatic_data_2019$PP_three_month_after_average=NA
climatic_data_2019$PP_three_month_after_accumulated=NA
climatic_data_2019$RH_three_month__after_average=NA
for(i in 1:21){
  climatic_data_2019$PP_three_month_after_average[i]=mean(climatic_data_2019$precip_month[c(i,i+1,i+2)])
  climatic_data_2019$PP_three_month_after_accumulated[i]=sum(climatic_data_2019$precip_month[c(i,i+1,i+2)])
  climatic_data_2019$RH_three_month_after_average[i]=mean(climatic_data_2019$mean_RH_month[c(i,i+1,i+2)])
  
}


woodstarch_all_species_climatic_comp=merge(woodstarch_all_species, climatic_data_2019,by = "date")



library(broom)
lm_precip_starch_per_sp=woodstarch_all_species_climatic_comp%>%group_by(species)%>%nest()%>%
  mutate(lm1=map(data, ~ lm(mean_starch_percentage~precip_month+I(precip_month^2), data=.x)),
         lm2=map(data, ~ lm(mean_starch_percentage~PP_three_month_prior_average+I(PP_three_month_prior_average^2), data=.x)),
         lm3=map(data, ~ lm(starch_woodcore_content_mg~PP_three_month_prior_average+I(PP_three_month_prior_average^2), data=.x)),
         lm4=map(data, ~ lm(starch_woodcore_content_mg~precip_month+I(precip_month^2), data=.x)),
         lm5=map(data, ~ lm(starch_woodcore_content_mg~PP_three_month_prior_accumulated+I(PP_three_month_prior_accumulated^2), data=.x)),
         lm6=map(data, ~ lm(starch_woodcore_content_mg~PP_three_month_accumulated+I(PP_three_month_accumulated^2), data=.x)),
         lm7=map(data, ~ lm(starch_woodcore_content_mg~precip_month, data=.x)),
         tidied1=map(lm1, tidy),
         tidied2=map(lm2, tidy),
         tidied3=map(lm3, tidy),
         tidied4=map(lm4, tidy),
         tidied5=map(lm5, tidy),
         tidied6=map(lm6, tidy),
         tidied7=map(lm7, tidy))
names(lm_precip_starch_per_sp$tidied1)=lm_precip_starch_per_sp$species
names(lm_precip_starch_per_sp$tidied2)=lm_precip_starch_per_sp$species
names(lm_precip_starch_per_sp$tidied3)=lm_precip_starch_per_sp$species
names(lm_precip_starch_per_sp$tidied4)=lm_precip_starch_per_sp$species
names(lm_precip_starch_per_sp$tidied5)=lm_precip_starch_per_sp$species
names(lm_precip_starch_per_sp$tidied6)=lm_precip_starch_per_sp$species
names(lm_precip_starch_per_sp$tidied7)=lm_precip_starch_per_sp$species
lm_precip_starch_per_sp$tidied1
lm_precip_starch_per_sp$tidied2
lm_precip_starch_per_sp$tidied3
lm_precip_starch_per_sp$tidied4
lm_precip_starch_per_sp$tidied5
lm_precip_starch_per_sp$tidied6
lm_precip_starch_per_sp$tidied7

woodstarch_all_species_climatic_comp$precip_month2=as.integer(as.character(woodstarch_all_species_climatic_comp$precip_month))

as.numeric(woodstarch_all_species_climatic_comp$precip_month)
as.factor(round(woodstarch_all_species_climatic_comp$precip_month,0))
as.integer(woodstarch_all_species_climatic_comp$precip_month)
str(woodstarch_all_species_climatic_comp)
str(monthly_means_growt$month_numeric)

woodstarch_all_species_climatic_comp$precip_month3=match(as.character(woodstarch_all_species_climatic_comp$precip_month2), c( "0", "49", "68", "108"))

woodstarch_all_species_climatic_comp$precip_month3[woodstarch_all_species_climatic_comp$precip_month3==1]=49
woodstarch_all_species_climatic_comp$precip_month3[woodstarch_all_species_climatic_comp$precip_month3==2]=0
woodstarch_all_species_climatic_comp$precip_month3[woodstarch_all_species_climatic_comp$precip_month3==3]=68
woodstarch_all_species_climatic_comp$precip_month3[woodstarch_all_species_climatic_comp$precip_month3==4]=108

woodstarch_all_species_climatic_comp$precip_month_factor=as.factor(woodstarch_all_species_climatic_comp$precip_month)
aov_precip_starch_wilcoxon=woodstarch_all_species_climatic_comp%>%group_by(species)%>%
  mapGroups(function(df){
    unique_ID=names(which(table(df$ID)>3))
    df=df%>%filter(ID%in%unique_ID)
    Wc_test=wilcox_test(starch_woodcore_content_mg~precip_month_factor, data=df, paired = T)
    return(Wc_test)
  })%>%ungroup()

woodstarch_parenchyma_storage_species=woodstarch_all_species_climatic_comp%>%filter(storage_strategy=="Parenchyma storage")

png("starch_content_and_monthly_precipitation_parenchyma_storing_sp.png", 
    width = 800, height = 600)
ggplot(woodstarch_parenchyma_storage_species, aes(x=precip_month, y=starch_woodcore_content_mg), color="Aquamarine")+
  geom_point(alpha=0.2, color="turquoise2")+
  geom_smooth(method="lm", formula=y~x, se=F, color="turquoise2")+
  xlab("Monthly precipitation (mm)")+ylab("starch (mg)")+
  facet_wrap(~species, scales = "free")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(2.5)),
        axis.title.y = element_text(size = rel(2.5)),
        strip.text = element_text(size = rel(3)),
        legend.text = element_text(size=rel(3)),
        legend.title = element_text(size = rel(3)),
        legend.key = element_rect(size=rel(3))
  )
dev.off()

woodstarch_Fiber_storage_species=woodstarch_all_species_climatic_comp%>%filter(storage_strategy=="Fiber storage")

png("starch_content_and_precipitation_fiber_storing_sp.png", 
    width = 400, height = 600)
ggplot(woodstarch_Fiber_storage_species, aes(x=precip_month, y=starch_woodcore_content_mg), color="orangered")+
  geom_point(alpha=0.2, color="orangered")+
  geom_smooth(method="lm", formula=y~x+I(x^2), se=F, color="orangered")+
  xlab("Monthly precipitation (mm)")+ylab("starch (mg)")+
  facet_wrap(~species, scales = "free")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(2.5)),
        axis.title.y = element_text(size = rel(2.5)),
        strip.text = element_text(size = rel(3)),
        legend.text = element_text(size=rel(3)),
        legend.title = element_text(size = rel(3)),
        legend.key = element_rect(size=rel(3))
  )
dev.off()




##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################
##################################################################################################################################################################################################################################

### RH


lm_RH_starch_per_sp=woodstarch_all_species_climatic_comp%>%group_by(species)%>%nest()%>%
  mutate(lm1=map(data, ~ lm(mean_starch_percentage~mean_RH_month+I(mean_RH_month^2), data=.x)),
         lm2=map(data, ~ lm(mean_starch_percentage~RH_three_month_average+I(RH_three_month_average^2), data=.x)),
         lm3=map(data, ~ lm(mean_starch_percentage~RH_three_month_prior_average+I(RH_three_month_prior_average^2), data=.x)),
         lm4=map(data, ~ lm(mean_starch_percentage~RH_three_month_after_average+I(RH_three_month_after_average^2), data=.x)),
         lm5=map(data, ~ lm(starch_woodcore_content_mg~mean_RH_month+I(mean_RH_month^2), data=.x)),
         lm6=map(data, ~ lm(starch_woodcore_content_mg~RH_three_month_average+I(RH_three_month_average^2), data=.x)),
         lm7=map(data, ~ lm(starch_woodcore_content_mg~RH_three_month_prior_average+I(RH_three_month_prior_average^2), data=.x)),
         lm8=map(data, ~ lm(starch_woodcore_content_mg~RH_three_month_after_average+I(RH_three_month_after_average^2), data=.x)),
         lm9=map(data, ~ lm(median_starch_percentage~mean_RH_month+I(mean_RH_month^2), data=.x)),
         
         tidied1=map(lm1, tidy),
         tidied2=map(lm2, tidy),
         tidied3=map(lm3, tidy),
         tidied4=map(lm4, tidy),
         tidied5=map(lm5, tidy),
         tidied6=map(lm6, tidy),
         tidied7=map(lm7, tidy),
         tidied8=map(lm8, tidy),
         tidied9=map(lm9, tidy))
names(lm_RH_starch_per_sp$tidied1)=lm_RH_starch_per_sp$species
names(lm_RH_starch_per_sp$tidied2)=lm_RH_starch_per_sp$species
names(lm_RH_starch_per_sp$tidied3)=lm_RH_starch_per_sp$species
names(lm_RH_starch_per_sp$tidied4)=lm_RH_starch_per_sp$species
names(lm_RH_starch_per_sp$tidied5)=lm_RH_starch_per_sp$species
names(lm_RH_starch_per_sp$tidied6)=lm_RH_starch_per_sp$species
names(lm_RH_starch_per_sp$tidied7)=lm_RH_starch_per_sp$species
names(lm_RH_starch_per_sp$tidied8)=lm_RH_starch_per_sp$species
names(lm_RH_starch_per_sp$tidied9)=lm_RH_starch_per_sp$species
lm_RH_starch_per_sp$tidied1
lm_RH_starch_per_sp$tidied2
lm_RH_starch_per_sp$tidied3
lm_RH_starch_per_sp$tidied4
lm_RH_starch_per_sp$tidied5
lm_RH_starch_per_sp$tidied6
lm_RH_starch_per_sp$tidied7
lm_RH_starch_per_sp$tidied8
lm_RH_starch_per_sp$tidied9


png("starch_content_and_RH_sp.png", 
    width = 1200, height = 600)
ggplot(woodstarch_all_species_climatic_comp, aes(x=RH_three_month_prior_average, y=starch_woodcore_content_mg, col=storage_strategy) )+
  geom_point(alpha=0.2)+
  geom_smooth(method="lm", formula=y~x+I(x^2), se=F)+
  xlab("Mean last three months RH (%)")+ylab("starch (mg)")+
  facet_wrap(~species, scales = "free")+
  theme_bw()+
  theme(axis.text.x=element_text(size = rel(3)),
        axis.text.y=element_text(size = rel(3)),
        axis.title.x = element_text(size = rel(2.5)),
        axis.title.y = element_text(size = rel(2.5)),
        strip.text = element_text(size = rel(3)),
        legend.text = element_text(size=rel(3)),
        legend.title = element_text(size = rel(3)),
        legend.position = "none",
        legend.key = element_rect(size=rel(3))
  )
dev.off()
