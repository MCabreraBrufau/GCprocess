#Quality control ECD vs FID detectors


#Author: Miguel Cabrera. 

#Description: This script uses the GC excel spreadsheet as input to calculate the calibration factor for each day of analysis as well as the stability of P5 standards. 

library(tidyverse)
library(readxl)
library(broom)
library(ggpmisc)
library(egg)


dropbox_root<- "C:/Users/Miguel/Dropbox/"
gcfolder<- paste0(dropbox_root,"RESTORE4Cs - Fieldwork/Data/GHG/GC data/")
plots_path<- paste0(dropbox_root,"GC data/plots/")


gcog<- read_xlsx(path = paste0(gcfolder,"DATA_R4Cs_formatted_ok_v20240615.xlsx"),range = cell_cols("A:I"))
names(gcog)
names(gcog)<- c("sampleid","seq_pos","date","batch","fileid","ch4","co2","n2o","obs")

gc_work<- gcog %>% 
  select(sampleid, seq_pos,date,batch,ch4,co2,n2o,obs) %>% 
  mutate(datef=as.Date.character(as.character(date),format="%Y%m%d")) %>% 
  filter(grepl('^P', sampleid)&datef>=as.Date("2024-03-19")) %>%  #keep only standards after introduction of P5 in each batch (remove all before 20240319)
  filter(!grepl("mix",sampleid)) %>% #remove Pmix standards
  mutate(sampleid=sub(",",".",sampleid),
         sampleid=sub("-","_",sampleid)) %>% 
  mutate(standardgas=sub("\\_.*", "", sampleid),#Extract PXX from sampleid
         vol= case_when(standardgas=="P20"~11.2,
                       TRUE~as.numeric(sub("P","", standardgas)))) %>% #Extract standard volume gas from PXX
  mutate(uniqueid=paste(rownames(.),standardgas,date,batch,sep = "_"))#add unique id with explicit info


#Check uniqueid is in fact unique
dim(gc_work)[1]==length(unique(gc_work$uniqueid))

#Check what dates have a complete cal curve integrated
dates_with_calcurve<- gc_work %>% 
  filter(batch==0)%>%
  group_by(date)%>%
  summarise(co2=mean(co2))%>%
  filter(!is.na(co2))%>%#drop dates with 1 or more NAs in co2 
  pull(date)#create vector

#Subset dates with complete cal curve
gc_complete<- gc_work %>% 
  filter(date%in%dates_with_calcurve) %>% 
  filter(date!=20240614&date!=20240615)





gc_complete%>%
  filter(date==20240410) %>% 
  filter(batch==0) %>% 
  ggplot(aes(x=vol, y=ch4))+
  geom_point()+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "n")))





#Calculate all the cal curves
cal<- gc_complete %>%
  filter(batch==0) %>% 
  # filter(vol!=11.2) %>% #NO P20s!
  # filter(vol!=0) %>%  #NO blanks!
  select(uniqueid,date,vol,co2,ch4,n2o)%>%
  pivot_longer(cols = c(co2,ch4,n2o), names_to = "gas",values_to = "area") %>% 
  nest_by(date,gas) %>% 
  mutate(mod = list(lm(area~vol, data=data))) %>%
  reframe(tidy(mod),glance(mod))
  
  
#Plot factors 
cal%>%
  filter(term=="vol")%>%
  # filter(r.squared>0.95) %>% 
ggplot(aes(x=as.Date(as.character(date),format="%Y%m%d"), y=estimate, col=gas))+
         geom_point()
  
##Plot R2
cal%>%
  filter(term=="vol")%>%
  ggplot(aes(x=as.Date(as.character(date),format="%Y%m%d"), y=r.squared, col=gas))+
  geom_point()

 
#We need to inspect all the calibration curves one by one!
 

# PLOT relative deviation of P5 areas for each of the gasses
gc_work %>% 
  filter(date!=20240614&date!=20240615) %>% #Exclude last days
  filter(vol==5) %>% 
  pivot_longer(cols = c(co2,ch4,n2o), names_to = "gas",values_to = "area") %>% 
  filter(!is.na(area)) %>% 
  group_by(gas) %>% 
  mutate(meanP5area=mean(area)) %>% 
  ggplot(aes(x=gas, y=area/meanP5area, col=gas))+geom_violin()+
  geom_jitter(width = 0.2,size = 1)



# ---- function to save a list of plots into pdf file ----

# gg_save_pdf = function(list, filename) {
#   pdf(filename,width = 14,height = 7)
#   for (p in list) {
#     print(p)
#   }
#   dev.off()
#   invisible(NULL)
# }
# 
# plt_list <- vector('list', length(unique(gc_complete$date)))
# 
# 
# 
# for(datecalcurve in unique(gc_complete$date)){
#   
#   calcurve<- gc_complete %>% filter(date==datecalcurve) %>% 
#     filter(batch==0) %>% 
#     filter(vol!=11.2) %>% 
#     filter(vol!=0)
# 
#   qual_inj<- gc_complete %>% 
#     filter(date==datecalcurve) %>% 
#     filter(!(uniqueid%in%calcurve$uniqueid))
#   
#   plt_CO2 <- ggplot(calcurve, aes(vol, co2))+
#     geom_point()+
#     stat_poly_line() +
#     stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
#     geom_point(data = qual_inj, aes(vol,co2, col="excluded"))+
#     theme_article()+
#     labs(col = "") +
#     theme(legend.position = "inside", legend.position.inside = c(0.8,0.1))+
#     xlab("Std gas vol (ml)")+
#     ylab("CO2 Area (FID)")+
#     ggtitle(paste0(datecalcurve," cal. plots"))
#             
#   plt_CH4 <- ggplot(calcurve, aes(vol, ch4))+
#     geom_point()+
#     stat_poly_line() +
#     stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
#     geom_point(data = qual_inj, aes(vol,ch4, col="excluded"))+
#   theme_article()+
#     labs(col = "") +
#     theme(legend.position = "inside", legend.position.inside = c(0.8,0.1))+
#     xlab("Std gas vol (ml)")+
#     ylab("CH4 Area (FID)")
#   
#   plt_N2O <- ggplot(calcurve, aes(vol, n2o))+
#     geom_point()+
#     stat_poly_line() +
#     stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
#     geom_point(data = qual_inj, aes(vol,n2o, col="excluded"))+
#   theme_article()+
#     labs(col = "") +
#     theme(legend.position = "inside", legend.position.inside = c(0.8,0.1))+
#     xlab("Std gas vol (ml)")+
#     ylab("N2O Area (ECD)")
#   
#   # plt <- ggarrange(plt_CO2, plt_CH4, plt_H2O, ncol = 1)
#   plt_list[[which(unique(gc_complete$date)==datecalcurve)]] <- ggarrange(plt_CO2, plt_CH4, plt_N2O, ncol = 3,nrow = 1)
#   
# }
# 
#   # Print pdf
#   setwd(plots_path)
#   gg_save_pdf(list = plt_list, filename = paste0("calcurves.pdf"))
# 
#   
#   
#   
#   
#   
#   gg_save_pdf = function(list, filename) {
#     pdf(filename,width = 14,height = 7)
#     for (p in list) {
#       print(p)
#     }
#     dev.off()
#     invisible(NULL)
#   }
#   
#   plt_list <- vector('list', length(unique(gc_complete$date)))
#   
rm(qual_inj,calcurve,plt_CH4,plt_CO2,plt_N2O,datecalcurve)
  gc_complete_long <- gc_complete %>% 
    pivot_longer(cols=c("co2","ch4","n2o"), names_to = "gas", values_to = "area") %>% 
    mutate(uniqueareaid=paste(date,batch,seq_pos,gas,sep="_"))
  
  #Outliers uniqueareaid:
  co2_out<- c("20240402_0_3_co2","20240404_0_3_co2","20240415_0_3_co2")
  ch4_out<- c("20240402_0_3_ch4","20240405_0_7_ch4","20240410_0_3_ch4","20240415_0_3_ch4")
  n2o_out<- c("20240321_0_5_n2o","20240404_0_5_n2o","20240506_0_3_n2o")
  
  plt_list <- vector('list', length(unique(gc_complete$date)))
  
  gg_save_pdf = function(list, filename) {
        pdf(filename,width = 14,height = 7)
        for (p in list) {
          print(p)
        }
        dev.off()
        invisible(NULL)
      }
  
  
  
  for(datecalcurve in unique(gc_complete$date)){
    
    calcurve<- gc_complete_long %>% filter(date==datecalcurve) %>% 
      filter(batch==0) %>% 
      filter(vol!=11.2) %>% 
      filter(vol!=0) %>% 
      filter(!(uniqueareaid)%in%c(co2_out,ch4_out,n2o_out))
    
    qual_inj<- gc_complete_long %>% 
      filter(date==datecalcurve) %>% 
      filter(!(uniqueareaid%in%calcurve$uniqueareaid))
    
    
    plt_CO2 <- ggplot(subset(calcurve, gas=="co2"), aes(vol, area))+
      geom_point()+
      geom_text(aes(label=uniqueareaid))+
      stat_poly_line() +
      stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
      geom_point(data = subset(qual_inj, gas=="co2"), aes(vol,area, col="excluded"))+
      geom_text(data=subset(qual_inj, gas=="co2"), aes(vol, area,label=uniqueareaid))+
      theme_article()+
      labs(col = "") +
      theme(legend.position = "inside", legend.position.inside = c(0.8,0.1))+
      xlab("Std gas vol (ml)")+
      ylab("CO2 Area (FID)")+
      ggtitle(paste0(datecalcurve," cal. plots"))
    
    plt_CH4 <- ggplot(subset(calcurve, gas=="ch4"), aes(vol, area))+
      geom_point()+
      geom_text(aes(label=uniqueareaid))+
      stat_poly_line() +
      stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
      geom_point(data = subset(qual_inj, gas=="ch4"), aes(vol,area, col="excluded"))+
      geom_text(data=subset(qual_inj, gas=="ch4"), aes(vol, area,label=uniqueareaid))+
      theme_article()+
      labs(col = "") +
      theme(legend.position = "inside", legend.position.inside = c(0.8,0.1))+
      xlab("Std gas vol (ml)")+
      ylab("CH4 Area (FID)")
    
    plt_N2O <- ggplot(subset(calcurve, gas=="n2o"), aes(vol, area))+
      geom_point()+
      geom_text(aes(label=uniqueareaid))+
      stat_poly_line() +
      stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
      geom_point(data = subset(qual_inj, gas=="n2o"), aes(vol,area, col="excluded"))+
      geom_text(data=subset(qual_inj, gas=="n2o"), aes(vol, area,label=uniqueareaid))+
      theme_article()+
      labs(col = "") +
      theme(legend.position = "inside", legend.position.inside = c(0.8,0.1))+
      xlab("Std gas vol (ml)")+
      ylab("N2O Area (ECD)")
    
    # plt <- ggarrange(plt_CO2, plt_CH4, plt_H2O, ncol = 1)
    plt_list[[which(unique(gc_complete$date)==datecalcurve)]] <- ggarrange(plt_CO2, plt_CH4, plt_N2O, ncol = 3,nrow = 1)
    
  }
  

  
  # Print pdf
  setwd(plots_path)
  gg_save_pdf(list = plt_list, filename = paste0("calcurves.pdf"))
  
  
  