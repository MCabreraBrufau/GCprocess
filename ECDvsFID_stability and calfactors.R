#Quality control ECD vs FID detectors


#Author: Miguel Cabrera. 

#Description: This script uses the GC excel spreadsheet as input to calculate the calibration factor for each day of analysis as well as the stability of P5 standards. 

library(tidyverse)
library(readxl)
library(broom)
library(ggpmisc)
library(egg)

#Directories for data and plots
dropbox_root<- "C:/Users/Miguel/Dropbox/"
gcfolder<- paste0(dropbox_root,"RESTORE4Cs - Fieldwork/Data/GHG/GC data/")
plots_path<- paste0(dropbox_root,"GC data/plots/")



# ---- function to save a list of plots into pdf file ----
gg_save_pdf = function(list, filename) {
  pdf(filename,width = 14,height = 7)
  for (p in list) {
    print(p)
  }
  dev.off()
  invisible(NULL)
}



#Data import and formatting

gcog<- read_xlsx(path = paste0(gcfolder,"DATA_R4Cs_formatted_ok_v20240618.xlsx"),range = cell_cols("A:I"))
names(gcog)<- c("sampleid","seq_pos","date","batch","fileid","ch4","co2","n2o","obs")


#Select only standard mix injections and formatting
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
  mutate(uniqueinj_id=paste(rownames(.),date,batch,seq_pos,sep="_"))#add unique id with explicit info


#Check uniqueid is in fact unique
dim(gc_work)[1]==length(unique(gc_work$uniqueinj_id))

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
  filter(date<20240614)#keep only dates before overnight ECD running test


#example plot of calcurve
gc_complete%>%
  filter(date==20240410) %>% 
  filter(batch==0) %>% 
  ggplot(aes(x=vol, y=ch4))+
  geom_point()+
  stat_poly_line() +
  stat_poly_eq(use_label(c("eq", "adj.R2", "n")))

#Pivot longer for calculations and add uniqueareaid
gc_complete_long <- gc_complete %>% 
  pivot_longer(cols=c("co2","ch4","n2o"), names_to = "gas", values_to = "area") %>% 
  mutate(uniqueareaid=paste(date,batch,seq_pos,gas,sep="_"))



#Outliers uniqueareaid: filled in manually after inspection of individual calplots 
co2_out<- c("20240402_0_3_co2","20240404_0_3_co2","20240415_0_3_co2")
ch4_out<- c("20240402_0_3_ch4","20240405_0_7_ch4","20240410_0_3_ch4","20240415_0_3_ch4")
n2o_out<- c("20240321_0_5_n2o","20240404_0_5_n2o","20240506_0_3_n2o")



#Calculate all the cal curves parameters
cal<- gc_complete_long %>%
  filter(batch==0) %>% 
  filter(vol!=11.2) %>% #NO P20s!
  filter(vol!=0) %>%  #NO blanks!
  filter(!(uniqueareaid)%in%c(co2_out,ch4_out,n2o_out)) %>% 
  select(uniqueinj_id,uniqueareaid,date,vol,gas,area)%>%
  nest_by(date,gas) %>% 
  mutate(mod = list(lm(area~vol, data=data))) %>%
  reframe(tidy(mod),glance(mod))

  
#Plot factors 
cal%>%
  filter(term=="vol")%>%
  # filter(r.squared>0.95) %>% 
ggplot(aes(x=as.Date(as.character(date),format="%Y%m%d"), y=estimate, col=gas))+
         geom_point()+
  facet_wrap(facets=.~gas, scales = "free")
  
##Plot R2
cal%>%
  filter(term=="vol")%>%
  ggplot(aes(x=as.Date(as.character(date),format="%Y%m%d"), y=r.squared, col=gas))+
  geom_point()

 
#We need to inspect all the calibration curves one by one! DONE, (outliers annotated in gas_out vectors)
 
#PLot every calcurve (after excluding outliers explictly, specified in "gas_out" vectors)

#Remove objects created in loop (for repeated executions of loop) 
rm(excluded,calcurve,plt_CH4,plt_CO2,plt_N2O,datecalcurve)

#create vector to store calcurve plots
plt_list <- vector('list', length(unique(gc_complete$date)))

#Loop over each date and plot the calcurve for all 3 ghgs
for(datecalcurve in unique(gc_complete$date)){
  
  calcurve<- gc_complete_long %>% 
    filter(date==datecalcurve) %>% 
    filter(batch==0) %>% 
    filter(vol!=11.2) %>% 
    filter(vol!=0) %>% 
    filter(!(uniqueareaid)%in%c(co2_out,ch4_out,n2o_out))
  
  excluded<- gc_complete_long %>% 
    filter(date==datecalcurve) %>% 
    filter(batch==0) %>% 
    filter(vol!=11.2) %>% 
    filter(vol!=0) %>% 
    filter(!(uniqueareaid%in%calcurve$uniqueareaid))
  
  
  plt_CO2 <- ggplot(subset(calcurve, gas=="co2"), aes(vol, area))+
    geom_point()+
    # geom_text(aes(label=uniqueareaid))+
    stat_poly_line() +
    stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
    geom_point(data = subset(excluded, gas=="co2"), aes(vol,area, col="excluded"))+
    # geom_text(data=subset(qual_inj, gas=="co2"), aes(vol, area,label=uniqueareaid))+
    theme_article()+
    labs(col = "") +
    theme(legend.position = "inside", legend.position.inside = c(0.8,0.1))+
    xlab("Std gas vol (ml)")+
    ylab("CO2 Area (FID)")+
    ggtitle(paste0(datecalcurve," cal. plots"))
  
  plt_CH4 <- ggplot(subset(calcurve, gas=="ch4"), aes(vol, area))+
    geom_point()+
    # geom_text(aes(label=uniqueareaid))+
    stat_poly_line() +
    stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
    geom_point(data = subset(excluded, gas=="ch4"), aes(vol,area, col="excluded"))+
    # geom_text(data=subset(qual_inj, gas=="ch4"), aes(vol, area,label=uniqueareaid))+
    theme_article()+
    labs(col = "") +
    theme(legend.position = "inside", legend.position.inside = c(0.8,0.1))+
    xlab("Std gas vol (ml)")+
    ylab("CH4 Area (FID)")
  
  plt_N2O <- ggplot(subset(calcurve, gas=="n2o"), aes(vol, area))+
    geom_point()+
    # geom_text(aes(label=uniqueareaid))+
    stat_poly_line() +
    stat_poly_eq(use_label(c("eq", "adj.R2", "n")))+
    geom_point(data = subset(excluded, gas=="n2o"), aes(vol,area, col="excluded"))+
    # geom_text(data=subset(qual_inj, gas=="n2o"), aes(vol, area,label=uniqueareaid))+
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
gg_save_pdf(list = plt_list, filename = paste0("calcurves_restore4c.pdf"))


rm(excluded,calcurve,plt_CH4,plt_CO2,plt_N2O,datecalcurve,plt_list)




#Calculate relative deviation of each Pinjection with respect to the daily calfactor

#Do all data in gc_complete_long have a calfactor?
unique(unique(gc_complete_long$date)%in%unique(cal$date))

date_factor<- cal %>% 
  filter(term=="vol") %>% 
  select(date, gas, estimate) %>% 
  rename(factor=estimate)

gc_complete_long_cal<- gc_complete_long %>% 
  merge.data.frame(date_factor, by=c("date","gas"), all = T) %>% 
  mutate(vol_est=area/factor)




#Timeseries of factor stability for each gas 
cal%>%
  filter(term=="vol")%>%
  group_by(gas) %>% 
  mutate(average_factor=mean(estimate)) %>% 
  # filter(r.squared>0.95) %>% 
  ggplot(aes(x=as.Date(as.character(date),format="%Y%m%d"), y=estimate/average_factor, col=gas))+
  geom_hline(yintercept = 1)+
  geom_point()+
  scale_x_date(name="Date of analysis", date_breaks = "1 month", date_minor_breaks = "7 days",date_labels = "%B")+
  scale_y_continuous(name="Rel. dev. from mean")+
  facet_grid(rows = vars(gas))+
  ggtitle("Relative deviation of daily calibration factor")

#PLot of relative calibration factor stability
cal%>%
  filter(term=="vol")%>%
  group_by(gas) %>% 
  mutate(average_factor=mean(estimate)) %>% 
  # filter(r.squared>0.95) %>% 
  ggplot(aes(x=gas, y=estimate/average_factor, fill=gas))+
  geom_violin()+
  geom_hline(yintercept = 1)+
  scale_y_continuous(name="Rel. dev. from mean")+
  geom_dotplot(binaxis= "y",
               stackdir = "center",binwidth=0.01,
               # dotsize = 1,
               fill = 1) +
  ggtitle("Relative deviation of daily calibration factor")

#Timeseries of P5 absolute areas 
gc_complete_long_cal%>%
  filter(vol==5) %>%
  # filter(r.squared>0.95) %>% 
  ggplot(aes(x=as.Date(as.character(date),format="%Y%m%d"), y=area, col=gas))+
  geom_hline(yintercept = 1)+
  geom_point()+
  scale_x_date(name="Date of analysis", date_breaks = "1 month", date_minor_breaks = "7 days",date_labels = "%B")+
  scale_y_continuous(name="Area of peak")+
  facet_grid(rows = vars(gas), scales="free")+
  ggtitle("Timeseries of P5 absolute areas ")


# PLOT absolute area of P5 injections
gc_complete_long_cal %>% 
  filter(vol==5) %>%
  group_by(gas) %>% 
  mutate(avg=mean(area,na.rm = T)) %>% 
  ggplot(aes(x=gas, y=area/avg, fill=gas))+geom_violin()+
  # geom_jitter(width = 0.2,size = 1)+
  geom_dotplot(binaxis= "y",
               stackdir = "center",binwidth=0.025,
               dotsize = 0.5,
               fill = 1) +
  geom_hline(yintercept = 1)+
  scale_y_continuous(name="Area of peak/mean area")+
  ggtitle("P5 deviation from mean (not calibrated with standard curves)")


# PLOT relative deviation of P5 areas for each of the gasses
#Comparison of GHG P5 injection quality
gc_complete_long_cal %>% 
  filter(vol==5) %>%
  ggplot(aes(x=gas, y=vol_est, fill=gas))+geom_violin()+
  # geom_jitter(width = 0.2,size = 1)+
  geom_hline(yintercept = 5)+
  geom_dotplot(binaxis= "y",
               stackdir = "center",binwidth=0.25,
               dotsize = 0.5,
               fill = 1) +
  scale_y_continuous(name="Estimated volume of std gas (ml)")+
  ggtitle("P5 deviation from daily calibration curve")




  
#Timeseries of measured P5 volume of gas
gc_complete_long_cal %>% 
  filter(vol==5) %>%
  ggplot(aes(x=datef, y=vol_est, col=gas, group=datef))+
  geom_hline(yintercept = 5)+
  # geom_boxplot()+
  geom_jitter(size=0.8, width=0.2)+
  scale_x_date(name="Date of analysis", date_breaks = "1 month", date_minor_breaks = "7 days",date_labels = "%B")+
  scale_y_continuous(name="Estimated volume of std gas (ml)")+
  facet_grid(rows = vars(gas))+
  ggtitle("P5 deviation from daily calibration curve")

#Timeseries of relative deviation of P5 measured vs expected volume
gc_complete_long_cal %>% 
  filter(vol==5) %>%
  ggplot(aes(x=datef, y=(vol_est/5)-1, col=gas, group=datef))+
  geom_hline(yintercept = 0)+
  # geom_boxplot()+
  geom_jitter(size=0.8, width=0.2)+
  scale_x_date(name="Date of analysis", date_breaks = "1 month", date_minor_breaks = "7 days",date_labels = "%B")+
  scale_y_continuous(name="Relative deviation from expected value")+
  facet_grid(rows = vars(gas))+
  ggtitle("P5 relative deviation from daily calibration curve")

  

#Expected vs measured volume across all days for the 3 GHG (scatter +jitter)
gc_complete_long_cal %>% 
  ggplot(aes(x=vol, y=vol_est,col=gas))+
  geom_abline(intercept = 0, slope = 1, linewidth=1)+
  geom_jitter(size=0.8, width=0.2)+
  facet_grid(cols=vars(gas))+
  theme_bw()+
  scale_y_continuous(name="Measured vol (ml)", breaks = c(0,1,2,4,6,8,10,12, 15,20,25))+
  scale_x_continuous(name="Actual vol injected (ml), with 0.2ml jitter", breaks = c(0,1,2,4,6,8,10,12))+
  ggtitle("Known vs measured volume in all standard gas injections")
  


#Expected vs measured volume across all days for the 3 GHG (scatter +jitter)
gc_complete_long_cal %>% 
  ggplot(aes(x=vol, y=vol_est,col=gas))+
  geom_abline(intercept = 0, slope = 1, linewidth=1)+
  geom_jitter(size=0.8, width=0.2)+
  facet_grid(rows=vars(gas))+
  theme_bw()+
  scale_y_continuous(name="Measured vol (ml)", breaks = c(0,1,2,4,6,8,10,12, 15,20,25))+
  scale_x_continuous(name="Actual vol injected (ml), with 0.2ml jitter", breaks = c(0,1,2,4,6,8,10,12))+
  ggtitle("Known vs measured volume (1:1 line)")

#Scatter plot (all data) area vs volume injected (with 0.2jitter)
gc_complete_long_cal %>% 
  ggplot(aes(x=vol, y=area,col=gas))+
  geom_jitter(size=0.8, width=0.2)+
  theme_bw()+
  scale_y_continuous(name="Area of peak")+
  scale_x_continuous(name="Actual vol injected (ml), with 0.2ml jitter", breaks = c(0,1,2,4,6,8,10,12))+
  ggtitle("Peak Area vs injected volume")+
  facet_grid(rows=vars(gas), scales = "free")



#Expected vs measured volume across all days for the 3 GHG (boxplot)
gc_complete_long_cal %>% 
  filter(vol!=3) %>% 
  ggplot(aes(x=vol, y=vol_est,col=gas, group=vol))+
  geom_boxplot()+
  geom_abline(intercept = 0, slope = 1)+
  facet_grid(cols=vars(gas))+
  scale_y_continuous(name="Measured vol (ml)", breaks = c(0,1,2,4,6,8,10,12, 15,20,25))+
  scale_x_continuous(name="Actual vol injected (ml), with 0.2ml jitter", breaks = c(0,1,2,4,6,8,10,12))+
  ggtitle("Known vs measured volume in all standard gas injections")

