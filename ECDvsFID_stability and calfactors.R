#Quality control ECD vs FID detectors


#Author: Miguel Cabrera. 

#Description: This script uses the GC excel spreadsheet as input to calculate the calibration factor for each day of analysis as well as the stability of P5 standards. 

library(tidyverse)
library(readxl)
library(broom)




dropbox_root<- "C:/Users/Miguel/Dropbox/RESTORE4Cs - Fieldwork/"
gcfolder<- "Data/GHG/GC data/"


gcog<- read_xlsx(path = paste0(dropbox_root,gcfolder,"DATA_R4Cs_formatted_ok_v20240615.xlsx"),range = cell_cols("A:I"))
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
  mutate(uniqueid=paste(rownames(gc_work),standardgas,date,batch,sep = "_"))#add unique id with explicit info


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



#Calculate all the cal curves
cal<- gc_complete %>%
  filter(batch==0) %>% 
  filter(vol!=11.2) %>% #NO P20s!
  filter(vol!=0) %>%  #NO blanks!
  select(uniqueid,date,vol,co2,ch4,n2o)%>%
  pivot_longer(cols = c(co2,ch4,n2o), names_to = "gas",values_to = "area") %>% 
  nest_by(date,gas) %>% 
  mutate(mod = list(lm(area~vol, data=data))) %>%
  reframe(tidy(mod),glance(mod))
  
  
#Plot factors 
cal%>%
  filter(term=="vol")%>%
  filter(r.squared>0.95) %>% 
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

