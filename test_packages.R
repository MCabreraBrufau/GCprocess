#Test to process files from GC#


#All chromatograms are batch-exported to csv format with the OpenChrom software. This first step produces 2 csv files per sample:

# - datetime_samplecode.csv   which contains the FID data (CO2 and CH4)
# - datetime_samplecode_REF1.csv  which contains the ECD data (N2O)

# All csv files are in the next folder:
csvfolder <- "C:/Users/Miguel/Desktop/GC_new/GCexportedfiles/"

#Test csv files in googledrive:

  #Thinkpad:
csvfolder<- "G:/Mi unidad/__RESTORE4Cs__/GC/Examples_csvformat/"





filelist<- list.files(path = csvfolder, full.names = T)
filelist

ecdfiles<- filelist[grepl(pattern = "REF1",x = filelist)]
fidfiles<- filelist[!grepl(pattern = "REF1",x = filelist)]



library(tidyverse)

#Following in this script are test to process these 2 types of data with different packages. I will try the following packages: 

# RpeakChrom
# chromatographR
# xcms

#None of these packages offers an acceptable workflow for our case, although they have useful functions and will be used as a guide for our custom pipeline. 



#Basis of peak indeces:
# https://www.agilent.com/cs/library/technicaloverviews/public/5990-7651EN.pdf

#Definitions: 

# signal: maximum value of chromatogram
# noise: sd of background noise calculated 60s before the peak or 30 seconds before and after the peak. Check for time window of other systems.



#What do we want: 

#Pre-processing: baseline correction
#Detection: peak detection (start&end RTs, RT of max)
#Area under peak
#General metrics: absolute noise (sd of TIC across RT windows, excluding the peak)

#Peak quality metrics: signal-to-noise ratio,  peak significance level , (sharpness, Gaussian similarity) maybe





#How are we going to achieve it: 
#split chromatograms to get 1 peak per vector: CO2, CH4 and N2O in separate vectors based on RTs (be generous): RT always in minutes (2nd column in openChrom exported csv files)

#Import calibration curve csv files as example to fit RT windows for each GHG
for (i in 1:length(fidfiles)){
  a<- read.csv(fidfiles[i])[,c(2,4)]
  b<- read.csv(ecdfiles[i])[,c(2,4)]
  names(a)<- c("RT", paste0("FID_","sample",i))
  names(b)<- c("RT", paste0("ECD_","sample",i))
  if(i==1){A<-merge.data.frame(a,b, by="RT")}else{A<-merge.data.frame(A,a, by="RT")%>%merge.data.frame(b,by="RT")}
}
rm(a,b,i)

A_long<-A%>%
  pivot_longer(-RT, names_to = "var",values_to = "value")%>%
  separate(var, into = c("sensor", "sampleno"))

#Set RT windows (in minutes) for each gas
co2RTs<- c(4.2,5.6)
ch4RTs<- c(2.65,3.2)
n2oRTs<- c(5.75,6.5)

#Plot RT windows to confirm:

#FID sensor
A_long %>% filter(sensor=="FID"&RT<n2oRTs&RT>2) %>% filter(RT%%0.005<.0002) %>% 
ggplot(aes(x=RT, y=value, group=paste(sensor,sampleno)))+
  geom_rect(aes(xmin=co2RTs[1], xmax=co2RTs[2],ymin=min(value), ymax=max(value), fill="CO2"), alpha=0.1)+
  geom_rect(aes(xmin=ch4RTs[1], xmax=ch4RTs[2],ymin=min(value), ymax=max(value), fill="CH4"), alpha=0.1)+
  geom_line()+
  facet_wrap(facets = ~sensor, scales = "free")


#ECD sensor
A_long %>% filter(sensor=="ECD"&RT>c(co2RTs)&RT<6.6)%>% filter(RT%%0.005<.0005) %>% 
  ggplot(aes(x=RT, y=value, group=paste(sensor,sampleno)))+
  geom_rect(aes(xmin=n2oRTs[1], xmax=n2oRTs[2],ymin=min(value), ymax=max(value), fill="N2O"), alpha=0.7)+
  # geom_vline(xintercept = c(n2oRTs))+
  geom_line()+
  facet_wrap(facets = ~sensor, scales = "free")


read.csv(fidfiles[1])[,c(2,4)] %>% 
  rename(RT=names(.)[1]) %>% 
  filter(between(RT, 2.65,3.2))%>%
  plot()

read.csv(ecdfiles[5])[,c(2,4)] %>% 
  rename(RT=names(.)[1]) %>% 
  filter(between(RT, 5.75,6.5))%>%
  plot()


#1. Baseline correction: 

#function ptw::baseline.corr(testb[, 2])
library(ptw)
# asysm: Trend estimation with asymmetric least squares 
baseline.corr(y = )# this function uses the asysm function to provide a baseline-corrected chromatogram

asysm(y,#data to correct, either as vector or as matrix with each chromatogram as a row
      lambda = 1e7, #smoothing parameter
      p=0.001, #asymetry parameter
      eps=1e-8, #numerical precission for convergence
      maxit = 25)#max number of iterations (warns when no convergence reached)


#Default seems to work relatively well for CH4 when RT windows are supplied
A%>%
  filter(between(RT,ch4RTs[1],ch4RTs[2]))%>%
ggplot(aes(x=RT))+
  geom_line(aes(y=FID_sample7, col="OG"))+
  geom_line(aes(y=baseline.corr(FID_sample7,lambda=1e10, p=0.0001), col="bs-corr"))+
  scale_y_continuous(limits = c(-1000,3e4))


#Baseline correction does not work for CO2 correctly, even after trying different parameters. Probably too asymetric a shape, the function tries to lower right hand tail of peak to zero, but this is not correct. 
A%>%
  filter(between(RT,co2RTs[1],co2RTs[2]))%>%
  ggplot(aes(x=RT))+
  geom_line(aes(y=FID_sample5, col="OG"))+
  geom_line(aes(y=baseline.corr(FID_sample5, lambda=1e12, p=0.00001), col="bs-corr"))+
  scale_y_continuous(limits = c(-10000,1e5))

#As the CO2 peak is always very clear and has a long baseline before the peak, we could do the baseline correction "manually", taking the mean value of the 0.3 minutes before the peak as the zero value, as in: 
A %>% 
  filter(between(RT,co2RTs[1],co2RTs[2]))%>%
  mutate(basecFID_sample7=FID_sample7-mean(FID_sample7[1:360])) %>% 
  ggplot(aes(x=RT))+
  geom_line(aes(y=FID_sample7, col="OG"))+
  geom_line(aes(y=basecFID_sample7, col="bs-corr"))+
  geom_rect(aes(xmin=co2RTs[1], xmax=4.5, ymin=-10000, ymax=10000, fill="bs-corr"), alpha=0.2)+
  # geom_line(aes(y=baseline.corr(FID_sample3, lambda=1e12, p=0.00001), col="bs-corr"))+
  scale_y_continuous(limits = c(-10000,1e6))




#N2O Baseline correction 
A%>%
  select(RT,ECD_sample3)%>%
  filter(between(RT,n2oRTs[1],n2oRTs[2]))%>%
  rename(OG=names(.)[2]) %>% 
  mutate(basecorr=baseline.corr(OG))%>%
  pivot_longer(cols = -RT, names_to = "chrom", values_to = "value") %>% 
  ggplot(aes(x=RT,y=value, col=chrom))+
  geom_line()+
  facet_grid(rows = vars(chrom),scales = "free")




data("gaschrom")

data(gaschrom)


#Peak detection: 
# find and integrate peaks using exponential-gaussian hybrid model (same approach as in chromatographR package)
#The "finding" part in the chromatographR package is too vague and produces 100s of peaks with very small RT windows instead of what we want: 1 wide peak (+ some "peaks" resulting from noise)




#####1. RpeakChrom pkg####

# 
if(!require(RpeakChrom)){
  install.packages('RpeakChrom')
}


test<-read.csv(filelist[2], header = F)
str(test)

######1.1 import######
#RpeakChrom requires only 2 columns in csv, RT and Intensity, without headers. 

testa<-RpeakChrom::readChrom(filepath = filelist[2], do.plot = T, t1=1.5e5, t2=3.5e5)
#Import function does not allow for headers in csv... causes all data to be imported as character
str(testa)



maxRT<-1539552

plot(testb)
plot(ptw::baseline.corr(testb[, 2]))



#try using an example file without headers
testb<-RpeakChrom::readChrom(filepath = filelist[2], do.plot = T, t1=1.65e5, t2=1.85e5)

str(testb)
#removing headers works


max(testb$V2)

RpeakChrom::processPeak(peak = testb,baseline = T,method ='direct',area = T, compound = "a", flow=F)


#something fails with the integration/modeling of the peak


parameters <- processPeak(peak, baseline = FALSE, flow = 0.1, method = "pvmg",
                          compound = "alanine", area = TRUE)

rm(testa, testb,test)


####2. chromatographR pkg####

# chromatographR is a package for the reproducible analysis of HPLC-DAD chromatographic data in R. It can also be used to analyze other "simple" chromatographic data like GC-FID, HPLC-UV, or HPLC-FD.

#Install:
if(!require(chromatographR)){
  if(!require(remotes)){
    install.packages("remotes") 
  }
  remotes::install_github("https://github.com/ethanbass/chromatographR/")
}





test<-chromatographR::read_chroms(paths = filelist[1],#path to files or folders containing files
                            find_files = F,#Logical. Set to TRUE (default) if you are providing the function with a folder or vector of folders containing the files. Otherwise, set toFALSE.
                            format_in = "other",#Format of files to be imported/converted. The current options are: chemstation_uv, chemstation_csv, masshunter_dad, shimadzu_fid, shimadzu_dad, chromeleon_uv, thermoraw, mzml, waters_arw, msd, csd, wsd, or other. 
                            #MIGUEL: current version of openchrom does not allow external packages to use their parsers, so we chose "other" and provide .csv files batch-exported from the openchrom software.
                            pattern =".csv",
                            parser = "openchrom",
                            format_out = "data.frame",
                            export = F,
                            )

sa<-Sa




read.csv(filelist[1])[,c(1,4)] %>% 
  filter(between(RT.milliseconds., 1.65e5,3.3e5))%>%
  plot()

test_plain<- read.csv(filelist[1])[,c(1,4)] %>% filter(between(RT.milliseconds., 1.65e5,3.3e5))
test_plain<- column_to_rownames(test_plain, "RT.milliseconds.")
test_matrix<- as.matrix(test_plain)


tpoints <- as.numeric(rownames(test_plain))
lambda <- 'TIC'

matplot(x = tpoints, y = test_plain[,lambda], 
        type = 'l', ylab = 'TIC', xlab = 'Time (milliseq)')
matplot(x = tpoints, y = ptw::baseline.corr(test_plain[,lambda], p = .0001, lambda = 1e9,maxit=50),
        type = 'l', add = TRUE, col='blue', lty = 3)







test_preprocessed<-chromatographR::preprocess(test_matrix, remove.time.baseline = T,dim1 = seq(1.65e5,1.85e5, by=50), dim2="TIC",spec.smooth = F,interpolate_rows = T,interpolate_cols = F)

read.csv(filelist[1])[,c(1,4)] %>% 
filter(between(RT.milliseconds., 1.65e5,1.85e5))%>%
  plot()

plot(test_preprocessed)

# find and integrate peaks using gaussian peak fitting
pks_gauss <- get_peaks(test_preprocessed,
                       lambdas = c("TIC"),
                       sd.max = 500,
                       fit = "gaussian")


# find and integrate peaks using exponential-gaussian hybrid model
pks_egh <- get_peaks(test_preprocessed, lambdas = c("TIC"), sd.max = 100,time.units = "ms", fit="egh")
 
# find and integrate peaks without modeling peak shape
pks_raw <- get_peaks(test_preprocessed, lambdas = c("TIC"), sd.max = 1000, fit="raw", )


read.csv(filelist[1])[,c(1,4)] %>% 
  filter(between(RT.milliseconds., 1.65e5,1.85e5))%>%
  mutate(der1=c(diff(TIC),NA)) %>% 
  ggplot(aes(x=RT.milliseconds., y=der1))+geom_line()+
  geom_line(aes(y=TIC), col="red")




read.csv(filelist[3])[,c(1,4)] %>% 
  ggplot(aes(x=RT.milliseconds., y=TIC))+geom_line()+
  geom_line(aes(y=TIC), col="red")


matplot(x = as.numeric(rownames(test_plain)), y = test_plain[,'TIC'], 
        type = 'l', ylab = 'TIC', xlab = 'Time (milliseq)')

matplot(x = as.numeric(rownames(test_plain)), y = ptw::baseline.corr(test_plain[,'TIC'], p = .001, lambda = 1e5,maxit=20),
        type = 'l', add = TRUE, col='blue', lty = 3)

plot(test_plain[,'TIC'])
plot(ptw::baseline.corr(test_plain[,'TIC'], p = .0001, lambda = 1e9,maxit=15))



read.csv(filelist[3])[,c(2,4)] %>% 
  rename("RT"=names(.)[1]) %>%  
  filter(between(RT, 5.7,6.5))%>%
  plot()

test_plain<- read.csv(filelist[3])[,c(2,4)] %>% 
  rename("RT"=names(.)[1]) %>%    filter(between(RT, 5.7,6.5))
test_plain<- column_to_rownames(test_plain, "RT")
test_matrix<- as.matrix(test_plain)

test_preprocessed<-chromatographR::preprocess(test_matrix, remove.time.baseline = T,dim1 = seq(5.7,6.5, by=0.001), dim2="TIC",spec.smooth = F,interpolate_rows = T,interpolate_cols = F,lambda=1e9)


plot(test_plain[,'TIC'])
plot(test_preprocessed)


# find and integrate peaks using exponential-gaussian hybrid model
pks_egh <- get_peaks(test_preprocessed, lambdas = c("TIC"), fit="egh",time.units = "min",sd.max = 1)
pks_gauss <- get_peaks(test_preprocessed,lambdas = c("TIC"), fit = "gaussian", time.units="min",sd.max=1e9)
pks_raw <- get_peaks(test_preprocessed, lambdas = c("TIC"), fit="raw", )



# find and integrate peaks using exponential-gaussian hybrid model
pks_egh <- get_peaks(test_preprocessed, lambdas = c("TIC"), sd.max = 100,time.units = "ms", fit="egh")

# find and integrate peaks without modeling peak shape
pks_raw <- get_peaks(test_preprocessed, lambdas = c("TIC"), sd.max = 1000, fit="raw", )


par(olpar)
par(oldpar)

pks_gauss[[1]]$TIC %>% 
  filter(FWHM==max(FWHM))


####3. xcmx pkg####

#the xcmx package belongs to the bioconductor universe, it was developed to work with LC-MS and GC-MS data, but we can make use of its functionalities for basic chromatographic processinga
#install:

if (!require("xcms", quietly = F)){
  if (!require("BiocManager", quietly = F)){
    install.packages("BiocManager")
  }
  BiocManager::install("xcms")
}

require(xcms)



