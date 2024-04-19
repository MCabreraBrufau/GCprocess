#TEst to import raw .dat files from GC


#All chromatograms are batch-exported to csv format with the OpenChrom software. This first step produces 2 csv files per sample:

# - datetime_samplecode.csv   which contains the FID data (CO2 and CH4)
# - datetime_samplecode_REF1.csv  which contains the ECD data (N2O)

# All csv files are in the next folder:

csvfolder <- "C:/Users/Miguel/Desktop/GC_new/GCexportedfiles/"

filelist<- list.files(path = csvfolder, full.names = F)


#Following in this script are test to process these 2 types of data with different packages. I will try the following packages: 

# RpeakChrom
# chromatographR
# xcms



#####1. RpeakChrom pkg####

# 
if(!require(RpeakChrom)){
  install.packages('RpeakChrom')
}


test<-read.csv("C:/Users/Miguel/Desktop/GC_new/chromatograma ejemplo.csv")
str(test)

######1.1 import######
testa<-RpeakChrom::readChrom(filepath = "C:/Users/Miguel/Desktop/GC_new/chromatograma ejemplo.csv", do.plot = T, t1=20, t2=500)
#Import function does not allow for headers in csv... causes all data to be imported as character
str(testa)



maxRT<-1539552
#try using an example file without headers
testb<-RpeakChrom::readChrom(filepath = "C:/Users/Miguel/Desktop/GC_new/chromatograma ejemplo_noheader.csv", do.plot = T, t1=1.65e5, t2=1.85e5)

str(testb)
#removing headers works


max(testb$V2)

processPeak(peak = testb,baseline = T,method ='direct',area = T, compound = 'text', 
flow = 'text2')
#something fails with the integration/modeling of the peak


parameters <- processPeak(peak, baseline = FALSE, flow = 0.1, method = "pvmg",
                          compound = "alanine", area = TRUE)


install.packages("chromConverter")
plot(peak)
plot(testb)
peak

####2. chromatographR pkg####

# chromatographR is a package for the reproducible analysis of HPLC-DAD chromatographic data in R. It can also be used to analyze other "simple" chromatographic data like GC-FID, HPLC-UV, or HPLC-FD.

#Install:
if(!require(chromatographR)){
  if(!require(remotes)){
    install.packages("remotes") 
  }
  remotes::install_github("https://github.com/ethanbass/chromatographR/")
}










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




