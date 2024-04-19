#TEst to import raw .dat files from GC


#All chromatograms are batch-exported to csv format with the OpenChrom software. This first step produces 2 csv files per sample:

# - datetime_samplecode.csv   which contains the FID data (CO2 and CH4)
# - datetime_samplecode_REF1.csv  which contains the ECD data (N2O)

# All csv files are in the next folder:

csvfolder <- "C:/Users/Miguel/Desktop/GC_new/GCexportedfiles/"

filelist<- list.files(path = csvfolder, full.names = F)


#Following in this script are test to process these 22 types of data with different packages.

#some are written in python, which requires the package reticulate to be run in R as an interface
install.packages('reticulate')

library(tidyverse)

#####1_RpeakChrom pkg####

# install.packages('RpeakChrom')

library(RpeakChrom)

test<-read.csv("C:/Users/Miguel/Desktop/GC_new/chromatograma ejemplo.csv")
str(test)

testch<-RpeakChrom::readChrom(filepath = "C:/Users/Miguel/Desktop/GC_new/chromatograma ejemplo.csv", do.plot = F, t1=10000, t2=50000)

str(testch)
str(test)

ggplot(test, aes(x=RT.milliseconds., y=TIC))+geom_line()

RpeakChrom::processPeak(peak = testch,baseline = T,method ='pvmg' ,compound = 's',flow = 5,area = T)

install.packages("chromConverter")

