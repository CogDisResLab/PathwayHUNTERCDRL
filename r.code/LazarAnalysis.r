#DrLazar dataset
#Proteomics dataset were we set to confirm
library(tidyverse)
library(readr)

#reading in Dr. Lazar dataset
Ldata<- read_csv("LazarData.csv")

#Using regex to extract GOIDs from column 2
library(dplyr)
GOID <- data.frame(do.call('rbind', strsplit(as.character(Ldata$Term),'~',fixed=TRUE)))
GOID <- data.frame(GOID$X1)
names(GOID)[names(GOID) == "GOID.X1"] <- "GOID"
myGOID <- data.frame(GOID)
myGOID <- data.frame(lapply(myGOID, as.character), stringsAsFactors=FALSE)

#Run PathwayHUNTER
PathwayHUNTER(myGOID,8,"LazarTest.csv",2)


LLdata<- read_csv("Lazar2.csv")

#Using regex to extract GOIDs from column 2
library(dplyr)
myid <- LLdata$`Biological Process`
GOID <- data.frame(myid)
names(GOID)[names(GOID) == "myid"] <- "GOID"
myGOID <- data.frame(GOID)
myGOID <- data.frame(lapply(myGOID, as.character), stringsAsFactors=FALSE)

#Run PathwayHUNTER
PathwayHUNTER(myGOID,8,"LazarTest2.csv",2)
