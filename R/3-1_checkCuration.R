setwd("/home/traschka/Projects/PhD/phd-project/zml-project/data/")
setwd("C://Users/traschka/Documents/Projects/PhD/zml-project/data/")

options(stringsAsFactors=FALSE)

rm(list=ls())

library("dplyr")

disCurated_tamara = read.csv("diseaseCuration_tamara_with_ids.tsv", sep = "\t")
disCurated_sumit = read.csv("diseaseCuration_sumit_with_labels.tsv", sep = "\t")

disCurated_sumit[ disCurated_sumit == "" ] = NA
disCurated_tamara[ disCurated_tamara == "" ] = NA

###---------------- get together files Sumit and Tamara ---------------------------
## check twice curated lines
curLinesTa = setdiff(c(1:NROW(disCurated_tamara)), which(is.na(disCurated_tamara$Curation.result)))
curLinesSu = setdiff(c(1:NROW(disCurated_sumit)), which(is.na(disCurated_sumit$Curation.result)))

twiceCur = intersect(curLinesTa, curLinesSu)

twiceData = rbind(disCurated_tamara[twiceCur,], disCurated_sumit[twiceCur,])

twiceData = twiceData %>% distinct() #remove exactly matching rows, but keep unique one'

twiceDataSub = subset(twiceData, select = c(acc_dataset, text_value))
whichLinesMultiple = which(duplicated(twiceDataSub))
whichLinesMultiple = unique(c(whichLinesMultiple, which(duplicated(twiceDataSub, fromLast = TRUE))))

doubleLines = twiceData[whichLinesMultiple,]
doubleLines = doubleLines[order(doubleLines$acc_dataset, doubleLines$text_value),] 
doubleLines

doubleLines_exp = unique(doubleLines$acc_dataset)
doubleLines_subTa = subset(disCurated_tamara, acc_dataset%in%doubleLines_exp)
doubleLines_subSu = subset(disCurated_sumit, acc_dataset%in%doubleLines_exp)

doubleLines_both = rbind(doubleLines_subSu, doubleLines_subTa)
doubleLines_both = doubleLines_both %>% distinct() #remove exactly matching rows, but keep unique one'
doubleLines_both = doubleLines_both[order(doubleLines_both$acc_dataset, doubleLines_both$text_value),] 


write.table(doubleLines_both, file="mismatchSumitTamara.tsv", quote=F, sep="\t", row.names = F)

