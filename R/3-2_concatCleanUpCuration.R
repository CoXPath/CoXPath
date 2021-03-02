#setwd("/home/traschka/Projects/PhD/phd-project/zml-project/data/")
#setwd("C://Users/traschka/Documents/Projects/PhD/zml-project/data/")
setwd("/home/rebecafig/Documents/fraunhofer/co-expression-vs-pathways/data/")

options(stringsAsFactors=FALSE)

rm(list=ls())

library("dplyr")

disCurated_tamara = read.csv("diseaseCuration_tamara_FINAL.tsv", sep = "\t")
disCurated_sumit = read.csv("diseaseCuration_sumit_FINAL.tsv", sep = "\t")

# only workaround, delete later!
#disCurated_tamara = disCurated_tamara[597:nrow(disCurated_tamara),]
#disCurated_sumit = disCurated_sumit[1:596,]

disCurated_sumit[ disCurated_sumit == "" ] = NA
disCurated_tamara[ disCurated_tamara == "" ] = NA

allCurated = rbind(disCurated_sumit, disCurated_tamara)

#write.table(allCurated, file="diseaseCuration_FINAL.tsv", quote=F, sep="\t", row.names = F)

###---------------- work with final Curation -------------------
rm(list=ls())
allCurated = read.csv(file="diseaseCuration_FINAL.tsv", sep="\t")
## this file includes now manual curations from evaluation of double mapped terms

###---------------- get datasets to delete ---------------------
datToDelete = subset(allCurated, subset = (text_index=="expTitle") & (Curation.result == "delete"))
datToDelete = unique(datToDelete$acc_dataset)

write(datToDelete, ncolumns = 1, file="datasetsToDelete.txt", )

allCurated = subset(allCurated, subset = !(acc_dataset%in%datToDelete))

###---------------- get mapped terms ---------------------------

mappedTerms = subset(allCurated, select = c("acc_dataset", "text_index", "text_value"))
mappedTerms$curatedTerm = NA
mappedTerms$curatedID = NA
for (i in 1:NROW(allCurated)) {
  curRes = allCurated$Curation.result[i]
  disTerm = c()
  disID = c()
  if(is.na(curRes)) {
    next
  }
  
  if(curRes=="correct") {
    disTerm = allCurated$predicted_label[i]
    disID = allCurated$predicted_id[i]
  } else if(curRes=="curated") {
    disTerm = allCurated$Curated.label[i]
    disID = allCurated$Curated.ID[i]
    if ( (is.na(disTerm) | is.na(disID))) {
      if (allCurated$Comment[i] != "no disease in title") {
        print(paste0("problem ", i))
      }
    }
  } else if(curRes=="delete") {
    disTerm = "delete"
    disID = NA
  } else if(curRes=="incorrect") {
    disTerm = "incorrect"
    disID = NA
  } else if(curRes=="normal") {
    disTerm = "normal"
    disID = NA
  } else {
    print(curRes)
    print("Problem")
  }
  mappedTerms$curatedTerm[i] = disTerm
  if (exists("disID")) {
    mappedTerms$curatedID[i] = disID
  }
}

###---------------- check mapped terms for duplicates ---------------------------

mappedTerms2 = subset(mappedTerms, subset = curatedTerm!="delete")
mappedTerms2 = subset(mappedTerms2, subset = curatedTerm!="incorrect")
mappedTerms2 = subset(mappedTerms2, select = c("text_value", "curatedTerm"))

uniqueRowTerms = mappedTerms2 %>% distinct()

whichTermsMultiple = which(duplicated(uniqueRowTerms$text_value))
whichTermsMultiple = unique(c(whichTermsMultiple, which(duplicated(uniqueRowTerms$text_value, fromLast = TRUE))))

doubleTerms = uniqueRowTerms[whichTermsMultiple,]
doubleTerms = doubleTerms[order(doubleTerms$text_value),] 
doubleTerms

# see reasoning in decisionsCuration.txt
subset(allCurated, text_value=="ependymoma") ### ok
subset(allCurated, text_value=="tumor") ### ok


###---------------- check mapped terms for other than DOID termIDs ---------------------------
# e.g. MONDO
mappedTerms3 = subset(mappedTerms, subset = curatedTerm!="delete")
mappedTerms3 = subset(mappedTerms3, subset = curatedTerm!="normal")

indxNonDOID = which(!grepl(pattern = "DOID", mappedTerms3$curatedID))

nonDOID = mappedTerms3[indxNonDOID,]
nrow(nonDOID)==0

###---------------- write one curation file ---------------------------
save(mappedTerms, file="diseaseCuration_mappedTerms_FINAL.RData")

