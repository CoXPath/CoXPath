#setwd("/home/traschka/Projects/PhD/phd-project/zml-project/data/")
#setwd("C://Users/traschka/Documents/Projects/PhD/zml-project/data/")
setwd("/home/rebecafig/Documents/fraunhofer/co-expression-vs-pathways/data/")

options(stringsAsFactors=FALSE)

rm(list=ls())

library("dplyr")

load("sampleDataClean.Rdata")
load("sampleDataALL_disease.Rdata")
load("diseaseCuration_mappedTerms_FINAL.RData")
### -> 51550 samples

###---------------- get only disease columns ---------------------------
whichColsDisease = grep(colnames(sampleDataALL_disease), pattern = "disease", ignore.case = TRUE)

# uninformative:
excludeDisease = c("diseased.free.survival..1..disease.", "diseased.free.survival..months.", "pre.transplant.disease", 
                   "clinical.status.at.last.follow.up.awd...alive.with.disease..dod...dead.from.disease..ned...no.evidence.of.disease",
                   "clinical.status.post.1st.line.chemotherapy.cr...complete.response..pr...partial.response..sd...stable.disease..p...progression",
                   "clinical.status.at.last.follow.up.awd...alive.with.disease..dod...dead.from.disease..ned...no.evidence.of.disease.1",
                   "clinical.status.post.1st.line.chemotherapy.cr...complete.response..pr...partial.response..sd...stable.disease..p...progression.1",
                   "dfs...disease.free.survival.days", "dfs...disease.free.survival.days.1", "dfs_event.disease.free.survival", 
                   "dfs_time.disease.free.survival.time..months", "dss_event.disease.specific.survival", "dss_time.disease.specific.survival.time..months",
                   "dfs_event.disease.free.survival.1", "dfs_time.disease.free.survival.time..months.1", "dss_event.disease.specific.survival.1",
                   "dss_time.disease.specific.survival.time..months.1", "disease.duration", "DISEASE.DURATION", "disease.free.survival.time.dfs",
                   "DISEASE.FREE.SURVIVAL.TIME.DFS", "disease.activity.score.das.28crp", "DISEASE.ACTIVITY.SCORE.DAS.28CRP", "disease.free.survival.time",
                   "DISEASE.FREE.SURVIVAL.TIME", "extent.of.disease", "EXTENT.OF.DISEASE", "DiseaseLocation", "residual.disease",
                   "code.disease.specific.survival", "disease.specific.survival.years", "DISEASE.SPECIFIC.SURVIVAL.YEARS", 
                   "Day.14.Minimal.Residual.Disease", "Day.28.Minimal.Residual.Disease", "Day.7.Minimal.Residual.Disease",
                   "Disease.Free.Survival", "Disease.Free.Survival.Status", "duration.of.disease..month.", "DURATION.OF.DISEASE..MONTH.", 
                   "disease.progression", "DISEASE.PROGRESSION", "Disease.free.interval..DFI..in.months", "DISEASE.FREE.INTERVAL..DFI..IN.MONTHS",
                   "diseased", "DISEASED")
excludeDisease_idx = which(colnames(sampleDataALL_disease)%in%excludeDisease)

keepCols = whichColsDisease[!(whichColsDisease%in%excludeDisease_idx)]

sampleDataALL_disease = subset(sampleDataALL_disease, select = c(sampleID, acc_sample, acc_dataset, keepCols))


###---------------- check for datasets to delete ---------------------------
whichDatDel = readLines("datasetsToDelete.txt")

sampleDataClean2 = subset(sampleDataClean, subset = !(acc_dataset%in%whichDatDel))
sampleDataALL_disease2 = subset(sampleDataALL_disease, subset = !(acc_dataset%in%whichDatDel))
### -> 46729 samples of 536 datasets

###---------------- prep for disease mapping ---------------------------
sampleDataClean2[ sampleDataClean2 == "NA" ] = NA
sampleDataALL_disease2[ sampleDataALL_disease2 == "NA" ] = NA
sampleDataClean2[ sampleDataClean2 == "none" ] = NA
sampleDataALL_disease2[ sampleDataALL_disease2 == "none" ] = NA
sampleDataClean2[ sampleDataClean2 == "not specified" ] = NA
sampleDataALL_disease2[ sampleDataALL_disease2 == "not specified" ] = NA
sampleDataClean2[ sampleDataClean2 == "NOT SPECIFIED" ] = NA
sampleDataALL_disease2[ sampleDataALL_disease2 == "NOT SPECIFIED" ] = NA
sampleDataClean2[ sampleDataClean2 == "N/A" ] = NA
sampleDataALL_disease2[ sampleDataALL_disease2 == "N/A" ] = NA

sampleDataClean2[ sampleDataClean2 == "?1" ] = 1
sampleDataALL_disease2[ sampleDataALL_disease2 == "?1" ] = 1

sampleDataClean2$diseaseTerm = NA
sampleDataClean2$diseaseID = NA

controlVals = c("control", "normal", "normal_culturecells", "Healthy control", 
                "Normal (no pulmonary hypertension)", "Healthy Control", "normal ",
                "healthy", "healthy control", "normal control", "pooled normal CD34+ cells",
                "normal CD3+ cells", "normal CD19+ cells", "normal bone marrow", "normal peripheral blood",
                "normal lung tissue", "normal frontal", "normal cerebellum", "Normal colonic mucosa",
                "Healthy donor", "healthy adults undergoing elective surgery (abdominoplasty)",
                "Normal", "adjacent normal tissue", "non-SS control", "healthy donor",
                "Neurologically normal control", "Healthy", "normal (healthy control)")
textExclude = c("1", "2", "3")


###---------------- choose the correct disease term for each sample ---------------------------
pb <- txtProgressBar(min = 0, max = NROW(sampleDataALL_disease2)-1, style = 3)

suppressMessages(
  for (i in 1:NROW(sampleDataALL_disease2)) {
    if (i==31359) {
      next
    }
    if (sampleDataALL_disease2$acc_dataset[i]=="E-TABM-945") {
      next
    }
    if (sampleDataALL_disease2$acc_dataset[i]=="E-GEOD-8581") {
      next
    }
    whichSample = sampleDataALL_disease2$sampleID[i]
    colsNotNA = colnames(sampleDataALL_disease2)[!is.na(sampleDataALL_disease2[i,])]
    colString = setdiff(colsNotNA, c("sampleID", "acc_sample", "acc_dataset"))
    textToMap = unique(as.character(sampleDataALL_disease2[i,colString]))
    textToMap = setdiff(textToMap, textExclude)
    whichCrazyText = which(startsWith(textToMap, "cell line"))
    if (length(whichCrazyText)!=0) {
      textToMap = textToMap[-whichCrazyText]
    }
    if(length(textToMap)==1) {
      if (textToMap%in%controlVals) {
        finalDisease_term = "normal"
        finalDisease_ID = NA
      } else {
        mappedTerms_sub = subset(mappedTerms, subset = acc_dataset==sampleDataALL_disease2$acc_dataset[i])
        whichRow = which(mappedTerms_sub$text_value == textToMap)
        finalDisease_term = mappedTerms_sub$curatedTerm[whichRow]
        finalDisease_ID = mappedTerms_sub$curatedID[whichRow]
        if (any(c(is.na(finalDisease_term), is.na(finalDisease_ID)))) {
          if (finalDisease_term=="delete") {
            next
          } else if (finalDisease_term=="incorrect") {
            next
          } else if (finalDisease_term!="normal") {
            print("stop1")
            print(i)
            break
          }
        }
      }
    } else {
      finalDisease_term = NA
      finalDisease_ID = NA
    }
    
    
    
    whichLineClean = which(sampleDataClean2$sampleID==whichSample)
    if(whichLineClean!=i) {
      print("mismatch i and whichLine")
    }
    
    sampleDataClean2$diseaseTerm[whichLineClean] = finalDisease_term
    sampleDataClean2$diseaseID[whichLineClean] = finalDisease_ID
    
    setTxtProgressBar(pb, i-1)
    
    rm(list=setdiff(ls(), c("pb", "i", "mappedTerms", "sampleDataALL_disease2", "sampleDataClean2",
                            "controlVals", "textExclude")))
  }
)

close(pb)

##get multiple columns if multiple terms/ids are present
sampleDataClean3 = sampleDataClean2
sampleDataClean3$diseaseTerm2 = NA
sampleDataClean3$diseaseID2 = NA
whichMultiple = grep(x=sampleDataClean3$diseaseTerm, pattern=";")
all(whichMultiple == grep(x=sampleDataClean3$diseaseID, pattern=";")) #all terms that have ; also have ID with ;
for (i in whichMultiple) {
  splitTerms = unlist(strsplit(sampleDataClean3$diseaseTerm[i], ";"))
  splitIDs = unlist(strsplit(sampleDataClean3$diseaseID[i], ";"))
  if (length(splitIDs)!=length(splitTerms)) {
    print("different number of terms and IDs")
    break
  }
  if (length(splitTerms)==2) {
    sampleDataClean3$diseaseTerm[i] = splitTerms[1]
    sampleDataClean3$diseaseTerm2[i] = splitTerms[2]
    sampleDataClean3$diseaseID[i] = splitIDs[1]
    sampleDataClean3$diseaseID2[i] = splitIDs[2]
  } else {
    break
  }
}

## write final metadata file
metadataFinal = sampleDataClean3
save(metadataFinal, file="metadata_final.Rdata")  ##### 2020-11-02
