setwd("/home/rebecafig/Documents/fraunhofer/co-expression-vs-pathways/data/")

library(dplyr)
library(data.table)
library(jsonlite)

options(stringsAsFactors=FALSE)

######---------remove interventional datasets and just keep patient data--------------

sampleDataALL_patient = data.frame()
removeColsALL_string = c()
removedExperiments = c()

for (i in 0:9) {
  ## load parts
  load(paste0("sampleData_Part", i, ".Rdata"))
  #sampleDataALL[[i+1]] = sampleData
  #colnamesALL = c(colnamesALL, colnames(sampleData))
  
  
  ### the used data should only include patient data, no interventional data or treated patients.
  ## Herefore, indicators for interventional data are needed. These are for examples: cellline (as these don't come from patients) or 
  ##  compound or dose (as this indicates, that something (a compound) was put on the samples with a specific dosage)
  
  ## exclude experiments, that have values in columns like:
  ##  cell.line, treatment, compound, dose, strain, antibody
  
  #print(nrow(sampleData))
  #print(length(unique(sampleData$acc_dataset)))
  sampleData_patient = sampleData

  
  ## get columns which match to the patterns
  #unique(sampleData_patient$TRANSFECTION)
  #unique(sampleData_patient$acc_dataset[!is.na(sampleData_patient$TRANSFECTION)])
  #colnames(sampleData_patient)[grep(colnames(sampleData_patient), pattern = "cell.*line.*", ignore.case = TRUE)]
  
  whichColsCellline = grep(colnames(sampleData_patient), pattern = "cell.*line", ignore.case = TRUE)
  whichColsTreatment = grep(colnames(sampleData_patient), pattern ="treat", ignore.case = TRUE) # treated.with and treatment
  whichColsCompound = grep(colnames(sampleData_patient), pattern ="compound", ignore.case = TRUE)
  whichColsDose = grep(colnames(sampleData_patient), pattern ="dose", ignore.case = TRUE)
  whichColsStrain = grep(colnames(sampleData_patient), pattern ="^strain", ignore.case = TRUE) # ^ to prevent from detecting "immunohistochemistry.straining.for.estrogen.receptor.positive.1"
  whichColsAntibody = grep(colnames(sampleData_patient), pattern ="^antibody", ignore.case = TRUE) # ^ to prevent from detecting  "immunohistochemistry.staining.with.ki67.antibody.positive.1"
  whichColsStimulus = grep(colnames(sampleData_patient), pattern ="stimu", ignore.case = TRUE) #stimulus and stimulated
  whichColsDecitabine = grep(colnames(sampleData_patient), pattern ="decitabine", ignore.case = TRUE) 
  #whichColsTherapy = grep(colnames(sampleData_patient), pattern ="therapy", ignore.case = TRUE)
  #whichColsExposure = grep(colnames(sampleData_patient), pattern ="^expos", ignore.case = TRUE) #exposed.to and exposure; ^ to prevent from detecting "Asbestos.Exposure"
  whichColsAllergen = grep(colnames(sampleData_patient), pattern ="allergen", ignore.case = TRUE)
  whichColsActivation = grep(colnames(sampleData_patient), pattern ="activa", ignore.case = TRUE) #activated and activation
  whichColsRespond = grep(colnames(sampleData_patient), pattern ="^respon", ignore.case = TRUE) #response and responder.status
  whichColsTransfection = grep(colnames(sampleData_patient), pattern ="transfect", ignore.case = TRUE) #transfection and transfected.with
  whichColsTrasnfection = grep(colnames(sampleData_patient), pattern ="trasnfect", ignore.case = TRUE) #trasnfection.protocol
  whichColsCulture = grep(colnames(sampleData_patient), pattern ="culture.", ignore.case = TRUE) #culture.condition and culture.substrate
  whichColsChlorambucil = grep(colnames(sampleData_patient), pattern ="chlorambucil", ignore.case = TRUE)
  whichColsIntervention = grep(colnames(sampleData_patient), pattern ="interven", ignore.case = TRUE) #intervention and interventional, intervened
  whichColsDoxycycline = grep(colnames(sampleData_patient), pattern ="doxycycline", ignore.case = TRUE)
  whichColsBreed = grep(colnames(sampleData_patient), pattern ="breed", ignore.case = TRUE)
  whichColsLineage = grep(colnames(sampleData_patient), pattern ="lineage", ignore.case = TRUE)
  whichColsPrednisolon = grep(colnames(sampleData_patient), pattern ="prednisolon", ignore.case = TRUE)
  whichColsTransduction = grep(colnames(sampleData_patient), pattern ="transduc", ignore.case = TRUE) #transduction and transduced
  whichColsVaccine = grep(colnames(sampleData_patient), pattern ="vaccin", ignore.case = TRUE) #vaccine and vaccination
  whichColsVector = grep(colnames(sampleData_patient), pattern ="vector", ignore.case = TRUE)
  
  removeCols = c(whichColsCellline, whichColsTreatment,    whichColsCompound,     whichColsDose,         whichColsStrain,
                 whichColsAntibody, whichColsStimulus,     whichColsDecitabine,   #whichColsTherapy,     whichColsExposure,
                 whichColsAllergen, whichColsActivation,   whichColsRespond,      whichColsTransfection, whichColsTrasnfection,
                 whichColsCulture,  whichColsChlorambucil, whichColsIntervention, whichColsDoxycycline,  whichColsBreed,
                 whichColsLineage,  whichColsPrednisolon,  whichColsTransduction, whichColsVaccine,      whichColsVector)
  #removeCols_string = colnames(sampleData_patient)[removeCols]
  
  ## remove experiments which are interventional
  removeExperiments = unique(sampleData_patient$acc_dataset[which((rowSums(!is.na(sampleData_patient[, removeCols])))>0)])
  removedExperiments = c(removedExperiments, removeExperiments)
  ## 1. check which samples have a value in one of the columns that indicate intervention ( code: !is.na(sampleData_patient[, removeCols]) )
  ## 2. check if there is at least 1 (>0) column that has a value ( code:e rowSums("1.") > 0 ) -> those samples need to be removed 
  ## 3. get the index of the samples of 2. ( code: which("2.") )
  ## 4. select the dataset accessions for those indices ( code: sampleData_patient$acc_dataset["3."] )
  ## 5. keep the dataset accession only once ( code: unique("4.") )
  
  sampleData_patient = subset(sampleData_patient, subset = !(acc_dataset%in%removeExperiments), select = -removeCols)
  
  ## those columns which don't have any entry now, are getting obsolete
  obsoleteColumns = which(colSums(!(is.na(sampleData_patient)))==0)
  sampleData_patient = subset(sampleData_patient, select = -obsoleteColumns)
  
  #colnames(sampleData_patient)
  
  sampleDataALL_patient = dplyr::bind_rows(sampleDataALL_patient, sampleData_patient)
  #removeColsALL_string = unique(c(removeColsALL_string, colnames(sampleData_patient)[removeCols], colnames(sampleData_patient)[obsoleteColumns]))
  rm(list=setdiff(ls(), c("sampleDataALL_patient", "removeColsALL_string", "i", "removedExperiments")))
}

write(removeColsALL_string, ncolumns = 1, file = paste0("removedColumns.csv")) ##### 2020-10-08
save(sampleDataALL_patient, file=paste0("sampleDataALL_patient.Rdata")) ##### 2020-10-08





######---------get disease information--------------
rm(list=ls())
load(file="sampleDataALL_patient.Rdata")

sampleDataALL_disease = sampleDataALL_patient

whichColsDisease = grep(colnames(sampleDataALL_disease), pattern = "disease", ignore.case = TRUE)

#colnames(sampleDataALL_patient)[grep(colnames(sampleDataALL_patient), pattern = "disease", ignore.case = TRUE)]
#unique(sampleData_Onlydisease$acc_dataset[!is.na(sampleData_Onlydisease$diseased.free.survival..1..disease.)])
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
usedColsALL = colnames(sampleDataALL_disease)[keepCols]

keepExperiments = unique(sampleDataALL_disease$acc_dataset[which((rowSums(!is.na(sampleDataALL_disease[, keepCols])))>0)])
## 1. check if samples have a value in any of the disease columns ( code: !is.na(sampleData_patient[, keepCols]) )
## 2. check if there is at least 1 (>0) column that has a value ( code: rowSums("1.") > 0 ) -> those samples need to be kept 
## 3. get the index of the samples of 2. ( code: which("2.") )
## 4. select the dataset accessions for those indices ( code: sampleData_patient$acc_dataset["3."] )
## 5. keep the dataset accession only once ( code: unique("4.") )

sampleDataALL_disease= subset(sampleDataALL_disease, subset = (acc_dataset%in%keepExperiments))

## those columns which don't have any entry now, are getting obsolete
obsoleteColumns = which(colSums(!(is.na(sampleDataALL_disease)))==0)
obsoleteColumnsALL = colnames(sampleDataALL_disease)[obsoleteColumns]
sampleDataALL_disease = subset(sampleDataALL_disease, select = -obsoleteColumns)


sampleData_Onlydisease = subset(sampleDataALL_patient, subset = (acc_dataset%in%keepExperiments), select = c(acc_sample, acc_dataset, keepCols))

save(sampleDataALL_disease, file="sampleDataALL_disease.Rdata") ##### 2020-10-08
save(sampleData_Onlydisease, file="sampleData_Onlydisease.Rdata") ##### 2020-10-08



####---add title of experiment---
sampleData_Onlydisease_withTitle = sampleData_Onlydisease

acc_list = unique(sampleData_Onlydisease$acc_dataset)
pb <- txtProgressBar(min = 0, max = length(acc_list)-1, style = 3)

suppressMessages(
  for (i in 1:length(acc_list)) {
    gse = acc_list[i]
    #print(gse)
    
    res = fromJSON(paste0("https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/", gse))
    expTitle = res$experiments$experiment$name
    
    sampleData_Onlydisease_withTitle$expTitle[sampleData_Onlydisease_withTitle$acc_dataset==gse] = expTitle
    
    setTxtProgressBar(pb, i-1)
    
  }
)
close(pb)
# ******* table writing termsmapped to entology. terms curated manually
save(sampleData_Onlydisease_withTitle, file="sampleData_Onlydisease_withTitle.Rdata") ##### 2020-10-08
write.table(sampleData_Onlydisease_withTitle, file="sampleData_Onlydisease_withTitle.tsv", 
            sep="\t", quote=F, row.names = F) ##### 2020-10-08

###---get term that need to be mapped
#termsForMapping = unique(unlist(sampleData_Onlydisease_withTitle[,-c(1:2)]))
#save(termsForMapping, file=paste0(Sys.Date(), "_termsForMapping.Rdata")) ##### 2020-10-08
#write(termsForMapping, file=paste0(Sys.Date(), "_termsForMapping.csv"), ncolumns = 1) ##### 2020-10-08

# disease is added with different script: diseaseCleanUp.R


### remove specific datasets
rm(list = ls())
load(file="sampleDataALL_disease.Rdata")

#remove dataset with agent
expAgent = unique(sampleDataALL_disease$acc_dataset[!is.na(sampleDataALL_disease$agent)])
expAgent = c(expAgent, unique(sampleDataALL_disease$acc_dataset[!is.na(sampleDataALL_disease$AGENT)]))
sampleDataALL_disease = subset(sampleDataALL_disease, subset = !(acc_dataset%in%expAgent))

## those columns which don't have any entry now, are getting obsolete
obsoleteColumnsALL = c()
obsoleteColumns = which(colSums(!(is.na(sampleDataALL_disease)))==0)
obsoleteColumnsALL = union(colnames(sampleDataALL_disease)[obsoleteColumns], obsoleteColumnsALL)
sampleDataALL_disease = subset(sampleDataALL_disease, select = -obsoleteColumns)

#remove dataset with organism != human
colsOrganismString = c("organism", "organism.1", "Organism", "Organism.1", "ORGANISM")
whichColsOrganism = which(colnames(sampleDataALL_disease)%in%colsOrganismString)
usedColsALL = c()
usedColsALL = union(colnames(sampleDataALL_disease)[whichColsOrganism], usedColsALL)

## normalize terms
organismALL = subset(sampleDataALL_disease, select = c(acc_sample, acc_dataset, whichColsOrganism))

organismALL[ organismALL == "Homo sapiens"] = "homo sapiens"
organismALL[ organismALL == "Homo Sapiens"] = "homo sapiens"
organismALL[ organismALL == "Homo sapiens "] = "homo sapiens"

organismALL$organismClean = NA
for (i in 1:NROW(organismALL)) {
  colsNotNA = colnames(organismALL)[!is.na(organismALL[i,])]
  colOrganism = setdiff(colsNotNA, c("acc_sample", "acc_dataset"))
  valsOrganism = unique(as.vector(as.matrix(organismALL[i, colOrganism])))
  if (length(valsOrganism)==0) {
    next
  } else if (length(valsOrganism)==1) {
    organismALL[i, "organismClean"] = tolower(valsOrganism)
  } else {
    print(i)
    break
  }
}

organismNotHuman_Samples = which(organismALL$organismClean != "homo sapiens")
organismALL = organismALL[-c(organismNotHuman_Samples),]
organismCleanDataframe = subset(organismALL, select = c(acc_sample, acc_dataset, organismClean))
save(organismCleanDataframe, file=paste0("organismClean.Rdata"))  ##### 2020-10-08

sampleDataALL_disease = sampleDataALL_disease[-c(organismNotHuman_Samples),]

## those columns which don't have any entry now, are getting obsolete
obsoleteColumns = which(colSums(!(is.na(sampleDataALL_disease)))==0)
obsoleteColumnsALL = union(colnames(sampleDataALL_disease)[obsoleteColumns], obsoleteColumnsALL)
sampleDataALL_disease = subset(sampleDataALL_disease, select = -obsoleteColumns)

save(sampleDataALL_disease, file="sampleDataALL_disease.Rdata")  ##### 2020-10-08
save(obsoleteColumnsALL, file="obsoleteCols.Rdata")  ##### 2020-10-08
save(usedColsALL, file="usedCols.Rdata")  ##### 2020-10-08

######---------map metadata--------------
rm(list = ls())

## add unique identifier
sampleDataALL_disease = data.frame(sampleID = 1:NROW(sampleDataALL_disease), sampleDataALL_disease)
save(sampleDataALL_disease, file="sampleDataALL_disease.Rdata")  ##### 2020-10-08

rm(list = ls())
load(file="sampleDataALL_disease.Rdata")
load(file="obsoleteCols.Rdata")
load(file="usedCols.Rdata")

## age
# get age columns and remove non related columns or uninformative age columns
whichColsAge = grep(colnames(sampleDataALL_disease), pattern = "age", ignore.case = TRUE)
whichColsStage = grep(colnames(sampleDataALL_disease), pattern = "stage", ignore.case = TRUE)
whichColsStorage = grep(colnames(sampleDataALL_disease), pattern = "storage", ignore.case = TRUE)
whichColsPassage = grep(colnames(sampleDataALL_disease), pattern = "passage", ignore.case = TRUE)
whichColsPercentage = grep(colnames(sampleDataALL_disease), pattern = "percentage", ignore.case = TRUE)

whichColsNotAge = c(whichColsStage, whichColsStorage, whichColsPassage, whichColsPercentage)
whichColsAge = setdiff(whichColsAge, whichColsNotAge)

usedColsALL = union(colnames(sampleDataALL_disease)[whichColsAge], usedColsALL)

#colnames(sampleDataALL_disease)[whichColsAge]
#unique(sampleDataALL_disease$acc_dataset[!is.na(sampleDataALL_disease$specimen.with.known.storage.state)])
age_varGroup1 = c("age", "age.1", "Age", "Age.1", "age.years", "age.years.1", "age.yrs", "AGE", "var_age", "AGE.YEARS",
                  "AGE.YRS", "age.in.years", "AGE.IN.YEARS", "age.yr", "AGE.YR", "age.y", "AGE.Y", "age..years.", "AGE..YEARS.",
                  "Age.years.")
age_varGroup2 = c("age.group", "age.group.1")
age_varGroup3 = c("age.at.diagnosis", "age.at.diagnosis.1", "age.at.diagnosis.years", "Age.At.Diagnosis", "age.at.death", 
                  "AGE.AT.DEATH", "AGE.AT.DIAGNOSIS.YEARS", "age.at.biopsy", "AGE.AT.BIOPSY", "age.at.surgery", "AGE.AT.SURGERY",
                  "Age.at.Diagnosis", "age.at.diagnosis..years.", "AGE.AT.DIAGNOSIS..YEARS.", "age.at.examination", 
                  "age.at.onset", "AGE.AT.EXAMINATION", "AGE.AT.ONSET", "age.on.study", "AGE.ON.STUDY", 
                  "age.at.diagonosis..years.", "AGE.AT.DIAGONOSIS..YEARS.")
age_varGroup4 = c("gestational.age", "Gestational.Age", "Gestation.Age", "gestational.age.1")

#colnames(ageALL)[!(colnames(ageALL)%in%c(age_varGroup1, age_varGroup2, age_varGroup3, age_varGroup4))]

#create clean age column
ageALL = subset(sampleDataALL_disease, select = c(sampleID, acc_sample, acc_dataset, whichColsAge))
ageALL = subset(ageALL, select = -c(AgeStatus, AgeStatus.1))

ageALL[ ageALL == "NA" ] = NA
ageALL[ ageALL == "not specified"] = NA
ageALL$ageClean = NA
ageALL$gestagenAgeClean = NA
ageALL$ageGroupClean = NA

for (i in 1:NROW(ageALL)) {
  # get which columns are not NA
  # i = 6945
  colsNotNA = colnames(ageALL)[!is.na(ageALL[i,])]
  colString = setdiff(colsNotNA, c("sampleID", "acc_sample", "acc_dataset"))
  if (length(colString)==0) {
    next
  } else {
    if (any(colString%in%age_varGroup2)) {
      colsNotGroup = setdiff(colString, age_varGroup2)
      #colsGroup = setdiff(colString, colsNotGroup)
      #valsGroup = unique(as.vector(as.matrix(ageALL[i, colsGroup])))
      if (length(colsNotGroup)>0) {
        vals = unique(as.vector(as.matrix(ageALL[i, colsNotGroup])))
        if (length(vals)==1) {
          ageALL[i, "ageGroupClean"] = vals
        } else {
          print(i)
          print(colsNotGroup)
          print(vals)
          print("Values in other than group columns are different.")
          break
        }
      } else {
        print(i)
        print("No other columns than group columns.")
        break
      }
    } else if (any(colString%in%age_varGroup4)) {
      colsNotGesta = setdiff(colString, age_varGroup4)
      colsGesta = setdiff(colString, colsNotGesta)
      valsGesta = unique(as.vector(as.matrix(ageALL[i, colsGesta])))
      if (length(colsNotGesta)>0) {
        vals = unique(as.vector(as.matrix(ageALL[i, colsNotGesta])))
        if (length(vals)==1) {
          if (vals==0) {
            ageALL[i, "ageClean"] = vals # = 0
            if (length(valsGesta)==1) {
              ageALL[i, "gestagenAgeClean"] = valsGesta
            } else {
              print(i)
              print("Gestagen column values are different.")
              break
            }
          } else {
            print(i)
            print("Age column not 0 when gestag age given.")
            break
          }
        } else {
          print(i)
          print(colsNotGroup)
          print(vals)
          print("Values in other than gesta columns are different.")
          break
        }
      } else if (length(valsGesta)==1) {
        ageALL[i, "ageClean"] = 0
        ageALL[i, "gestagenAgeClean"] = valsGesta
      } else {
        print(i)
        print("Some problem with gestagen column.")
        break
      }
    } else if (any(colString%in%age_varGroup3)) {
      colsNotAt = setdiff(colString, age_varGroup3)
      colsAt = setdiff(colString, colsNotAt)
      valsAt = unique(as.vector(as.matrix(ageALL[i, colsAt])))
      if (length(colsNotAt)>0) {
        break
      } else if (length(valsAt)==1) {
        ageALL[i, "ageClean"] = valsAt
      } else if (any(grep("examination", colsAt))) {
        ageALL[i, "ageClean"] = ageALL[i, "AGE.AT.EXAMINATION"]
      } else {
        print(i)
        print(colsAt)
        print("Different values for age at sth (not examination).")
        break
      }
    } else if (any(colString%in%age_varGroup1)) {
      vals = unique(as.vector(as.matrix(ageALL[i, colString])))
      vals = unique(trimws(unlist(strsplit(vals, ";"))))
      if (length(vals)==1) {
        ageALL[i, "ageClean"] = vals
      } else {
        checkYearsString = any(grepl(x = vals, pattern = " years"))
        if (checkYearsString) {
          vals = unlist(strsplit(vals, " years"))
          vals = unique(vals)
          if (length(vals)==1) {
            ageALL[i, "ageClean"] = vals
          } else {
            print(i)
            print(colString)
            print(vals)
            print(sampleDataALL_disease[i,"acc_dataset"])
            print("Values in age columns are different after check for years string.")
            break
          }
        } else {
          print(i)
          print(colString)
          print(vals)
          print(sampleDataALL_disease[i,"acc_dataset"])
          print("Values in age columns are different.")
          break
        }
      }
    } else {
      print(i)
      break
    }
  }
}

ageCleanDataframe = subset(ageALL, select = c(sampleID, acc_sample, acc_dataset, ageClean, gestagenAgeClean, ageGroupClean))
save(ageCleanDataframe, file="ageClean.Rdata")  ##### 2020-10-08
save(obsoleteColumnsALL, file="obsoleteCols.Rdata")  ##### 2020-10-08
save(usedColsALL, file="usedCols.Rdata")  ##### 2020-10-08


## gender
# get gender/sex columns
whichColsGender = grep(colnames(sampleDataALL_disease), pattern = "gender", ignore.case = TRUE)
whichColsSex = grep(colnames(sampleDataALL_disease), pattern = "sex", ignore.case = TRUE)
whichColsGender = union(whichColsGender, whichColsSex)
#colnames(sampleDataALL_disease)[whichColsGender]
#unique(sampleDataALL_disease$acc_dataset[!is.na(sampleDataALL_disease$female.gender)])
#unique(sampleDataALL_disease$var_sex)

usedColsALL = union(colnames(sampleDataALL_disease)[whichColsGender], usedColsALL)

genderALL = subset(sampleDataALL_disease, select = c(sampleID, acc_sample, acc_dataset, whichColsGender))

#create clean gender column
## set all unavailable data to real NA's
genderALL[ genderALL == "NA" ] = NA
genderALL[ genderALL == "not available" ] = NA
genderALL[ genderALL == "unknown" ] = NA
genderALL[ genderALL == "unknown sex" ] = NA
genderALL[ genderALL == "unknown_sex" ] = NA
genderALL[ genderALL == "Not available" ] = NA
genderALL[ genderALL == "not specified" ] = NA
genderALL[ genderALL == "n/a" ] = NA

## normalize terms
normGender = c("female", "male")

genderALL[ genderALL == "Female" ] = "female"
genderALL[ genderALL == "Male" ] = "male"

genderALL[ genderALL == "F" ] = "female"
genderALL[ genderALL == "M" ] = "male"

genderALL[ genderALL == "female " ] = "female"
genderALL[ genderALL == "male " ] = "male"

genderALL[ genderALL == "Female" ] = "female"
genderALL[ genderALL == "Male " ] = "male"

genderALL[ genderALL == "femaleale" ] = "female"

## get clean column
genderALL$genderClean = NA
for (i in 1:NROW(genderALL)) {
  colsNotNA = colnames(genderALL)[!is.na(genderALL[i,])]
  colGender = setdiff(colsNotNA, c("sampleID", "acc_sample", "acc_dataset"))
  valsGender = unique(as.vector(as.matrix(genderALL[i, colGender])))
  if (genderALL$acc_dataset[i] == "E-GEOD-38627") {
    if (length(valsGender)==0) {
      genderALL[i, "genderClean"] = "male"
    } else if (valsGender==1) {
      genderALL[i, "genderClean"] = "female"
    } else {
      break
    }
  } else {
    if (length(valsGender)==0) {
      next
    } else if (length(valsGender)==1) {
      genderALL[i, "genderClean"] = valsGender
    } else {
      print(i)
      break
    }
  }
}

genderCleanDataframe = subset(genderALL, select = c(sampleID, acc_sample, acc_dataset, genderClean))
save(genderCleanDataframe, file="genderClean.Rdata")  ##### 2020-10-08
save(obsoleteColumnsALL, file="obsoleteCols.Rdata")  ##### 2020-10-08
save(usedColsALL, file="usedCols.Rdata")  ##### 2020-10-08


## ethnecity
# get ethnic / race columns
whichColsEthnic = grep(colnames(sampleDataALL_disease), pattern = "ethn", ignore.case = TRUE)
whichColsRace = grep(colnames(sampleDataALL_disease), pattern = "race", ignore.case = TRUE)
whichColsEthnic = union(whichColsEthnic, whichColsRace)
#colnames(sampleDataALL_disease)[whichColsEthnic]
#unique(sampleDataALL_disease$acc_dataset[!is.na(sampleDataALL_disease$ethnic.group)])
#unique(sampleDataALL_disease$ethnic.group)
#unique(as.vector(as.matrix(subset(sampleDataALL_disease, select=whichColsEthnic))))
#unique(as.vector(as.matrix(subset(ethnicALL, select=-c(1,2)))))

usedColsALL = union(colnames(sampleDataALL_disease)[whichColsEthnic], usedColsALL)

ethnicALL = subset(sampleDataALL_disease, select = c(sampleID, acc_sample, acc_dataset, whichColsEthnic))

#create clean gender column
## set all unavailable data to real NA's
ethnicALL[ ethnicALL == "NA" ] = NA
ethnicALL[ ethnicALL == "Unknown" ] = NA
ethnicALL[ ethnicALL == "?NA" ] = NA
ethnicALL[ ethnicALL == "not available" ] = NA
ethnicALL[ ethnicALL == "unknown" ] = NA
ethnicALL[ ethnicALL == "race unknown" ] = NA
ethnicALL[ ethnicALL == "not specified" ] = NA

## normalize terms
#normEthnic = ??
ethnicALL[ ethnicALL == "W" ] = "white"
ethnicALL[ ethnicALL == "B" ] = "black"
ethnicALL[ ethnicALL == "Afr" ] = "african"
ethnicALL[ ethnicALL == "Eur" ] = "european"
ethnicALL[ ethnicALL == "other (not caucasian, black, or hispanic)" ] = "other"
ethnicALL[ ethnicALL == "Arab" ] = "arabian"
ethnicALL[ ethnicALL == "caucasuian" ] = "caucasian"

## get clean column
ethnicALL$ethnicClean = NA
for (i in 1:NROW(ethnicALL)) {
  colsNotNA = colnames(ethnicALL)[!is.na(ethnicALL[i,])]
  colEthnic = setdiff(colsNotNA, c("sampleID", "acc_sample", "acc_dataset"))
  valsEthnic = unique(as.vector(as.matrix(ethnicALL[i, colEthnic])))
  if (length(valsEthnic)==0) {
    next
  } else if (length(valsEthnic)==1) {
    ethnicALL[i, "ethnicClean"] = tolower(valsEthnic)
  } else {
    print(i)
    break
  }
}

#unique(ethnicALL$ethnicClean)

ethnicCleanDataframe = subset(ethnicALL, select = c(sampleID, acc_sample, acc_dataset, ethnicClean))
save(ethnicCleanDataframe, file="ethnicClean.Rdata")  ##### 2020-10-08
save(obsoleteColumnsALL, file="obsoleteCols.Rdata")  ##### 2020-10-08
save(usedColsALL, file="usedCols.Rdata")  ##### 2020-10-08


## cell type
# get cell type columns
whichColsCelltype = grep(colnames(sampleDataALL_disease), pattern = "cell.*ty", ignore.case = TRUE)
whichColsProperty = grep(colnames(sampleDataALL_disease), pattern = "property", ignore.case = TRUE)
whichColsCelltype = setdiff(whichColsCelltype, whichColsProperty)
#test = unique(as.vector(as.matrix(subset(celltypeALL, select=-c(1,2)))))
#write(test, ncolumns = 1, file="../../../../../Desktop/colsTemp.txt")

usedColsALL = union(colnames(sampleDataALL_disease)[whichColsCelltype], usedColsALL)

celltypeALL = subset(sampleDataALL_disease, select = c(sampleID, acc_sample, acc_dataset, whichColsCelltype))

#create clean celltype column
## set all unavailable data to real NA's
celltypeALL[ celltypeALL == "NA" ] = NA
celltypeALL[ celltypeALL == "Not applicable" ] = NA
celltypeALL[ celltypeALL == "not specified" ] = NA

## get clean column
celltypeALL$celltypeClean = NA
celltypeALL$celltypeClean2 = NA
celltypeALL$celltypeClean3 = NA
for (i in 1:NROW(celltypeALL)) {
  colsNotNA = colnames(celltypeALL)[!is.na(celltypeALL[i,])]
  colCelltype = setdiff(colsNotNA, c("sampleID", "acc_sample", "acc_dataset"))
  valsCelltype = unique(as.vector(as.matrix(celltypeALL[i, colCelltype])))
  if (length(valsCelltype)==0) {
    next
  } else {
    valsCelltype = unique(tolower(trimws(unlist(strsplit(valsCelltype, ";")))))
    if (length(valsCelltype)==1) {
      celltypeALL[i, "celltypeClean"] = tolower(valsCelltype)
    } else if (length(valsCelltype==3)) {
      celltypeALL[i, "celltypeClean"] = tolower(valsCelltype[1])
      celltypeALL[i, "celltypeClean2"] = tolower(valsCelltype[2])
      celltypeALL[i, "celltypeClean3"] = tolower(valsCelltype[3])
    } else {
      break
    }
  }
}

#unique(celltypeALL$celltypeClean)
#test = subset(celltypeALL, acc_dataset==celltypeALL$acc_dataset[i])

celltypeCleanDataframe = subset(celltypeALL, select = c(sampleID, acc_sample, acc_dataset, celltypeClean, celltypeClean2, celltypeClean3))
save(celltypeCleanDataframe, file="celltypeClean.Rdata")  ##### 2020-10-08
save(obsoleteColumnsALL, file= "obsoleteCols.Rdata")  ##### 2020-10-08
save(usedColsALL, file="usedCols.Rdata")  ##### 2020-10-08


## organism part
# get organsim part (and organism)
whichColsOrganismPart = grep(colnames(sampleDataALL_disease), pattern = "organism", ignore.case = TRUE)
colsOrganismString = c("organism", "organism.1", "Organism", "Organism.1", "ORGANISM")
excludeOrganism = which(colnames(sampleDataALL_disease)%in%colsOrganismString)
whichColsOrganismPart = setdiff(whichColsOrganismPart, excludeOrganism)
whichColsOrganismStatus = grep(colnames(sampleDataALL_disease), pattern = "status", ignore.case = TRUE)
whichColsOrganismPart = setdiff(whichColsOrganismPart, whichColsOrganismStatus)

#colnames(sampleDataALL_disease)[whichColsEthnic]
#unique(sampleDataALL_disease$acc_dataset[!is.na(sampleDataALL_disease$ethnic.group)])
#unique(sampleDataALL_disease$ethnic.group)
#unique(as.vector(as.matrix(subset(sampleDataALL_disease, select=whichColsEthnic))))
#unique(as.vector(as.matrix(subset(ethnicALL, select=-c(1,2)))))

usedColsALL = union(colnames(sampleDataALL_disease)[whichColsOrganismPart], usedColsALL)

organismPartALL = subset(sampleDataALL_disease, select = c(sampleID, acc_sample, acc_dataset, whichColsOrganismPart))
#test = unique(as.vector(as.matrix(subset(organismPartALL, select=-c(1,2)))))
#write(test, ncolumns = 1, file="../../../../../Desktop/colsTemp.txt")

#create clean organismPart column
## set all unavailable data to real NA's
organismPartALL[ organismPartALL == "NA" ] = NA
organismPartALL[ organismPartALL == "--" ] = NA
organismPartALL[ organismPartALL == "Other" ] = NA
organismPartALL[ organismPartALL == "other" ] = NA
organismPartALL[ organismPartALL == "not specified" ] = NA
organismPartALL[ organismPartALL == "n/a" ] = NA
organismPartALL[ organismPartALL == "no data" ] = NA

## get clean column
#datasetsSplitStr = c("E-TABM-389")
#datasetsNormedVals = c("E-GEOD-11755")

organismPartALL$organismPartClean = NA
organismPartALL$organismPartClean2 = NA
organismPartALL$organismPartClean3 = NA
organismPartALL$organismPartClean4 = NA
for (i in 1:NROW(organismPartALL)) {
  colsNotNA = colnames(organismPartALL)[!is.na(organismPartALL[i,])]
  colOrganismPart = setdiff(colsNotNA, c("sampleID", "acc_sample", "acc_dataset"))
  valsOrganismPart = unique(as.vector(as.matrix(organismPartALL[i, colOrganismPart])))
  valsOrganismPart = unique(trimws(tolower(valsOrganismPart)))
  if (length(valsOrganismPart)==0) {
    next
  } else {
    valsOrganismPart = unique(trimws(unlist(strsplit(valsOrganismPart, ";"))))
    if (length(valsOrganismPart)==1) {
      organismPartALL[i, "organismPartClean"] = tolower(valsOrganismPart)
    } else if (length(valsOrganismPart)==2) {
      organismPartALL[i, "organismPartClean"] = tolower(valsOrganismPart[1])
      organismPartALL[i, "organismPartClean2"] = tolower(valsOrganismPart[2])
    } else if (length(valsOrganismPart)==4) {
      organismPartALL[i, "organismPartClean"] = tolower(valsOrganismPart[1])
      organismPartALL[i, "organismPartClean2"] = tolower(valsOrganismPart[2])
      organismPartALL[i, "organismPartClean3"] = tolower(valsOrganismPart[3])
      organismPartALL[i, "organismPartClean4"] = tolower(valsOrganismPart[4])
    } else {
      break
    }
  }
}

#unique(organismPartALL$organismPartClean)

organismPartCleanDataframe = subset(organismPartALL, select = c(sampleID, acc_sample, acc_dataset, organismPartClean, organismPartClean2, organismPartClean3, organismPartClean4))
save(organismPartCleanDataframe, file="organismPartClean.Rdata")  ##### 2020-10-08
save(obsoleteColumnsALL, file="obsoleteCols.Rdata")  ##### 2020-10-08
save(usedColsALL, file="usedCols.Rdata") ##### 2020-10-08


## batch
# get batch
whichColsBatch = grep(colnames(sampleDataALL_disease), pattern = "batch", ignore.case = TRUE)

#colnames(sampleDataALL_disease)[whichColsBatch]
#unique(sampleDataALL_disease$acc_dataset[!is.na(sampleDataALL_disease$ethnic.group)])
#unique(sampleDataALL_disease$ethnic.group)
#unique(as.vector(as.matrix(subset(sampleDataALL_disease, select=whichColsEthnic))))
#unique(as.vector(as.matrix(subset(ethnicALL, select=-c(1,2)))))

usedColsALL = union(colnames(sampleDataALL_disease)[whichColsBatch], usedColsALL)

batchALL = subset(sampleDataALL_disease, select = c(sampleID, acc_sample, acc_dataset, whichColsBatch))

#create clean batch column
## set all unavailable data to real NA's
batchALL[ batchALL == "not specified" ] = NA

batchALL$batchClean = NA
for (i in 1:NROW(batchALL)) {
  colsNotNA = colnames(batchALL)[!is.na(batchALL[i,])]
  colBatch = setdiff(colsNotNA, c("sampleID", "acc_sample", "acc_dataset"))
  valsBatch = unique(as.vector(as.matrix(batchALL[i, colBatch])))
  valsBatch = unique(trimws(tolower(valsBatch)))
  if (length(valsBatch)==0) {
    next
  } else if (length(valsBatch)==1) {
    batchALL[i, "batchClean"] = tolower(valsBatch)
  } else {
    print(i)
    break
  }
}

batchCleanDataframe = subset(batchALL, select = c(sampleID, acc_sample, acc_dataset, batchClean))
save(batchCleanDataframe, file="batchClean.Rdata")  ##### 2020-10-08
save(obsoleteColumnsALL, file="obsoleteCols.Rdata")  ##### 2020-10-08
save(usedColsALL, file="usedCols.Rdata")  ##### 2020-10-08



## tissue
# get tissue (type)
whichColsTissue = grep(colnames(sampleDataALL_disease), pattern = "tissue", ignore.case = TRUE)
colsCollectionString = c("tissue.collection")
excludeTissue = which(colnames(sampleDataALL_disease)%in%colsCollectionString)
whichColsTissue = setdiff(whichColsTissue, excludeTissue)

#colnames(sampleDataALL_disease)[whichColsEthnic]
#unique(sampleDataALL_disease$acc_dataset[!is.na(sampleDataALL_disease$ethnic.group)])
#unique(sampleDataALL_disease$ethnic.group)
#unique(as.vector(as.matrix(subset(sampleDataALL_disease, select=whichColsEthnic))))
#unique(as.vector(as.matrix(subset(ethnicALL, select=-c(1,2)))))

usedColsALL = union(colnames(sampleDataALL_disease)[whichColsTissue], usedColsALL)

tissueALL = subset(sampleDataALL_disease, select = c(sampleID, acc_sample, acc_dataset, whichColsTissue))

#create clean tissue column
## set all unavailable data to real NA's
tissueALL[ tissueALL == "NOT SPECIFIED" ] = NA

tissueALL$tissueClean = NA
tissueALL$tissueClean2 = NA
for (i in 1:NROW(tissueALL)) {
  colsNotNA = colnames(tissueALL)[!is.na(tissueALL[i,])]
  colTissue = setdiff(colsNotNA, c("sampleID", "acc_sample", "acc_dataset"))
  valsTissue = unique(as.vector(as.matrix(tissueALL[i, colTissue])))
  valsTissue = unique(trimws(tolower(valsTissue)))
  if (length(valsTissue)==0) {
    next
  } else {
    valsTissue = unique(trimws(unlist(strsplit(valsTissue, ";"))))
    if (length(valsTissue)==1) {
      tissueALL[i, "tissueClean"] = tolower(valsTissue)
    } else if (length(valsTissue)==2) {
      tissueALL[i, "tissueClean"] = tolower(valsTissue[1])
      tissueALL[i, "tissueClean2"] = tolower(valsTissue[2])
    } else {
      print(i)
      break
    }
  }
}

tissueCleanDataframe = subset(tissueALL, select = c(sampleID, acc_sample, acc_dataset, tissueClean, tissueClean2))
save(tissueCleanDataframe, file="tissueClean.Rdata")  ##### 2020-10-08
save(obsoleteColumnsALL, file="obsoleteCols.Rdata")  ##### 2020-10-08
save(usedColsALL, file="usedCols.Rdata")  ##### 2020-10-08


######## write all metadata together
rm(list=ls())
load(file="ageClean.Rdata")
load(file="genderClean.Rdata")
load(file="ethnicClean.Rdata")
load(file="celltypeClean.Rdata")
load(file="organismPartClean.Rdata")
load(file="batchClean.Rdata")
load(file="tissueClean.Rdata")

sampleDataClean = merge(ageCleanDataframe, genderCleanDataframe[,-c(2,3)], by = "sampleID")
sampleDataClean = merge(sampleDataClean, ethnicCleanDataframe[,-c(2,3)], by = "sampleID")
sampleDataClean = merge(sampleDataClean, celltypeCleanDataframe[,-c(2,3)], by = "sampleID")
sampleDataClean = merge(sampleDataClean, organismPartCleanDataframe[,-c(2,3)], by = "sampleID")
sampleDataClean = merge(sampleDataClean, batchCleanDataframe[,-c(2,3)], by = "sampleID")
sampleDataClean = merge(sampleDataClean, tissueCleanDataframe[,-c(2,3)], by = "sampleID")

save(sampleDataClean, file="sampleDataClean.Rdata")   ##### 2020-10-08

#############################################################################################################
# 
# # check which columns are not used until now
# load(file="2020-09-07_sampleDataALL_disease.Rdata")
# load(file="2020-09-07_obsoleteCols.Rdata")
# load(file="2020-09-07_usedCols.Rdata")
# 
# setdiff(colnames(sampleDataALL_disease), union(obsoleteColumnsALL, usedColsALL))[201:300]
# 
# 
# colnames(sampleDataALL_disease)[grep(colnames(sampleDataALL_disease), pattern = "cell.*ty", ignore.case = TRUE)]
# unique(sampleDataALL_disease$smoking)
