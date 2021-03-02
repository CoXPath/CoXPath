library("ArrayExpress")
library("pd.hg.u133.plus.2")
library("oligo")

## load all data
load(file="metadata_final.Rdata")

datasets = unique(metadataFinal$acc_dataset)

datasetsCausingProblems = readLines("datasetsCausingProblems.txt")
pb <- txtProgressBar(min = 0, max = length(datasets)-1, style = 3)

suppressMessages(
  for (i in 1:length(datasets)) {
    if (datasets[i]%in%datasetsCausingProblems) {
      next
    }
    fname_exprsSet = paste0("~/export", datasets[i], ".Rdata")
    if (file.exists(fname_exprsSet)) {
      next
    }
    rawset = ArrayExpress(datasets[i])
    if (!is.null(rawset)) {
      if (length(rawset)>1) {
        indx = which(names(rawset)=="A-AFFY-44")
        rawData = rawset[[indx]]
        AEsetnorm = oligo::rma(rawData)
        save(AEsetnorm, file = fname_exprsSet)
      } else {
        AEsetnorm = oligo::rma(rawset)
        save(AEsetnorm, file = fname_exprsSet)
      }
    } else {
      datasetsCausingProblems = c(datasetsCausingProblems, datasets[i])
    }
    setTxtProgressBar(pb, i-1)
    rm(list=setdiff(ls(), c("datasets", "datasetsCausingProblems", "pb", "i")))
    filesToRemove = paste0(tempdir(), "/", list.files(tempdir()))
    filesToRemove = setdiff(filesToRemove, paste0(tempdir(), "/", "A-AFFY-44.adf.txt"))
    filesToRemove = setdiff(filesToRemove, '/tmp/RtmpVpoO8X/rs-graphics-1b408c82-c603-41ef-b5a3-26723f71f4eb')
    file.remove(filesToRemove)
    gc()
  }
)
close(pb)


write(datasetsCausingProblems, ncolumns = 1, file="datasetsCausingProblems.txt")
#write(datasetsCausingProblems, ncolumns = 1, file="C://Users/traschka/Documents/Projects/PhD/zml-project/data/datasetsCausingProblems.txt")

datasets[i]
i

#####------------- remove metadata of experiments that couldn't be loaded -----------------------
## load all data
rm(list=ls())
load(file="metadata_final.Rdata")

datasetsALL = unique(metadataFinal$acc_dataset)
#datasetsLoaded = list.files("Z:/NO BACKUP/exprsSets_rma/")
datasetsLoaded = list.files("~/export")
datasetsLoaded = unlist(strsplit(datasetsLoaded, ".Rdata"))

missingDatasets = setdiff(datasetsALL, datasetsLoaded)
datasetsCausingProblems = readLines("datasetsCausingProblems.txt")
datasetsToDelete = readLines("datasetsToDelete.txt")

missingDatasets = setdiff(missingDatasets, datasetsCausingProblems)
missingDatasets = setdiff(missingDatasets, datasetsToDelete)
length(missingDatasets)==0

metadataFinal = subset(metadataFinal, acc_dataset%in%datasetsLoaded)
save(metadataFinal, file="metadataFinal_afterDataLoading.RData")
#28543 samples, 476 experiments
