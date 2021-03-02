# path
setwd("/home/traschka/Projects/PhD/co-expression-vs-pathways/data/data_for_coexp_network_construction")

options(stringsAsFactors=FALSE)
rm(list=ls())
library("ArrayExpress")
library("sva")
library("jsonlite")
library("data.table")
library("dplyr")
library("tidyr")
library("genefilter")
library("hgu133plus2.db")

#copied from group definition script
groupIDs = c('1612','9119','9256','9538','1324','3069','3908','0070004','9952','2394',
             '0050745','2841','1040','8893','5520','0040085','1793','1909','5603','3083',
             '0050908','263','4450','9970','10652','0050866','9778','3910','169','0060074',
             '2377','8567','0050746','10283','5419','1115','657','0050749','0050750','5844',
             '3070','4074','3247','769','3382','1967','1596','9261','0050902','9074','289',
             '0050909','3234','7502','3312','1470','8778','9253','11335','3310','10223','9744',
             '0080199','1883','normal')

skippedIDs = c('263', 'normal')

for(ID in groupIDs) {

  print(paste0("Next group: ", ID))

  if (ID%in%skippedIDs) {
    next
  }

  #datasets = c("E-MTAB-3513", "E-MTAB-6844", "E-MTAB-5933", "E-MTAB-3516", "E-MTAB-5724", "E-MTAB-5279")
  datasets = readLines(paste0(ID, "/datasets.txt"))

  #load the correct metadata file (the one where you have also defined your per disease datasets from)
  metadataFinal = read.table(file = "../metadata_FINAL.tsv", sep = '\t', header = TRUE, quote = "")

  #load the metadata file that has only samples for a specific disease
  metadata = read.table(file = paste0(ID,"/metadata.tsv"), sep = '\t', header = TRUE, quote = "")

  source("../../scripts/getCelFile.R")

  datadir = paste0("./", ID, "/")

  if(file.exists(paste0(datadir, "data_annotated.RData"))) {
    next
  }

  ###### Download data
  datasetFiles = paste0(datadir, datasets, ".Rdata")


  if (!(all(file.exists(paste0(datadir, datasets, ".Rdata"))))) {
    print("Start downloading datasets!")

    pb <- txtProgressBar(min = 0, max = length(datasets)-1, style = 3)

    suppressMessages(
      for (i in 1:length(datasets)) {
        #file.exists(paste0(datadir, datasets[7], ".Rdata"))
        fname_exprsSet = paste0(datadir, datasets[i], ".Rdata")
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
        }
        setTxtProgressBar(pb, i-1)
      }
    )

    close(pb)
    print("Finished downloading datasets!")

  } else {
    print("Download was already done!")
  }

  datasets = datasets[file.exists(paste0(datadir, datasets, ".Rdata"))]
  metadataFinal = subset(metadataFinal, acc_dataset%in%datasets)
  metadata = subset(metadata, acc_dataset%in%datasets)


  ###### Merge data
  if (!file.exists(paste0(datadir, "exprsALL.RData"))) {

    getSampleIDs = function() {
      tryCatch({
        sampleIDs = sampleData_dataset$sampleID[sapply(sampleNames, FUN = grep, x = sampleData_dataset$acc_sample)]
        return(sampleIDs)
      }, error = function(msg) {
        sampleData_dataset$celFiles = getCelFileNamesFromAE(datasets[i])
        sampleNames = colnames(exprsData)
        sampleIDs = sampleData_dataset$sampleID[match(sampleNames,sampleData_dataset$celFiles)]
        return(sampleIDs)
      })
    }

    print("Start merging datasets!")

    pb <- txtProgressBar(min = 0, max = length(datasets)-1, style = 3)

    exprsALL = matrix(data = , nrow = 54675, ncol = 0)

    suppressMessages(
      for (i in 1:length(datasets)) {
        fname_exprsSet = paste0(datadir, datasets[i], ".Rdata")
        load(fname_exprsSet) #AEsetnorm
        exprsData = exprs(AEsetnorm)
        sampleNames = colnames(exprsData)
        sampleNames = unlist(strsplit(sampleNames, ".CEL.gz|.cel.gz|.CEL|.cel|_.*"))
        sampleData_dataset = subset(metadataFinal, subset=acc_dataset==datasets[i])
        sampleIDs = getSampleIDs()
        if (length(sampleIDs)!=length(sampleNames)) {
          sampleData_dataset$celFiles = getCelFileNamesFromAE(datasets[i])
          sampleNames = colnames(exprsData)
          sampleIDs = sampleData_dataset$sampleID[match(sampleNames,sampleData_dataset$celFiles)]
        }
        colnames(exprsData) = sampleIDs
        if (all(rownames(exprsALL)==rownames(exprsData))) {
          exprsALL = cbind(exprsALL, exprsData)
        } else {
          print("Error: unequal probe identifier")
          break
        }
        setTxtProgressBar(pb, i-1)
      }
    )
    save(exprsALL, file=paste0(datadir, "exprsALL.RData"))
    print("Finished merging datasets!")
  } else {
    print("Merging was already done!")
    load(file=paste0(datadir, "exprsALL.RData"))
  }

  dim(exprsALL)

  metadataFinal = subset(metadataFinal, sampleID%in%colnames(exprsALL))
  metadata = subset(metadata, sampleID%in%colnames(exprsALL))


  ###### Batch correct data
  if (!file.exists(paste0(datadir, "data_batchCorrected.RData"))) {
    print("Start batch correcting data!")
    metadata_subset = subset(metadataFinal, acc_dataset%in%datasets)
    modcombat = model.matrix(~1, data=metadata_subset)
    #combat_edata = ComBat(dat=exprsALL, batch=metadata_subset$acc_dataset, mod=modcombat,
    #                      par.prior=TRUE, prior.plots=FALSE)
    combat_edata = ComBat(dat=exprsALL, batch=as.factor(as.character(metadata_subset$acc_dataset)),
                          mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
    save(combat_edata, file=paste0(datadir, "data_batchCorrected.RData"))
    print("Finished batch correcting data!")
  } else {
    print("Batch correction was already done!")
    load(file=paste0(datadir, "data_batchCorrected.RData"))
  }


  ###### Filter for needed samples only
  if (!file.exists(paste0(datadir, "data_batchCorrected_filtered.RData"))) {
    print("Start filtering datasets for needed samples!")
    sampleIDs = metadata$sampleID
    edata_filtered = subset(combat_edata, select = colnames(combat_edata)%in%sampleIDs)
    save(edata_filtered, file=paste0(datadir, "data_batchCorrected_filtered.RData"))
    print("Finished filtering datasets for needed samples!")
  } else {
    print("Filtering datasets for needed samples was already done!")
    load(file=paste0(datadir, "data_batchCorrected_filtered.RData"))
  }


  ###### Probe filtering and annotation
  if (!file.exists(paste0(datadir, "data_annotated.RData"))) {
    print("Start probe filtering and annotation datasets!")
    eset_combat = ExpressionSet(edata_filtered, annotation = "hgu133plus2")
    eset_combat_filtered = featureFilter(eset_combat, require.entrez=TRUE, remove.dupEntrez=TRUE,
                                         feature.exclude="^AFFX")

    mapGeneNames_probe2HGNC = function(probeIDs) {
      genes_HGNC = unlist(mapIds(hgu133plus2.db,
                                 keys=probeIDs,
                                 column="SYMBOL",
                                 keytype="PROBEID",
                                 multiVals="first"))
      return(genes_HGNC)
    }

    combat_edata_filtered = exprs(eset_combat_filtered)
    
    symbols = mapGeneNames_probe2HGNC(rownames(combat_edata_filtered))
    combat_edata_annotated = combat_edata_filtered[!is.na(symbols),]
    symbols_filtered = symbols[!is.na(symbols)]
    rownames(combat_edata_annotated) = symbols_filtered

    save(combat_edata_annotated, file=paste0(datadir, "data_annotated.RData"))
    print("Finished probe filtering and annotation datasets!")
    print("Everything is prepared for WGCNA.")
  } else {
    print("Probe filtering and annotation was already done!")
    print("Everything is prepared for WGCNA. Nothing to do anymore.")
  }



}
