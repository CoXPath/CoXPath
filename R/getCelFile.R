getCelFileNamesFromAE = function(dataset) {
  source("../../scripts/helperFunctions.R")
  res = fromJSON(paste0("https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/", dataset, "/samples"))
  samples = res$experiment$sample
  files = samples$file
  if(length(unique(files[[1]]$type))==length(files[[1]]$type)) {   # only try to find CEL file, if keys are unique, otherwise skip. I have to get the data for some of the datasets manually anyway.
    if(checkExists_files(files)) {
      fileNA = which(unlist(lapply(files, FUN = is.null)))
      fileNotNA = which(!unlist(lapply(files, FUN = is.null)))
      if (length(fileNA)>0) {
        files_woNA = files[-fileNA]
      } else {
        files_woNA = files
      }
      files_parsed_woNA = rbindlist(lapply(files_woNA, FUN = function(x) {
        x %>%
          dplyr::select(c("type", "name")) %>%
          spread(key = type, value = name)
      }), fill = TRUE)
      files_parsed = data.frame(matrix(vector(), NROW(sampleData_dataset$acc_sample), NCOL(files_parsed_woNA)),
                                stringsAsFactors=F)
      colnames(files_parsed) = colnames(files_parsed_woNA)
      files_parsed[fileNotNA,] = files_parsed_woNA
      whichColCEL = names(files_parsed)[which(endsWith(as.character(files_parsed[1,]), ".CEL"))]
      if (length(whichColCEL)==0) {
        whichColCEL = names(files_parsed)[which(endsWith(as.character(files_parsed[1,]), ".cel"))]
      }
      if (length(whichColCEL)!=0) {
        cel_file = files_parsed[, whichColCEL]
      } else {
        cel_file = rep(NA, length(sampleData_gse$acc_sample))
      }
    } else {
      print("stop2")
      break
    }
  }
  return(cel_file)
}


# source("../scripts/helperFunctions.R")
# setwd("/home/traschka/Projects/PhD/phd-project/zml-project/data/")
# setwd("C://Users/traschka/Documents/Projects/PhD/zml-project/data/")
# source("../scripts/helperFunctions.R")
# 
# ###---------------------- get sample data from AE
# library(ArrayExpress)
# library(jsonlite)
# library(tidyr)
# library(dplyr)
# library(data.table)
# 
# load("2020-10-08_sampleDataClean.Rdata")
# datasetsClean = unique(sampleDataClean$acc_dataset)
# 
# options(stringsAsFactors=FALSE)
# celFilesALL = data.frame()
# #nbr_samples = 0
# #part = 9
# 
# pb <- txtProgressBar(min = 0, max = length(datasetsClean)-1, style = 3)
# 
# suppressMessages(
#   for (i in 1:length(datasetsClean)) {
#     gse = datasetsClean[i]
#     print(gse)
#     
#     sampleData_gse = subset(sampleDataClean, acc_dataset==gse)
# 
#     res = fromJSON(paste0("https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/", gse, "/samples"))
#     samples = res$experiment$sample
#     if(is.null(samples)) {
#       print(paste0("This dataset (", gse, ") does not contain any sample."))
#       next
#     }
#     if (gse=="E-GEOD-15735") {
#       samples = samples[-10,] # only no variables for the 10th line, all samples are multiple times included, therefore just excluded this line
#     }
#     
#     # raw file name
#     files = samples$file
#     noFiles = c("E-GEOD-52847", "E-MTAB-1328", "E-GEOD-11956", "E-GEOD-25639", "E-AFMX-11", "E-GEOD-21687",
#                 "E-GEOD-13319")
#     # E-GEOD-25639, E-AFMX-11, E-GEOD-21687, E-GEOD-13319: different numbers of samples due to mouse data -> will be excluded completly anyway..
#     if (gse%in%noFiles) {
#       cel_file = rep(NA, length(sampleData_gse$acc_sample))
#     } else if(length(unique(files[[1]]$type))==length(files[[1]]$type)) {   # only try to find CEL file, if keys are unique, otherwise skip. I have to get the data for some of the datasets manually anyway.
#       if(checkExists_files(files)) {
#         fileNA = which(unlist(lapply(files, FUN = is.null)))
#         fileNotNA = which(!unlist(lapply(files, FUN = is.null)))
#         if (length(fileNA)>0) {
#           files_woNA = files[-fileNA]
#         } else {
#           files_woNA = files
#         }
#         files_parsed_woNA = rbindlist(lapply(files_woNA, FUN = function(x) {
#           x %>%
#             select(c("type", "name")) %>%
#             spread(key = type, value = name)
#         }), fill = TRUE)
#         files_parsed = data.frame(matrix(vector(), length(sampleData_gse$acc_sample), NCOL(files_parsed_woNA)),
#                                   stringsAsFactors=F)
#         colnames(files_parsed) = colnames(files_parsed_woNA)
#         files_parsed[fileNotNA,] = files_parsed_woNA
#         whichColCEL = names(files_parsed)[which(endsWith(as.character(files_parsed[1,]), ".CEL"))]
#         if (length(whichColCEL)==0) {
#           whichColCEL = names(files_parsed)[which(endsWith(as.character(files_parsed[1,]), ".cel"))]
#         }
#         if (length(whichColCEL)!=0) {
#           cel_file = files_parsed[, whichColCEL]
#         } else {
#           cel_file = rep(NA, length(sampleData_gse$acc_sample))
#         }
#       } else {
#         print("stop2")
#         break
#         files_parsed = rbindlist(lapply(files, FUN = function(x) {
#           x %>%
#             select(c("type", "name")) %>%
#             spread(key = type, value = name)
#         }), fill = TRUE)
#         whichColCEL = names(files_parsed)[which(endsWith(as.character(files_parsed[1,]), ".CEL"))]
#         if (length(whichColCEL)!=0) {
#           whichColCEL = names(files_parsed)[which(endsWith(as.character(files_parsed[1,]), ".cel"))]
#         }
#         if (length(whichColCEL)!=0) {
#           cel_file = files_parsed[, get(whichColCEL)]
#         } else {
#           print("Problem with cel-file at :")
#           print(gse)
#           cel_file = rep(NA, length(sampleData_gse$acc_sample))
#         }
#       }
#     } else {
#       cel_file = rep(NA, length(sampleData_gse$acc_sample))
#     }
#     
#     cel_file_df = data.frame(sampleData_gse$acc_sample, cel_file)
#     
#     celFilesALL = dplyr::bind_rows(celFilesALL, cel_file_df)
#     
#     rm(list=setdiff(ls(), c("datasetsClean", "celFilesALL", "sampleDataClean", "checkExists_files", "i", "pb")))
#     
#     setTxtProgressBar(pb, i-1)
#   }
# )
# 
# write.table(celFilesALL, file=paste0(Sys.Date(), "_celFiles.tsv"), sep="\t", quote=F, row.names=F, col.names=T)
# save(celFilesALL, file=paste0(Sys.Date(), "_celFiles.Rdata"))
# 
