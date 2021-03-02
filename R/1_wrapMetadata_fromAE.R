# path
setwd("/home/rebecafig/Documents/fraunhofer/co-expression-vs-pathways/data/")
source("../scripts/helperFunctions.R")

###---------------------- load GSE list
acc_list = as.vector(as.matrix(read.table("accession_gpl570ArrayExpress.csv")))
excAcc = read.csv("errorDatasets.csv", header = F)$V2

###---------------------- get sample data from AE
library(ArrayExpress)
library(jsonlite)
library(tidyr)
library(dplyr)
library(data.table)

options(stringsAsFactors=FALSE)
sampleData = data.frame()
nbr_samples = 0
part = 0

#rdm_gse = "E-GEOD-28792"    #"E-GEOD-12345"     #"E-MEXP-1621"    #"E-GEOD-68694"    #acc_list[sample(0:length(acc_list), 1)]
pb <- txtProgressBar(min = 0, max = length(acc_list)-1, style = 3)

# Iter over all geo accession numbers
suppressMessages(
#  for (i in 1:length(acc_list)) { # part 0
#  for (i in 171:length(acc_list)) { # part 1
#  for (i in 547:length(acc_list)) { # part 2
#  for (i in 1076:length(acc_list)) { # part 3
#  for (i in 1550:length(acc_list)) { # part 4
#  for (i in 2089:length(acc_list)) { # part 5
#  for (i in 2614:length(acc_list)) { # part 6
#  for (i in 3121:length(acc_list)) { # part 7
#  for (i in 3662:length(acc_list)) { # part 8
  for (i in 4219:length(acc_list)) { # part 9
    metadata = c()
    gse = acc_list[i]
    print(gse)
    if (gse%in%excAcc) {
      next
    }
    
    #rdm_gse = acc_list[sample(0:length(acc_list), 1)]
    #print(rdm_gse)
    
    # Get data from the API
    Sys.sleep(.5)
    print("api")
    res = fromJSON(paste0("https://www.ebi.ac.uk/arrayexpress/json/v3/experiments/", gse, "/samples"))
    print("api done")
    samples = res$experiment$sample
    if(is.null(samples)) {
      print(paste0("This dataset (", gse, ") does not contain any sample."))
      next
    }
    if (gse=="E-GEOD-15735") {
      samples = samples[-10,] # only no variables for the 10th line, all samples are multiple times included, therefore just excluded this line
    }
    
    
    # accession
    if (all(startsWith(as.character(samples$assay[,1]), "GSM"))) {
      sample_acc = samples$assay$name
    } else {
      if(any(samples$characteristic[[1]]$category=="GEO sample accession")) {
        acc = lapply(X = samples$characteristic, FUN = function(x) {
          gsm = x %>% filter(category=="GEO sample accession") %>% select(value)
        })
        sample_acc = unlist(acc)
      } else {
        #print(paste0("This dataset (", rdm_gse, ") has no GSM accession in the assay field."))
        gsmWhere = grep("GSM", samples)
        if (length(gsmWhere)==0) {
          #print("No GSM accession could be found.")
          # just use the names of the samples.
          sample_acc = samples$assay$name
        } else {
          sapply(gsmWhere, function(x) {
            length(grep("GSM", unlist(samples[gsmWhere])))
          })
          
          
          gsms_vector = unlist(samples[gsmWhere[1]])
          names(gsms_vector) = c()
          gsmForRegExCheck = gsms_vector[grep("GSM", gsms_vector)[1]]
          if (grepl("GSE[0-9]*GSM[0-9]*", gsmForRegExCheck)) {
            #print("GSM accession could be found in other fields.")
            # in case it looks like: "GSE12345GSM309986 something"but also "GSE8059GSM198943"
            regex = "(?<=.)(?=(GSE[0-9]*)|(GSM[0-9]*)| )"
            acc = lapply(X = gsms_vector, FUN = function(x) {
              splt = strsplit(x, split=regex, perl = TRUE)[[1]]
              gsm = splt[grep("GSM", splt)]
            })
            acc[sapply(acc, function(x) {length(x)==0})] <- NA
            sample_acc = unlist(acc)
          } else if (grepl("E-GEOD-[0-9]*_GSM[0-9]*.CEL", gsmForRegExCheck)) {
            # in case it looks like: "E-GEOD-19665_GSM490987.CEL"
            regex = "(?<=.)(?=(E-GEOD-[0-9]*_)|(GSM[0-9]*)|.CEL)"
            acc = lapply(X = gsms_vector, FUN = function(x) {
              splt = strsplit(x, split=regex, perl = TRUE)[[1]]
              gsm = splt[grep("GSM", splt)]
            })
            acc[sapply(acc, function(x) {length(x)==0})] <- NA
            sample_acc = unlist(acc)
          } else if (grepl("GSM[0-9]*", gsmForRegExCheck[1])){
            # in case it looks like: "GSM318121_Wh" or "GSM107522 extract"
            acc = lapply(X = gsms_vector, FUN = function(x) {
              splt = strsplit(x, split="_| ", perl = TRUE)[[1]]
              gsm = splt[grep("GSM", splt)]
            })
            sample_acc = unlist(acc)
          } else {
            print("Problem with gsm vector.")
            print(gsmForRegExCheck)
            break
          }
        }
      }
    }
    sample_acc = as.character(sample_acc)
    
    # characteristics
    characteristics = samples$characteristic
    if (checkSharedKeys_char(characteristics)) {
      charCollapseDF = matrix(NA, ncol = length(unique(tolower(characteristics[[1]]$category))), nrow = length(samples$assay$name))
      colnames(charCollapseDF) = unique(tolower(characteristics[[1]]$category))
      for (j in unique(tolower(characteristics[[1]]$category))) {
        charCollapse = lapply(characteristics, FUN = function(x) {
          charFilter = x %>% filter(tolower(category)==j)
          charValues = charFilter$value[!is.na(charFilter$value)]
          if (length(charValues)>1) {
            charValues = paste(trimws(charFilter$value), collapse=";")
            #print(charValues)
          }
          return(charValues)
        })
        charCollapse[sapply(charCollapse, FUN = function(x) {
          length(x) == 0
        })] = NA
        charCollapseDF[,j] = unlist(charCollapse)
      }
      #colnames(charCollapseDF) = unique(tolower(characteristics[[1]]$category))
      characteristics_parsed = as.data.frame(charCollapseDF)
    } else {
      characteristics_parsed = rbindlist(lapply(characteristics, FUN = function(x) {
        x %>%
          select(c(category, value)) %>%
          spread(key = category, value = value)
      }), fill = TRUE)
    }

    
    # # sample source
    # sample_source = samples$source$comment # only samples$source??
    # if (is.null(dim(sample_source))) {
    #   source_parsed = rbindlist(lapply(sample_source, FUN = function(x) {
    #     x %>%
    #       spread(key = name, value = value)
    #   }))
    # } else {
    #   sample_source = tibble::rowid_to_column(sample_source)
    #   source_parsed = sample_source %>% 
    #     spread(key = name, value = value) %>%
    #     select(Sample_source_name)
    # }
    # 
    
    # variables
    variables = samples$variable
    if (!is.null(variables)) {
      if (length(unique(variables[[1]]$name))==length(variables[[1]]$name)) { #no duplicated names
        variables_parsed = rbindlist(lapply(variables, FUN = function(x) {
          x %>%
            select(c(name, value)) %>%
            spread(key = name, value = value)
        }))
        variables_parsed = sapply(variables_parsed, as.character)
        names(variables_parsed) = paste0("var_", names(variables_parsed))
      } else { # duplicated names
        varCollapseDF = matrix(NA, ncol = length(unique(tolower(variables[[1]]$name))), nrow = length(samples$assay$name))
        colnames(varCollapseDF) = unique(tolower(variables[[1]]$name))
        for (j in unique(tolower(variables[[1]]$name))) {
          varCollapse = lapply(variables, FUN = function(x) {
            varFilter = x %>% filter(tolower(name)==j)
            varFilter$value[varFilter$value=="not specified"] = NA
            varValues = varFilter$value[!is.na(varFilter$value)]
            if (length(varValues)>1) {
              varValues = paste(varFilter$value, collapse=";")
            }
            return(varValues)
          })
          varCollapse[sapply(varCollapse, FUN = function(x) {
            length(x) == 0
          })] = NA
          varCollapseDF[,j] = unlist(varCollapse)
        }
        colnames(varCollapseDF) = paste0("var_", colnames(varCollapseDF))
        variables_parsed = as.data.frame(varCollapseDF)
        variables_parsed = sapply(variables_parsed, as.character)
      }
    } else {
      variables_parsed = rep(NA, length(sample_acc))
    }
    
    ### put all together
    metadata = data.frame(sample_acc, gse, characteristics_parsed, variables_parsed) #source_parsed,
    names(metadata)[1:2] = c("acc_sample", "acc_dataset")
    metadata %>% mutate_all(as.character)
    
    #start = i + 1
    
    nbr_samples = nbr_samples + NROW(metadata)
    
    ## add to overall list
    sampleData = dplyr::bind_rows(sampleData, metadata)
    rm(list=setdiff(ls(), c("acc_list", "excAcc", "sampleData", "pb", "i", "part", "checkExists_files", "checkSharedKeys_char", "nbr_samples")))
    
    setTxtProgressBar(pb, i-1)
    
    print("part")
    print(part)
    if(nbr_samples>20000) {
      write.table(sampleData, file=paste0("sampleData_Part", part, ".tsv"), sep="\t", quote=F, row.names=F, col.names=T)
      save(sampleData, file=paste0("sampleData_Part", part, ".Rdata"))
      part = part + 1
      print(i)
      print(nbr_samples)
      nbr_samples = 0
      rm(sampleData)
      sampleData = data.frame()
      break
    }
  }
)
close(pb)

#### 2020-10-08

# #####---------------------###
# 
# colnamesALL = c()
# #sampleDataALL = list()
# for (i in 0:9) {
#   load(paste0("2020-10-08_sampleData_Part", i, ".Rdata"))
#   #sampleDataALL[[i+1]] = sampleData
#   colnamesALL = c(colnamesALL, colnames(sampleData))
# }
# 
# write(unique(colnamesALL), file="colnamesALL.csv", ncolumns = 1)
# 
# #d = adist(colnamesALL)   # find distance matrix 
# #hc = hclust(as.dist(d))                # apply hirarchical clustering 
# #plot(hc)     

