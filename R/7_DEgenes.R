setwd("/home/traschka/Projects/PhD/CoExprVsPathways/data_for_coexp_network_construction/")

options(stringsAsFactors = FALSE)

library(limma)

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
  
  datadir = paste0("./", ID, "/")
  
  if(file.exists(paste0(datadir, "DEgenes.tsv"))) {
    next
  }
  
  #load the metadata file
  metadataFinal = read.table(file = "metadata_FINAL.tsv", sep = '\t', header = TRUE, quote = "")
  
  #load the merged dataset for the disease
  load(file=paste0(datadir, "data_batchCorrected.RData")) #combat_edata
  sampleIDs = colnames(combat_edata)
  
  #annotate the merged dataset
  ###### Probe filtering and annotation
  library("genefilter")
  library("hgu133plus2.db")
  eset_combat = ExpressionSet(combat_edata, annotation = "hgu133plus2")
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
    
  
  #filter metadata for needed datasets
  metadata = subset(metadataFinal, sampleID%in%sampleIDs)
  metadata_sub = subset(metadata, select=c("sampleID", "diseaseTerm", "diseaseID"))
  ## order metadata in same way as data
  metadata_sub = metadata_sub[order(match(metadata_sub$sampleID,colnames(combat_edata_annotated))),]
  ## get only "number" of disease ID to make easier comparison later
  metadata_sub$diseaseID = substr(metadata_sub$diseaseID, 37, 1000)
  ## encode "normal" with ID="0000"
  metadata_sub$diseaseID[metadata_sub$diseaseTerm=="normal"] = "0000"
  metadata_sub$diseaseID[metadata_sub$diseaseTerm==""] = "exclude"
  all(metadata_sub$sampleID==colnames(combat_edata_annotated))
  
  #define groups for differential analysis
  ## one group for disease and one for normal
  if ("0000"%in%unique(metadata_sub$diseaseID)) {
    diseases = setdiff(unique(metadata_sub$diseaseID), "0000")
    group = factor(metadata_sub$diseaseID, levels = c("0000", diseases))
    design1 = model.matrix(~group)
  } else {
    print("No normal data")
    next
  }
  
  #do differential analysis
  fit1 = lmFit(combat_edata_annotated, design=design1)
  efit1 <- eBayes(fit1)
  top.table1 <- topTable(efit1, n = Inf, coef=paste0("group",ID))
  top.table1 = subset(top.table1, adj.P.Val<0.05)
  top.table1$HGNC = rownames(top.table1)
  print(dim(top.table1))
  if (NROW(top.table1)==0) {
    print("No DE genes")
  } else {
    write.table(top.table1, file=paste0(datadir, "DEgenes.tsv"), sep="\t", quote=F, row.names=F, col.names=T)
  }

  rm(list=setdiff(ls(), c("groupIDs", "ID", "skippedIDs", "metadataFinal")))

}