setwd("/home/traschka/Projects/PhD/CoExprVsPathways/data_for_coexp_network_construction/")

library(WGCNA)

options(stringsAsFactors = FALSE)

#groupIDs = c('1612', '8692', '9256', '9538', '1324', '3908', '0070004', '707', '3068', '9952', '2394',
#             '3459', '2841', '1040', '8893', '5520', '0040085', '1909', '1793', '4450', '5603', '3083',
#             '0050908', '263', '9970', '10652', '3069', '0050866', '9778', '3910',  '169', '0060074', 
#             '2377', '8567', '10283', '5419', '1115', '657', '0050749', '3070', '5844', '4074', '3247',
#             '1192', '3620', '3382', '1967', '1596', '0050902', '9261', '9074', '289', '7502', '3312', 
#             '8778', '1470', '9253', '3310', '11335', '10223', '0080199', '9744', '1883', 'normal')

groupIDs = c('1612','9119','9256','9538','1324','3069','3908','0070004','9952','2394',
             '0050745','2841','1040','8893','5520','0040085','1793','1909','5603','3083',
              '0050908','263','4450','9970','10652','0050866','9778','3910','169','0060074',
              '2377','8567','0050746','10283','5419','1115','657','0050749','0050750','5844',
              '3070','4074','3247','769','3382','1967','1596','9261','0050902','9074','289',
              '0050909','3234','7502','3312','1470','8778','9253','11335','3310','10223','9744',
              '0080199','1883','normal')

skippedIDs = c('263')

for(ID in groupIDs) {
  
  print(paste0("Next group: ", ID))
  
  if (ID%in%skippedIDs) {
    next
  }
  
  # for ID in groupIDs:
  #phenotype = ID
  #phenotype = ""
  
  datadir = paste0("./", ID, "/")
  
  if(file.exists(paste0(datadir, "coexp_network_edges.tsv"))) {
    next
  }
  
  data = get(load(paste0(datadir, "data_annotated.RData"))) # or however the files will be named
  powers = c(seq(from = 2, to=12, by=2))                                              
  sft = pickSoftThreshold(t(data), powerVector = powers, verbose = 5)                 #
  power = sft$powerEstimate                                                           #
  if (is.na(power)) {                                                                 #
    if ((max(sft$fitIndices[,2]))>0.7) {                                              #
      power = sft$fitIndices[which(sft$fitIndices[,2]==max(sft$fitIndices[,2])),1]    #
    } else {
      if (ID=="9538") {
        power = 30
      } else if (ID=="9952") {
        power = 30
      } else if (ID=="1909") {
        power = 30
      } else if (ID=="0050745") {
        power = 30
      } else {
        powers = c(seq(from = 12, to=30, by=2))                                              
        sft = pickSoftThreshold(t(data), powerVector = powers, verbose = 5)                 #
        power = sft$powerEstimate
      }
    }
  }
  print(power)
  if (!file.exists(paste0(datadir, "coexp_network.rda"))) {
    net = blockwiseModules(t(data), power = power, maxBlockSize = nrow(data),
                           minModuleSize = 20, verbose = 4, #numericLabels = TRUE, 
                           saveTOMs = TRUE, loadTOM = FALSE,
                           saveTOMFileBase = paste0(datadir, "coexp_network"))
    save(net, file = paste0(datadir, "coexp_network.rda"))
  } else {
    load(file=paste0(datadir, "coexp_network.rda"))
  }
  
  # TOM <- TOMsimilarityFromExpr(datExpr, power = power)
  load(paste0(datadir, "coexp_network-block.1.RData"))
  #add names to TOM
  
  vis <- exportNetworkToVisANT(TOM,
                               #file = paste0(datadir, "VisANTInput-", ID),
                               weighted = TRUE,
                               threshold = 0,
                               probeToGene = data.frame(c(seq(from = 1, to=length(rownames(data)))), rownames(data)) )
  # cyt <- exportNetworkToCytoscape(TOM,
  #                                 edgeFile = paste("CytoscapeInput-edges-", ID, ".txt", sep=""),
  #                                 nodeFile = paste("CytoscapeInput-nodes-", ID, ".txt", sep=""),
  #                                 weighted = TRUE,
  #                                 threshold = 0.02,
  #                                 nodeNames = rownames(data),
  #                                 nodeAttr = net$colors);
  # cyt$edgeData

  threshold = quantile(vis$weight,0.99)
  print(threshold)
  vis2 = subset(vis, weight>=threshold)
  
  write.table(net$colors, file=paste0(datadir, "coexp_network_modules.tsv"), sep="\t", quote=F, row.names=T, col.names=F)
  write.table(vis2, file=paste0(datadir, "coexp_network_edges.tsv"), sep="\t", quote=F, row.names=F, col.names=T)
  
  rm(list=setdiff(ls(), c("groupIDs", "ID", "skippedIDs")))
}