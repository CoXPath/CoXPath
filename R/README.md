# R Scripts to systematically generate co-expression networks from ArrayExpress

## Scripts
- *1_wrapMetadata_fromAE.R*: get all the metadata from ArrayExpress API for a given list of datasets and store it in a dataframe
- *2_wrapMetadata_fromAE_cleanUp.R*: remove non-patient data, get clean version of metadata by summarizing some columns to a single one (e.g. disease or celltype)
- *3-1_checkCuration.R*: after manual curation check consistency of terms
- *3-2_concatCleanUpCuration.R*: process curation to get mapping between original and DOID terms
- *3-3_diseaseCleanUp.R*: map original disease terms in metadata to DOID terms and remove non-patient data
- *4_downloadData.R*: download remaining (mapable) datasets
- *5_prepDataForWGCNA.R*: preprocess downloaded datasets individually, concatenate disease data, batch correct data, and annotate probes
- *6_getNetworkFromExpData.R*: run WGCNA for each disease
- *7_DEgenes.R*: get DE genes for each disease

### helper functions
- *getCelFile.R*: get the .cel-file name for each sample, if available
- *helperFunctions.R*: checks shared keys in metadata JSON and checks if files are available
