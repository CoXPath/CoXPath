checkSharedKeys_char <- function(characteristics) {
  ## This function checks if there are shared keys in the given list. 
  ### TRUE: shared keys in at least one sample
  ### FALSE: no shared keys for any sample
  return(!all(sapply(characteristics, FUN = function(x) {length(unique(tolower(x$category)))==length(tolower(x$category))})))
}

checkExists_files <- function(files) {
  ## This function checks if there are NULL fields in the given list. 
  ### TRUE: NULL fields in at least one sample
  ### FALSE: no NULL fields for any sample
  return(!all(unlist(lapply(files, FUN = is.null))))
}