
create_feature_count_table = function(feature_count_folder, suffix = ".Aligned.sorted.out.bam",prefix = ""){
  library(data.table)
  # feature counts gives you feature coutns and summariues, given a folder name get a list of files that don't contain
  # the word summary
  feature_count_files = grep(list.files(path=feature_count_folder,full.names = TRUE), pattern='summary', inv=T, value=T)
  # from the list, read them all
  l <- lapply(feature_count_files, fread)
  # make it into a data table
  wide_feature_counts = setDT(unlist(l, recursive = FALSE), check.names = TRUE)[]
  # there's only a few columns we really want from this table it's going to be
  # Geneid and the things that labeled by their file name so we do two things
  # first when we do this setDT unlist thing names become X.path.to. so we'lre gooing to take the feature counts files and replace the / with .
  selector_list = c("Geneid",names(wide_feature_counts)[grepl("^X",names(wide_feature_counts))])
# take only those columns
  wide_feature_counts = wide_feature_counts[,.SD, .SDcols = selector_list]
  # ----remove unnecessary information from the column names
  if(suffix != ""){
    colnames(wide_feature_counts) = gsub(suffix,"",colnames(wide_feature_counts))
  }
  if(prefix != ""){
    colnames(wide_feature_counts) = gsub(prefix,"",colnames(wide_feature_counts))
  }
  return(wide_feature_counts)
}
