#this function takes the total rna tables produced by featureCounts, and gives a reasonable output data frame and
#metadata frame
make_deseq_dfs = function(total_table, grep_pattern = "", leave_out = "", base_grep = "", contrast_grep = ""){

  if(grep_pattern == ""){
    grep_pattern = paste0(colnames(total_table[,2:length(total_table)]),collapse = "|")
  }
  #grep pattern is being used to select small parts of this overall
  total_table = as.data.table(total_table, keep.rownames = TRUE)

  if(leave_out == ""){
    conv_df = as.data.frame(total_table[,round(.SD), 
                                        .SDcols = grep(grep_pattern,names(total_table))])
  }else{
    conv_df = as.data.frame(total_table[,round(.SD), 
                                        .SDcols = grep(grep_pattern,names(total_table))])
    drop_col = paste0(leave_out_sample,".aligned.sorted.out.bam")
    print(paste("Dropping column:",drop_col ))
    conv_df = conv_df[ , !(names(conv_df) %in% drop_col)]
  }

  if("geneid" %in% names(total_table)){
    rownames(conv_df) = total_table$geneid
  }else if("Geneid" %in% names(total_table)){
    rownames(conv_df) = total_table$Geneid
  }else if("rn" %in% names(total_table)){
    rownames(conv_df) = total_table$rn
  }
  else{
    rownames(conv_df) = total_table$gene
  }
  coldata = as.data.table(names(conv_df))
  if(base_grep == "" & contrast_grep == ""){
    coldata[grep("a.a",V1), cond := "base"]
    coldata[grep("b.a",V1), cond := "during"]
    coldata[grep("c.a",V1), cond := "post"]
    coldata[grep("wt",V1), cell := "wt"]
    coldata[grep("f",V1), cell := "f210i"]
  }else if(base_grep != ""){
    coldata[grep(base_grep,V1), cond := "base"]
  }else if(contrast_grep != ""){
    coldata[grep(contrast_grep,V1), cond := "contrast"]
  }
  coldata[is.na(cond), cond := "contrast"]
  
  coldata = as.data.frame(coldata[,2:ncol(coldata)])
  rownames(coldata) = names(conv_df)
  morphed = list(conv_df,coldata)
  names(morphed) = c("conv_df","coldata")
  print(names(conv_df))

  return(morphed)
}


# this function takes a deseq metadata df, the name you want to call the baseline and contrast conditions and returns the metadata table with
# a new column called 'comparison' which is a factor with baseline as the first level
rename_relevel_for_deseq = function(coldata, baseName = "", contrastName = ""){
  coldata$comparison  = factor(ifelse(coldata$cond == "base", baseName, contrastName), levels = c(baseName, contrastName))
  return(coldata)
}
