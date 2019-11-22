#this function takes the total rna tables produced by featureCounts, and gives a reasonable output data frame and
#metadata frame
make_deseq_dfs = function(total_table, grep_pattern = "_new_readcount", leave_out = ""){
  #grep pattern is being used to select small parts of this overall
  if(leave_out == ""){
    conv_df = as.data.frame(total_table[,round(.SD), .SDcols = grep(grep_pattern,names(total_table))])
  }else{
    conv_df = as.data.frame(total_table[,round(.SD), .SDcols = grep(grep_pattern,names(total_table))])
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
  coldata[grep("a.a",V1), cond := "base"]
  coldata[grep("b.a",V1), cond := "during"]
  coldata[grep("c.a",V1), cond := "post"]
  coldata[grep("wt",V1), cell := "wt"]
  coldata[grep("f",V1), cell := "f210i"]

  coldata = as.data.frame(coldata[,2:ncol(coldata)])
  rownames(coldata) = names(conv_df)
  morphed = list(conv_df,coldata)
  names(morphed) = c("conv_df","coldata")
  print(names(conv_df))

  return(morphed)
}
