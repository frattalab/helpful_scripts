
DESeqDataSetFromIRFinder <-  function(filePaths,designMatrix,designFormula){
    res=c()
    libsz=c()
    spl=c()
    irtest=read.table(filePaths[1])
    if (irtest[1,1]=="Chr"){irtest=irtest[-1,]}
    irnames=unname(apply(as.matrix(irtest),1,FUN=function(x){return(paste0(x[4],"/",x[1],":",x[2],"-",x[3],":",x[6]))}))
    n=1
    for (i in filePaths){
        print(paste0("processing file ",n," at ",i))
        irtab=read.table(i)
        if (irtab[1,1]=="Chr"){irtab=irtab[-1,]}
        #rn=unname(apply(irtab,1,FUN=function(x){return(paste0(x[4],"/",x[1],":",x[2],"-",x[3],":",x[6]))}))
        #row.names(irtab)=rn
        #tmp1=round(as.numeric(as.vector(irtab[irnames,9])))
        #tmp2=as.numeric(as.vector(irtab[irnames,19]))
        tmp1=as.numeric(as.vector(irtab[,9]))
        tmp2=as.numeric(as.vector(irtab[,19]))
        tmp3=tmp1+tmp2
        tmp4=as.numeric(as.vector(irtab[,17]))
        tmp5=as.numeric(as.vector(irtab[,18]))
        tmp6=pmax(tmp4,tmp5, na.rm=T)
        res=cbind(res,tmp1)
        libsz=cbind(libsz,tmp2)
        spl=cbind(spl,tmp6)
        n=n+1
    }
    res.rd=round(res)
    libsz.rd=round(libsz)
    spl.rd=round(spl)
    colnames(res.rd)=paste("intronDepth",as.vector(designMatrix[,1]),sep=".")
    rownames(res.rd)=irnames
    colnames(libsz.rd)=paste("totalSplice",as.vector(designMatrix[,1]),sep=".")
    rownames(libsz.rd)=irnames
    colnames(spl.rd)=paste("maxSplice",as.vector(designMatrix[,1]),sep=".")
    rownames(spl.rd)=irnames
    
    ir=c(rep("IR",dim(designMatrix)[1]),rep("Splice",dim(designMatrix)[1]))
    group=rbind(designMatrix,designMatrix)
    group$IRFinder=ir
    group$IRFinder=factor(group$IRFinder,levels=c("Splice","IR"))
    
    #counts.IRFinder=cbind(res.rd,libsz.rd)
    counts.IRFinder=cbind(res.rd,spl.rd)
    
    dd = DESeqDataSetFromMatrix(countData = counts.IRFinder, colData = group, design = designFormula)
    sizeFactors(dd)=rep(1,dim(group)[1])
    rownames(dd)=irnames
    final=list(dd,res,libsz,spl)
    names(final)=c("DESeq2Object","IntronDepth","SpliceDepth","MaxSplice")
    return(final)
}

make_deframe_irfinder <-  function(ir_top_level, 
                                 dir_type = "IRFinder-IR-nondir.txt",
                                 base_grep = "Control",
                                 contrast_grep = "TDP43KD", 
                                 baseName = "control",
                                 contrastName = 'TDPKD'){
    paths <- ir_top_level %>% 
        list.files(,full.names = TRUE) %>% 
        file.path(., dir_type) %>% 
        as.vector()
    
    SampleNames <-  paths %>% 
        gsub(paste0("/",dir_type),"",.) %>% basename()
    
    experiment <-  as.data.frame(SampleNames) 
    
    experiment <- experiment %>% 
        mutate(Condition = as.factor(ifelse(grepl(contrast_grep,SampleNames),contrastName,baseName)))
    experiment$Condition = factor(experiment$Condition, levels = c(baseName,contrastName))
    rownames(experiment)=NULL
    
    morphed = list(paths,experiment)
    
    names(morphed) = c("paths","metadata")
    return(morphed)
}


run_deseq_ir = function(paths, experiment){
    metaList=DESeqDataSetFromIRFinder(filePaths=paths, designMatrix=experiment, designFormula=~1)
    dds = metaList$DESeq2Object                       # Extract DESeq2 Object with normalization factors ready

    design(dds) = ~Condition + Condition:IRFinder     # Build a formula of GLM. Read below for more details. 
    dds = DESeq(dds)                                  # Estimate parameters and fit to model
    
    return(dds)
    
}


calculate_ir = function(dds, result_name){
    
    intron_result = results(dds, name = result_name)
    IR_vs_splice=2^intron_result$log2FoldChange
    IRratio = IR_vs_splice/(1+IR_vs_splice)
    return(IRratio)
}

return_formated_results = function(dds){
    
    result_names = resultsNames(dds)
    ir_contrast = result_names[2]
    ir_cond1 = result_names[3]
    ir_cond2 = result_names[4]
    message(glue::glue("IR Condition 1 {ir_cond1}"))
    message(glue::glue("IR Condition 2 {ir_cond2}"))
    
    ir_ratio_cond1 = calculate_ir(dds, ir_cond1)
    column_name_c1 = gsub(".IRFinderIR","",gsub("Condition","",ir_cond1))
    

    ir_ratio_cond2 = calculate_ir(dds, ir_cond2)
    column_name_c2 = gsub(".IRFinderIR","",gsub("Condition","",ir_cond2))
    
    
    
    res.diff = results(dds, contrast=list(ir_cond2,ir_cond1))  
    res.diff <- res.diff %>% 
        as.data.frame() %>% 
        rownames_to_column('event') %>% 
        separate(event, into = c("symbol","gene","type","location"),sep = "\\/") %>% 
        separate(location, into = c("chr","start","end","strand")) %>% 
        mutate(paste_me_into_igv = paste0(chr,":",start,"-",end))
    
    res.diff <- res.diff %>% 
        mutate(!!column_name_c1 := ir_ratio_cond1) %>% 
        mutate(!!column_name_c2 := ir_ratio_cond2) %>% 
        as.data.table()
        
    return(res.diff)
    
}


