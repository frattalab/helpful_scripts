library("pheatmap")

select_n_top_sig <- function(standard_deseq_output, n_genes = 35,padj_cut = 0.1,direction = "up"){
    if(direction == 'up'){
        #I'm going to select the geneids for the top n upregulated genes
        the_genes <- standard_deseq_output$results_table %>% 
            filter(padj < padj_cut) %>% 
            slice_max(log2FoldChange, n = n_genes) %>% 
            pull(Geneid)
    }
    if(direction == "down"){
    
        #I'm going to select the geneids for the bottom n upregulated genes
        the_genes <- standard_deseq_output$results_table %>% 
            filter(padj < padj_cut) %>% 
            slice_min(log2FoldChange, n = n_genes) %>% 
            pull(Geneid)
    }
    
    
    return(the_genes)

    
}


generate_pheatmaps <-  function(standard_deseq_output, gene_ids){
    normalized_deseq_obj <- normTransform(standard_deseq_output$deseq_obj) #normalize the counts
    
    #we're taking the colData from the DESeq object (which is metadata)
    
    
    ph_annotation <- as.data.frame(colData(standard_deseq_output$deseq_obj)) %>% 
        dplyr::select(comparison) #if you just check what this: colData(one_hour_dds$deseq_obj) you see that we have 2 additional columns which we don't want on a heatmap
    

    #I'm going to select the geneids for the top n upregulated genes
    norm_counts <- assay(normalized_deseq_obj) %>% #Normalized counts are stored in the "assay" slot on this object
        as.data.frame() %>% #it's not a data frame so I have to make it one
        rownames_to_column('Geneid') %>% #our old friend
        filter(Geneid %in% gene_ids) %>% #select these genes
        left_join(standard_deseq_output$results_table %>% #left join from the results table to get the gene_name
                      dplyr::select(Geneid,gene_name)) %>% #notice that i'm selecting inside the left_join
        mutate(gene_name = make.unique(gene_name)) %>% #this is to save me from having to do this incase the gene_name are like repeated
        column_to_rownames('gene_name') %>% #the heatmap function pheatmap uses rownames for plotting
        dplyr::select(-Geneid) #get rid of Geneid column because we want a full numeric data.frame
    
    
    the_heatmap = pheatmap(norm_counts, 
                       annotation_col=ph_annotation)
    

    return(the_heatmap)
    
}

top_expressed_heatmaps <- function(standard_deseq_output, n_genes = 35,padj_cut = 0.1){
    top_n_up <- select_n_top_sig(standard_deseq_output, n_genes = n_genes, 
                                 padj_cut = padj_cut,
                                 direction = 'up')
    
    top_n_down <- select_n_top_sig(standard_deseq_output, n_genes = n_genes, 
                                 padj_cut = padj_cut,
                                 direction = 'down')
    
    up_pheatmap <- generate_pheatmaps(standard_deseq_output, top_n_up)
    down_pheatmap <- generate_pheatmaps(standard_deseq_output, top_n_down)
    
    pheatmap_list = list(up_pheatmap, down_pheatmap)
    names(pheatmap_list) = c("up","down")
    
    return(pheatmap_list)
    
}