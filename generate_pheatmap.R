library("pheatmap")



generate_pheatmaps <-  function(standard_deseq_output, n_genes = 35){
    normalized_deseq_obj <- normTransform(standard_deseq_output$deseq_obj) #normalize the counts
    
    #we're taking the colData from the DESeq object (which is metadata)
    
    
    ph_annotation <- as.data.frame(colData(standard_deseq_output$deseq_obj)) %>% 
        select(comparison) #if you just check what this: colData(one_hour_dds$deseq_obj) you see that we have 2 additional columns which we don't want on a heatmap
    
    #I'm going to select the geneids for the top n upregulated genes
    top_n_up <- standard_deseq_output$results_table %>% 
        filter(padj < 0.1) %>% 
        slice_max(log2FoldChange, n = n_genes) %>% 
        pull(Geneid)
    
    #I'm going to select the geneids for the bottom n upregulated genes
    bottom_n_down <- standard_deseq_output$results_table %>% 
        filter(padj < 0.1) %>% 
        slice_min(log2FoldChange, n = n_genes) %>% 
        pull(Geneid)
    
    #I'm going to select the geneids for the top n upregulated genes
    norm_up_counts <- assay(normalized_deseq_obj) %>% #Normalized counts are stored in the "assay" slot on this object
        as.data.frame() %>% #it's not a data frame so I have to make it one
        rownames_to_column('Geneid') %>% #our old friend
        filter(Geneid %in% top_n_up) %>% #select these genes
        left_join(standard_deseq_output$results_table %>% #left join from the results table to get the gene_name
                      select(Geneid,gene_name)) %>% #notice that i'm selecting inside the left_join
        mutate(gene_name = make.unique(gene_name)) %>% #this is to save me from having to do this incase the gene_name are like repeated
        column_to_rownames('gene_name') %>% #the heatmap function pheatmap uses rownames for plotting
        dplyr::select(-Geneid) #get rid of Geneid column because we want a full numeric data.frame
    
    
    #I'm going to select the geneids for the top n upregulated genes
    norm_down_counts <- assay(normalized_deseq_obj) %>% #Normalized counts are stored in the "assay" slot on this object
        as.data.frame() %>% #it's not a data frame so I have to make it one
        rownames_to_column('Geneid') %>% #our old friend
        filter(Geneid %in% bottom_n_down) %>% #select these genes
        left_join(standard_deseq_output$results_table %>% #left join from the results table to get the gene_name
                      select(Geneid,gene_name)) %>% #notice that i'm selecting inside the left_join
        mutate(gene_name = make.unique(gene_name)) %>% #this is to save me from having to do this incase the gene_name are like repeated
        column_to_rownames('gene_name') %>% #the heatmap function pheatmap uses rownames for plotting
        dplyr::select(-Geneid) #get rid of Geneid column because we want a full numeric data.frame
    
    heatmap_up = pheatmap(norm_up_counts, 
                    annotation_col=ph_annotation)
    
    heatmap_down = pheatmap(norm_down_counts, 
                       annotation_col=ph_annotation)
    
    
    return_list = list(heatmap_up, heatmap_down)
    
    names(return_list) = c("upregulated","downregulated")

    return(return_list)
    
}