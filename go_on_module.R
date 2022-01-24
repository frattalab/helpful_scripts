# ------Probabaly should put here but this function takes a cemitool object and 
go_on_module = function(cem_object, species = "hsa", module = "M1", de_dataframe = ""){
    module_table = module_genes(cem_object)
    if(!(species %in% c("hsa","mmu"))){
        print("Species has to be 'hsa' or 'mmu'")
        print(paste("Species you entered:",species))
        return(1)
    }
    if(species == "hsa"){
        gene_mart_map = annotables::grch38
        
    }else if(species == "mmu"){
        gene_mart_map = annotables::grcm38
    }

        module_table = module_table %>% left_join(gene_mart_map %>% dplyr::select(entrez ,symbol),
                                              by = c("genes" = "symbol")) %>% 
        as.data.table()
    
    mgenes = module_table[modules == module & !is.na(entrez),entrez]

    if(de_dataframe != ""){
        
        mgenes_fold = de_dataframe[ncbi_gene_id %in% mgenes,log2FoldChange]
        
        names(mgenes_fold) = de_dataframe[ncbi_gene_id %in% mgenes,ncbi_gene_id]
        
    }else(
        mgenes_fold = "No enrichment provided"
    )

    if(species == "hsa"){
        goenrich = clusterProfiler::enrichGO(mgenes,OrgDb = 'org.Hs.eg.db',ont = "ALL")
        
        go_module = clusterProfiler::setReadable(goenrich,'org.Hs.eg.db', 'ENTREZID')
        
    }else if(species == "mmu"){
        goenrich = clusterProfiler::enrichGO(mgenes,OrgDb = 'org.Mm.eg.db',ont = "ALL")
        
        go_module = clusterProfiler::setReadable(goenrich,'org.Mm.eg.db', 'ENTREZID')
    }
    temp = list(go_module,mgenes_fold)
    names(temp) = c("go_module","foldChange")
    return(temp)
}


get_mapped_module = function(cem_object, species = "hsa", module = "M1"){
    module_table = module_genes(cem_object)
    if(!(species %in% c("hsa","mmu"))){
        print("Species has to be 'hsa' or 'mmu'")
        print(paste("Species you entered:",species))
        return(1)
    }
    if(species == "hsa"){
        gene_mart_map = as.data.table(clean_names(fread("/Users/annaleigh/Documents/data/reference_genomes/human_gene_mart_export.txt")))
        
    }else if(species == "mmu"){
        gene_mart_map = as.data.table(clean_names(fread("/Users/annaleigh/Documents/data/reference_genomes/mouse_gene_mart_export.txt")))
    }
    module_table = module_table %>% left_join(gene_mart_map,
                                              by = c("genes" = "gene_stable_id_version")) %>% 
        as.data.table()
    mgenes = module_table[modules == module & !is.na(ncbi_gene_id),ncbi_gene_id]
    module_table = module_table[modules == module]
    return(module_table)
}
