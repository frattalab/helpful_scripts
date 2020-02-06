
draw_sample_heatmap_from_deseq = function(deseq_object, normalization = "vst"){
# this is all essentially copy paste from the DESEQ2 vignette
    library(colorRamps)
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    if(normalization == "poisson"){
        library("PoiClaClu")
        poisd <- PoissonDistance(t(counts(deseq_object,normalized = TRUE)))
        
        samplePoisDistMatrix <- as.matrix( poisd$dd )
        rownames(samplePoisDistMatrix) <-colnames(deseq_object)
        colnames(samplePoisDistMatrix) <-colnames(deseq_object)
        
        simple_heat = pheatmap(samplePoisDistMatrix,
                 clustering_distance_rows = poisd$dd,
                 clustering_distance_cols = poisd$dd,
                 col = colors)
        return(simple_heat)
        
    }else if(normalization == "rlog"){
        transformed_dds <- rlog(deseq_object, blind = FALSE)
    }else if(normalization == "vst"){
        transformed_dds <- vst(deseq_object, blind = FALSE)
    }
        
        sampleDists <- dist(t(assay(transformed_dds)))

        
        sampleDistMatrix <- as.matrix( sampleDists )

        rownames(sampleDistMatrix) <- NULL

        
        colnames(sampleDistMatrix) <- colnames(deseq_object)

        simple_heat = pheatmap(sampleDistMatrix,
                 clustering_distance_rows = sampleDists,
                 clustering_distance_cols = sampleDists,
                 col = colors)
        
        
        
    
    return(simple_heat)
}