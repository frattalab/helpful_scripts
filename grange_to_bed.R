grange_to_bed = function(gr, write_spot, names_col = "",scores_col = ""){
    # allow there to either be empty names or if you fill someting out
    if(names_col == ""){
        write_names = c(rep(".", length(gr)))
        # then the name is taken fromt he metadata
    }else(
        write_names = gr@elementMetadata %>% 
            as.data.frame() %>% dplyr::select(!!names_col)
    )
    # similar concept with scores
    # allow there to either be empty names or if you fill someting out
    if(scores_col == ""){
        write_scores = c(rep(".", length(gr)))
        # then the name is taken fromt he metadata
    }else(
        write_scores = gr@elementMetadata %>% 
            as.data.frame() %>% dplyr::select(!!scores_col)
    )
    df <- data.frame(seqnames=seqnames(gr),
                     starts=start(gr)-1,
                     ends=end(gr),
                     names=write_names,
                     scores=write_scores,
                     strands=strand(gr))
    write.table(df, file=write_spot, quote=F, sep="\t", row.names=F, col.names=F)
    return(df)
}