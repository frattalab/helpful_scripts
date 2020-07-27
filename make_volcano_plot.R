

make_volcano_plot = function(results_table, logFoldCut = 7, padjCut = 0.1,nonSigAlpha = 0.3){
    #make a plot
    results_edit = copy(results_table)
    results_edit = results_edit %>% 
        mutate(plotlabel = ifelse(abs(log2FoldChange) > logFoldCut & padj < padjCut,gene_name,"")) %>% 
        mutate(plotAlpha = ifelse( padj < 0.1,1,nonSigAlpha)) 
    
    
    volc = ggplot(results_edit,aes(x = log2FoldChange, y = -log10(padj))) + 
        geom_point(aes(alpha = plotAlpha),show_guide  = F ) + 
        geom_hline(yintercept = -log10(0.1),linetype="dotted") + 
        geom_text_repel(aes(label = plotlabel)) + 
        ggpubr::theme_classic2() + 
        ylab(bquote('-Log'[10]~ 'Adjusted p-value')) + 
        xlab(bquote('Log'[2]~ 'Fold Change')) 
    
        return(volc)
}
