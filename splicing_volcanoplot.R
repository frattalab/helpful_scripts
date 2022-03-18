splicing_dots_tables_function <- function(input_splicing, list_gene = "",psi_cutoff = 0,
                                          probability_changing_cutoff = 0) {

    if(psi_cutoff != 0 ){
      psi_genes <- input_splicing %>% 
        filter(mean_dpsi_per_lsv_junction >= psi_cutoff) %>% 
        filter(probability_changing >= 0.9)
    }
  
    if(probability_changing_cutoff != 0){
      psi_genes = psi_genes %>% 
        filter(probability_changing >= probability_changing_cutoff)
    }
    
    if(exists('psi_genes')){

      psi_genes = psi_genes %>% pull(gene_name)
      list_gene = unique(c(list_gene,psi_genes))
    }
    splicing_dots_tables <- input_splicing %>% 
      mutate(junction_name = case_when(gene_name %in% list_gene &
                                         probability_changing > 0.9  &
                                         mean_dpsi_per_lsv_junction > 0 ~ gene_name,
                                       T ~ "")) %>%
      mutate(`Novel Junction` = de_novo_junctions == 1) %>%
      mutate(log10_test_stat = -log10(1 - probability_changing)) %>%
      mutate(log10_test_stat = ifelse(is.infinite(log10_test_stat), 4.5, log10_test_stat)) %>%
      mutate(graph_alpha = ifelse(probability_changing > 0.9, 1, 0.2)) %>%
      mutate(label_junction = case_when(gene_name %in% list_gene &
                                          probability_changing > 0.9  &
                                          mean_dpsi_per_lsv_junction > 0 ~ junction_name,
                                        T ~ ""))


fig1 <- ggplot() +
  geom_point(data = splicing_dots_tables %>% dplyr::filter(de_novo_junctions == 0),
             aes(x = mean_dpsi_per_lsv_junction, y = log10_test_stat,
                 alpha = graph_alpha, fill = "Annotated Junction"), pch = 21, size = 2) +
  geom_point(data = splicing_dots_tables %>% dplyr::filter(de_novo_junctions != 0),
             aes(x = mean_dpsi_per_lsv_junction, y = log10_test_stat,
                 alpha = graph_alpha, fill = "Novel Junction"), pch = 21, size = 2) +
  geom_text_repel(data = splicing_dots_tables[19 != ""],
                  aes(x = mean_dpsi_per_lsv_junction, y =log10_test_stat,
                      label = label_junction,
                      color = as.character(de_novo_junctions)), point.padding = 0.3,
                  nudge_y = 0.2,
                  min.segment.length = 0,
                  box.padding  = 2,
                  max.overlaps = Inf,
                  size=4, show.legend = F) +
  geom_hline(yintercept = -log10(1 - .9)) +
  geom_vline(xintercept = -0,linetype="dotted") +
  scale_fill_manual(name = "",
                    breaks = c("Annotated Junction", "Novel Junction"),
                    values = c("Annotated Junction" = "#648FFF", "Novel Junction" = "#fe6101") ) +
  scale_color_manual(name = "",
                     breaks = c("0", "1"),
                     values = c("0" = "#648FFF", "1" = "#fe6101") ) +
  guides(alpha = "none", size = "none") +
  theme(legend.position = 'top') +
  ggpubr::theme_pubr() +
  xlab("delta PSI") +
  ylab(expression(paste("-Lo", g[10], " Test Statistic"))) +
  theme(text = element_text(size = 24)) +
  theme(legend.text=element_text(size=22)) +
  scale_x_continuous(labels = scales::percent, limits = c(-1,1))
  return(fig1)
}
