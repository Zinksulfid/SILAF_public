#############################
## GSEA plotting functions ##
## by Ilian Atanassov #######
#############################

## I used this post as a basis https://stackoverflow.com/questions/35942317/how-to-plot-two-geoms-side-by-side-for-same-categorical-in-r-errorbarjitter

## ---- functions ----

kegg_ids_tab_to_long <- function(kegg_id_plus_entrez, kegg_id){
  
  gsub(paste0(kegg_id, ".*", kegg_id, "\t"), "", kegg_id_plus_entrez)
  
}

add_kegg_id <- function(entrez_ids, kegg_id){
  
  entrez_ids %>% 
    mutate(kegg_id = kegg_id)
}

## ---- color scheme -----

# 11-class RdBu color scheme https://colorbrewer2.org/#type=diverging&scheme=RdBu&n=11, data range from -1 to 1, 102 ranges and colors, gets unique colors.

breaks_ <- c(-Inf ,seq(from = -1, to = 1, by = 0.02), Inf)

cut_ranges_for_color <- cut(seq(from = -10, to = 10, by = 0.001), 
                            breaks = breaks_, 
                            include.lowest=TRUE, 
                            right=FALSE)


color_scheme <- c(colorRampPalette(colors = c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7"))(length(unique(cut_ranges_for_color))/2),
                  colorRampPalette(colors = c("#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"))(length(unique(cut_ranges_for_color))/2))

tibble(logFC_cut = unique(cut_ranges_for_color),
       colors_ = color_scheme) -> logFC_cut


plot_gsea <- function(entrez_gene_kegg_annot = NULL, gsea_result = NULL, quant_data = NULL, FDR_set = 0.05, pValue = NULL){
  
  quant_data <- read_tsv(quant_data)
  
  
  if(is.null(pValue)){
  gsea_result <- read_tsv(gsea_result) %>%
    filter(FDR < FDR_set)} else {
  gsea_result <- read_tsv(gsea_result) %>%
    filter(pValue < 0.05)    
    }
  
  entrez_gene_kegg_annot <- read_lines(entrez_gene_kegg_annot)

    ## ---- wrangle ----
    
    ## entrez id to kegg set mapping
    
    entrez_gene_kegg_annot %>%
    map(str_split,pattern = "\\\t") %>% 
    map(unlist) %>% 
    map(1) -> kegg_ids_list
    
    entrez_gene_kegg_annot %>%
    map2(.x = .,
         .y = kegg_ids_list,
         .f = kegg_ids_tab_to_long) %>% 
    map(as_tibble) %>% 
    map(separate_rows, value, sep = "\\\t") %>% 
    map2(.x = .,
         .y = kegg_ids_list,
         .f =  add_kegg_id) %>% 
    bind_rows() %>% 
    dplyr::rename(geneSet = kegg_id,
                  Id = value) %>% 
    dplyr::select(geneSet, Id) -> entrez_gene_kegg_annot_long
    
    entrez_gene_kegg_annot_long %>% 
    group_by(geneSet) %>% 
    dplyr::summarise(SizeID = paste(Id, collapse = ";")) -> geneSet_SizeID
    
    ## add entrez ids corresponding to size of gsea results
    
    gsea_result %>% 
    mutate(FDR = format(FDR,digits = 2,scientific = TRUE)) %>% 
    left_join(geneSet_SizeID) -> gsea_result_final
    
    quant_data %>% 
    dplyr::rename(logFC = score) %>% 
    arrange(logFC) %>% 
    mutate(logFC_rank = 1:nrow(.)) -> quant_data_final
    
    gsea_result_final %>% 
    left_join(geneSet_SizeID) %>% 
    dplyr::select(description, SizeID) %>%  
    separate_rows(SizeID , sep = ";",convert = TRUE) %>% 
    dplyr::rename(entrezgene = SizeID ) %>% 
    left_join(quant_data_final[,c("entrezgene", "logFC", "logFC_rank")]) %>% 
    filter(!is.na(logFC)) -> gsea_data_long
    
    ## ---- y axis pathway order using normalizedEnrichmentScore ----
    
    gsea_result_final %>% 
    arrange(normalizedEnrichmentScore) %>% 
    dplyr::select(description) %>% 
    unlist(use.names = FALSE) %>% 
    rev()-> pathway_order
    
    ## ---- rank vs logFC plot----
    
    quant_data_final %>% 
    mutate(logFC_cut = cut(logFC,
                           breaks = breaks_,
                           include.lowest=TRUE,
                           right=FALSE)) %>% 
    left_join(logFC_cut) %>%
    ggplot(aes(x = logFC_rank, y = logFC, fill = colors_, color = colors_))+
    geom_bar(stat = "identity")+
    scale_fill_identity()+
    scale_color_identity()+
    scale_x_continuous(expand = c(0, 0))+
    theme_classic()+
    theme(axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = 0.5))+
    labs(x = NULL, title = "KEGG")-> gsea_rank_vs_logFC
    
    ## ---- prepare data for plot ----
    
    quant_data_final %>% 
    mutate(logFC_cut = cut(logFC,
                           breaks = breaks_,
                           include.lowest=TRUE,
                           right=FALSE)) %>% 
    left_join(logFC_cut) %>% 
    right_join(unique(gsea_data_long[,c("entrezgene", "description")])) %>%
    mutate(description = factor(description, levels = pathway_order)) %>% 
    mutate(description_num = as.numeric(description)) %>% 
    mutate(description_num_y_min = description_num-0.2) %>% 
    mutate(description_num_y_max = description_num+0.2) -> data_for_plot
    
    ## ---- table with description and description position number ----
    
    data_for_plot %>% 
    dplyr::select(description,description_num) %>% 
    arrange(description_num) %>% 
    unique() -> y_axis_annot
    
    ## ---- plot results ----
    
    data_for_plot %>% 
    ggplot(aes(x = logFC_rank, y = description_num, color = colors_, color = colors_ ))+
    geom_linerange(aes(ymin = description_num_y_min, ymax = description_num_y_max))+
    scale_fill_identity()+
    scale_color_identity()+
    scale_x_continuous(expand = c(0, 0))+
    scale_y_continuous(breaks = y_axis_annot$description_num,labels = as.character(y_axis_annot$description))+
    theme_classic()+
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank())+
    labs(x = NULL, y = NULL) -> gsea_categories

    p <- plot_grid(gsea_rank_vs_logFC, gsea_categories, ncol = 1, align = "v",rel_heights = c(1,2))
    return(p)
}
