plot_id_symbol_melting <- function (id_symbol,
                            additional_peptides = NULL,
                              e_set,
                              x_col = "temperature",
                              facet_col = NULL,
                              dot_color = 'black',
                              add_splines = FALSE,
                              spline_color = 'steelblue',
                              hide_lines = TRUE,
                              x_label = "Temperature (C)",
                              y_label = "Fraction Non-Denatured",
                              y_limits = c(-0.2, 1.2),
                              custom_theme = NULL, 
                            label_proteins = FALSE,
                            highlight_additional = FALSE,
                            subset_cell_line = NULL) {
  
  library(RColorBrewer)
  
  logical_info <- grepl(id_symbol, fData(e_set)$id) | grepl(additional_peptides, fData(e_set)$Peptide_AA)
  id_e_set <- e_set[logical_info, ]
  
  if (nrow(id_e_set) == 0) stop("No data found for this gene symbol of interest.")
  
  data <- id_e_set %>%
    tidy_e_set() %>%
    set_names(gsub("^gene$", "id_symbol", colnames(.))) %>%
    left_join(fData(id_e_set), by = "id")
  
  data$x_col <- data[[x_col]]
  
  if(length(subset_cell_line) > 0){
    
    data <- data[unique(data$sample_name) %in% subset_cell_line,]
    
  }
  
  p <- data %>%
    ggplot(aes(x = x_col,
               y = value,
               group = interaction(peptide.x, sample_name))) +
    ylim(y_limits)
  
  if (!hide_lines) {
    p <- p + geom_line(alpha = .7)
  }
  
  if (add_splines & !highlight_additional) {
    p <- p + geom_smooth(method = "lm",
                         formula = 'y ~ splines::ns(x, df = 4)',
                         se = FALSE,
                         alpha = .4,
                         size = .7,
                         color = spline_color)
  }
  
  if (!is.null(facet_col)) {
    p <- p + facet_wrap(facets = facet_col)
  }
  
  if (!is.null(custom_theme)) {
    p <- p + custom_theme()
  }
  if (label_proteins) {
    p <- p + geom_point(aes(color = protein_ids.x))
  } else {
    p <- p + geom_point(color = dot_color)
  }
  if (!highlight_additional & label_proteins) {
    p <- p + aes(color = protein_ids.x)
  } 
  if (highlight_additional){
    color_names <- data$protein_ids.x
    colors <- ifelse(data$Peptide_AA.y %in% additional_peptides, "Red", "Grey30")
    names(colors) <- color_names
    sizes <- ifelse(data$Peptide_AA.y %in% additional_peptides, 1.5, .5)
    alphas <- ifelse(data$Peptide_AA.y %in% additional_peptides, .9, .1)
    names(sizes) <- color_names
    names(alphas) <- color_names
    p <- p + aes(color = protein_ids.x, size = sizes, alpha = alphas) + scale_color_manual(name = "Peptide", values = colors)
  }
  if (add_splines & highlight_additional) {
    p <- p + geom_line(stat = "smooth",
                       method = "lm",
                         formula = 'y ~ splines::ns(x, df = 4)',
                         se = FALSE,
                         alpha = .3,
                         size = 1) 
    
    
  }
  
  p <- p +
    xlab(as.character(x_label)) + ylab(as.character(y_label)) +
    ggtitle(as.character(id_symbol)) +
    theme(legend.position = "bottom")
  
  print(p)
  
  return(p)
}
