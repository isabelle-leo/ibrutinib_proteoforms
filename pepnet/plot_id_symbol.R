plot_id_symbol <- function (id_symbol,
                            additional_peptides,
                              e_set,
                              x_col = "sample_name",
                              facet_col = NULL,
                              dot_color = 'black',
                              add_splines = FALSE,
                              spline_color = 'steelblue',
                              hide_lines = FALSE,
                              x_label = "Sample name",
                              y_label = "log2 protein level",
                              y_limits = c(-0.2, 1.2),
                              custom_theme = NULL) {
  
  logical_info <- grepl(id_symbol, fData(e_set)$id) | grepl(additional_peptides, fData(e_set)$Peptide_AA)
  id_e_set <- e_set[logical_info, ]
  
  if (nrow(id_e_set) == 0) stop("No data found for this gene symbol of interest.")
  
  data <- id_e_set %>%
    tidy_e_set() %>%
    set_names(gsub("^gene$", "id", colnames(.))) %>%
    left_join(fData(id_e_set), by = "id")
  
  data$x_col <- data[[x_col]]
  
  p <- data %>%
    ggplot(aes(x = x_col,
               y = value,
               group = interaction(id, sample_name))) +
    geom_point(color = dot_color) +
    ylim(y_limits)
  
  if (!hide_lines) {
    p <- p + geom_line(alpha = .7)
  }
  
  if (add_splines) {
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
  
  p <- p +
    xlab(x_label) + ylab(y_label) +
    ggtitle(id) +
    theme(legend.position = "bottom")
  
  print(p)
  
  return(p)
}
