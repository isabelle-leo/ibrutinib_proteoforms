#' Retrieve peptide community memberships
#'
#' @param graphs \code{list} of graphs
#'
#' @return \code{data.frame}
#' 
#' @importFrom magrittr "%>%"
#' @importFrom igraph V
#' 
#' @export
#'
#' @examples
get_peptide_memberships_legacy <- function (graphs) {
  lapply(X = names(graphs),
         FUN = function (ioi) {
           g <- graphs[[ioi]]
           if (length(V(g)) == 0) return(NA)
           data.frame(ioi = ioi,
                      peptide = V(g)$name,
                      membership = V(g)$membership,
                      stringsAsFactors = FALSE)
         }) %>%
    .[!is.na(.)] %>%
    do.call(rbind, .)
}

get_peptide_memberships <- function (graphs) {
  processed_graphs <- lapply(names(graphs), function (ioi) {
    g <- graphs[[ioi]]
    if (!inherits(g, "igraph") || length(V(g)) == 0) return(NA)
    data.frame(ioi = ioi,
               peptide = V(g)$name,
               membership = V(g)$membership,
               stringsAsFactors = FALSE)
  })
  
  # Filter out NA values and empty data frames
  valid_indexes <- sapply(processed_graphs, function(df) {
    !is.na(df) && (is.data.frame(df) && nrow(df) > 0)
  })
  valid_graphs <- processed_graphs[valid_indexes]
  
  # Merge data frames
  do.call(rbind, valid_graphs)
}


