#' Generate color vector from a vector of values
#'
#' @param colors the color palette to use
#' @param vec the value vactor
#' @param gray gray tone for values exceeding the color palette
#'
#' @return the color vector
#' @export
#'
#' @examples
get_color_vector <- function (colors,
                              vec,
                              gray = "#494e57") {
  
  n <- length(unique(vec))
  
  suppl_needed <- max(n - length(colors), 0)
  col_suppl <- c(colors[1:min(n, length(colors))], rep(gray, suppl_needed))
  
  # count appearances of each group
  vec_accum <- table(vec) %>% sort(decreasing = TRUE)
  
  names(col_suppl) <- names(vec_accum) %>% as.character()
  
  # add a gray missing value (-1)
  col_suppl <- c(col_suppl, c("-1" = gray))
  
  return(col_suppl)
}
