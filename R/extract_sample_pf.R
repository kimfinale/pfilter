#' Estimate parameter values using the maximum likelihood approach
#'
#' The \code{maxlik()} estimates the parameter values via the maximum likelihood approach
#' @param pf Output from a particle filter
#' @param str Variable name characters
#' @param days Days, starting from the end, over which values of the variable is averaged
#' @param ids Particle IDs. It's here to consistently select the same particles
#' @export
#' @examples
#' sample <- extract_sample_pf(pf, str = "Rt", days = 1, ids = sample(1:npart, 200))
extract_sample_pf <- function(pf, str = "Rt", days = 1, ids = NULL){
  df <- as.data.frame(sapply(pf, function(x) x[, str]))
  nr <- nrow(df)
  df_sub <- df[(nr - days + 1):nr, ids]
  if(days > 1) {
    vec <- colMeans(df_sub)
  }
  else{
    vec <- unlist(df_sub)
  }
  return (vec)
}

