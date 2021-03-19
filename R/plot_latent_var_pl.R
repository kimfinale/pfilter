#' This function plots latent variables
#'
#' The \code{plot_latent_var_pl()} estimates posterior distribution using the particle filtering method
#' population for both sexes and incidence rate
#' @param sim_res Simulation results
#' @param data data to compare
#' @param data_type data based on the date of infection, symptom onset, or confirmation
#' @export
#' @import tidyverse
#' @examples
#' plot_latent_var_pl(sim_res = pf)
plot_latent_var_pl <- function(sim_res = pf_result,
                               data = Rt_data,
                               data_type = c("infection", "symptom onset", "confirmation")) {
  dtype <- match.arg(data_type)
  col_index <- 8 #dR, i.e., daily confirmation
  var <- "CE"
  col_nm <- "daily_infect"

  if (dtype == "infection") {
    var <- "CE"
    col_index <- 6
    col_nm <- "daily_infect"
  }
  else if (dtype == "symptom onset") {
    col_index <- 7
    var <- "CI"
    col_nm <- "daily_symptom"
  }
  else if (dtype == "confirmation") {
    col_index <- 8
    var <- "CR"
    col_nm <- "daily_confirm"
  }
  dd <- r[col_index, ]
  ddf <- data.frame(matrix(unlist(dd), ncol = length(dd), byrow = F))
  pr <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  dR_quantile <- as.data.frame(t(apply(ddf, 1, function(x) {quantile(x, pr)})))
  dR_mean <- data.frame(mean = rowMeans(ddf))
  df <- cbind(dR_quantile, dR_mean, data[, c("t", col_nm)])
  suppressMessages(library(tidyverse))
  df$date <- seq(as.Date("2020-01-20"), length.out = nrow(df), by = "day")

  plt <- ggplot(df, aes(x = date)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "steelblue", alpha = 0.3) +
    geom_ribbon(aes(ymin = `25%`, ymax = `75%`), fill = "steelblue", alpha = 0.6) +
    geom_line(aes(y = `50%`), color = "steelblue", size = 1, linetype = "dotted") +
    geom_line(aes(y = mean), color = "steelblue", size = 1) +
    geom_point(aes(y = eval(parse(text = col_nm))), color = "darkred", size = 1) +
    labs(title = paste0("Daily ", dtype, " via particle filtering with backward sampling"), y = paste0("Daily ", dtype), x = "") +
    scale_x_date(date_breaks = "2 month", date_minor_breaks = "1 month", limits = c(as.Date("2020-01-20"), max(df$date))) +
    theme(axis.text.x = element_text(angle = 0, vjust = 1.0, hjust = 0.5))

  return (plt)

}
