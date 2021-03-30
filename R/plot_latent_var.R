#' This function plots latent variables
#'
#' The \code{plot_latent_var()} estimates posterior distribution using the particle filtering method
#' population for both sexes and incidence rate
#' @param sim_res Simulation results
#' @param data data to compare
#' @param data_type data based on the date of infection, symptom onset, or confirmation
#' @export
#' @import tidyverse
#' @examples
#' pf <- particle_filter(data_type = "infection", npart = 1e3, dt = 0.1); plot_latent_var(sim_res = pf)
plot_latent_var <- function(sim_res = pf_result,
                            data = Rt_data,
                            data_type = c("infection", "symptom onset", "confirmation")) {
  dtype <- match.arg(data_type)
  var <- "CE"
  col_nm <- "daily_infected"

  if (dtype == "infection") {
    var <- "CE"
    col_nm <- "daily_infected"
  }
  else if (dtype == "symptom onset") {
    var <- "CI"
    col_nm <- "daily_symptom_onset"
  }
  else if (dtype == "confirmation") {
    var <- "CR"
    col_nm <- "daily_confirmed"
  }

  var_quantile <- as.data.frame(t(apply(sim_res$latent_var_filtered[, , var], 2, function(x) {quantile(x, pr)})))
  var_mean <- data.frame(mean = rowMeans(t(sim_res$latent_var_filtered[, , var])))
  df <- cbind(var_quantile, var_mean, data[, c("t", col_nm)])
  df$date <- seq(as.Date("2020-01-20"), length.out = nrow(df), by = "day")
  suppressMessages(library(tidyverse))
  plt <- ggplot(df, aes(x = date)) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), fill = "steelblue", alpha = 0.3) +
    geom_ribbon(aes(ymin = `25%`, ymax = `75%`), fill = "steelblue", alpha = 0.6) +
    geom_line(aes(y = `50%`), color = "steelblue", size = 1, linetype = "dotted") +
    geom_line(aes(y = mean), color = "steelblue", size = 1) +
    geom_point(aes(y = eval(parse(text = col_nm))), color = "darkred", size = 1) +
    labs(title = paste0("Daily ", dtype, " via particle filtering"), y = paste0("Daily ", dtype), x = "") +
    scale_x_date(date_breaks = "2 month", date_minor_breaks = "1 month", limits = c(as.Date("2020-01-20"), max(df$date)))

  return (plt)
}
