#' This function plots projection of daily confirmed cases assuming the last Rt
#' estimates (by default) stay the same over a defined period (14 days by default)
#'
#' The \code{plot_projection()} estimates posterior distribution using the particle filtering method
#' population for both sexes and incidence rate
#' @param pf Particle filtering resutls
#' @param dat data to compare (column plot)
#' @param days Reference Rt estimates, 1 means the latest estimate
#' @export
#' @import ggplot2
#' @examples
#' res <- plot_projection(); res$plot
case_projection_plot <- function(dat = NULL,
                                 smpl_daily_confirmed = NULL,
                                 proj = NULL,
                                 probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {

  df1 <-
    as.data.frame(t(apply(smpl_daily_confirmed, 1, function(x) quantile(x, probs))))
  df2 <-
    as.data.frame(t(apply(proj, 1, function(x) quantile(x, probs))))
  df2 <- rbind(df1[nrow(df1),], df2)
  df1$date <- dat$date
  df2$date <- seq(as.Date(max(dat$date)), by = "day", length.out = nrow(df2))

  plt <- ggplot(df1) +
    geom_col(data = dat, aes(x = date, y = daily_confirmed),
             fill = "grey50", alpha = 0.3) +
    geom_ribbon(aes(x = date, ymin = `2.5%`, ymax = `97.5%`),
                fill = "steelblue", alpha = 0.3) +
    geom_ribbon(aes(x = date, ymin = `25%`, ymax = `75%`),
                fill = "steelblue", alpha = 0.7) +
    geom_line(aes(x = date, y = `50%`), color = "steelblue", size = 1) +
    geom_ribbon(data = df2, aes(x = date, ymin = `2.5%`, ymax = `97.5%`),
                fill = "indianred", alpha = 0.3) +
    geom_ribbon(data = df2, aes(x = date, ymin = `25%`, ymax = `75%`),
                fill = "indianred", alpha = 0.7) +
    geom_line(data = df2, aes(x = date, y = `50%`), color = "indianred", size = 1) +
    labs(title = "", x = "", y = "") +
    scale_x_date(date_breaks = "1 month")

  out <- list(fit = df1, proj = df2, plot = plt)
}
