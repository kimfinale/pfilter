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
case_projection2 <- function(dat = NULL,
                             smpl_states = NULL,
                             smpl_Rt = NULL,
                             smpl_daily_confirmed = NULL,
                             Rt_sim = NULL,
                             tend = 14,
                             dt = 0.1){

  # library(pfilter)

  fit_confirm <- smpl_daily_confirmed
  # add 1 to put the last data points
  proj_confirm <- data.frame(matrix(NA, nrow = tend + 1, ncol = nsample))
  proj_confirm[1, ] <- tail(fit_confirm, 1)

  # library(dplyr)

  y <- smpl_states

  if(!is.null(Rt_sim)){
    samp_Rt <- (Rt_sim / mean(samp_Rt)) * samp_Rt
  }

  for (i in 1:tend) {
    y <- process_model(params = theta,
                       y = y,
                       tbegin = 0,
                       tend = 1,
                       dt = dt,
                       beta = samp_Rt / R0_dur(params = theta),
                       stoch = FALSE)
    names(y) <- str[1:nstates]
    proj_confirm[i+1, ] <- y[,"CR"]
  }

  pr <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  fit_confirm_qt <- as.data.frame(t(apply(fit_confirm, 1, function(x) quantile(x, pr))))
  proj_confirm_qt <- as.data.frame(t(apply(proj_confirm, 1, function(x) quantile(x, pr))))

  sub1 <- data.frame(matrix(NA,
                            nrow = nrow(proj_confirm_qt) - 1,
                            ncol = ncol(proj_confirm_qt)))

  names(sub1) <- names(fit_confirm_qt)

  df1 <- rbind(fit_confirm_qt, sub1)
  # -1 to match the end of the fit
  sub2 <- data.frame(matrix(NA,
                            nrow = nrow(fit_confirm_qt) - 1,
                            ncol = ncol(fit_confirm_qt)))
  names(sub2) <- names(fit_confirm_qt)

  df2 <- rbind(sub2, proj_confirm_qt)

  df1$date <- seq(min(dat$date), by = "day", length.out = nrow(df1))
  df2$date <- df1$date

  library(ggplot2)
  plt <- ggplot(df1, aes(x = date)) +
    geom_col(data = dat, aes(x = date, y = daily_confirmed),
             fill = "grey50", alpha = 0.3) +
    geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`),
                fill = "steelblue", alpha = 0.3) +
    geom_ribbon(aes(ymin = `25%`, ymax = `75%`),
                fill = "steelblue", alpha = 0.7) +
    geom_line(aes(y = `50%`), color = "steelblue", size = 1) +
    geom_ribbon(data = df2, aes(x = date, ymin = `2.5%`, ymax = `97.5%`),
                fill = "indianred", alpha = 0.3) +
    geom_ribbon(data = df2, aes(x = date, ymin = `25%`, ymax = `75%`),
                fill = "indianred", alpha = 0.7) +
    geom_line(data = df2, aes(x = date, y = `50%`), color = "indianred", size = 1) +
    labs(title = "Daily confirmed case", x = "", y = "") +
    scale_x_date(date_breaks = "1 month")

  out <- list(fit = df1, projection = df2, plot = plt)

  return (out)
}
