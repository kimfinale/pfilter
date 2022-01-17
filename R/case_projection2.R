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
                             smpl_last_states = NULL,
                             smpl_Rt = NULL,
                             smpl_daily_confirmed = NULL,
                             Rt_sim = NULL,
                             tend = 14,
                             dt = 0.1){

  # library(pfilter)
  if (is.null(dat)) {
    dat <- data.table::fread("daily_sim/dat.csv")
  }
  if (is.null(smpl_last_states)) {
    smpl_last_states <- data.table::fread("daily_sim/smpl_last_states.csv")
  }
  if (is.null(smpl_Rt)) {
    smpl_Rt <- data.table::fread("daily_sim/smpl_Rt.csv")
    smpl_Rt <- unlist(smpl_Rt[nrow(smpl_Rt),])
  }
  if (is.null(smpl_daily_confirmed)) {
    smpl_daily_confirmed <- data.table::fread("daily_sim/smpl_daily_confirmed.csv")
  }
  # add 1 to put the last data points
  proj_confirm <-
    data.frame(matrix(NA, nrow = tend + 1, ncol = ncol(smpl_daily_confirmed)))
  proj_confirm[1, ] <- tail(smpl_daily_confirmed, 1)

  # library(dplyr)
  y <- as.data.frame(smpl_last_states)
  nms <- names(y)

  if (!is.null(Rt_sim)) {
    smpl_Rt <- (Rt_sim / mean(smpl_Rt)) * smpl_Rt
  }
  # beta <- as.matrix(smpl_Rt / R0_dur(params = theta), nrow = length(smpl_Rt))
  beta <- smpl_Rt / R0_dur(params = theta)

  for (i in 1:tend) {
    y <- process_model(params = theta,
                       y = y,
                       tbegin = 0,
                       tend = 1,
                       dt = dt,
                       beta = beta,
                       stoch = FALSE)
    names(y) <- nms
    proj_confirm[i+1, ] <- y[,"CR"]
  }

  pr <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  smpl_daily_confirmed_qt <-
    as.data.frame(t(apply(smpl_daily_confirmed, 1, function(x) quantile(x, pr))))
  proj_confirm_qt <-
    as.data.frame(t(apply(proj_confirm, 1, function(x) quantile(x, pr))))

  sub1 <- data.frame(matrix(NA,
                            nrow = nrow(proj_confirm_qt) - 1,
                            ncol = ncol(proj_confirm_qt)))

  names(sub1) <- names(smpl_daily_confirmed_qt)

  df1 <- rbind(smpl_daily_confirmed_qt, sub1)
  # -1 to match the end of the fit
  sub2 <- data.frame(matrix(NA,
                            nrow = nrow(smpl_daily_confirmed_qt) - 1,
                            ncol = ncol(smpl_daily_confirmed_qt)))
  names(sub2) <- names(smpl_daily_confirmed_qt)

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
