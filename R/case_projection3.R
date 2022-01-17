#' This function plots projection of daily confirmed cases assuming the last Rt
#' estimates (by default) stay the same over a defined period (14 days by default)
#'
#' The \code{plot_projection()} estimates posterior distribution using the particle filtering method
#' population for both sexes and incidence rate
#' @param smpl_last_states PF results of the state variables at t=t_N
#' @param smpl_Rt PF results for the Rt at t = t_N
#' @param Rt_sim Reference Rt estimates, 1 means the latest estimate
#' @param nsample Reference Rt estimates, 1 means the latest estimate
#' @param tend Reference Rt estimates, 1 means the latest estimate
#' @param dt Reference Rt estimates, 1 means the latest estimate
#' @param probs Reference Rt estimates, 1 means the latest estimate
#' @export
#' @import ggplot2
#' @examples
#' res <- plot_projection(); res$plot
project_case <- function(smpl_last_states = NULL,
                         smpl_Rt = NULL,
                         Rt_sim = NULL,
                         nsample = 200,
                         tend = 14,
                         dt = 0.1,
                         probs = c(0.025, 0.25, 0.5, 0.75, 0.975)){

  if (is.null(smpl_last_states)) {
    smpl_last_states <- data.table::fread("daily_sim/smpl_last_states.csv")
  }
  if (is.null(smpl_Rt)) {
    smpl_Rt <- data.table::fread("daily_sim/smpl_Rt.csv")
    smpl_Rt <- unlist(smpl_Rt[nrow(smpl_Rt),])
  }

  proj_confirm <-
    data.frame(matrix(NA, nrow = tend, ncol = 200))

  y <- as.data.frame(smpl_last_states)
  nms <- names(y)

  if (!is.null(Rt_sim)) {
    smpl_Rt <- (Rt_sim / mean(smpl_Rt)) * smpl_Rt
  }
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
    proj_confirm[i, ] <- y[,"CR"]
  }

  proj_confirm_qt <-
    as.data.frame(t(apply(proj_confirm, 1, function(x) quantile(x, probs))))

  return (proj_confirm_qt)
}
