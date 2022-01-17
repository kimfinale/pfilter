#' Models the process (infection transmission, recovery, etc.).
#' State variables and beta are given as a vector in which the number of elements are the number of particles
#'
#' The \code{process_model()} does SE1E2IR model integration using Euler method
#' @param tstart start time for simulation with default value of 0
#' @param tend end time with default value of 1 (i.e., one day process of SEI1I2R)
#' @param dt A step size
#' @param params Parameters
#' @param y A vector of a vector of state variables
#' @param beta A vector of beta
#' @param VOC If TRUE, variants of concern (VOC) are also included
#' @export
#' @import data.table
#' @examples
#' pop <- process_model_variant(params, y, tbegin = 0, tend = 1,)
# Process model for simulation --------------------------------------------
process_model_variant <- function (t = NULL,
                                   params = NULL,
                                   y = NULL,
                                   tbegin = 0,
                                   tend = 1,
                                   dt = 0.2,
                                   beta = NULL) {

  # daily_* values are updated every process_model_step
  y1 <- y[["Alpha"]]
  y1[, c("CE", "CI", "CR")] <- 0
  S <- y1[, "S"]
  E <- y1[, "E"]
  P <- y1[, "P"]
  A <- y1[, "A"]
  I <- y1[, "I"]
  R <- y1[, "R"]
  daily_infected <- y1[, "CE"]
  daily_symptom_onset <- y1[, "CI"]
  daily_confirmed <- y1[, "CR"]
  pop <- S + E + P + A + I + R
  durP <- (1 / params[["delta"]] - 1 / params[["epsilon"]]) # residence time in P
  zeta <- 1 / durP # rate from P to I
  beta_alpha <- beta[["Alpha"]]

  delta_start <- params[["Delta_start"]]
  omicron_start <- params[["Omicron_start"]]

  if (t > delta_start) { # make the condition for t
    y2 <- y[["Delta"]]
    y2[, c("CE", "CI", "CR")] <- 0
    E_delta <- y2[, "E"] # no S for the variant model
    P_delta <- y2[, "P"]
    A_delta <- y2[, "A"]
    I_delta <- y2[, "I"]
    R_delta <- y2[, "R"]
    daily_infected_delta <- y2[, "CE"]
    daily_symptom_onset_delta <- y2[, "CI"]
    daily_confirmed_delta <- y2[, "CR"]
    pop_delta <- E_delta + P_delta + A_delta + I_delta + R_delta
    pop <- pop + pop_delta
    beta_delta <- beta[["Delta"]]
  }
  if (t > omicron_start) { # make the condition for t
      y3 <- y[["Omicron"]]
      y3[, c("CE", "CI", "CR")] <- 0
      E_omicron <- y3[, "E"] # no S for the variant model
      P_omicron <- y3[, "P"]
      A_omicron <- y3[, "A"]
      I_omicron <- y3[, "I"]
      R_omicron <- y3[, "R"]
      daily_infected_omicron <- y3[, "CE"]
      daily_symptom_onset_omicron <- y3[, "CI"]
      daily_confirmed_omicron <- y3[, "CR"]
      pop_omicron <- E_omicron + P_omicron + A_omicron + I_omicron + R_omicron
      pop <- pop + pop_omicron
      beta_omicron <- beta[["Omicron"]]
  }

  for (ii in seq((tbegin + dt), tend, dt)) {
    foi <- beta_alpha * (params[["bp"]] * P + params[["ba"]] * A + I)
    # message(paste('length of S =',length(S)))
    S_to_E <- foi * dt
    E_to_P <- E * params[["epsilon"]] * dt
    P_to_A <- P * params[["fa"]] * zeta * dt
    P_to_I <- P * (1 - params[["fa"]]) * zeta * dt
    I_to_R <- I * params[["gamma"]]* dt
    A_to_R <- A * params[["gamma"]] * dt

    # Process model for SEPAIR Alpha
    S <- S - S_to_E
    E <- E + S_to_E - E_to_P
    P <- P + E_to_P - P_to_A - P_to_I
    A <- A + P_to_A - A_to_R
    I <- I + P_to_I - I_to_R
    R <- R + I_to_R + A_to_R

    daily_infected <- daily_infected + S_to_E
    daily_symptom_onset <- daily_symptom_onset + P_to_I
    daily_confirmed <- daily_confirmed + A_to_R + I_to_R

    # Delta variant
    if (t > delta_start) {
      foi_delta <- beta_delta * (params[["bp"]] * P_delta +
                                   params[["ba"]] * A_delta + I_delta)

      S_to_E_delta <- foi_delta * dt
      E_to_P_delta <- E_delta * params[["epsilon"]] * dt
      P_to_A_delta <- P_delta * params[["fa"]] * zeta * dt
      P_to_I_delta <- P_delta * (1 - params[["fa"]]) * zeta * dt
      I_to_R_delta <- I_delta * params[["gamma"]] * dt
      A_to_R_delta <- A_delta * params[["gamma"]] * dt

      S <- S - S_to_E_delta
      E_delta <- E_delta + S_to_E_delta - E_to_P_delta
      P_delta <- P_delta + E_to_P_delta - P_to_A_delta - P_to_I_delta
      A_delta <- A_delta + P_to_A_delta - A_to_R_delta
      I_delta <- I_delta + P_to_I_delta - I_to_R_delta
      R_delta <- R_delta + I_to_R_delta + A_to_R_delta

      daily_infected_delta <- daily_infected_delta + S_to_E_delta
      daily_symptom_onset_delta <- daily_symptom_onset_delta + P_to_I_delta
      daily_confirmed_delta <- daily_confirmed_delta + A_to_R_delta + I_to_R_delta
    }

    # Omicrion variant
    if (t > omicron_start) {
      foi_omicron <- beta_omicron * (params[["bp"]] * P_omicron +
                          params[["ba"]] * A_omicron + I_omicron)

      # message("t = ", t, ", beta_omicron = ", beta_omicron[1],
      #         ", I_omicron = ", I_omicron[1])

      S_to_E_omicron <- foi_omicron * dt
      E_to_P_omicron <- E_omicron * params[["epsilon"]] * dt
      P_to_A_omicron <- P_omicron * params[["fa"]] * zeta * dt
      P_to_I_omicron <- P_omicron * (1 - params[["fa"]]) * zeta * dt
      I_to_R_omicron <- I_omicron * params[["gamma"]] * dt
      A_to_R_omicron <- A_omicron * params[["gamma"]] * dt

      S <- S - S_to_E_omicron
      E_omicron <- E_omicron + S_to_E_omicron - E_to_P_omicron
      P_omicron <- P_omicron + E_to_P_omicron - P_to_A_omicron - P_to_I_omicron
      A_omicron <- A_omicron + P_to_A_omicron - A_to_R_omicron
      I_omicron <- I_omicron + P_to_I_omicron - I_to_R_omicron
      R_omicron <- R_omicron + I_to_R_omicron + A_to_R_omicron

      daily_infected_omicron <- daily_infected_omicron + S_to_E_omicron
      daily_symptom_onset_omicron <- daily_symptom_onset_omicron + P_to_I_omicron
      daily_confirmed_omicron <- daily_confirmed_omicron + A_to_R_omicron + I_to_R_omicron
    }
  }

  y1[, "S"] <- S
  y1[, "E"] <- E
  y1[, "P"] <- P
  y1[, "A"] <- A
  y1[, "I"] <- I
  y1[, "R"] <- R
  y1[, "CE"] <- daily_infected
  y1[, "CI"] <- daily_symptom_onset
  y1[, "CR"] <- daily_confirmed

  ylist <- list(Alpha = y1)

  if (t > delta_start) {
    y2[, "E"] <- E_delta
    y2[, "P"] <- P_delta
    y2[, "A"] <- A_delta
    y2[, "I"] <- I_delta
    y2[, "R"] <- R_delta
    y2[, "CE"] <- daily_infected_delta
    y2[, "CI"] <- daily_symptom_onset_delta
    y2[, "CR"] <- daily_confirmed_delta

    ylist[["Delta"]] <- y2
  }
  if (t > omicron_start) {
    y3[, "E"] <- E_omicron
    y3[, "P"] <- P_omicron
    y3[, "A"] <- A_omicron
    y3[, "I"] <- I_omicron
    y3[, "R"] <- R_omicron
    y3[, "CE"] <- daily_infected_omicron
    y3[, "CI"] <- daily_symptom_onset_omicron
    y3[, "CR"] <- daily_confirmed_omicron

    ylist[["Omicron"]] <- y3
  }

  return(ylist)
}
