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
#' @export
#' @import data.table
#' @examples
#' pop <- process_model(params, y, tbegin = 0, tend = 1,)
# Process model for simulation --------------------------------------------
process_model_sev <- function (params = NULL,
                           y = NULL,
                           tbegin = 0,
                           tend = 1,
                           dt = 0.2,
                           beta = NULL,
                           fsev = NULL, # prop of becomeing severe cases
                           stoch = TRUE) {

  # daily_* values are updated every process_model_step
  y[, c("CE", "CI", "CR", "CS")] <- 0

  S <- y[, "S"]
  E <- y[, "E"]
  P <- y[, "P"]
  A <- y[, "A"]
  I <- y[, "I"]
  R <- y[, "R"]
  daily_infected <- y[, "CE"]
  daily_symptom_onset <- y[, "CI"]
  daily_confirmed <- y[, "CR"]
  daily_severe <- y[, "CS"]

  pop <- S + E + P + A + I + R
  durP <- (1 / params[["delta"]] - 1 / params[["epsilon"]]) # residence time in P
  zeta <- 1 / durP # rate from P to I

  for (ii in seq((tbegin + dt), tend, dt)) {
    foi <- beta * (params[["bp"]] * P + params[["ba"]] * A + I)
    # message(paste('length of S =',length(S)))
    if (stoch){
      S_to_E <- rbinom(length(S), S, 1 - exp( - foi / S * dt))
      E_to_P <- rbinom(length(E), E, 1 - exp( - params[["epsilon"]] * dt))
      P_to_AI <- rbinom(length(P), P, 1 - exp( - zeta * dt))
      P_to_A <- rbinom(length(P_to_AI), P_to_AI, params[["fa"]])
      P_to_I <- P_to_AI - P_to_A
      I_to_R <- rbinom(length(I), I, 1 - exp( - params[["gamma"]] * dt))
      A_to_R <- rbinom(length(A), A, 1 - exp( - params[["gamma"]] * dt))
    } else {
      S_to_E <- foi * dt
      E_to_P <- E * params[["epsilon"]] * dt
      P_to_A <- P * params[["fa"]] * zeta * dt
      P_to_I <- P * (1 - params[["fa"]]) * zeta * dt
      I_to_R <- I * params[["gamma"]]* dt
      A_to_R <- A * params[["gamma"]] * dt
      # 1/eta = mean days before symptoms become severe after confirmation
      # 1/gamma + 1/eta ~ 7 days
      # fsev = fraction of cases that become severe
      # R represents confirmed cases, not removed ones in this case
      R_to_Sev <- R * params[["eta"]] * fsev * dt
    }

    # Process model for SEPIR
    S <- S - S_to_E
    E <- E + S_to_E - E_to_P
    P <- P + E_to_P - P_to_A - P_to_I
    A <- A + P_to_A - A_to_R
    I <- I + P_to_I - I_to_R
    R <- R + I_to_R + A_to_R

    daily_infected <- daily_infected + S_to_E
    daily_symptom_onset <- daily_symptom_onset + P_to_I
    daily_confirmed <- daily_confirmed + A_to_R + I_to_R
    daily_severe <- daily_severe + R_to_Sev

  }

  y[, "S"] <- S
  y[, "E"] <- E
  y[, "P"] <- P
  y[, "A"] <- A
  y[, "I"] <- I
  y[, "R"] <- R
  y[, "CE"] <- daily_infected
  y[, "CI"] <- daily_symptom_onset
  y[, "CR"] <- daily_confirmed
  y[, "CS"] <- daily_severe

  return(y)
}
