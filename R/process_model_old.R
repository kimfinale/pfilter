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
process_model_old <- function (params = NULL, y = NULL, tbegin = 0, tend = 1, dt = 0.2, beta = NULL) {
  S <- y[, "S"]
  E <- y[, "E"]
  P <- y[, "P"]
  I <- y[, "I"]
  R <- y[, "R"]
  daily_infect <- rep(0, nrow(y))
  daily_symptom <- rep(0, nrow(y))
  daily_confirm <- rep(0, nrow(y))

  pop <- S + E + P + I + R
  zeta <- 1/(1/params[["sigma"]] - 1/params[["epsilon"]]) # rate from P to I

  for (ii in seq((tbegin + dt), tend, dt)) {
    foi <- beta * I / pop
    S_to_E <- S * foi * dt
    E_to_P <- E * params[["epsilon"]] * dt
    P_to_I <- P * zeta * dt
    I_to_R <- I * params[["gamma"]] * dt

    # Process model for SEPIR
    S <- S - S_to_E
    E <- E + S_to_E - E_to_P
    P <- P + E_to_P - P_to_I
    I <- I + P_to_I - I_to_R
    R <- R + I_to_R

    daily_infect <- daily_infect + S_to_E
    daily_symptom <- daily_symptom + P_to_I
    daily_confirm <- daily_confirm + I_to_R

  }

  y[, "S"] <- S
  y[, "E"] <- E
  y[, "P"] <- P
  y[, "I"] <- I
  y[, "R"] <- R
  y[, "dE"] <- daily_infect
  y[, "dI"] <- daily_symptom
  y[, "dR"] <- daily_confirm


  # nm <- c("S", "E1", "E2", "I", "R", "CE1", "CE2", "CI")
  # for (s in nm) y[, s] <- eval(parse(text = s))

  return(y)
}
