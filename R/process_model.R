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
process_model <- function (params = NULL, y = NULL, tbegin = 0, tend = 1, dt = 0.2, beta = NULL) {
  S <- y[, "S"]
  E <- y[, "E"]
  P <- y[, "P"]
  I <- y[, "I"]
  R <- y[, "R"]

  daily_infected <- rep(0, nrow(y))
  daily_symptom_onset <- rep(0, nrow(y))
  daily_confirmed <- rep(0, nrow(y))

  # yini <- data.frame(S = S, E = E, P = P, I = I, R = R, CE = daily_infect, CI = daily_symptom, CR = daily_confirm)
  # out <- apply(yini, 1, run_step, times = c(tbegin, tend), params = as.list(params))
  # outdf <- data.frame(matrix(unlist(out), nrow = length(out), byrow = T))
  # outdf <- outdf[, -1] # remove the time column
  # names(outdf) <- names(y)

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

    daily_infected <- daily_infected + S_to_E
    daily_symptom_onset <- daily_symptom_onset + P_to_I
    daily_confirmed <- daily_confirmed + I_to_R

  }

  y[, "S"] <- S
  y[, "E"] <- E
  y[, "P"] <- P
  y[, "I"] <- I
  y[, "R"] <- R
  y[, "CE"] <- daily_infected
  y[, "CI"] <- daily_symptom_onset
  y[, "CR"] <- daily_confirmed


  # nm <- c("S", "E1", "E2", "I", "R", "CE1", "CE2", "CI")
  # for (s in nm) y[, s] <- eval(parse(text = s))

  return(y)
}
