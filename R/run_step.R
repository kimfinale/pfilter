#' Runs a step (day in this case)
#'
#' The \code{step_run} returns the rate of change of the SE1E2IR model that can be integrated using deSolve package
#' @param y A vector of state variables
#' @param times Times over which ODE is integrated (1 day in our case)
#' @param params A vector of parameter values (i.e., sigma, epsilon, gamma, R)
#' @import deSolve
#' @export
#' @examples
#' res <- sepir(t, y, param)
run_step <- function(y, times, params){
  params$presymp_infect <- FALSE
  out <- as.data.frame(ode(y = y, func = sepir, times = times, parms = params))
  tail(out, 1)
}
