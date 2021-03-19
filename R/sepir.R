#' Implements the rate of change of the SE1E2IR model for integration using deSolve
#'
#' The \code{sepir} returns the rate of change of the SE1E2IR model that can be integrated using deSolve package
#' @param t A vector of times for integration
#' @param y A vector of state variables
#' @param param A vector of parameter values (i.e., sigma, epsilon, gamma, R)
#' @export
#' @examples
#' res <- sepir(t, y, param)

sepir <- function (t, y, param) {
  sigma <- param$sigma
  epsilon <- param$epsilon
  gamma <- param$gamma
  beta <- param$R[floor(t)+1]*gamma
  presymp_infect <- param$presymp_infect
  if (presymp_infect) {
    beta <- param$R[floor(t)+1]*(1/(1/sigma - 1/epsilon + 1/gamma))
  }

  S <- y["S"]
  E <- y["E"]
  P <- y["P"]
  I <- y["I"]
  R <- y["R"]
  ## cumulative values
  CE <- y["CE"]
  CI <- y["CI"]
  CR <- y["CR"]

  N <- S + E + P + I + R

  foi <- beta * I / N
  if (presymp_infect) {
    foi <- beta*(P + I)/N #force of infection
  }
  zeta <- 1/(1/sigma - 1/epsilon) # rate from P to I

  dSdt = - S*foi
  dEdt = + S*foi - epsilon*E
  dPdt = + epsilon*E - zeta*P
  dIdt = + zeta*P - gamma*I
  dRdt = + gamma*I

  dCEdt = + S*foi
  dCIdt = + zeta*P
  dCRdt = + gamma *I # It is redundant, but I keep it to make it consitent across other files


  dydt <- c(dSdt, dEdt, dPdt, dIdt, dRdt, dCEdt, dCIdt, dCRdt)

  list(dydt)
}
