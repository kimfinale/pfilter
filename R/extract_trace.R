#' Extract results from particle filtering and backward sampling (trace)
#'
#' The \code{extract_trace()} extracts trace
#' @param theta Parameter values
#' @param y lstate variables
#' @param data data to calculate likelihood against
#' @param data_type infection, symptom onset, or confirmation
#' @export
#' @examples
#' res <- extract_trace()
extract_trace <- function (params = NULL,
                       y = NULL,
                       beta0 = NULL,
                       data = NULL,
                       data_type = c("infection", "symptom onset", "confirmation"),
                       rep = 1,
                       npart = 100,
                       tend = 200,
                       dt = 0.2,
                       error_pdf = c("pois", "negbin", "binom"),
                       negbin_size = 5,
                       binom_prob = 0.8,
                       systematic_resampling = FALSE) {

  type <- match.arg(data_type)
  error_pdf <- match.arg(error_pdf)

  nstatevar <- length(y)
  if(is.data.frame(y)) {
    nstatevar <- ncol(y)
  }

  res <- particle_filter(params = params,
                         beta0 = beta0,
                         y = y, data = data, data_type = type,
                         npart = npart, tend = tend, dt = dt,
                         error_pdf = error_pdf,
                         negbin_size = negbin_size,
                         binom_prob = binom_prob,
                         systematic_resampling = systematic_resampling)

 # variables to return = latent state variable + beta
  output <- data.frame(matrix(NA, nrow = tend, ncol = nstatevar + 1))
  names(output) <- c(names(y), "Rt")

  output[, "S"] <- res$trace[["S"]]
  output[, "E"] <- res$trace[["E"]]
  output[, "P"] <- res$trace[["P"]]
  output[, "A"] <- res$trace[["A"]]
  output[, "I"] <- res$trace[["I"]]
  output[, "R"] <- res$trace[["R"]]
  output[, "CE"] <- res$trace[["CE"]]
  output[, "CI"] <- res$trace[["CI"]]
  output[, "CR"] <- res$trace[["CR"]]

  durP <- (1 / params[["delta"]] - 1 / params[["epsilon"]])
  durI <- (1 / params[["gamma"]])
  fa <- params[["frac_a"]]
  rhoa <- params[["rho_a"]]
  rhop <- params[["rho_p"]]
  R0_dur <- ((1 - fa) + fa * rhoa) * durI + rhop * durP
  # beta <- params["beta"]
  beta <- params[["R0"]] / R0_dur

  ##
  output[, "Rt"] <- res$trace[["beta"]] * R0_dur * res$trace[["S"]] / 50000030
  # output[, "Rt"] <- res$trace[["beta"]] * R0_dur
  return (output)
}
