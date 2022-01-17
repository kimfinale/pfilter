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
extract_trace_variant <- function (params = NULL,
                       y = NULL,
                       data = NULL,
                       data_type = c("infection", "symptom onset", "confirmation"),
                       npart = 10000,
                       tend = 200,
                       dt = 0.2,
                       error_pdf = c("pois", "negbin", "binom"),
                       negbin_size = 15,
                       binom_prob = 0.8) {

  type <- match.arg(data_type)
  error_pdf <- match.arg(error_pdf)

  res <- particle_filter_variant(params = params,
                         y = y, data = data, data_type = type,
                         npart = npart, tend = tend, dt = dt,
                         error_pdf = error_pdf,
                         negbin_size = negbin_size,
                         binom_prob = binom_prob)

 # variables to return = latent state variable + beta
  y1 <- y[["Alpha"]]

  nstatevar <- length(y1)
  if(is.data.frame(y1)) {
    nstatevar <- ncol(y1)
  }

  output <- data.frame(matrix(NA, nrow = tend, ncol = nstatevar + 1))
  names(output) <- c(names(y1), "Rt")

  for(nm in names(y1)) {
    output[, nm] <- res$trace[["Alpha"]][[nm]]
  }
  output[, "Rt"] <- res$trace[["Alpha"]][["beta"]] * R0_dur(params = params)
  output_list <- list(Alpha = output)
  # Delta variant
  if (length(res$trace) > 1) {
    y2 <- y[["Delta"]]
    nstatevar_delta <- length(y2)
    if(is.data.frame(y2)) {
      nstatevar_delta <- ncol(y2)
    }
    output_delta <- data.frame(matrix(NA, nrow = tend, ncol = nstatevar_delta + 1))
    names(output_delta) <- c(names(y2), "Rt")

    for(nm in names(y2)) {
      output_delta[, nm] <- res$trace[["Delta"]][[nm]]
    }
    output_delta[, "Rt"] <-
      res$trace[["Delta"]][["beta"]] * R0_dur(params = params)
    output_list$Delta <- output_delta
  }
  # Delta variant
  if (length(res$trace) > 2) {
    y3 <- y[["Omicron"]]
    nstatevar_omicron <- length(y3)
    if(is.data.frame(y3)) {
      nstatevar_omicron <- ncol(y3)
    }
    output_omicron <-
      data.frame(matrix(NA, nrow = tend, ncol = nstatevar_omicron + 1))
    names(output_omicron) <- c(names(y3), "Rt")

    for(nm in names(y3)) {
      output_omicron[, nm] <- res$trace[["Omicron"]][[nm]]
    }
    output_omicron[, "Rt"] <-
      res$trace[["Omicron"]][["beta"]] * R0_dur(params = params)
    output_list$Omicron <- output_omicron
  }
  return (output_list)
}
