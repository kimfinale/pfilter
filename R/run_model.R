#' Run particle filtering models and sample trajectories of variables
#'
#' The \code{run_model()} runs particle filtering models
#' @param theta Parameter values
#' @param y lstate variables
#' @param data data to calculate likelihood against
#' @param data_type infection, symptom onset, or confirmation
#' @export
#' @examples
#' res <- run_model()
# Likelihood calc for SMC --------------------------------------------
run_model <- function (params = theta,
                       y = y0,
                       beta0 = NULL,
                       data = Rt_data,
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

  # output <- list()
  output <- replicate((length(y)+1), matrix(NA, nrow = tend, ncol = rep), simplify = FALSE)
  names(output) <-  c(names(y), "Rt")

  for (i in seq_len(rep)) {
    cat( "i =", i, "\n" )
    res <- particle_filter(params = params,
                           beta0 = beta0,
                           y = y, data = data, data_type = type,
                           npart = npart, tend = tend, dt = dt,
                           error_pdf = error_pdf,
                           negbin_size = negbin_size,
                           binom_prob = binom_prob,
                           systematic_resampling = systematic_resampling)
    # invisible(lapply(names(y), function(x) output[[x]][,i] <- res$trace[[x]]))
    output[["S"]][, i] <- res$trace[["S"]]
    output[["E"]][, i] <- res$trace[["E"]]
    output[["P"]][, i] <- res$trace[["P"]]
    output[["I"]][, i] <- res$trace[["I"]]
    output[["R"]][, i] <- res$trace[["R"]]
    output[["CE"]][, i] <- res$trace[["CE"]]
    output[["CI"]][, i] <- res$trace[["CI"]]
    output[["CR"]][, i] <- res$trace[["CR"]]
    output[["Rt"]][, i] <- res$trace[["beta"]]/params[["gamma"]]
    # output$latent_var_filtered[[i]] <- res$latent_var_filtered
    # output$beta_filtered[[i]] <- res$beta_vol
    # output$W[[i]] <- res$W
    # output$A[[i]] <- res$A
  }
  return (output)
  # saveRDS(output, file = paste0("outputs/fit_", filename, ".rds"))
}
