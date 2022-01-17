#' Assign weights (i.e., likelihoods) to each of the value
#'
#' The \code{assign_weights()} calculate the likelihoods assuming case incidences are Poisson-distributed random numbers
#' @param var A vector of state variables from simulation
#' @param t current time
#' @param data data on the times of the infection transmission process (e.g., infection, symptom onset, or confirmation)
#' @param data_type c("infection", "symptom onset", "confirmation")
#' @export
#' @examples
#'  wt <- assign_weights_variant(var=latent_var, t=t, data=data, data_type=type)
assign_weights_variant <-
  function (params,
            var,
            t,
            data,
            data_type = c("infection", "symptom onset", "confirmation"),
            error_pdf = c("pois", "negbin", "binom"),
            negbin_size = 5,
            binom_prob = 0.8) {

  type <- match.arg(data_type)
  error_pdf <- match.arg(error_pdf)

  delta_start <- params[["Delta_start"]]
  omicron_start <- params[["Omicron_start"]]

  var_alpha <- var[["Alpha"]]
  dat_alpha <- data[["Alpha"]]

  if (type == "infection") {
    case_expected <- var_alpha[, "CE"]
    case_data <- round(unlist(dat_alpha[t, "daily_infected"]))
  }
  else if (type == "symptom onset") {
    case_expected <- var_alpha[, "CI"]
    case_data <- round(unlist(dat_alpha[t, "daily_symptom_onset"]))
  }
  else if (type == "confirmation") {
    case_expected <- var_alpha[, "CR"]
    case_data <- round(unlist(dat_alpha[t, "daily_confirmed"]))
  }

  if (!is.na(case_data)) {
    expected_val <- pmax(0, case_expected)# case_expected is a vector of length npart
    if (error_pdf == "pois"){
      log_lik <- dpois(round(case_data), lambda = expected_val, log = T)
    }
    else if (error_pdf == "negbin") {
      log_lik <- dnbinom(round(case_data),
                         size = negbin_size, mu = expected_val, log = T)
    }
    else if (error_pdf == "binom") {
      log_lik <- dbinom(round(case_data),
                        size = round(expected_val), prob = binom_prob, log = T)
    }
  }
  else {
    log_lik <- -Inf
  }

  # repeat for the Delta variant
  if (t > delta_start) {
    var_delta <- var[["Delta"]] # returns dataframe unlike latent_var, 3D array
    dat_delta <- data[["Delta"]]
    if (type == "infection") {
      case_expected_delta <- var_delta[, "CE"]
      case_data_delta <- round(unlist(dat_delta[t, "daily_infected"]))
    }
    else if (type == "symptom onset") {
      case_expected_delta <- var_delta[, "CI"]
      case_data_delta <- round(unlist(dat_delta[t, "daily_symptom_onset"]))
    }
    else if (type == "confirmation") {
      case_expected_delta <- var_delta[, "CR"]
      case_data_delta <- round(unlist(dat_delta[t, "daily_confirmed"]))
    }
    if (!is.na(case_data_delta)) {# case_expected is a vector of length npart
      expected_val_delta <- pmax(0, case_expected_delta)
      if (error_pdf == "pois"){
        log_lik_delta <- dpois(round(case_data_delta),
                               lambda = expected_val_delta, log = T)
      }
      else if (error_pdf == "negbin") {
        log_lik_delta <- dnbinom(round(case_data_delta),
                  size = negbin_size, mu = expected_val_delta, log = T)
      }
      else if (error_pdf == "binom") {
        log_lik_delta <- dbinom(round(case_data_delta),
                  size = round(expected_val_delta), prob = binom_prob, log = T)
      }
    }
    else {
      log_lik_delta <- -Inf
    }
    log_lik <- log_lik + log_lik_delta # update the likelihood
  }

  if (t > omicron_start) {
    var_omicron <- var[["Omicron"]] # returns dataframe unlike latent_var, 3D array
    dat_omicron <- data[["Omicron"]]

    if (type == "infection") {
      case_expected_omicron <- var_omicron[, "CE"]
      case_data_omicron <- round(unlist(dat_omicron[t, "daily_infected"]))
    }
    else if (type == "symptom onset") {
      case_expected_omicron <- var_omicron[, "CI"]
      case_data_omicron <- round(unlist(dat_omicron[t, "daily_symptom_onset"]))
    }
    else if (type == "confirmation") {
      case_expected_omicron <- var_omicron[, "CR"]
      case_data_omicron <- round(unlist(dat_omicron[t, "daily_confirmed"]))
    }
#
    # message("t = ", t, ", dat = ", case_data_omicron,  ", model = ",
    #         case_expected_omicron[1:5])
    # repeat for the Delta variant
    if (!is.na(case_data_omicron)) {# case_expected is a vector of length npart
      expected_val_omicron <- pmax(0, case_expected_omicron)
      if (error_pdf == "pois"){
        log_lik_omicron <- dpois(round(case_data_omicron),
                               lambda = expected_val_omicron, log = T)
      }
      else if (error_pdf == "negbin") {
        log_lik_omicron <- dnbinom(round(case_data_omicron),
                                 size = negbin_size,
                                 mu = expected_val_omicron, log = T)
      }
      else if (error_pdf == "binom") {
        log_lik_omicron <- dbinom(round(case_data_omicron),
                                size = round(expected_val_omicron),
                                prob = binom_prob, log = T)
      }
    }
    else {
      log_lik_omicron <- -Inf
    }

    log_lik <- log_lik + log_lik_omicron # update the likelihood
  }

  # cat("sum(is.na(log_lik)) =", sum(is.na(log_lik)), ", sum(log_lik) =", sum(log_lik),"\n")
  return (exp(log_lik)) # convert to normal probability
}

