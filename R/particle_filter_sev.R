#' This function implements forward filtering-backward sampling algorithm
#'
#' The \code{particle_filter_sv()} estimates posterior distribution of Rt and
#' severity parameter (i.e., prop of confirmed cases that go on to the severe
#' disease
#' population for both sexes and incidence rate
#' @param params Parameters
#' @param y A vector of state variables (e.g., c(S,E,P,I,R))
#' @param data data to calculate the likelihoods against
#' @param data_type data based on the date of infection, symptom onset, or confirmation
#' @param npart number of particles
#' @param dt time step size for the Euler method
#' @export
#' @import data.table
#' @examples
#' sample <- particle_filter()
particle_filter_sev <- function (params = theta,
                            y = y0,
                            beta0 = NULL,
                            data = NULL,
                            data_type = c("infection", "symptom onset", "confirmation"),
                            npart = 1000,
                            tend = 200,
                            dt = 0.2,
                            error_pdf = c("pois", "negbin", "binom"),
                            negbin_size = 5,
                            binom_prob = 0.8,
                            systematic_resampling = FALSE,
                            filter_traj = TRUE,
                            stoch = TRUE) {

  # Assumptions - using daily growth rate
  type <- match.arg(data_type)
  error_pdf <- match.arg(error_pdf)
  nstatevar <- length(y)
  if(is.data.frame(y)){
    nstatevar <- ncol(y)
  }
  latent_var <- array(0,
                  dim = c(npart, tend, nstatevar),
                  dimnames = list(NULL, NULL, names(y)))
  # Add initial values
  if(is.data.frame(y)) { ## output matrix from other particle filtering simulations
    for (i in 1:ncol(y)) {
      latent_var[, 1, i] <- y[, i]
    }
  }
  else if (is.vector(y)) {
    for (nm in names(y)) {
      latent_var[, 1, nm] <- y[[nm]]
    }
  }

  # compute the duration of infection that accounts for relative infectiousness
  # proportion of asymptomatic and symptomatic
  dur <- R0_dur(params)
  beta <- params[["R0"]] / dur
  beta_vol <-
    matrix(rnorm(npart*tend, mean=0, sd=params[["betavol"]]), nrow=tend)
  beta_vol[1,] <- beta * exp(beta_vol[1,]) # initial beta

  wt <- matrix(NA, nrow = npart, ncol = tend) # weight (likelihood)
  wt[, 1] <- 1 / npart  # initial weights
  W <- matrix(NA, nrow = npart, ncol = tend) # normalized weights
  A <- matrix(NA, nrow = npart, ncol = tend) # Resample according to the normalized weight

  # begin particle loop
  for (t in 2:tend) {
    # cat("t =", t, "\n")
    beta_vol[t, ] <- beta_vol[t-1, ] * exp(beta_vol[t, ])
    # run process model
    latent_var[, t, ] <- process_model_sev(params = params,
                                       y = latent_var[, t-1, ],
                                       tbegin = t-1,
                                       tend = t,
                                       dt = dt,
                                       beta = beta_vol[t,],
                                       fsev = fsev[t,],
                                       stoch = stoch)
    # calculate weights (likelihood)
    wt[, t] <- assign_weights(var = latent_var, t = t, data = data,
                              data_type = type, error_pdf = error_pdf,
                              negbin_size = negbin_size,
                              binom_prob = binom_prob)

    wt[, t] <- assign_weights(var = latent_var, t = t, data = data,
                              data_type = type, error_pdf = error_pdf,
                              negbin_size = negbin_size,
                              binom_prob = binom_prob)
    # normalize particle weights
    # cat("sum(is.na(wt[, t])) =", sum(is.na(wt[, t])), ", sum(wt[, t] =", sum(wt[, t]),"\n")
    W[, t] <- wt[, t] / sum(wt[, t])
    # resample particles by sampling parent particles according to weights
    A[, t] <- sample(1:npart, prob = W[1:npart, t], replace = T)
    if (systematic_resampling) {
      A[, t] <- systematic_resampling(W[1:npart, t])
    }
    # Resample particles for corresponding variables
    latent_var[, t,] <- latent_var[A[, t], t,]
    beta_vol[t,] <- beta_vol[t, A[, t]] # needed for random walk on beta
  } # end particle loop

  # Marginal likelihoods
  lik_values <- rep(NA, tend)
  for (t in 1:tend) {
    lik_values[t] <- log(sum(wt[1:npart, t])) # log-likelihoods
  }
  # averaged log likelihoods log(L/(npart^tend))
  loglik <- - tend * log(npart) + sum(lik_values)

  traj <- replicate(nstatevar + 1, rep(NA, tend), simplify = F)
  names(traj) <- c(names(y), "beta")
  if (filter_traj) {
    # Latent variables sampled from the filtered distribution (backward smoothing)
    loc <- rep(NA, tend)
    # loc_tend <-  sample(1:npart, prob = W[1:npart, tend], replace = T)
    loc[tend] <-  sample.int(npart, size = 1, prob = W[, tend], replace = T)
    # particle for the last time step
    traj[["beta"]][tend] <- beta_vol[tend, loc[tend]]
    for (nm in names(y)) {
      traj[[nm]][tend] <- latent_var[loc[tend], tend, nm]
    }
    # update backward
    for (i in seq(tend, 2, -1)) {
      loc[i-1] <- A[loc[i], i]
      traj[["beta"]][i-1] <- beta_vol[i-1, loc[i-1]]
      for (nm in names(y)) {
        traj[[nm]][i-1] <- latent_var[loc[i-1], i-1, nm]
      }
    }
  }


  # DEBUG  plot(Rep_traj[,1]-C_traj[,1])
  return (list(trace = traj, lik_marginal = lik_values,
               lik_overall_average = loglik, latent_var_filtered = latent_var,
              beta_filtered = beta_vol,
              W = W, A = A))
}
