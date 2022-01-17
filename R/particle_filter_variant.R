#' This function implements forward filtering-backward sampling algorithm
#'
#' The \code{particle_filter_variant()} estimates posterior distribution of Rt of
#' variants
#' @param params Parameters
#' @param y A vector of state variables (e.g., c(S,E,P,I,R))
#' @param data data to calculate the likelihoods against
#' @param data_type data based on the date of infection, symptom onset, or confirmation
#' @param npart number of particles
#' @param dt time step size for the Euler method
#' @export
#' @import data.table
#' @examples
#' sample <- particle_filter_variant()
particle_filter_variant <- function (params = theta,
                            y = y0,
                            data = NULL,
                            data_type = c("infection", "symptom onset", "confirmation"),
                            npart = 1000,
                            tend = 200,
                            dt = 0.2,
                            error_pdf = c("pois", "negbin", "binom"),
                            negbin_size = 5,
                            binom_prob = 0.8,
                            filter_traj = TRUE) {

  y1 <- y[["Alpha"]]

  type <- match.arg(data_type)
  error_pdf <- match.arg(error_pdf)

  nstatevar <- length(y1)
  if(is.data.frame(y1)){
    nstatevar <- ncol(y1)
  }

  latent_var <- array(0,
                  dim = c(npart, tend, nstatevar),
                  dimnames = list(NULL, NULL, names(y1)))

  # Add initial values
  if(is.data.frame(y1)) { ## output matrix from other particle filtering simulations
    for (i in 1:ncol(y1)) {
      latent_var[, 1, i] <- y1[, i]
    }
  }
  else if (is.atomic(y1)) {
    for (nm in names(y1)) {
      latent_var[, 1, nm] <- y1[[nm]]
    }
  }

  delta_start <- params[["Delta_start"]]
  omicron_start <- params[["Omicron_start"]]

  # Delta
  y2 <- y[["Delta"]] # no S compartment
  nstatevar_delta <- length(y2)
  if(is.data.frame(y2)){
    nstatevar_delta  <- ncol(y2)
  }

  latent_var_delta <- array(0,
                            dim = c(npart, tend, nstatevar_delta), # no S compartment
                            dimnames = list(NULL, NULL, names(y2)))
  if(is.data.frame(y2)) { ## output matrix from other particle filtering simulations
    for (i in 1:ncol(y2)) {
      latent_var[, 1, i] <- y2[, i]
    }
  }
  else if (is.atomic(y2)) {
    for (nm in names(y2)) {
      latent_var_delta[, 1, nm] <- y2[[nm]]
    }
  }

  # Omicron
  y3 <- y[["Omicron"]] # no S compartment
  nstatevar_omicron <- length(y3)
  if(is.data.frame(y3)){
    nstatevar_omicron  <- ncol(y3)
  }

  latent_var_omicron <- array(0,
                            dim = c(npart, tend, nstatevar_omicron), # no S compartment
                            dimnames = list(NULL, NULL, names(y3)))
  if(is.data.frame(y3)) { ## output matrix from other particle filtering simulations
    for (i in 1:ncol(y3)) {
      latent_var_omicron[, 1, i] <- y3[, i]
    }
  }
  else if (is.atomic(y3)) {
    for (nm in names(y3)) {
      latent_var_omicron[, 1, nm] <- y3[[nm]]
    }
  }

  # to check if the error caused by the latent_var was set to zero
  # when the Omicron variable started
  for (nm in names(y3)) {
    for (tt in 2:omicron_start+1){
      latent_var_omicron[, tt, nm] <- y3[[nm]]
    }
  }

  # compute the duration of infection that accounts for relative infectiousness
  # proportion of asymptomatic and symptomatic
  dur <- R0_dur(params)
  beta0 <- params[["R0"]] / dur
  beta_vol <-
    matrix(rnorm(npart*tend, mean=0, sd=params[["betavol"]]), nrow=tend)
  beta_vol[1,] <- beta0 * exp(beta_vol[1,]) # initial beta

  # Delta
  beta_vol_delta <-
    matrix(rnorm(npart*tend, mean=0, sd=params[["betavol"]]), nrow=tend)
  beta_vol_delta[1,] <- beta0 * exp(beta_vol_delta[1,]) # initial beta
  # Omicron
  beta_vol_omicron <-
    matrix(rnorm(npart*tend, mean=0, sd=params[["betavol"]]), nrow=tend)
  beta_vol_omicron[1,] <- beta0 * exp(beta_vol_omicron[1,]) # initial beta

  # update the beta_vol_omicron to be initial values when the Omicron starts
  for (tt in 2:omicron_start){
    beta_vol_omicron[tt, ] <- beta_vol_omicron[1,] # initial beta
  }

  wt <- matrix(NA, nrow = npart, ncol = tend) # weight (likelihood)
  wt[, 1] <- 1 / npart  # initial weights
  W <- matrix(NA, nrow = npart, ncol = tend) # normalized weights
  A <- matrix(NA, nrow = npart, ncol = tend) # Resample according to the normalized weight

  # begin particle loop
  for (t in 2:tend) {
    # cat("t =", t, "\n")
    beta_vol[t,] <- beta_vol[t-1,] * exp(beta_vol[t,])
    varlist <- list(Alpha = latent_var[, t-1,])
    betalist <- list(Alpha = beta_vol[t,])
    if (t > delta_start) {
      beta_vol_delta[t,] <- beta_vol_delta[t-1,] * exp(beta_vol_delta[t,])
      varlist[["Delta"]] <- latent_var_delta[, t-1,]
      betalist[["Delta"]] <- beta_vol_delta[t,]
    }
    if (t > omicron_start) {
      beta_vol_omicron[t,] <- beta_vol_omicron[t-1,] * exp(beta_vol_omicron[t,])
      varlist[["Omicron"]] <- latent_var_omicron[, t-1,]
      betalist[["Omicron"]]<- beta_vol_omicron[t,]
    }
    # run process model
    varlist <- process_model_variant(t = t,
                                     params = params,
                                     y = varlist,
                                     tbegin = t-1,
                                     tend = t,
                                     dt = dt,
                                     beta = betalist)

    # calculate weights (likelihood)
    wt[, t] <- assign_weights_variant(params = params,
                                      var = varlist, t = t, data = data,
                                      data_type = type, error_pdf = error_pdf,
                                      negbin_size = negbin_size,
                                      binom_prob = binom_prob)
    # normalize particle weights
    # cat("sum(is.na(wt[, t])) =", sum(is.na(wt[, t])), ", sum(wt[, t] =", sum(wt[, t]),"\n")
    W[, t] <- wt[, t] / sum(wt[, t])
    # resample particles by sampling parent particles according to weights
    A[, t] <- sample(1:npart, prob = W[1:npart, t], replace = T)

    latent_var[, t,] <- varlist[["Alpha"]]
    # update variables
    latent_var[, t,] <- latent_var[A[, t], t,]
    beta_vol[t,] <- beta_vol[t, A[, t]] # needed for random walk on beta

    if (t > delta_start) {
      latent_var_delta[, t,] <- varlist[["Delta"]]
      latent_var_delta[, t,] <- latent_var_delta[A[, t], t,]
      beta_vol_delta[t,] <- beta_vol_delta[t, A[, t]] # needed for random walk on beta
    }
    if (t > omicron_start) {
      latent_var_omicron[, t,] <- varlist[["Omicron"]]
      latent_var_omicron[, t,] <- latent_var_omicron[A[, t], t,]
      beta_vol_omicron[t,] <- beta_vol_omicron[t, A[, t]] # needed for random walk on beta
    }
  } # end particle loop

  # Marginal likelihoods
  lik_values <- rep(NA, tend)
  for (t in 1:tend) {
    lik_values[t] <- log(sum(wt[1:npart, t])) # log-likelihoods
  }
  # averaged log likelihoods log(L/(npart^tend))
  loglik <- - tend * log(npart) + sum(lik_values)

  traj <- replicate(nstatevar + 1, rep(NA, tend), simplify = F)
  names(traj) <- c(names(y1), "beta")
  if (t > delta_start) {
    traj_delta <- replicate(nstatevar_delta + 1, rep(NA, tend), simplify = F) # no S variable
    names(traj_delta) <- c(names(y2), "beta")
  }
  if (t > omicron_start) {
    traj_omicron <- replicate(nstatevar_omicron + 1, rep(NA, tend), simplify = F) # no S variable
    names(traj_omicron) <- c(names(y3), "beta")
  }

  if (filter_traj) {
    # Latent variables sampled from the filtered distribution (backward sampling)
    loc <- rep(NA, tend)
    loc[tend] <-  sample.int(npart, size = 1, prob = W[, tend], replace = T)
    # particle for the last time step
    traj[["beta"]][tend] <- beta_vol[tend, loc[tend]]

    for (nm in names(y1)) {
      traj[[nm]][tend] <- latent_var[loc[tend], tend, nm]
    }
    if (t > delta_start) {
      traj_delta[["beta"]][tend] <- beta_vol_delta[tend, loc[tend]]
      for (nm in names(y2)) {
        traj_delta[[nm]][tend] <- latent_var_delta[loc[tend], tend, nm]
      }
    }
    if (t > omicron_start) {
      traj_omicron[["beta"]][tend] <- beta_vol_omicron[tend, loc[tend]]
      for (nm in names(y3)) {
        traj_omicron[[nm]][tend] <- latent_var_omicron[loc[tend], tend, nm]
      }
    }
    # update backward
    for (i in seq(tend, 2, -1)) {
      loc[i-1] <- A[loc[i], i]
      traj[["beta"]][i-1] <- beta_vol[i-1, loc[i-1]]
      for (nm in names(y1)) {
        traj[[nm]][i-1] <- latent_var[loc[i-1], i-1, nm]
      }
      if (t > delta_start) {
        traj_delta[["beta"]][i-1] <- beta_vol_delta[i-1, loc[i-1]]
        for (nm in names(y2)) {
          traj_delta[[nm]][i-1] <- latent_var_delta[loc[i-1], i-1, nm]
        }
      }
      if (t > omicron_start) {
        traj_omicron[["beta"]][i-1] <- beta_vol_omicron[i-1, loc[i-1]]
        for (nm in names(y3)) {
          traj_omicron[[nm]][i-1] <- latent_var_omicron[loc[i-1], i-1, nm]
        }
      }
    }
  }

  trajlist <- list(Alpha = traj)
  if (t > delta_start) {
    trajlist[["Delta"]] = traj_delta
  }
  if (t > omicron_start) {
    trajlist[["Omicron"]] = traj_omicron
  }

  return (list(trace = trajlist, lik_marginal = lik_values,
               lik_overall_average = loglik,
               latent_var_filtered = latent_var,
               beta_filtered = beta_vol,
               W = W, A = A))
}
