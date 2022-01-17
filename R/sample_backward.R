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
sample_backward <- function (var1 = NULL,
                             var2 = NULL,
                             weight = NULL,
                             sample = NULL,
                             npart = 1000,
                             tend = 200) {

# Latent variables sampled from the filtered distribution (backward sampling)
loc <- rep(NA, tend)
# loc_tend <-  sample(1:npart, prob = W[1:npart, tend], replace = T)
loc[tend] <-  sample.int(npart, size = 1, prob = weight[, tend], replace = T)
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



  return (traj)
}
