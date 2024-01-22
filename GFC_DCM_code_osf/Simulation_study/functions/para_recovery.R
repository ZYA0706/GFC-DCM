# compute the bias, MAE, MSE and RMSE of estimates
para_recovery <- function(est, real) {
  valid_indices <- is.finite(est)
  valid_indices <- valid_indices[valid_indices!=0]
  if (any(!valid_indices)) {
    warning("Estimates contain NA or infinite values. These values will be skipped in the calculations.")
  }
  est <- c(est);real <- c(real)
  bias <- mean(est[valid_indices] - real[valid_indices])
  MAE <- mean(abs(est[valid_indices] - real[valid_indices]) )
  MSE <- mean((est[valid_indices] - real[valid_indices])^2)
  RMSE <- sqrt(MSE)
  inds <- c(bias,MAE,MSE,RMSE)
  # inds <- data.frame(bias = bias, MAE = MAE, MSE = MSE, RMSE = RMSE,REL = REL)
  return(inds)
}
