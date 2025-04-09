#===============================================================================
# Function to calculate the effective sample size in a multilevel model
#===============================================================================

# Note that this function requires the packages "lme4" and "Matrix"

get_neff <- function(model) {
  # Extract random effects structure and residual variance
  reVar <- VarCorr(model)
  sigma2 <- sigma(model)^2
  
  # Extract the grouping structure and number of observations per subject
  N <- length(unique(model@frame$id))  # Number of subjects
  n <- nrow(model@frame) / N  # Number of observations per subject
  
  # Extract subject-specific Z matrix
  Zt <- as.matrix(getME(model, "Zt")[1:2, 1:n])  # Transposed random-effects design matrix (Z')
  Z <- t(Zt)
  
  X_control <- model.matrix(model)[1:n,]
  X_treat <- model.matrix(model)[(n+1):(2*n),]
  
  # Construct D for one subject
  D <- as.matrix(bdiag(lapply(reVar, function(x) as(bdiag(x), "sparseMatrix"))))
  
  # Compute V for the smaller Z and X
  V <- Z %*% D %*% t(Z) + sigma2 * Diagonal(nrow(Z))  # V = Z D Z' + σ² I
  V_inv <- solve(V)
  
  # Compute W (diagonal matrix of V)
  W <- Diagonal(x = diag(V))
  W_inv <- solve(W)
  
  # Compute sum_mat and sum_mat_indep
  sum_mat_control <- t(X_control) %*% V_inv %*% X_control
  sum_mat_treat <- t(X_treat) %*% V_inv %*% X_treat
  sum_mat <- sum_mat_control + sum_mat_treat
  
  sum_mat_control_indep <- t(X_control) %*% W_inv %*% X_control
  sum_mat_treat_indep <- t(X_treat) %*% W_inv %*% X_treat
  sum_mat_indep <- sum_mat_control_indep + sum_mat_treat_indep
  
  # Compute var_beta_hat and var_betahat_indep
  var_beta_hat <- solve(sum_mat)
  var_betahat_indep <- solve(sum_mat_indep)
  
  # Compute weight w and effective sample size N_eff
  w <- var_betahat_indep / var_beta_hat
  N_eff <- w * N * n
  
  return(N_eff[4, 4])
}
