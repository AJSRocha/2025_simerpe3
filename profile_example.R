# fix parameter on par list
# add map(list(logalpha = factor(NA)) )
# add map = map to the obj() call

library(TMB)

# Compile & load model
compile("model.cpp")
dyn.load(dynlib("model"))

# Simulate data
set.seed(123)
y <- rnorm(100, mean = 5, sd = 2)
data <- list(y = y)

# Define fixed grid of sigma values (log scale)
sigma_vals <- seq(log(0.5), log(4), length.out = 100)
profile_nll <- numeric(length(sigma_vals))

# Profile loop: fix sigma, estimate mu
for (i in seq_along(sigma_vals)) {
  fixed_sigma <- sigma_vals[i]
  
  obj <- MakeADFun(
    data = data,
    parameters = list(mu = 0, sigma = fixed_sigma),
    map = list(sigma = factor(NA)),  # fix sigma
    DLL = "model"
  )
  
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  profile_nll[i] <- opt$objective
}

# Compute threshold for 95% CI
nll_min <- min(profile_nll)
threshold <- nll_min + qchisq(0.95, df = 1) / 2  # ≈ nll_min + 1.92

# Interpolate to find bounds
above <- profile_nll > threshold
ci_indices <- which(!above)
lower_idx <- min(ci_indices)
upper_idx <- max(ci_indices)

# Add a simple linear interpolation for better precision
approx_left <- approx(profile_nll[1:lower_idx], exp(sigma_vals[1:lower_idx]), xout = threshold)$y
approx_right <- approx(profile_nll[upper_idx:length(profile_nll)],
                       exp(sigma_vals[upper_idx:length(profile_nll)]),
                       xout = threshold)$y

cat("95% profile likelihood CI for sigma:\n")
cat(sprintf("Lower bound: %.3f\n", approx_left))
cat(sprintf("Upper bound: %.3f\n", approx_right))

# Plot
plot(exp(sigma_vals), profile_nll, type = "l", lwd = 2,
     xlab = "sigma", ylab = "Negative Log-Likelihood",
     main = "Profile Likelihood for sigma")
abline(h = threshold, col = "red", lty = 2)
abline(v = approx_left, col = "blue", lty = 3)
abline(v = approx_right, col = "blue", lty = 3)
legend("topright", legend = c("Threshold", "95% CI bounds"),
       col = c("red", "blue"), lty = c(2, 3))

# --- 1. Estimate both parameters (for Hessian CI) ---
obj_full <- MakeADFun(data = data, parameters = parameters, DLL = "model")
opt_full <- nlminb(obj_full$par, obj_full$fn, obj_full$gr)
rep <- sdreport(obj_full)

# Extract MLE and SE of log(sigma)
log_sigma_hat <- rep$value["sigma"]
log_sigma_se  <- rep$sd["sigma"]

# Wald (Hessian) 95% CI on log scale
wald_log_ci <- log_sigma_hat + c(-1.96, 1.96) * log_sigma_se
wald_ci <- exp(wald_log_ci)

cat("Wald (Hessian-based) 95% CI for sigma:\n")
cat(sprintf("Lower: %.3f, Upper: %.3f\n", wald_ci[1], wald_ci[2]))

# --- 2. Profile likelihood for sigma ---
sigma_vals <- seq(log(0.5), log(4), length.out = 100)
profile_nll <- numeric(length(sigma_vals))

for (i in seq_along(sigma_vals)) {
  fixed_sigma <- sigma_vals[i]
  
  obj <- MakeADFun(
    data = data,
    parameters = list(mu = 0, sigma = fixed_sigma),
    map = list(sigma = factor(NA)),
    DLL = "model"
  )
  
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  profile_nll[i] <- opt$objective
}

# Compute threshold for 95% profile CI
nll_min <- min(profile_nll)
threshold <- nll_min + qchisq(0.95, df = 1) / 2  # ≈ +1.92

# Interpolate profile CI bounds
ci_indices <- which(profile_nll <= threshold)
lower_idx <- min(ci_indices)
upper_idx <- max(ci_indices)

approx_left <- approx(profile_nll[1:lower_idx], exp(sigma_vals[1:lower_idx]), xout = threshold)$y
approx_right <- approx(profile_nll[upper_idx:length(profile_nll)],
                       exp(sigma_vals[upper_idx:length(profile_nll)]),
                       xout = threshold)$y

cat("Profile likelihood 95% CI for sigma:\n")
cat(sprintf("Lower: %.3f, Upper: %.3f\n", approx_left, approx_right))

# --- 3. Plot both CIs on profile likelihood curve ---
plot(exp(sigma_vals), profile_nll, type = "l", lwd = 2,
     xlab = "sigma (original scale)", ylab = "Negative Log-Likelihood",
     main = "Profile Likelihood with CI Comparison")

# Add profile CI
abline(h = threshold, col = "red", lty = 2)
abline(v = approx_left, col = "blue", lty = 3)
abline(v = approx_right, col = "blue", lty = 3)

# Add Hessian CI
abline(v = wald_ci[1], col = "darkgreen", lty = 4)
abline(v = wald_ci[2], col = "darkgreen", lty = 4)

# Legend
legend("topright", legend = c("Profile Threshold", "Profile CI", "Wald (Hessian) CI"),
       col = c("red", "blue", "darkgreen"), lty = c(2, 3, 4), lwd = 2)


#TODO

opt <- optim(
  par = obj$par,
  fn = obj$fn,
  gr = obj$gr,
  method = "CG",
  control = list(maxit = 1000)
)