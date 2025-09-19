
optimx_out = spg1
use_units = F

estimate_catdyn_biomass <- function(par, effort, catch = cat_df$Data$`MIS+OTB-S`$obscat.kg, nstep = nT, units = list(), 
                                    use_units = FALSE, model_type = 0, distribution = "apnormal") {
  
  # EXPONENTIATE parameters if needed

  # Extract parameters (CatDyn ordering)
  M     <- par[["logM"]] %>% exp
  N0    <- par[["logN0_scaled"]] %>% exp
  k     <- par[["logK"]] %>% exp
  alpha <- par[["logalpha"]] %>% exp
  beta  <- par[["logbeta"]] %>% exp
  
  # Biomass vector
  B <- numeric(nT)
  B[1] <- N0 * nmult
  
  # Optional: Unit scaling
  if (use_units && length(units) >= 4) {
    unit_weight <- switch(units[[3]],
                          "g"  = 1e-6,
                          "kg" = 0.001,
                          1)
    unit_scale <- switch(units[[4]],
                         "bill" = 1e9,
                         "mill" = 1e6,
                         "thou" = 1e3,
                         "hund" = 1e2,
                         1)
  } else {
    unit_weight <- 1
    unit_scale  <- 1000
  }
  
  for (t in 2:nstep) {
    F <- NA
    # Update biomass with a delay-difference-like model
    Bt_1 <- B[t - 1]
    Ct_1 <- nmult * catch[t - 1]
    B[t] <- Bt_1 + k * (alpha * Bt_1^beta * exp(-M) - Ct_1)  # Simplified DD model
  }
  
  # Exploitation rate = Catch / Biomass
  U_obs <- catch / B
  
  # Optional: Estimate instantaneous fishing mortality F from U
  estimate_F <- function(U, M) {
    if (U <= 0) return(0)
    uniroot(function(F) U - (F / (F + M)) * (1 - exp(-F - M)),
            c(1e-8, 10), extendInt = "yes")$root
  }
  F_est <- sapply(U_obs, estimate_F, M = M)
  
  # Return results
  results <- data.frame(
    time = seq_len(nstep),
    biomass = B * unit_weight * unit_scale,
    catch = catch,
    effort = effort,
    exploitation_rate = U_obs,
    F = F_est
  )
  
  return(results)
}

library(msm)  # for deltamethod

estimate_biomass_uncertainty <- function(par, cov_mat, timing, weight_kg_mean, weight_kg_sd) {
  # Assume: par is on log scale, so we exponentiate
  par_exp <- exp(par)
  names(par_exp) <- names(par)
  
  M  <- par_exp["M"]
  N0 <- par_exp["N0"]
  P1 <- par_exp["P1"]
  
  # Create storage vectors
  n <- length(timing)
  N_thou       <- numeric(n)
  N_thou_se    <- numeric(n)
  B_ton        <- numeric(n)
  B_ton_se     <- numeric(n)
  
  for (i in seq_along(timing)) {
    # Predict abundance in thousands
    Pt <- timing[i] * P1
    Nt <- N0 * exp(-M) + Pt * exp(-M)
    N_thou[i] <- Nt
    
    # Delta method SE for Nt
    # Define function g: N = N0 * exp(-M) + P * exp(-M)
    # Variables: M (x1), N0 (x2), P (x3)
    g <- ~ x2 * exp(-x1) + x3 * exp(-x1)
    
    mean_vals <- c(M, N0, Pt)
    cov_sub <- cov_mat[c("M", "N0", "P1"), c("M", "N0", "P1")]
    
    # SE for abundance
    N_se <- deltamethod(g, mean = mean_vals, cov = cov_sub)
    N_thou_se[i] <- N_se
    
    # Biomass in tonnes
    B_ton[i] <- 1000 * Nt * weight_kg_mean * 0.001
    
    # SE for biomass using error propagation
    B_ton_se[i] <- sqrt(
      (1000 * N_se)^2 * (weight_kg_mean * 0.001)^2 +
        (1000 * Nt)^2 * (weight_kg_sd * 0.001)^2
    )
  }
  
  return(data.frame(
    timestep = seq_along(timing),
    N_thou = N_thou,
    N_thou_se = N_thou_se,
    B_ton = B_ton,
    B_ton_se = B_ton_se
  ))
}

