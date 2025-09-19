
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


