

x = cobaia_ref
z = obj$report()$nstep
method = "spg"
multi = T
mbw.sd = predicos
nmult = 1000
dates = cat_df$Properties$Dates
df = cat_df
p=15

# check 1
ggplot() +
  geom_point(aes(x = 1:180,
                 y = obj$report()$nstep)) +
  geom_point(aes(x = 1:180,
                 y = annual_biomass$N.thou), col = 'red')

# check 2 - thousands of individuals to tons
ggplot() +
  geom_point(aes(x = 1:180,
                 y = annual_biomass$N.thou * nmult /1000 *
              cat_df$Data$`MIS+OTB-S`$obsmbw.kg)) +

  geom_point(aes(x = 1:180,
                 y = annual_biomass$B.ton), col = 'red')

# rtmb
ggplot() +
  geom_point(aes(x = 1:180,
                 y = obj$report()$nstep * nmult /1000 *
                   cat_df$Data$`MIS+OTB-S`$obsmbw.kg)) +
  
  geom_point(aes(x = 1:180,
                 y = annual_biomass$B.ton), col = 'red')




effort_df = cat_df$Data$`MIS+OTB-S` %>% 
  transmute(time = time.step,
            fleet = "MIS+OTB-S",
            effort = obseff.trips)

rep = cobaia_ref[["spg"]]


manual_biomass_estimate <- function(report_out, effort_df, method = "spg") {
  # Extract from report
  rep <- report_out[[method]]
  bt.par <- rep$bt.par
  dates <- rep$Dates
  fleets <- gsub("^P\\d+\\.", "", grep("^P", names(bt.par), value = TRUE))
  p <- rep$Type  # number of time blocks
  
  # Initial abundance and natural mortality
  N0 <- bt.par["N0"]
  M <- bt.par["M"]
  
  # Get catchability (P1.FleetA, P2.FleetA, etc.)
  P_vals <- bt.par[grep("^P\\d+\\.", names(bt.par))]
  names(P_vals) <- gsub("^P\\d+\\.", "", names(P_vals))  # simplify names to just fleet
  
  # Time structure
  # time_seq <- seq(as.Date(dates["ts.start"]), as.Date(dates["ts.end"]), by = "week")
  n_time <- max(effort_df$time)
  
  period_indices <- ceiling((1:n_time) / (n_time / p))
  biomass <- numeric(n_time)
  biomass[1] <- N0
  
  # Main loop
  for (t in 2:n_time) {
    F_total <- 0
    
    # Effort for previous time step (t-1)
    day_effort <- effort_df[effort_df$time == (t - 1), ]
    
    for (fleet in fleets) {
      fleet_effort <- sum(day_effort$effort[day_effort$fleet == fleet])
      period <- period_indices[t - 1]
      P_name <- paste0("P", period, ".", fleet)
      q <- bt.par[P_name]
      F_total <- F_total + q * fleet_effort
    }
    
    Z <- M + F_total
    biomass[t] <- biomass[t - 1] * exp(-Z)
  }
  
  # Return results
  return(data.frame(
    timestep = 1:n_time,
    biomass = unlist(biomass)
  ))
}

res = manual_biomass_estimate(cobaia_ref, effort_df)



# scaler =    Thou.scaler <- 1e+06 * (x$Data$Properties$Units[4] == 
#                                       "bill") + 1000 * (x$Data$Properties$Units[4] == 
  #                                                         "mill") + 1 * (x$Data$Properties$Units[4] == "thou") + 
  # 0.1 * (x$Data$Properties$Units[4] == "hund")



function (x, df, dates, 
          method, multi, mbw.sd) 

PopDyn  = data.frame(M = x[[method]]$bt.par$M, 
                     SE.M = x[[method]]$bt.stdev["M"],
                     N0 = nmult * x[[method]]$bt.par$N0,
                     SE.N0 = nmult * x[[method]]$bt.stdev["N0"])

# se houver NAs nos stdev, corrigir aqui. Ver CatdynBSD()
if (is.na(PopDyn[2])) {
  PopDyn[2] <- PopDyn[1] * mean(unlist(x[[method]]$bt.stdev[
    which(!is.na(x[[method]]$bt.stdev))])/
                                  unlist(x[[method]]$bt.par[
                                    which(!is.na(x[[method]]$bt.stdev))]))
}

if (is.na(PopDyn[4])) {
  PopDyn[4] <- PopDyn[3] * mean(unlist(x[[method]]$bt.stdev[
    which(!is.na(x[[method]]$bt.stdev))])/
                                  unlist(x[[method]]$bt.par[
                                    which(!is.na(x[[method]]$bt.stdev))]))
}

    Perts = data.frame(Pest = 
                         unlist(x[[method]]$bt.par[
                           grep("P",names(x[[method]]$bt.par))]) * nmult, 
                        SE.Pest = 
                         unlist(x[[method]]$bt.stdev[
                           grep("P.",names(x[[method]]$bt.par))]) * nmult, 
                        tsteps = 
                         x[[method]]$Dates[-c(1,length(x[[method]]$Dates))])
    
# se houver NAs nos stdev, corrigir aqui. Ver CatdynBSD() 
    if (any(is.na(Perts$SE.Pest))) {
      Perts$SE.Pest[which(is.na(Perts$SE.Pest))] =
        Perts$Pest[which(is.na(Perts$SE.Pest))] * 
        mean(Perts$SE.Pest[which(!is.na(Perts$SE.Pest))]/
               Perts$Pest[which(!is.na(Perts$SE.Pest))])
    }

    mt = x[[method]]$Type # nr de perturbacoes
    
    # if (length(x$Data$Properties$Fleets$Fleet) == 1) {

      
      Timing = matrix(0, nrow = 12 * mt, ncol = 1)
      for (j in 1:mt) {
        Timing[Perts$tsteps[j]:(j*12)] = 1 # .
      }  
      
      
      Cov.Mat <- cor2cov(cor.mat = x[[method]]$Cor[c(1:(mt + 2)),
                                                   c(1:(mt + 2))],
                         sd = c(PopDyn$SE.M,
                                PopDyn$SE.N0,
                                Perts$SE.Pest))
      
      yr1 <- as.numeric(format(as.Date(dates[1]), 
                               "%Y"))
      yr2 <- as.numeric(format(as.Date(dates[2]), 
                               "%Y"))
      # class(x) = "catdyn"
      # x$Model = x
      # x$Data = df
      
      PredStock <- data.frame(
        Year = sort(rep(yr1:yr2,12)),
        Month = c("Jan", "Feb", "Mar", "Apr",
                  "May", "Jun", "Jul", "Aug", 
                  "Sep", "Oct", "Nov", "Dec"), 
        TimeStep = 1:((yr2 - yr1 + 1) * 12), 
        Mmw.kg = df$Data[['MIS+OTB-S']]$obsmbw.kg,
        SDmw.kg = predicos$se.fit, 
        N.thou = z,
        N.thou.SE = 0, 
        B.ton.SE = 0) %>% 
        mutate(B.ton = z * nmult /1000 * Mmw.kg,) # 1000 = ton scaler
      
      for (i in 1:length(z)) {
        m = ceiling(i/12)
        
        matN =  matrix(
          c(Cov.Mat[1, 1], 0,
            Cov.Mat[1, 2+m], 0,
            PredStock$N.thou.SE[(m-1)*12]^2, 0,
            Cov.Mat[2+m, 1], 0,
            Cov.Mat[2+m, 2+m]), 3, 3)
        
        if(m==1){cova = Cov.Mat[c(1:3), c(1:3)]
        }else{cova = matN}
        
        PredStock$N.thou.SE[i] =
          deltamethod(g = list(~x2 * exp(-x1) + x3 * exp(-x1)),
                      mean = c(PopDyn$M,
                               PopDyn$N0,
                               Timing[i] * Perts$Pest[m]),
                      cov = cova)
                      
        PredStock$B.ton.SE[i] <- 
          sqrt((nmult * PredStock$N.thou.SE[i] )^2 * 
                      (PredStock$Mmw.kg[i] * .001)^2 + 
                 (nmult * PredStock$N.thou[i] )^2 * 
                 (PredStock$SDmw.kg[i] * .001)^2)
      }
      
      PredStock$B.ton/annual_biomass$B.ton
      PredStock$B.ton.SE /annual_biomass$B.ton.SE     
      PredStock$N.thou.SE/annual_biomass$N.thou.SE
###############
    
  }
  if (mt == 4) {
    P1 <- Thou.scaler * unlist(x[[i]]$Model[[method[i]]]$bt.par[3])
    P1.SE <- Thou.scaler * unlist(x[[i]]$Model[[method[i]]]$bt.stdev[3])
    P1.back <- x[[i]]$Model[[method[i]]]$Dates[2] - 
      x[[i]]$Model[[method[i]]]$Dates[1] + 1
    P2 <- Thou.scaler * unlist(x[[i]]$Model[[method[i]]]$bt.par[4])
    P2.SE <- Thou.scaler * unlist(x[[i]]$Model[[method[i]]]$bt.stdev[4])
    P2.back <- x[[i]]$Model[[method[i]]]$Dates[3] - 
      x[[i]]$Model[[method[i]]]$Dates[1] + 1
    P3 <- Thou.scaler * unlist(x[[i]]$Model[[method[i]]]$bt.par[5])
    P3.SE <- Thou.scaler * unlist(x[[i]]$Model[[method[i]]]$bt.stdev[5])
    P3.back <- x[[i]]$Model[[method[i]]]$Dates[4] - 
      x[[i]]$Model[[method[i]]]$Dates[1] + 1
    P4 <- Thou.scaler * unlist(x[[i]]$Model[[method[i]]]$bt.par[6])
    P4.SE <- Thou.scaler * unlist(x[[i]]$Model[[method[i]]]$bt.stdev[6])
    P4.back <- x[[i]]$Model[[method[i]]]$Dates[5] - 
      x[[i]]$Model[[method[i]]]$Dates[1] + 1
    cov <- cor2cov(x[[i]]$Model[[method[i]]]$Cor[1:6, 
                                                 1:6], sd = c(SE.M, SE.N0, P1.SE, P2.SE, 
                                                              P3.SE, P4.SE))
    PredStock$N0Tot.thou[i] <- N0 + P1 * exp(M * 
                                               P1.back) + P2 * exp(M * P2.back) + P3 * 
      exp(M * P3.back) + P4 * exp(M * P4.back)
    form <- sprintf("~x2+x3*exp(x1*%i)+x4*exp(x1*%i)+x5*exp(x1*%i)+x6*exp(x1*%i)", 
                    P1.back, P2.back, P3.back, P4.back)
    PredStock$N0Tot.thou.SE[i] <- deltamethod(g = as.formula(form), 
                                              mean = c(M, N0, P1, P2, P3, P4), cov = cov)
    PredStock$B0Tot.ton[i] <- PredStock$N0Tot.thou[i] * 
      mbw.sd[i, 2]
    PredStock$B0Tot.ton.SE[i] <- sqrt((PredStock$N0Tot.thou.SE[i])^2 * 
                                        (mbw.sd[i, 2])^2 + (PredStock$N0Tot.thou[i])^2 * 
                                        (mbw.sd[i, 3])^2)    
    
    