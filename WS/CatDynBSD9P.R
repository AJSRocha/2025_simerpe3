CatDynBSD9P <- function (x, method, multi, mbw.sd)
{
  if (!multi)
  {
    stop("This function only works for multi-annual depletion models")
  }
  if (multi) {
        if (class(x) != "catdyn") 
            {
            stop("In multi-annual models 'x' must be a single object of class 'catdyn' run at monthly time steps")
            }
        if (x$Data$Properties$Units[3] == "ind") 
            {
            stop("This function is used to calculate standard deviation of annual biomass when the catch is recorded in weight")
            }
        if (x$Data$Properties$Units["Time Step"] != "month") 
            {
            stop("In multi-annual models 'x' must be a single object of class 'catdyn' run at monthly time steps")
            }
        if (length(x$Data$Properties$Fleets$Fleet) != 1)
            {
            stop("This function only works for single-fleet depletion models")
            }
        if (abs(x$Model[[method]]$Type) != 9) 
            {
            stop("This function only works with 9 years of data")
            }
        if (length(mbw.sd) != 12 & length(mbw.sd) != 12 * abs(x$Model[[method]]$Type)) 
            {
            stop("mbw.sd must be a vector of length 12 (monthly mean weight) or 12*number of years (in kg)")
            }
        Thou.scaler <- 1000000 * (x$Data$Properties$Units[4] ==  "bill") + 
                       1000 * (x$Data$Properties$Units[4] == "mill") + 
                       1 * (x$Data$Properties$Units[4] == "thou") + 
                       0.1 * (x$Data$Properties$Units[4] == "hund")
        PopDyn <- data.frame(M = x$Model[[method]]$bt.par$M,
                             SE.M = x$Model[[method]]$bt.stdev["M"], 
                             N0 = Thou.scaler * x$Model[[method]]$bt.par$N0, 
                             SE.N0 = Thou.scaler * x$Model[[method]]$bt.stdev["N0"])
        if (is.na(PopDyn[2])) 
            {
            PopDyn[2] <- PopDyn[1] * mean(unlist(x$Model[[method]]$bt.stdev[which(!is.na(x$Model[[method]]$bt.stdev))])/unlist(x$Model[[method]]$bt.par[which(!is.na(x$Model[[method]]$bt.stdev))]))
            }
        if (is.na(PopDyn[4])) 
            {
            PopDyn[4] <- PopDyn[3] * mean(unlist(x$Model[[method]]$bt.stdev[which(!is.na(x$Model[[method]]$bt.stdev))])/unlist(x$Model[[method]]$bt.par[which(!is.na(x$Model[[method]]$bt.stdev))]))
            }
        Perts <- data.frame(Pest = unlist(x$Model[[method]]$bt.par[3:(abs(x$Model[[method]]$Type)+2)]) * Thou.scaler, 
                            SE.Pest = unlist(x$Model[[method]]$bt.stdev[3:(abs(x$Model[[method]]$Type)+2)]) * Thou.scaler, 
                            tsteps = x$Model[[method]]$Dates[grep("ts.P",names(x$Model[[method]]$Dates))])
        if (any(is.na(Perts$SE.Pest)))
            {
            Perts$SE.Pest[which(is.na(Perts$SE.Pest))] <- Perts$Pest[which(is.na(Perts$SE.Pest))] *
                mean(Perts$SE.Pest[which(!is.na(Perts$SE.Pest))]/Perts$Pest[which(!is.na(Perts$SE.Pest))])
            }
        if (length(x$Data$Properties$Fleets$Fleet) == 1) 
            {
            mt <- abs(x$Model[[method]]$Type)
            Timing <- matrix(0, 12 * mt, 1)
            Timing[1:12] <- ifelse(row(Timing)[1:12] >= Perts$tsteps[1],1, 0)
            Timing[13:24] <- ifelse(row(Timing)[13:24] >= Perts$tsteps[2],1, 0)
            Timing[25:36] <- ifelse(row(Timing)[25:36] >= Perts$tsteps[3],1, 0)
            Timing[37:48] <- ifelse(row(Timing)[37:48] >= Perts$tsteps[4],1, 0)
            Timing[49:60] <- ifelse(row(Timing)[49:60] >= Perts$tsteps[5],1, 0)
            Timing[61:72] <- ifelse(row(Timing)[61:72] >= Perts$tsteps[6],1, 0)
            Timing[73:84] <- ifelse(row(Timing)[73:84] >= Perts$tsteps[7],1, 0)
            Timing[85:96] <- ifelse(row(Timing)[85:96] >= Perts$tsteps[8],1, 0)
            Timing[97:108] <- ifelse(row(Timing)[97:108] >= Perts$tsteps[9],1, 0)
            #Timing[109:120] <- ifelse(row(Timing)[109:120] >= Perts$tsteps[10], 1, 0)
            #Timing[121:132] <- ifelse(row(Timing)[121:132] >= Perts$tsteps[11], 1, 0)
            # Timing[133:144] <- ifelse(row(Timing)[133:144] >= Perts$tsteps[12], 1, 0)
            # Timing[145:156] <- ifelse(row(Timing)[145:156] >= Perts$tsteps[13], 1, 0)
            # Timing[157:168] <- ifelse(row(Timing)[157:168] >= Perts$tsteps[14], 1, 0)
            # Timing[169:180] <- ifelse(row(Timing)[169:180] >= Perts$tsteps[15], 1, 0)
            fleet1 <- x$Data$Properties$Fleets[1, 1]
            Cov.Mat <- cor2cov(cor.mat = x$Model[[method]]$Cor[c(1:(mt + 2)), c(1:(mt + 2))], 
                               sd = c(PopDyn$SE.M, PopDyn$SE.N0,Perts$SE.Pest))
            yr1 <- as.numeric(format(as.Date(x$Data$Properties$Dates[1]),"%Y"))
            yr2 <- as.numeric(format(as.Date(x$Data$Properties$Dates[2]),"%Y"))
            z <- CatDynPred(x, method)
            PredStock <- data.frame(Year = sort(rep(yr1:yr2,12)), 
                                    Month = c("Jan", "Feb", "Mar","Apr", "May", "Jun", "Jul","Aug", "Sep", "Oct", "Nov","Dec"), 
                                    TimeStep = 1:((yr2 - yr1 + 1) * 12), 
                                    Mmw.kg = x$Data$Data[[fleet1]]$obsmbw.kg,  
                                    SDmw.kg = mbw.sd, N.thou = z$Model$Results[, 10], 
                                    N.thou.SE = 0, 
                                    B.ton = z$Model$Results[,11], B.ton.SE = 0)
            for (m in 1:12) {
                PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * exp(-x1) + x3 * exp(-x1)), 
                                                      mean = c(PopDyn$M, PopDyn$N0, Timing[m] * Perts$Pest[1]), 
                                                      cov = Cov.Mat[c(1:3), c(1:3)])
                PredStock$B.ton.SE[m]  <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 0.001)^2)
            }
            for (m in 13:24) {
                PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * exp(-x1) + x3 * exp(-x1)), 
                                                      mean = c(PopDyn$M, PredStock$N.thou[12], Timing[m] * Perts$Pest[2]),
                                                      cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1,4], 0, PredStock$N.thou.SE[12]^2, 0, Cov.Mat[4,1], 0, Cov.Mat[4, 4]), 3, 3))
                PredStock$B.ton.SE[m]  <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 0.001)^2)
            }
            for (m in 25:36) {
                PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * exp(-x1) + x3 * exp(-x1)), 
                                                      mean = c(PopDyn$M, PredStock$N.thou[24], Timing[m] * Perts$Pest[3]),
                                                      cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1,5], 0, PredStock$N.thou.SE[24]^2, 0, Cov.Mat[5,1], 0, Cov.Mat[5, 5]), 3, 3))
                PredStock$B.ton.SE[m]  <- sqrt((1000 * PredStock$N.thou.SE[m])^2 *(PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 0.001)^2)
            }
            for (m in 37:48) {
                PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * exp(-x1) + x3 * exp(-x1)), 
                                                      mean = c(PopDyn$M,PredStock$N.thou[36], Timing[m] * Perts$Pest[4]),
                                                      cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1,6], 0, PredStock$N.thou.SE[36]^2, 0, Cov.Mat[6,1], 0, Cov.Mat[6, 6]), 3, 3))
                PredStock$B.ton.SE[m]  <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 0.001)^2)
            }
            for (m in 49:60) {
                 PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * exp(-x1) + x3 * exp(-x1)), 
                                                       mean = c(PopDyn$M, PredStock$N.thou[48], Timing[m] * Perts$Pest[5]),
                                                       cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 7], 0, PredStock$N.thou.SE[48]^2, 0, Cov.Mat[7,1], 0, Cov.Mat[7, 7]), 3, 3))
                 PredStock$B.ton.SE[m]  <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 0.001)^2)
                }
            for (m in 61:72) {
                PredStock$N.thou.SE[m]  <- deltamethod(g = list(~x2 * exp(-x1) + x3 * exp(-x1)), 
                                                       mean = c(PopDyn$M, PredStock$N.thou[60], Timing[m] * Perts$Pest[6]),
                                                       cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 8], 0, PredStock$N.thou.SE[60]^2, 0, Cov.Mat[8,1], 0, Cov.Mat[8, 8]), 3, 3))
                PredStock$B.ton.SE[m]   <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 0.001)^2)
            }
            for (m in 73:84) {
                PredStock$N.thou.SE[m]  <- deltamethod(g = list(~x2 * exp(-x1) + x3 * exp(-x1)), 
                                                       mean = c(PopDyn$M, PredStock$N.thou[72], Timing[m] * Perts$Pest[7]),
                                                       cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 9], 0, PredStock$N.thou.SE[72]^2, 0, Cov.Mat[9,1], 0, Cov.Mat[9, 9]), 3, 3))
                PredStock$B.ton.SE[m]   <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 0.001)^2)
            }
            for (m in 85:96) {
                PredStock$N.thou.SE[m]  <- deltamethod(g = list(~x2 * exp(-x1) + x3 * exp(-x1)), 
                                                       mean = c(PopDyn$M, PredStock$N.thou[84], Timing[m] * Perts$Pest[8]),
                                                       cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 10], 0, PredStock$N.thou.SE[84]^2, 0, Cov.Mat[10,1], 0, Cov.Mat[10, 10]), 3, 3))
                PredStock$B.ton.SE[m]   <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 0.001)^2)
            }
            for (m in 97:108) {
                PredStock$N.thou.SE[m]  <- deltamethod(g = list(~x2 * exp(-x1) + x3 * exp(-x1)), 
                                                       mean = c(PopDyn$M, PredStock$N.thou[96], Timing[m] * Perts$Pest[9]),
                                                       cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 11], 0, PredStock$N.thou.SE[96]^2, 0, Cov.Mat[11,1], 0, Cov.Mat[11, 11]), 3, 3))
                PredStock$B.ton.SE[m]   <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 0.001)^2)
            }
         return(PredStock)
        }
    }
}
