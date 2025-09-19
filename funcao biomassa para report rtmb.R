

x = cobaia
method = "spg"
multi = T
mbw.sd = predicos
nmult = 1000
dates = cat_df$Properties$Dates
df = cat_df
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
  PopDyn[2] <- PopDyn[1] * mean(unlist(x[[method]]$bt.stdev[which(!is.na(x[[method]]$bt.stdev))])/
                                  unlist(x[[method]]$bt.par[which(!is.na(x[[method]]$bt.stdev))]))
}
if (is.na(PopDyn[4])) {
  PopDyn[4] <- PopDyn[3] * mean(unlist(x[[method]]$bt.stdev[which(!is.na(x[[method]]$bt.stdev))])/
                                  unlist(x[[method]]$bt.par[which(!is.na(x[[method]]$bt.stdev))]))
}

    Perts = data.frame(Pest = unlist(x[[method]]$bt.par[grep("P", names(x[[method]]$bt.par))]) * nmult, 
                        SE.Pest = unlist(x[[method]]$bt.stdev[grep("P.",names(x[[method]]$bt.par))]) * nmult, 
                        tsteps = x[[method]]$Dates[-c(1,length(x[[method]]$Dates))])
# se houver NAs nos stdev, corrigir aqui. Ver CatdynBSD() 
    if (any(is.na(Perts$SE.Pest))) {
      Perts$SE.Pest[which(is.na(Perts$SE.Pest))] <- Perts$Pest[which(is.na(Perts$SE.Pest))] * 
        mean(Perts$SE.Pest[which(!is.na(Perts$SE.Pest))]/Perts$Pest[which(!is.na(Perts$SE.Pest))])
    }
    
   
    
    mt = x[[method]]$Type # nr de perturbacoes
    
     

    if (length(x$Data$Properties$Fleets$Fleet) == 1) {
      mt <- x[[method]]$Type
      
      Timing_2 = matrix(0, nrow = 12 * mt, ncol = 1)
      for (j in 1:mt) {
        Timing_2[Perts$tsteps[j]:(j*12)] = 1 # .
      }  
      
      
      Cov.Mat <- cor2cov(cor.mat = x[[method]]$Cor[c(1:(mt + 2)),
                                                   c(1:(mt + 2))],
                         sd = c(PopDyn$SE.M, PopDyn$SE.N0,
                                Perts$SE.Pest))
      
      
      
      yr1 <- as.numeric(format(as.Date(dates[1]), 
                               "%Y"))
      yr2 <- as.numeric(format(as.Date(dates[2]), 
                               "%Y"))
      class(x) = "catdyn"
      x$Model = x
      x$Data = df
      z <- CatDynPred(x, method)
      

      
      PredStock <- data.frame(Year = sort(rep(yr1:yr2, 
                                              12)), Month = c("Jan", "Feb", "Mar", "Apr", 
                                                              "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", 
                                                              "Dec"), TimeStep = 1:((yr2 - yr1 + 1) * 12), 
                              Mmw.kg = x$Data$Data[[fleet1]]$obsmbw.kg, SDmw.kg = mbw.sd, 
                              N.thou = z$Model$Results[, 10], N.thou.SE = 0, 
                              B.ton = z$Model$Results[, 11], B.ton.SE = 0)
      for (m in 1:12) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PopDyn$N0, Timing[m] * Perts$Pest[1]), cov = Cov.Mat[c(1:3), 
                                                                                                                                                  c(1:3)])
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 13:24) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[12], Timing[m] * Perts$Pest[2]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       4], 0, PredStock$N.thou.SE[12]^2, 0, Cov.Mat[4, 
                                                                                                                                    1], 0, Cov.Mat[4, 4]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 25:36) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[24], Timing[m] * Perts$Pest[3]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       5], 0, PredStock$N.thou.SE[24]^2, 0, Cov.Mat[5, 
                                                                                                                                    1], 0, Cov.Mat[5, 5]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 37:48) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[36], Timing[m] * Perts$Pest[4]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       6], 0, PredStock$N.thou.SE[36]^2, 0, Cov.Mat[6, 
                                                                                                                                    1], 0, Cov.Mat[6, 6]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 49:60) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[48], Timing[m] * Perts$Pest[5]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       7], 0, PredStock$N.thou.SE[48]^2, 0, Cov.Mat[7, 
                                                                                                                                    1], 0, Cov.Mat[7, 7]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 61:72) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[60], Timing[m] * Perts$Pest[6]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       8], 0, PredStock$N.thou.SE[60]^2, 0, Cov.Mat[8, 
                                                                                                                                    1], 0, Cov.Mat[8, 8]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 73:84) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[72], Timing[m] * Perts$Pest[7]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       9], 0, PredStock$N.thou.SE[72]^2, 0, Cov.Mat[9, 
                                                                                                                                    1], 0, Cov.Mat[9, 9]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 85:96) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[84], Timing[m] * Perts$Pest[8]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       10], 0, PredStock$N.thou.SE[84]^2, 0, Cov.Mat[10, 
                                                                                                                                     1], 0, Cov.Mat[10, 10]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 97:108) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[96], Timing[m] * Perts$Pest[9]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       11], 0, PredStock$N.thou.SE[96]^2, 0, Cov.Mat[11, 
                                                                                                                                     1], 0, Cov.Mat[11, 11]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 109:120) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[108], Timing[m] * Perts$Pest[10]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       12], 0, PredStock$N.thou.SE[108]^2, 0, Cov.Mat[12, 
                                                                                                                                      1], 0, Cov.Mat[12, 12]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 121:132) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[120], Timing[m] * Perts$Pest[11]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       13], 0, PredStock$N.thou.SE[120]^2, 0, Cov.Mat[13, 
                                                                                                                                      1], 0, Cov.Mat[13, 13]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 133:144) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[132], Timing[m] * Perts$Pest[12]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       14], 0, PredStock$N.thou.SE[132]^2, 0, Cov.Mat[14, 
                                                                                                                                      1], 0, Cov.Mat[14, 14]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 145:156) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[144], Timing[m] * Perts$Pest[13]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       15], 0, PredStock$N.thou.SE[144]^2, 0, Cov.Mat[15, 
                                                                                                                                      1], 0, Cov.Mat[15, 15]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 157:168) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[156], Timing[m] * Perts$Pest[14]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       16], 0, PredStock$N.thou.SE[156]^2, 0, Cov.Mat[16, 
                                                                                                                                      1], 0, Cov.Mat[16, 16]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
      }
      for (m in 169:180) {
        PredStock$N.thou.SE[m] <- deltamethod(g = list(~x2 * 
                                                         exp(-x1) + x3 * exp(-x1)), mean = c(PopDyn$M, 
                                                                                             PredStock$N.thou[168], Timing[m] * Perts$Pest[15]), 
                                              cov = matrix(c(Cov.Mat[1, 1], 0, Cov.Mat[1, 
                                                                                       17], 0, PredStock$N.thou.SE[168]^2, 0, Cov.Mat[17, 
                                                                                                                                      1], 0, Cov.Mat[17, 17]), 3, 3))
        PredStock$B.ton.SE[m] <- sqrt((1000 * PredStock$N.thou.SE[m])^2 * 
                                        (PredStock$Mmw.kg[m] * 0.001)^2 + (1000 * 
                                                                             PredStock$N.thou[m])^2 * (PredStock$SDmw.kg[m] * 
                                                                                                         0.001)^2)
    
      }
    }
    
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