

dat = list()
dat$Ct = as.vector(tibec.14.22.pg.2$Data$Pelagic$obscat.thou)         # observed catch
dat$Et = as.vector(tibec.14.22.pg.2$Data$Pelagic$obseff.ntrips)       # effort
dat$u = c(12, 24, 36, 47, 60, 72, 84, 96, 107)                        # recruitment pulse months
nT = length(dat$Ct)
nR = length(dat$u)

I = matrix(0, nrow = nT, ncol = nR)
for (j in 1:nR) {
  I[dat$u[j]:nT, j] <- 1
}

dat$I = I

par_cd = tibec.14.22.pg.apn.1.fit$Model$CG$bt.par

# Initial parameter values
par = list(
  Rt_scaled     = c(par_cd$P1.Pelagic, par_cd$P2.Pelagic, par_cd$P3.Pelagic, par_cd$P4.Pelagic,
                    par_cd$P5.Pelagic, par_cd$P6.Pelagic, par_cd$P7.Pelagic, par_cd$P8.Pelagic,
                    par_cd$P9.Pelagic),
  logalpha      = log(par_cd$alpha.Pelagic),
  logbeta       = log(par_cd$beta.Pelagic),
  logK          = log(par_cd$k.Pelagic),
  logN0_scaled  = log(par_cd$N0),
  logM          = log(par_cd$M),
  logsdCt       = log(0.25 * sd(dat$Ct[dat$Ct > 0]))  # consistent with CatDyn
)

attach(dat)
attach(par)

# initialize joint negative loglikelihood function


  # Extract parameters
  Ct     = OBS(Ct)
  Et     = OBS(Et)
  I_mat  = OBS(I)
  alpha  = exp(logalpha)
  beta   = exp(logbeta)
  K      = exp(logK)
  N0     = exp(logN0_scaled)
  M      = exp(logM)
  sdCt   = 0.7032531#exp(logsdCt)
  Rt     = Rt_scaled
  
  
  # Initialize predicted catch and biomass
  nT     = length(Ct)
  catch = rep(0, nT) * 1   # Ensures catch becomes an AD type
  Bt = rep(0, nT) * 1     # Ensures Bt becomes an AD type
  # Bt[1]  = N0  # start biomass at N0
  m = exp(-M/2)
  # Forward simulation of catch and biomass
  # for (t in 1:nT) {
  #   # Predict catch at time t
  #   catch[t] = K * exp(-M/2) *  Et[t]^alpha * Bt[t]^beta
  #   
  #   # Update biomass if not last timestep
  #   if (t < nT) {
  #     # Rt_t = sum(I_mat[t + 1, ] * Rt)
  #     for(j in 1:nR){R_t = sum(I_mat[t+1,j] * Rt[j]) * exp(-M *(t-j))}
  #     Bt[t + 1] = Bt[t] * exp(-M) - catch[t] * exp(-M / 2) + R_t
  #   }
  # }
  # Bt
  # Negative log-likelihood
  # jnll = -sum(dnorm(Ct, mean = catch, sd = sdCt, log = TRUE))
  
  for(mes in 1:length(Ct)){
    # core
      if(mes == 1){
      Bt[mes]  = N0*exp(-M)
      
    }else{
      
           # M-driven decay
      
      part1 = N0*exp(-M*mes)
      
      
      part2 = 0
      
      for(mes_ant in 1:(mes-1)){
        part2 = part2 + catch[mes_ant]*exp(-M*(mes-mes_ant-1))}
      
      
      # Recruitment pulses
      
      part3 = 0
          part3 = sum(I_mat[mes,] * Rt)
    
      part4 = 0
      
      # Assemble
      Bt[mes] = part1 - m * part2 + part3 - part4}
    
      pred =  (K*Et[mes]^alpha) * m * ((Bt[mes])^beta)
      catch[mes] = pred
   }
  
  

catch_cat = unlist(tibec.14.22.pg.apn.1.fit.pred.CG$Model$Results['Predicted.Catch.thou'])
catch


ggplot() + 
  geom_line(aes(x = 1:108,
                y = catch_cat)) + 
  # geom_line(aes(x = 1:108,
                # y = rnorm(108,catch,sdCt)), col = 'red') + 
  geom_line(aes(x = 1:108,
                y = catch), col = 'green') + 
  theme_bw()


testador(par = c(log(par_cd$M),
                 log(par_cd$N0),
                 log(par_cd$P1.Pelagic),
                 log(par_cd$P2.Pelagic),
                 log(par_cd$P3.Pelagic),
                 log(par_cd$P4.Pelagic),
                 log(par_cd$P5.Pelagic),
                 log(par_cd$P6.Pelagic),
                 log(par_cd$P7.Pelagic),
                 log(par_cd$P8.Pelagic),
                 log(par_cd$P9.Pelagic),
                 log(par_cd$k.Pelagic),
                 log(par_cd$alpha.Pelagic),
                 log(par_cd$beta.Pelagic),
                 0.25*sd(tibec.14.22.pg.2$Data$Pelagic$obscat.thou[tibec.14.22.pg.2$Data$Pelagic$obscat.thou>0])^2,
                 0.25*sd(log(tibec.14.22.pg.2$Data$Pelagic$obscat.thou[tibec.14.22.pg.2$Data$Pelagic$obscat.thou>0]))^2),
         dates = tibec.14.22.pg.dates.1,
         obscat1 = Ct,
         obseff1 = Et,
         obsmbm1 = tibec.14.22.pg.1[,5],
         output = 'no',
         partial = F)

# funcao extraida
# 
testador = function (par, dates, obscat1, obseff1, obsmbm1, distr, properties, 
                     output, pp = 9, partial) 
{
  ts.start <- head(dates, 1)
  ts.P1 <- dates[2]
  ts.P2 <- dates[3]
  ts.P3 <- dates[4]
  ts.P4 <- dates[5]
  ts.P5 <- dates[6]
  ts.P6 <- dates[7]
  ts.P7 <- dates[8]
  ts.P8 <- dates[9]
  ts.P9 <- dates[10]
  ts.end <- tail(dates, 1)
  sealen <- ts.end - ts.start + 1 #length of series
  nstep <- vector("numeric", sealen) #init vector
  mccum <- vector("numeric", sealen)
  effeff1 <- vector("numeric", sealen)
  effn1 <- vector("numeric", sealen)
  predcat1 <- vector("numeric", sealen)
  
  #indicator vectors
  ind.P1 <- ifelse(1:sealen < (ts.P1 - ts.start), 0, 1)
  ind.P2 <- ifelse(1:sealen < (ts.P2 - ts.start), 0, 1)
  ind.P3 <- ifelse(1:sealen < (ts.P3 - ts.start), 0, 1)
  ind.P4 <- ifelse(1:sealen < (ts.P4 - ts.start), 0, 1)
  ind.P5 <- ifelse(1:sealen < (ts.P5 - ts.start), 0, 1)
  ind.P6 <- ifelse(1:sealen < (ts.P6 - ts.start), 0, 1)
  ind.P7 <- ifelse(1:sealen < (ts.P7 - ts.start), 0, 1)
  ind.P8 <- ifelse(1:sealen < (ts.P8 - ts.start), 0, 1)
  ind.P9 <- ifelse(1:sealen < (ts.P9 - ts.start), 0, 1)
  
  logM <- par[1]
  logN0 <- par[2]
  logP1 <- par[3]
  logP2 <- par[4]
  logP3 <- par[5]
  logP4 <- par[6]
  logP5 <- par[7]
  logP6 <- par[8]
  logP7 <- par[9]
  logP8 <- par[10]
  logP9 <- par[11]
  logscale <- par[12]
  logalpha <- par[13]
  logbeta <- par[14]
  mccum[1] <- 0
  nstep[1] <- exp(logN0) * exp(-exp(logM))
  for (i in 2:sealen) {
    mccum[i] <- obscat1[i - 1] + mccum[i - 1] * exp(-exp(logM))
    nstep[i] <- exp(logN0) * exp(-exp(logM) * i) + ind.P1[i] * 
      exp(logP1) * exp(-exp(logM) * (i - (ts.P1 - ts.start) + 
                                       1)) + ind.P2[i] * exp(logP2) * exp(-exp(logM) * (i - 
                                                                                          (ts.P2 - ts.start) + 1)) + ind.P3[i] * exp(logP3) * 
      exp(-exp(logM) * (i - (ts.P3 - ts.start) + 1)) + 
      ind.P4[i] * exp(logP4) * exp(-exp(logM) * (i - (ts.P4 - 
                                                        ts.start) + 1)) + ind.P5[i] * exp(logP5) * exp(-exp(logM) * 
                                                                                                         (i - (ts.P5 - ts.start) + 1)) + ind.P6[i] * exp(logP6) * 
      exp(-exp(logM) * (i - (ts.P6 - ts.start) + 1)) + 
      ind.P7[i] * exp(logP7) * exp(-exp(logM) * (i - (ts.P7 - 
                                                        ts.start) + 1)) + ind.P8[i] * exp(logP8) * exp(-exp(logM) * 
                                                                                                         (i - (ts.P8 - ts.start) + 1)) + ind.P9[i] * exp(logP9) * 
      exp(-exp(logM) * (i - (ts.P9 - ts.start) + 1)) - 
      mccum[i] * exp(-exp(logM)/2)
  }
  effeff1 <- obseff1^(exp(logalpha))
  effn1 <- nstep^(exp(logbeta))
  predcat1 <- exp(logscale) * (effeff1 * effn1) * exp(-exp(logM)/2)
  
  
  
  Likel <- .CatDynLik1F(obscat1, predcat1, distr, par)
  if (output == "predict") {
    catdynexp <- .CatDynExp1F.Res(properties, nstep, obsmbm1, 
                                  pp, dates, distr, par, Likel, ts.start, ts.end, obseff1, 
                                  obscat1, predcat1, partial)
    class(catdynexp) <- "CatDynExp"
    return(catdynexp)
  }
  else {
    if (distr == "apnormal" | distr == "aplnormal") {
      negsup <- ((sealen - 2)/2) * log(sum(Likel[["Likelihood"]]))
    }
    else {
      negsup <- -sum(Likel[["Likelihood"]])
    }
    return(negsup)
  }
}