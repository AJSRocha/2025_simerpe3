#getAnywhere(".CDMN9P") - sem transito
#getAnywhere(".CDMNT9P") - com transito

function (par, dates, obscat1, obseff1, obsmbm1, distr, properties, 
          output, pp = -1, partial) 
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
  sealen <- ts.end - ts.start + 1
  nstep <- vector("numeric", sealen)
  mccum <- vector("numeric", sealen)
  effeff1 <- vector("numeric", sealen)
  effn1 <- vector("numeric", sealen)
  predcat1 <- vector("numeric", sealen)
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

CDMN1P1P
function (par, dates, obscat1, obseff1, obsmbm1, obscat2, obseff2, 
          obsmbm2, distr, properties, output, pp = c(1, 1), partial) 
{
  ts.start <- head(dates, 1)
  ts1.P1 <- dates[2]
  ts2.P1 <- dates[3]
  ts.end <- tail(dates, 1)
  sealen <- ts.end - ts.start + 1
  nstep <- vector("numeric", sealen)
  mccum <- vector("numeric", sealen)
  effeff1 <- vector("numeric", sealen)
  effn1 <- vector("numeric", sealen)
  predcat1 <- vector("numeric", sealen)
  ind1.P1 <- ifelse(1:sealen < (ts1.P1 - ts.start + 1), 0, 
                    1)
  effeff2 <- vector("numeric", sealen)
  effn2 <- vector("numeric", sealen)
  predcat2 <- vector("numeric", sealen)
  ind2.P1 <- ifelse(1:sealen < (ts2.P1 - ts.start + 1), 0, 
                    1)
  logM <- par[1]
  logN0 <- par[2]
  logP1F1 <- par[3]
  logscale1 <- par[4]
  logalpha1 <- par[5]
  logbeta1 <- par[6]
  logP1F2 <- par[7]
  logscale2 <- par[8]
  logalpha2 <- par[9]
  logbeta2 <- par[10]
  mccum[1] <- 0
  nstep[1] <- exp(logN0) * exp(-exp(logM))
  for (i in 2:sealen) {
    mccum[i] <- obscat1[i - 1] + obscat2[i - 1] + mccum[i - 
                                                          1] * exp(-exp(logM))
    nstep[i] <- exp(logN0) * exp(-exp(logM) * i) + ind1.P1[i] * 
      exp(logP1F1) * exp(-exp(logM) * (i - (ts1.P1 - ts.start + 
                                              1))) + ind2.P1[i] * exp(logP1F2) * exp(-exp(logM) * 
                                                                                       (i - (ts2.P1 - ts.start + 1))) - mccum[i] * exp(-exp(logM)/2)
  }
  effeff1 <- obseff1^(exp(logalpha1))
  effn1 <- nstep^(exp(logbeta1))
  predcat1 <- exp(logscale1) * (effeff1 * effn1) * exp(-exp(logM)/2)
  effeff2 <- obseff2^(exp(logalpha2))
  effn2 <- nstep^(exp(logbeta2))
  predcat2 <- exp(logscale2) * (effeff2 * effn2) * exp(-exp(logM)/2)
  Likel <- .CatDynLik2F(obscat1, predcat1, obscat2, predcat2, 
                        distr, par)
  if (output == "predict") {
    catdynexp <- .CatDynExp2F.Res(properties, nstep, obsmbm1, 
                                  obsmbm2, pp, dates, distr, par, Likel, ts.start, 
                                  ts.end, obseff1, obscat1, predcat1, obseff2, obscat2, 
                                  predcat2)
    class(catdynexp) <- "CatDynExp"
    return(catdynexp)
  }
  else {
    if (distr[1] == "apnormal" | distr[1] == "aplnormal" & 
        distr[2] != "apnormal" & distr[2] != "aplnormal") {
      negsup <- ((sealen - 2)/2) * log(sum(Likel[["Likelihood"]][, 
                                                                 1])) - sum(Likel[["Likelihood"]][, 2])
    }
    else if (distr[2] == "apnormal" | distr[2] == "aplnormal" & 
             distr[1] != "apnormal" & distr[1] != "aplnormal") {
      negsup <- ((sealen - 2)/2) * log(sum(Likel[["Likelihood"]][, 
                                                                 2])) - sum(Likel[["Likelihood"]][, 1])
    }
    else if (distr[1] == "apnormal" | distr[1] == "aplnormal" & 
             distr[2] == "apnormal" | distr[2] == "aplnormal") {
      negsup <- ((sealen - 2)/2) * (log(sum(Likel[["Likelihood"]][, 
                                                                  1])) + log(sum(Likel[["Likelihood"]][, 2])))
    }
    else {
      negsup <- -sum(Likel[["Likelihood"]][, 1]) - sum(Likel[["Likelihood"]][, 
                                                                             2])
    }
    return(negsup)
  }
}