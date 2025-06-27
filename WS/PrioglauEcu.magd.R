################################################################################
#                                                                              #
#             Stock assessment of Prionace glauca (blue shark) caught          #
#                  as by-catch in Ecuadorian pelagic fisheries                 #
#                with multi-annual generalized depletion models                #
#                     January 2014 to December 2022                            #
#                    Example for U. do Algarve course                          #
#                       Ruben H. Roa-Ureta                                     #
#                                                                              #
################################################################################
#
wdir <- "/home/ruben/Workplace/Consultancy/PortugalUdoAlgarve/Course2/Day2/"
setwd(wdir)
options(max.print=12000,scipen=6)
library(CatDyn)
source("CatDynBSD9P.R")
#
################################################################################
# 2.5. Prionace glauca                                                         #
################################################################################
#
tibec.14.22.pg.1 <- read.csv("WS/tibec.14.22.pg.1.csv",header=TRUE)
#
################################################################################
# 2.5.1. Exploratory analysis for initial values                               #
################################################################################
#
tibec.14.22.pg.2 <- as.CatDynData(x=tibec.14.22.pg.1,
                                  step="month",
                                  fleet.name=c("Pelagic"),
                                  coleff=4,
                                  colcat=3,
                                  colmbw=5,
                                  unitseff="ntrips",
                                  unitscat="kg",
                                  unitsmbw="kg",
                                  nmult="thou",
                                  season.dates=c("2014-01-01","2022-12-31"))

ggplot(tibec.14.22.pg.1) +
  geom_line(aes(x = Month,
                y = Effort.ntrips))




#
plot(x=tibec.14.22.pg.2,mark=TRUE,offset=c(9,10),hem="N")
#
par(mfrow=c(2,1))
plot(tibec.14.22.pg.2$Data$Pelagic$time.step,
     tibec.14.22.pg.2$Data$Pelagic$spikecat,
     type="n")
axis(side=3,at=seq(1,108,12),lab=2014:2022)
text(tibec.14.22.pg.2$Data$Pelagic$time.step,
     tibec.14.22.pg.2$Data$Pelagic$spikecat,
     tibec.14.22.pg.2$Data$Pelagic$time.step,
     cex=0.5)
abline(v=seq(1,108,12))
plot(tibec.14.22.pg.2$Data$Pelagic$time.step,
     tibec.14.22.pg.2$Data$Pelagic$obsmbw.kg,
     type="n")
axis(side=3,at=seq(1,108,12),lab=2014:2022)
text(tibec.14.22.pg.2$Data$Pelagic$time.step,
     tibec.14.22.pg.2$Data$Pelagic$obsmbw.kg,
     tibec.14.22.pg.2$Data$Pelagic$time.step,
     cex=0.5)
abline(v=seq(1,108,12))
#
P1.Pelagic  <- 12  #2014
P2.Pelagic  <- 24  #2015
P3.Pelagic  <- 36  #2016
P4.Pelagic  <- 47  #2017 
P5.Pelagic  <- 60  #2018 
P6.Pelagic  <- 72  #2019
P7.Pelagic  <- 84  #2020 
P8.Pelagic  <- 96  #2021
P9.Pelagic  <- 107 #2022
#
tibec.14.22.pg.dates.1  <- c(head(tibec.14.22.pg.2$Data[[1]]$time.step,1),
                             P1.Pelagic,
                             P2.Pelagic,
                             P3.Pelagic,
                             P4.Pelagic,
                             P5.Pelagic,
                             P6.Pelagic,
                             P7.Pelagic,
                             P8.Pelagic,
                             P9.Pelagic,
                             tail(tibec.14.22.pg.2$Data[[1]]$time.step,1))
#
max.age   <- 25
time.step <- "month"
M.Hoenig(max.age,time.step)
#   M.pred.month M.pred.SE.month
# 1   0.01490822      0.00279294
#
x <- data.frame(key=sort(rep(1:9,12)),x=tibec.14.22.pg.2$Data$Pelagic$obscat.thou)
aggregate(x$x,
          list(x$key),
          sum)
#   Group.1        x
# 1       1 22.86496
# 2       2 19.54268
# 3       3 20.06188
# 4       4 22.80183
# 5       5 23.71732
# 6       6 14.50843
# 7       7 30.00649
# 8       8 44.04814
# 9       9 49.93188
rm(x)
#
M.ini             <- 0.015
N0.ini            <- 500
P1.ini.Pelagic    <- 100 #12  2014
P2.ini.Pelagic    <- 50  #24  2015
P3.ini.Pelagic    <- 50  #36  2016
P4.ini.Pelagic    <- 50  #47  2017
P5.ini.Pelagic    <- 100 #60  2018
P6.ini.Pelagic    <- 60  #72  2019
P7.ini.Pelagic    <- 350 #84  2020
P8.ini.Pelagic    <- 150 #96  2021
P9.ini.Pelagic    <- 150 #107 2022
k.ini.Pelagic     <- 2e-6
alpha.ini.Pelagic <- 1
beta.ini.Pelagic  <- 1
#
psi.ini.Pelagic    <- 0.25*sd(tibec.14.22.pg.2$Data$Pelagic$obscat.thou[tibec.14.22.pg.2$Data$Pelagic$obscat.thou>0])^2
psilog.ini.Pelagic <- 0.25*sd(log(tibec.14.22.pg.2$Data$Pelagic$obscat.thou[tibec.14.22.pg.2$Data$Pelagic$obscat.thou>0]))^2
#
tibec.14.22.pg.pars.ini <- log(c(M.ini,
                                 N0.ini,
                                 P1.ini.Pelagic,
                                 P2.ini.Pelagic,
                                 P3.ini.Pelagic,
                                 P4.ini.Pelagic,
                                 P5.ini.Pelagic,
                                 P6.ini.Pelagic,
                                 P7.ini.Pelagic,
                                 P8.ini.Pelagic,
                                 P9.ini.Pelagic,
                                 k.ini.Pelagic,
                                 alpha.ini.Pelagic,
                                 beta.ini.Pelagic,
                                 psi.ini.Pelagic,
                                 psilog.ini.Pelagic))
#
tibec.14.22.pg.ini.apn <- catdynexp(x=tibec.14.22.pg.2,                                 
                                    p=9,
                                    par=tibec.14.22.pg.pars.ini[1:14],
                                    dates=tibec.14.22.pg.dates.1,
                                    distr="apnormal")
plot(x=tibec.14.22.pg.ini.apn,
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.25,
     Biom.ypos=0.75,
     Cat.tstep=12,
     Cat.xpos=0.25,
     Cat.ypos=0.65)
#
################################################################################
# 2.5.2. Statistical optimization with 3 numerical methods and 4 likelihood f. #
################################################################################
#
start.time.apn.fit       <- Sys.time()
tibec.14.22.pg.apn.1.fit <- CatDynFit(x=tibec.14.22.pg.2,
                                      p=9,
                                      par=tibec.14.22.pg.pars.ini[1:14],
                                      dates=tibec.14.22.pg.dates.1,
                                      distr="apnormal",
                                      method=c("CG","spg","Nelder-Mead"),
                                      itnmax=50000)
end.time.apn.fit         <- Sys.time()
comp.time.pg.apn.1.fit   <- end.time.apn.fit - start.time.apn.fit
#
save.image("PrioglauEcu.magd.RData")
#
tibec.14.22.pg.apn.1.fit.pred.CG <- CatDynPred(x=tibec.14.22.pg.apn.1.fit,method="CG")
plot(x=tibec.14.22.pg.apn.1.fit.pred.CG,                                        ###01 
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.55,
     Biom.ypos=0.8,
     Cat.tstep=12,
     Cat.xpos=0.55,
     Cat.ypos=0.7)
#
tibec.14.22.pg.apn.1.fit.pred.spg <- CatDynPred(x=tibec.14.22.pg.apn.1.fit,method="spg")
plot(x=tibec.14.22.pg.apn.1.fit.pred.spg,                                       ###02 
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.55,
     Biom.ypos=0.8,
     Cat.tstep=12,
     Cat.xpos=0.55,
     Cat.ypos=0.7)
#
tibec.14.22.pg.apn.1.fit.pred.NM <- CatDynPred(x=tibec.14.22.pg.apn.1.fit,method="Nelder-Mead")
plot(x=tibec.14.22.pg.apn.1.fit.pred.NM,                                        ###03 
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.55,
     Biom.ypos=0.8,
     Cat.tstep=12,
     Cat.xpos=0.55,
     Cat.ypos=0.7)
#
start.time.apln.fit       <- Sys.time()
tibec.14.22.pg.apln.1.fit <- CatDynFit(x=tibec.14.22.pg.2,
                                       p=9,
                                       par=tibec.14.22.pg.pars.ini[1:14],
                                       dates=tibec.14.22.pg.dates.1,
                                       distr="aplnormal",
                                       method=c("CG","spg","Nelder-Mead"),#c("spg","Nelder-Mead"),
                                       itnmax=50000)
end.time.apln.fit         <- Sys.time()
comp.time.apln.1.fit      <- end.time.apln.fit - start.time.apln.fit
#
save.image("PrioglauEcu.magd.RData")
#
# tibec.14.22.pg.apln.1.fit.pred.CG <- CatDynPred(x=tibec.14.22.pg.apln.1.fit,method="CG")
# plot(x=tibec.14.22.pg.apln.1.fit.pred.CG,                                       
#      leg.pos="topleft",
#      Biom.tstep=1,
#      Biom.xpos=0.55,
#      Biom.ypos=0.8,
#      Cat.tstep=12,
#      Cat.xpos=0.55,
#      Cat.ypos=0.7)
#
tibec.14.22.pg.apln.1.fit.pred.spg <- CatDynPred(x=tibec.14.22.pg.apln.1.fit,method="spg")
plot(x=tibec.14.22.pg.apln.1.fit.pred.spg,                                      ###04
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.55,
     Biom.ypos=0.8,
     Cat.tstep=12,
     Cat.xpos=0.55,
     Cat.ypos=0.7)
#
tibec.14.22.pg.apln.1.fit.pred.NM <- CatDynPred(x=tibec.14.22.pg.apln.1.fit,method="Nelder-Mead")
plot(x=tibec.14.22.pg.apln.1.fit.pred.NM,                                       ###05 4813
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.55,
     Biom.ypos=0.8,
     Cat.tstep=12,
     Cat.xpos=0.55,
     Cat.ypos=0.7)
#
start.time.n.fit       <- Sys.time()
tibec.14.22.pg.n.1.fit <- CatDynFit(x=tibec.14.22.pg.2,
                                    p=9,
                                    par=tibec.14.22.pg.pars.ini[1:15],
                                    dates=tibec.14.22.pg.dates.1,
                                    distr="normal",
                                    method=c("CG","spg","Nelder-Mead"),
                                    itnmax=50000)
end.time.n.fit         <- Sys.time()
comp.time.pg.n.1.fit   <- end.time.n.fit - start.time.n.fit
#
save.image("PrioglauEcu.magd.RData")
#
tibec.14.22.pg.n.1.fit.pred.CG <- CatDynPred(x=tibec.14.22.pg.n.1.fit,method="CG")
plot(x=tibec.14.22.pg.n.1.fit.pred.CG,                                          ###06 
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.55,
     Biom.ypos=0.8,
     Cat.tstep=12,
     Cat.xpos=0.55,
     Cat.ypos=0.7)
#
tibec.14.22.pg.n.1.fit.pred.spg <- CatDynPred(x=tibec.14.22.pg.n.1.fit,method="spg")
plot(x=tibec.14.22.pg.n.1.fit.pred.spg,                                         ###07 
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.55,
     Biom.ypos=0.8,
     Cat.tstep=12,
     Cat.xpos=0.55,
     Cat.ypos=0.7)
#
tibec.14.22.pg.n.1.fit.pred.NM <- CatDynPred(x=tibec.14.22.pg.n.1.fit,method="Nelder-Mead")
plot(x=tibec.14.22.pg.n.1.fit.pred.NM,                                          ###08 
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.55,
     Biom.ypos=0.8,
     Cat.tstep=12,
     Cat.xpos=0.55,
     Cat.ypos=0.7)
#
start.time.ln.fit       <- Sys.time()
tibec.14.22.pg.ln.1.fit <- CatDynFit(x=tibec.14.22.pg.2,
                                     p=9,
                                     par=tibec.14.22.pg.pars.ini[c(1:14,16)],
                                     dates=tibec.14.22.pg.dates.1,
                                     distr="lognormal",
                                     method=c("spg","Nelder-Mead"),#c("CG","spg","Nelder-Mead"),
                                     itnmax=50000)
end.time.ln.fit         <- Sys.time()
comp.time.pg.ln.1.fit   <- end.time.ln.fit - start.time.ln.fit
#
save.image("PrioglauEcu.magd.RData")
#
# tibec.14.22.pg.ln.1.fit.pred.CG <- CatDynPred(x=tibec.14.22.pg.ln.1.fit,method="CG")
# plot(x=tibec.14.22.pg.ln.1.fit.pred.CG,                                         
#      leg.pos="topleft",
#      Biom.tstep=1,
#      Biom.xpos=0.55,
#      Biom.ypos=0.8,
#      Cat.tstep=12,
#      Cat.xpos=0.55,
#      Cat.ypos=0.7)
#
tibec.14.22.pg.ln.1.fit.pred.spg <- CatDynPred(x=tibec.14.22.pg.ln.1.fit,method="spg")
plot(x=tibec.14.22.pg.ln.1.fit.pred.spg,                                        ###09 
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.55,
     Biom.ypos=0.8,
     Cat.tstep=12,
     Cat.xpos=0.55,
     Cat.ypos=0.7)
#
tibec.14.22.pg.ln.1.fit.pred.NM <- CatDynPred(x=tibec.14.22.pg.ln.1.fit,method="Nelder-Mead")
plot(x=tibec.14.22.pg.ln.1.fit.pred.NM,                                         ###10 
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.55,
     Biom.ypos=0.8,
     Cat.tstep=12,
     Cat.xpos=0.55,
     Cat.ypos=0.7)
#
################################################################################
# 2.5.3. Examining variants                                                    #
################################################################################
#
x <- list(tibec.14.22.pg.apn.1.fit,                                             ###01
          tibec.14.22.pg.apn.1.fit,                                             ###02
          tibec.14.22.pg.apn.1.fit,                                             ###03
          tibec.14.22.pg.apln.1.fit,                                            ###04
          tibec.14.22.pg.apln.1.fit,                                            ###05
          tibec.14.22.pg.n.1.fit,                                               ###06
          tibec.14.22.pg.n.1.fit,                                               ###07
          tibec.14.22.pg.n.1.fit,                                               ###08
          tibec.14.22.pg.ln.1.fit,                                              ###09
          tibec.14.22.pg.ln.1.fit)                                              ###10
y <- c("CG",
       "spg",
       "Nelder-Mead",
       "spg",
       "Nelder-Mead",
       "CG",
       "spg",
       "Nelder-Mead",
       "spg",
       "Nelder-Mead")
#
# Gradients
tibec.14.22.pg.summary.0  <- CatDynSum(x=x,
                                       season="2014-2022",
                                       method=y)
# |gradients| < 1
tibec.14.22.pg.summary.1 <- tibec.14.22.pg.summary.0[tibec.14.22.pg.summary.0$Max.Abs.Grads.<1,]
#
# Correlations
CatDynCor(x=x[as.numeric(row.names(tibec.14.22.pg.summary.1))],
          ttl=paste(tibec.14.22.pg.summary.1$Model,
                    tibec.14.22.pg.summary.1$Distribution,
                    tibec.14.22.pg.summary.1$Method,
                    as.numeric(row.names(tibec.14.22.pg.summary.1))),
          method=y[as.numeric(row.names(tibec.14.22.pg.summary.1))],
          arr=c(3,2))
# Variant ###02 has the best good correlations, followed by variant ###01
#
# Estimates and CVs
CatDynPar(x=tibec.14.22.pg.apn.1.fit,method="CG")                               ###01
# cccccccParameter  Timing        Estimates CVpCent
# 1      M.1/month           0.014532852318    13.4
# 2        N0.thou         371.179681761009    31.9
# 3  Rec.thou.2015  2015-0 120.788282070255    42.1
# 4  Rec.thou.2016  2016-0  45.880927158144    80.3
# 5  Rec.thou.2017  2017-0  81.549894404828    47.4
# 6  Rec.thou.2017 2017-10 110.034917804256    38.3
# 7  Rec.thou.2019  2019-0  85.412910587074    40.4
# 8  Rec.thou.2020  2020-0  61.583821222408    65.8
# 9  Rec.thou.2021  2021-0 230.414084925536    34.6
# 10 Rec.thou.2022  2022-0 203.786025680308    30.2
# 11 Rec.thou.2022 2022-10  47.218045468544    68.7
# 12    k.1/ntrips           0.000002321036     1.8
# 13         alpha           0.134053988069    65.4
# 14          beta           2.134324864119     0.5
CatDynPar(x=tibec.14.22.pg.apn.1.fit,method="spg")                              ###02
#        Parameter  Timing        Estimates CVpCent
# 1      M.1/month           0.009440680017    24.2
# 2        N0.thou         252.016765093671    13.3
# 3  Rec.thou.2015  2015-0  71.948807845750    34.8
# 4  Rec.thou.2016  2016-0  30.323969224374    66.1
# 5  Rec.thou.2017  2017-0  47.034996146035    44.8
# 6  Rec.thou.2017 2017-10  66.396331021551    30.1
# 7  Rec.thou.2019  2019-0  44.422151572425    49.3
# 8  Rec.thou.2020  2020-0  37.792705942800    57.3
# 9  Rec.thou.2021  2021-0 144.548094207353    15.2
# 10 Rec.thou.2022  2022-0 129.050108293689    12.9
# 11 Rec.thou.2022 2022-10  27.171730106148    79.1
# 12    k.1/ntrips           0.000002014773     3.9
# 13         alpha           0.086789278592    49.6
# 14          beta           2.366809834336     0.4
# Both variants have very good CVs
#
# Fit to data, biomass prediction, and residual analysis
plot(x=tibec.14.22.pg.apn.1.fit.pred.CG,                                        ###01 54075
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.55,
     Biom.ypos=0.8,
     Cat.tstep=12,
     Cat.xpos=0.55,
     Cat.ypos=0.7)
plot(x=tibec.14.22.pg.apn.1.fit.pred.spg,                                       ###02 36202
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.55,
     Biom.ypos=0.8,
     Cat.tstep=12,
     Cat.xpos=0.55,
     Cat.ypos=0.7)
#
tibec.14.22.pg.summary.0$Sel.Model[2] <- 1
#
save.image("PrioglauEcu.magd.RData")
#
################################################################################
# 2.5.4. Standard error of monthly biomass                                     #
################################################################################
#
x                        <- data.frame(key=rep(1:12,9),x=tibec.14.22.pg.2$Data$Pelagic$obsmbw.kg)
y                        <- aggregate(x$x,list(x$key),sd)
tibec.14.22.pg.apn.1.bsd <- CatDynBSD9P(x=tibec.14.22.pg.apn.1.fit,
                                        method="spg",
                                        multi=TRUE,
                                        mbw.sd=y$x)
#
tibec.14.22.pg.apn.1.bsd$B.CV.pCent   <- 100*tibec.14.22.pg.apn.1.bsd$B.ton.SE/tibec.14.22.pg.apn.1.bsd$B.ton
tibec.14.22.pg.apn.1.bsd$Landings.ton <- 1e-3*tibec.14.22.pg.2$Data$Pelagic$obscat.kg
#
aggregate(tibec.14.22.pg.apn.1.bsd$B.CV.pCent,list(tibec.14.22.pg.apn.1.bsd$Month),mean) #Jan, Apr, Jun, Aug, Oct, Oct, Nov, Dec
aggregate(tibec.14.22.pg.apn.1.bsd$B.CV.pCent,list(tibec.14.22.pg.apn.1.bsd$Month),min)  #Jan, Apr, Jun, Aug, Oct, Oct, Nov, Dec
aggregate(tibec.14.22.pg.apn.1.bsd$B.CV.pCent,list(tibec.14.22.pg.apn.1.bsd$Month),max)  #Jan, Apr, Jun, Aug, Oct, Oct, Nov, Dec
#
tibec.14.22.pg.apn.1.bsd.jun <- tibec.14.22.pg.apn.1.bsd[tibec.14.22.pg.apn.1.bsd$Month=="Jun",c(1,8:9)]
#     Year    B.ton  B.ton.SE
# 6   2014 32640.47  7004.202
# 18  2015 36906.74  7926.002
# 30  2016 27294.44  7053.478
# 42  2017 33458.28  8560.491
# 54  2018 27444.64  7614.874
# 66  2019 27271.33  7986.130
# 78  2020 32217.06  9551.042
# 90  2021 41626.95 10286.433
# 102 2022 33381.62  8942.782
tibec.97.22.pg.landings <- aggregate(1e-3*tibec$Catch.kg[tibec$Species=="Prionace.glauca"],
                                     list(tibec$Year[tibec$Species=="Prionace.glauca"]),
                                     sum)
names(tibec.97.22.pg.landings) <- c("Year","Landings.ton")
#
save.image("PrioglauEcu.magd.RData")
#
################################################################################
# 2.5.5. Presentation of results, management tools                             #
################################################################################
#
# Figure 1
tiff("EcuPela.Prionace.glauca.fit.tiff",width=18,height=18,units="cm",type="cairo",res=350)
plot(x=tibec.14.22.pg.apn.1.fit.pred.spg,                                       ###02 36202
     leg.pos="topleft",
     Biom.tstep=1,
     Biom.xpos=0.55,
     Biom.ypos=0.8,
     Cat.tstep=12,
     Cat.xpos=0.55,
     Cat.ypos=0.7)
dev.off()
#
# Table 1
tibec.14.22.pg.apn.1.spg.pars <- CatDynPar(x=tibec.14.22.pg.apn.1.fit,method="spg") 
#
# Figure 2
piv  <- 1
x1   <- tibec.14.22.pg.apn.1.bsd$B.ton
x2   <- tibec.14.22.pg.apn.1.bsd$B.ton.SE
x3   <- tibec.14.22.pg.apn.1.bsd$Landings.ton
z1   <- tibec.14.22.pg.apn.1.fit$Model$spg$bt.par[3:(abs(tibec.14.22.pg.apn.1.fit$Model$spg$Type)+2)]
z2   <- tibec.14.22.pg.apn.1.fit$Model$spg$bt.stdev[3:(abs(tibec.14.22.pg.apn.1.fit$Model$spg$Type)+2)]
y1   <- tibec.14.22.pg.apn.1.fit.pred.spg$Model$Results$"Observed.F.1/month"
y2   <- tibec.14.22.pg.apn.1.fit.pred.spg$Model$Results$"Predicted.F.1/month"
t    <- tibec.14.22.pg.apn.1.fit.pred.spg$Model$Results$Period.month
M    <- as.numeric(tibec.14.22.pg.apn.1.fit$Model$spg$bt.par$M)
M.sd <- as.numeric(tibec.14.22.pg.apn.1.fit$Model$spg$bt.stdev[1])
y3   <- y1/(y1 + M) # Observed instantaneous exploitation rate F_obs/(F_obs+M)
y4   <- y2/(y2 + M) # Predicted instantaneous exploitation rate F_pred/(F_pred+M)
tiff("EcuPela.Prionace.glauca.Status.tiff",width=20,height=20,units="cm",type="cairo",res=350)
par(mfrow=c(2,2),mar=c(4,4,3,1),oma=c(1,1,1,1))
# Biomass and landings
plot(t,
     x1,
     type="l",
     xlab=paste("Month of ",fyr,"to ", lyr),
     ylab="Biomass (tonnes)",
     ylim=c(0,1.25*max(x1)),
     col="green4")
axis(side=3,at=seq(1,12*nyr,12),lab=fyr:lyr)
polygon(x=c(t[1]:tail(t,1),tail(t,1):t[1]),
        y=c(x1-piv*x2,rev(x1+piv*x2)),
        col=gray(0.15,0.45), border=NA, fillOddEven=TRUE)
lines(x=t[1]:tail(t,1),y=x3,type="h",lwd=5,col="blue")
legend("topright",lty=1,col="green4",bty="n",
       legend=expression(Monthly~population~biomass%+-%1~standard~deviations))
segments(x0=5,x1=5,y0=2000,y1=4000,col="blue",lwd=5)
text(x=7,y=1750,lab="Monthly biomass",adj=0)
# Recruitment
plot(fyr:lyr,
     z1,
     type="l",
     xlab="Year",
     ylab="Recruitment (thousand sharks)",
     ylim=c(0,1.25*max(z1)),
     col="green4")
polygon(x=c(fyr:lyr,lyr:fyr),
        y=c(z1-piv*z2,rev(z1+piv*z2)),
        col=gray(0.15,0.45), border=NA, fillOddEven=TRUE)
legend("topright",lty=1,col="green4",bty="n",
       legend=expression(Annual~recruitment%+-%1~standard~deviations))
# Fishing mortality and natural mortality
plot(t,
     y1,
     type="l",
     lty=1,
     xlab=paste("Month of ",fyr,"to ", lyr),
     ylab="Fishing Mortality (1/month)",
     ylim=c(0,0.02))
lines(t,y2,lty=2)
legend("topright",lty=c(1,2),legend=c("Observed","Predicted"),bty="n")
polygon(x=c(head(t,1),tail(t,1),tail(t,1),head(t,1)),
        y=c(M + piv*M.sd, M + piv*M.sd, M - piv*M.sd, M - piv*M.sd),
        col=gray(0.21,0.45),
        border=NA)
lines(x=c(head(t,1),tail(t,1)),
      y=c(M,M),col="green4")
legend("bottomleft",lty=1,col="green4",bty="n",
       legend=expression(Natural~mortality~rate%+-%1~standard~deviations))
axis(side=3,at=seq(1,12*nyr,12),lab=fyr:lyr)
# Instantaneous exploitation rate
plot(t,
     y3,
     type="l",
     ylim=c(0,1),
     xlab=paste("Month of ",fyr,"to ", lyr),
     ylab="Instantaneous exploitation rate")
lines(t,y4,type="l",ylim=c(0,1),lty=2)
abline(h=0.4,lwd=2)
legend("topright",lty=c(1,2),legend=c("Observed","Predicted"),bty="n")
axis(side=3,at=seq(1,12*nyr,12),lab=fyr:lyr)
dev.off()
#
save.image("PrioglauEcu.magd.RData")
#
################################################################################
# 2.5.6. Pella-Tomlinson surplus production model                              #
################################################################################
#
# Toy model for initial values
fyr.landings    <- 1997
B0.ini          <- 40000
K.ini           <- 40000
p.ini           <- 2
r.ini           <- 0.5
Stock.pg.pt     <- rep(0,lyr-fyr.landings+1)
Stock.pg.pt[1]  <- B0.ini
for(i in 2:(lyr-fyr.landings+1))
{
  Stock.pg.pt[i] <- Stock.pg.pt[i-1] + (r.ini/p.ini)*Stock.pg.pt[i-1]*(1-(Stock.pg.pt[i-1]/K.ini)^(p.ini-1))-tibec.97.22.pg.landings$Landings.ton[i-1];
} 
#
par(mfrow=c(1,1),mai=c(1,1,1,1))
plot(fyr.landings:lyr,Stock.pg.pt,type="l",ylim=c(0,50000))
points(fyr.landings:lyr,c(rep(NA,fyr-fyr.landings),tibec.14.22.pg.apn.1.bsd.jun$B.ton))
lines(fyr.landings:lyr,tibec.97.22.pg.landings$Landings.ton,type="h",lwd=20,col="blue")
#
(log.B0.ini <- log(B0.ini))
# [1] 10.59663
(log.K.ini  <- log(K.ini))
# [1] 10.59663
(log.p.ini  <- log(p.ini))
# [1] 0.6931472
(log.r.ini  <- log(r.ini))
# [1] -0.6931472
#
save.image("PrioglauEcu.magd.RData")
#
tibec.97.22.pg.landings$Landings.ton
#35.73069 12.46824 15.07896 41.43262 50.46805 78.43099 294.69662 52.00875 110.45789 489.93647 646.00733 832.16444 439.48788 1918.56360 1698.71698 1728.78014 3366.09372 2674.07745 2396.98280 2078.99992 2360.09016 2496.00873 1358.10350 2955.77170 4309.36228 4813.19810
tibec.14.22.pg.apn.1.bsd.jun
#     Year    B.ton  B.ton.SE
# 6   2014 32640.47  7004.202
# 18  2015 36906.74  7926.002
# 30  2016 27294.44  7053.478
# 42  2017 33458.28  8560.491
# 54  2018 27444.64  7614.874
# 66  2019 27271.33  7986.130
# 78  2020 32217.06  9551.042
# 90  2021 41626.95 10286.433
# 102 2022 33381.62  8942.782
#
# Collecting results from ADMB
# 4p variant did not converge succesfully
maxgrad.admb.pt.pgec.3p        <- read.table("marlikptpgec3p.par",comment.char="",nrows=1)$V16
#[1] 0.00009773193
lik.admb.pt.pgec.3p            <- read.table("marlikptpgec3p.par",comment.char="",nrows=1)$V11
npar.admb.pt.pgec.3p           <- read.table("marlikptpgec3p.par",comment.char="",nrows=1)$V6
AIC.admb.pt.pgec.3p            <- 2*npar.admb.pt.pgec.3p - 2*lik.admb.pt.pgec.3p
#[1] -175.8333
#
pars.admb.pt.pgec.3p           <- read.table("marlikptpgec3p.std",skip=4,header=FALSE)[1:3,3:4]
names(pars.admb.pt.pgec.3p)    <- c("MLE","SD")
rownames(pars.admb.pt.pgec.3p) <- c("K","p","r")
K3p.pgec                       <- round(pars.admb.pt.pgec.3p["K","MLE"])
p3p.pgec                       <- pars.admb.pt.pgec.3p["p","MLE"]
r3p.pgec                       <- pars.admb.pt.pgec.3p["r","MLE"]
Stock.admb.pt.pgec.3p          <- cbind(fyr.landings:(lyr+1),read.table("marlikptpgec3p.std",skip=7,header=FALSE)[,3:4])
names(Stock.admb.pt.pgec.3p)   <- c("Year","Biomass.ton","SD.Biomass.ton")
Stock.admb.pt.pgec.3p$B.pcCV   <- 100*Stock.admb.pt.pgec.3p$SD.Biomass.ton/Stock.admb.pt.pgec.3p$Biomass.ton
#
# Derived parameters and their SE
pt.cor.pgec.3p                   <- read.table("marlikptpgec3p.cor",skip=5,header=FALSE,fill=TRUE,)[1:3,8:10]
pt.cor.pgec.3p[1,2]              <- pt.cor.pgec.3p[2,1]
pt.cor.pgec.3p[1,3]              <- pt.cor.pgec.3p[3,1]
pt.cor.pgec.3p[2,3]              <- pt.cor.pgec.3p[3,2]
row.names(pt.cor.pgec.3p)        <- c("K","p","r")
names(pt.cor.pgec.3p)            <- c("K","p","r")
pt.sd.pgec.3p                    <- read.table("marlikptpgec3p.cor",skip=5,header=FALSE,fill=TRUE)[1:3,4]
pt.cov.pgec.3p                   <- cor2cov(cor.mat=as.matrix(pt.cor.pgec.3p),sd=as.vector(pt.sd.pgec.3p))
#
MSY.pgec.3p          <- round(K3p.pgec*r3p.pgec*(p3p.pgec-1)/(p3p.pgec^(p3p.pgec/(p3p.pgec-1))))
#[1] 2035
MSY.SE.pgec.3p       <- round(deltamethod(g=list(~x1*x3*(x2-1)/(x2^(x2/(x2-1)))),
                                          mean=c(K3p.pgec,p3p.pgec,r3p.pgec),
                                          cov=pt.cov.pgec.3p,
                                          ses=TRUE))
#[1] 7348
Mean.Catch.pgec.t    <- mean(tibec.97.22.pg.landings$Landings.ton)
#[1] 1432.812
#
B.MSY.pgec.3p        <- round(K3p.pgec*(p3p.pgec^(1/(1-p3p.pgec))))
#[1] 16659
B.MSY.SE.pgec.3p     <- round(deltamethod(g=list(~x1*x2^(1/(1-x2))),
                                          mean=c(K3p.pgec,p3p.pgec),
                                          cov=pt.cov.pgec.3p[c(1,2),c(1,2)],
                                          ses=TRUE))
#[1] 8914
#
F.MSY.pgec.3p        <- MSY.SE.pgec.3p/B.MSY.pgec.3p
#[1] 0.1832359
F.MSY.SE.pgec.3p     <- deltamethod(g=list(~x1/x2),
                                    mean=c(MSY.SE.pgec.3p,B.MSY.pgec.3p),
                                    cov=matrix(c(MSY.SE.pgec.3p^2,0,0,B.MSY.SE.pgec.3p^2),2,2),
                                    ses=TRUE)
#[1] 0.1874385
#
Mean.Catch.pgec.t <- mean(tibec.97.22.pg.landings$Landings.ton)
#[1] 103.0185
SD.Catch.pgec.t   <- sd(tibec.97.22.pg.landings$Landings.ton)
#[1] 271.626
#
# Biomass dynamics plot
pgec.popdyn.3p <- data.frame(Year=fyr.landings:(lyr+1),
                             Catch.ton=c(tibec.97.22.pg.landings$Landings.ton,NA),
                             Catdyn.Biom.ton=c(rep(NA,fyr-fyr.landings),tibec.14.22.pg.apn.1.bsd.jun$B.ton,NA),
                             Catdyn.Biom.SE.ton=c(rep(NA,fyr-fyr.landings),tibec.14.22.pg.apn.1.bsd.jun$B.ton.SE,NA),
                             Catdyn.Biom.pcCV=NA,
                             PT.Biom.ton=Stock.admb.pt.pgec.3p$Biomass.ton,
                             PT.Biom.SE.ton=Stock.admb.pt.pgec.3p$SD.Biomass.ton,
                             PT.Biom.pcCV=Stock.admb.pt.pgec.3p$B.pcCV)
#
pgec.popdyn.3p$Catdyn.Biom.pcCV <- 100*pgec.popdyn.3p$Catdyn.Biom.SE.ton/pgec.popdyn.3p$Catdyn.Biom.ton
#
tiff("EcuPela.Prionace.glauca.BiomDyn.tiff",width=25,height=25,units="cm",type="cairo",res=350)
par(mfrow=c(1,1),mar=c(4,4.5,1,1),oma=c(4,4,1,1),bty="o")
plot(x=fyr.landings:(lyr+1),
     y=pgec.popdyn.3p$Catdyn.Biom.ton,
     pch=19,
     ylim=c(0,1.25*max(pgec.popdyn.3p$Catdyn.Biom.ton,na.rm=TRUE)),
     cex=2.0/pgec.popdyn.3p$Catdyn.Biom.pcCV/min(1/pgec.popdyn.3p$Catdyn.Biom.pcCV,na.rm=TRUE),
     cex.axis=1.25,
     cex.lab=1.25,
     col="green4",
     type="b",
     lwd=2,
     ylab="",
     xlab="",
     main=expression(paste(italic(Prionace~glauca),"- Ecuadorian pelagic fisheries"),sep=""))
polygon(x=c(fyr.landings:(lyr+1),(lyr+1):fyr.landings),
        y=c(pgec.popdyn.3p$PT.Biom.ton-piv*pgec.popdyn.3p$PT.Biom.SE.ton,
            rev(pgec.popdyn.3p$PT.Biom.ton+piv*pgec.popdyn.3p$PT.Biom.SE.ton)),
        col=gray(0.15,0.45), 
        border=NA, 
        fillOddEven=TRUE)
lines(fyr.landings:(lyr+1),
      pgec.popdyn.3p$PT.Biom.ton,
      lwd=2,
      type="b",
      pch=19)
lines(fyr.landings:(lyr+1),
      pgec.popdyn.3p$Catch.ton,
      lwd=2,
      type="b",
      pch=19,
      col="blue")
legend("left",
       legend=c("K=38,375 (8,284)","p=1.421 (1.621)","r=0.412 (0.571)"))
legend("topright",
       pch=c(19,19,19),
       col=c("green4","black","blue"),
       lty=c(1,1,1),
       lwd=c(1,2,2),
       legend=c("Depletion model biomass (Size of dot proportional to precision)",
                expression(Surplus~production~biomass%+-%1~standard~errors),
                "Annual landings"))
#
mtext(side=2,outer=TRUE,text="Biomass (tonnes)",cex=2)
mtext(side=1,outer=TRUE,text="Year",cex=2)
dev.off()
#
save.image("PrioglauEcu.magd.RData")
#
################################################################################
# 2.5.7. Phase diagram                                                         #
################################################################################
#
pgec.popdyn.3p$B.BMSY      <- NA
pgec.popdyn.3p$F.FMSY      <- NA
pgec.popdyn.3p$B.BMSY.SE   <- NA
pgec.popdyn.3p$F.FMSY.SE   <- NA
pgec.popdyn.3p$MSY.Status  <- NA
#
pgec.popdyn.3p$B.BMSY  <- (pgec.popdyn.3p$PT.Biom.ton)/(B.MSY.pgec.3p)
pgec.popdyn.3p$F.FMSY  <- (pgec.popdyn.3p$Catch.ton/pgec.popdyn.3p$PT.Biom.ton)/(MSY.pgec.3p/B.MSY.pgec.3p)
# Standard errors of B/B_MSY and F/F_MSY
for(i in 1:dim(pgec.popdyn.3p)[1])
{
  pgec.popdyn.3p$B.BMSY.SE[i] <- deltamethod(g=list(~x1/x2),
                                             mean=c(pgec.popdyn.3p$PT.Biom.ton[i],B.MSY.pgec.3p),
                                             cov=matrix(c(pgec.popdyn.3p$PT.Biom.SE.ton[i]^2,0,0,B.MSY.SE.pgec.3p^2),2,2),
                                             ses=TRUE)
}
for(i in 1:dim(pgec.popdyn.3p)[1])
{
  yield <- pgec.popdyn.3p$Catch.ton[i]
  pgec.popdyn.3p$F.FMSY.SE[i] <- deltamethod(g=list(~(yield/x1)/(x2/x3)),
                                             mean=c(pgec.popdyn.3p$PT.Biom.ton[i],MSY.pgec.3p,B.MSY.pgec.3p),
                                             cov=matrix(c(pgec.popdyn.3p$PT.Biom.SE.ton[i]^2,0,0,0,MSY.SE.pgec.3p^2,0,0,0,B.MSY.SE.pgec.3p^2),3,3),
                                             ses=TRUE)
  
}
# Stock status for state counts barplot                                                    
pgec.popdyn.3p$MSY.Status <- ifelse(pgec.popdyn.3p$B.BMSY > 1 & pgec.popdyn.3p$F.FMSY < 1,1,
                                    ifelse(pgec.popdyn.3p$B.BMSY < 1 & pgec.popdyn.3p$F.FMSY < 1,2,
                                           ifelse(pgec.popdyn.3p$B.BMSY > 1 & pgec.popdyn.3p$F.FMSY > 1,3,4)))
# Kobe plot
tiff("EcuPela.Prionace.glauca.Kobe.tiff",width=22,height=20,units="cm",res=350,type="cairo")
par(mfrow=c(1,1),mar=c(4,5,2,1))
xlim  <- c(0,2.5)
ylim  <- c(0,2.5)
plot(pgec.popdyn.3p$B.BMSY[1:(dim(pgec.popdyn.3p)[1]-1)],
     pgec.popdyn.3p$F.FMSY[1:(dim(pgec.popdyn.3p)[1]-1)],
     type="b",
     xlim=xlim,
     ylim=ylim,
     xlab=expression(B/B[MSY]),
     ylab=expression(F/F[MSY]),
     main=expression(paste(italic(Prionace~glauca)," - Ecuadorian pelagic fisheries")),
     cex.axis=2,
     cex.lab=2,
     cex.main=2)
abline(h=1)
abline(v=1)
rect(xleft=-1,xright=1,ybottom=-1,ytop=1,col="yellow")
rect(xleft=-1,xright=1,ybottom=1,ytop=ylim[2]+1,col="red")
rect(xleft=1,xright=xlim[2]+1,ybottom=1,ytop=ylim[2]+1,col="orange1")
rect(xleft=1,xright=xlim[2]+1,ybottom=-1,ytop=1,col="green")
points(pgec.popdyn.3p$B.BMSY,
       pgec.popdyn.3p$F.FMSY,
       pch=19,
       cex=3,
       col="blue")
lines(pgec.popdyn.3p$B.BMSY,
      pgec.popdyn.3p$F.FMSY,
      col="blue",
      lwd=2)
text(x=0,y=ylim[2],lab="Overfished",cex=1.5,adj=0)
text(x=xlim[2],y=ylim[2],lab="Overfishing",cex=1.5,adj=1)
text(x=0,y=0,lab="Rebuilding",cex=1.5,adj=0)
text(x=xlim[2],y=0,lab="Sustainable",cex=1.5,adj=1)
segments(x0=pgec.popdyn.3p$B.BMSY[dim(pgec.popdyn.3p)[1]-1]-piv*pgec.popdyn.3p$B.BMSY.SE[dim(pgec.popdyn.3p)[1]-1],
         x1=pgec.popdyn.3p$B.BMSY[dim(pgec.popdyn.3p)[1]-1]+piv*pgec.popdyn.3p$B.BMSY.SE[dim(pgec.popdyn.3p)[1]-1],
         y0=pgec.popdyn.3p$F.FMSY[dim(pgec.popdyn.3p)[1]-1],
         y1=pgec.popdyn.3p$F.FMSY[dim(pgec.popdyn.3p)[1]-1],
         lwd=4,
         col="black")
segments(x0=pgec.popdyn.3p$B.BMSY[dim(pgec.popdyn.3p)[1]-1]+piv*pgec.popdyn.3p$B.BMSY.SE[dim(pgec.popdyn.3p)[1]-1],
         x1=pgec.popdyn.3p$B.BMSY[dim(pgec.popdyn.3p)[1]-1]+piv*pgec.popdyn.3p$B.BMSY.SE[dim(pgec.popdyn.3p)[1]-1],
         y0=pgec.popdyn.3p$F.FMSY[dim(pgec.popdyn.3p)[1]-1]-0.05,
         y1=pgec.popdyn.3p$F.FMSY[dim(pgec.popdyn.3p)[1]-1]+0.05,
         lwd=4,
         col="black")
segments(x0=pgec.popdyn.3p$B.BMSY[dim(pgec.popdyn.3p)[1]-1]-piv*pgec.popdyn.3p$B.BMSY.SE[dim(pgec.popdyn.3p)[1]-1],
         x1=pgec.popdyn.3p$B.BMSY[dim(pgec.popdyn.3p)[1]-1]-piv*pgec.popdyn.3p$B.BMSY.SE[dim(pgec.popdyn.3p)[1]-1],
         y0=pgec.popdyn.3p$F.FMSY[dim(pgec.popdyn.3p)[1]-1]-0.05,
         y1=pgec.popdyn.3p$F.FMSY[dim(pgec.popdyn.3p)[1]-1]+0.05,
         lwd=4,
         col="black")
segments(x0=pgec.popdyn.3p$B.BMSY[dim(pgec.popdyn.3p)[1]-1],
         x1=pgec.popdyn.3p$B.BMSY[dim(pgec.popdyn.3p)[1]-1],
         y0=pgec.popdyn.3p$F.FMSY[dim(pgec.popdyn.3p)[1]-1]-piv*pgec.popdyn.3p$F.FMSY.SE[dim(pgec.popdyn.3p)[1]-1],
         y1=pgec.popdyn.3p$F.FMSY[dim(pgec.popdyn.3p)[1]-1]+piv*pgec.popdyn.3p$F.FMSY.SE[dim(pgec.popdyn.3p)[1]-1],
         lwd=4,
         col="black")
segments(x0=pgec.popdyn.3p$B.BMSY[dim(pgec.popdyn.3p)[1]-1]-0.05,
         x1=pgec.popdyn.3p$B.BMSY[dim(pgec.popdyn.3p)[1]-1]+0.05,
         y0=pgec.popdyn.3p$F.FMSY[dim(pgec.popdyn.3p)[1]-1]-piv*pgec.popdyn.3p$F.FMSY.SE[dim(pgec.popdyn.3p)[1]-1],
         y1=pgec.popdyn.3p$F.FMSY[dim(pgec.popdyn.3p)[1]-1]-piv*pgec.popdyn.3p$F.FMSY.SE[dim(pgec.popdyn.3p)[1]-1],
         lwd=4,
         col="black")
segments(x0=pgec.popdyn.3p$B.BMSY[dim(pgec.popdyn.3p)[1]-1]-0.05,
         x1=pgec.popdyn.3p$B.BMSY[dim(pgec.popdyn.3p)[1]-1]+0.05,
         y0=pgec.popdyn.3p$F.FMSY[dim(pgec.popdyn.3p)[1]-1]+piv*pgec.popdyn.3p$F.FMSY.SE[dim(pgec.popdyn.3p)[1]-1],
         y1=pgec.popdyn.3p$F.FMSY[dim(pgec.popdyn.3p)[1]-1]+piv*pgec.popdyn.3p$F.FMSY.SE[dim(pgec.popdyn.3p)[1]-1],
         lwd=4,
         col="black")
text(pgec.popdyn.3p$B.BMSY[1:(dim(pgec.popdyn.3p)[1]-1)],
     pgec.popdyn.3p$F.FMSY[1:(dim(pgec.popdyn.3p)[1]-1)],
     lab=substring(pgec.popdyn.3p$Year[1:(dim(pgec.popdyn.3p)[1]-1)],3,4),
     col="yellow")
par(new=TRUE,mai=c(4.5,4.5,2,2))
barplot(c(sum(pgec.popdyn.3p$MSY.Status==1,na.rm=TRUE),
          sum(pgec.popdyn.3p$MSY.Status==2,na.rm=TRUE),
          sum(pgec.popdyn.3p$MSY.Status==3,na.rm=TRUE),
          sum(pgec.popdyn.3p$MSY.Status==4,na.rm=TRUE)),
        width=c(0.5,0.5,0.5,0.5),
        col=c("green","yellow","orange1","red"),
        names.arg=c(paste(round(100*sum(pgec.popdyn.3p$MSY.Status==1,na.rm=TRUE)/sum(!is.na(pgec.popdyn.3p$MSY.Status))),"%",sep=""),
                    paste(round(100*sum(pgec.popdyn.3p$MSY.Status==2,na.rm=TRUE)/sum(!is.na(pgec.popdyn.3p$MSY.Status))),"%",sep=""),
                    paste(round(100*sum(pgec.popdyn.3p$MSY.Status==3,na.rm=TRUE)/sum(!is.na(pgec.popdyn.3p$MSY.Status))),"%",sep=""),
                    paste(round(100*sum(pgec.popdyn.3p$MSY.Status==4,na.rm=TRUE)/sum(!is.na(pgec.popdyn.3p$MSY.Status))),"%",sep="")),
        cex.names=1.0)
dev.off()
#
save.image("PrioglauEcu.magd.RData")
#
##########################   END OF SCRIPT    ##################################