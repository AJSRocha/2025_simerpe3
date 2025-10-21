
cat_df = as.CatDynData(x=df_effort %>%
                         filter(as.numeric(
                           as.character(year_sale)) %in% c(2006:2020)),
                       step="month",
                       fleet.name="MIS+OTB-S",
                       coleff=6,
                       colcat=5,
                       colmbw=9,
                       unitseff="trips",
                       unitscat="kg",
                       unitsmbw="kg",
                       nmult="thou",
                       season.dates=c(as.Date("2006-01-01"),
                                      last_date_of_week(2020, 52)-1))

# par = list(
#   logRt_scaled = log(c(10000,10000,10000,
#                        10000,10000,10000,
#                        10000,10000,10000,
#                        10000,10000,10000,
#                        10000,10000,10000)*1.5),
#   logalpha      = log(1.2),
#   logbeta       = log(0.6),
#   logK          = log(.000035),
#   logN0_scaled  = log(10000),
#   logM          = log(0.05),
#   logsdCt       = log(0.25 * sd(df_effort$catch[df_effort$catch>0]))  # consistent with CatDyn
# )


# Hipotese 1

par = list(
  logRt_scaled = log(c(5700,40000,3500,
                       21000,2700,9700,
                       55000,3400,16500,
                       10700,14500,2500,
                       22700,12500,34400)),
  logalpha      = log(1.14),
  logbeta       = log(0.6),
  logK          = log(.000034),
  logN0_scaled  = log(28000),
  logM          = log(0.1),
  logsdCt       = log(0.25 * sd(df_effort$catch[df_effort$catch>0]))  # consistent with CatDyn
)


# Hipotese 2

par = list(
  logRt_scaled = log(c(5700,40000,3500,
                       21000,2700,9700,
                       55000,3400,16500,
                       10700,14500,2500,
                       22700,12500,34400)),
  logalpha      = log(1.14),
  logbeta       = log(0.6),
  logK          = log(.000034),
  logN0_scaled  = log(25000),
  logM          = log(0.1),
  logsdCt       = log(0.25 * sd(df_effort$catch[df_effort$catch>0]))  # consistent with CatDyn
)


indice_manual =
  list(
    10,11,12,
    12,12,10,11,12,
    11,10,12,8,10,
    12,11)

for(i in 0:(length(indice_manual)-1)){
  indice_manual[[i+1]] = 12*i + indice_manual[[i+1]]
}

fit_null =
  trialer(cat_df,
          p = 15,
          M = exp(par$logM),
          N0.ini = exp(par$logN0_scaled),
          P = indice_manual,
          P.ini  = as.list(exp(par$logRt_scaled)),
          k.ini = exp(par$logK),
          alpha.ini = exp(par$logalpha),
          beta.ini  = exp(par$logbeta),
          distr = distribuicoes[5],
          method = optimizadores[2],
          itnmax = 10000,
          disp = list(100))

plotador(cat_df, fit_null, pre =F, post1 = T, post2 = T)

fit_null$fit$Model$spg
catdynpred <- CatDynPred(x=fit_null$fit,method="spg")
cat_pars_1 = fit_null$fit$Model$spg$bt.par
cat_pars_1

annual_biomass =
  CatDynBSD(fit_null$fit,
            method = names(fit_null$fit$Model),
            multi = T,
            mbw.sd = predicos$se.fit) %>% 
  mutate(TimeStep = 180,
         x =seq(2006,2020+11/12,1/12)) 

catch_cat = unlist(CatDynPred(x=fit_null$fit,method="spg")$Model$Results['Predicted.Catch.thou'])


ggplot() + 
  geom_line(aes(x = annual_biomass$x,
                y = annual_biomass$B.ton,
                group = 1),
            size = 1) +
  geom_ribbon(aes(x = annual_biomass$x,
                  y = annual_biomass$B.ton,
                  ymin= annual_biomass$B.ton- 2*annual_biomass$B.ton.SE,
                  ymax= annual_biomass$B.ton+ 2*annual_biomass$B.ton.SE),
              alpha=0.2) +
  
  # coord_cartesian(ylim = c(0, 300000), xlim = c(1995,2024)) +
  theme_bw()

results = fit_null$pred$Model$Results %>% 
  as.data.frame() %>% 
  mutate(x = 1:180)

natural_mortality = fit_null$fit$Model$spg$bt.par$M
natural_mortality_sd = fit_null$fit$Model$spg$bt.stdev[['M']]

results %>% 
  ggplot() + 
  geom_line(aes(x = Period.month,
                y = `Observed.F.1/month`),
            col = 'tomato',
            size = 1) +
  geom_line(aes(x = Period.month,
                y = `Predicted.F.1/month`),
            col = 'darkred',
            size = 1,
            linetype = 2) + 
  geom_hline(yintercept = natural_mortality,
             col = 'darkgreen',
             size = 1,
             linetype = 1) +
  geom_hline(yintercept = natural_mortality + 2*natural_mortality,
             col = 'darkgreen',
             size = 1,
             linetype = 2) +
  geom_hline(yintercept = natural_mortality - 2*natural_mortality,
             col = 'darkgreen',
             size = 1,
             linetype = 2) +
  theme_bw()


results %>% 
  ggplot() + 
  geom_line(aes(x = Period.month,
                y = `Observed.Explotrate`),
            col = 'tomato',
            size = 1) +
  
  geom_line(aes(x = Period.month,
                y = `Predicted.Explotrate`),
            col = 'darkred',
            size = 1,
            linetype = 2) + 
  geom_hline(yintercept = 0.4) + 
  # geom_line(aes(x = Period.month,
  #               y = `Observed.F.1/month`/(`Observed.F.1/month`+ natural_mortality)),
  #           col = 'purple',
  #           size = 1) +
  # geom_hline(yintercept = 0.4,
  #            col = 'darkgreen',
  #            size = 1,
  #            linetype = 1) +
  theme_bw()

