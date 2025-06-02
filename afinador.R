cat_df = as.CatDynData(x=df_effort %>% 
                         filter(as.numeric(
                           as.character(year_sale)) %in% c(2006:2023)),
                       step="month",
                       fleet.name="MIS+OTB-S",
                       coleff=6,
                       colcat=5,
                       colmbw=9,
                       unitseff="trips",
                       unitscat="kg",
                       unitsmbw="kg",
                       nmult="mill",
                       season.dates=c(as.Date("2006-01-01"),
                                      last_date_of_week(2023, 52)-1))

View(cat_df$Data$`MIS+OTB-S`)

fit_null =
  trialer(cat_df,
          p = 18,
          M = 0.05976,
          N0.ini = 14020,
          P = indice_manual,
          P.ini  = list(
            5744,40315,3534,
            21558,2766,9735,
            55453,3438,16525,
            10746,14552,2524,
            22789,12567,34463,
            20000,20000,20000),
          k.ini = 0.00003546,
          alpha.ini = 1.1455,
          beta.ini  = 0.6134,
          distr = 'apnormal',
          method = 'spg',
          itnmax = 10000,
          disp = list(100))



res = fit_null$fit


annual_biomass =
  CatDynBSD(res,
            method = names(res$Model),
            multi = T,
            mbw.sd = predicos$se.fit) %>% 
  mutate(TimeStep = 216,
         x =seq(2006,2023+11/12,1/12))


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
  mutate(x = 1:216)

natural_mortality = res$Model$spg$bt.par$M
natural_mortality_sd = res$Model$spg$bt.stdev[['M']]

# Fishing Mortality
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

# Exploitation Rate
# Exploitation Rate

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


results %>% 
  summarise(pilas = `Observed.F.1/month`/(`Observed.F.1/month`+natural_mortality),
            expo = Observed.Explotrate) %>% 
  mutate(teste = pilas/expo)

