# funcoes para fittar catdyn

trialer = function(data, p, M, N0.ini, P.ini, k.ini,
                   alpha.ini, beta.ini, P,
                   distr, method, itnmax, disp = list()){
  
  # psi.ini   = 0.33*sd(log(data$Data$`Polyvalent-S`$obscat.thou))^2
  
  
  if(p>0){
    pars.ini = log(c(M,
                     N0.ini,
                     unlist(P.ini), # estimativa de amplitude da perturbacao
                     k.ini,
                     alpha.ini,
                     beta.ini))
    
    dates = c(head(data$Data[[1]]$time.step,1),
              unlist(P), #estimativa do timing da perturbacao
              tail(data$Data[[1]]$time.step,1))}
  else{
    pars.ini = log(c(M,
                     N0.ini,
                     # unlist(P.ini), # estimativa de amplitude da perturbacao
                     k.ini,
                     alpha.ini,
                     beta.ini))
    
    
    # negative binomial, normal, lognormal, gamma, robust lognormal or gumbel distribution  
    
    
    
    dates = c(head(data$Data[[1]]$time.step,1),
              # unlist(P), #estimativa do timing da perturbacao
              tail(data$Data[[1]]$time.step,1))
  }
  
  if(distr %in% c('negbin','normal','lognormal','gamma')){
    pars.ini = c(pars.ini, log(unlist(disp)))
  }  
  
  res = list()
  
  res$pre_fit = catdynexp(x=data,
                          p=p,
                          par=pars.ini,
                          dates=dates,
                          distr=distr)
  
  res$fit = CatDynFit_2(x = data,
                      p = p,
                      par = pars.ini,
                      dates = dates,
                      distr = distr,
                      method = method,
                      itnmax = itnmax)
  
  res$pred = CatDynPred(res$fit,method)
  
  return(res)
}

plotador = function(data, model, pre = T, post1 = T, post2 = T){
  if(pre){
    plot.CatDynData(data,
                    mark = T,
                    offset = c(0,1,10),
                    hem = 'N')}
  if(post1){
    plot(x=model$pre_fit,
         leg.pos="topright",
         Biom.tstep=7,
         Cat.tstep=120,
         Biom.xpos=0.4,
         Biom.ypos=0,
         Cat.xpos=0.4,
         Cat.ypos=0.1)}
  
  if(post2){
    plot(x=model$pred,
         leg.pos="topright",
         Biom.tstep=7,
         Cat.tstep=10,
         Biom.xpos=0.18,
         Biom.ypos=0.1,
         Cat.xpos=0.18,
         Cat.ypos=0.2)}
  
}

# funcao que tapa buracos de defeso
defeso = function(df){
  semanas =c(1:52)[!(c(1:52) %in% unique(df$week))]
  ref = df_effort %>% filter(year_sale == 2000)
  
  for(i in semanas){
    linha = data.frame("year_sale" = df$year_sale[1],
                       "week" = i,
                       "catch" = 0,
                       "effort" = 0,
                       "res" = NA,
                       "res.se" = NA,
                       "mbw" = ref[ref$week == i,]$mbw,
                       "se" = ref[ref$week == i,]$se,
                       "mbw_rand" = ref[ref$week == i,]$mbw_rand)
    
    df = rbind(data.frame(df), linha) %>% 
      arrange(week)
  }
  return(df)}

#Funcao que determina ultimo dia da semana
last_date_of_week = function(year, week){strptime(paste(year, week, 1), format = "%Y %W %u")}

# funcoes para simular catdyn
catch_module = function(index,M){ #index é o t na formulacao original
  res = 0
  for(i in 2:index-1){
    res = res + df$Ct[i]*exp(-M*(index-i-1))
  }
  return(res)
}

recruit_module = function(index, M){ #index é o t na formulacao original
  contador = 1
  res = 0
  while(perturbacoes$timing[contador] <= index){
    res = res + perturbacoes$R[contador] * exp(-M*(index-perturbacoes$timing[contador]))
    if(contador == nrow(perturbacoes)){break()}
    contador = contador + 1
  }
  return(res)
}

simulador = function(data, k, alpha,
                     beta, M, N0){
  res = c()
  for(i in 1:nrow(df)){
    if(i == 1){# t0
      res[i] = k * df$obseff.trips[i]^ alpha * exp(-M/2) * (
        N0 * exp(-M*i) -
          exp(-M/2) * (0) # catch
        + 0# recrutamento
      )^beta 
    }
    else{
      res[i] = k * df$obseff.trips[i]^ alpha * exp(-M/2) * (
        N0 * exp(-M*i) -
          exp(-M/2) * catch_module(i,M) # catch
        + recruit_module(i,M)# recrutamento
      )^beta
    }
  }
  
  return(res)
}


## funcao que corrige delta.glm

custom_delta_glm = function(input.data){
  
  input.data$year <- as.factor(input.data$year)
  input.data$fishing.season <- as.factor(input.data$fishing.season)
  input.data$rectangle <- as.factor(input.data$rectangle)
  input.data$power.class <- as.factor(input.data$power.class)
  input.data$factor.month = input.data$month
  # input.data[input.data$month == 1, ]$factor.month <- "01"
  # input.data[input.data$month == 2, ]$factor.month <- "02"
  # input.data[input.data$month == 3, ]$factor.month <- "03"
  # input.data[input.data$month == 4, ]$factor.month <- "04"
  # input.data[input.data$month == 5, ]$factor.month <- "05"
  # input.data[input.data$month == 6, ]$factor.month <- "06"
  # input.data[input.data$month == 7, ]$factor.month <- "07"
  # input.data[input.data$month == 8, ]$factor.month <- "08"
  # input.data[input.data$month == 9, ]$factor.month <- "09"
  # input.data[input.data$month == 10, ]$factor.month <- "10"
  # input.data[input.data$month == 11, ]$factor.month <- "11"
  # input.data[input.data$month == 12, ]$factor.month <- "12"
  input.data <- data.frame(input.data$year, input.data$fishing.season, 
                           input.data$factor.month, input.data$rectangle, input.data$power.class, 
                           input.data$lpue)
  colnames(input.data) <- c("year", "fishing.season", "month", 
                            "rectangle", "power.class", "lpue")
  input.data$presence <- 1
  input.data[input.data$lpue == 0, ]$presence <- 0
  binomial.glm <- glm(presence ~ fishing.season + month + rectangle + 
                        power.class, family = "binomial", data = input.data)
  binomial.summary <- summary(binomial.glm)
  binomial.residuals <- residuals(binomial.glm)
  binomial.fit <- fitted(binomial.glm)
  positive.input.data <- input.data[input.data$lpue > 0, ]
  gaussian.glm <- glm(log(lpue) ~ fishing.season + month + 
                        rectangle + power.class, family = "gaussian", data = positive.input.data)
  gaussian.summary <- summary(gaussian.glm)
  gaussian.residuals <- residuals(gaussian.glm)
  gaussian.fit <- fitted(gaussian.glm)
  positive.input.data <- input.data[input.data$lpue > 0, ]
  positive.input.data$year <- as.factor(as.character(positive.input.data$year))
  positive.input.data$fishing.season <- as.factor(as.character(positive.input.data$fishing.season))
  positive.input.data$month <- as.factor(as.character(positive.input.data$month))
  positive.input.data$rectangle <- as.factor(as.character(positive.input.data$rectangle))
  positive.input.data$power.class <- as.factor(as.character(positive.input.data$power.class))
  l.fishing.season <- length(levels(as.factor(positive.input.data$fishing.season)))
  l.month <- length(levels(as.factor(positive.input.data$month)))
  l.rectangle <- length(levels(as.factor(positive.input.data$rectangle)))
  l.power.class <- length(levels(as.factor(positive.input.data$power.class)))
  l.total <- l.fishing.season * l.month * l.rectangle * l.power.class
  predicted.lpue <- matrix(NA, nrow = l.total, ncol = 4)
  predicted.lpue[, 4] <- rep(levels(positive.input.data$power.class), 
                             (l.total/l.power.class))
  for (k in 1:l.fishing.season) {
    for (j in 1:l.rectangle) {
      for (m in 1:l.month) {
        start.fishing.season <- k * l.rectangle * l.month * 
          l.power.class - l.rectangle * l.month * l.power.class
        end.fishing.season <- k * l.rectangle * l.month * 
          l.power.class
        start.rectangle <- start.fishing.season + j * 
          l.month * l.power.class - l.month * l.power.class
        end.rectangle <- start.fishing.season + j * l.month * 
          l.power.class
        start.month <- start.rectangle + m * l.power.class - 
          l.power.class + 1
        end.month <- start.rectangle + m * l.power.class
        predicted.lpue[start.month:end.month, 3] <- rep(levels(positive.input.data$month)[m], 
                                                        l.power.class)
      }
      predicted.lpue[(start.rectangle + 1):end.rectangle, 
                     2] <- rep(levels(positive.input.data$rectangle)[j], 
                               l.month * l.power.class)
    }
    predicted.lpue[(start.fishing.season + 1):end.fishing.season, 
                   1] <- rep(levels(positive.input.data$fishing.season)[k], 
                             l.rectangle * l.month * l.power.class)
  }
  predicted.lpue <- as.data.frame(predicted.lpue)
  colnames(predicted.lpue) <- c("fishing.season", "rectangle", 
                                "month", "power.class")
  binomial.glm.prediction <- predict.glm(binomial.glm, predicted.lpue[, 
                                                                      1:4], type = "r", se.fit = T)
  gaussian.glm.prediction <- predict.glm(gaussian.glm, predicted.lpue[, 
                                                                      1:4], type = "r", se.fit = T)
  binomial.glm.prediction.fit <- binomial.glm.prediction$fit
  gaussian.glm.prediction.fit <- gaussian.glm.prediction$fit
  gaussian.glm.prediction.var <- (gaussian.glm.prediction$se.fit)^2
  gaussian.prediction <- exp(gaussian.glm.prediction.fit + 
                               gaussian.glm.prediction.var/2)
  predicted.lpue$st.lpue <- binomial.glm.prediction.fit * gaussian.prediction
  delta.outputs <- list(binomial.glm, binomial.summary, binomial.residuals, 
                        binomial.fit, gaussian.glm, gaussian.summary, gaussian.residuals, 
                        gaussian.fit, predicted.lpue)
  names(delta.outputs) <- c("binomial.glm", "binomial.summary", 
                            "binomial.residuals", "binomial.fit", "gaussian.glm", 
                            "gaussian.summary", "gaussian.residuals", "gaussian.fit", 
                            "predicted.lpue")
  
  
  return(delta.outputs) } 

