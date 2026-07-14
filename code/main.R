# title: "Escenarios de estrés en portafolios de inversiones: Un enfoque empírico"
# author: Stefan Bolta, FRM.
# affiliation: Superintendencia de Bancos
# date: "Junio 30, 2025"
# subtitle: Superintendencia de Bancos, Departamento de Estudios Económicos
# abstract: "Esta investigación presenta un proceso de construcción de escenarios de estrés 
# siguiendo principios puramente empíricos. La aplicación se desarrolla utilizando la técnica de meta-labeling 
# para identificar períodos con características similares en la serie de tasas de interés y, posteriormente, 
# analizar el comportamiento del índice benchmark durante esos períodos. Los datos se aproximan a la función de densidad 
# de probabilidad empírica mediante la solución cerrada de la familia de distribuciones Pearson Tipo IV. A partir de esto, 
# se generan los escenarios por medio de Simulación de Monte Carlo (SMC) utilizando datos sintéticos y se interpretan los resultados."  
#
# NOTA: Versión modular de MR_EscenariosEstres.R — ejecutar con el directorio de
# trabajo en la raíz del repositorio. Las funciones auxiliares se cargan desde
# code/functions.R; el resto del pipeline es idéntico al script original.

# Cargar librerias ----
require(quantmod)
require(PerformanceAnalytics)
require(roll)
require(pastecs)
require(stargazer)
require(knitr)
require(fAssets)
require(readxl)
require(dplyr)
require(PearsonDS)
require(bizdays)
require(fitdistrplus)
require(lubridate)
require(kableExtra)
require(cubature)

# Parámetros ----
# Valores ajustables centralizados. Cada uno sustituye 1:1 el literal usado en el
# script original, por lo que el comportamiento es idéntico por construcción.

# La simulación walk-forward OOS es computacionalmente costosa; por defecto se
# utilizan los resultados precomputados en empirical_VaR_rollingOOS_p*.csv
# (leídos en la sección OOS). Poner TRUE para re-ejecutar la simulación.
RUN_WALK_FORWARD <- FALSE

# Archivos / fuentes de datos
URL_GOBIX      <- "https://www.bvrd.com.do/indice/Data/GobixDataIRP.csv"
FILE_REGIMENES <- "Regimenes de TPM.xlsx"
FILE_VAR_P95   <- "empirical_VaR_rollingOOS_p95.csv"
FILE_VAR_P975  <- "empirical_VaR_rollingOOS_p975.csv"
FILE_VAR_P99   <- "empirical_VaR_rollingOOS_p99.csv"

# Ventanas y fechas
RANGO_GOBIX_2023  <- "2022-12-30/2023-12-15"
RANGO_GOBIX_2024  <- "2023-12-30/2024-12-15"
RANGO_GOBIX_OOS   <- "2022-11-22/2024-12-30"
RANGO_OOS         <- "2019::2022"
RANGO_ITS         <- "2014::2022"
RANGO_TRAIN       <- "2014::2018"
RANGO_BIII_STRESS <- "2021-12-01/2022"
FECHA_CORTE_ITS   <- "/2023-01-01"
FECHA_OOS_INICIO  <- "2023-03-01"
FECHA_OOS_FIN     <- "2024-12-30"
PERIODO_CENSURADO <- "2022-11-15/2023-03-01"

# Simulación / estimación
ROLL_WIDTH <- 65      # ventana rolling de la desviación estándar diaria
N_BOOT     <- 10000   # réplicas bootstrap en descdist / descdist_custom
N_SIM      <- 5000    # número de simulaciones Monte Carlo (ventana fija)
T_SIM      <- 250     # horizonte (días) de la simulación de ventana fija
# Nota: las semillas están fijadas dentro de las funciones (set.seed(1)) y en el
# ANEXO-01, igual que en el script original.

source("code/functions.R")

# Descarga GOBIXDR diario ----
gobix <- read.csv(url(URL_GOBIX)) # descarga
gobix <- gobix[,1:2] # se queda con dos primeras columnas
colnames(gobix)[1] <- "fecha" # renombra la columna 1
gobix$fecha <- lubridate::mdy(gobix[,1]) # la convierte en formato fechas mes-dia-año
gobix <- xts::as.xts(gobix[,-1], order.by=gobix$fecha) # convierte el formato en XTS
colnames(gobix)[1] <- "Close.GOBIX" # renombra la serie (indice al cierre del dia)

# VaR. Aplica el ejercicio del VaR historico visto en clase.
gobix$var.diaria <- CalculateReturns(Cl(gobix), method = "discrete") 
gobix$daily.vol <- diff(log(Cl(gobix)))
gobix$rolling.sd <- roll_sd(gobix$daily.vol, width = ROLL_WIDTH)*sqrt(252)
gobix <- na.omit(gobix)
gobix$drawdown <- Drawdowns(gobix[,'daily.vol'])

gobix.2023 <- gobix[RANGO_GOBIX_2023]
gobix.2024 <- gobix[RANGO_GOBIX_2024]
gobix.oos <- gobix[RANGO_GOBIX_OOS]

gobix.completo <- gobix

# OOS
oos <- gobix[RANGO_OOS] 

# Traininig: In the Sample
gobix <- gobix[RANGO_ITS] 
gobix.train <- gobix[RANGO_TRAIN] 


gobix$yesterday.ret <- data.table::shift(as.numeric(gobix[,'daily.vol',drop=FALSE]), n=1, type = "lag")
gobix <- na.omit(gobix)
# Asignacion regimenes ----
Regimenes_PM <- read_excel(FILE_REGIMENES)
Regimenes_PM <- as.data.frame(Regimenes_PM) %>% dplyr::filter(ANO >= 2014)

regimes.df <- as.data.frame(na.omit(Regimenes_PM))
regimes.df$TPM_move <- NA
for(i in 2:nrow(regimes.df)){
  if(regimes.df[i,'TPM'] > regimes.df[i-1,'TPM']){
    regimes.df[i,'TPM_move'] <- "HIGHER"
  } else if(regimes.df[i,'TPM'] < regimes.df[i-1,'TPM']){
    regimes.df[i,'TPM_move'] <- "LOWER"
  } else{
    regimes.df[i,'TPM_move'] <- "FLAT"
  }
}

regimes.df$CONDITIONAL <- paste0(regimes.df$ESTADO,"-",
                                 regimes.df$TENDENCIA,"-",
                                 regimes.df$TPM_move)

regimes.df <- na.omit(regimes.df)

regimes.list <- list()
for(i in 1:length(unique(regimes.df$CONDITIONAL))){
  regimes.list[[i]] <- regimes.df %>% dplyr::filter(CONDITIONAL == unique(regimes.df$CONDITIONAL)[i])
}

higher <- regimes.df %>% dplyr::filter(TPM_move == 'HIGHER')
flat <- regimes.df %>% dplyr::filter(TPM_move == 'FLAT')
lower <- regimes.df %>% dplyr::filter(TPM_move == 'LOWER')

gobix.df <- as.data.frame(gobix.completo)
gobix.df$YM <- as.yearmon(index(gobix.completo))

# MAPEAR DESDE EL CONTROL FEATURE (TASA MENSUAL) A LA SERIE DIARIA DEL GOBIX
periodo.lower <- as.yearmon(paste0(lower$ANO,"-",lower$MES))
periodo.flat <- as.yearmon(paste0(flat$ANO,"-",flat$MES))
periodo.higher <- as.yearmon(paste0(higher$ANO,"-",higher$MES))
# FALTA HACERLO EN LOOP PARA FILTRAR TODAS LAS CONFIGURACIONES
regime.lower <- list()
for(i in 1:length(periodo.lower)){
  regime.lower[[i]] <- gobix.df %>% dplyr::filter(YM == periodo.lower[i])
}

regime.lower <- do.call(rbind,regime.lower)
regime.lower <- regime.lower[,-ncol(regime.lower)]
regime.lower <- as.xts(regime.lower, order.by = ymd(rownames(regime.lower)))
regime.lower_full <- regime.lower
regime.lower <- regime.lower[FECHA_CORTE_ITS] # ITS

regime.flat <- list()
for(i in 1:length(periodo.flat)){
  regime.flat[[i]] <- gobix.df %>% dplyr::filter(YM == periodo.flat[i])
}

regime.flat <- do.call(rbind,regime.flat)
regime.flat <- regime.flat[,-ncol(regime.flat)]
regime.flat <- as.xts(regime.flat, order.by = ymd(rownames(regime.flat)))
regime.flat_full <- regime.flat
regime.flat <- regime.flat[FECHA_CORTE_ITS] # ITS

regime.higher <- list()
for(i in 1:length(periodo.higher)){
  regime.higher[[i]] <- gobix.df %>% dplyr::filter(YM == periodo.higher[i])
}

regime.higher <- do.call(rbind,regime.higher)
regime.higher <- regime.higher[,-ncol(regime.higher)]
regime.higher <- as.xts(regime.higher, order.by = ymd(rownames(regime.higher)))
regime.higher_full <- regime.higher
regime.higher <- regime.higher[FECHA_CORTE_ITS] # ITS

baselIII.12m.stress <- gobix[RANGO_BIII_STRESS] 

# CALIBRACION ----

# resumen autocorrelacion
autocorrelation.table <- cbind(table.Autocorrelation(gobix[,'daily.vol',drop=FALSE]), 
                               table.Autocorrelation(regime.higher[,'daily.vol',drop=FALSE]),
                               table.Autocorrelation(regime.flat[,'daily.vol',drop=FALSE]),
                               table.Autocorrelation(regime.lower[,'daily.vol',drop=FALSE]))

colnames(autocorrelation.table) <- c("todos","R(i=1)","R(i=2)","R(i=3)")

# estadisticas sumarizadas por regimen
unconditional.summary.stats <- stat.desc(gobix[,'daily.vol',drop=FALSE], basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)
regime.lower.summary.stats <- stat.desc(na.omit(regime.lower[,'daily.vol',drop=FALSE]), basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)
regime.flat.summary.stats <- stat.desc(na.omit(regime.flat[,'daily.vol',drop=FALSE]), basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)
regime.higher.summary.stats <- stat.desc(na.omit(regime.higher[,'daily.vol',drop=FALSE]), basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)
regime.stress.summary.stats <- stat.desc(na.omit(baselIII.12m.stress["2022-01-01/2022-06-30",'daily.vol',drop=FALSE]), basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)

summary.stats <- cbind(unconditional.summary.stats,
                       regime.lower.summary.stats,
                       regime.flat.summary.stats, 
                       regime.higher.summary.stats,
                       regime.stress.summary.stats)

summary.stats <- round(summary.stats,4)
  colnames(summary.stats) <- c('todos','reduccion','sin_cambios','incremento','estres')

summary.stats_short <- summary.stats[c('mean','std.dev','skewness','kurtosis'),]

# simulacion
regime.lower.sim <- descdist(as.numeric(regime.lower[,'daily.vol',drop=FALSE]), discrete=FALSE, boot=N_BOOT, graph = FALSE)
regime.flat.sim <- descdist(as.numeric(regime.flat[,'daily.vol',drop=FALSE]), discrete=FALSE, boot=N_BOOT, graph = FALSE)
regime.higher.sim <- descdist(as.numeric(regime.higher[,'daily.vol',drop=FALSE]), discrete=FALSE, boot=N_BOOT, graph = FALSE)
regime.all.sim <- descdist(as.numeric(gobix[,'daily.vol',drop=FALSE]), discrete=FALSE, boot=N_BOOT, graph = FALSE)
regime.stress.sim <- descdist(as.numeric(baselIII.12m.stress[,'daily.vol',drop=FALSE]), discrete=FALSE, boot=N_BOOT, graph = FALSE)

sim.list <- list(regime.all.sim, 
                 regime.lower.sim,
                 regime.flat.sim,
                 regime.higher.sim,
                 regime.stress.sim)

sim.params <- lapply(sim.list, FUN = function(x){extract_params(x)})
summary.stats_sim <- do.call(cbind, sim.params)
colnames(summary.stats_sim) <- c('todos','R(i=3)','R(i=2)','R(i=1)','BIII')

var <- summary.stats_sim['sd',]^2
summary.stats_sim <- rbind(summary.stats_sim,var)

df <- head(gobix.2023[,'Close.GOBIX'],1)

# Run only one piece 2023-2024
# At each month t, revise the Ri and simulate from correct distribution (updated)

# Later: Write a function that updates the parameters
forecast.estres <- list()
forecast.sin_cambios <- list()
forecast.incremento <- list()
forecast.reduccion <- list()

stress_mean <- list()
stress_variance <- list()
stress_skew <- list()
stress_kurt <- list()

sin.cambios_mean <- list()
sin.cambios_variance <- list()
sin.cambios_skew <- list()
sin.cambios_kurt <- list()

higher_mean <- list()
higher_variance <- list()
higher_skew <- list()
higher_kurt <- list()

lower_mean <- list()
lower_variance <- list()
lower_skew <- list()
lower_kurt <- list()

oos.start <- which(index(gobix.completo) == FECHA_OOS_INICIO)
oos.end <- which(index(gobix.completo) == FECHA_OOS_FIN)

# manage exclusions. 
censored.period <- PERIODO_CENSURADO
censored.dates <- index(gobix.completo[censored.period,'daily.vol'])
oos_sim.dates <- index(ExcludeDates(gobix.oos[,'daily.vol'], exclude = censored.dates))

for(i in oos.start:oos.end){
  
  date <- index(gobix.completo)[i]
  date.limit <- paste0(index(gobix.completo)[1],"/",date)
  #print(date)
  
  stress_mean[[i]] <- mean(ExcludeDates(baselIII.12m.stress[date.limit,'daily.vol'], exclude = censored.dates))
  stress_variance[[i]] <- sd(ExcludeDates(baselIII.12m.stress[date.limit,'daily.vol'], exclude = censored.dates))^2
  stress_skew[[i]] <- moments::skewness(ExcludeDates(baselIII.12m.stress[date.limit,'daily.vol'], exclude = censored.dates))
  stress_kurt[[i]] <- moments::kurtosis(ExcludeDates(baselIII.12m.stress[date.limit,'daily.vol'], exclude = censored.dates))
  
  sin.cambios_mean[[i]] <- mean(ExcludeDates(regime.flat_full[date.limit,'daily.vol'], exclude = censored.dates))
  sin.cambios_variance[[i]] <- sd(ExcludeDates(regime.flat_full[date.limit,'daily.vol'], exclude = censored.dates))^2
  sin.cambios_skew[[i]] <- moments::skewness(ExcludeDates(regime.flat_full[date.limit,'daily.vol'], exclude = censored.dates))
  sin.cambios_kurt[[i]] <- moments::kurtosis(ExcludeDates(regime.flat_full[date.limit,'daily.vol'], exclude = censored.dates))
  
  higher_mean[[i]] <- mean(ExcludeDates(regime.higher_full[date.limit,'daily.vol'], exclude = censored.dates))
  higher_variance[[i]] <- sd(ExcludeDates(regime.higher_full[date.limit,'daily.vol'], exclude = censored.dates))^2
  higher_skew[[i]] <- moments::skewness(ExcludeDates(regime.higher_full[date.limit,'daily.vol'], exclude = censored.dates))
  higher_kurt[[i]] <- moments::kurtosis(ExcludeDates(regime.higher_full[date.limit,'daily.vol'], exclude = censored.dates))
  
  lower_mean[[i]] <- mean(ExcludeDates(regime.lower_full[date.limit,'daily.vol'], exclude = censored.dates))
  lower_variance[[i]] <- sd(ExcludeDates(regime.lower_full[date.limit,'daily.vol'], exclude = censored.dates))^2
  lower_skew[[i]] <- moments::skewness(ExcludeDates(regime.lower_full[date.limit,'daily.vol'], exclude = censored.dates))
  lower_kurt[[i]] <- moments::kurtosis(ExcludeDates(regime.lower_full[date.limit,'daily.vol'], exclude = censored.dates))
  
  # RESULTADOS EXPORTADOS a los archivos:
  # empirical_VaR_rollingOOS_p95.csv, empirical_VaR_rollingOOS_p975.csv, empirical_VaR_rollingOOS_p99.csv
  # Por defecto (RUN_WALK_FORWARD = FALSE) NO se re-ejecuta esta simulación:
  # los CSV precomputados se leen más abajo en la sección OOS.
  if (RUN_WALK_FORWARD) {
    forecast.estres[[i]] <- empirical.simulation.simple(data = gobix.completo[i,],
                                                        n.sim = 5000,
                                                        t = 30,
                                                        p = 0.01,
                                                        mean = stress_mean[[i]],
                                                        variance = stress_variance[[i]],
                                                        skew = stress_skew[[i]],
                                                        kurt = stress_kurt[[i]])

    forecast.sin_cambios[[i]] <- empirical.simulation.simple(data = gobix.completo[i,],
                                                            n.sim = 5000,
                                                            t = 30,
                                                            p = 0.01,
                                                            mean = sin.cambios_mean[[i]],
                                                            variance = sin.cambios_variance[[i]],
                                                            skew = sin.cambios_skew[[i]],
                                                            kurt = sin.cambios_kurt[[i]])

    forecast.incremento[[i]] <- empirical.simulation.simple(data = gobix.completo[i,],
                                                            n.sim = 5000,
                                                            t = 30,
                                                            p = 0.01,
                                                            mean = higher_mean[[i]],
                                                            variance = higher_variance[[i]],
                                                            skew = higher_skew[[i]],
                                                            kurt = higher_kurt[[i]])

    forecast.reduccion[[i]] <- empirical.simulation.simple(data = gobix.completo[i,],
                                                            n.sim = 5000,
                                                            t = 30,
                                                            p = 0.01,
                                                            mean = lower_mean[[i]],
                                                            variance = lower_variance[[i]],
                                                            skew = lower_skew[[i]],
                                                            kurt = lower_kurt[[i]])
  }
}

# EVOLUCION DE STD.DEV
stress_variance.df <- do.call(rbind,stress_variance)
sin.cambios_variance.df <- do.call(rbind,sin.cambios_variance)
highers_variance.df <- do.call(rbind,higher_variance)
lower_variance.df <- do.call(rbind,lower_variance)

variance.df <- data.frame(lower_variance.df, 
                          sin.cambios_variance.df, 
                          highers_variance.df, 
                          stress_variance.df)

colnames(variance.df) <- c("lower","sin_cambios","higher","BIII.stress")

variance.df <- as.data.frame(apply(variance.df, MARGIN = 2, FUN = function(x){ sqrt(x)*sqrt(250) }))
  variance.df <- variance.df[-1,]
  variance.df$date <- oos_sim.dates

# EVOLUCION DE SKEW
stress_skew.df <- do.call(rbind,stress_skew)
sin.cambios_skew.df <- do.call(rbind,sin.cambios_skew)
highers_skew.df <- do.call(rbind,higher_skew)
lower_skew.df <- do.call(rbind,lower_skew)

skew.df <- data.frame(lower_skew.df, 
                      sin.cambios_skew.df, 
                      highers_skew.df, 
                      stress_skew.df)

  colnames(skew.df) <- c("lower","sin_cambios","higher","BIII.stress")
  skew.df <- skew.df[-1,]
  skew.df$date <- oos_sim.dates

# EVOLUCION DE KURT
stress_kurt.df <- do.call(rbind,stress_kurt)
sin.cambios_kurt.df <- do.call(rbind,sin.cambios_kurt)
highers_kurt.df <- do.call(rbind,higher_kurt)
lower_kurt.df <- do.call(rbind,lower_kurt)

kurt.df <- data.frame(lower_kurt.df, 
                      sin.cambios_kurt.df, 
                      highers_kurt.df, 
                      stress_kurt.df)

  colnames(kurt.df) <- c("lower","sin_cambios","higher","BIII.stress")
  kurt.df <- kurt.df[-1,]
  kurt.df$date <- oos_sim.dates

# Estimacion Pearson IV ----
a <- round(ajustar_pearson_iv(na.omit(gobix[,'daily.vol',drop=FALSE])),4)
b <- round(ajustar_pearson_iv(na.omit(regime.lower[,'daily.vol',drop=FALSE])),4)
c <- round(ajustar_pearson_iv(na.omit(regime.flat[,'daily.vol',drop=FALSE])),4)
d <- round(ajustar_pearson_iv(na.omit(regime.higher[,'daily.vol',drop=FALSE])),4)
e <- round(ajustar_pearson_iv(na.omit(baselIII.12m.stress[,'daily.vol',drop=FALSE])),4)

estimation.results <- cbind(a,b,c,d,e)
  colnames(estimation.results) <- c("todos","Ri=3","Ri=2","Ri=1","BIII")
  
# OOS ----
# FORECAST OOS
# carga los objetos pre-procesados.
fwd_p95.df <- read.csv(FILE_VAR_P95)
fwd_p975.df <- read.csv(FILE_VAR_P975)
fwd_p99.df <- read.csv(FILE_VAR_P99)

fwd_oos.all <- list(fwd_p95.df, fwd_p975.df, fwd_p99.df)

fwd.df <- lapply(fwd_oos.all, FUN = function(x){apply(x[,2:ncol(x)], MARGIN = 2, FUN = function(x){quantile(x,0.5)})})
fwd.df <- do.call(rbind,fwd.df)
colnames(fwd.df) <- c('BIII','R(i=1)','R(i=2)','R(i=3)')
rownames(fwd.df) <- c('VaR.950','VaR.975','VaR.990')

fwd_all.df1 <- as.data.frame(fwd_oos.all[[3]])
fwd_all.df1$dates <- lubridate::ymd(fwd_all.df1$X)

selected.range <- paste0(index(df),"/",index(df)+360)
oos.data <- gobix[selected.range][,1]

sim.stress <- empirical.simulation(data = df, oos.series = NULL, n.sim = N_SIM, t = T_SIM, mean=summary.stats_sim['mean','BIII'], variance=summary.stats_sim['var','BIII'], skew=summary.stats_sim['skewness','BIII'], kurt=summary.stats_sim['kurtosis','BIII'], main="FHS Stress Times")
sim.higher <- empirical.simulation(data = df,  oos.series = NULL, n.sim = N_SIM, t = T_SIM, mean=summary.stats_sim['mean','R(i=1)'], variance=summary.stats_sim['var','R(i=1)'], skew=summary.stats_sim['skewness','R(i=1)'], kurt=summary.stats_sim['kurtosis','R(i=1)'], main="FHS Higher Rates")
sim.lower <- empirical.simulation(data = df,  oos.series = NULL, n.sim = N_SIM, t = T_SIM, mean=summary.stats_sim['mean','R(i=3)'], variance=summary.stats_sim['var','R(i=3)'], skew=summary.stats_sim['skewness','R(i=3)'], kurt=summary.stats_sim['kurtosis','R(i=3)'], main="FHS Lower Rates")
sim.flat <- empirical.simulation(data = df,  oos.series = NULL, n.sim = N_SIM, t = T_SIM, mean=summary.stats_sim['mean','R(i=2)'], variance=summary.stats_sim['var','R(i=2)'], skew=summary.stats_sim['skewness','R(i=2)'], kurt=summary.stats_sim['kurtosis','R(i=2)'], main="FHS No Changes in Rates")

x.higher <- quantile(na.omit(as.numeric(sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),])),0.50)/100-1
x.lower <- quantile(na.omit(as.numeric(sim.lower[[1]][[1]][nrow(sim.lower[[1]][[1]]),])),0.50)/100-1
x.flat <- quantile(na.omit(as.numeric(sim.flat[[1]][[1]][nrow(sim.flat[[1]][[1]]),])),0.50)/100-1
x.stress <- quantile(na.omit(as.numeric(sim.stress[[1]][[1]][nrow(sim.stress[[1]][[1]]),])),0.50)/100-1

higher.q05 <- quantile(na.omit(as.numeric(sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),])),0.05)/100-1
lower.q05 <- quantile(na.omit(as.numeric(sim.lower[[1]][[1]][nrow(sim.lower[[1]][[1]]),])),0.05)/100-1
flat.q05 <- quantile(na.omit(as.numeric(sim.flat[[1]][[1]][nrow(sim.flat[[1]][[1]]),])),0.05)/100-1
stress.q05 <- quantile(na.omit(as.numeric(sim.stress[[1]][[1]][nrow(sim.stress[[1]][[1]]),])),0.05)/100-1

higher.q025 <- quantile(na.omit(as.numeric(sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),])),0.025)/100-1
lower.q025 <- quantile(na.omit(as.numeric(sim.lower[[1]][[1]][nrow(sim.lower[[1]][[1]]),])),0.025)/100-1
flat.q025 <- quantile(na.omit(as.numeric(sim.flat[[1]][[1]][nrow(sim.flat[[1]][[1]]),])),0.025)/100-1
stress.q025 <- quantile(na.omit(as.numeric(sim.stress[[1]][[1]][nrow(sim.stress[[1]][[1]]),])),0.025)/100-1

higher.q01 <- quantile(na.omit(as.numeric(sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),])),0.01)/100-1
lower.q01 <- quantile(na.omit(as.numeric(sim.lower[[1]][[1]][nrow(sim.lower[[1]][[1]]),])),0.01)/100-1
flat.q01 <- quantile(na.omit(as.numeric(sim.flat[[1]][[1]][nrow(sim.flat[[1]][[1]]),])),0.01)/100-1
stress.q01 <- quantile(na.omit(as.numeric(sim.stress[[1]][[1]][nrow(sim.stress[[1]][[1]]),])),0.01)/100-1

# EXTRACT PCENTILES

P01.sim <- data.frame(sim.stress[[1]][[1]][,'Q010'],
                      sim.higher[[1]][[1]][,'Q010'],
                      sim.flat[[1]][[1]][,'Q010'],
                      sim.lower[[1]][[1]][,'Q010'])
  colnames(P01.sim) <- c("BaselIII.stress","higher.rate","no.change","lower.rate")
  P01.sim <- as.xts(P01.sim, order.by = ymd(rownames(P01.sim)))


P025.sim <- data.frame(sim.stress[[1]][[1]][,'Q025'],
                       sim.higher[[1]][[1]][,'Q025'],
                       sim.flat[[1]][[1]][,'Q025'],
                       sim.lower[[1]][[1]][,'Q025'])
  colnames(P025.sim) <- c("BaselIII.stress","higher.rate","no.change","lower.rate")
  P025.sim <- as.xts(P025.sim, order.by = ymd(rownames(P025.sim)))

P05.sim <- data.frame(sim.stress[[1]][[1]][,'Q050'],
                      sim.higher[[1]][[1]][,'Q050'],
                      sim.flat[[1]][[1]][,'Q050'],
                      sim.lower[[1]][[1]][,'Q050'])
  colnames(P05.sim) <- c("BaselIII.stress","higher.rate","no.change","lower.rate")
  P05.sim <- as.xts(P05.sim, order.by = ymd(rownames(P05.sim)))

# merge
gobix.2023$indexed <- cumprod(1+gobix.2023$var.diaria)*100

# normalizacion
P01.sim <- na.omit(merge(P01.sim,gobix.2023$indexed))
  P01.sim_df <- as.data.frame(P01.sim)
  P01.sim_df <- as.data.frame(apply(P01.sim_df, MARGIN = 2, FUN = function(x){(x-100)/100}))
  P01.sim_df$dates <- lubridate::ymd(rownames(P01.sim_df))

# tabla
sim.resume <- as.data.frame(rbind(tail(P05.sim,1), tail(P025.sim,1), tail(P01.sim[,-ncol(P01.sim)],1)))
  sim.resume <- round(sim.resume, 2)
  sim.resume <- apply(sim.resume, FUN = function(x){ x-100}, MARGIN = 2)
  rownames(sim.resume) <- c("VaR.950","VaR.975","VaR.990")

# RESUMIR EN UNA TABLA PARA CADA PARAMETRO (mu,sigma,kurt,skew)
sim_lower.regime.mean <- apply(sim.lower[[1]][[1]], MARGIN = 2, FUN = function(x){mean(diff(log(x)))*sqrt(250)})
sim_flat.regime.mean <- apply(sim.flat[[1]][[1]], MARGIN = 2, FUN = function(x){mean(diff(log(x)))*sqrt(250)})
sim_higher.regime.mean <- apply(sim.higher[[1]][[1]], MARGIN = 2, FUN = function(x){mean(diff(log(x)))*sqrt(250)})
sim_stress.regime.mean <- apply(sim.stress[[1]][[1]], MARGIN = 2, FUN = function(x){mean(diff(log(x)))*sqrt(250)})

sim_lower.regime.sd <- apply(sim.lower[[1]][[1]], MARGIN = 2, FUN = function(x){sd(diff(log(x)))*sqrt(250)})
sim_flat.regime.sd <- apply(sim.flat[[1]][[1]], MARGIN = 2, FUN = function(x){sd(diff(log(x)))*sqrt(250)})
sim_higher.regime.sd <- apply(sim.higher[[1]][[1]], MARGIN = 2, FUN = function(x){sd(diff(log(x)))*sqrt(250)})
sim_stress.regime.sd <- apply(sim.stress[[1]][[1]], MARGIN = 2, FUN = function(x){sd(diff(log(x)))*sqrt(250)})

sim_lower.regime.skew <- apply(sim.lower[[1]][[1]], MARGIN = 2, FUN = function(x){moments::skewness(diff(log(x)))})
sim_flat.regime.skew <- apply(sim.flat[[1]][[1]], MARGIN = 2, FUN = function(x){moments::skewness(diff(log(x)))})
sim_higher.regime.skew <- apply(sim.higher[[1]][[1]], MARGIN = 2, FUN = function(x){moments::skewness(diff(log(x)))})
sim_stress.regime.skew <- apply(sim.stress[[1]][[1]], MARGIN = 2, FUN = function(x){moments::skewness(diff(log(x)))})

sim_lower.regime.kurt <- apply(sim.lower[[1]][[1]], MARGIN = 2, FUN = function(x){moments::kurtosis(diff(log(x)))})
sim_flat.regime.kurt <- apply(sim.flat[[1]][[1]], MARGIN = 2, FUN = function(x){moments::kurtosis(diff(log(x)))})
sim_higher.regime.kurt <- apply(sim.higher[[1]][[1]], MARGIN = 2, FUN = function(x){moments::kurtosis(diff(log(x)))})
sim_stress.regime.kurt <- apply(sim.stress[[1]][[1]], MARGIN = 2, FUN = function(x){moments::kurtosis(diff(log(x)))})

# compute the distribution of parameter uncertainty
sim.regime.kurt <- list(sim_lower.regime.kurt,
                        sim_flat.regime.kurt,
                        sim_higher.regime.kurt,
                        sim_stress.regime.kurt)

sim.regime.skew <- list(sim_lower.regime.skew,
                        sim_flat.regime.skew,
                        sim_higher.regime.skew,
                        sim_stress.regime.skew)

sim.regime.sd <- list(sim_lower.regime.sd,
                      sim_flat.regime.sd,
                      sim_higher.regime.sd,
                      sim_stress.regime.sd)

sim.regime.mean <- list(sim_lower.regime.mean,
                        sim_flat.regime.mean,
                        sim_higher.regime.mean,
                        sim_stress.regime.mean)


fourth_moment.sim <- sapply(sim.regime.kurt, FUN = function(x){compute_stats(object = x)})
  colnames(fourth_moment.sim) <- c("R_(i=3)","R_(i=2)","R_(i=1)","BIII")

third_moment.sim <- sapply(sim.regime.skew, FUN = function(x){compute_stats(object = x)})
  colnames(third_moment.sim) <- c("R_(i=3)","R_(i=2)","R_(i=1)","BIII")

second_moment.sim <- sapply(sim.regime.sd, FUN = function(x){compute_stats(object = x)})
  colnames(second_moment.sim) <- c("R_(i=3)","R_(i=2)","R_(i=1)","BIII")

first_moment.sim <- sapply(sim.regime.mean, FUN = function(x){compute_stats(object = x)})
  colnames(first_moment.sim) <- c("R_(i=3)","R_(i=2)","R_(i=1)","BIII")


# GRAFICOS ----
  
# graph-01: Serie de tiempo x regimen
par(mfrow=c(1,1), mar=c(1,1,1,1))
plot.1 <- plot(gobix[,'var.diaria'], main=" Retornos diarios GOBIXDR segregado por regímen R(i) ", type="l", xlab="T", grid.col=NA, lwd=0.85, col="grey20")
plot.1 <- addLegend("bottomleft", legend.names=c('reducción','sin cambios','incremento'), pch=16, col=c("forestgreen","blue","red"), bty="n", y.intersp = 0.75)
plot.1 <- lines(regime.lower[,'var.diaria'], on=NA, col="forestgreen", main="R(i=3,t:T)")
plot.1 <- lines(regime.flat[,'var.diaria'], on=NA, col="blue", main="R(i=2,t:T)")
plot.1 <- lines(regime.higher[,'var.diaria'], on=NA, col="red", main="R(i=1,t:T)")
plot.1

# graph-02: PDF historica x regimen
par(mfrow=c(1,1), mar=c(2,2,1,2))
plot.pdf(gobix[,'daily.vol',drop=FALSE], breaks=64, col='grey80', title="Distribución retornos por regímen R(i,t:T)")
lines(density(regime.flat[,'daily.vol'], bw=sd(regime.flat[,'daily.vol'])), lwd=2, col="blue")
lines(density(regime.lower[,'daily.vol'], bw=sd(regime.lower[,'daily.vol'])), lwd=2, col="forestgreen")
lines(density(regime.higher[,'daily.vol'], bw=sd(regime.higher[,'daily.vol'])), lwd=2, col="red")
legend("topleft", legend = c('R(i=1)',
                             'R(i=2)',
                             'R(i=3)',
                             'incondicional',
                             'realizado'), lwd=2, col = c('red',
                                                          'blue',
                                                          'forestgreen',
                                                          'grey50',
                                                          'grey'), lty=c(1,1,1,2,1), bty="n", cex=0.70)

#
lower.rate <- na.omit(as.numeric(sim.lower[[1]][[1]][nrow(sim.lower[[1]][[1]]),]))
unchanged.rate <- na.omit(as.numeric(sim.flat[[1]][[1]][nrow(sim.flat[[1]][[1]]),]))
higher.rate <- na.omit(as.numeric(sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),]))
baselIII.stress <- na.omit(as.numeric(sim.stress[[1]][[1]][nrow(sim.stress[[1]][[1]]),]))

data.box <- data.frame(lower.rate, unchanged.rate, higher.rate, baselIII.stress)

# transform from index to pct
sim.stress_250d <- sim.stress[[1]][[1]][nrow(sim.stress[[1]][[1]]),]
sim.stress_250d <- as.numeric(sim.stress_250d) / 100 -1

sim.higher_250d <- sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),]
sim.higher_250d <- as.numeric(sim.higher_250d) / 100 -1

sim.flat_250d <- sim.flat[[1]][[1]][nrow(sim.flat[[1]][[1]]),]
sim.flat_250d <- as.numeric(sim.flat_250d) / 100 -1

sim.lower_250d <- sim.lower[[1]][[1]][nrow(sim.lower[[1]][[1]]),]
sim.lower_250d <- as.numeric(sim.lower_250d) / 100 -1

# graph-03: CDF simulacion ventana fija OOS
par(mfrow=c(1,1), mar=c(2,3,3,3))
plot(ecdf(na.omit(as.numeric(sim.higher_250d))),
     xlab="",
     ylab="Densidad cumulativa",
     ylim=c(0,0.10),
     xlim=c(-0.60,0.20),
     main="Estimación VaR ventana fija - por estado R(i)",
     col="red", lwd=1.5)
lines(ecdf(na.omit(as.numeric(sim.lower_250d))), col="forestgreen", lwd=1.5)
lines(ecdf(na.omit(as.numeric(sim.flat_250d))), col="blue", lwd=1.5)
lines(ecdf(na.omit(as.numeric(sim.stress_250d))), col="red4", lwd=1.5)
abline(h=0.50, col='grey', lty=2, lwd=0.75)
mtext("Nivel cumulativo esperado en P(T) (inicial P(t) = 100)",  side=3, cex = 0.80)
mtext("(truncado)",  side=3, adj = 0, cex = 0.80)
mtext("T: 1A futuro",  side=3, adj = 1, cex = 0.80)

points(x = stress.q05, y = 0.05, pch=16, col="red4")
points(x = higher.q05, y = 0.05, pch=16, col="red")
points(x = flat.q05, y = 0.05, pch=16, col="blue")
points(x = lower.q05, y = 0.05, pch=16, col="green4")

points(x = stress.q025, y = 0.025, pch=17, col="red4", bg="white", lwd=1)
points(x = higher.q025, y = 0.025, pch=17, col="red", bg="white", lwd=1)
points(x = flat.q025, y = 0.025, pch=17, col="blue", bg="white", lwd=1)
points(x = lower.q025, y = 0.025, pch=17, col="green4", bg="white", lwd=1)

points(x = stress.q01, y = 0.01, pch=15, col="red4", bg="white", lwd=1)
points(x = higher.q01, y = 0.01, pch=15, col="red", bg="white", lwd=1)
points(x = flat.q01, y = 0.01, pch=15, col="blue", bg="white", lwd=1)
points(x = lower.q01, y = 0.01, pch=15, col="green4", bg="white", lwd=1)

text(x=as.numeric(stress.q01)-0.14, y=0.01, labels=paste0(round(as.numeric(stress.q01)*100,2),"%"), pos=4, cex=1, col = "red4")
text(x=as.numeric(stress.q025)-0.12, y=0.025, labels=paste0(round(as.numeric(stress.q025)*100,2),"%"), pos=4, cex=1, col = "red4")
text(x=as.numeric(stress.q05)-0.14, y=0.05, labels=paste0(round(as.numeric(stress.q05)*100,2),"%"), pos=4, cex=1, col = "red4")

text(x=as.numeric(higher.q01)+0.02, y=0.01, labels=paste0(round(as.numeric(higher.q01)*100,2),"%"), pos=4, cex=1, col = "red")
text(x=as.numeric(higher.q025)+0.01, y=0.025, labels=paste0(round(as.numeric(higher.q025)*100,2),"%"), pos=4, cex=1, col = "red")
text(x=as.numeric(higher.q05)+0.01, y=0.05, labels=paste0(round(as.numeric(higher.q05)*100,2),"%"), pos=4, cex=1, col = "red")

text(x=as.numeric(lower.q01)+0.01, y=0.01, labels=paste0(round(as.numeric(lower.q01)*100,2),"%"), pos=4, cex=1, col = "green4")
text(x=as.numeric(lower.q025)+0.01, y=0.025, labels=paste0(round(as.numeric(lower.q025)*100,2),"%"), pos=4, cex=1, col = "green4")
text(x=as.numeric(lower.q05)+0.01, y=0.05, labels=paste0(round(as.numeric(lower.q05)*100,2),"%"), pos=4, cex=1, col = "green4")

text(x=as.numeric(flat.q01)-0.13, y=0.01, labels=paste0(round(as.numeric(flat.q01)*100,2),"%"), pos=4, cex=1, col = "blue")
text(x=as.numeric(flat.q025)-0.09, y=0.025, labels=paste0(round(as.numeric(flat.q025)*100,2),"%"), pos=4, cex=1, col = "blue")
text(x=as.numeric(flat.q05)-0.12, y=0.05, labels=paste0(round(as.numeric(flat.q05)*100,2),"%"), pos=4, cex=1, col = "blue")

legend("topleft", legend=c('R(i=3): reducción',
                           'R(i=2): sin cambios',
                           'R(i=1): incremento',
                           'Basel_III estres'),
       col=c("forestgreen","blue","red","red4"), text.col=c("forestgreen","blue","red","red4"), lty=1, border.col = "white", border.lwd = 0, cex = 0.90)
box(col="grey")

legend("topright", legend=c('P.950',
                            'P.975',
                            'P.990'),
       col=c("grey60"), text.col = "grey60", pch = c(16,17,15), lwd=0, cex=1, lty=1, bty="n", y.intersp = 1.10)
box(col="grey")
par(mfrow=c(1,1))

# graph-04: PDF std dev estimacion ventana fija
par(mfrow=c(1,1), mar=c(2,3,2.5,3))
hist(unlist(sim_lower.regime.sd), breaks=64, xlim=c(0,0.25), ylim=c(0,1100), col="white", xlab="", ylab="", main="GOBIXDR - Distribución marginal condicional", border = "white")
hist(sim_lower.regime.sd, breaks=64, add=TRUE, col="green4")
rug(as.numeric(sim_lower.regime.sd), ticksize = 0.02, lwd=0.5, col="green4")
hist(sim_flat.regime.sd, breaks=64, add=TRUE, col="blue")
rug(as.numeric(sim_flat.regime.sd), ticksize = 0.02, lwd=0.5, col="blue")
hist(sim_higher.regime.sd, breaks=64, add=TRUE, col="red")
rug(as.numeric(sim_higher.regime.sd), ticksize = 0.02, lwd=0.5, col="red")
hist(sim_stress.regime.sd, breaks=64, add=TRUE, col="red4")
rug(as.numeric(sim_stress.regime.sd), ticksize = 0.02, lwd=0.5, col="red4")
legend("topright", legend=c('R(i=3)',
                            'R(i=2)',
                            'R(i=1)',
                            'Basel_III estres'),
       col=c("green4","blue","red","red4"), text.col=c("green4","blue","red","red4"), pch=15,  border.col = "white", bty="n", cex=1)
mtext("Desviación estándar anualizada",  side=3, cex = 0.90)  
box(col="grey")

# graph-05: CDF simulacion walk-forward
par(mfrow=c(1,1), mar=c(2,3,3,3))

lower.q05 <- quantile(na.omit(as.numeric(fwd_all.df1$Q01.reduccion)),0.05)
lower.q025 <- quantile(na.omit(as.numeric(fwd_all.df1$Q01.reduccion)),0.025)
lower.q01 <- quantile(na.omit(as.numeric(fwd_all.df1$Q01.reduccion)),0.01)

flat.q05 <- quantile(na.omit(as.numeric(fwd_all.df1$Q01.sin_cambio)),0.05)
flat.q025 <- quantile(na.omit(as.numeric(fwd_all.df1$Q01.sin_cambio)),0.025)
flat.q01 <- quantile(na.omit(as.numeric(fwd_all.df1$Q01.sin_cambio)),0.01)

higher.q05 <- quantile(na.omit(as.numeric(fwd_all.df1$Q01.incremento)),0.05)
higher.q025 <- quantile(na.omit(as.numeric(fwd_all.df1$Q01.incremento)),0.025)
higher.q01 <- quantile(na.omit(as.numeric(fwd_all.df1$Q01.incremento)),0.01)

stress.q05 <- quantile(na.omit(as.numeric(fwd_all.df1$Q01.estres)),0.05)
stress.q025 <- quantile(na.omit(as.numeric(fwd_all.df1$Q01.estres)),0.025)
stress.q01 <- quantile(na.omit(as.numeric(fwd_all.df1$Q01.estres)),0.01)

plot(ecdf(na.omit(as.numeric(fwd_all.df1$Q01.estres))),
     xlab="Nivel de pérdida esperado en t",
     ylab="Densidad cumulativa",
     ylim=c(0,0.2),
     xlim=c(-0.25,0),
     main="Estimación VaR.99 walk-forward - por estado R(i)",
     col="red4", lwd=1)
lines(ecdf(na.omit(as.numeric(fwd_all.df1$Q01.reduccion))), col="green4", lwd=1)
lines(ecdf(na.omit(as.numeric(fwd_all.df1$Q01.sin_cambio))), col="blue", lwd=1)
lines(ecdf(na.omit(as.numeric(fwd_all.df1$Q01.incremento))), col="red", lwd=1)
mtext("CDF en intervalo (t,T)",  side=3, cex = 0.80)
mtext("(truncado)",  side=3, adj = 0, cex = 0.80)
mtext("(t:T) = 18 meses",  side=3, adj = 1, cex = 0.80)

text(x=round(as.numeric(stress.q01),4)-0.04, y=0.01, labels=paste0(round(as.numeric(stress.q01)*100,2),"%"), pos=4, cex=1, col = "red4")
text(x=round(as.numeric(stress.q025),4)-0.035, y=0.025, labels=paste0(round(as.numeric(stress.q025)*100,2),"%"), pos=4, cex=1, col = "red4")
text(x=round(as.numeric(stress.q05),4)-0.035, y=0.05, labels=paste0(round(as.numeric(stress.q05)*100,2),"%"), pos=4, cex=1, col = "red4")

text(x=round(as.numeric(higher.q01),4)+0.001, y=0.01, labels=paste0(round(as.numeric(higher.q01)*100,2),"%"), pos=4, cex=1, col = "red")
text(x=round(as.numeric(higher.q025),4)+0.001, y=0.025, labels=paste0(round(as.numeric(higher.q025)*100,2),"%"), pos=4, cex=1, col = "red")
text(x=round(as.numeric(higher.q05),4)+0.005, y=0.05, labels=paste0(round(as.numeric(higher.q05)*100,2),"%"), pos=4, cex=1, col = "red")

text(x=round(as.numeric(lower.q01),4)+0.001, y=0.01, labels=paste0(round(as.numeric(lower.q01)*100,2),"%"), pos=4, cex=1, col = "green4")
text(x=round(as.numeric(lower.q025),4)+0.001, y=0.025, labels=paste0(round(as.numeric(lower.q025)*100,2),"%"), pos=4, cex=1, col = "green4")
text(x=round(as.numeric(lower.q05),4)+0.001, y=0.05, labels=paste0(round(as.numeric(lower.q05)*100,2),"%"), pos=4, cex=1, col = "green4")

text(x=round(as.numeric(flat.q01),4)-0.04, y=0.01, labels=paste0(round(as.numeric(flat.q01)*100,2),"%"), pos=4, cex=1, col = "blue")
text(x=round(as.numeric(flat.q025),4)-0.035, y=0.025, labels=paste0(round(as.numeric(flat.q025)*100,2),"%"), pos=4, cex=1, col = "blue")
text(x=round(as.numeric(flat.q05),4)-0.035, y=0.05, labels=paste0(round(as.numeric(flat.q05)*100,2),"%"), pos=4, cex=1, col = "blue")

legend("topleft", legend=c('R(i=3)',
                           'R(i=2)',
                           'R(i=1)',
                           'Basel_III.stress'),
       col=c("green4","blue","red","red4"), text.col=c("green4","blue","red","red4"), pch=16, border.col = "white", bty="n")
box(col="grey")
par(mfrow=c(1,1), mar=c(4,4,4,4))

# tabla 01: probabilidades de transicion ----
require(kableExtra)
require(markovchain)

mcFit.regime_transitions <- markovchainFit(data=Regimenes_PM$RATE)
regime_transition.probs <- as.data.frame(mcFit.regime_transitions$estimate@transitionMatrix)
regime_transition.probs <- round(regime_transition.probs,4)

  colnames(regime_transition.probs) <- c("R(i=3,t+1)","R(i=2,t+1)","R(i=1,t+1)")
  rownames(regime_transition.probs) <- c("R(i=3,t)","R(i=2,t)","R(i=1,t)")

knitr::kable(regime_transition.probs, caption="Probabilidad de transición entre regímenes TPM") %>% kable_styling(latex_options = "hold", position = "center")

# tabla 02: resumen parametros ----
colnames(summary.stats_short) <- c("todos","R(i=3)","R(i=2)","R(i=1)","BIII.estres")
knitr::kable(summary.stats_short[,-ncol(summary.stats_short)], caption="Resumen parámetros: período entrenamiento") %>%
  kable_styling(latex_options = "HOLD_position")
# tabla 03: autocorrelacion ----
knitr::kable(autocorrelation.table, caption="Ljung-Box Q") %>% kable_styling(latex_options = "hold", position = "center")
# tabla 04: parametros estimados ----
knitr::kable(estimation.results, caption="Parámetros estimados") %>%
  kable_styling(latex_options = "HOLD_position")
# tabla 05: parametros simulacion ----
summary.stats_sim[1,] <- summary.stats_sim[1,]*256
summary.stats_sim[2,] <- summary.stats_sim[2,]*sqrt(256)
rownames(summary.stats_sim)[1:2] <- c("media (ann.)","std.dev (ann.)")
summary.stats_sim <- round(summary.stats_sim,4)
knitr::kable(summary.stats_sim[-nrow(summary.stats_sim),], caption="Parámetros estimados") %>%
  kable_styling(latex_options = "HOLD_position")
# tabla 06: resumen walk-forward -----
fwd.df.summary.stats <- lapply(fwd_oos.all, FUN = function(x){apply(x, MARGIN = 2, FUN = function(x){stat.desc(x, basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)})})
fwd.df.summary.stats <- lapply(fwd.df.summary.stats, FUN = function(x){round(x[-(1:3),],5)})
knitr::kable(round(fwd.df,4), caption="Resumen simulación walk-forward - percentil 50") %>%
  kable_styling(latex_options = "HOLD_position")
# ANEXO-01: ----
# Ejemplo de aplicacion VaR desde el benchmark al portafolio particular ---
df.returns <- gobix[,'daily.vol']
colnames(df.returns) <- "gobix.rets"

set.seed(1)
sim.portfolio_vol.index <- sample(x = df.returns$gobix.rets, size = nrow(df.returns), replace = FALSE) # simulate the base level
sim.portfolio_vol.random <- rnorm(n=nrow(df.returns), mean=0, sd=0.002) # generate the noise
sim.portfolio_vol <- sim.portfolio_vol.index+(sim.portfolio_vol.random) # base + noise
sim.portfolio_vol <- sim.portfolio_vol*0.70 # scale the vol down

df.returns$portfolio.rets <- sim.portfolio_vol
beta.estimation <- lm(coredata(df.returns$portfolio.rets)~coredata(df.returns$gobix.rets))

regression.result <- summary(beta.estimation)

par(mar=c(4,4,4,4), oma=c(0,0,0,0))
plot(coredata(df.returns$gobix.rets), coredata(df.returns$portfolio.rets),
     xlim=c(-0.03,0.03), ylim=c(-0.03,0.03),
     main="Retornos diarios\n portafolio simulado vs. benchmark",
     ylab="Portafolio", xlab="GOBIXDR", col="blue", pch=16, cex=0.60, cex.main=0.70, cex.sub=0.75, cex.lab=0.70, cex.axis=0.75)
abline(h=mean(df.returns$portfolio.rets), v=mean(df.returns$gobix.rets), col="grey40", lty=2)
abline(beta.estimation,  col="red", lwd=1.75)
box(col="grey")

# adjust
sim.resume.adjusted <- apply(sim.resume, MARGIN = 2, FUN = function(x){x*regression.result$coefficients[2,1]*(sd(coredata(df.returns$portfolio.rets)) / sd(coredata(df.returns$gobix.rets)))})

sim.resume.adjusted <- round(sim.resume.adjusted,2)
colnames(sim.resume.adjusted) <- c("BIII.estres","R(i=1)","R(i=2)","R(i=3)")
colnames(sim.resume) <- c("BIII.estres","R(i=1)","R(i=2)","R(i=3)")

# tabla benchmark 
knitr::kable(sim.resume, caption="Simulación estrés benchmark") %>% kable_styling(latex_options = "hold", position = "center")

# tabla ajuste
knitr::kable(sim.resume.adjusted, caption="Estimación estrés portafolio ajustado") %>% kable_styling(latex_options = "hold", position = "center")

# ANEXO-02: graph ----
par(mfrow=c(2,2), mar=c(4,4,4,4))
regime.lower.sim <- descdist_custom(as.numeric(regime.lower[,'daily.vol',drop=FALSE]), 
                                    discrete=FALSE, 
                                    boot=N_BOOT, 
                                    title = "TPM: Reducción",
                                    subtitle = "R(i=3)", 
                                    label.align = 0.5)

regime.flat.sim <- descdist_custom(as.numeric(regime.flat[,'daily.vol',drop=FALSE]), 
                                   discrete=FALSE, 
                                   boot=N_BOOT, 
                                   title = "TPM: Sin cambios",
                                   subtitle = "R(i=2)",
                                   label.align = 0.5)

regime.higher.sim <- descdist_custom(as.numeric(regime.higher[,'daily.vol',drop=FALSE]), 
                                     discrete=FALSE, 
                                     boot=N_BOOT, 
                                     title = "TPM: Incremento",
                                     subtitle = "R(i=1)",
                                     label.align = 0.5)

regime.flat.sim <- descdist_custom(as.numeric(baselIII.12m.stress[,'daily.vol',drop=FALSE]), 
                                   discrete=FALSE, 
                                   boot=N_BOOT, 
                                   title = "Basel III",
                                   subtitle = "BIII",
                                   label.align = 0.5)

par(mfrow=c(1,1))
# ANEXO-02: tablas ----
# INCERTIDUBRE PARAMETROS: RESUME EN UNA TABLA PARA CADA PARAMETRO (mu,sigma,kurt,skew)
sim_lower.regime.mean <- apply(sim.lower[[1]][[1]], MARGIN = 2, FUN = function(x){mean(diff(log(x)))*sqrt(250)})
sim_flat.regime.mean <- apply(sim.flat[[1]][[1]], MARGIN = 2, FUN = function(x){mean(diff(log(x)))*sqrt(250)})
sim_higher.regime.mean <- apply(sim.higher[[1]][[1]], MARGIN = 2, FUN = function(x){mean(diff(log(x)))*sqrt(250)})
sim_stress.regime.mean <- apply(sim.stress[[1]][[1]], MARGIN = 2, FUN = function(x){mean(diff(log(x)))*sqrt(250)})

sim_lower.regime.sd <- apply(sim.lower[[1]][[1]], MARGIN = 2, FUN = function(x){sd(diff(log(x)))*sqrt(250)})
sim_flat.regime.sd <- apply(sim.flat[[1]][[1]], MARGIN = 2, FUN = function(x){sd(diff(log(x)))*sqrt(250)})
sim_higher.regime.sd <- apply(sim.higher[[1]][[1]], MARGIN = 2, FUN = function(x){sd(diff(log(x)))*sqrt(250)})
sim_stress.regime.sd <- apply(sim.stress[[1]][[1]], MARGIN = 2, FUN = function(x){sd(diff(log(x)))*sqrt(250)})

sim_lower.regime.skew <- apply(sim.lower[[1]][[1]], MARGIN = 2, FUN = function(x){moments::skewness(diff(log(x)))})
sim_flat.regime.skew <- apply(sim.flat[[1]][[1]], MARGIN = 2, FUN = function(x){moments::skewness(diff(log(x)))})
sim_higher.regime.skew <- apply(sim.higher[[1]][[1]], MARGIN = 2, FUN = function(x){moments::skewness(diff(log(x)))})
sim_stress.regime.skew <- apply(sim.stress[[1]][[1]], MARGIN = 2, FUN = function(x){moments::skewness(diff(log(x)))})

sim_lower.regime.kurt <- apply(sim.lower[[1]][[1]], MARGIN = 2, FUN = function(x){moments::kurtosis(diff(log(x)))})
sim_flat.regime.kurt <- apply(sim.flat[[1]][[1]], MARGIN = 2, FUN = function(x){moments::kurtosis(diff(log(x)))})
sim_higher.regime.kurt <- apply(sim.higher[[1]][[1]], MARGIN = 2, FUN = function(x){moments::kurtosis(diff(log(x)))})
sim_stress.regime.kurt <- apply(sim.stress[[1]][[1]], MARGIN = 2, FUN = function(x){moments::kurtosis(diff(log(x)))})

# compute the distribution of parameter uncertainty
sim.regime.kurt <- list(sim_lower.regime.kurt,
                        sim_flat.regime.kurt,
                        sim_higher.regime.kurt,
                        sim_stress.regime.kurt)

sim.regime.skew <- list(sim_lower.regime.skew,
                        sim_flat.regime.skew,
                        sim_higher.regime.skew,
                        sim_stress.regime.skew)

sim.regime.sd <- list(sim_lower.regime.sd,
                      sim_flat.regime.sd,
                      sim_higher.regime.sd,
                      sim_stress.regime.sd)

sim.regime.mean <- list(sim_lower.regime.mean,
                        sim_flat.regime.mean,
                        sim_higher.regime.mean,
                        sim_stress.regime.mean)

fourth_moment.sim <- sapply(sim.regime.kurt, FUN = function(x){compute_stats(object = x)})
  colnames(fourth_moment.sim) <- c("R(i=1)","R(i=2)","R(i=3)","R(i=4)")

third_moment.sim <- sapply(sim.regime.skew, FUN = function(x){compute_stats(object = x)})
  colnames(third_moment.sim) <- c("R(i=1)","R(i=2)","R(i=3)","R(i=4)")

second_moment.sim <- sapply(sim.regime.sd, FUN = function(x){compute_stats(object = x)})
  colnames(second_moment.sim) <- c("R(i=1)","R(i=2)","R(i=3)","R(i=4)")

first_moment.sim <- sapply(sim.regime.mean, FUN = function(x){compute_stats(object = x)})
  colnames(first_moment.sim) <- c("R(i=1)","R(i=2)","R(i=3)","R(i=4)")

knitr::kable(first_moment.sim, caption="Distribución: media")

knitr::kable(second_moment.sim, caption="Distribución: std dev")

knitr::kable(third_moment.sim, caption="Distribución: sesgo")

knitr::kable(fourth_moment.sim, caption="Distribución: curtosis")

# WALK FORWARD
fwd.df.summary.stats <- lapply(fwd_oos.all, FUN = function(x){apply(x[,2:ncol(x)], MARGIN = 2, FUN = function(x){stat.desc(x, basic=TRUE, desc=TRUE, norm=TRUE, p=0.95)})})
fwd.df.summary.stats <- lapply(fwd.df.summary.stats, FUN = function(x){round(x[-c(1:3,7,10,11,12,14,20),],5)})

colnames(fwd.df.summary.stats[[1]]) <- c('BIII','R(i=1)','R(i=2)', 'R(i=3)')
colnames(fwd.df.summary.stats[[2]]) <- c('BIII','R(i=1)','R(i=2)', 'R(i=3)')
colnames(fwd.df.summary.stats[[3]]) <- c('BIII','R(i=1)','R(i=2)', 'R(i=3)')

knitr::kable(fwd.df.summary.stats[[1]], caption="Estadísticas sumarizadas - VaR.95") %>%
  kable_styling(latex_options = "HOLD_position")

knitr::kable(fwd.df.summary.stats[[2]], caption="Estadísticas sumarizadas - VaR.975") %>%
  kable_styling(latex_options = "HOLD_position")

knitr::kable(fwd.df.summary.stats[[3]], caption="Estadísticas sumarizadas - VaR.99") %>%
  kable_styling(latex_options = "HOLD_position")
# ANEXO-02: LPM-UPM ----
lower.quantiles <- as.data.frame(sim.lower[[1]][[1]][,c('Q010','Q025','Q050','Q950','Q975','Q990')])
lower.quantiles <- apply(lower.quantiles, MARGIN = 2, FUN = function(x){x-100})

flat.quantiles <- as.data.frame(sim.flat[[1]][[1]][,c('Q010','Q025','Q050','Q950','Q975','Q990')])
flat.quantiles <- apply(flat.quantiles, MARGIN = 2, FUN = function(x){x-100})

higher.quantiles <- as.data.frame(sim.higher[[1]][[1]][,c('Q010','Q025','Q050','Q950','Q975','Q990')])
higher.quantiles <- apply(higher.quantiles, MARGIN = 2, FUN = function(x){x-100})

par(mfrow=c(3,1), mar=c(4,4,4,4), oma=c(0,0,0,0))
plot(lower.quantiles[,'Q990'], lower.quantiles[,'Q010'], 
     ylim = rev(range(lower.quantiles[,'Q010'])),
     type="l", lwd=2, col="green4", xlab="nivel cola derecha (%)", ylab="nivel cola izquierda (%)", main="R(i=3)",
     cex=1.25, cex.axis=1.25, cex.lab=1, cex.main=1.25)
lines(lower.quantiles[,'Q975'], lower.quantiles[,'Q025'], lwd=2, col="green2")
lines(lower.quantiles[,'Q950'], lower.quantiles[,'Q050'], lwd=2, col="green")
legend("topleft", legend = c("[P.010,P.990]","[P.025,P.975]","[P.050,P.950]"), lwd=2,
       col = c("green4","green2","green"), text.col = c("green4","green2","green"), bty = "n", cex=1, y.intersp = 1.25)
box(col="grey")

plot(flat.quantiles[,'Q990'], flat.quantiles[,'Q010'],type="l", 
     ylim = rev(range(flat.quantiles[,'Q010'])),
     col="blue4", lwd=2, xlab="nivel cola derecha (%)", ylab="nivel cola izquierda (%)", main="R(i=2)",
     cex=1.25, cex.axis=1.25, cex.lab=1, cex.main=1.25)
lines(flat.quantiles[,'Q975'], flat.quantiles[,'Q025'], lwd=2, col="blue2")
lines(flat.quantiles[,'Q950'], flat.quantiles[,'Q050'], lwd=2, col="blue")
box(col="grey")

plot(higher.quantiles[,'Q990'], higher.quantiles[,'Q010'], type="l", 
     ylim = rev(range(higher.quantiles[,'Q010'])),
     col="red4", lwd=2, xlab="nivel cola derecha (%)", ylab="nivel cola izquierda (%)", main="R(i=1)",
     cex=1.25, cex.axis=1.25, cex.lab=1, cex.main=1.25)
lines(higher.quantiles[,'Q975'], higher.quantiles[,'Q025'], lwd=2, col="red2")
lines(higher.quantiles[,'Q950'],higher.quantiles[,'Q050'], lwd=2, col="red")
box(col="grey")
# ANEXO-02: evolucion parametros ----

par(mfrow=c(3,1), mar=c(3,4,3,4), oma=c(0,0,0,0))
plot(variance.df$lower~variance.df$date, type="l", ylim=c(0,0.20), ylab="std dev anualizada", xlab="", col="green4", main="Desviación estándar")
  lines(variance.df$sin_cambios~variance.df$date, col="blue", lwd=1.5, lty=1)
  lines(variance.df$higher~variance.df$date, col="red", lwd=1.5, lty=1)
  lines(variance.df$BIII.stress~variance.df$date, col="red4", lwd=1.5, lty=1)
  legend("topleft", legend = c("reducción","sin.cambios","incremento","BIII.stress"), lwd=2, lty=1, col=c("green4","blue","red","red4"), bty="n", y.intersp = 0.8)
  box(col="grey")
plot(skew.df$lower~skew.df$date, type="l", ylim=c(-2,2), ylab="skew (sesgo)", xlab="", col="green4", main="Sesgo")
  lines(skew.df$sin_cambios~skew.df$date, col="blue", lwd=1.5, lty=1)
  lines(skew.df$higher~skew.df$date, col="red", lwd=1.5, lty=1)
  lines(skew.df$BIII.stress~skew.df$date, col="red4", lwd=1.5, lty=1)
  box(col="grey")
plot(kurt.df$lower~kurt.df$date, type="l", ylim=c(0,20), ylab="kurtosis (curtosis)", xlab="", col="green4", main="Curtosis")
  lines(kurt.df$sin_cambios~kurt.df$date, col="blue", lwd=1.5, lty=1)
  lines(kurt.df$higher~kurt.df$date, col="red", lwd=1.5, lty=1)
  lines(kurt.df$BIII.stress~kurt.df$date, col="red4", lwd=1.5, lty=1)
  box(col="grey")
par(mfrow=c(1,1), mar=c(4,4,4,4))
