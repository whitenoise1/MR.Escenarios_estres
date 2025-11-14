# MR.Escenarios_estres
Código fuente investigacion: ["Escenarios de estrés en portafolios de inversiones: Un enfoque empírico"](https://sb.gob.do/publicaciones/publicaciones-tecnicas/escenarios-de-estres-en-portafolios-de-inversiones-un-enfoque-empirico/).

Aquí, se presenta un caso de aplicación de pruebas de estrés en portafolios de bonos gubernamentales, utilizando una serie univariada del índice precio-retorno como benchmark de referencia. El enfoque adoptado se enmarca dentro de la taxonomía de métodos no paramétricos. La construcción de escenarios se basó en fundamentos teóricos, a partir de la información públicamente disponible y siguiendo principios estrictamente empíricos. Se aplicó la técnica de meta-labeling para identificar períodos con características similares y definir los regímenes de mercado. A partir de las características identificadas, se procedió a generar choques condicionales a cada régimen mediante la Simulación de Monte Carlo (SMC). 

* Regimenes de TPM.xlsx: Contiene la serie de tiempo de la Tasa de Pólitica Monetaria (TPM). Fuente: Banco Central RD.
* empirical_VaR_rollingOOS_p95.csv: Es el resultado de la simulación (VaR 95%) walk forward (out-of-sample) presentado.
* empirical_VaR_rollingOOS_p975.csv: Es el resultado de la simulación (VaR 97.5%) walk forward (out-of-sample) presentado.
* empirical_VaR_rollingOOS_p99.csv: Es el resultado de la simulación (VaR 99%) walk forward (out-of-sample) presentado.
* MR_EscenariosEstres.R: Código fuente completo.
* EscenariosEstres_empirico.RData: Contiene los resultados completos del proceso.

# Librerias requeridas
```{r}
require(quantmod)
require(PerformanceAnalytics)
require(roll)
require(pastecs)
require(stargazer)
require(knitr)
require(fAssets)
require(readxl)
require(dplyr)
require(lubridate)
require(bizdays)
require(PearsonDS)
require(fitdistrplus)
require(kableExtra)
require(cubature)
```
