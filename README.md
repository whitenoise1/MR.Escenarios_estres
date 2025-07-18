# MR.Escenarios_estres
Código fuente de la investigacion: "Escenarios de estrés en portafolios de inversiones: Un enfoque empírico".

Se presenta un proceso de construcción de escenarios de estrés siguiendo principios puramente empíricos. La aplicación se desarrolla utilizando la técnica de meta-labeling para identificar períodos con características similares en la serie de tasas de interés y, posteriormente, analizar el comportamiento del índice benchmark durante esos períodos. Los datos se aproximan a la función de densidad de probabilidad empírica mediante la solución cerrada de la familia de distribuciones Pearson Tipo IV. A partir de esto, se generan los escenarios por medio de Simulación de Monte Carlo (SMC) utilizando datos sintéticos y se interpretan los resultados.

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
