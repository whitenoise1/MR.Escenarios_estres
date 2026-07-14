# MR.Escenarios_estres

Código fuente que replica la investigación: ["Escenarios de estrés en portafolios de inversiones: Un enfoque empírico"](https://sb.gob.do/publicaciones/publicaciones-tecnicas/escenarios-de-estres-en-portafolios-de-inversiones-un-enfoque-empirico/), Superintendencia de Bancos de la República Dominicana (2025).

**Palabras clave:** portafolio de inversiones, riesgos financieros, prueba de estrés, simulación de Monte-Carlo.
**Clasificación JEL:** G00, G10, G17, G20, G21.

## Resumen

En este repositorio se presenta un caso de aplicación de pruebas de estrés en portafolios de bonos gubernamentales, utilizando una serie univariada del índice precio-retorno (GOBIXDR) como benchmark de referencia. El enfoque adoptado se enmarca dentro de la taxonomía de métodos no paramétricos. La construcción de escenarios se basó en fundamentos teóricos, a partir de la información públicamente disponible y siguiendo principios estrictamente empíricos. Se aplicó la técnica de meta-labeling para identificar períodos con características similares y definir los regímenes de mercado. A partir de las características identificadas, se procedió a generar choques condicionales a cada régimen mediante la Simulación de Monte Carlo (SMC).

## Contenido del repositorio

| Archivo | Descripción |
|---|---|
| `MR_EscenariosEstres.R` | Código fuente completo de la investigación (punto de entrada principal). |
| `Regimenes de TPM.xlsx` | Serie de tiempo de la Tasa de Política Monetaria (TPM). **Insumo requerido** por el script. Fuente: Banco Central de la República Dominicana. |
| `empirical_VaR_rollingOOS_p95.csv` | Resultado precalculado de la simulación walk-forward (out-of-sample), VaR 95%. |
| `empirical_VaR_rollingOOS_p975.csv` | Resultado precalculado de la simulación walk-forward (out-of-sample), VaR 97.5%. |
| `empirical_VaR_rollingOOS_p99.csv` | Resultado precalculado de la simulación walk-forward (out-of-sample), VaR 99%. |
| `duracion_GOBIX_2015-2025.csv` | Serie de duración del GOBIXDR, insumo del análisis de duración presentado en el paper (no es leído por el script principal). |
| `EscenariosEstres_empirico.RData` | Workspace de R con las funciones, el análisis y los resultados completos del proceso. |
| `paper_EscenariosEstres_final.pdf` | Versión final del documento publicado. |
| `refactor/` | Versión modular y limpia del script (`functions.R` + `main.R`), con resultados idénticos por construcción. Ver `refactor/README.md`. |

## Cómo ejecutar

1. **Requisitos:** R ≥ 4.0 y conexión a internet — el script descarga la serie del GOBIXDR directamente desde la BVRD (`https://www.bvrd.com.do/indice/Data/GobixDataIRP.csv`).
2. **Instalar los paquetes requeridos:**

```r
install.packages(c("quantmod", "PerformanceAnalytics", "roll", "pastecs",
                   "stargazer", "knitr", "fAssets", "readxl", "dplyr",
                   "lubridate", "bizdays", "PearsonDS", "fitdistrplus",
                   "kableExtra", "cubature"))
```

3. **Fijar el directorio de trabajo en la raíz del repositorio** (el script lee `Regimenes de TPM.xlsx` y los `.csv` por ruta relativa):

```r
setwd("ruta/al/repositorio/MR.Escenarios_estres")
source("MR_EscenariosEstres.R")
```

También puede ejecutarse por secciones desde RStudio (el script está organizado con encabezados `# Sección ----`).

**Nota sobre tiempos de ejecución:** la simulación walk-forward (out-of-sample) es computacionalmente costosa; por esa razón sus resultados se distribuyen precalculados en los tres archivos `empirical_VaR_rollingOOS_p*.csv`, que el script lee directamente. Alternativamente, `EscenariosEstres_empirico.RData` puede cargarse con `load()` para explorar todos los objetos y resultados sin re-ejecutar el análisis.

## Principales funciones

* **`empirical.simulation()`**: ejecuta la simulación Pearson a partir del *fit* a la distribución empírica observada en la muestra. Retorna el resultado de la distribución completa.
* **`empirical.simulation.simple()`**: versión abreviada que ejecuta el mismo proceso, pero retorna únicamente el *point-estimate* del resultado.

## Licencia

Este proyecto se distribuye bajo la licencia [GNU GPL v3](LICENSE).
