# refactor/ — versión modular de MR_EscenariosEstres.R

Esta carpeta contiene una versión modular y limpia del script `MR_EscenariosEstres.R`:

- **`functions.R`** — todas las funciones auxiliares, copiadas verbatim del script original.
- **`main.R`** — el pipeline de análisis completo, en el mismo orden que el original, con un bloque de `# Parámetros ----` al inicio y `source("refactor/functions.R")`. La simulación walk-forward (costosa) está desactivada por defecto (`RUN_WALK_FORWARD <- FALSE`); en su lugar se leen los archivos precomputados `empirical_VaR_rollingOOS_p*.csv`.

El script original `MR_EscenariosEstres.R` (en la raíz del repositorio) sigue siendo la versión canónica que replica el paper.

## Uso

Ejecutar `main.R` con el directorio de trabajo en la **raíz del repositorio** (no dentro de `refactor/`), para que las lecturas relativas de `Regimenes de TPM.xlsx` y de los CSV precomputados funcionen sin cambios:

```r
# setwd("<raíz del repositorio>")
source("refactor/main.R")
```

Los resultados son idénticos a los del script original por construcción (refactorización que preserva la estructura: solo se extrajeron funciones, se centralizaron parámetros con sustituciones 1:1 y se eliminó código muerto comentado).
