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
require(lubridate)
require(PearsonDS)
require(lubridate)
require(bizdays)
require(fitdistrplus)
require(lubridate)
require(bizdays)
require(kableExtra)
require(cubature)

# Funciones auxiliares ----
descdist_custom <- function(data, discrete = FALSE, boot = NULL, method = "unbiased", 
                            graph = TRUE, print = TRUE, obs.col = "red", obs.pch = 16, 
                            boot.col = "grey40", title, subtitle, label.align) 
{
  if (missing(data) || !is.vector(data, mode = "numeric")) 
    stop("data must be a numeric vector")
  if (length(data) < 4) 
    stop("data must be a numeric vector containing at least four values")
  moment <- function(data, k) {
    m1 <- mean(data)
    return(sum((data - m1)^k)/length(data))
  }
  if (method == "unbiased") {
    skewness <- function(data) {
      sd <- sqrt(moment(data, 2))
      n <- length(data)
      gamma1 <- moment(data, 3)/sd^3
      unbiased.skewness <- sqrt(n * (n - 1)) * gamma1/(n - 
                                                         2)
      return(unbiased.skewness)
    }
    kurtosis <- function(data) {
      n <- length(data)
      var <- moment(data, 2)
      gamma2 <- moment(data, 4)/var^2
      unbiased.kurtosis <- (n - 1)/((n - 2) * (n - 3)) * 
        ((n + 1) * gamma2 - 3 * (n - 1)) + 3
      return(unbiased.kurtosis)
    }
    standdev <- function(data) {
      sd(data)
    }
  }
  else if (method == "sample") {
    skewness <- function(data) {
      sd <- sqrt(moment(data, 2))
      return(moment(data, 3)/sd^3)
    }
    kurtosis <- function(data) {
      var <- moment(data, 2)
      return(moment(data, 4)/var^2)
    }
    standdev <- function(data) {
      sqrt(moment(data, 2))
    }
  }
  else stop("The only possible value for the argument method are 'unbiased' or 'sample'")
  res <- list(min = min(data), max = max(data), median = median(data), 
              mean = mean(data), sd = standdev(data), skewness = skewness(data), 
              kurtosis = kurtosis(data), method = method)
  skewdata <- res$skewness
  kurtdata <- res$kurtosis
  if (graph) {
    if (!is.null(boot)) {
      if (!is.numeric(boot) || boot < 10) {
        stop("boot must be NULL or a integer above 10")
      }
      n <- length(data)
      databoot <- matrix(sample(data, size = n * boot, 
                                replace = TRUE), nrow = n, ncol = boot)
      s2boot <- sapply(1:boot, function(iter) skewness(databoot[, 
                                                                iter])^2)
      kurtboot <- sapply(1:boot, function(iter) kurtosis(databoot[, 
                                                                  iter]))
      kurtmax <- max(10, ceiling(max(kurtboot)))
      xmax <- max(4, ceiling(max(s2boot)))
    }
    else {
      kurtmax <- max(10, ceiling(kurtdata))
      xmax <- max(4, ceiling(skewdata^2))
    }
    ymax <- kurtmax - 1
    plot(skewdata^2, kurtmax - kurtdata, pch = obs.pch, xlim = c(0, xmax), ylim = c(0, ymax), yaxt = "n", xlab = "square of skewness", 
         ylab = "kurtosis", main = title, col=scales::alpha(boot.col, 0.25))
    mtext(subtitle, side = 3, cex=0.90, adj = 0.5, padj = -0.5)
    box(col="grey")
    yax <- as.character(kurtmax - 0:ymax)
    axis(side = 2, at = 0:ymax, labels = yax, col="grey")
    axis(side = 1, col="grey")
    if (!discrete) {
      p <- exp(-100)
      lq <- seq(-100, 100, 0.1)
      q <- exp(lq)
      s2a <- (4 * (q - p)^2 * (p + q + 1))/((p + q + 2)^2 * 
                                              p * q)
      ya <- kurtmax - (3 * (p + q + 1) * (p * q * (p + 
                                                     q - 6) + 2 * (p + q)^2)/(p * q * (p + q + 2) * 
                                                                                (p + q + 3)))
      p <- exp(100)
      lq <- seq(-100, 100, 0.1)
      q <- exp(lq)
      s2b <- (4 * (q - p)^2 * (p + q + 1))/((p + q + 2)^2 * 
                                              p * q)
      yb <- kurtmax - (3 * (p + q + 1) * (p * q * (p + 
                                                     q - 6) + 2 * (p + q)^2)/(p * q * (p + q + 2) * 
                                                                                (p + q + 3)))
      s2 <- c(s2a, s2b)
      y <- c(ya, yb)
      polygon(s2, y, col = "lightgrey", border = "lightgrey")
      lshape <- seq(-100, 100, 0.1)
      shape <- exp(lshape)
      s2 <- 4/shape
      y <- kurtmax - (3 + 6/shape)
      lines(s2[s2 <= xmax], y[s2 <= xmax], lty = 2)
      lshape <- seq(-100, 100, 0.1)
      shape <- exp(lshape)
      es2 <- exp(shape^2)
      s2 <- (es2 + 2)^2 * (es2 - 1)
      y <- kurtmax - (es2^4 + 2 * es2^3 + 3 * es2^2 - 3)
      lines(s2[s2 <= xmax], y[s2 <= xmax], lty = 3)
      legend(xmax * label.align, ymax * 0.90, legend = "Distributions",
             bty = "n", cex = 0.8)
      legend(xmax * label.align, ymax * 0.86, pch = 8, legend = "normal",
             bty = "n", cex = 0.8, col = "blue")
      legend(xmax * label.align, ymax * 0.82, pch = 2, legend = "uniform",
             bty = "n", cex = 0.8, col = "blue")
      legend(xmax * label.align, ymax * 0.78, pch = 7, legend = "exponential",
             bty = "n", cex = 0.8, col = "blue")
      legend(xmax * label.align, ymax * 0.74, pch = 3, legend = "logistic",
             bty = "n", cex = 0.8, col = "blue")
      legend(xmax * label.align, ymax * 0.62, fill = "grey80",
             legend = "beta", bty = "n", cex = 0.8)
      legend(xmax * label.align, ymax * 0.70, lty = 3, legend = "lognormal",
             bty = "n", cex = 0.8)
      legend(xmax * label.align, ymax * 0.66, lty = 2, legend = "gamma",
             bty = "n", cex = 0.8)
      legend(xmax * label.align, ymax * 1.02, pch = obs.pch, legend = "data", 
             bty = "n", cex = 0.8, pt.cex = 1.2, col = obs.col)
      if (!is.null(boot)) {
        legend(xmax * label.align, ymax * 0.98, pch = 1, legend = "bootstrap", 
               bty = "n", cex = 0.8, col = boot.col)
      }
    }
    else {
      p <- exp(-10)
      lr <- seq(-100, 100, 0.1)
      r <- exp(lr)
      s2a <- (2 - p)^2/(r * (1 - p))
      ya <- kurtmax - (3 + 6/r + p^2/(r * (1 - p)))
      p <- 1 - exp(-10)
      lr <- seq(100, -100, -0.1)
      r <- exp(lr)
      s2b <- (2 - p)^2/(r * (1 - p))
      yb <- kurtmax - (3 + 6/r + p^2/(r * (1 - p)))
      s2 <- c(s2a, s2b)
      y <- c(ya, yb)
      polygon(s2, y, col = "grey80", border = "grey80")
      legend(xmax * 0.55, ymax * 1.03, legend = "Theoretical", 
             bty = "n", cex = 0.8)
      legend(xmax * 0.6, ymax * 0.98, pch = 8, legend = "normal", 
             bty = "n", cex = 0.8)
      legend(xmax * 0.6, ymax * 0.94, fill = "grey80", 
             legend = "negative binomial", bty = "n", cex = 0.8)
      legend(xmax * 0.6, ymax * 0.9, lty = 2, legend = "Poisson", 
             bty = "n", cex = 0.8)
      legend(xmax * 0.55, ymax * 0.85, legend = "Empirical", 
             bty = "n", cex = 0.8)
      legend(xmax * 0.6, ymax * 0.81, pch = obs.pch, legend = "data", 
             bty = "n", cex = 0.8, pt.cex = 1.2, col = obs.col)
      if (!is.null(boot)) {
        legend(xmax * 0.6, ymax * 0.77, pch = 1, legend = "bootstrap", 
               bty = "n", cex = 0.8, col = boot.col)
      }
      llambda <- seq(-100, 100, 0.1)
      lambda <- exp(llambda)
      s2 <- 1/lambda
      y <- kurtmax - (3 + 1/lambda)
      lines(s2[s2 <= xmax], y[s2 <= xmax], lty = 2)
    }
    if (!is.null(boot)) {
      points(s2boot, kurtmax - kurtboot, pch = 1, col = scales::alpha(boot.col, 0.25), 
             cex = 0.5)
    }
    points(skewness(data)^2, kurtmax - kurtosis(data), pch = obs.pch, 
           cex = 2, col = obs.col)
    points(0, kurtmax - 3, pch = 8, cex = 1.5, lwd = 1.25, col = "blue")
    if (!discrete) {
      points(0, kurtmax - 9/5, pch = 2, cex = 1.5, lwd = 1.25, col = "blue")
      points(2^2, kurtmax - 9, pch = 7, cex = 1.5, lwd = 1.25, col = "blue")
      points(0, kurtmax - 4.2, pch = 3, cex = 1.5, lwd = 1.25, col = "blue")
    }
  }
  if (!print) 
    invisible(structure(res, class = "descdist"))
  else structure(res, class = "descdist")
}

empirical.simulation <- function(data, 
                                 n.sim, 
                                 t, 
                                 oos.series,
                                 mean, 
                                 variance, 
                                 skew, 
                                 kurt, 
                                 main){
  # data has to be XTS [OHLC] object
  if(is.xts(data) == FALSE){
    print("Error: data field is required in XTS format [OHLC].")
  } else if(is.xts(data) == TRUE){
    
    require(lubridate)
    require(bizdays)
    require(xts)
    require(quantmod)
    require(PearsonDS)
    
    set.seed(1)
    days.ahead = t
    n.sim <- n.sim
    emprical.sim <- matrix(nrow = days.ahead, ncol=n.sim)
    for(i in 1:n.sim){
      df <- rpearson(n=days.ahead, moments=c(mean=mean, variance=variance, skewness=skew, kurtosis=kurt))
      emprical.sim[,i] <- df
    }
    
    sim <- as.data.frame(emprical.sim)
    
    sim$Q050 <- NA
    sim$Q025 <- NA
    sim$Q010 <- NA
    sim$Q950 <- NA
    sim$Q975 <- NA
    sim$Q990 <- NA
    
    sim$Q010 <- apply(sim[,2:(ncol(sim)-6)], FUN = function(x){quantile(na.omit(x),0.010)}, MARGIN = 1)
    sim$Q025 <- apply(sim[,2:(ncol(sim)-6)], FUN = function(x){quantile(na.omit(x),0.025)}, MARGIN = 1)
    sim$Q050 <- apply(sim[,2:(ncol(sim)-6)], FUN = function(x){quantile(na.omit(x),0.050)}, MARGIN = 1)
    sim$Q950 <- apply(sim[,2:(ncol(sim)-6)], FUN = function(x){quantile(na.omit(x),0.950)}, MARGIN = 1)
    sim$Q975 <- apply(sim[,2:(ncol(sim)-6)], FUN = function(x){quantile(na.omit(x),0.975)}, MARGIN = 1)
    sim$Q990 <- apply(sim[,2:(ncol(sim)-6)], FUN = function(x){quantile(na.omit(x),0.990)}, MARGIN = 1)
    
    sim.std <- sd(unlist(sim[,2:(ncol(sim)-12)]))
    
    # adjust for business days
    business.calendar <- create.calendar('my_calendar', weekdays = c('saturday','sunday'))
    sim.dates <- bizdays::offset(ymd(last(index(data))+1), 1:days.ahead, cal = business.calendar)
    
    df <- cbind(as.data.frame(sim.dates),sim)
    df <- as.xts(df[,-1], order.by=ymd(df$sim.dates))
    df <- as.numeric(tail(Cl(data),1)) * cumprod(1 + df)
    df$empty <- NA
    
    df.combined <- rbind(Cl(data),df$empty)
    df.combined <- merge(df.combined,df[,1:(ncol(df)-1)]) # do not merge the empty column
    df.combined[1,2:ncol(df.combined)] <- df.combined[1,1]
    df.combined <- df.combined[,-1]
    
    df.combined$Q010 <- apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){quantile(na.omit(x),0.010)}, MARGIN = 1)
    df.combined$Q025 <- apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){quantile(na.omit(x),0.025)}, MARGIN = 1)
    df.combined$Q050 <- apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){quantile(na.omit(x),0.050)}, MARGIN = 1)
    df.combined$Q950 <- apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){quantile(na.omit(x),0.950)}, MARGIN = 1)
    df.combined$Q975 <- apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){quantile(na.omit(x),0.975)}, MARGIN = 1)
    df.combined$Q990 <- apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){quantile(na.omit(x),0.990)}, MARGIN = 1)
    df.combined$mean <- as.numeric(apply(df.combined[,2:(ncol(df.combined)-7)], FUN = function(x){mean(na.omit(x))}, MARGIN = 1))
    
    last.day <- index(tail(Cl(data),1))
    last.price <- as.numeric(tail(Cl(data),1))
    
    z <- as.data.frame(merge(Cl(data),df.combined[,2:ncol(df.combined)]))
    
    # below here, it breaks when the dates are not aligned
    quantInv <- function(distr, value) ecdf(distr)(value)
    level.prob <- quantInv(na.omit(as.numeric(df.combined[last.day])),last.price)
    quantile.series <- apply(z[,2:ncol(df.combined)], MARGIN = 1, FUN = function(x) { quantInv(x,z[,1]) })
    quantile.series <- diag(quantile.series)
    
    results <- list(df.combined, list(level.prob,last.day,last.price), quantile.series, data, sim.std)
    
    
    par(mfrow=c(1,1), mar=c(5,5,5,5))
    sim.plot <- plot(df.combined[,1], ylim=c(min(na.omit(df.combined[,2:ncol(df.combined)])),max(na.omit(df.combined[,2:ncol(df.combined)]))*1.1), col="white", main=main, grid.col=NA)
    sim.plot <- lines(df.combined[,2:(ncol(df.combined)-6)], col=scales::alpha("grey",0.3), on=1, lty=2, lwd=0.75)
    sim.plot <- lines(df.combined[,'Q010'], col="red", on=1, lty=1, lwd=2)
    sim.plot <- lines(df.combined[,'Q025'], col="red", on=1, lty=1, lwd=1)
    sim.plot <- lines(df.combined[,'Q050'], col="red", on=1, lty=2, lwd=0.5)
    sim.plot <- lines(df.combined[,'mean'], col="grey40", on=1, lty=2, lwd=1.75)
    sim.plot <- lines(df.combined[,'Q950'], col="blue", on=1, lty=2, lwd=0.5)
    sim.plot <- lines(df.combined[,'Q975'], col="blue", on=1, lty=1, lwd=1)
    sim.plot <- lines(df.combined[,'Q990'], col="blue", on=1, lty=1, lwd=2)
    #sim.plot <- lines(Cl(oos), col="black", on=1, lty=1, lwd=1.5)
    sim.plot <- points(tail(df.combined[,'mean'],1), col="red", pch=16)
    
    output <- list(results,
                   sim.plot)
    
    #return(output)
    
    if(is.null(oos.series) == TRUE){
      
      return(output)
      
    } else if(is.xts(oos.series) == TRUE){
      
      par(mfrow=c(1,1), mar=c(5,5,5,5))
      
      sim.plot <- plot(df.combined[,1], ylim=c(min(na.omit(df.combined[,2:ncol(df.combined)])),max(na.omit(df.combined[,2:ncol(df.combined)]))*1.1), col="white", main=main, grid.col=NA)
      sim.plot <- lines(df.combined[,2:(ncol(df.combined)-6)], col=scales::alpha("grey",0.3), on=1, lty=2, lwd=0.75)
      sim.plot <- lines(df.combined[,'Q010'], col="red", on=1, lty=1, lwd=2)
      sim.plot <- lines(df.combined[,'Q025'], col="red", on=1, lty=1, lwd=1)
      sim.plot <- lines(df.combined[,'Q050'], col="red", on=1, lty=2, lwd=0.5)
      sim.plot <- lines(df.combined[,'mean'], col="grey40", on=1, lty=2, lwd=1.75)
      sim.plot <- lines(df.combined[,'Q950'], col="blue", on=1, lty=2, lwd=0.5)
      sim.plot <- lines(df.combined[,'Q975'], col="blue", on=1, lty=1, lwd=1)
      sim.plot <- lines(df.combined[,'Q990'], col="blue", on=1, lty=1, lwd=2)
      sim.plot <- lines(oos.series[,1], col="purple", on=1, lty=1, lwd=2)
      sim.plot <- points(tail(df.combined[,'mean'],1), col="red", pch=16)
      
      output <- list(results,
                     sim.plot)
      
      return(output)
    }
  }
}

empirical.simulation.simple <- function(data, 
                                        n.sim, 
                                        t, 
                                        p,
                                        mean, 
                                        variance, 
                                        skew, 
                                        kurt){
  # data has to be XTS [OHLC] object
  if(is.xts(data) == FALSE){
    print("Error: data field is required in XTS format [OHLC].")
  } else if(is.xts(data) == TRUE){
    
    require(lubridate)
    require(bizdays)
    require(xts)
    require(quantmod)
    require(PearsonDS)
    
    set.seed(1)
    days.ahead = t
    n.sim <- n.sim
    emprical.sim <- matrix(nrow = days.ahead, ncol=n.sim)
    
    for(i in 1:n.sim){
      df <- rpearson(n=days.ahead, moments=c(mean=mean, variance=variance, skewness=skew, kurtosis=kurt))
      emprical.sim[,i] <- df
    }
    
    sim <- as.data.frame(emprical.sim)
    
    sim$quantile <- NA
   
    sim.cumulative <- apply(sim[,1:(ncol(sim)-1)], FUN = function(x){cumprod(1+x)}, MARGIN = 2)
    sim$quantile <- apply(sim.cumulative, FUN = function(x){quantile(na.omit(x),p)}, MARGIN = 1)

    # adjust for business days
    business.calendar <- create.calendar('my_calendar', weekdays = c('saturday','sunday'))
    sim.dates <- bizdays::offset(ymd(last(index(data))+1), 1:days.ahead, cal = business.calendar)
    
    df <- cbind(as.data.frame(sim.dates),sim[,'quantile'])
    df$empty <- NA
    
    df <- as.xts(df[,-1], order.by=ymd(df$sim.dates))
    df <- as.numeric(tail(Cl(data),1)) * df
    
    df.combined <- rbind(Cl(data),df$empty)
    df.combined <- merge(df.combined,df[,1:(ncol(df)-1)]) # do not merge the empty column
    df.combined[1,2:ncol(df.combined)] <- df.combined[1,1]
    df.combined <- df.combined[,-1]
    colnames(df.combined) <- "quantile"
    #quantile(df.combined[nrow(df.combined),], 0.50)
    
    VaR.forecast <- tail(df.combined[,'quantile'],1)
    
    colnames(VaR.forecast) <- paste0("Q",p)
    
    return(VaR.forecast)
    
  }
}

extract_params <- function(list){
  a <- list[c(4:7)]
  df <- do.call(rbind,a)
  return(df)
}

ExcludeDates <- function(x, exclude) {
  idx <- index(x)
  x[idx[!format(idx, "%Y-%m-%d") %in% paste(exclude)]]
}

compute_stats <- function(object){
  
  min <- min(as.numeric(object))
  max <- max(as.numeric(object))
  p.25 <- as.numeric(quantile(as.numeric(object),0.25))
  mean <- mean(as.numeric(object))
  p.50 <- as.numeric(quantile(as.numeric(object),0.50))
  p.75 <- as.numeric(quantile(as.numeric(object),0.75))
  stdev <- sd(as.numeric(object))
  skew <- moments::skewness(as.numeric(object))
  kurt <- moments::kurtosis(as.numeric(object))
  
  result <- round(data.frame(min, max, p.25, mean, p.50, p.75, stdev, skew, kurt),4)
  colnames(result) <- c("min","max","p.25","mean","p.50","p.75","stdev","skewness","kurtosis")
  return(result)
  
}

ajustar_pearson_iv <- function(x, verbose = TRUE) {
  if (!requireNamespace("moments", quietly = TRUE)) {
    install.packages("moments")
  }
  if (!requireNamespace("cubature", quietly = TRUE)) {
    install.packages("cubature")
  }
  
  library(moments)
  library(cubature)
  
  # Momentos empiricos
  g1 <- skewness(x)
  g2 <- kurtosis(x)
  
  # Verificacion condicion Pearson IV
  if (g2 <= 1.5) {
    stop("g2 <= 1.5: La distribucion no puede ser modelada como Pearson Tipo IV.")
  }
  
  # Parámetros base
  mu <- mean(x)
  beta <- sd(x)
  nu <- 2.5 + 3 / (g2 - 1.5)
  lambda <- g1 / sqrt(g2 - 1.5)
  eta <- lambda * nu
  
  # Densidad sin normalizar
  pearsonIV_density <- function(x_val) {
    (1 + ((x_val - mu)^2 / beta^2))^(-nu) * exp(-eta * atan((x_val - mu) / beta))
  }
  
  # Estimacion numerica de K
  integrand <- function(x_vec) pearsonIV_density(x_vec)
  integral_result <- cubature::adaptIntegrate(integrand, lowerLimit = -Inf, upperLimit = Inf)
  K <- 1 / integral_result$integral
  
  # Output
  resultado <- list(
    mu = mu,
    beta = beta,
    g1 = g1,
    g2 = g2,
    nu = nu,
    eta = eta,
    lambda = lambda,
    K = K
    #densidad = function(x_input) K * pearsonIV_density(x_input)
  )
  
  resultado <- do.call(rbind,resultado)
  
  # if (verbose) {
  #   cat("Resultados de Pearson Tipo IV:\n")
  #   cat("Media (mu):", mu, "\n")
  #   cat("Desviación (beta):", beta, "\n")
  #   cat("Sesgo (g1):", g1, "\n")
  #   cat("Curtosis (g2):", g2, "\n")
  #   cat("nu:", nu, "\n")
  #   cat("lambda:", lambda, "\n")
  #   cat("eta:", eta, "\n")
  #   cat("Constante de normalizacion (K):", K, "\n")
  # }
  
  return(resultado)
}

# Descarga GOBIXDR diario ----
gobix <- read.csv(url("https://www.bvrd.com.do/indice/Data/GobixDataIRP.csv")) # descarga
gobix <- gobix[,1:2] # se queda con dos primeras columnas
colnames(gobix)[1] <- "fecha" # renombra la columna 1
gobix$fecha <- lubridate::mdy(gobix[,1]) # la convierte en formato fechas mes-dia-año
gobix <- xts::as.xts(gobix[,-1], order.by=gobix$fecha) # convierte el formato en XTS
colnames(gobix)[1] <- "Close.GOBIX" # renombra la serie (indice al cierre del dia)

# VaR. Aplica el ejercicio del VaR historico visto en clase.
gobix$var.diaria <- CalculateReturns(Cl(gobix), method = "discrete") 
gobix$daily.vol <- diff(log(Cl(gobix)))
gobix$rolling.sd <- roll_sd(gobix$daily.vol, width = 65)*sqrt(252)
gobix <- na.omit(gobix)
gobix$drawdown <- Drawdowns(gobix[,'daily.vol'])

gobix.2023 <- gobix["2022-12-30/2023-12-15"]
gobix.2024 <- gobix["2023-12-30/2024-12-15"]
gobix.oos <- gobix["2022-11-22/2024-12-30"]

gobix.completo <- gobix

# OOS
oos <- gobix["2019::2022"] 

# Traininig: In the Sample
gobix <- gobix["2014::2022"] 
gobix.train <- gobix["2014::2018"] 


gobix$yesterday.ret <- data.table::shift(as.numeric(gobix[,'daily.vol',drop=FALSE]), n=1, type = "lag")
gobix <- na.omit(gobix)
# Asignacion regimenes ----
Regimenes_PM <- read_excel("Regimenes de TPM.xlsx")
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
#periodo.lower <- paste0(lower$ANO,"-",lower$MES)
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
regime.lower <- regime.lower["/2023-01-01"] # ITS

regime.flat <- list()
for(i in 1:length(periodo.flat)){
  regime.flat[[i]] <- gobix.df %>% dplyr::filter(YM == periodo.flat[i])
}

regime.flat <- do.call(rbind,regime.flat)
regime.flat <- regime.flat[,-ncol(regime.flat)]
regime.flat <- as.xts(regime.flat, order.by = ymd(rownames(regime.flat)))
regime.flat_full <- regime.flat
regime.flat <- regime.flat["/2023-01-01"] # ITS

regime.higher <- list()
for(i in 1:length(periodo.higher)){
  regime.higher[[i]] <- gobix.df %>% dplyr::filter(YM == periodo.higher[i])
}

regime.higher <- do.call(rbind,regime.higher)
regime.higher <- regime.higher[,-ncol(regime.higher)]
regime.higher <- as.xts(regime.higher, order.by = ymd(rownames(regime.higher)))
regime.higher_full <- regime.higher
regime.higher <- regime.higher["/2023-01-01"] # ITS

baselIII.12m.stress <- gobix["2021-12-01/2022"] 

# CALIBRACION ----

# resumen autocorrelacion
autocorrelation.table <- cbind(table.Autocorrelation(gobix[,'daily.vol',drop=FALSE]), 
                               table.Autocorrelation(regime.lower[,'daily.vol',drop=FALSE]),
                               table.Autocorrelation(regime.flat[,'daily.vol',drop=FALSE]),
                               table.Autocorrelation(regime.higher[,'daily.vol',drop=FALSE]))

  colnames(autocorrelation.table) <- c('todos', 'Ri=2', 'Ri=3', 'Ri=1')

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
regime.lower.sim <- descdist(as.numeric(regime.lower[,'daily.vol',drop=FALSE]), discrete=FALSE, boot=10000, graph = FALSE)
regime.flat.sim <- descdist(as.numeric(regime.flat[,'daily.vol',drop=FALSE]), discrete=FALSE, boot=10000, graph = FALSE)
regime.higher.sim <- descdist(as.numeric(regime.higher[,'daily.vol',drop=FALSE]), discrete=FALSE, boot=10000, graph = FALSE)
regime.all.sim <- descdist(as.numeric(gobix[,'daily.vol',drop=FALSE]), discrete=FALSE, boot=10000, graph = FALSE)
regime.stress.sim <- descdist(as.numeric(baselIII.12m.stress[,'daily.vol',drop=FALSE]), discrete=FALSE, boot=10000, graph = FALSE)

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
# summary.stats_sim[1,] <- summary.stats_sim[1,]*256
# summary.stats_sim[2,] <- summary.stats_sim[2,]*sqrt(256)
# summary.stats_sim <- round(summary.stats_sim,4)

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

oos.start <- which(index(gobix.completo) == "2023-03-01")
oos.end <- which(index(gobix.completo) == "2024-12-30")

# manage exclusions. 
censored.period <- paste0("2022-11-15/2023-03-01")
censored.dates <- index(gobix.completo[censored.period,'daily.vol'])

# test
# ExcludeDates(baselIII.12m.stress[date.limit,'daily.vol'], exclude = censored.dates)

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
  # forecast.estres[[i]] <- empirical.simulation.simple(data = gobix.completo[i,],
  #                                                     n.sim = 5000,
  #                                                     t = 30,
  #                                                     p = 0.01,
  #                                                     mean = stress_mean[[i]],
  #                                                     variance = stress_variance[[i]],
  #                                                     skew = stress_skew[[i]],
  #                                                     kurt = stress_kurt[[i]])
  # 
  # forecast.sin_cambios[[i]] <- empirical.simulation.simple(data = gobix.completo[i,], 
  #                                                         n.sim = 5000, 
  #                                                         t = 30, 
  #                                                         p = 0.01,
  #                                                         mean = sin.cambios_mean[[i]], 
  #                                                         variance = sin.cambios_variance[[i]], 
  #                                                         skew = sin.cambios_skew[[i]], 
  #                                                         kurt = sin.cambios_kurt[[i]])
  # 
  # forecast.incremento[[i]] <- empirical.simulation.simple(data = gobix.completo[i,], 
  #                                                         n.sim = 5000, 
  #                                                         t = 30, 
  #                                                         p = 0.01,
  #                                                         mean = higher_mean[[i]], 
  #                                                         variance = higher_variance[[i]], 
  #                                                         skew = higher_skew[[i]], 
  #                                                         kurt = higher_kurt[[i]])
  # 
  # forecast.reduccion[[i]] <- empirical.simulation.simple(data = gobix.completo[i,], 
  #                                                         n.sim = 5000, 
  #                                                         t = 30, 
  #                                                         p = 0.01,
  #                                                         mean = lower_mean[[i]], 
  #                                                         variance = lower_variance[[i]], 
  #                                                         skew = lower_skew[[i]], 
  #                                                         kurt = lower_kurt[[i]])
  
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

# Estimacion Pearson IV
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
fwd_p95.df <- read.csv("empirical_VaR_rollingOOS_p95.csv")
fwd_p975.df <- read.csv("empirical_VaR_rollingOOS_p975.csv")
fwd_p99.df <- read.csv("empirical_VaR_rollingOOS_p99.csv")

fwd_oos.all <- list(fwd_p95.df, fwd_p975.df, fwd_p99.df)

fwd.df <- lapply(fwd_oos.all, FUN = function(x){apply(x[,2:ncol(x)], MARGIN = 2, FUN = function(x){quantile(x,0.5)})})
fwd.df <- do.call(rbind,fwd.df)
colnames(fwd.df) <- c('BIII','R(i=1)','R(i=2)','R(i=3)')
rownames(fwd.df) <- c('VaR.950','VaR.975','VaR.990')

fwd_all.df1 <- as.data.frame(fwd_oos.all[[3]])
fwd_all.df1$dates <- lubridate::ymd(fwd_all.df1$X)
#fwd_all.df1 <- as.xts(as.data.frame(fwd_oos.all[[3]]), order.by = ymd(rownames(as.data.frame(fwd_oos.all[[3]]))))

selected.range <- paste0(index(df),"/",index(df)+360)
oos.data <- gobix[selected.range][,1]

sim.stress <- empirical.simulation(data = df, oos.series = NULL, n.sim = 5000, t = 250, mean=summary.stats_sim['mean','BIII'], variance=summary.stats_sim['var','BIII'], skew=summary.stats_sim['skewness','BIII'], kurt=summary.stats_sim['kurtosis','BIII'], main="FHS Stress Times")
sim.higher <- empirical.simulation(data = df,  oos.series = NULL, n.sim = 5000, t = 250, mean=summary.stats_sim['mean','R(i=1)'], variance=summary.stats_sim['var','R(i=1)'], skew=summary.stats_sim['skewness','R(i=1)'], kurt=summary.stats_sim['kurtosis','R(i=1)'], main="FHS Higher Rates")
sim.lower <- empirical.simulation(data = df,  oos.series = NULL, n.sim = 5000, t = 250, mean=summary.stats_sim['mean','R(i=3)'], variance=summary.stats_sim['var','R(i=3)'], skew=summary.stats_sim['skewness','R(i=3)'], kurt=summary.stats_sim['kurtosis','R(i=3)'], main="FHS Lower Rates")
sim.flat <- empirical.simulation(data = df,  oos.series = NULL, n.sim = 5000, t = 250, mean=summary.stats_sim['mean','R(i=2)'], variance=summary.stats_sim['var','R(i=2)'], skew=summary.stats_sim['skewness','R(i=2)'], kurt=summary.stats_sim['kurtosis','R(i=2)'], main="FHS No Changes in Rates")
# quantile(na.omit(as.numeric(sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),])),0.050)
# quantile(na.omit(as.numeric(sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),])),0.025)
# quantile(na.omit(as.numeric(sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),])),0.001)

x.higher <- quantile(na.omit(as.numeric(sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),])),0.50)
x.lower <- quantile(na.omit(as.numeric(sim.lower[[1]][[1]][nrow(sim.lower[[1]][[1]]),])),0.50)
x.flat <- quantile(na.omit(as.numeric(sim.flat[[1]][[1]][nrow(sim.flat[[1]][[1]]),])),0.50)
x.stress <- quantile(na.omit(as.numeric(sim.stress[[1]][[1]][nrow(sim.stress[[1]][[1]]),])),0.50)

higher.q05 <- quantile(na.omit(as.numeric(sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),])),0.05)
lower.q05 <- quantile(na.omit(as.numeric(sim.lower[[1]][[1]][nrow(sim.lower[[1]][[1]]),])),0.05)
flat.q05 <- quantile(na.omit(as.numeric(sim.flat[[1]][[1]][nrow(sim.flat[[1]][[1]]),])),0.05)
stress.q05 <- quantile(na.omit(as.numeric(sim.stress[[1]][[1]][nrow(sim.stress[[1]][[1]]),])),0.05)

higher.q025 <- quantile(na.omit(as.numeric(sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),])),0.025)
lower.q025 <- quantile(na.omit(as.numeric(sim.lower[[1]][[1]][nrow(sim.lower[[1]][[1]]),])),0.025)
flat.q025 <- quantile(na.omit(as.numeric(sim.flat[[1]][[1]][nrow(sim.flat[[1]][[1]]),])),0.025)
stress.q025 <- quantile(na.omit(as.numeric(sim.stress[[1]][[1]][nrow(sim.stress[[1]][[1]]),])),0.025)

higher.q01 <- quantile(na.omit(as.numeric(sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),])),0.01)
lower.q01 <- quantile(na.omit(as.numeric(sim.lower[[1]][[1]][nrow(sim.lower[[1]][[1]]),])),0.01)
flat.q01 <- quantile(na.omit(as.numeric(sim.flat[[1]][[1]][nrow(sim.flat[[1]][[1]]),])),0.01)
stress.q01 <- quantile(na.omit(as.numeric(sim.stress[[1]][[1]][nrow(sim.stress[[1]][[1]]),])),0.01)

# EXTRACT PCENTILES
# sim.stress[[1]][[1]][,'Q010']
# sim.higher[[1]][[1]][,'Q010']
# sim.flat[[1]][[1]][,'Q010']
# sim.lower[[1]][[1]][,'Q010']

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

# graph-03: CDF simulacion ventana fija OOS
par(mfrow=c(1,1), mar=c(2,3,3,3))
plot(ecdf(na.omit(as.numeric(sim.higher[[1]][[1]][nrow(sim.higher[[1]][[1]]),]))),
     xlab="",
     ylab="Densidad cumulativa",
     ylim=c(0,0.10),
     xlim=c(40,120),
     main="Estimación VaR ventana fija - por estado R(i)",
     col="red", lwd=1.5)
lines(ecdf(na.omit(as.numeric(sim.lower[[1]][[1]][nrow(sim.lower[[1]][[1]]),]))), col="forestgreen", lwd=1.5)
lines(ecdf(na.omit(as.numeric(sim.flat[[1]][[1]][nrow(sim.flat[[1]][[1]]),]))), col="blue", lwd=1.5)
lines(ecdf(na.omit(as.numeric(sim.stress[[1]][[1]][nrow(sim.stress[[1]][[1]]),]))), col="red4", lwd=1.5)
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

text(x=round(as.numeric(stress.q01),0)-11, y=0.01, labels=paste0(round(as.numeric(stress.q01)-100,1),"%"), pos=4, cex=1, col = "red4")
text(x=round(as.numeric(stress.q025),0)-10, y=0.025, labels=paste0(round(as.numeric(stress.q025)-100,1),"%"), pos=4, cex=1, col = "red4")
text(x=round(as.numeric(stress.q05),0)-11, y=0.05, labels=paste0(round(as.numeric(stress.q05)-100,1),"%"), pos=4, cex=1, col = "red4")

text(x=round(as.numeric(higher.q01),0)+1, y=0.01, labels=paste0(round(as.numeric(higher.q01)-100,1),"%"), pos=4, cex=1, col = "red")
text(x=round(as.numeric(higher.q025),0)+1, y=0.025, labels=paste0(round(as.numeric(higher.q025)-100,1),"%"), pos=4, cex=1, col = "red")
text(x=round(as.numeric(higher.q05),0)+1, y=0.05, labels=paste0(round(as.numeric(higher.q05)-100,1),"%"), pos=4, cex=1, col = "red")

text(x=round(as.numeric(lower.q01),0)+1, y=0.01, labels=paste0(round(as.numeric(lower.q01),1)-1*100,"%"), pos=4, cex=1, col = "green4")
text(x=round(as.numeric(lower.q025),0)+1, y=0.025, labels=paste0(round(as.numeric(lower.q025)-100,1),"%"), pos=4, cex=1, col = "green4")
text(x=round(as.numeric(lower.q05),0)+1, y=0.05, labels=paste0(round(as.numeric(lower.q05)-100,1),"%"), pos=4, cex=1, col = "green4")

text(x=round(as.numeric(flat.q01),0)-11, y=0.01, labels=paste0(round(as.numeric(flat.q01)-100,1),"%"), pos=4, cex=1, col = "blue")
text(x=round(as.numeric(flat.q025),0)-10, y=0.025, labels=paste0(round(as.numeric(flat.q025)-100,1),"%"), pos=4, cex=1, col = "blue")
text(x=round(as.numeric(flat.q05),0)-9, y=0.05, labels=paste0(round(as.numeric(flat.q05)-100,1),"%"), pos=4, cex=1, col = "blue")

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
#abline(h=0.50, col='grey', lty=2, lwd=0.75)
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
# regression.result$coefficients[2,1]
# regression.result$coefficients[2,2]

# CONT'par(mar=c(2,4,2,4), oma=c(0,0,0,0))
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
                                    boot=10000, 
                                    title = "TPM: Reducción",
                                    subtitle = "R(i=3)", 
                                    label.align = 0.5)

regime.flat.sim <- descdist_custom(as.numeric(regime.flat[,'daily.vol',drop=FALSE]), 
                                   discrete=FALSE, 
                                   boot=10000, 
                                   title = "TPM: Sin cambios",
                                   subtitle = "R(i=2)",
                                   label.align = 0.5)

regime.higher.sim <- descdist_custom(as.numeric(regime.higher[,'daily.vol',drop=FALSE]), 
                                     discrete=FALSE, 
                                     boot=10000, 
                                     title = "TPM: Incremento",
                                     subtitle = "R(i=1)",
                                     label.align = 0.5)

regime.flat.sim <- descdist_custom(as.numeric(baselIII.12m.stress[,'daily.vol',drop=FALSE]), 
                                   discrete=FALSE, 
                                   boot=10000, 
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
par(mfrow=c(3,1), mar=c(4,4,4,4), oma=c(0,0,0,0))
lower.quantiles <- as.data.frame(sim.lower[[1]][[1]][,c('Q010','Q025','Q050','Q950','Q975','Q990')])
lower.quantiles <- apply(lower.quantiles, MARGIN = 2, FUN = function(x){x-100})
plot(lower.quantiles[,'Q990'], lower.quantiles[,'Q010'], 
     ylim = rev(range(lower.quantiles[,'Q010'])),
     type="l", lwd=2, col="green4", xlab="nivel cola derecha (%)", ylab="nivel cola izquierda (%)", main="R(i=3)",
     cex=1.25, cex.axis=1.25, cex.lab=1, cex.main=1.25)
lines(lower.quantiles[,'Q975'], lower.quantiles[,'Q025'], lwd=2, col="green2")
lines(lower.quantiles[,'Q950'], lower.quantiles[,'Q050'], lwd=2, col="green")
box(col="grey")

flat.quantiles <- as.data.frame(sim.flat[[1]][[1]][,c('Q010','Q025','Q050','Q950','Q975','Q990')])
flat.quantiles <- apply(flat.quantiles, MARGIN = 2, FUN = function(x){x-100})
plot(flat.quantiles[,'Q990'], flat.quantiles[,'Q010'],type="l", 
     ylim = rev(range(flat.quantiles[,'Q010'])),
     col="blue4", lwd=2, xlab="nivel cola derecha (%)", ylab="nivel cola izquierda (%)", main="R(i=2)",
     cex=1.25, cex.axis=1.25, cex.lab=1, cex.main=1.25)
lines(flat.quantiles[,'Q975'], flat.quantiles[,'Q025'], lwd=2, col="blue2")
lines(flat.quantiles[,'Q950'], flat.quantiles[,'Q050'], lwd=2, col="blue")
legend("topleft", legend = c("[P.010,P.990]","[P.025,P.975]","[P.050,P.950]"), lwd=2,
       col = c("blue4","blue2","blue"), text.col = c("blue4","blue2","blue"), bty = "n", cex=1, y.intersp = 1.25)
box(col="grey")

higher.quantiles <- as.data.frame(sim.higher[[1]][[1]][,c('Q010','Q025','Q050','Q950','Q975','Q990')])
higher.quantiles <- apply(higher.quantiles, MARGIN = 2, FUN = function(x){x-100})
plot(higher.quantiles[,'Q990'], higher.quantiles[,'Q010'], type="l", 
     ylim = rev(range(higher.quantiles[,'Q010'])),
     col="red4", lwd=2, xlab="nivel cola derecha (%)", ylab="nivel cola izquierda (%)", main="R(i=1)",
     cex=1.25, cex.axis=1.25, cex.lab=1, cex.main=1.25)
lines(higher.quantiles[,'Q975'], higher.quantiles[,'Q025'], lwd=2, col="red2")
lines(higher.quantiles[,'Q950'],higher.quantiles[,'Q050'], lwd=2, col="red")
box(col="grey")