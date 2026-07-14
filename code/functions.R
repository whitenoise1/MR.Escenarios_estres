# Funciones auxiliares de MR_EscenariosEstres.R (versión modular).
# title: "Escenarios de estrés en portafolios de inversiones: Un enfoque empírico"
# author: Stefan Bolta, FRM.
# Todas las funciones se copian VERBATIM del script original.

# ---- descdist_custom: gráfico Cullen & Frey (personalización de fitdistrplus::descdist) ----
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

# ---- empirical.simulation: simulación Monte Carlo Pearson IV (ventana fija, con gráficos) ----
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

# ---- empirical.simulation.simple: pronóstico VaR al cuantil p ----
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

# ---- extract_params: extrae momentos de un objeto descdist ----
extract_params <- function(list){
  a <- list[c(4:7)]
  df <- do.call(rbind,a)
  return(df)
}

# ---- ExcludeDates: excluye fechas de un objeto xts ----
ExcludeDates <- function(x, exclude) {
  idx <- index(x)
  x[idx[!format(idx, "%Y-%m-%d") %in% paste(exclude)]]
}

# ---- compute_stats: estadísticas descriptivas resumidas ----
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

# ---- ajustar_pearson_iv: solución cerrada Pearson Tipo IV ----
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
