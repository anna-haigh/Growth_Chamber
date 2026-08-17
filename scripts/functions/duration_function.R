duration_tol <- function(fv_fm, time, individual, data = NULL, print_graph = TRUE, graph_dir) {
  
  require(car)
  
  #Assign columns to objects within the function
  
  fv_fm <- eval(substitute(fv_fm), data, parent.frame())
  time <- eval(substitute(time), data, parent.frame())
  individual <- eval(substitute(individual), data, parent.frame())
  
  data <- data.frame(individual, fv_fm, time)  # new data frame to perform the analysis
  
  # Objects to store results across all individuals
  
  ID = Fv_Fm_Max = D50 = D50.lci = D50.uci = Dcrit = Dcrit.lci = Dcrit.uci = D15 = D15.lci = D15.uci = R2 = c()
  
  id_levels <- unique(data$individual)
  
  # Loop over individuals to run the model
  
  for (t in 1:length(id_levels)) {
    
    print(paste("Starting fitting for: ", unique(data$individual)[t], sep = ""))
    
    ptm <- proc.time()
    
    datum2 <- data[data$individual == id_levels[t], ]
    
    datum2 <- droplevels(datum2)
    
    IT <- min(datum2$time)
    
    # Perform linear logit regression to extract initial parameters for theta2 and theta3
    
    cof <- numeric()
    
    cof[1] <- lm(logit(fv_fm) ~ time, data = datum2)$coefficients[1]  # theta 2
    cof[2] <- lm(logit(fv_fm) ~ time, data = datum2)$coefficients[2]  # theta 3
    
    # Fit the logistic decay function of FV/Fm as a function of temperature
    
    hotmod <- tryCatch({
      nls(fv_fm ~ theta1 / (1 + exp(-(theta2 + theta3 * time))),
          start = list(theta1 = 0.8, theta2 = cof[1], theta3 = cof[2]),
          data = datum2,
          trace = F,
          control = list(maxiter = 1000, tol = 1e-3))
    }, error = function(e) {
      message(paste0("  Fit FAILED for ", unique(datum2$individual), ": ", conditionMessage(e)))
      return(NULL)
    })
    
    if (is.null(hotmod)) {
      # Record the ID so it's still in the output, mark everything else NA, and skip
      # bootstrapping / plotting entirely for this individual.
      ID[t] <- as.character(unique(datum2$individual))
      Fv_Fm_Max[t] <- NA
      D50[t] <- D50.lci[t] <- D50.uci[t] <- NA
      Dcrit[t] <- Dcrit.lci[t] <- Dcrit.uci[t] <- NA
      D15[t] <- D15.lci[t] <- D15.uci[t] <- NA
      R2[t] <- NA
      next
    }
    
    # Extract parameters from the non linear model
    
    phia <- coef(hotmod)[1]  # fmax, maximum fluorescense
    Fv_Fm_Max[t] <- phia
    phib <- coef(hotmod)[2]
    phic <- coef(hotmod)[3]
    
    x <- seq(from = min(datum2$time), to = max(datum2$time))
    lx <- length(x)
    y <- phia / (1 + exp(-(phib + phic * x)))
    predict <- data.frame(x, y)
    
    # calculates half of the initial Fv/Fm (D50 target)
    initial <- subset(datum2, datum2$time == IT, select = fv_fm)
    half <- mean(initial[, 1]) / 2
    
    # calculates 85% of the initial Fv/Fm (15% reduction)
    fifteen <- mean(initial[, 1]) * 0.85
    
    # Estimate parameters with 95 Confidence Interval
    predict.boot <- matrix(NA, nrow = 41, ncol = 500)
    D50.half.bs = Dcrit.bs = D15.bs = c()
    
    for (k in 1:100) {  # start loop for bootstrap estimation
      srows <- sample(x = 1:nrow(datum2), size = nrow(datum2), replace = TRUE)
      
      cof2 <- coef(lm(logit(fv_fm) ~ time, data = datum2[srows, ]))
      
      hotmod2 <- tryCatch({
        nls(fv_fm ~ theta1 / (1 + exp(-(theta2 + theta3 * time))),
            start = list(theta1 = .8, theta2 = cof2[1], theta3 = cof2[2]),
            data = datum2[srows, ],
            trace = F,
            control = list(maxiter = 1000, tol = 1e-3))
      }, error = function(e) NULL)
      
      if (!is.null(hotmod2)) {
        phia2 <- coef(hotmod2)[1]
        phib2 <- coef(hotmod2)[2]
        phic2 <- coef(hotmod2)[3]
        
        x2 <- seq(from = min(datum2$time), to = max(datum2$time), length.out = 41)
        y2 <- phia2 / (1 + exp(-(phib2 + phic2 * x2)))
        
        predict.boot[, k] <- y2
        D50.fl1 <- (-log((phia2 / half) - 1) - phib2) / phic2  # estimate D50
        D50.half.bs[k] <- D50.fl1  # estimate D50
        
        D15.fl1 <- (-log((phia2 / fifteen) - 1) - phib2) / phic2  # estimate D15
        D15.bs[k] <- D15.fl1  # estimate T15
        
        predict <- data.frame(x2, y2)  # create the prediction data frame
        df1 <- cbind(predict[-1, ], predict[-nrow(predict), ])
        df1 <- df1[, c(3, 1, 4, 2)]
        df1$slp <- as.vector(apply(df1, 1, function(x) summary(lm((x[3:4]) ~ x[1:2]))[[4]][[2]]))
        max.slp <- round(min(df1$slp), 3)
        slp.at.Dcrit <- max.slp * .15
        b4ct <- which(df1[, 1] < D50.fl1)
        
        if (length(b4ct) > 0) {
          fvfv.at.Dcrit <- df1[which(abs(df1[b4ct, ]$slp - slp.at.Dcrit) == min(abs(df1[b4ct, ]$slp - slp.at.Dcrit))), ][1, 3]
          Dcrit.bs[k] <- round((-log((phia / fvfv.at.Dcrit) - 1) - phib) / phic, 2)
        } else {
          Dcrit.bs[k] <- NA
        }
      } else {
        predict.boot[, k] <- NA
        D50.half.bs[k] <- NA
        Dcrit.bs[k] <- NA
        D15.bs[k] <- NA
      }
    }  # finish loop for bootstrap estimation
    
    # print time ellapsed in bootstrapping
    print(proc.time() - ptm); print(paste(unique(datum2$individual), ": initial fv/fm > 0.7", sep = ""))
    
    predict.boot <- apply(predict.boot, 2, as.numeric)
    D50.half.bs <- as.numeric(D50.half.bs)
    Dcrit.bs <- as.numeric(Dcrit.bs)
    D15.bs <- as.numeric(D15.bs)
    
    safe_quantile <- function(v) {
      v <- v[!is.na(v)]
      if (length(v) == 0) {
        return(c(`2.5%` = NA_real_, `97.5%` = NA_real_))
      }
      quantile(v, c(0.025, 0.975))
    }
    
    D50.half.CI <- safe_quantile(D50.half.bs)
    Dcrit.CI <- safe_quantile(Dcrit.bs)
    D15.CI <- safe_quantile(D15.bs)
    
    ID[t] <- as.character(unique(datum2$individual))  # get individual name
    D50[t] <- if (all(is.na(D50.half.bs))) NA else round(mean(D50.half.bs, na.rm = TRUE), 2)
    D50.lci[t] <- round(D50.half.CI[1], 2)  # D50 lower confidence interval
    D50.uci[t] <- round(D50.half.CI[2], 2)  # D50 upper confidence interval
    
    Dcrit[t] <- if (all(is.na(Dcrit.bs))) NA else round(mean(Dcrit.bs, na.rm = TRUE), 2)
    Dcrit.lci[t] <- round(Dcrit.CI[1], 2)  # Dcrit lower confidence interval
    Dcrit.uci[t] <- round(Dcrit.CI[2], 2)  # Dcrit upper confidence interval
    
    D15[t] <- if (all(is.na(D15.bs))) NA else round(mean(D15.bs, na.rm = TRUE), 2)
    D15.lci[t] <- round(D15.CI[1], 2)  # D15 lower confidence interval
    D15.uci[t] <- round(D15.CI[2], 2)  # D15 upper confidence interval
    
    # Estimate non linear model fit by plotting real Fv/Fm to predicted Fv/Fm
    x3 <- datum2$time
    y3 <- phia / (1 + exp(-(phib + phic * x3)))
    model.fit <- lm(datum2$fv_fm ~ y3)
    R2[t] <- summary(model.fit)$r.squared
    
    # Plot the results for a given individual
    if (!print_graph) {  # start plotting action 1
      par(mfrow = c(1, 2), oma = c(2, 2, 2, 2), mar = c(2, 2, 2, 2))
      
      plot(datum2$fv_fm ~ datum2$time, ylim = c(0, 0.88), xlim = c(0, 7), type = "n")
      text(x = 55, y = 0.88, labels = bquote("F"[v] * "/" * "F"[m] == ""))
      text(x = 61, y = 0.88, labels = paste(round(phia, digits = 3)))
      text(x = 55.8, y = 0.82, labels = bquote("D"[50] == ""))
      text(x = 61.1, y = 0.82, labels = paste(round(D50[t], digits = 2)))
      text(x = 55.8, y = 0.76, labels = bquote("D"[crit] == ""))
      text(x = 61.1, y = 0.76, labels = paste(round(Dcrit[t], digits = 2)))
      text(x = 55.8, y = 0.70, labels = bquote("D"[15] == ""))
      text(x = 61.1, y = 0.70, labels = paste(round(D15[t], digits = 2)))
      mtext(text = bquote(Time ~ "(hours)"), side = 1, line = 2, cex = 1)
      mtext(text = bquote("F"[v] * "/" * "F"[m]), side = 2, line = 2, cex = 1)
      mtext(text = unique(datum2$individual), side = 3, outer = T)
      lines(x, y, col = "#B22222", lwd = 2.5)
      abline(v = round(mean(D50.half.bs, na.rm = TRUE), 2), lty = 2, lwd = 1.5, col = "purple")
      abline(v = round(mean(Dcrit.bs, na.rm = TRUE), 2), lty = 2, lwd = 1.5, col = "orange")
      abline(v = round(mean(D15.bs, na.rm = TRUE), 2), lty = 2, lwd = 1.5, col = "steelblue")
      points(datum2$fv_fm ~ datum2$time, pch = 16, col = "#00000095", cex = 1.5)
      
      plot(datum2$fv_fm ~ y3, ylim = c(0, 0.88), xlim = c(0, 0.88), type = "n")
      text(x = 0.03, y = 0.88, labels = bquote(R^2 == ""))
      text(x = 0.13, y = 0.872, labels = paste(round(summary(model.fit)$r.squared, digits = 3)))
      abline(a = 0, b = 1, col = "#808080", lwd = 2.5)
      abline(model.fit, col = "#B22222", lwd = 2.5)
      mtext(text = bquote(Model ~ "F"[v] * "/" * "F"[m]), side = 1, line = 2, cex = 1)
      mtext(text = bquote(Data ~ "F"[v] * "/" * "F"[m]), side = 2, line = 2, cex = 1)
      points(datum2$fv_fm ~ y3, pch = 16, col = "#00000095", cex = 1.5)
    } else {  # start plotting action 2
      jpeg(filename = paste(graph_dir,  # directory
                            as.character(unique(datum2$individual)),  # individual_id
                            ".jpg", sep = ""), width = 5500, height = 2500, units = "px", res = 550)
      par(mfrow = c(1, 2), oma = c(2, 2, 2, 2), mar = c(2, 2, 2, 2))
      
      plot(datum2$fv_fm ~ datum2$time, ylim = c(0, 0.88), xlim = c(0, 7), type = "n")
      text(x = 55, y = 0.88, labels = bquote("F"[v] * "/" * "F"[m] == ""))
      text(x = 61, y = 0.88, labels = paste(round(phia, digits = 3)))
      text(x = 55.8, y = 0.82, labels = bquote("D"[50] == ""))
      text(x = 61.1, y = 0.82, labels = paste(round(D50[t], digits = 2)))
      text(x = 55.8, y = 0.76, labels = bquote("D"[crit] == ""))
      text(x = 61.1, y = 0.76, labels = paste(round(Dcrit[t], digits = 2)))
      text(x = 55.8, y = 0.70, labels = bquote("D"[15] == ""))
      text(x = 61.1, y = 0.70, labels = paste(round(D15[t], digits = 2)))
      mtext(text = bquote(Time ~ "(hours)"), side = 1, line = 2, cex = 1)
      mtext(text = bquote("F"[v] * "/" * "F"[m]), side = 2, line = 2, cex = 1)
      mtext(text = unique(datum2$individual), side = 3, outer = T)
      lines(x, y, col = "#B22222", lwd = 2.5)
      abline(v = round(mean(D50.half.bs, na.rm = TRUE), 2), lty = 2, lwd = 1.5, col = "purple")
      abline(v = round(mean(Dcrit.bs, na.rm = TRUE), 2), lty = 2, lwd = 1.5, col = "orange")
      abline(v = round(mean(D15.bs, na.rm = TRUE), 2), lty = 2, lwd = 1.5, col = "steelblue")
      points(datum2$fv_fm ~ datum2$time, pch = 16, col = "#00000095", cex = 1.5)
      
      plot(datum2$fv_fm ~ y3, ylim = c(0, 0.88), xlim = c(0, 0.88), type = "n")
      text(x = 0.03, y = 0.88, labels = bquote(R^2 == ""))
      text(x = 0.13, y = 0.872, labels = paste(round(summary(model.fit)$r.squared, digits = 3)))
      abline(a = 0, b = 1, col = "#808080", lwd = 2.5)
      abline(model.fit, col = "#B22222", lwd = 2.5)
      mtext(text = bquote(Model ~ "F"[v] * "/" * "F"[m]), side = 1, line = 2, cex = 1)
      mtext(text = bquote(Data ~ "F"[v] * "/" * "F"[m]), side = 2, line = 2, cex = 1)
      points(datum2$fv_fm ~ y3, pch = 16, col = "#00000095", cex = 1.5)
      dev.off()
    }  # close plotting action 1 & 2
  }  # close individual loop
  
  # Get the output into a single dataframe 
  dur_output <- data.frame(
    ID = ID,
    Fv_Fm_Max = Fv_Fm_Max,
    D50 = D50, D50.uci = D50.uci, D50.lci = D50.lci,  # D50
    Dcrit = Dcrit, Dcrit.uci = Dcrit.uci, Dcrit.lci = Dcrit.lci,  # Dcrit
    D15 = D15, D15.uci = D15.uci, D15.lci = D15.lci,  # D15
    R2 = R2
  )
  
  return(dur_output)
}