#Check the effects of single variables  (TChl and Temp) on mean size and size variance
Chl_b <- quantile(Data$Ln_TChl, probs = c(0.01, 0.99))
N     <- 100 #Number of discreted points
TChls <- seq_range(Chl_b, n = N, pretty = T)

B <- Spectra_keras(Temp  = median(Data$Temp), 
                 Ln_TChl = TChls)

d <- Mean_VAR_size(B$CDF) 
d$CWM$TChl <- exp(d$CWM$TChl)
d$VAR$TChl <- exp(d$VAR$TChl)

CWM_chl_plot <- ggplot(data = d$CWM, mapping = aes(x = TChl, y = CWM)) +
    geom_smooth() +
    scale_x_log10() +
    xlab('') +
    ylab('Community mean Ln ESD (µm)') +
    theme_light()

VAR_chl_plot <- ggplot(data = d$VAR, mapping = aes(x = TChl, y = VAR)) +
    geom_smooth() +
    scale_x_log10() +
    xlab('Chl (µg/L)') +
    ylab('Variance of Ln ESD') +
    theme_light()


Temp_b <- quantile(Data$Temp, probs = c(0.01, 0.99))
Temps  <- seq_range(Temp_b, n = N, pretty = T)
G <- Spectra_keras(Temp  = Temps, 
                 Ln_TChl = median(Data$Ln_TChl))
dd <- Mean_VAR_size(G$CDF, varname = 'Temp')

CWM_temp_plot <- ggplot(data = dd$CWM, mapping = aes(x = Temp, y = CWM)) +
    geom_smooth() +
    xlab('') +
    ylab('') +
    theme_light()

VAR_temp_plot <- ggplot(data = dd$VAR, mapping = aes(x = Temp, y = VAR)) +
    geom_smooth() +
    xlab('Temperature (ºC)') +
    ylab('') +
    theme_light()

pdf("Environmental_controls_mean_var.pdf", width = 9, height = 9)
ggarrange(CWM_chl_plot, 
          CWM_temp_plot,
          VAR_chl_plot,
          VAR_temp_plot,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
dev.off()