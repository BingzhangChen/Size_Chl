#Check the effects of single variables  (TChl and Temp) on mean size and size variance
N      <- 20 #Number of discreted points
Chl_b  <- quantile(Data$Ln_TChl, probs = c(0.01, 0.99))
TChls  <- seq_range(Chl_b, n = N, pretty = T)
Temp_b <- quantile(Data$Temp, probs = c(0.01, 0.99))
Temps  <- seq_range(Temp_b, n = N, pretty = T)
NO3_b  <- quantile(Data$NO3, probs = c(0.01, 0.99))
NO3s   <- seq_range(NO3_b, n = N, pretty = T)

B <- Spectra_keras(Temp  = Temps, 
                 Ln_TChl = TChls,
                   NO3   = NO3s,
                 cal.density = F)
newD <- B$CDF[,1:(ncol(B$CDF)-1)]
SI <- array(NA, dim = c(Nrep, nrow(B$CDF), 2))

for (i in 1:nrow(B$CDF)){
	for (j in 1:Nrep){
      CDF0 <- B$CDF[i,]$data[[1]]
      dat <- CDF0 %>%
        select(-Size)
    
      S <- Size_spectra_slope_intercept(Cp = dat[,j][[1]], TChl = exp(B$CDF[i,]$Ln_TChl))
      SI[j,i,1] <- S$Intercept
      SI[j,i,2] <- S$Slope
    }
}

#Compute averages
SI_avg <- apply(SI, c(2,3), mean)
newD$Intercept <- SI_avg[,1]
newD$Slope       <- SI_avg[,2]

pdf("Interaction_Slope_Intercept_4panel.pdf", 
    width = 10, height = 9)
par(font.lab  = 1,
    family    = "serif",
    mgp       = c(2.2,1,0),
    oma       = c(3,2,2,2),
    mar       = c(2.5,2.5,1,2),
    mfrow     = c(2,2)
   )
oriChl <- c(0.1, 0.2, 0.5, 1, 2, 5)
NO3_names <- c('Low', 'High')

#plot effects of TChl and Temp on regression slope and intercept at two levels of NO3
for (i in 1:2){
  if (i == 1){
    m <- 1
  }else if (i == 2){
    m <- length(NO3s)
  }else{
    stop('Incorrect i!')
  }
  d1 <- newD %>%
    filter(NO3 == NO3s[m])
  
  V <- matrix(d1$Slope, nr = length(TChls), nc = length(Temps))
  image2D(V, TChls, Temps,
          zlim = range(newD$Slope),
          xaxt = 'n',
          xlab = '',
          ylab = '')
  axis(1, at = log(oriChl), labels = oriChl)
  mtext(paste0(LETTERS[i],') Slope of size spectra, ', NO3_names[i],' NO3'), adj=0)
  
  V <- matrix(d1$Intercept, nr = length(TChls), nc = length(Temps))
  image2D(V, TChls, Temps,
          zlim = range(newD$Intercept),
          xaxt = 'n',
          xlab = '',
          ylab = '')
  axis(1, at = log(oriChl), labels = oriChl)
  mtext(paste0(LETTERS[i+1],') Intercept of size spectra, ',
               NO3_names[i],' NO3'), adj = 0)
}
mtext('Chl (µg/L)', side = 1, outer=T, adj = .5)
mtext('Temperature (ºC)', side = 2, outer=T, adj = .5)
dev.off()
