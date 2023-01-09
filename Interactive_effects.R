#Check the effects of single variables  (TChl and Temp) on mean size and size variance
Chl_b  <- quantile(Data$Ln_TChl, probs = c(0.01, 0.99))
N      <- 100 #Number of discreted points
TChls  <- seq_range(Chl_b, n = N, pretty = T)
Temp_b <- quantile(Data$Temp, probs = c(0.01, 0.99))
Temps  <- seq_range(Temp_b, n = N, pretty = T)

B <- Spectra_keras(Temp  = Temps, 
                 Ln_TChl = TChls)

Mean_VAR_size_2D <- function(dat){
  #Interactive effects of 'Chl' and 'Temp'
  
  #Compute mean and variance of ln ESD
  #Mean log ESD
  Nrep <- ncol(dat$data[[1]])-1
  Mean_LnESD <- matrix(NA, 
                       nr = nrow(dat), 
                       nc = Nrep)
  
  #Variance of log ESD
  VAR_LnESD <- Mean_LnESD
 
  for (i in 1:nrow(dat)){
    #Extract the data
    d <- dat[i,]$data[[1]]
    
    for (k in 1:Nrep){
       Z <- d[,k][[1]]
       Z <- moments(as.double(Z))
      Mean_LnESD[i,k] <- Z$CWM
       VAR_LnESD[i,k] <- Z$VAR     
    }
  }
  
  #Compute average mean size and size variance at each combination of Temp and Chl
  MEAN <- apply(Mean_LnESD, 1, mean)
  MEAN <- matrix(MEAN, nr = length(TChls), nc = length(Temps))
  VAR  <- apply( VAR_LnESD, 1, mean)
  VAR  <- matrix(VAR,  nr = length(TChls), nc = length(Temps))
    
  return(list(CWM = MEAN, VAR = VAR))
}

d <- Mean_VAR_size_2D(B$CDF) 

pdf("Environmental_controls_interaction.pdf", width = 5, height = 9)
par(font.lab  = 1,
    family    = "serif",
    mgp       = c(2.2,1,0),
    oma       = c(3,2,2,2),
    mar       = c(4,4,1,.2),
    mfrow     = c(2,1)
   )
oriChl <- c(0.1, 0.2, 0.5, 1, 2, 5)
image2D(d$CWM, TChls, Temps,
        xaxt = 'n',
        xlab = 'Chl (µg/L)',
        ylab = 'Temperature (ºC)')
axis(1, at = log(oriChl), labels = oriChl)
mtext('A) Mean Ln ESD', adj=0)

image2D(d$VAR, TChls, Temps,
        xaxt = 'n',
        xlab = 'Chl (µg/L)',
        ylab = 'Temperature (ºC)')
axis(1, at = log(oriChl), labels = oriChl)
mtext('B) Variance of Ln ESD', adj=0)
dev.off()
