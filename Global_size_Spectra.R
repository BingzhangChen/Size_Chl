load('Global_grid_SFChl.Rdata')

Global_grid <- Global_grid %>%
  select(Longitude, Latitude, DOY, Temp, NO3, TChl) %>%
  mutate(CWM = NA) %>%
  mutate(VAR = NA) %>%
  mutate(Intercept = NA) %>%
  mutate(Slope = NA)

#Run prediction by season
for (k in 2:length(DOYs)){

  #Extract non-NA data in each month
  wx <- which(!is.na(Global_grid$Temp) & Global_grid$DOY == DOYs[k])
  
  #Extrapolate global surface Chl at a 1x1 grid every month
  B <- Spectra_keras(Temp    =     Global_grid[wx,]$Temp, 
                     Ln_TChl = log(Global_grid[wx,]$TChl),
                     NO3     = sqrt(sqrt(Global_grid[wx,]$NO3)),
                     cal.density = F,
                     expand  = F)

  #For each grid, collect the size spectra of each month and compute the mean size and size variance
  d <- Mean_VAR_size(dat = B$CDF, keep.rep = F) 

  #Compute slope and intercept of each size spectra
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
  Global_grid[wx,]$Intercept <- SI_avg[,1]
  Global_grid[wx,]$Slope     <- SI_avg[,2]
  Global_grid[wx,]$CWM       <- d$CWM
  Global_grid[wx,]$VAR       <- d$VAR
  save(Global_grid, file = 'Global_CWM_VAR_slope_wNO3.Rdata')
}


#Calculate the annual mean CWM and VAR
CWM <- array(Global_grid$CWM, dim = c(length(Lons), length(Lats), length(DOYs)))
VAR <- array(Global_grid$VAR, dim = c(length(Lons), length(Lats), length(DOYs))) 

meanCWM     <- apply(CWM, c(1,2), function(x)mean(x, na.rm=T))
CWM_latmean <- apply(meanCWM, 2,  function(x)mean(x, na.rm=T))

meanVAR     <- apply(VAR, c(1,2), function(x)mean(x, na.rm=T))
VAR_latmean <- apply(meanVAR, 2,  function(x)mean(x, na.rm=T))

pdf("Global_CWM_VAR_wNO3.pdf", width = 8, height = 9)
par(font.lab  = 1,
    family    = "serif",
    mgp       = c(2.2,1,0),
    oma       = c(3,2,2,2),
    mar       = c(4,4,1,.2)
   )
nf <- layout(matrix(c(1,1,2, 3,3,4), 2, 3, byrow = TRUE))

#Plot mean ESD
z <- exp(meanCWM)
CWMmax <- quantile(z, probs = 0.975, na.rm=T)
z[z > CWMmax] <- CWMmax
image2D(z, Lons, Lats,
        col = jet.colors(18), 
        frame= F,
        xlab = "", 
        ylab = "Latitude (ºN)",
        cex.lab = 1.2,
        cex.axis= 1.2,
        main = '')
mtext('A) Mean ESD (µm)', adj=0)

plot(exp(CWM_latmean), Lats, 
     xlab = 'Mean ESD (µm)',
     ylab = '',
     cex.lab = 1.2,
     cex.axis= 1.2,
     type = 'l'
     )
mtext('B)', adj=0)

z <- meanVAR
VARmax <- quantile(z, probs = 0.975, na.rm=T)
z[z > VARmax] <- VARmax
image2D(z, Lons, Lats,
        col = jet.colors(18), 
        frame= F,
        xlab = "Longitude (ºE)", 
        ylab = "Latitude (ºN)",
     cex.lab = 1.2,
     cex.axis= 1.2,
        main = '')
mtext(expression(paste('C) Size variance ((ln µm)'^2*')')), adj=0)
plot(VAR_latmean, Lats, 
     xlab = 'Variance of Ln ESD',
     ylab = '',
     cex.lab = 1.2,
     cex.axis= 1.2,
     type = 'l'
     )
mtext('D)', adj=0)
dev.off()


#Obtain size spectra for 11.5 W 56.5 N
LON1 <- -11.5
LAT1 <- 56.5
wx   <- which.min((Lons-LON1)**2)
wy   <- which.min((Lats-LAT1)**2)

cff <- which(Global_grid$Longitude == Lons[wx] & Global_grid$Latitude  == Lats[wy])



#Plot CWM in Jan.
wr <- which(Global_grid$DOY == DOYs[1])



