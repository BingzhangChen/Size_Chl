ESM_Chl_file <- 'chl_Omon_CESM2-WACCM_ssp585_r1i1p1f1_gr_201501-210012.nc'

Time <- ncread(ESM_Chl_file, 'time') #Days since 0001-01-01

#Convert Time to year and DOY
newTime <- as.Date(Time, origin = "0001-01-01")
Years   <- as.numeric(strftime(newTime, format='%Y'))
Months  <- as.numeric(strftime(newTime, format='%m'))

ESM_Lat <- ncread(ESM_Chl_file, 'lat')
ESM_Lon <- ncread(ESM_Chl_file, 'lon')
Z       <- ncread(ESM_Chl_file, 'lev_partial')

#Only extract surface Chl
ESM_Chl <- ncread(ESM_Chl_file, 'chl', 
                  start = c(1, 1, 1, 1),
                  count = c(length(ESM_Lon), length(ESM_Lat), 1, length(Time))) 
ESM_Chl <- ESM_Chl*1e6 #Correct the unit of Chl to µg/L
   
ESM_sst_file <- 'tos_Omon_CESM2-WACCM_ssp585_r1i1p1f1_gr_201501-210012.nc'
ESM_SST      <- ncread(ESM_sst_file, 'tos')
Time_sst     <- ncread(ESM_sst_file, 'time') #The same with Time
ESM_Lat_sst  <- ncread(ESM_sst_file, 'lat') #The same with ESM_Lat

ESM_no3_file <- 'no3os_Omon_CESM2-WACCM_ssp585_r1i1p1f1_gr_201501-210012.nc'
ESM_NO3      <- ncread(ESM_no3_file, 'no3os')
ESM_NO3      <- ESM_NO3*1e3 #Convert the unit of NO3 to mmol m^-3
Time_NO3     <- ncread(ESM_no3_file, 'time') #The same with Time
ESM_Lat_no3  <- ncread(ESM_no3_file, 'lat') #The same with ESM_Lat

#Calculate the difference between the average between 2015 and 2020 and that between 2095 and 2100
#First, extract the indices
w1 <- which(Years >= 2015 & Years <= 2019)
w2 <- which(Years >= 2094 & Years <= 2098)
w  <- c(w1, w2)

ESM_grid <- expand.grid(ESM_Lon, ESM_Lat, newTime[w])
colnames(ESM_grid) <- c('Lon', 'Lat', 'Date')
ESM_grid$Temp <- as.vector(ESM_SST[,,w])
ESM_grid$TChl <- as.vector(ESM_Chl[,,w])
ESM_grid$NO3  <- as.vector(ESM_NO3[,,w])
ESM_grid$TChl[ESM_grid$TChl < 0 &!is.na( ESM_grid$TChl)] <- NA
ESM_grid$NO3[ ESM_grid$NO3  < 0 &!is.na( ESM_grid$NO3)]  <- NA
ESM_grid$CWM  <- NA
ESM_grid$VAR  <- NA

#Run prediction by season
DATES <- unique(ESM_grid$Date)
for (k in 1:length(w)){

  #Extract non-NA data in each month
  wx <- which(!is.na(ESM_grid$TChl) & ESM_grid$Date == DATES[k])
  
  system.time(
    B <- Spectra_keras(Temp    =     ESM_grid[wx,]$Temp, 
                       Ln_TChl = log(ESM_grid[wx,]$TChl),
                       NO3     = sqrt(sqrt(ESM_grid[wx,]$NO3)),
                       cal.density = F,
                       expand  = F)
  )
  #For each grid, collect the size spectra of each month and compute the mean size and size variance
  system.time(
    d <- Mean_VAR_size(dat = B$CDF, keep.rep = F) 
  )
  ESM_grid[wx,]$CWM <- d$CWM
  ESM_grid[wx,]$VAR <- d$VAR

  if (k %% 10 == 0) save(ESM_grid, file = 'ESM_Mean_VAR_NO3.Rdata')
}
save(ESM_grid, file = 'ESM_Mean_VAR_NO3.Rdata')

load('ESM_Mean_VAR_NO3.Rdata')

#Calculate the climatological annual mean between 2015 and 2019
ESM_grid <- ESM_grid %>%
  mutate(Year = as.numeric(strftime(Date, format='%Y')))

ESM1 <- ESM_grid %>%
  filter(Year <= 2019) %>%
  group_by(Lat, Lon) %>%
  summarize(CWM = mean(CWM, na.rm = T),
            VAR = mean(VAR, na.rm = T),
            SST = mean(Temp,na.rm = T),
            CHL = mean(TChl,na.rm = T),
            NO3 = mean(NO3, na.rm = T)
            )

ESM2 <- ESM_grid %>%
  filter(Year >= 2094) %>%
  group_by(Lat, Lon) %>%
  summarize(CWM = mean(CWM, na.rm = T),
            VAR = mean(VAR, na.rm = T),
            SST = mean(Temp,na.rm = T),
            CHL = mean(TChl,na.rm = T),
            NO3 = mean(NO3, na.rm = T)
            )

#Transform the dataset
CWM1 <- matrix(ESM1$CWM, 
               nr = length(ESM_Lon),
               nc = length(ESM_Lat))  
CWM2 <- matrix(ESM2$CWM, 
               nr = length(ESM_Lon),
               nc = length(ESM_Lat))  
VAR1 <- matrix(ESM1$VAR, 
               nr = length(ESM_Lon),
               nc = length(ESM_Lat))  
VAR2 <- matrix(ESM2$VAR, 
               nr = length(ESM_Lon),
               nc = length(ESM_Lat))  
SST1 <- matrix(ESM1$SST, 
               nr = length(ESM_Lon),
               nc = length(ESM_Lat))  
SST2 <- matrix(ESM2$SST, 
               nr = length(ESM_Lon),
               nc = length(ESM_Lat))  
CHL1 <- matrix(ESM1$CHL, 
               nr = length(ESM_Lon),
               nc = length(ESM_Lat))  
CHL2 <- matrix(ESM2$CHL, 
               nr = length(ESM_Lon),
               nc = length(ESM_Lat))  
NO31 <- matrix(ESM1$NO3, 
               nr = length(ESM_Lon),
               nc = length(ESM_Lat))  
NO32 <- matrix(ESM2$NO3, 
               nr = length(ESM_Lon),
               nc = length(ESM_Lat))  

#Calculate changes in % compared to baseline
CWM_change_perc <- (exp(CWM2) - exp(CWM1))/exp(CWM1)
VAR_change_perc <- (VAR2 - VAR1)/VAR1
SST_change      <- SST2 - SST1
CHL_change_perc <- (CHL2 - CHL1)/CHL1
NO3_change_perc <- (NO32 - NO31)/NO31

#Keep consistent with previous plotting (by putting Atlantic in the central position)
w1 <- which(ESM_Lon >  180)
w2 <- which(ESM_Lon <= 180)

#Change ESM_lon
ESM_Lon1 <- ESM_Lon[c(w1, w2)]
ESM_Lon1[ESM_Lon1 > 180] <- ESM_Lon1[ESM_Lon1 > 180] - 360

CWM_change_perc1 <- CWM_change_perc[ c(w1, w2),]
VAR_change_perc1 <- VAR_change_perc[ c(w1, w2),]
CHL_change_perc1 <- CHL_change_perc[ c(w1, w2),]
NO3_change_perc1 <- NO3_change_perc[ c(w1, w2),]
SST_change1      <- SST_change[ c(w1, w2),]
CWM1             <- CWM1[ c(w1, w2),]

#Check the transformation is correct
image2D(exp(CWM1), ESM_Lon1, ESM_Lat,
        col = jet.colors(18), 
        frame= F,
        xlab = "", 
        ylab = "Latitude (ºN)",
        cex.lab = 1.2,
        cex.axis= 1.2,
        main = '')
mtext('Mean ESD (µm) 2015-2019', adj=0)

#Plot the differences of mean size and size diversity
pdf("Diff_ESM_CWM_VAR_NO3.pdf", width = 12, height = 9)
par(font.lab  = 1,
    family    = "serif",
    mfrow     = c(3,2),
    mgp       = c(2.2,1,0),
    oma       = c(3,2,2,2),
    mar       = c(2.5,3.1,1,.5)
   )

#Plot changes in mean ESD
CWM_range <- quantile(CWM_change_perc1, probs = c(0.025,0.975), na.rm=T)
C1 <- CWM_range[1]
z <- CWM_change_perc1
z[z <  C1] <- C1
z[z > -C1] <- -C1
image2D(z, ESM_Lon1, ESM_Lat,
        col = jet.colors(18), 
        frame= F,
        xlab = "", 
        ylab = "",
        cex.lab = 1.2,
        cex.axis= 1.2,
        main = '')
mtext('A) Percentage in changes in mean size', adj=0)

#Plot changes in mean ESD
VAR_range <- quantile(VAR_change_perc1, probs = c(0.025,0.975), na.rm=T)
C1 <- VAR_range[1]
z  <- VAR_change_perc1
z[z <  C1] <- C1
z[z > -C1] <- -C1
image2D(z, ESM_Lon1, ESM_Lat,
        col = jet.colors(18), 
        frame= F,
        xlab = "", 
        ylab = "",
        cex.lab = 1.2,
        cex.axis= 1.2,
        main = '')
mtext('B) Percentage in changes in phytoplankton size diversity', adj=0)

#Plot changes in mean SST
SST_range <- quantile(SST_change, probs = c(0.025,0.975), na.rm=T)
C1 <- max(abs(SST_range))
z  <- SST_change1
z[z < -C1] <- -C1
z[z > C1]  <-  C1
image2D(z, ESM_Lon1, ESM_Lat,
        col = jet.colors(18), 
        frame= F,
        xlab = "", 
        ylab = "",
        cex.lab = 1.2,
        cex.axis= 1.2,
        main = '')
mtext('C) Changes in sea surface temperature (ºC)', adj=0)

#Plot changes in mean Chl
CHL_range  <- quantile(CHL_change_perc1, probs = c(0.025,0.975), na.rm=T)
C1         <- max(abs(CHL_range))
z          <- CHL_change_perc1
z[z < -C1] <- -C1
z[z > C1]  <-  C1

image2D(z, ESM_Lon1, ESM_Lat,
        col = jet.colors(18), 
        frame= F,
        xlab = "", 
        ylab = "",
        cex.lab = 1.2,
        cex.axis= 1.2,
        main = '')
mtext('D) Percentage in changes in surface Chlorophyll', adj=0)

#Plot changes in mean NO3
NO3_range  <- quantile(NO3_change_perc1, probs = c(0.025,0.975), na.rm=T)
C1         <- max(abs(NO3_range))
z          <- NO3_change_perc1
z[z < -C1] <- -C1
z[z > C1]  <-  C1

image2D(z, ESM_Lon1, ESM_Lat,
        col = jet.colors(18), 
        frame= F,
        xlab = "", 
        ylab = "",
        cex.lab = 1.2,
        cex.axis= 1.2,
        main = '')
mtext('E) Percentage in changes in surface nitrate', adj=0)

mtext('Longitude (ºE)', adj=0.5, side=1, outer = T)
mtext('Latitude (ºN)',  adj=0.5, side=2, outer = T)

dev.off()