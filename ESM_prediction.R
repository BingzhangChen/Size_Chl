ESM_Chl_file <- 'chl_Omon_CESM2-WACCM_ssp585_r1i1p1f1_gr_201501-210012.nc'

Time <- ncread(ESM_Chl_file, 'time') #Days since 0001-01-01

#Convert Time to year and DOY
newTime <- as.Date(Time, origin = "0001-01-01")
Years <- as.numeric(strftime(newTime, format='%Y'))
Months <- as.numeric(strftime(newTime, format='%m'))

ESM_Lat <- ncread(ESM_Chl_file, 'lat')
ESM_Lon <- ncread(ESM_Chl_file, 'lon')
Z <- ncread(ESM_Chl_file, 'lev_partial')

#Only extract surface Chl
ESM_Chl <- ncread(ESM_Chl_file, 'chl', 
                  start = c(1, 1, 1, 1),
                  count = c(length(ESM_Lon), length(ESM_Lat), 1, length(Time))) 
ESM_Chl <- ESM_Chl*1e6 #Correct the unit of Chl to µg/L
   
ESM_sst_file <- 'tos_Omon_CESM2-WACCM_ssp585_r1i1p1f1_gr_201501-210012.nc'
ESM_SST      <- ncread(ESM_sst_file, 'tos')
#Time_sst     <- ncread(ESM_sst_file, 'time') #The same with Time
#ESM_Lat_sst  <- ncread(ESM_sst_file, 'lat') #The same with ESM_Lat

#Calculate the difference between the average between 2015 and 2020 and that between 2095 and 2100
#First, extract the indices
w1 <- which(Years >= 2015 & Years <= 2019)
w2 <- which(Years >= 2094 & Years <= 2098)
w  <- c(w1, w2)

ESM_grid <- expand.grid(ESM_Lon, ESM_Lat, newTime[w])
colnames(ESM_grid) <- c('Lon', 'Lat', 'Date')
ESM_grid$Temp <- as.vector(ESM_SST[,,w])
ESM_grid$TChl <- as.vector(ESM_Chl[,,w])
ESM_grid$TChl[ESM_grid$TChl < 0 &!is.na( ESM_grid$TChl)] <- NA
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
                       cal.density = F,
                       expand  = F)
  )
  #For each grid, collect the size spectra of each month and compute the mean size and size variance
  system.time(
    d <- Mean_VAR_size(dat = B$CDF, keep.rep = F) 
  )
  ESM_grid[wx,]$CWM <- d$CWM
  ESM_grid[wx,]$VAR <- d$VAR

  save(ESM_grid, file = 'ESM_Mean_VAR.Rdata')
}
