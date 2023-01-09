ESM_Chl_file <- 'chl_Omon_CESM2-WACCM_ssp585_r1i1p1f1_gr_201501-210012.nc'

Time <- ncread(ESM_Chl_file, 'time') #Days since 0001-01-01
ESM_Lat <- ncread(ESM_Chl_file, 'lat')
ESM_Lon <- ncread(ESM_Chl_file, 'lon')
ESM_Chl <- ncread(ESM_Chl_file, 'chl') 

ESM_sst_file <- 'tos_Omon_CESM2-WACCM_ssp585_r1i1p1f1_gr_201501-210012.nc'
ESM_SST <- ncread(ESM_sst_file, 'tos')
Time_sst <- ncread(ESM_sst_file, 'time')
ESM_Lat_sst <- ncread(ESM_sst_file, 'lat')
#Calculate the difference between the average between 2015 and 2020 and that between 2095 and 2100
