#Construct an example size spectrum
Lmin    <- 0.2
Lmax    <- 200
N       <- 100 #Number of size classes
Newsize <- exp(seq(log(Lmin), log(Lmax), length.out=N))

#The interval of the Newsize vector
dsize   <- diff(Newsize)

#Midpoint of Newsize for computing density profile
Mid_size <- numeric(length(dsize))

for (i in 1:length(Mid_size)) Mid_size[i] <- sqrt(Newsize[i]* 
                                                  Newsize[i+1])

#Convert Mid_size  to biovolume
Vol <- pi/6*Mid_size^3

#Convert biovolume to carbon (pgC cell-1)
Size_class_C <- 0.216 *Vol^0.939
 
#Run the model with TChl, Size,  temperature, and NO3
xname <- c('label', 'Size', 'Ln_TChl', 'Temp', 'NO3')
Data  <- Full_data %>% 
  select(all_of(xname))

fitted_models <- function(Nrep = 10){
  models <- vector(mode = "list", length = Nrep)
  for (i in 1:Nrep){
    models[[i]] <- Fitted_model(TRAIN_DATA = Data)
  }
  return(models)
}

Nrep <- 10
models <- fitted_models()

#Repeat the model simulation for 10 times and compute the average
moments <- function(Cp){
  #Cp: cumulative probability (vector), has to monotonically increasing
  #output: mean and variance of this commulative probability
  
  N <- length(Cp) 
  D <- diff(Cp)/dsize
  
  CWM <- sum(D * log(Mid_size) * dsize) / 
          sum(D                * dsize)
  
  #Variance of log ESD
  VAR <- sum(D * (log(Mid_size) - CWM)**2 * dsize)/ 
             sum(D * dsize)
  
  return(list(CWM = CWM, VAR = VAR))
 
}
Spectra_keras <- function(Temp, Ln_TChl, NO3, cal.density = T, expand = T){
  #This function computes the size spectrum, mean size and size variance based on Temperature, Chl, and NO3
  
  if (expand) { #Need to construct combinations of Temp and Chl and add Size
    
    newX <- expand.grid(
                       Size    = Newsize,
                       Ln_TChl = Ln_TChl, 
                       Temp    = Temp,
                       NO3     = NO3
                     )
  
  }else{ #Take original inputs of Temp, NO3, Ln_TChl to construct predictive dataframe
    G  <- data.frame(Temp = Temp, Ln_TChl = Ln_TChl, NO3 = NO3) %>%
      mutate(ID = row.names(.))
    newX <- matrix(rep(t(G), N), ncol=ncol(G), byrow=TRUE) %>%
      as.data.frame() 
    colnames(newX) <- names(G)

    for (i in 1:N){
      newX[((i-1)*nrow(G) + 1):(i*nrow(G)), 'Size'] <- Newsize[i]
    }
  }
  
  newX <- newX %>%
    mutate(Temp = as.double(Temp)) %>%
    mutate(Ln_TChl = as.double(Ln_TChl)) %>%
    mutate(NO3 = as.double(NO3))
 
  #A new dataframe holding the values of cumulative probability
  newY <- matrix(NA, nc = Nrep, nr = nrow(newX)) %>%
    as.data.frame() %>%
    mutate(Size    = newX$Size) %>%
    mutate(Temp    = newX$Temp) %>%
    mutate(NO3     = newX$NO3) %>%
    mutate(Ln_TChl = newX$Ln_TChl)
  
  if (!expand){
    newY <- newY %>%
      mutate(ID = newX$ID)
    newX <- newX %>% 
      select(-ID) 
  }
  
  for (i in 1:Nrep){
    model <- models[[i]]
    
    #predict the cumulative probability based on each unique combination of TChl and Temp
    newY[, i] <- model %>% 
         predict(newX) %>%
         sigmoidrev()
  }
  
  #Convert newY to a column-list
  if (!expand){
  newY <- newY %>%
    group_by(ID, Temp, Ln_TChl, NO3) %>%
    nest()
  }else{
   newY <- newY %>%
    group_by(Temp, Ln_TChl, NO3) %>%
    nest()   
  }
  
  if (cal.density){
    
    #Compute average and sd of density of each unique combination of environmental factors
    if (expand) { #Need to construct combinations of Temp and Chl and add Size
      newX1 <- expand.grid(
                    Size    = Mid_size,
                    Ln_TChl = Ln_TChl, 
                    Temp    = Temp,
                    NO3     = NO3
                  )
    } else{
      G  <- data.frame(Temp = Temp, Ln_TChl = Ln_TChl, NO3 = NO3)
      newX <- matrix(rep(t(G), N), ncol=ncol(G), byrow=TRUE) %>%
              as.data.frame()
      colnames(newX) <- names(G)
  
      for (i in 1:N){
        newX[((i-1)*nrow(G) + 1):(i*nrow(G)), 'Size'] <- Mid_size[i]
      }       
    }
   
    newD <- matrix(NA, nc = Nrep, nr = nrow(newX1)) %>%
      as.data.frame() %>%
      mutate(Size    = newX1$Size) %>%
      mutate(Temp    = newX1$Temp) %>%
      mutate(NO3     = newX1$NO3)  %>%
      mutate(Ln_TChl = newX1$Ln_TChl) %>%
      group_by(Temp, Ln_TChl, NO3) %>%
      nest()
  
    for (i in 1:nrow(newD)){
      y <- newY[i,]$data[[1]]
      y <- y[, 1:Nrep]
      z <- apply(y, 2, function(x)diff(x)/dsize)   
      newD[i,]$data[[1]][,1:Nrep] <- z
    }
    return(list(CDF = newY,
                PDF = newD ) )
  }else{
    return(list(CDF = newY) )
  }

}

#Compute mean and variance of ln ESD
Mean_VAR_size <- function(dat, keep.rep = T){
  #varname can be either 'Chl' or 'Temp' or 'NO3'
  
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
  
  if (keep.rep){
    #Keep the replications of calculated Mean and Var of size
    MEAN <- as.data.frame(Mean_LnESD) %>%
      reshape(direction = 'long', 
              varying = list(1:ncol(.)),
              v.names = 'CWM') 
      
    VAR <- as.data.frame(VAR_LnESD) %>%
      reshape(direction = 'long', 
              varying = list(1:ncol(.)),
              v.names = 'VAR')

  }else{
    #Simply calculate the mean CWM and VAR
    MEAN <- apply(Mean_LnESD, 1, function(x)mean(x, na.rm=T))
    VAR  <- apply( VAR_LnESD, 1, function(x)mean(x, na.rm=T))
  }
  return(list(CWM = MEAN, VAR = VAR))
}

#Compute mean and variance of ln ESD
Size_spectra_slope_intercept <- function(Cp, TChl, DSize = dsize){
  #Cp: cumulative probability  
  #Mean log ESD
  N    <- length(Cp) 
  D    <- diff(Cp)/DSize #Density probability
  
  #Compute total carbon biomass (mgC/m3)
  theta  <- 1/50 #gChl:gC
  TCarbon <- TChl/theta #total carbon: mgC/m3
  
  #Carbon in each size class
  Carbon <- TCarbon /sum(D*DSize) * (D*DSize)

  #Compute abundance (cells/mL)
  Abun <- Carbon/Size_class_C/1e3 
  
  #Run linear regression of abundance against size
  Dat <- data.frame(Size = log(Size_class_C), Abun = Abun ) %>%
    filter(Abun > 0)
  A_Lm <- lm(log(Abun) ~ Size, data = Dat)
  
# pdf('test.pdf', width = 4, height = 4, page = 'a4')
# plot(log(Size_class_C), log(Abun))
# abline(A_Lm)
# dev.off()
  return(list(Slope = coef(summary(A_Lm))[2,1], 
                    Intercept = coef(summary(A_Lm))[1,1]))
}

