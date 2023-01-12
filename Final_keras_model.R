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

#Run the model with TChl, Size, and temperature
xname <- c('label', 'Size', 'Ln_TChl', 'Temp')
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
Spectra_keras <- function(Temp, Ln_TChl, cal.density = T, expand = T){
  #This function computes the size spectrum, mean size and size variance based on Temperature, Chl, and NO3
  
  if (expand) { #Need to construct combinations of Temp and Chl and add Size
    
    newX <- expand.grid(
                       Size    = Newsize,
                       Ln_TChl = Ln_TChl, 
                       Temp    = Temp
                     )
  
  }else{ #Take original inputs of Temp and Ln_TChl to construct predictive dataframe
    G  <- data.frame(Temp = Temp, Ln_TChl = Ln_TChl) %>%
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
    mutate(Ln_TChl = as.double(Ln_TChl)) 
 
  #A new dataframe holding the values of cumulative probability
  newY <- matrix(NA, nc = Nrep, nr = nrow(newX)) %>%
    as.data.frame() %>%
    mutate(Size    = newX$Size) %>%
    mutate(Temp    = newX$Temp) %>%
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
  newY <- newY %>%
    group_by(ID, Temp, Ln_TChl) %>%
    nest()
  
  if (cal.density){
    
    #Compute average and sd of density of each unique combination of environmental factors
    if (expand) { #Need to construct combinations of Temp and Chl and add Size
      newX1 <- expand.grid(
                    Size    = Mid_size,
                    Ln_TChl = Ln_TChl, 
                    Temp    = Temp
                  )
    } else{
      G  <- data.frame(Temp = Temp, Ln_TChl = Ln_TChl)
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
      mutate(Ln_TChl = newX1$Ln_TChl) %>%
      group_by(Temp, Ln_TChl) %>%
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
Mean_VAR_size <- function(dat, keep.rep = T, varname = 'Chl'){
  #varname can be either 'Chl' or 'Temp'
  
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
    
    if (varname == 'Chl'){
      MEAN <- MEAN %>%
        mutate(TChl = TChls[id])
        
      VAR <- VAR %>%
        mutate(TChl = TChls[id])
    }else if (varname == 'Temp'){
       MEAN <- MEAN %>%
        mutate(Temp = Temps[id])
        
      VAR <- VAR %>%
        mutate(Temp = Temps[id])   
    }else{
      #do nothing
    }
  }else{
    #Simply calculate the mean CWM and VAR
    MEAN <- apply(Mean_LnESD, 1, function(x)mean(x, na.rm=T))
    VAR  <- apply( VAR_LnESD, 1, function(x)mean(x, na.rm=T))
  }
 
  return(list(CWM = MEAN, VAR = VAR))
}
