#Load libraries
library(tidyverse)
library(keras)
library(tensorflow)
library(tidymodels)
library(tfdatasets)
library(modelr)
library(ggpubr)
library(plot3D)
library(ncdf4)

#Function reading data from netcdf files
ncread  <- function(file, VAR, start = NA, count = NA){
  if(file.exists(file)){
    nc    <- nc_open(file)
  }else{
    stop(paste(file,'does not exist!')   )
  }
  data    <- ncvar_get(nc, VAR, start = start, count = count)
  nc_close(nc)
  return(data)
}

jet.colors   <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                   "#7FFF7F", "yellow", "#FF7F00", "red",
                                   "#7F0000")) 

options(tibble.width = Inf)
#Set random seed
set.seed(1)

#Load relevant datasets
load('train_test_full.Rdata')

#Convert Perc back to [-inf, inf] using the sigmoid curve
sigmoid <- function(x)-log(1/x - 1)

#back converting to percentage
sigmoidrev <- function(x)1/(1+exp(-x))

#Transform the predictors
train_df <- train_df %>% 
  mutate(Ln_Zb   = log(abs(Zb))) %>%
  mutate(Ln_MLD  = log(MLD))     %>%
  mutate(NO3     = sqrt(sqrt(NO3))) %>%
  mutate(Ln_TChl = log(TChl))    %>%
  mutate(Ln_PAR  = log(PAR))     %>%
  mutate(Ln_dFe  = log(dFe))     %>%
  select(-Zb, -MLD, -TChl, -PAR, -dFe)
  
test_df <- test_df %>% 
  mutate(Ln_Zb   = log(abs(Zb))) %>%
  mutate(Ln_MLD  = log(MLD))     %>%
  mutate(NO3     = sqrt(sqrt(NO3))) %>%
  mutate(Ln_TChl = log(TChl))    %>%
  mutate(Ln_PAR  = log(PAR))     %>%
  mutate(Ln_dFe  = log(dFe))     %>%
  select(-Zb, -MLD, -TChl, -PAR, -dFe)

#Merge the full dataset
Full_data <- full_join(train_df, test_df)

#Function depending on the dataset and number of units
SFChl_keras_model <- function(train_data = train_df, m = ncol(train_data)-1){
  
  #Normalize features
  spec <- feature_spec(train_data, label ~ . ) %>% 
    tfdatasets::step_numeric_column(all_numeric(), 
                                    normalizer_fn = scaler_standard()) %>% 
    fit()
  
  #wrap the model building code into a function
  build_model <- function() {
    input <- layer_input_from_dataset(train_data %>%
                                        select(-label))
    
    output <- input %>% 
      layer_dense_features(dense_features(spec)) %>% 
      layer_dense(units = m, activation = "relu") %>%
      layer_dense(units = m, activation = "relu") %>%
      layer_dense(units = 1) 
    
    model <- keras_model(input, output)
    
    #Compile the code
    model %>% 
      compile(
        loss = "mse",
        optimizer = "rmsprop",
        metrics = list("mean_absolute_error")
      )
    
    model
  }
  
  return(build_model())
}

# The patience parameter is the amount of epochs to check for improvement.
early_stop <- callback_early_stopping(monitor = "val_loss", patience = 20)
 
Fitted_model <- function(TRAIN_DATA){
  
  m     <- ncol(TRAIN_DATA) - 1
  model <- SFChl_keras_model(train_data = TRAIN_DATA, m = m)
  
  #Train the model
  history <- model %>% fit(
    x = TRAIN_DATA %>% select(-label),
    y = TRAIN_DATA$label,
    epochs = 500,
    validation_split = .2,
    batch_size = 16,
    verbose = 0,
    callbacks = list(early_stop)
  )
  
  #Obtain the best epochs
  best_epoch <- length(history$metrics$loss)
  
  #Train the final model with all data
  model <- SFChl_keras_model(train_data = TRAIN_DATA, m = m)
  model %>% fit(
    x = TRAIN_DATA %>% select(-label),
    y = TRAIN_DATA$label,
    epochs = best_epoch,
    batch_size = 16,
    verbose = 0
  )
  return(model)
}

