#This script is the main script that provides the flow of the size-fractionated Chl analysis
#Load data and packages
source('Model_prep.R')

#Check input fields
#source('Exploration.R')

#Comparing different models using 10-fold cross validation
#source('varImp_keras.R') #Results are that the models of TChl, size, Temp, NO3 is good enough

#Final model for prediction
source('Final_keras_model.R')

#Plot an example phytoplankton size spectrum
#source('Example_spectrum.R')

#Plot the partial effects of each individual predictors on mean size and size variance
#source('Single_variable_effect_spectra.R')

#Plot the partial effects of each individual predictors on mean size and size variance
#source('Interactive_effects.R')

# Plot the global distributions of means size and size variance
#source('GLobal_size_Spectra.R')

# Predict the future changes of phytoplankton size spectra
source('ESM_prediction.R')
