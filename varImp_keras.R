Xnames <- list(c('label', 'Size', 'Ln_TChl'),
               c('label', 'Size', 'Ln_TChl', 'Temp'),
               c('label', 'Size', 'Ln_TChl', 'NO3'),
               c('label', 'Size', 'Ln_TChl', 'Temp', 'NO3'),
               c('label', 'Size', 'Ln_TChl', 'Temp', 'Ln_dFe'),
               c('label', 'Size', 'Ln_TChl', 'Temp', 'NO3', 'Ln_dFe')
              )

#k-fold validation
N_fold     <- 10L
Nrep       <- 10L #Number of repetitions
MAE_rm_var <- matrix(NA, nc = length(Xnames) + 1, nr = N_fold * Nrep)
COR_rm_var <- matrix(NA, nc = length(Xnames) + 1, nr = N_fold * Nrep)

Compute_MAE_COR <- function(TRAIN_DATA, TEST_DATA){
  model <- Fitted_model(TRAIN_DATA = TRAIN_DATA)
 
  #Evaluate the test data
  pred.keras <- model %>% 
    predict(TEST_DATA %>% 
              select(names(TRAIN_DATA)) %>%
              select(-label)) 
  
  return(list(MAE = mean(abs(TEST_DATA$label - pred.keras)),
              COR =  cor(    TEST_DATA$label,  pred.keras) ))
}

mm <- as.integer(nrow(Full_data)/N_fold)
for (i in 1:N_fold){
  index    <- (mm*(i-1) + 1):(mm*i)
  test_df  <- Full_data[index, ]
  train_df <- Full_data[-index,]
 
  for (k in 1:length(Xnames)){
 
    varnames <- Xnames[k][[1]]
    train <- train_df %>% select(all_of(varnames))
    test  <-  test_df %>% select(all_of(varnames))
  
    #Mean MAE
    for (j in 1:Nrep){
      
      cff <- Compute_MAE_COR(TRAIN_DATA=train, TEST_DATA=test)
      
      print(paste('i =',i, 'k =', k, 'j = ', j, 'MAE =', round(cff$MAE, 2)))
      print(paste('i =',i, 'k =', k, 'j = ', j, 'COR =', round(cff$COR, 2)))
      
      L <- (i - 1)*Nrep + j
      MAE_rm_var[L, k] <- cff$MAE
    
      #correlation
      COR_rm_var[L, k] <- cff$COR
    
    }
  }
  #Compute the MAE and COR for the full model 
  for (j in 1:Nrep){
    cff <- Compute_MAE_COR(TRAIN_DATA=train_df, TEST_DATA=test_df)
      L <- (i - 1)*Nrep + j
    MAE_rm_var[L, ncol(MAE_rm_var)] <- cff$MAE
    COR_rm_var[L, ncol(COR_rm_var)] <- cff$COR
  }
}

#save MAE_rm_var and COR_rm_var
save(MAE_rm_var, COR_rm_var, file = 'MAE_COR.Rdata')

#Conclusion: Model 4 or 2 is good enough
pdf("Keras_model_comparisons.pdf", width = 8, height = 10)
par(font.lab  = 1,
    family    = "serif",
    mgp       = c(2.2,1,0),
    oma       = c(3,2,2,2),
    mar       = c(4,4,2,.2),
    mfrow     = c(2,2)
   )
#Run a RBD for MAE and COR
MAE <- as.data.frame(t(MAE_rm_var)) %>%
  reshape(direction = 'long', 
          varying = list(1:ncol(.)),
          v.names = 'MAE') %>%
  mutate(id = as.factor(id)) %>%
  mutate(time = as.factor(time))

RBD.MAE <- aov(MAE ~ id + time, data = MAE)
boxplot(MAE ~ id, data = MAE,
        xlab = 'Model ID')
mtext('A)', adj = 0)

plot(TukeyHSD(RBD.MAE, which = 'id'),
     las = 1)
mtext('B)', adj = 0)

COR <- as.data.frame(t(COR_rm_var)) %>%
  reshape(direction = 'long', 
          varying = list(1:ncol(.)),
          v.names = 'COR') %>%
  mutate(id = as.factor(id)) %>%
  mutate(time = as.factor(time))

RBD.COR <- aov(COR ~ id + time, data = COR)
boxplot(COR ~ id, data = COR,
        xlab = 'Model ID',
        ylab = 'Correlation')
mtext('C)', adj = 0)

plot(TukeyHSD(RBD.COR, which = 'id'),
     las = 1)
mtext('D)', adj = 0)

dev.off()
