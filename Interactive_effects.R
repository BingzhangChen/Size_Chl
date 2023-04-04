#Check the effects of single variables  (TChl and Temp) on mean size and size variance
N      <- 20 #Number of discreted points
Chl_b  <- quantile(Data$Ln_TChl, probs = c(0.01, 0.99))
TChls  <- seq_range(Chl_b, n = N, pretty = T)
Temp_b <- quantile(Data$Temp, probs = c(0.01, 0.99))
Temps  <- seq_range(Temp_b, n = N, pretty = T)
NO3_b  <- quantile(Data$NO3, probs = c(0.01, 0.99))
NO3s   <- seq_range(NO3_b, n = N, pretty = T)

B <- Spectra_keras(Temp  = Temps, 
                 Ln_TChl = TChls,
                   NO3   = NO3s,
                 cal.density = F)

d  <- Mean_VAR_size(B$CDF, keep.rep = F) #just averages over each condition
dr <- Mean_VAR_size(B$CDF, keep.rep = T) #with replication

#Run ANOVA to check the percentages of mean size and size variance explained by environmental factors
#Construct new data frame
newD <- B$CDF[,1:(ncol(B$CDF)-1)]
newD$MEAN <- d$CWM
newD$VAR  <- d$VAR

#Add replication into newD
newDr <- dr$CWM
newDr$VAR <- dr$VAR$VAR
newDr$Temp <- NA
newDr$Ln_TChl <- NA
newDr$NO3 <- NA
for (i in 1:Nrep){
  newDr[newDr$time == i, ]$Temp    <- newD$Temp
  newDr[newDr$time == i, ]$NO3     <- newD$NO3
  newDr[newDr$time == i, ]$Ln_TChl <- newD$Ln_TChl
}
newDr <- newDr %>%
  mutate(Temp = as.factor(Temp)) %>%
  mutate(NO3 = as.factor(NO3)) %>%
  mutate(Ln_TChl = as.factor(Ln_TChl))

CWM.aov <- aov(CWM ~ Ln_TChl + Temp + NO3 + 
               Ln_TChl:Temp + Ln_TChl:NO3 + Temp:NO3,
             data = newDr)
summary(CWM.aov)

par(mfrow=c(2,2))
plot(CWM.aov)

VAR.aov <- aov(VAR ~ Ln_TChl + Temp + NO3 + 
               Ln_TChl:Temp + Ln_TChl:NO3 + Temp:NO3,
             data = newDr)
summary(VAR.aov)
  
pdf("Environmental_controls_interaction_4panel.pdf", 
    width = 10, height = 9)
par(font.lab  = 1,
    family    = "serif",
    mgp       = c(2.2,1,0),
    oma       = c(3,2,2,2),
    mar       = c(2.5,2.5,1,.2),
    mfrow     = c(2,2)
   )
oriChl <- c(0.1, 0.2, 0.5, 1, 2, 5)
NO3_names <- c('Low', 'High')

#plot effects of TChl and Temp on mean size and size variance at two levels of NO3
for (i in 1:2){
  if (i == 1){
    m <- 1
  }else if (i == 2){
    m <- length(NO3s)
  }else{
    stop('Incorrect i!')
  }
  d1 <- newD %>%
    filter(NO3 == NO3s[m])
  
  V <- matrix(d1$MEAN, nr = length(TChls), nc = length(Temps))
  image2D(V, TChls, Temps,
          zlim = range(newD$MEAN),
          xaxt = 'n',
          xlab = '',
          ylab = '')
  axis(1, at = log(oriChl), labels = oriChl)
  mtext(paste0(LETTERS[i],') Mean Ln ESD, ', NO3_names[i],' NO3'), adj=0)
  
  V <- matrix(d1$VAR, nr = length(TChls), nc = length(Temps))
  image2D(V, TChls, Temps,
          zlim = range(newD$VAR),
          xaxt = 'n',
          xlab = '',
          ylab = '')
  axis(1, at = log(oriChl), labels = oriChl)
  mtext(paste0(LETTERS[i+1],') Variance of Ln ESD, ',
               NO3_names[i],' NO3'), adj = 0)
}
mtext('Chl (µg/L)', side = 1, outer=T, adj = .5)
mtext('Temperature (ºC)', side = 2, outer=T, adj = .5)
dev.off()
