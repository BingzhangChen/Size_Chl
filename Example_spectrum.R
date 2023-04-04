library(grid)
#Predict the cumulative size proportion at different Chl, Temp, and size
B <- Spectra_keras(Temp    = median(Data$Temp), 
                   Ln_TChl = median(Data$Ln_TChl),
                   NO3     = median(Data$NO3)
                   ) 

#Try ggplot and smoothing
CDF0 <- B$CDF$data[[1]]

#Compute median and 2.5% and 97.5% quantiles
# Create a text
#grob <- grobTree(textGrob("A", x=0.01,  y=0.98, hjust=0,
#                          gp=gpar(col=1, fontsize=12, fontface="bold")))

cdf1 <- as.data.frame(CDF0) %>%
  reshape(direction = 'long', 
          varying = list(1:(ncol(.)-1)),
          v.names = 'CDF') %>%
  select(-id) 

gamma1 <- gamm(CDF ~ s(Size), family = gaussian(), data = cdf1,
               random = list(time=~1))
  ggplot(mapping = aes(x = Size, y = CDF)) +
    geom_smooth() +
    scale_x_log10() +
    xlab('') +
    ylab('Cumulative probability') +
    theme_light()

#check whether the sum of probabilities equals to 1
#sum(newD[,10] * dsize)

#Try ggplot and smoothing
PDF0 <- B$PDF$data[[1]]

pdf_plot <- as.data.frame(PDF0) %>%
  reshape(direction = 'long', 
          varying = list(1:(ncol(.)-1)),
          v.names = 'PDF') %>%
  select(-id, -time) %>%
  ggplot(mapping = aes(x = Size, y = PDF)) +
    geom_smooth() +
    scale_x_log10() +
    coord_cartesian(ylim = c(0, .20)) +  #Do not drop data
    xlab('Size (µm)') +
    ylab('Density') +
    theme_light()

#arrange them in a single pdf file
pdf("Keras_size_spectrum_example.pdf", width = 5, height = 9)
ggarrange(cdf_plot, pdf_plot, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)
dev.off()
