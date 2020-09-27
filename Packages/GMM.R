#update.packages(c("mixtools","sBIC","mclust"))

library("sBIC")
#library(mclust, quietly=TRUE)
library(mixtools)
library(ggplot2)
library(tidyverse)

Names <- c("0min1","0min2","0min0","5min","10min")

Input_ParentPath = "/Volumes/WWT/Final-Version/1802/Distance/BetweenSignal/"
Output_ParentPath = "/Volumes/WWT/Final-Version/1802/Distance/BetweenSignal/"
xlim = c(8,20)

for(i in 1:5)
{
  tmp = read.table(paste(Input_ParentPath, Names[i],".bed",sep = '') , sep = '\t', header = F, stringsAsFactors=FALSE)
  observations <- tibble(value = log(tmp$V2,2))
  pdf(file= paste(Output_ParentPath,Names[i],".pdf",sep = ''))
  x = observations$value
  my_mix <- normalmixEM(x,k=3)
  
  
  Plot <- ggplot(observations, aes(x = value)) +
    geom_histogram(binwidth = 0.05) +
    xlab("distence log2 value between signal")+
    xlim(xlim)+
    ggtitle("Distance between Signal")+
    mapply(
      function(mean, sd, lambda, n, binwidth) {
        stat_function(
          fun = function(x) {
            (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
          }
        )
      },
      mean = my_mix[["mu"]], #mean
      sd = my_mix[["sigma"]], #standard deviation
      lambda = my_mix[["lambda"]], #amplitude
      n = length(observations$value), #sample size
      binwidth = 0.06 #binwidth used for histogram
    )
  print(Plot)
  dev.off()
}
