

###### reference example: https://www.jianshu.com/p/710feceee42f ########
FitCurve <- function(arg1,window,overlaplength,Signal_column) # window must bigger than overlaplength
{
  Allceiling=nrow(arg1)
  Allcontrol = rep(0,Allceiling)
  
  for (i in 1:24) 
  {
    
    
    if(grepl("chr", arg1[1,1]))
    {
      Select <- which(arg1[,1]==paste("chr",i,sep = ""))
      #chr.name = sub("chr", "Chr", chr.name)
    }
    else
    {
      Select <- which(arg1[,1]==paste("Chr",i,sep = ""))
    }
    
    print(paste("Chr",i,sep = ""))
    
    
    Chosen <- arg1[Select,]
    Chosen = Chosen[-1,]
    ceiling = nrow(Chosen)
    control = rep(0,ceiling)
    transition1=c()
    transition2=c()
    weight=seq(from=0,to=1,by=(1/overlaplength))
    for(i in seq(from = 1, to = ceiling,by = window))
    {
      if(i==1)
      {
        GenomicPosition <- (Chosen[i:(i+window-1+overlaplength),2]+Chosen[i:(i+window-1+overlaplength),3])/2  #Left overlapped overlaplength bins as transition period
        SignalStrenth <- Chosen[i:(i+window-1+overlaplength),Signal_column]
        IZ_Chr <- data.frame(GP=GenomicPosition,Signal=SignalStrenth)
        
        #Loess_R <- loess(Signal~GP,IZ_Chr,parametric=IZ_Chr$Signal,family = c("gaussian"))
        Loess_R <- loess(Signal~GP,IZ_Chr,span=0.75, degree=2)
        
        
        smooth_vals <- predict(Loess_R, IZ_Chr$GP)
        control[1:window] = smooth_vals
        transition1 = smooth_vals[(window+1):(window+overlaplength)]
      }
      else
      {
        if( (i+window+overlaplength-1) < ceiling)
        {
          GenomicPosition <- (Chosen[(i-overlaplength):(i+window-1+overlaplength),2]+Chosen[(i-overlaplength):(i+window-1+overlaplength),3])/2
          SignalStrenth <- Chosen[(i-overlaplength):(i+window-1+overlaplength),Signal_column]
          IZ_Chr <- data.frame(GP=GenomicPosition,Signal=SignalStrenth)
          #Loess_R <- loess(Signal~GP,IZ_Chr,parametric=IZ_Chr$Signal,family = c("gaussian", "symmetric"))
          Loess_R <- loess(Signal~GP,IZ_Chr,span=0.75, degree=2)
          smooth_vals <- predict(Loess_R, IZ_Chr$GP)
          
          transition2 = smooth_vals[1+overlaplength:(overlaplength*2-1)]
          
          # start = transition1[1]
          # end = transition2[overlaplength]
          # step = (end-start)/overlaplength
          # control[i:(i+overlaplength-1)] = seq(from=start,to=end,by=step)
          
          ControlTmp =c()
          for (k in 1:overlaplength) {
            ControlTmp[k] = transition1[k]*(1-weight[k]) + transition2[k]*weight[k]
          }
          
          control[i:(i+overlaplength-1)] = ControlTmp
          
          control[(i+overlaplength):(i+window-1)] = smooth_vals[(1+overlaplength*2):(window+overlaplength)]  
          transition1 = smooth_vals[(window+overlaplength+1):(window+overlaplength*2)]
        }
      }
    }
    
    Allcontrol[Select] = control
  
  }
  
  Result <- data.frame(arg1[1:3],fitcurve=Allcontrol)
  return(Result)
}


RFD <- read.table('/Users/wwang/Desktop/Final-Version/RFD_InWindow/input/EDC_OKseq_RFD_RawData.bedgraph',header = F,sep='\t',stringsAsFactors=FALSE)
RFD[RFD$V4==-2,4] <- NA
RFD[RFD$V1=="ChrX",1] <- "Chr23"
RFD[RFD$V1=="ChrY",1] <- "Chr24"
Fit <- FitCurve(RFD,20,10,4)
Fit[is.na(Fit$fitcurve),4] <- 0

write.table(Fit, "/Users/wwang/Desktop/Final-Version/RFD_InWindow/Smoothing_RFD.bedgraph",quote = F,row.names = F,col.names = F,sep = "\t")




