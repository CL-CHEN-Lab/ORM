library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(RColorBrewer)
library(scales)
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)




#######   Input dataframe is the LOESS fitting output  ########
GetInitialZone <- function(Fit,SN,Sign_SN) 
{
  arg2 = data.frame(Fit,SampleNumber=SN)
  IZ = data.frame()
  Fit[which(Fit[,4]<0),4] <- 0
  for(i in 1:23)
  {
    options(scipen = 300)
    chromNames <-  read.table("~/Desktop/enrichment/ChrLength.txt",header=FALSE,sep="\t",comment.char = "#",stringsAsFactors = FALSE)
    chr.name <- chromNames[i,1]
    print(chr.name)
    Chosen <- which(arg2==chr.name)
    
    
    if (length(Chosen)) {
      
      #Read data
      if(grepl("chr", arg2[1,1]))
      {
        chr.name = sub("chr", "Chr", chr.name)
      }
      
      Data.chr <- arg2[arg2[,1]==chr.name,]	
      
      S <- which(Data.chr[,6]=="Peak")
      
      for(t in S)
      {
        
        m=1
        while(Data.chr[(t-m),6]!="Valley" & Data.chr[(t+m),6]!="Valley")
        {
          m=m+1
          if(m>100000)
          {
            print("err")
            break
          }
        }
        
        PeakPosition = floor( (Data.chr[t,2] + Data.chr[t,3])/2)
          
        Delta_Hight=0
        
        if(m!=1)
        {
          Hight1 = Data.chr[t,4] - Data.chr[(t-m),4]
          Hight2 = Data.chr[t,4] - Data.chr[(t+m),4]
          Delta_Hight = min(Hight1,Hight2)
        }
        
        SUM = Data.chr[t,5]
        N=0
        Start=t
        End=t
        Boundry = FALSE
        
        if(Data.chr[t,7]!=0)
        {
          while (SUM < 0.4 ) 
          {
            N=N+1
            #print(paste(N,Start,End,SUM))
            if(Start-1 == 0 & End+1 != nrow(Data.chr))
            {
              End = End+1;
              Boundry = TRUE;
            }
            
            if(Start-1 != 0 & End+1 == nrow(Data.chr))
            {
              Start = Start-1;
              Boundry = TRUE;
            }
            
            if((Data.chr[Start,6]=="Valley") & (Data.chr[End,6]!="Valley") )
            {
              End = End+1;
              Boundry = TRUE;
            }
            
            if((Data.chr[Start,6]!="Valley") & (Data.chr[End,6]=="Valley") )
            {
              Start = Start-1;
              Boundry = TRUE;
            }
            
            if(!Boundry)
            {
              if((Data.chr[Start-1,4] > Data.chr[End+1,4]))
              {
                Start = Start-1
              }
              else
              {
                End = End+1
              }
            }
            
            SUM = sum(Data.chr[Start:End,5])
          }
          
          
          if((End-Start)>0)
          {
            samplenumber=0
            SampleSet <- Data.chr[Start:End , 8]
            for (tmpsign in Sign_SN) {
              if(any(str_detect(SampleSet,tmpsign)))
              {
                samplenumber=samplenumber+1
              }
            }
            
            Chr_IZ = cbind(Data.chr[Start,1:2],Data.chr[End,3],Data.chr[t,4],samplenumber,Delta_Hight,PeakPosition);
            colnames(Chr_IZ)<-c("Chr","start","end","PeakFireEfficiency","ContainSample","Delta_Hight","PeakPosition")
            IZ = rbind(IZ,Chr_IZ)
          }
        }
          
        
      }
    }
  }
  return(IZ)
}

FinalIZ <- GetInitialZone(All0min,SN, Sign)




#######   Input dataframe is the fitting output  ########
GetInitialZone_byBin <- function(Fit) 
{
  Fit[which(Fit[,4]<0),4] <- 0
  Result <- list()
  for(i in 1:23)
  {
    options(scipen = 300)
    chromNames <-  read.table("~/Desktop/enrichment/ChrLength.txt",header=FALSE,sep="\t",comment.char = "#",stringsAsFactors = FALSE)
    chr.name <- chromNames[i,1]
    print(chr.name)
    Chosen <- which(Fit==chr.name)
    
    if (length(Chosen)) {
      
      #Read data
      if(grepl("chr", Fit[1,1]))
      {
        chr.name = sub("chr", "Chr", chr.name)
      }
      
      Data.chr <- Fit[Fit[,1]==chr.name,]	
      Data.chr <- data.frame(Data.chr,Calling=rep(0,nrow(Data.chr)))
      loop = nrow(Data.chr)-1
      Sign = rep("NA",loop+1)
      Window=50;
      
      for(i in 2:loop)
      {
        Middle=Data.chr$V4[i]
        Previous=Data.chr$V4[i-1]
        Latter=Data.chr$V4[i+1]
        
        if(Previous < Middle & Latter > Middle)
        {
          Sign[i]="Up"
        }
        
        if(Previous <= Middle & Latter <= Middle)
        {
          Sign[i]="Peak"
        }
        
        if(Previous >= Middle & Latter >= Middle)
        {
          Sign[i]="Valley"
        }
        
        if(Previous > Middle & Latter < Middle)
        {
          Sign[i]="Down"
        }
      }
      
    
      
      Sign[length(Sign)] <-"Valley"
      Sign[1]<-"Valley"
      Data.chr <- data.frame(Data.chr,Sign=Sign)
      Select <- which(Sign=="Peak")
      options(scipen = 200)
      
      Data.chr <- data.frame(Data.chr,IZ=rep(0,nrow(Data.chr)))
      
      Output <- matrix(NA,nrow=length(Select),(Window*2+1))
      N=1
      
      for(i in Select)
      {
        j=1
        
        while(Sign[(i-j)]!="Valley" & Sign[(i+j)]!="Valley")
        {
          j=j+1
          if(j>100000)
          {
            print("err")
            break
          }
        }
        
        
        if(j!=1)
        {
          SignalSUM = sum(Data.chr[(i-j):(i+j),4])
          Out = rep(NA, Window*2+1)
          Data.chr[(i-j):(i+j),5] <- Data.chr[(i-j):(i+j),4]/SignalSUM
          
          SUM = Data.chr[i,5]
          N=1
          t=i
          Start=t
          End=t
          Boundry=FALSE
          
          while (SUM < 0.4 ) 
          {
            
            if(Start-1 == 0 & End+1 != nrow(Data.chr))
            {
              End = End+1;
              Boundry = TRUE;
            }
            
            if(Start-1 != 0 & End+1 == nrow(Data.chr))
            {
              Start = Start-1;
              Boundry = TRUE;
            }
            
            if((Data.chr[Start,6]=="Valley") & (Data.chr[End,6]!="Valley") )
            {
              End = End+1;
              Boundry = TRUE;
            }
            
            if((Data.chr[Start,6]!="Valley") & (Data.chr[End,6]=="Valley") )
            {
              Start = Start-1;
              Boundry = TRUE;
            }
            
            if(!Boundry)
            {
              if((Data.chr[Start-1,4] > Data.chr[End+1,4]))
              {
                Start = Start-1
              }
              else
              {
                End = End+1
              }
            }
            
            SUM = sum(Data.chr[Start:End,5])
          }
          
          Data.chr[Start:End,]$IZ=1
          
        }
      }
    }
    Result <-rbind(Result,Data.chr)
  }
  
  return(Result)
}




#####  Generate IZ-Bin #########
All0min <- read.table("/Users/wwang/Desktop/Final-Version/Final_Initialzone/CallPeak/Treat_Peak/Sliding/Fitting/All0min.bedgraph")
S1708 <- read.table("/Users/wwang/Desktop/Final-Version/Final_Initialzone/CallPeak/Treat_Peak/Sliding/Fitting/1708.bedgraph")
S1802 <- read.table("/Users/wwang/Desktop/Final-Version/Final_Initialzone/CallPeak/Treat_Peak/Sliding/Fitting/1802.bedgraph")
S1807 <- read.table("/Users/wwang/Desktop/Final-Version/Final_Initialzone/CallPeak/Treat_Peak/Sliding/Fitting/1807.bedgraph")
S1905 <- read.table("/Users/wwang/Desktop/Final-Version/Final_Initialzone/CallPeak/Treat_Peak/Sliding/Fitting/1905.bedgraph")


Sample <- list(S1708,S1802,S1807,S1905,All0min)
Names <- c("1708","1802","1807","1905","All0min")
for (i in 5:5) 
{
  Sample[[i]][which(Sample[[i]][,4]<0),4] <- 0
  
  Out_table = GetInitialZone_byBin(Sample[[i]])
  write.table(Out_table,paste("/Users/wwang/Desktop/Final-Version/Final_Initialzone/IZ_Update/SignalRatio/Sliding/",Names[i],"_SignalRatio.bed",sep=""),quote = FALSE,row.names = FALSE,col.names = FALSE, sep="\t")  
  
}





######  Get replicate IZ ####
#
# Names <- c("1708","1802","1807","1905")
# for(k in 1:4)
# {
#   Path = paste("/Users/wwang/Desktop/Final-Version/Final_Initialzone/IZ_Update/SignalRatio/Sliding/",Names[k],"_SignalRatio.bed",sep="")
#   arg2 = read.table(Path)
#   IZ = list()
#   print(Names[k])
#   for(i in 1:23)
#   {
#     options(scipen = 300)
#     chromNames <-  read.table("~/Desktop/enrichment/ChrLength.txt",header=FALSE,sep="\t",comment.char = "#",stringsAsFactors = FALSE)
#     chr.name <- chromNames[i,1]
#     print(chr.name)
#     Chosen <- which(arg2==chr.name)
# 
#     if (length(Chosen)) {
# 
#       #Read data
#       if(grepl("chr", arg2[1,1]))
#       {
#         chr.name = sub("chr", "Chr", chr.name)
#       }
# 
#       Data.chr <- arg2[arg2[,1]==chr.name,]
#       colnames(Data.chr) <- c("chr","start","end","FireEfficiency","FirePercentage","Sign","IZ")
#       FireEfficiency = Data.chr[which(Data.chr$Sign=="Peak" & Data.chr$IZ!=0),4]
# 
#       S <- which(Data.chr$IZ==1)
#       Data.chr <- Data.chr[S,]
# 
#       Chr_IZ <- Data.chr%>%makeGRangesFromDataFrame()%>%GenomicRanges::reduce(min.gapwidth=10)%>%as.data.frame()
#       Chr_IZ <- data.frame(Chr_IZ,FireEfficiency=FireEfficiency)
#       IZ = rbind(IZ,Chr_IZ)
# 
#     }
#   }
# 
#   write.table(IZ,paste("/Users/wwang/Desktop/Final-Version/Final_Initialzone/IZ_Update/IZ_Replicate/",Names[k],"_Update_IZ_Sliding.bed",sep=""),quote = FALSE,row.names = FALSE,col.names = FALSE, sep="\t")
# 
# }




####   Call All0min IZ ####
S1708 <- read.table("/Users/wwang/Desktop/Final-Version/Final_Initialzone/IZ_Update/SignalRatio/Sliding/1708_SignalRatio.bed")
S1802 <- read.table("/Users/wwang/Desktop/Final-Version/Final_Initialzone/IZ_Update/SignalRatio/Sliding/1802_SignalRatio.bed")
S1807 <- read.table("/Users/wwang/Desktop/Final-Version/Final_Initialzone/IZ_Update/SignalRatio/Sliding/1807_SignalRatio.bed")
S1905 <- read.table("/Users/wwang/Desktop/Final-Version/Final_Initialzone/IZ_Update/SignalRatio/Sliding/1905_SignalRatio.bed")

IZ <- list(S1708,S1802,S1807,S1905)
Sign <- c("a","b","c","d")

SampleNumber <- rep(0,nrow(IZ[[1]]))

for(i in 1:4)
{
  Tmp <- IZ[[i]]
  Tmp[,8] <- "0"
  Tmp[which(IZ[[i]]$V7 > 0),8] <- Sign[i]
  SampleNumber <- paste(SampleNumber,Tmp$V8,sep = "")
}

library(stringr)
SN <- str_replace_all(SampleNumber, "0", "")
All0min <- read.table("/Users/wwang/Desktop/Final-Version/Final_Initialzone/IZ_Update/SignalRatio/Sliding/All0min_SignalRatio.bed")
FinalIZ <- GetInitialZone(All0min,SN, Sign)
colnames(FinalIZ)<-c("Chr","start","end","PeakFireEfficiency","ContainSample","PeakPosition")
write.table(FinalIZ,paste("/Users/wwang/Desktop/Final-Version/Final_Initialzone/IZ_Update/IZ/All0min.bed",sep=""),quote = FALSE,row.names = FALSE,col.names = FALSE, sep="\t") 

