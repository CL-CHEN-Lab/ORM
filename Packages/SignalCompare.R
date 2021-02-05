
#Names = c("0min","20min","30min","45min","60min","90min")
#Names = c("0min0","0min1","0min2","5min","10min")

Names = c("1807","1802")

for(i in 1:length(Names))
{
  In_Green = read.table(paste("/Volumes/WWT/Final-Version/StrangeHotDot/Step3_HotDotFilter/FilterTXT/In/",Names[i],"_Green.txt",sep=""))
  Out_Green = read.table(paste("/Volumes/WWT/Final-Version/StrangeHotDot/Step3_HotDotFilter/FilterTXT/Out/",Names[i],"_Green.txt",sep=""))
  In_Red = read.table(paste("/Volumes/WWT/Final-Version/StrangeHotDot/Step3_HotDotFilter/FilterTXT/In/",Names[i],"_Red.txt",sep=""))
  Out_Red = read.table(paste("/Volumes/WWT/Final-Version/StrangeHotDot/Step3_HotDotFilter/FilterTXT/Out/",Names[i],"_Red.txt",sep=""))
  
  In_Green = data.frame(SNR = In_Green$V1,Type=rep("Green_InHotDotRegion",nrow(In_Green)))
  Out_Green = data.frame(SNR = Out_Green$V1,Type=rep("Green_OutHotDotRegion",nrow(Out_Green)))
  
  In_Red = data.frame(SNR = In_Red$V1,Type=rep("Red_InHotDotRegion",nrow(In_Red)))
  Out_Red = data.frame(SNR = Out_Red$V1,Type=rep("Red_OutHotDotRegion",nrow(Out_Red)))
  
  AllData = rbind(In_Red,Out_Red)
  library(ggplot2)
  ggplot(AllData, aes(SNR,color = Type)) + geom_density()+xlim(0,50)
  ggsave(paste("/Volumes/WWT/Final-Version/StrangeHotDot/Step3_HotDotFilter/FilterTXT/Plots/Red/",Names[i],"_Red_line.pdf",sep = ""))
  

  AllData = rbind(In_Green,Out_Green)
  cname = 
  library(ggplot2)
  ggplot(AllData, aes(SNR,color = Type)) + geom_density()+xlim(0,500)
  ggsave(paste("/Volumes/WWT/Final-Version/StrangeHotDot/Step3_HotDotFilter/FilterTXT/Plots/Green/",Names[i],"_Green_line.pdf",sep = ""))
}