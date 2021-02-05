
#############  1802 solo signal  ###############
T0min <- read.table("/Users/wwang/Desktop/Final-Version/Seg_AddSolo/Raw/Bed/1802/1802_0min.bed")
T5min <- read.table("/Users/wwang/Desktop/Final-Version/Seg_AddSolo/Raw/Bed/1802/5min.bed")
T10min <- read.table("/Users/wwang/Desktop/Final-Version/Seg_AddSolo/Raw/Bed/1802/10min.bed")

Solo_0min = T0min[which(T0min$V5==1),]
Solo_5min = T0min[which(T5min$V5==1),]
Solo_10min = T0min[which(T10min$V5==1),]

write.table(Solo_0min,"/Users/wwang/Desktop/Final-Version/Strange_Region/SoloEnrich/1802_0min.bed",quote = FALSE,row.names = FALSE,col.names = FALSE, sep="\t") 
write.table(Solo_5min,"/Users/wwang/Desktop/Final-Version/Strange_Region/SoloEnrich/1802_5min.bed",quote = FALSE,row.names = FALSE,col.names = FALSE, sep="\t") 
write.table(Solo_10min,"/Users/wwang/Desktop/Final-Version/Strange_Region/SoloEnrich/1802_10min.bed",quote = FALSE,row.names = FALSE,col.names = FALSE, sep="\t")


#############  1807 solo signal  ###############
T0min <- read.table("/Users/wwang/Desktop/Final-Version/Seg_AddSolo/Raw/Bed/1807/1807_0min.bed")
T20min <- read.table("/Users/wwang/Desktop/Final-Version/Seg_AddSolo/Raw/Bed/1807/20min.bed")
T30min <- read.table("/Users/wwang/Desktop/Final-Version/Seg_AddSolo/Raw/Bed/1807/30min.bed")
T45min <- read.table("/Users/wwang/Desktop/Final-Version/Seg_AddSolo/Raw/Bed/1807/45min.bed")
T60min <- read.table("/Users/wwang/Desktop/Final-Version/Seg_AddSolo/Raw/Bed/1807/60min.bed")
T90min <- read.table("/Users/wwang/Desktop/Final-Version/Seg_AddSolo/Raw/Bed/1807/90min.bed")


Solo_0min = T0min[which(T0min$V5==1),]
Solo_20min = T20min[which(T20min$V5==1),]
Solo_30min = T30min[which(T30min$V5==1),]
Solo_45min = T45min[which(T45min$V5==1),]
Solo_60min = T60min[which(T60min$V5==1),]
Solo_90min = T90min[which(T90min$V5==1),]


write.table(Solo_0min,"/Users/wwang/Desktop/Final-Version/Strange_Region/SoloEnrich/1807_0min.bed",quote = FALSE,row.names = FALSE,col.names = FALSE, sep="\t") 
write.table(Solo_20min,"/Users/wwang/Desktop/Final-Version/Strange_Region/SoloEnrich/1807_20min.bed",quote = FALSE,row.names = FALSE,col.names = FALSE, sep="\t") 
write.table(Solo_30min,"/Users/wwang/Desktop/Final-Version/Strange_Region/SoloEnrich/1807_30min.bed",quote = FALSE,row.names = FALSE,col.names = FALSE, sep="\t")
write.table(Solo_45min,"/Users/wwang/Desktop/Final-Version/Strange_Region/SoloEnrich/1807_45min.bed",quote = FALSE,row.names = FALSE,col.names = FALSE, sep="\t") 
write.table(Solo_60min,"/Users/wwang/Desktop/Final-Version/Strange_Region/SoloEnrich/1807_60min.bed",quote = FALSE,row.names = FALSE,col.names = FALSE, sep="\t") 
write.table(Solo_90min,"/Users/wwang/Desktop/Final-Version/Strange_Region/SoloEnrich/1807_90min.bed",quote = FALSE,row.names = FALSE,col.names = FALSE, sep="\t")

#############  Example to set the bin size as 300 bases ###############

Strange_Region <- read.table("/Users/wwang/Desktop/Final-Version/Strange_Region/Annotation/Strange_Region.bed")
Mid = (Strange_Region$V3+Strange_Region$V2)/2

for(i in 1 : length(Mid)-1)
{
  print(Mid[i+1]-Mid[i])
}

#############  Then merge all solo segments and use bedtools map All_Solo.bed to Genome_300b.bed to Get Map.bed  ####################

#############  Then run java GetMappedCount to get Count_300.bed  ####################

library("ggplot2")
Count<-read.table("/Users/wwang/Desktop/Final-Version/Strange_Region/SoloEnrich/Count_300.bedgraph")
Count1 <- Count[which(Count$V4>20),]
colnames(Count1) <- c("Chr","Start","End","SignalNumber")

ggplot(Count1, aes(log(SignalNumber,2)))+ geom_histogram(binwidth = 0.1)+geom_vline(xintercept=log(20,2),color = "orange", size=1.5)

write.table(Count1,"/Users/wwang/Desktop/Final-Version/Strange_Region/SoloEnrich/CountBigger20_300.bed",quote = FALSE,row.names = FALSE,col.names = FALSE, sep="\t") 



