
EnrichmentPlot <-  function(arg1,arg2,Plotname,colum,resolution,Window,total_Tittle,SampleName,REFpath)
{
  
  chromNames <-  read.table(REFpath,header=FALSE,sep="\t",comment.char = "#",stringsAsFactors = FALSE)
  Gap <- c() 
  
  arg1.Center <- matrix(NA,nrow=nrow(arg1),ncol=Window*2+1)
  
  N <- 0
  
  for (i in 1:23)   {
    
    chr.name <- chromNames[i,1]
    print(chr.name)
    arg1[,1] <- gsub('Chr', 'chr', arg1[,1])
    arg1[,1] <- gsub('23', 'X', arg1[,1])
    arg1[,1] <- gsub('24', 'Y', arg1[,1])
    
    
    Chosen <- which(arg1[,1]==chr.name)

    
    if (length(Chosen)) {
      
      
      Data.chr <- arg2[arg2[,1]==chr.name,]		
      
      #Compute
      for (j in 1:length(Chosen)) {
        N <- N+1
        
        Center <- (arg1[Chosen[j],2]+arg1[Chosen[j],3])/2
        
        
        if (Center < 1000) {
          Select <- 0
        } else {  
          Select <- which(Data.chr[,2] <= Center & Data.chr[,3] >= Center)
        }
        
        if (length(Select)==1) {
          if (Select > Window) {
            Data.Chosen <- Data.chr[(Select-Window):(Select+Window),colum]
          } else {
            Data.Chosen <- c(rep(NA,Window-Select+1),Data.chr[1:(Select+Window),colum])
          }
          
          tmpline = as.numeric(as.character(Data.Chosen)) #Remove NA element
          #tmpline = as.numeric(as.character(replace(Data.Chosen, Data.Chosen == 0, NA)))  # Remove 0 element
          arg1.Center[N,] <- tmpline
          
          #if(length(which(tmpline>1)) > 0)
          #{
          # print(tmpline)
          #}
          
          
        }
      }
    }
  }
  

  title = paste(Plotname," record number:",nrow(arg1.Center))
  #pdf(file=PlotOutput, bg="transparent")
  
  Draw <- colMeans(arg1.Center,na.rm=T)
  
  plot(seq(Window*(-1),Window,by=resolution), Draw, type="b",main = title,cex=0.1,ylim = c(-0.4,0.4),xlab = "Bin around Tpeak(unit 1kb)",ylab ="FDI_RFD")
  
  Out_table <- data.frame(Xvalue = seq(Window*(-1),Window,by=resolution),Yvalue = Draw)
  
  #write.table(Out_table,"~/Desktop/TpeakAround/SignalNumber.txt",quote = FALSE,row.names = FALSE,col.names = T, sep="\t")  
  
  
  mtext(total_Tittle, outer = TRUE, cex = 1.5)
  #dev.off()
  return(Out_table)
}





Step = 0.25
S50_Range = c()

for(i in 1:4)
{
  S50_Range[i] <- paste(Step * (i-1) , "~" , Step * i,sep="")
}




#AllTpeak_OKseq <- read.table("/Users/wwang/Desktop/Final-Version/Tpeak/Step0.25/OKseq_ClassifyByRawS50/AllTpeak.bed")
#TpeakName = "OKseq"
AllTpeak_OKseq <- read.table("/Users/wwang/Desktop/Final-Version/Tpeak/Step0.25/UCSC_Encode_Repliseq_hg19_Peak/All_Tpeak")
TpeakName = "UCSC"

Early_OKseq_S <- AllTpeak_OKseq[which(AllTpeak_OKseq$V5<= Step),]
MidEarly_OKseq_S <- AllTpeak_OKseq[which(AllTpeak_OKseq$V5 <= Step*2 & AllTpeak_OKseq$V5 > Step),]
MidLate_OKseq_S <- AllTpeak_OKseq[which(AllTpeak_OKseq$V5 <= Step*3 & AllTpeak_OKseq$V5 > Step*2),]
Late_OKseq_S <- AllTpeak_OKseq[which(AllTpeak_OKseq$V5 > Step*3),]
S50Classification <- list(Early_OKseq_S, MidEarly_OKseq_S, MidLate_OKseq_S, Late_OKseq_S)


path = "/Users/wwang/Desktop/ORM/FDI_RFD/Test/Output.bed"

Sum = read.table(path)


par(mfrow=c(1,4),oma = c(0, 0, 2, 0))


for(i in 1:4)
{
  Tittle = paste(" S50:" , S50_Range[i],"Record number:", nrow(S50Classification[[i]]))
  EnrichmentPlot(S50Classification[[i]],Sum ,Tittle,7  ,1,200, paste("ORM_FDI_RFD in AllFDI in 1kb around early" ,TpeakName,"Tpeak"), S50_Range[i],"~/Desktop/enrichment/ChrLength_hg19.txt")
}



