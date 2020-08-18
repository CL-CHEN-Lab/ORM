# ORM
Bioinformatic tools for the analysis of Optical Replication Mapping (ORM) data to study genome-wide DNA replication program at the single molecule level.

1. AllRawDataRefining_1G2R.jar

Option: 
-B: BNX parent directory
-X: XMP parent directory
-R: RCmap parent directory
-Q: QCmap parent directory
-O: Output directory 
-SampleName: All SampleName seperate by "," 

Notice: 
1.Pleasw pay attention all input files must have same name beside suffix, including .bnx, .xmp, .qcmap, .rcmap
2.There are 2 channels in bionano, here the channel 1 is for green signals for mapping and channel 2 is red signals for labelling

Example:
java -jar AllRawDataRefining_1G2R.jar 
     -B /Volumes/WWT/Final-Version/1905/BNX/ 
     -X /Volumes/WWT/Final-Version/1905/XMP/ 
     -R /Volumes/WWT/Final-Version/1905/RCmap/ 
     -Q /Volumes/WWT/Final-Version/1905/QCmap/ 
     -O /Volumes/WWT/Final-Version/1905/ 
     -SampleName 1905async,1905.FC0
     
