# Optical Replication Mapping (ORM) analysis toolkit 

Bioinformatic tools for the analysis of Optical Replication Mapping (ORM) data to study genome-wide DNA replication program at the single molecule level.


# Introduction

Origin usage in eukaryotic cells has been shown to be remarkably dynamic and flexible: the choice of origins to be activated varies from cell to cell within a population. However, studying the degree of such stochastic variation is a challenging issue, since so far the current high throughput methods can only provide a population-averaged picture.

In order to map replication origins at single molecule level at a high-throughput manner along the human genome, in collaborated with N. Rhind (UMass Med., USA) and D.M. Gilbert (Florida Univ., USA), we have developed a new approach to map early-firing human replication origins by Optical Replication Mapping (ORM). We combined in vivo labeling of origins labelled by fluorescent nucleotides with Bionano visualization and mapping of long individual DNA molecules. Bionano analysis involves fluorescently labeling deproteinated DNA at specific loci (e.g. the sites of the Nt.BspQI restriction enzyme, or more recent developed Direct Label and Stain – DLS labeling), stretching individual molecules in nanochannels and imaging both the replication incorporation tracks and the restriction sites at about 1 kb resolution. The pattern of fluorescent restriction sites allows the molecules and their associated replication tracks to be mapped to the human genome. We identified the fluorescent-nucleotide incorporation segments, called Optical Replication Mapping (ORM) tracks, and mapped them to the genome to get the location of early firing origins within each individual molecule.


# The main function of the toolkit

Common usage (Could be applied on any DNA signal labeling analysis):
      - Calibrate and calculate the labeling signal position	
      - Create GTF like file to visualize signal by IGV
      - Generate all DNA fiber coordinates in bed file format
      - Calculate labeling signal number in Sliding or Adjacent windows, and do normalization by mapped fiber depth in each bin
      - Using Gaussian mixed model to analyze the distance distribution between labelling signals to identify proper cutoff for the signal segementation
      - Perform the segmentation to get the segments with cluster of adjacent signals (in our case, ORM tracks, i.e. DNA replication origins or on-going replication forks), in bed file format

Specific usage for ORM segmentation analysis:
      •	Add replication Fork Direction Index (FDI) for illustrating signal polarity into the ORM segmentation bed file

Specific usage for DNA replication analysis:
      •	Add replication timing (RT) obtained by Repli-Seq or other related techiniques into the ORM bed file
      •	Add replication fork directional (RFD) data obtained by OK-Seq into the ORM bed file
      •	Calculate FDI_RFD and draw enrichment plot around given regions (R)
      •	Based on normalized signal numbers in sliding or adjacent bins using LOESS fit to generate smoothing fire efficiency curves and capture initial zone with peak calling (R)
      •	Draw labeling signaling enrichment in given regions(R)

The majority of the toolkit was developped using JAVA and some fonctions are in R, specified as (R). See the User Manual and our preprint manuscript for detail.


# Authors:

Weitao WANG: weitao.want@curie.fr

Chunlong CHEN: chunlong.chen@curie.fr

Don't hesitate to contact the authors or open an issue for any questions.



# Reference:

Wang W, Klein K, Proesmans K, Yang H, Marchal C, Zhu X, Borrman T, Hastie A, Weng Z, Bechhoefer J*, Chen CL*, Gilbert DM* and Rhind N*. 2020. Genome-Wide Mapping of Human DNA Replication by Optical Replication Mapping Supports a Stochastic Model of Eukaryotic Replication Timing. bioRxiv 2020.08.24.263459. https://www.biorxiv.org/content/10.1101/2020.08.24.263459v1 





