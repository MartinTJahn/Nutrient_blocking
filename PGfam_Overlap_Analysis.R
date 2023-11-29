# Code to calculate PGfam inhibitory score

#'Preparation: 
#'Define Genomes of Interest as set in PATRIC (n=50 + Pathogens)
#'Download Genome Metadata from Patric > data/PATRIC_metadata.xlsx (add species identifiers, matching Experimental data sheets)
#'Prepare Experimental data sheets (Community_set =  Members of community matching Metadata sheet, also includes growth data) 
#'Download pan genome table from PATRIC > PATRIC_PGfams.csv 


# Set env

setwd("Path to Nutrient_blocking" )
library(plyr)
library(readxl)
library(tidyverse)
library(stringr)
library(data.table)  
library(ggpubr)

# Set vars

PATHOGEN = "Klebsiella" # Can be [Salmonella,Klebsiella] # Pathogen to analyse 

# Set paths

if (PATHOGEN == "Salmonella"){
  pathoLink = "data/annotations/Salmonella_216597.6.PATRIC.pathogen.tab.tsv" 
  comLink = "data/Comunities_v2_Salmonella_EB.xlsx"
  print("Processing Salmonella")
} else if (PATHOGEN == "Klebsiella"){
  pathoLink = "data/annotations/Klebsiella_1162296.3.PATRIC.pathogen.tab.tsv"
  comLink = "data/Comunities_v2_KlebsiellaSet.xlsx"
  print("Processing Klebsiella")
} else {
  print("PATHOGEN not in the list")
}

genomeMetadata = "data/PATRIC_metadata.xlsx" # Genome metadata from PATRIC (thank you!)
pgfamMetadata = "data/PATRIC_PGfams.csv" # PGfam medatata for Set from PATRIC

# Load

CommunitySets <- read_excel(comLink, sheet = "Formated") # Load community sets
genome_meta <- read_excel(genomeMetadata) # Load genome info table
pan = read.csv(pgfamMetadata) # Load cluster metadata for set
colnames(pan) = paste0("PGfam_",colnames(pan))

# Download annotations from PATRIC (optional, see data/annotations for original files)

## Community

genome_meta_Com50 = subset(genome_meta , !genome_meta$StrainID %in% c("EcoHS","EcoMG1655", "EcoZ1269","EcoZ1331")) # focus on EcoIAI1 as community member

#for (i in 1:nrow(genome_meta_Com50)){
#  print(paste0("Downloaded:   ",i, "  of  ",nrow(genome_meta_Com50) ))
#  url = paste0("ftp://ftp.bvbrc.org/genomes/", genomeID, "/", genomeID, "*.features.tab")
#  genomeID = genome_meta_Com50$`Genome ID`[i]
#  download.file(url, destfile = paste0("data/annotations/", genomeID,".features.tab"), quiet = F, mode = "wb",method = "wget")
#} 

## Pathogens

#Purl.Salmo = paste0("ftp://ftp.bvbrc.org/genomes/", "216597.6", "/", "216597.6", "*.features.tab")
#download.file(Purl.Salmo, destfile = "data/annotations/Salmonella_216597.6.PATRIC.pathogen.tab.tsv", quiet = F, mode = "wb",method = "wget")
#Purl.Salmo = NULL
#Purl.Klebs = paste0("ftp://ftp.bvbrc.org/genomes/", "1162296.3", "/", "1162296.3", "*.features.tab")
#download.file(Purl.Klebs, destfile = "data/annotations/Klebsiella_1162296.3.PATRIC.pathogen.tab.tsv", quiet = F, mode = "wb",method = "wget")
#Purl.Klebs  = NULL


# Combine Community annotation from PATRIC
annotation <- data.frame() 
for (i in list.files(pattern = ".features.tab$",path = "data/annotations/", recursive = TRUE, full.names = T)) {
  print(i)
  tmp = as.data.frame(fread(i, select = c("genome_id","genome_name","feature_type","patric_id","annotation","pgfam_id"), colClasses = 'character'))
  annotation <- rbind(annotation, tmp ) 
}
  annotation = annotation %>% filter(genome_id != "genome_id") # remove headers
  length(unique(annotation$genome_id)) # sanity check: # genomes in dataset n = 50


# Cluster stats for Community set 
table(annotation$feature_type) # print available annotation types
annotation <- subset(annotation, feature_type == "CDS") # focus on CDS for analyses 
table(annotation$pgfam_id != "") # PGFams annotation rate: 153,771 of 170,334 = 90.28 % (for PATRIC genome annotation versions as of Sept 2023, see: data/annotations/)

mean(pan$PGfam_Genomes) # PGfams covered by on average 2.44 genomes in set
ggplot(pan, aes(x = PGfam_Genomes)) + geom_histogram(size = 100, bins = 50, binwidth = .5, )  + theme_bw() + geom_vline(aes(xintercept=mean(pan$PGfam_Genomes)), linetype="dashed") # Plot distrition, genomes per pgfam

# Clean & merge data
annotation_select <- subset(annotation, feature_type == "CDS" & pgfam_id != "") # subset CDS with PGFam annotation
genome_meta_red = genome_meta[,c(1:5)] # focus on essentials 
table(genome_meta$`Genome Name` %in% annotation_select$genome_name) # sanity check if annotation files represent all set genomes
annotation_select_meta = merge(genome_meta_red , annotation_select, by.x = "Genome Name", by.y = "genome_name") # Map annotation to genome meta
unique(annotation_select_meta$`Genome Name`) # double check strain list


# Load pathogen PGFAM annotation
PathogenPGF = read.csv(pathoLink, sep="\t")

PathogenPGF <- subset(PathogenPGF, feature_type == "CDS") # subset pathogen CDS with PGfam annotation
table(PathogenPGF$pgfam_id != "") / sum(table(PathogenPGF$pgfam_id == "")) * 100 # Annotation Rate Pathogen

PathogenPGF <- subset(PathogenPGF, feature_type == "CDS" & pgfam_id != "") # subset pathogen CDS with PGfam annotation
PathogenPGFList = unique(PathogenPGF$pgfam_id) # pathogen PGfam set
table(PathogenPGFList %in% annotation_select$pgfam_id) # matched PGfams by any of the 50 symbiont strains


# Filter releavant communities to test
SET = subset(CommunitySets , CommunitySets$paper == "T") 

# Calculate Pathogen community overlap for each community

foo <- data.frame()
com.annot = NULL

for (i in 1:nrow(SET)){
  # gather community features
  strain.list = strsplit(SET$Community_set, split = " ")[[i]] 
  nstrains = length(strain.list)
  CommunityID = as.character(SET[i,"CommunityID"])
  d1 = as.numeric(SET[i,"Day1"]) # d1 pathogen CFU/ml
  d2 = as.numeric(SET[i,"Day2"]) # d2 pathogen CFU/ml
  ecoCom = as.character(SET[i,"Eco_community"])# does community contain E. coli?
  # gather community annotations
  com.annot = subset(annotation_select_meta , annotation_select_meta$StrainNr %in% strain.list) # focal strain annotation
  ComPGFAM = unique(com.annot$pgfam_id) # list unique community PGfam-clusters
  nCluster = length(ComPGFAM) # count uniqe community PGfam-clusters
  
  # calculate pathogen - community PGfams overlap
  PathoOverlap = PathogenPGFList %in% ComPGFAM # TRUE = Overlap, FALSE = no overlap
  OverlapCluster = table(PathogenPGFList %in% ComPGFAM)["TRUE"] # covered
  DifferCluster = table(PathogenPGFList %in% ComPGFAM)["FALSE"] # not covered
  OverlapProp = OverlapCluster / (OverlapCluster + DifferCluster)
  
  # populate output table
    
  foo[i,"CommunityID"] = CommunityID
  foo[i,"nstrains"] = nstrains
  foo[i,"nCluster"] = nCluster
  foo[i,"OverlapProp"] = OverlapProp
  foo[i,"ecoCom"] = ecoCom
  
  foo[i,"d1"] = d1
  foo[i,"d2"] = d2
  }

# classify pathogen itself as 100 % overlap  
foo$OverlapProp[foo$CommunityID == "SL1344 itself"] <- 1
foo$OverlapProp[foo$CommunityID == "Klebsiella itself"] <- 1


# Plot 

# diversity vs overlap - continuous 
foo %>% filter(CommunityID != "SL1344 itself") %>% filter(CommunityID != "Klebsiella itself") %>%
ggscatter( y = "OverlapProp", x = "nstrains", color = "ecoCom", 
          palette = c("black","green"),
          ylab = "Cluster Overlap (%)", xlab = "Strains (n)")  

# diversity vs overlap - categorical 
df_filtered <- foo %>%
  filter(CommunityID != "SL1344 itself") %>%
  filter(CommunityID != "Klebsiella itself") %>%
  filter(nstrains %in% c(1, 2, 3, 5, 9, 10, 49, 50)) %>%
  mutate(nstrains = factor(nstrains, levels = c(1, 2, 3, 5, 9, 10, 49, 50)))

ggscatter(df_filtered, y = "OverlapProp", x = "nstrains", color = "ecoCom", 
          palette = c("black", "green"),
          ylab = "Cluster Overlap (%)", xlab = "Strains (n)") +
  coord_cartesian(ylim = c(0, 0.8)) 

# overlap vs inhibition
foo %>% 
  filter(ecoCom == "T") %>% 
  ggscatter( x = "OverlapProp", y = "d2", 
            size = 3,
             xlab = "Cluster Overlap (%)", ylab = "Pathogen growth day 2 [CFU/ml]")  + 
  coord_cartesian(ylim = c(10^5, 10^9),xlim=c(0.64,0.72)) +
  scale_y_continuous(trans = "log10", breaks = scales::trans_breaks("log10", function(x) 10^x), labels = scales::trans_format("log10", scales::math_format(10^.x))) 

# strains vs diversity 
foo  %>% filter(CommunityID != "SL1344 itself") %>% filter(CommunityID != "Klebsiella itself") %>%
  ggscatter( x = "nstrains", y = "nCluster", color = "ecoCom", 
             palette = c("black","green"),
             xlab = "n strains", ylab = "cluster diversity") 
