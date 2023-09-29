# Code to simulate community signatures in silico & to calculate functional overlap to a given pathogen 

#'Preparation: 
#'Define Genomes of Interest as set in PATRIC
#'Download Genome Metadata from Patric > data/PATRIC_metadata.xlsx (add species identifiers, matching Experimental data sheets)

#'Application
#'This code can be adapted to run BIOLOG- based predictions


# Set variables
## Pathogen
PATHOGEN = "Eco" # Either:  "Salmo" or "Klebs" or "Eco" 
## Community set 
COMMUNITYSET = "all" # "all", or "TOP10


# Set env
setwd("~/Downloads/Nutrient_blocking-main/")
library(plyr)
library(readxl)
library(tidyverse)
library(stringr)


# Set data 
## Metadata
genomeMetadata = "data/PATRIC_metadata.xlsx" # Genome metadata from PATRIC (thank you!)
genome_meta <- read_excel(genomeMetadata) # Load genome info table
genome_meta = subset(genome_meta , !genome_meta$StrainID %in% c("EcoHS","EcoMG1655", "EcoZ1269","EcoZ1331")) # focus on EcoIAI1 as community member


## Load Community 
annotation <- data.frame() #  init
for (i in list.files(pattern = ".features.tab$",path = "data/annotations", recursive = TRUE, full.names = T)) {
  print(i)
  tmp <- read.csv(i, sep= "\t", stringsAsFactors = FALSE )
  tmp$genome_id = as.character(tmp$genome_id)
  annotation= rbind.fill(annotation , tmp) 
  }
  tmp = NULL
  # sanity checks 
  length(unique(annotation$genome_id)) # n genomes in community dataset = 50
  # ifs match 
  down = unique(as.character(as.numeric(annotation$genome_id)))
  should = unique(as.character(as.numeric(genome_meta$`Genome ID`)))
  table(should %in% down) # # match
  
  
# clean data & merge 
annotation_select <- subset(annotation, feature_type == "CDS" ) # focus on CDS
set = unique(annotation_select[,c("genome_id","genome_name")]) # check set assayed

annotation_select = unique(annotation_select[,c("genome_id","pgfam_id")]) # select relevant columns
annotation_select = subset(annotation_select , annotation_select$pgfam_id != "") # focus on genes with PGfam annotation 

# Possible strain space 

strainUniverse = unique(annotation_select$genome_id) # list ofstrains used for sims 

# optional subset for top 10's
if (COMMUNITYSET == "TOP10"){
  strainUniverse = c("555970.3","610130.3","1423823.4","349741.6","272621.13","565042.3","585034.5","435590.9","411460.6","1423799.3")
}


# Load Pathogen 

if (PATHOGEN == "Salmo"){
  Patho = read.csv("data/annotations/Salmonella_216597.6.PATRIC.pathogen.tab.tsv", sep = "\t", header = T)# Salmonella 
  Patho <- subset(Patho, feature_type == "CDS" )
  Patho = unique(Patho$pgfam_id)
  print("Simulating Salmonella")
} else if (PATHOGEN == "Klebs"){
  Patho = read.csv("data/annotations/Klebsiella_1162296.3.PATRIC.pathogen.tab.tsv", sep = "\t", header = T)# Klebsiella
  Patho <- subset(Patho, feature_type == "CDS" )
  Patho = unique(Patho$pgfam_id)
  print("Simulating Klebsiella")
} else if (PATHOGEN == "Eco"){
  Patho = read.csv("data/annotations/Escherichia_coli_19Y000018.PATRIC.pathogen.tab.tsv", sep = "\t", header = T) # Escherichia_coli_19Y000018
  Patho <- subset(Patho, Feature.Type == "CDS" )
  Patho = unique(Patho$PATRIC.cross.genus.families..PGfams.)
  print("Simulating E coli")
} else {
  print("PATHOGEN not in the list")
}


###################################
# Simulate E.coli containing communities 
###################################

strainUniverse2 = strainUniverse[strainUniverse != "585034.5"] # remove Eco 585034.5 from strain universe and add it later everywhere
foo <- data.frame()
com.annot = NULL
PossibleComs = NULL
counter = 1
samples = 100000000000000 # max picking trials
nCOMs = 1000 #  unique communities to be simulated   


for (i in 1:samples){
  # pick a random list of strains
  strain.list = sample(x = strainUniverse2 , size = sample(c(1,2,4,9),1, replace = F), replace = F) # pick random community of n= 1,2,4,9 + Eco
  strain.list = c(strain.list, "585034.5") # add Eco 
  nstrains = length(strain.list)
  com.annot = subset(annotation_select , annotation_select$genome_id %in% strain.list) # focal strain annotation
  nPGFAM = length(unique(com.annot$pgfam_id))
  CommunityComp = paste(sort(strain.list), collapse = ", ")
  
  # Calculate overlap pathogen - community PGfams overlap
  ComPGFAM = unique(com.annot$pgfam_id) # list community PGfam repertoire
  
  PathoOverlap = Patho %in% ComPGFAM # TRUE = Overlap, FALSE = no overlap
  OverlapCluster = table(Patho %in% ComPGFAM)["TRUE"] # covered
  DifferCluster = table(Patho %in% ComPGFAM)["FALSE"] # not covered
  
  OverlapProp = OverlapCluster / (OverlapCluster + DifferCluster)
  
  # Feed summary table
  foo[i,"CommunityMembers"] = CommunityComp
  foo[i,"nstrains"] = nstrains
  foo[i,"nPGFAM"] = nPGFAM
  foo[i,"OverlapProp"] = OverlapProp
  foo[i,"Eco_present"] = any(c("331112.6","511145.12","562.68767","562.6877","585034.5") %in% strain.list)
  #print(i/samples*100)
  uniqueSamples = nrow(unique(foo))
  print(uniqueSamples)
  
  # if condition with break
  if(uniqueSamples >= nCOMs ) {
    break
  }
    }
EcoComs = unique(foo)


###################################
# Simulate n random communities 
###################################

strainUniverse2 = strainUniverse
foo <- data.frame()
com.annot = NULL
PossibleComs = NULL
counter = 1
samples = 100000000000000 # max picking trials
nCOMs = 100  #  unique communities to be simulated   

for (i in 1:samples){
  # pick a random list of strains
  strain.list = sample(x = strainUniverse2 , size = sample(c(2,3,5,10),1, replace = F), replace = F) # pick random community of n= 2,3,5,10 
  nstrains = length(strain.list)
  com.annot = subset(annotation_select , annotation_select$genome_id %in% strain.list) # focal strain annotation
  nPGFAM = length(unique(com.annot$pgfam_id))
  CommunityID = paste0("Com_",i)
  CommunityComp = paste(sort(strain.list), collapse = ", ")
  # Calculate overlap pathogen - community PGfams overlap
  ComPGFAM = unique(com.annot$pgfam_id) # list community Pgfam repertoire
  
  PathoOverlap = Patho %in% ComPGFAM # TRUE = Overlap, FALSE = no overlap
  OverlapCluster = table(Patho %in% ComPGFAM)["TRUE"] # covered
  DifferCluster = table(Patho %in% ComPGFAM)["FALSE"] # not covered
  OverlapProp = OverlapCluster / (OverlapCluster + DifferCluster)
  
  # Feed summary table
  foo[i,"CommunityMembers"] = CommunityComp
  foo[i,"nstrains"] = nstrains
  foo[i,"nPGFAM"] = nPGFAM
  foo[i,"OverlapProp"] = OverlapProp
  foo[i,"Eco_present"] = any(c("331112.6","511145.12","562.68767","562.6877","585034.5") %in% strain.list)
  #print(i/samples*100)
  uniqueSamples = nrow(unique(foo))
  print(uniqueSamples)
  
  # if condition with break
  if(uniqueSamples >= nCOMs) {
    break
  }
}
AllComs = unique(foo)
