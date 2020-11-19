
#'**Community diversity and composition at the haplotype, 3% and 5% Clustering levels for each of the eight taxonomic orders studied**'#

#This script gets all scripts of diversity using 8 arthropods order at multi-hierarchical levels 
#We used 8 groups: Arachnda, Coleoptera, Collembola, Diptera, Hemiptera, Hymenoptera, Lepidoptera, and Myriapoda. 

library(stats)
library(base)
library(dplyr)
library(tidyr)
library(knitr)
library(PMCMR)
library(rcdd) 
library(vegan)
library(betapart) 
library(stringr)


#'**We get the diversity analysis at haplotypes, CL 3%, and CL 5% in Diptera**'#
# Source to scripts Diptera and data
source("Diptera/Com_matrixes_data_exploration_Diptera_h_Site.R") 
source("Diptera/Com_matrixes_data_exploration_Diptera_3P_Site.R") 
source("Diptera/Com_matrixes_data_exploration_Diptera_5P_Site.R") 

#'**We get the diversity analysis at haplotypes, CL 3%, and CL 5% in Collembola**'#
# Source to scripts Collembola and data
source("Collembola/Com_matrixes_data_exploration_Collembola_h_Site.R") 
source("Collembola/Com_matrixes_data_exploration_Collembola_3P_Site.R") 
source("Collembola/Com_matrixes_data_exploration_Collembola_5P_Site.R") 

#'**We get the diversity analysis at haplotypes, CL 3%, and CL 5% in Arachnida**'#
# Source to scripts Arachnida and data
source("Arachnida/Com_matrixes_data_exploration_Arachnida_h_Site.R") 
source("Arachnida/Com_matrixes_data_exploration_Arachnida_3P_Site.R") 
source("Arachnida/Com_matrixes_data_exploration_Arachnida_5P_Site.R")

#'**We get the diversity analysis at haplotypes, CL 3%, and CL 5% in Hemiptera**'#
# Source to scripts Hemiptera and data
source("Hemiptera/Com_matrixes_data_exploration_Hemiptera_h_Site.R") 
source("Hemiptera/Com_matrixes_data_exploration_Hemiptera_3P_Site.R") 
source("Hemiptera/Com_matrixes_data_exploration_Hemiptera_5P_Site.R") 

#'**We get the diversity analysis at haplotypes, CL 3%, and CL 5% in Hymenoptera**'#
# Source to scripts Hymenoptera and data
source("Hymenoptera/Com_matrixes_data_exploration_Hymenoptera_h_Site.R") 
source("Hymenoptera/Com_matrixes_data_exploration_Hymenoptera_3P_Site.R") 
source("Hymenoptera/Com_matrixes_data_exploration_Hymenoptera_5P_Site.R") 

#'**We get the diversity analysis at haplotypes, CL 3%, and CL 5% in Coleoptera**'#
# Source to scripts Coleoptera and data
source("Coleoptera/Com_matrixes_data_exploration_Coleoptera_h_Site.R") 
source("Coleoptera/Com_matrixes_data_exploration_Coleoptera_3P_Site.R") 
source("Coleoptera/Com_matrixes_data_exploration_Coleoptera_5P_Site.R") 

#'**We get the diversity analysis at haplotypes, CL 3%, and CL 5% in Myriapoda**'#
# Source to scripts Myriapoda and data
source("Myriapoda/Com_matrixes_data_exploration_Myriapoda_h_Site.R") 
source("Myriapoda/Com_matrixes_data_exploration_Myriapoda_3P_Site.R") 
source("Myriapoda/Com_matrixes_data_exploration_Myriapoda_5P_Site.R")

#'**We get the diversity analysis at haplotypes, CL 3%, and CL 5% in Lepidoptera#'**
# Source to scripts Lepidoptera and data
source("Lepidoptera/Com_matrixes_data_exploration_Lepidoptera_h_Site.R") 
source("Lepidoptera/Com_matrixes_data_exploration_Lepidoptera_3P_Site.R") 
source("Lepidoptera/Com_matrixes_data_exploration_Lepidoptera_5P_Site.R") 




####################################END######################################################################



