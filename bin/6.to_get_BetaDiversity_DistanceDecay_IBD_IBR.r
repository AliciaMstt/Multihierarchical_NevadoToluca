#**Beta diversity versus Geografic distantes (IBD) from haplotype levels to 7.5% Clustering levels**

#####This script gets all scripts of *multilevel_distance_decay_Sitegeomatrix_IBD* using two arthropods order at multi-hierarchical levels. We used 2 groups: Collembola and Diptera. 
library(stats)
library(base)
library(dplyr)
library(dplyr)
library(tidyr)
library(knitr)
library(PMCMR)
library(vegan)
library(betapart) 
library(stringr)
library(permute)
library(lattice)
library(ecodist)
library(Imap)

#**Showed scripts of the analyzes in Figure5_DistanceDecay_Flat_Altitude.png**
##**We get the *Beta diversity* vs *Flat* and *Altitude* distances analysis from haplotype levels to 7.5% CL at large scale**
##Source to scripts Collembola and Diptera genetic data, and Conductance Flat and altitude matrix.
source("Collembola/multilevel_distance_decay_Collembola_IBResistFlat_Site.R") 
source("Collembola/multilevel_distance_decay_Collembola_IBResist3000_Site.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistFlat_Site.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAltB_Site.R")
#

#**Showed scripts of the analyzes in Figure6_DistanceDecay_FinerScaleFlat.png**
##**We get the *Beta diversity* vs *Flat* distances analysis from haplotype levels to 7.5% CL at finer scale flat**
##Source to scripts Collembola and Diptera genetic data, and Conductance Flat matrix at finer scale.
source("Collembola/multilevel_distance_decay_Collembola_IBResistFlat_East.R") 
source("Collembola/multilevel_distance_decay_Collembola_IBResistFlat_West.R") 
source("Diptera/multilevel_distance_decay_Diptera_IBResistFlat_East.R") 
source("Diptera/multilevel_distance_decay_Diptera_IBResistFlat_West.R") 
#

#**Showed scripts of the analyzes in FigureS6_DistanceDecay_FinerScaleAltitudes.png**
##**We get the *Beta diversity* vs *Altitude* distances analysis from haplotype levels to 7.5% CL at finer scale flat**
##Source to scripts Collembola and Diptera genetic data, and Conductance Altitude matrix at finer scale.
source("Collembola/multilevel_distance_decay_Collembola_IBResist3000_East.R")
source("Collembola/multilevel_distance_decay_Collembola_IBResist3000_West.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAltB_East.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAltB_West.R")
#                

#**Showed scripts of Results in Table S3, S4 and S5 of Distance decay in Collembola and Diptera. Also, values of fractal pattern by a log–log Pearson correlation in Table S6** 
##**We get the *Beta diversity* vs *Flat*, *Altitude*, *slope*, and *vegetation types* distances analysis from haplotype levels to 7.5% CL at large scale *no significant* differences**

##Source to scripts Collembola genetic data, and geographic matrix.
source("Collembola/multilevel_distance_decay_Collembola_IBDist_Site.R") 
source("Collembola/multilevel_distance_decay_Collembola_IBDist_East.R") 
source("Collembola/multilevel_distance_decay_Collembola_IBDist_West.R") 

#**We get the Beta diversity vs Geographic distances analysis from haplotype levels to 7.5% CL**
##Source to scripts Diptera, genetic data, and geographic matrix.
source("Diptera/multilevel_distance_decay_Diptera_IBDist_Site.R") 
source("Diptera/multilevel_distance_decay_Diptera_IBDist_East.R") 
source("Diptera/multilevel_distance_decay_Diptera_IBDist_West.R") 

#**END**
