#**Beta diversity versus Geografic distantes (IBD) from haplotype levels to lineages 7.5%**

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
##**We get the *Beta diversity* vs *Flat* and *Altitude* distances analysis from haplotype levels to lineages 7.5% at large scale**
##Source to scripts Collembola and Diptera genetic data, and Conductance Flat and altitude matrix.
source("Collembola/multilevel_distance_decay_Collembola_IBResistFlat_Site.R") 
source("Collembola/multilevel_distance_decay_Collembola_IBResist3000_Site.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistFlat_Site.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAltB_Site.R")
#

#**Showed scripts of the analyzes in Figure6_DistanceDecay_FinerScaleFlat.png**
##**We get the *Beta diversity* vs *Flat* distances analysis from haplotype levels to lineages 7.5% at finer scale flat**
##Source to scripts Collembola and Diptera genetic data, and Conductance Flat matrix at finer scale.
source("Collembola/multilevel_distance_decay_Collembola_IBResistFlat_East.R") 
source("Collembola/multilevel_distance_decay_Collembola_IBResistFlat_West.R") 
source("Diptera/multilevel_distance_decay_Diptera_IBResistFlat_East.R") 
source("Diptera/multilevel_distance_decay_Diptera_IBResistFlat_West.R") 
#

#**Showed scripts of the analyzes in FigureS6_DistanceDecay_FinerScaleAltitudes.png**
##**We get the *Beta diversity* vs *Altitude* distances analysis from haplotype levels to lineages 7.5% at finer scale flat**
##Source to scripts Collembola and Diptera genetic data, and Conductance Altitude matrix at finer scale.
source("Collembola/multilevel_distance_decay_Collembola_IBResist3000_East.R")
source("Collembola/multilevel_distance_decay_Collembola_IBResist3000_West.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAltB_East.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAltB_West.R")
#                

#**Showed scripts of the analyzes in FigureS7_DistanceDecay_Sites.png**
##**We get the *Beta diversity* vs *Flat* distances analysis from haplotype levels to lineages 7.5% at large scale**
##Source to scripts Arachnida, Coleoptera, Hemiptera and Hymenoptera genetic data.

source("Arachnida/multilevel_distance_decay_Arachnida_IBResistFlat_Site.R")
source("Coleoptera/multilevel_distance_decay_Coleoptera_IBResistFlat_Site.R")
source("Hemiptera/multilevel_distance_decay_Hemiptera_IBResistFlat_Site.R")
source("Hymenoptera/multilevel_distance_decay_Hymenoptera_IBResistFlat_Site.R")
#

#**Showed scripts of Results in Table S3, S4 and S5 of Distance decay in Collembola and Diptera. Also, values of fractal pattern by a log–log Pearson correlation in Table S6** 
##**We get the *Beta diversity* vs *Flat*, *Altitude*, *slope*, and *vegetation types* distances analysis from haplotype levels to lineages 7.5% at large scale**

#**Vegetation type**
source("Collembola/multilevel_IBR_CollembolaResistance_VegRcl_A_g.R") 
source("Collembola/multilevel_IBR_CollembolaResistance_VegRcl_B_d.R")
source("Collembola/multilevel_IBR_CollembolaResistance_VegRcl_C_c.R") 
source("Collembola/multilevel_IBR_CollembolaResistance_VegRcl_D_a.R")
source("Collembola/multilevel_IBR_CollembolaResistance_VegRcl_E_b.R")
source("Collembola/multilevel_IBR_CollembolaResistance_VegRcl_F_e.R") 
source("Collembola/multilevel_IBR_CollembolaResistance_VegRcl_G_f.R")
source("Collembola/multilevel_IBR_CollembolaResistance_VegRcl_H_h.R") 
source("Collembola/multilevel_IBR_CollembolaResistance_VegRcl_I_i.R")  
source("Collembola/multilevel_IBR_CollembolaResistance_VegRcl_J.R") 

source("Diptera/multilevel_IBR_DipteraResistance_VegRcl_A_g.R")
source("Diptera/multilevel_IBR_DipteraResistance_VegRcl_B_d.R")
source("Diptera/multilevel_IBR_DipteraResistance_VegRcl_C_c.R") 
source("Diptera/multilevel_IBR_DipteraResistance_VegRcl_D_a.R")
source("Diptera/multilevel_IBR_DipteraResistance_VegRcl_E_b.R")
source("Diptera/multilevel_IBR_DipteraResistance_VegRcl_F_e.R") 
source("Diptera/multilevel_IBR_DipteraResistance_VegRcl_G_f.R")
source("Diptera/multilevel_IBR_DipteraResistance_VegRcl_H_h.R") 
source("Diptera/multilevel_IBR_DipteraResistance_VegRcl_I_i.R")  
source("Diptera/multilevel_IBR_DipteraResistance_VegRcl_J.R") 

#**Elevation**
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt2000.R") 
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt2700.R")
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt2800.R")
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt2900.R")
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt3100.R")
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt3200.R")
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt3300.R")                  
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt3400.R")
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt3500.R")
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt3600.R")
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt3700.R")
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt_A.R")
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt_B.R")
source("Collembola/multilevel_distance_decay_Collembola_IBResistAlt_C.R")

source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt2000.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt2700.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt2800.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt2900.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt3000.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt3100.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt3200.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt3300.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt3400.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt3500.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt3600.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt3700.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt_A.R")
source("Diptera/multilevel_distance_decay_Diptera_IBResistAlt_C.R")

#**Slope**
source("Collembola/multilevel_IBR_CollembolaResistance_Slope_A.R")
source("Collembola/multilevel_IBR_CollembolaResistance_Slope_B.R")
source("Collembola/multilevel_IBR_CollembolaResistance_Slope_C.R")
source("Collembola/multilevel_IBR_CollembolaResistance_Slope_D.R")

source("Diptera/multilevel_IBR_DipteraResistance_Slope_A.R")
source("Diptera/multilevel_IBR_DipteraResistance_Slope_B.R")
source("Diptera/multilevel_IBR_DipteraResistance_Slope_C.R")
source("Diptera/multilevel_IBR_DipteraResistance_Slope_D.R")

#**We get the Beta diversity vs Geographic distances analysis from haplotype levels to lineages 7.5%**
##Source to scripts Collembola and Diptera genetic data, and geographic matrix.
source("Collembola/multilevel_distance_decay_Collembola_IBDist_Site.R") 
source("Collembola/multilevel_distance_decay_Collembola_IBDist_East.R") 
source("Collembola/multilevel_distance_decay_Collembola_IBDist_West.R") 

source("Diptera/multilevel_distance_decay_Diptera_IBDist_Site.R") 
source("Diptera/multilevel_distance_decay_Diptera_IBDist_East.R") 
source("Diptera/multilevel_distance_decay_Diptera_IBDist_West.R") 

#**END**
