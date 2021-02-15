#**Defining lineages at different clustering levels for Diptera order (date set "b"): STEP 2** 

#**This script gets Specie delimitation using **Arachnida order**. For this, the ASV filtered dataset was used to generate an UPGMA tree with corrected distances under a F84 model, and based on this tree, all haplotypes were nested into clusterin levels (CLs) following the genetic similarity at different thresholds (0.5%, 1.5%, 3%, 5% and 7.5%), plus an additional threshold corresponding to the result of a specie delimitation analyses conducted with the generalized mixed Yule-coalescent (GMYC).

#**PART 1. Distance Matrix and UPGMA**
library(ape)
library(phangorn)
library(seqinr)
library(stringr)

Diptera1016 <- read.FASTA("../../../../genetic/Data_in/Diptera/7Diptera1016_zotus_Nostops.fasta")  #leer fasta
dist.dna(Diptera1016, model = "F84", as.matrix = T)->dist_Diptera1016  #calcular matriz distancias por pares
Diptera1016_UPGMA <- upgma(dist_Diptera1016)  #hacer Arbol
plot(Diptera1016_UPGMA)

#PART 3: To obtain subtrees at multihierarchical levels.
#This funcion are in bin. NOTA Nos muestra un warning, aunque sigue estando bien.

source("../../../NODE.MIN_Get_trees_multiple.R")

#**Diptera b**
#We used the result of the threshold time: e.g. "0.01258205" in the next analisys of each lineages.

#threshold time Diptera_a: -0.009679657
#threshold time Diptera_b: -0.00786858
#Sum: (0.009679657 + 0.00786858)/2= 0.0087741185

read.table ("../../../../genetic/Data_in/Diptera/7Diptera1016_zotus_Nostops.fasta", header=FALSE)->Diptera1016.fasta
#read.tree("Diptera1016_All_s2_sinSTOPS_revised.newick")->Diptera1016_SHORT2.tree
GMYC.thershold<-0.0087741185


NODE.MIN(Diptera1016.fasta,Diptera1016_UPGMA,0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.075, 0.075,0.075, GMYC.thershold, print.subtrees = "NO", print.subtrees.fasta = "NO")

#**END**




