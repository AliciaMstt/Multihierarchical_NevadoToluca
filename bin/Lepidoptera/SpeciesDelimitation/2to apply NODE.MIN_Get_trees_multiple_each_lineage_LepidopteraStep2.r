#**Defining lineages at different clustering levels for Lepidoptera order: STEP 2** 

#**This script gets Specie delimitation using **Lepidoptera order**. For this, the ASV filtered dataset was used to generate an UPGMA tree with corrected distances under a F84 model, and based on this tree, all haplotypes were nested into clusterin levels (CLs) following the genetic similarity at different thresholds (0.5%, 1.5%, 3%, 5% and 7.5%), plus an additional threshold corresponding to the result of a specie delimitation analyses conducted with the generalized mixed Yule-coalescent (GMYC).

#**PART 1. Distance Matrix and UPGMA**
library(ape)
library(phangorn)
library(seqinr)
library(stringr)

Lepidoptera <- read.FASTA("../../../genetic/Data_in/Lepidoptera/5Lepidoptera60_zotus_Nostops.fasta")  #leer fasta
dist.dna(Lepidoptera, model = "F84", as.matrix = T)->dist_Lepidoptera  #calcular matriz distancias por pares
Lepidoptera_UPGMA <- upgma(dist_Lepidoptera)  #hacer arbol
plot(Lepidoptera_UPGMA)

#PART 3: To obtain subtrees at multihierarchical levels.
#This funcion are in bin. NOTA Nos muestra un warning, aunque sigue estando bien.

source("../../NODE.MIN_Get_trees_multiple.R")

#**Arachnida**
#We used the result of the threshold time: e.g. "0.002997824 in the next analisys of each lineages.
read.table ("../../../genetic/Data_in/Lepidoptera/5Lepidoptera60_zotus_Nostops.fasta", header=FALSE)->Lepidoptera.fasta
#read.tree("Lepidoptera_All_s2_sinSTOPS_revised.newick")->Lepidoptera_SHORT2.tree
GMYC.thershold<-0.002997824

NODE.MIN(Lepidoptera.fasta,Lepidoptera_UPGMA,0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.075, 0.075, 0.075, GMYC.thershold, print.subtrees = "NO", print.subtrees.fasta = "NO")

#**END**




