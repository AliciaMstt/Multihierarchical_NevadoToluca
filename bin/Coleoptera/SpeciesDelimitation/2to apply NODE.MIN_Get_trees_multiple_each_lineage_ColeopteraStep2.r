#**Defining lineages at different clustering levels for Coleoptera order: STEP 2** 

#**This script gets Specie delimitation using **Coleoptera order**. For this, the ASV filtered dataset was used to generate an UPGMA tree with corrected distances under a F84 model, and based on this tree, all haplotypes were nested into clusterin levels (CLs) following the genetic similarity at different thresholds (0.5%, 1.5%, 3%, 5% and 7.5%), plus an additional threshold corresponding to the result of a specie delimitation analyses conducted with the generalized mixed Yule-coalescent (GMYC).

#**PART 1. Distance Matrix and UPGMA**
library(ape)
library(phangorn)
library(seqinr)
library(stringr)

Coleoptera <- read.FASTA("../../../genetic/Data_in/Coleoptera/6Coleoptera209_zotus_Nostops.fasta")  #leer fasta
dist.dna(Coleoptera, model = "F84", as.matrix = T)->dist_Coleoptera  #calcular matriz distancias por pares
Coleoptera_UPGMA <- upgma(dist_Coleoptera)  #hacer árbol
plot(Coleoptera_UPGMA)

#PART 3: To obtain subtrees at multihierarchical levels.
#This funcion are in bin. NOTA Nos muestra un warning, aunque sigue estando bien.

source("../../NODE.MIN_Get_trees_multiple.R")

#**Coleoptera**

read.table ("../../../genetic/Data_in/Coleoptera/6Coleoptera209_zotus_Nostops.fasta", header=FALSE)->Coleoptera.fasta
#read.tree("Coleoptera_All_s2_sinSTOPS_revised.newick")->Coleoptera_SHORT2.tree
GMYC.thershold<-0.007255866   

NODE.MIN(Coleoptera.fasta,Coleoptera_UPGMA,0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.075, 0.075, 0.075, GMYC.thershold, print.subtrees = "NO", print.subtrees.fasta = "NO")

#**END**

