#**Defining Diptera lineage at different clustering levels: STEP 1** 

#**This script gets Specie delimitation using **Diptera order with two data set: "a" and "b"**. For this, the ASV filtered dataset was used to generate an UPGMA tree with corrected distances under a F84 model, and based on this tree, all haplotypes were nested into clusterin levels (CLs) following the genetic similarity at different thresholds (0.5%, 1.5%, 3%, 5% and 7.5%), plus an additional threshold corresponding to the result of a specie delimitation analyses conducted with the generalized mixed Yule-coalescent (GMYC).

#PART 1. Distance matrix and UPGMA
library(igraph)
library(ape)
library(phangorn)
library(permute)
library(geiger)
library(lattice)
library(stats)
library(vegan)
library(base)
library(MASS)
library(paran)
library(gtools)
library(splits)
library(seqinr)
library(spam)
library(dotCall64)
library(grid)
library(maps)
library(fields)

#**Diptera a**
Diptera1015 <- read.FASTA("../genetic/Data_in/Diptera/7Diptera1015_zotus_Nostops.fasta")  #leer fasta
dist.dna(Diptera1015, model = "F84", as.matrix = T)->dist_Diptera1015  #calcular matriz distancias por pares
Diptera1015_UPGMA <- upgma(dist_Diptera1015)  #hacer árbol
plot(Diptera1015_UPGMA)

#PART 2: GMYC

#Those funsions are in bin
source("gmyc.pkg.0.9.6.R")
source("Powell_supplemental_script.R")

#**single gmyc threshold**
Diptera1015_UPGMA.GMYC.singleThreshold<-gmyc.edit(Diptera1015_UPGMA, method="s") #hacer GMYC
#resumen de los resultados
summary(Diptera1015_UPGMA.GMYC.singleThreshold) #Mostrar resumen GMYC

sink("Diptera/SpeciesDelimitation/Diptera_b/Diptera1015_UPGMA.GMYC.singleThreshold.log", type=c("output","message"), split=TRUE) #guardar resultado del GMYC (resumen)
#MIT.test.s o MIT.test.m para ver todos los datos
summary(Diptera1015_UPGMA.GMYC.singleThreshold)
sink()


#**Diptera b**
Diptera1016 <- read.FASTA("../genetic/Data_in/Diptera/7Diptera1016_zotus_Nostops.fasta")  #leer fasta
dist.dna(Diptera1016, model = "F84", as.matrix = T)->dist_Diptera1016  #calcular matriz distancias por pares
Diptera1016_UPGMA <- upgma(dist_Diptera1016)  #hacer árbol
plot(Diptera1016_UPGMA)

#PART 2: GMYC

#Those funsions are in bin
source("gmyc.pkg.0.9.6.R")
source("Powell_supplemental_script.R")

#**single gmyc threshold**
Diptera1016_UPGMA.GMYC.singleThreshold<-gmyc.edit(Diptera1016_UPGMA, method="s") #hacer GMYC
#resumen de los resultados
summary(Diptera1016_UPGMA.GMYC.singleThreshold) #Mostrar resumen GMYC

sink("Diptera/SpeciesDelimitation/Diptera_a/Diptera1016_UPGMA.GMYC.singleThreshold.log", type=c("output","message"), split=TRUE) #guardar resultado del GMYC (resumen)
#MIT.test.s o MIT.test.m para ver todos los datos
summary(Diptera1016_UPGMA.GMYC.singleThreshold)
sink()

#**END**
