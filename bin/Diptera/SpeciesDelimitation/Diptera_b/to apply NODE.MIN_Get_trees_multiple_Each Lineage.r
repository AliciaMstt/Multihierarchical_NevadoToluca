setwd("~/Documents/NancyTesis/GMYCyLimites/7Diptera2031readsRealigned1b")

#PARTE 1. Matriz de distancias y UPGMA

#abrir paquetes necesarios ya instalados 
library(ape)
library(phangorn)
library(seqinr) 

Diptera1016 <- read.FASTA("7Diptera1016.fasta")  #leer fasta
dist.dna(Diptera1016, model = "F84", as.matrix = T)->dist_Diptera1016  #calcular matriz distancias por pares
Diptera1016_UPGMA <- upgma(dist_Diptera1016)  #hacer Arbol
plot(Diptera1016_UPGMA)

#Parte 3: Correr funci?n para obtener subtrees a multiples niveles


library("stringr")
#Funcion que he puesto en el directorio de trabajo
source("NODE.MIN_Get_trees_multiple.R")


##################### Diptera1016 ################################

read.table ("7Diptera1016.fasta", header=FALSE)->Diptera1016.fasta
#read.tree("Diptera1016_All_s2_sinSTOPS_revised.newick")->Diptera1016_SHORT2.tree
GMYC.thershold<-0.0087741185


NODE.MIN(Diptera1016.fasta,Diptera1016_UPGMA,0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.04, 0.05, 0.06, 0.075, 0.075,0.075, GMYC.thershold, print.subtrees = "NO", print.subtrees.fasta = "NO")





