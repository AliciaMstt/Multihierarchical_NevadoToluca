#'**INITIAL STEPS**

#'In excel remove the simbol #' from the names of the table and rename the samples acccording to the code used in the gradient e.g. GRA_S10_D_F_A10



library(stats)
library(base)
library(dplyr)
library(dplyr)
library(tidyr)
library ("knitr")
library(PMCMR)
library("vegan")
library("betapart") 
library("stringr")
library(permute) #Me pide que carge permute y lattice en el Rstudio que estoy usando
library(lattice)
library(ggsignif)
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(ggsignif)
library(broom)

#'**TABLES AND COMMUNITY MATRIXES**'# 

################################################################################################################################################################################

###########'open table with names including Region and habitat parameters
s2_raw_all <- read.table("../genetic/Data_in/Arachnida/s2_raw_all_Arachnida_threshold.txt", header=TRUE)
dim(s2_raw_all)

###########'remove additional columns and leave only names (of haplotipes), samples and taxa (and threshold in this case)
s2_raw_all[,c(1:66,68)]->s2_raw
dim(s2_raw) #'49 samples = 48 plus 1 neg (the second neg from DOM_REPS is not there because all 0)
colnames(s2_raw)

###########'Applying the conservative threshold (this is a binary column)
s2_raw[which(s2_raw$conservative_threshold == "1"),]->s2_raw_threshold 
s2_raw_threshold [,1:66]->s2_raw_threshold #'remove threshold col
dim(s2_raw_threshold)
colnames(s2_raw_threshold)


################################################################################################################################################################################

#'loop to create the new matrix combining haplotype by otu pertenency, i.e. submatrix by limit
#'
#'*Arachnida*
unique (s2_raw_threshold$limite0.03)->levels_limite0.03

data.frame()->s2_raw_Arachnida_limite0.03

for (i in 1:length(unique (s2_raw_threshold$limite0.03)))
{
  levels_limite0.03[i]->level0.03
  s2_raw_threshold[which(s2_raw_threshold$limite0.03==level0.03),]->subcom_level_names0.03
  subcom_level_names0.03[,c(2:52)]->subcom_level0.03  #delete names, level and also the negative column
  colSums(subcom_level0.03)->sum0.03
  as.data.frame(sum0.03)->sum0.03
  t(sum0.03)->sum0.03
  row.names(sum0.03)<-subcom_level_names0.03[1,1] #keep the name of the first haplotype
  rbind(s2_raw_Arachnida_limite0.03,sum0.03)->s2_raw_Arachnida_limite0.03
}

#'delete the ocurrences with less than 4 reads by library (same criteria than denoising)
#s2_raw_Arachnida_limite0.03->s2_f4_with_abundance_Arachnida_limite0.03 
#s2_f4_with_abundance_Arachnida_limite0.03[s2_f4_with_abundance_Arachnida_limite0.03<4]<-0  #'2 warning corresponding wiht the columms of the names and taxa

#'transform in present/absence table
s2_raw_Arachnida_limite0.03->s2_raw_Arachnida_limite0.03
s2_raw_Arachnida_limite0.03[s2_raw_Arachnida_limite0.03>1]<-1 #'transform in present/absence table 

#'transform in present/absence table
#s2_f4_with_abundance->s2_f4 #NOTA_Nancy: Tengo un subset de Coleoptear
#s2_f4[s2_f4>1]<-1 #'2 warning corresponding wiht the columms of the names and taxa

#'checking if there is any row with no presence
s2_raw_Arachnida_limite0.03[,1:51]->data0.03 #Nota_Nancy: Modifique el numero: 2:50 por 1:51, aunque no funcionó. Volví a la version de 2:50
rowSums(data0.03)
length(which(rowSums(data0.03)!=0))
length(which(rowSums(data0.03)==0))


#'Community matrixes (samples in rows and h in cols). 
#'*Arachnida*
t(s2_raw_Arachnida_limite0.03)->t_s2_f4_Arachnida_limite0.03 #'trasp
t_s2_f4_Arachnida_limite0.03[1:51,]->community_Arachnida_limite0.03 #NOTA_Nancy: Este numero es importante. Colocar exactamente el numero de "s2_f4[,2:52]->data".
colnames(t_s2_f4_Arachnida_limite0.03)<-community_Arachnida_limite0.03[1,] 
as.data.frame(community_Arachnida_limite0.03)->community_Arachnida0.03 #'trasp including col and row names
#community_Acari[-49,]->community_Arachnida #'removing neg
dim(community_Arachnida0.03)
community_Arachnida0.03[order(row.names(community_Arachnida0.03)),]->community_Arachnida0.03 #'order samples
write.table (community_Arachnida0.03, file="../genetic/Data_out/Arachnida/Arachnida3P/community_Arachnida0.03.txt") #'this is necessary for the format, not able to solve in other way
read.table ("../genetic/Data_out/Arachnida/Arachnida3P/community_Arachnida0.03.txt")->community_Arachnida0.03


#'submatrixes by MOUNTAIN. NOTA_Nancy: Quiero hacer tablas que incluyan datos de las tres montañas en la estación de seca.
dim(community_Arachnida0.03)
community_Arachnida0.03[which(str_extract (row.names(community_Arachnida0.03), "_M_") %in% "_M_"),]->community_Arachnida0.03_Mountain
dim(community_Arachnida0.03_Mountain)
community_Arachnida0.03_Mountain[,which(colSums(community_Arachnida0.03_Mountain)!=0)]->community_Arachnida0.03_Mountain #'to remove no data colums
dim(community_Arachnida0.03_Mountain)


#'submatrixes by SITE in Nevado Toluca. NOTA_Nancy: Quiero hacer tablas que incluyan datos de una montañas en la epoca de lluvida con localidades dentro del Nevado de Toluca.
dim(community_Arachnida0.03)
community_Arachnida0.03[which(str_extract (row.names(community_Arachnida0.03), "_NTO_") %in% "_NTO_"),]->community_Arachnida0.03_Site
dim(community_Arachnida0.03_Site)
community_Arachnida0.03_Site[,which(colSums(community_Arachnida0.03_Site)!=0)]->community_Arachnida0.03_Site #'to remove no data colums
dim(community_Arachnida0.03_Site)


####################################################

#'Generating a general table with names and habitat parameters.
#BY MOUNTAIN
row.names(community_Arachnida0.03_Mountain)->sample_names_Mountain0.03
as.data.frame(sample_names_Mountain0.03)->sample_names_Mountain0.03
sample_names_Mountain0.03 %>% separate(sample_names_Mountain0.03, c("Conservation","Region","Mountain","Site","ID"), sep="_",remove=FALSE)->general_sample_RegionMountain0.03
general_sample_RegionMountain0.03
general_sample_RegionMountain0.03 %>% unite(RegionMountain, Region, Mountain, sep="_",remove=FALSE)->general_sample_RegionMountain0.03 #'generating a variable combining layer and habitat
general_sample_RegionMountain0.03
write.table(general_sample_RegionMountain0.03, file="../genetic/Data_out/Arachnida/Arachnida3P/general_sample_RegionMountain0.03.txt") #'this is the only way I found to be able to work later
read.table("../genetic/Data_out/Arachnida/Arachnida3P/general_sample_RegionMountain0.03.txt",header=TRUE)->general_sample_RegionMountain0.03


#BY SITE
#'Generating a general table with names and habitat parameters NOTA_Nancy: POR EL MOMENTO NO para este analisis porque no tengo habitats. Mas adelante cuando decida ver todo por montañas o localida
row.names(community_Arachnida0.03_Site)->sample_names_Mountain1_0.03
as.data.frame(sample_names_Mountain1_0.03)->sample_names_Mountain1_0.03
sample_names_Mountain1_0.03 %>% separate(sample_names_Mountain1_0.03, c("Conservation","Mountain1","Site","ID"), sep="_",remove=FALSE)->general_sample_Mountain1Site0.03
general_sample_Mountain1Site0.03
general_sample_Mountain1Site0.03 %>% unite(Mountain1andSite, Mountain1,Site, sep="_",remove=FALSE)->general_sample_Mountain1Site0.03 #'generating a variable combining layer and habitat
general_sample_Mountain1Site0.03
write.table(general_sample_Mountain1Site0.03, file="../genetic/Data_out/Arachnida/Arachnida3P/general_sample_Mountain1Site0.03.txt") #'this is the only way I found to be able to work later
read.table("../genetic/Data_out/Arachnida/Arachnida3P/general_sample_Mountain1Site0.03.txt",header=TRUE)->general_sample_Mountain1Site


####################################################
####################################################

#'**HAPLOTYPE RICHNESS TABLES, PLOTS AND ANALYSES**'# BY MOUNTAIN

#'*Arachnida*

as.matrix(community_Arachnida0.03_Mountain)->community_Arachnida0.03_Mountain
row.names(community_Arachnida0.03_Mountain)->sample_names_Mountain0.03
dim(community_Arachnida0.03_Mountain)->dims_Mountain0.03
dims_Mountain0.03
.rowSums (community_Arachnida0.03_Mountain,dims_Mountain0.03[1],dims_Mountain0.03[2])->sample_richness_Mountain0.03 #'summatory by rows
rbind(sample_names_Mountain0.03,sample_richness_Mountain0.03)->richness_Mountain0.03
t(richness_Mountain0.03)->richness_Mountain0.03
colnames(richness_Mountain0.03)<-c("sample_names_Mountain0.03","sample_richness_Mountain0.03")
richness_Mountain0.03
as.data.frame(richness_Mountain0.03)->richness_Mountain0.03


#'**HAPLOTYPE RICHNESS TABLES, PLOTS AND ANALYSES**'# by SITES

#'*Arachnida*

as.matrix(community_Arachnida0.03_Site)->community_Arachnida0.03_Site
row.names(community_Arachnida0.03_Site)->sample_names_Site0.03
dim(community_Arachnida0.03_Site)->dims_Site0.03
dims_Site0.03
.rowSums (community_Arachnida0.03_Site,dims_Site0.03[1],dims_Site0.03[2])->sample_richness_Site0.03 #'summatory by rows
rbind(sample_names_Site0.03,sample_richness_Site0.03)->richness_Site0.03
t(richness_Site0.03)->richness_Site0.03
colnames(richness_Site0.03)<-c("sample_names_Site0.03","sample_richness_Site0.03")
richness_Site0.03
as.data.frame(richness_Site0.03)->richness_Site0.03

##'Generating variables with MOUNTAIN. En la tabla rishnees. Generar variable de montana y sitio. 
richness_Mountain0.03 %>% separate(sample_names_Mountain0.03, c("Conservation","Region","Mountain","Site","ID"), sep="_",remove=FALSE)->richness_Mountain0.03
richness_Mountain0.03
richness_Mountain0.03%>% unite(RegionMountain, Region, Mountain, sep="_",remove=FALSE)->richness_Mountain0.03 #'generating a variable combining layer and habitat
richness_Mountain0.03

##'Generating variables with SITE. En la tabla rishnees. Generar variable de montana y sitio. 
richness_Site0.03 %>% separate(sample_names_Site0.03, c("Conservation","Mountain1","Site","ID"), sep="_",remove=FALSE)->richness_Site0.03
richness_Site0.03
richness_Site0.03 %>% unite(Mountain1Site, Mountain1, Site, sep="_",remove=FALSE)->richness_Site0.03 #'generating a variable combining layer and habitat
richness_Site0.03

#By MOUNTAIN
write.table(richness_Mountain0.03, file="../genetic/Data_out/Arachnida/Arachnida3P/richness_Mountain_Arachnida0.03.txt") #'this is the only way I found to be able to work later
read.table("../genetic/Data_out/Arachnida/Arachnida3P/richness_Mountain_Arachnida0.03.txt",header=TRUE)->richness_Mountain0.03

#BY SITE
write.table(richness_Site0.03, file="../genetic/Data_out/Arachnida/Arachnida3P/richness_Site_Arachnida0.03.txt") #'this is the only way I found to be able to work later
read.table("../genetic/Data_out/Arachnida/Arachnida3P/richness_Site_Arachnida0.03.txt",header=TRUE)->richness_Site0.03

##'General plot of richness by sample in MOUNTAIN
#Mycolors<-c("plum", "hotpink", "orange")
barplot(richness_Mountain0.03$sample_richness_Mountain,col=richness_Mountain0.03$RegionMountain,names.arg= richness_Mountain0.03$sample_names_Mountain,las=2, cex.names=0.5, ylab="richness_Mountain0.03", main="H richness_Mountain_Arachnida_0.03")
richness_Mountain0.03 %>% group_by(RegionMountain) %>% summarise(mean(sample_richness_Mountain0.03))

##'General plot of richness by sample in SITE
barplot(richness_Site0.03$sample_richness_Site,col=richness_Site0.03$Mountain1Site,names.arg= richness_Site0.03$sample_names_Site,las=2,cex.names=0.5, ylab="richness_Site0.03", main="H richness_Site_Arachnida_0.03")
richness_Site0.03 %>% group_by(Mountain1Site) %>% summarise(mean(sample_richness_Site0.03))

##'Global richness by MOUNTAIN.  
plot(richness_Mountain0.03$RegionMountain,richness_Mountain0.03$sample_richness_Mountain,ylab="richness_Mountain0.03", main="H richness_Mountain Arachnida_0.03")
kruskal.test(sample_richness_Mountain0.03 ~ RegionMountain, data = richness_Mountain0.03)
posthoc.kruskal.nemenyi.test(x=richness_Mountain0.03$sample_richness_Mountain, g=richness_Mountain0.03$RegionMountain, method="Bonferroni")




##'Global richness by SITE.  
plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site,ylab="richness_Site0.03", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Arachnida_0.03")
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site, g=richness_Site0.03$Mountain1Site, method="Bonferroni")

# Comparison of each group against  

ggplot(richness_Site0.03, aes(x=Mountain1Site, y=sample_richness_Site0.03))+ 
  geom_boxplot() +
  ylim(0,30) +
  geom_signif(comparisons = list(c("NTO_AAB", "NTO_SJH"), c("NTO_AAB", "NTO_TLC"),
                                 c("NTO_ASB", "NTO_SJH"), c("NTO_ASB", "NTO_TLC")), y_position = c(20, 23, 26, 29), textsize = 3) +
  theme_classic()                               
                                 

##################################################


#'**HAPLOTYPE OCURRENCE TABLES AND SINGLETONS**'


#'*Arachnida*

#'General
#'
##'Singletons by MOUNTAIN
##'summatory by rows
colnames(community_Arachnida0.03_Mountain)->h_names_Mountain0.03
dim (community_Arachnida0.03_Mountain)->dims_Mountain0.03
.colSums (community_Arachnida0.03_Mountain,dims_Mountain0.03[1],dims_Mountain0.03[2])->h_ocurrence_Mountain0.03 #'summatory by cols

rbind(h_names_Mountain0.03,h_ocurrence_Mountain0.03)->h_ocurrence_Mountain0.03
t(h_ocurrence_Mountain0.03)->h_ocurrence_Mountain0.03
colnames(h_ocurrence_Mountain0.03)<-c("h_names_Mountain0.03","h_ocurrence_Mountain0.03")
dim(h_ocurrence_Mountain0.03)
write.table(h_ocurrence_Mountain0.03, file="../genetic/Data_out/Arachnida/Arachnida3P/h_ocurrence_Mountain_Arachnida0.03.txt") #'this is the only way I found to be able to work later
read.table("../genetic/Data_out/Arachnida/Arachnida3P/h_ocurrence_Mountain_Arachnida0.03.txt",header=TRUE)->h_ocurrence_Mountain0.03

#' percentege of singletons by sample
h_ocurrence_Mountain0.03
which(h_ocurrence_Mountain0.03$h_ocurrence_Mountain==1)->singletons_Mountain0.03
length(singletons_Mountain0.03)
length(singletons_Mountain0.03)/length(h_ocurrence_Mountain0.03$h_ocurrence_Mountain)*100

#' percentege of h in more than 2 samples
which(h_ocurrence_Mountain0.03$h_ocurrence_Mountain>2)->more2_Mountain0.03
more2_Mountain0.03
length(more2_Mountain0.03)
length(more2_Mountain0.03)/length(h_ocurrence_Mountain0.03$h_ocurrence_Mountain)*100


##'Singletons by SITE
#'summatory by rows
colnames(community_Arachnida0.03_Site)->h_names_Site0.03
dim (community_Arachnida0.03_Site)->dims_Site0.03
.colSums (community_Arachnida0.03_Site,dims_Site0.03[1],dims_Site0.03[2])->h_ocurrence_Site0.03 #'summatory by cols

rbind(h_names_Site0.03,h_ocurrence_Site0.03)->h_ocurrence_Site0.03
t(h_ocurrence_Site0.03)->h_ocurrence_Site0.03
colnames(h_ocurrence_Site0.03)<-c("h_names_Site0.03","h_ocurrence_Site0.03")
dim(h_ocurrence_Site0.03)
write.table(h_ocurrence_Site0.03, file="../genetic/Data_out/Arachnida/Arachnida3P/h_ocurrence_Site_Arachnida0.03.txt") #'this is the only way I found to be able to work later
read.table("../genetic/Data_out/Arachnida/Arachnida3P/h_ocurrence_Site_Arachnida0.03.txt",header=TRUE)->h_ocurrence_Site0.03

#' percentege of singletons by sample
h_ocurrence_Site0.03
which(h_ocurrence_Site0.03$h_ocurrence_Site==1)->singletons_Site0.03
length(singletons_Site0.03)
length(singletons_Site0.03)/length(h_ocurrence_Site0.03$h_ocurrence_Site)*100

#' percentege of h in more than 2 samples
which(h_ocurrence_Site0.03$h_ocurrence_Site>2)->more2_Site0.03
more2_Site0.03
length(more2_Site0.03)
length(more2_Site0.03)/length(h_ocurrence_Site0.03$h_ocurrence_Site)*100

#*********************
#*********************

#'number of singletons by MOUNTAIN
#community_Arachnida_Mountain
#singletons
community_Arachnida0.03_Mountain[,singletons_Mountain0.03]->community_Arachnida0.03_singletons_Mountain
row.names(community_Arachnida0.03_singletons_Mountain)->sample_names_Mountain0.03
dim (community_Arachnida0.03_singletons_Mountain)->dims_Mountain0.03
dims_Mountain0.03
.rowSums (community_Arachnida0.03_singletons_Mountain,dims_Mountain0.03[1],dims_Mountain0.03[2])->sample_richness_singletons_Mountain0.03 #'summatory by rows
rbind(sample_names_Mountain0.03,sample_richness_singletons_Mountain0.03)->richness_singletons_Mountain0.03
t(richness_singletons_Mountain0.03)->richness_singletons_Mountain0.03
colnames(richness_singletons_Mountain0.03)<-c("sample_names_Mountain0.03","sample_richness_singletons_Mountain0.03")
richness_singletons_Mountain0.03
as.data.frame(richness_singletons_Mountain0.03)->richness_singletons_Mountain0.03


#'number of singletons by SITE
#community_Arachnida_Site
#singletons
community_Arachnida0.03_Site[,singletons_Site0.03]->community_Arachnida0.03_singletons_Site
row.names(community_Arachnida0.03_singletons_Site)->sample_names_Site0.03
dim (community_Arachnida0.03_singletons_Site)->dims_Site0.03
dims_Site0.03
.rowSums (community_Arachnida0.03_singletons_Site,dims_Site0.03[1],dims_Site0.03[2])->sample_richness_singletons_Site0.03 #'summatory by rows
rbind(sample_names_Site0.03,sample_richness_singletons_Site0.03)->richness_singletons_Site0.03
t(richness_singletons_Site0.03)->richness_singletons_Site0.03
colnames(richness_singletons_Site0.03)<-c("sample_names_Site0.03","sample_richness_singletons_Site0.03")
richness_singletons_Site0.03
as.data.frame(richness_singletons_Site0.03)->richness_singletons_Site0.03


####################################################
####################################################

#'**ACCUMULATION CURVES AND EXTRAPOLATED RICHNESS**'#
#'
#'*Arachnida* by MOUNTAIN
#'General
specaccum(community_Arachnida0.03_Mountain,"random", permutations=1000)->cum_Mountain0.03
plot(cum_Mountain0.03, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, main="h_Arachnida_Mountain_0.03")
specpool(community_Arachnida0.03_Mountain)->specpool_Mountain0.03
specpool_Mountain0.03$Species/specpool_Mountain0.03$chao*100


#'*Arachnida* by SITE
#'General
specaccum(community_Arachnida0.03_Site,"random", permutations=1000)->cum_Site0.03
plot(cum_Site0.03, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,100), main="h_Arachnida_Site_0.03")
specpool(community_Arachnida0.03_Site)->specpool_Site0.03
specpool_Site0.03$Species/specpool_Site0.03$chao*100



####################################################
####################################################

#'**BETADIVERSITY ORDINATIONS**'# by MOUNTAIN
#'
#'
#'*Arachnida*
#'
#'beta general

beta.multi(community_Arachnida0.03_Mountain, index.family="sorensen")

#'turnover by pairs, nmds, anosim
beta.pair(community_Arachnida0.03_Mountain, index.family="sorensen")->beta.pair_Arachnida_M0.03  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair_Arachnida_M0.03$beta.sim)->MDSbetasim_Arachnida_M0.03 #'NMDS
plot (MDSbetasim_Arachnida_M0.03, main="Arachnida_Mountain_0.03") 
x<- MDSbetasim_Arachnida_M0.03$points[,1]
y<- MDSbetasim_Arachnida_M0.03$points[,2]
text(x, y, pos = 1, cex=0.7, labels = row.names (community_Arachnida0.03_Mountain))

#PLOTS BY MOUNTAIN SIN OUTLIER

###Hasta aqui modificado Nancy. Nota_Nancy: Modificar para colembolos
plot (MDSbetasim_Arachnida_M0.03, main="Arachnida_Mountain_0.03")
with(general_sample_RegionMountain0.03,ordispider(MDSbetasim_Arachnida_M0.03, Mountain, label=T, col="green"))

#Anosim
anosim(beta.pair_Arachnida_M0.03$beta.sim, general_sample_RegionMountain0.03$Mountain, permutations=999)

#'**BETADIVERSITY ORDINATIONS**'# by Sites
#'
#'
#'*Arachnida*
#'
#'beta general

beta.multi(community_Arachnida0.03_Site, index.family="sorensen")

#'turnover by pairs, nmds, anosim
beta.pair(community_Arachnida0.03_Site, index.family="sorensen")->beta.pair_Arachnida_S0.03  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair_Arachnida_S0.03$beta.sim)->MDSbetasim_Arachnida_S0.03 #'NMDS
plot (MDSbetasim_Arachnida_S0.03, main="Arachnida_Site_0.03") 
x<- MDSbetasim_Arachnida_S0.03$points[,1]
y<- MDSbetasim_Arachnida_S0.03$points[,2]
text(x, y, pos = 1, cex=0.7, labels = row.names (community_Arachnida0.03_Site))

plot (MDSbetasim_Arachnida_S0.03, main="Arachnida_Site_0.03")
with(general_sample_Mountain1Site0.03,ordispider(MDSbetasim_Arachnida_S0.03, Site, label=T, col="blue"))

plot (MDSbetasim_Arachnida_S0.03, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex.axis=1.4, cex=1.2, cex.lab=1.4, main="Arachnida_Site_0.03")
with(general_sample_Mountain1Site0.03,ordispider(MDSbetasim_Arachnida_S0.03, Site, label=T, cex.lab=0.9, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0")))

#Anosim
anosim(beta.pair_Arachnida_S0.03$beta.sim, general_sample_Mountain1Site0.03$Site, permutations=999)

plot (MDSbetasim_Arachnida_S0.03, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5, main="Arachnida_Site_0.03")
with(general_sample_Mountain1Site0.03,ordispider(MDSbetasim_Arachnida_S0.03, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
legend("bottomright", inset=.02, legend=c("R= 0.452 ***"),
       col=c("black"), box.lty=0, cex=1.3)

##############################################END############################################################3


