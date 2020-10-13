#'**INITIAL STEPS**'#

#'**Community diversity and composition at the haplotype level of the Arachnida order**'#

#'**Community diversity and composition of Arachnida at the Clustering level 3% lineages**#
#'In excel remove the simbol #' from the names of the table and rename the samples acccording to the code used in the gradient e.g. GRA_S10_D_F_A10

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

#'**TABLES AND COMMUNITY MATRIXES**'# 
################################################################################################################################################################################

###########'open table with names including Region and Site parameters
s2_raw_all <- read.table("../genetic/Data_in/Arachnida/s2_raw_all_Arachnida_threshold.txt", header=TRUE)
dim(s2_raw_all)

#'Table Haplotipos'
###########'remove additional columns and leave only names (of haplotipes), samples and taxa (and threshold in this case)
s2_raw_all[,c(1:52,68)]->s2_raw 
dim(s2_raw) #'51 samples = 51 plus 1 neg (the second neg from DOM_REPS is not there because all 0)
colnames(s2_raw)

###########'Applying the conservative threshold (this is a binary column)
s2_raw[which(s2_raw$conservative_threshold == "1"),]->s2_raw_threshold 
s2_raw_threshold [,1:52]->s2_raw_threshold #'remove threshold col
dim(s2_raw_threshold)
colnames(s2_raw_threshold)

################################################################################################################################################################################
#'transform in present/absence table
s2_raw_threshold->s2_f4_h #NOTA_Nancy: Tengo un subset de Arachnida
s2_f4_h[s2_f4_h>1]<-1 #'2 warning corresponding wiht the columms of the names and taxa

#'checking if there is any row with no presence
s2_f4_h[,2:52]->data_h 
rowSums(data_h)
length(which(rowSums(data_h)!=0))
length(which(rowSums(data_h)==0))

#'*Arachnida*
t(s2_f4_h)->t_s2_f4_h #'trasp
t_s2_f4_h[2:52,]->community_Arachnida_h #NOTA_Nancy: Este numero es importante. Colocar exactamente el numero de "s2_f4[,2:52]->data".
colnames(community_Arachnida_h)<-t_s2_f4_h[1,]
as.data.frame(community_Arachnida_h)->community_Arachnida_h #'trasp including col and row names
#community_Arachnida[-49,]->community_Arachnida #'removing neg
dim(community_Arachnida_h)
community_Arachnida_h[order(row.names(community_Arachnida_h)),]->community_Arachnida_h #'order samples
write.table (community_Arachnida_h, file="../genetic/Data_out/Arachnida/Arachnida_Haplotypes/community_Arachnida_h.txt") #'this is necessary for the format, not able to solve in other way
read.table ("../genetic/Data_out/Arachnida/Arachnida_Haplotypes/community_Arachnida_h.txt")->community_Arachnida_h

#'submatrixes by SITE in Nevado Toluca.
dim(community_Arachnida_h)
community_Arachnida_h[which(str_extract (row.names(community_Arachnida_h), "_NTO_") %in% "_NTO_"),]->community_Arachnida_Site
dim(community_Arachnida_Site)
community_Arachnida_Site[,which(colSums(community_Arachnida_Site)!=0)]->community_Arachnida_Site #'to remove no data colums
dim(community_Arachnida_Site)

####################################################
#BY SITE
#'Generating a general table with names and habitat parameters
row.names(community_Arachnida_Site)->sample_names_Mountain1
as.data.frame(sample_names_Mountain1)->sample_names_Mountain1
sample_names_Mountain1 %>% separate(sample_names_Mountain1, c("Conservation","Mountain1","Site","ID"), sep="_",remove=FALSE)->general_sample_Mountain1Site_h
general_sample_Mountain1Site_h
general_sample_Mountain1Site_h %>% unite(Mountain1andSite, Mountain1,Site, sep="_",remove=FALSE)->general_sample_Mountain1Site_h #'generating a variable combining layer and habitat
general_sample_Mountain1Site_h
write.table(general_sample_Mountain1Site_h, file="../genetic/Data_out/Arachnida/Arachnida_Haplotypes/general_sample_Mountain1Site_h.txt") #'this is the only way I found to be able to work later
read.table("../genetic/Data_out/Arachnida/Arachnida_Haplotypes/general_sample_Mountain1Site_h.txt",header=TRUE)->general_sample_Mountain1Site_h

####################################################
####################################################
#'**HAPLOTYPE RICHNESS TABLES, PLOTS AND ANALYSES by SITES**'# 
#'*Arachnida*
as.matrix(community_Arachnida_Site)->community_Arachnida_Site
row.names(community_Arachnida_Site)->sample_names_Site
dim(community_Arachnida_Site)->dims_Site
dims_Site
.rowSums (community_Arachnida_Site,dims_Site[1],dims_Site[2])->sample_richness_Site #'summatory by rows
rbind(sample_names_Site,sample_richness_Site)->richness_Site
t(richness_Site)->richness_Site
colnames(richness_Site)<-c("sample_names_Site","sample_richness_Site")
richness_Site
as.data.frame(richness_Site)->richness_Site

##'Generating variables with SITE. En la tabla rishnees. Generar variable de montana y sitio. 
richness_Site %>% separate(sample_names_Site, c("Conservation","Mountain1","Site","ID"), sep="_",remove=FALSE)->richness_Site
richness_Site
richness_Site %>% unite(Mountain1Site, Mountain1, Site, sep="_",remove=FALSE)->richness_Site #'generating a variable combining layer and habitat
richness_Site

#BY SITE
write.table(richness_Site, file="../genetic/Data_out/Arachnida/Arachnida_Haplotypes/richness_Site_Arachnida_h.txt") #'this is the only way I found to be able to work later
read.table("../genetic/Data_out/Arachnida/Arachnida_Haplotypes/richness_Site_Arachnida_h.txt",header=TRUE)->richness_Site

##'General plot of richness by sample in SITE
barplot(richness_Site$sample_richness_Site,col=richness_Site$Mountain1Site,names.arg= richness_Site$sample_names_Site, las=2,cex.names=0.5, ylab="richness_Site", main="H richness_Site Arachnida")
richness_Site %>% group_by(Mountain1Site) %>% summarise(mean(sample_richness_Site))

##'Global richness by SITE.  
plot(richness_Site$Mountain1Site,richness_Site$sample_richness_Site,ylab="richness_Site", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Arachnida")
kruskal.test(sample_richness_Site ~ Mountain1Site, data = richness_Site)
posthoc.kruskal.nemenyi.test(x=richness_Site$sample_richness_Site, g=richness_Site$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
text(x=c(1,2,3,4), y=(28), labels=c("a","a","b","b"), cex=1.4)
text(x=4.5, y=29.5, labels="***", cex=2)

####################################################
#'**HAPLOTYPE OCURRENCE TABLES AND SINGLETONS by SITE**'
#'*Arachnida*
#'General
##'Singletons by SITE
#'summatory by rows
colnames(community_Arachnida_Site)->h_names_Site
dim (community_Arachnida_Site)->dims_Site
.colSums (community_Arachnida_Site,dims_Site[1],dims_Site[2])->h_ocurrence_Site #'summatory by cols

rbind(h_names_Site,h_ocurrence_Site)->h_ocurrence_Site
t(h_ocurrence_Site)->h_ocurrence_Site
colnames(h_ocurrence_Site)<-c("h_names_Site","h_ocurrence_Site")
dim(h_ocurrence_Site)
write.table(h_ocurrence_Site, file="../genetic/Data_out/Arachnida/Arachnida_Haplotypes/h_ocurrence_Site_Arachnida.txt") #'this is the only way I found to be able to work later
read.table("../genetic/Data_out/Arachnida/Arachnida_Haplotypes/h_ocurrence_Site_Arachnida.txt",header=TRUE)->h_ocurrence_Site

#' percentege of singletons by sample
h_ocurrence_Site
which(h_ocurrence_Site$h_ocurrence_Site==1)->singletons_Site
length(singletons_Site)
length(singletons_Site)/length(h_ocurrence_Site$h_ocurrence_Site)*100

#' percentege of h in more than 2 samples
which(h_ocurrence_Site$h_ocurrence_Site>2)->more2_Site
more2_Site
length(more2_Site)
length(more2_Site)/length(h_ocurrence_Site$h_ocurrence_Site)*100

#'number of singletons by SITE
#community_Arachnida_Site
#singletons
community_Arachnida_Site[,singletons_Site]->community_Arachnida_singletons_Site
row.names(community_Arachnida_singletons_Site)->sample_names_Site
dim (community_Arachnida_singletons_Site)->dims_Site
dims_Site
.rowSums (community_Arachnida_singletons_Site,dims_Site[1],dims_Site[2])->sample_richness_singletons_Site #'summatory by rows
rbind(sample_names_Site,sample_richness_singletons_Site)->richness_singletons_Site
t(richness_singletons_Site)->richness_singletons_Site
colnames(richness_singletons_Site)<-c("sample_names_Site","sample_richness_singletons_Site")
richness_singletons_Site
as.data.frame(richness_singletons_Site)->richness_singletons_Site

####################################################
#'*Arachnida* by SITE
#'General
specaccum(community_Arachnida_Site,"random", permutations=1000)->cum_Site
plot(cum_Site, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, main="h_Arachnida_Site")
specpool(community_Arachnida_Site)->specpool_Site
specpool_Site$Species/specpool_Site$chao*100

####################################################
#'**BETADIVERSITY ORDINATIONS**'# by Sites
#'*Arachnida*
#'beta general
beta.multi(community_Arachnida_Site, index.family="sorensen")

#'turnover by pairs, nmds, anosim
beta.pair(community_Arachnida_Site, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim #'NMDS
plot (MDSbetasim, main="Arachnida_Site") 
x<- MDSbetasim$points[,1]
y<- MDSbetasim$points[,2]
text(x, y, pos = 1, cex=0.7, labels = row.names (community_Arachnida_Site))

plot (MDSbetasim, main="Arachnida_Site")
with(general_sample_Mountain1Site_h,ordispider(MDSbetasim, Site, label=T, col="blue"))

plot (MDSbetasim, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex.axis=1.4, cex=1.2, cex.lab=1.4, main="Arachnida_Site_h")
with(general_sample_Mountain1Site_h,ordispider(MDSbetasim, Site, label=T, cex.lab=0.9, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0")))

#Anosim
anosim(beta.pair$beta.sim, general_sample_Mountain1Site_h$Site, permutations=999)

plot (MDSbetasim, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5, main="Arachnida_Site_h")
with(general_sample_Mountain1Site_h,ordispider(MDSbetasim, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))

#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.30, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.569 ***", cex=2)

##################################### E N D ##################################################

