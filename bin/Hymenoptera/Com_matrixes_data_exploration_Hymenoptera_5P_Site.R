#**INITIAL STEPS**

#**Community diversity and composition at the 3% Clustering level of the Hymenoptera order**

##In excel remove the simbol # from the names of the table and rename the samples acccording to the code used in the gradient e.g. GRA_S10_D_F_A10
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
library(permute) 
library(lattice)

#**TABLES AND COMMUNITY MATRIXES** 
#####open table with names including Region and habitat parameters
s2_raw_all <- read.table("../genetic/Data_in/Hymenoptera/s2_raw_all_Hymenoptera_threshold.txt", header=TRUE)
dim(s2_raw_all)

#####remove additional columns and leave only names (of haplotipes), samples and taxa (and threshold in this case)
s2_raw_all[,c(1:65,67)]->s2_raw
dim(s2_raw) ##49 samples = 48 plus 1 neg (the second neg from DOM_REPS is not there because all 0)
colnames(s2_raw)

#####Applying the conservative threshold (this is a binary column)
s2_raw[which(s2_raw$conservative_threshold == "1"),]->s2_raw_threshold 
s2_raw_threshold [,1:65]->s2_raw_threshold ##remove threshold col
dim(s2_raw_threshold)
colnames(s2_raw_threshold)

####loop to create the new matrix combining haplotype by otu pertenency, i.e. submatrix by limit
##**Hymenoptera**
unique (s2_raw_threshold$limite0.05)->levels_limite0.05

data.frame()->s2_raw_Hymenoptera_limite0.05

for (i in 1:length(unique (s2_raw_threshold$limite0.05)))
{
  levels_limite0.05[i]->level0.05
  s2_raw_threshold[which(s2_raw_threshold$limite0.05==level0.05),]->subcom_level_names
  subcom_level_names[,c(2:51)]->subcom_level0.05  #delete names, level and also the negative column
  colSums(subcom_level0.05)->suml0.05
  as.data.frame(suml0.05)->suml0.05
  t(suml0.05)->suml0.05
  row.names(suml0.05)<-subcom_level_names[1,1] #keep the name of the first haplotype
  rbind(s2_raw_Hymenoptera_limite0.05,suml0.05)->s2_raw_Hymenoptera_limite0.05
}

##**transform in present/absence table**
s2_raw_Hymenoptera_limite0.05->s2_raw_Hymenoptera_limite0.05
s2_raw_Hymenoptera_limite0.05[s2_raw_Hymenoptera_limite0.05>1]<-1 ##transform in present/absence table 

##**checking if there is any row with no presence**
s2_raw_Hymenoptera_limite0.05[,1:50]->data0.05
rowSums(data0.05)
length(which(rowSums(data0.05)!=0))
length(which(rowSums(data0.05)==0))

##**Community matrixes (samples in rows and h in cols)**
##**Hymenoptera**
t(s2_raw_Hymenoptera_limite0.05)->t_s2_f4_Hymenoptera_limite0.05 ##trasp
t_s2_f4_Hymenoptera_limite0.05[1:50,]->community_Hymenoptera_limite0.05 #NOTA_Nancy: Este numero es importante. Colocar exactamente el numero de "s2_f4[,2:52]->data0.05".
colnames(t_s2_f4_Hymenoptera_limite0.05)<-community_Hymenoptera_limite0.05[1,] 
as.data.frame(community_Hymenoptera_limite0.05)->community_Hymenoptera0.05 ##trasp including col and row names
####community_Acari[-49,]->community_Hymenoptera ##removing neg
dim(community_Hymenoptera0.05)
community_Hymenoptera0.05[order(row.names(community_Hymenoptera0.05)),]->community_Hymenoptera0.05 ##order samples
write.table (community_Hymenoptera0.05, file="../genetic/Data_out/Hymenoptera/Hymenoptera5P/community_Hymenoptera0.05.txt") ##this is necessary for the format, not able to solve in other way
read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera5P/community_Hymenoptera0.05.txt")->community_Hymenoptera0.05

##**submatrixes by SITE in Nevado Toluca**
dim(community_Hymenoptera0.05)
community_Hymenoptera0.05[which(str_extract (row.names(community_Hymenoptera0.05), "_NTO_") %in% "_NTO_"),]->community_Hymenoptera_Site0.05
dim(community_Hymenoptera_Site0.05)
community_Hymenoptera_Site0.05[,which(colSums(community_Hymenoptera_Site0.05)!=0)]->community_Hymenoptera_Site0.05 ##to remove no data colums
dim(community_Hymenoptera_Site0.05)
write.table (community_Hymenoptera_Site0.05, file="../genetic/Data_out/Hymenoptera/Hymenoptera5P/community_Hymenoptera_Site0.05.txt") ##this is necessary for the format, not able to solve in other way
read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera5P/community_Hymenoptera_Site0.05.txt")->community_Hymenoptera_Site0.05

##**Generating a general table with names and habitat parameters.**
##**BY SITE**
row.names(community_Hymenoptera_Site0.05)->sample_names_Mountain1_0.05
as.data.frame(sample_names_Mountain1_0.05)->sample_names_Mountain1_0.05
sample_names_Mountain1_0.05 %>% separate(sample_names_Mountain1_0.05, c("Conservation","Mountain1","Site","ID"), sep="_",remove=FALSE)->general_sample_Mountain1Site0.05
general_sample_Mountain1Site0.05
general_sample_Mountain1Site0.05 %>% unite(Mountain1andSite, Mountain1,Site, sep="_",remove=FALSE)->general_sample_Mountain1Site0.05 ##generating a variable combining layer and habitat
general_sample_Mountain1Site0.05
write.table(general_sample_Mountain1Site0.05, file="../genetic/Data_out/Hymenoptera/Hymenoptera5P/general_sample_Mountain1Site0.05.txt") ##this is the only way I found to be able to work later
read.table("../genetic/Data_out/Hymenoptera/Hymenoptera5P/general_sample_Mountain1Site0.05.txt",header=TRUE)->general_sample_Mountain1Site0.05
#
#**HAPLOTYPE RICHNESS TABLES, PLOTS AND ANALYSES by SITES** 
##**Hymenoptera**
as.matrix(community_Hymenoptera_Site0.05)->community_Hymenoptera_Site0.05
row.names(community_Hymenoptera_Site0.05)->sample_names_Site0.05
dim(community_Hymenoptera_Site0.05)->dims_Site0.05
dims_Site0.05
.rowSums (community_Hymenoptera_Site0.05,dims_Site0.05[1],dims_Site0.05[2])->sample_richness_Site0.05 ##summatory by rows
rbind(sample_names_Site0.05,sample_richness_Site0.05)->richness_Site0.05
t(richness_Site0.05)->richness_Site0.05
colnames(richness_Site0.05)<-c("sample_names_Site0.05","sample_richness_Site0.05")
richness_Site0.05
as.data.frame(richness_Site0.05)->richness_Site0.05

##**Generating variables with SITE. En la tabla rishnees. Generar variable de montana y sitio.** 
richness_Site0.05 %>% separate(sample_names_Site0.05, c("Conservation","Mountain1","Site","ID"), sep="_",remove=FALSE)->richness_Site0.05
richness_Site0.05
richness_Site0.05 %>% unite(Mountain1Site, Mountain1, Site, sep="_",remove=FALSE)->richness_Site0.05 ##generating a variable combining layer and habitat
richness_Site0.05

##**Generating variables with in total SITE_C.**
richness_Site0.05 %>% separate(sample_names_Site0.05, c("Conservation","Mountain1","Site","ID"), sep="_",remove=FALSE)->richness_SiteC
richness_SiteC
richness_SiteC %>% unite(ConservationMountain1, Conservation, Mountain1, sep="_",remove=FALSE)->richness_SiteC ##generating a variable combining layer and habitat
richness_SiteC

##**BY SITE**
write.table(richness_Site0.05, file="../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_Site0.05_Hymenoptera.txt") ##this is the only way I found to be able to work later
read.table("../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_Site0.05_Hymenoptera.txt",header=TRUE)->richness_Site0.05

##**BY SITE_C general**
write.table(richness_SiteC, file="../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_SiteC_Hymenoptera.txt") ##this is the only way I found to be able to work later
read.table("../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_SiteC_Hymenoptera.txt",header=TRUE)->richness_SiteC

##**General plot of richness by sample in SITE**
barplot(richness_Site0.05$sample_richness_Site0.05,col=richness_Site0.05$Mountain1Site,names.arg= richness_Site0.05$sample_names_Site0.05,las=2,cex.names=0.5, ylab="richness_Site0.05", main="H richness_Site Hymenoptera_0.05")
richness_Site0.05 %>% group_by(Mountain1Site) %>% summarise(mean(sample_richness_Site0.05))

##**min, max, ds Summarise**
richness_Site0.05 %>% group_by(Mountain1Site) %>% summarise(min(sample_richness_Site0.05))
richness_Site0.05 %>% group_by(Mountain1Site) %>% summarise(max(sample_richness_Site0.05))
richness_Site0.05 %>% group_by(Mountain1Site) %>% summarise(sd(sample_richness_Site0.05))

##**General mean, min, max, ds by sample in SITE_richness_SiteC**
richness_SiteC %>% group_by(ConservationMountain1) %>% summarise(mean(sample_richness_Site0.05))
richness_SiteC %>% group_by(ConservationMountain1) %>% summarise(min(sample_richness_Site0.05))
richness_SiteC %>% group_by(ConservationMountain1) %>% summarise(max(sample_richness_Site0.05))
richness_SiteC %>% group_by(ConservationMountain1) %>% summarise(sd(sample_richness_Site0.05))

##**Global richness by SITE.**  
plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05,ylab="richness_Site0.05", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Hymenoptera_0.05")
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
text(x=4.5, y=29, labels="ns", cex=2)
mtext(c("Hymenoptera"), side = 4, col = "black", line = 1, cex = 2)
#

#**HAPLOTYPE OCURRENCE TABLES AND SINGLETONS by SITE**
##**Hymenoptera**
##General
##Singletons 
##summatory by rows
colnames(community_Hymenoptera_Site0.05)->h_names_Site0.05
dim (community_Hymenoptera_Site0.05)->dims_Site0.05
.colSums (community_Hymenoptera_Site0.05,dims_Site0.05[1],dims_Site0.05[2])->h_ocurrence_Site0.05 ##summatory by cols

rbind(h_names_Site0.05,h_ocurrence_Site0.05)->h_ocurrence_Site0.05
t(h_ocurrence_Site0.05)->h_ocurrence_Site0.05
colnames(h_ocurrence_Site0.05)<-c("h_names_Site0.05","h_ocurrence_Site0.05")
dim(h_ocurrence_Site0.05)
write.table(h_ocurrence_Site0.05, file="../genetic/Data_out/Hymenoptera/Hymenoptera5P/h_ocurrence_Site0.05_Hymenoptera.txt") ##this is the only way I found to be able to work later
read.table("../genetic/Data_out/Hymenoptera/Hymenoptera5P/h_ocurrence_Site0.05_Hymenoptera.txt",header=TRUE)->h_ocurrence_Site0.05

##**percentege of singletons by sample**
h_ocurrence_Site0.05
which(h_ocurrence_Site0.05$h_ocurrence_Site0.05==1)->singletons_Site0.05
length(singletons_Site0.05)
length(singletons_Site0.05)/length(h_ocurrence_Site0.05$h_ocurrence_Site0.05)*100

##**percentege of h in more than 2 samples**
which(h_ocurrence_Site0.05$h_ocurrence_Site0.05>2)->more2_Site0.05
more2_Site0.05
length(more2_Site0.05)
length(more2_Site0.05)/length(h_ocurrence_Site0.05$h_ocurrence_Site0.05)*100

##**number of singletons by SITE**
##**community_Hymenoptera_Site0.05**
##singletons
community_Hymenoptera_Site0.05[,singletons_Site0.05]->community_Hymenoptera_singletons_Site0.05
row.names(community_Hymenoptera_singletons_Site0.05)->sample_names_Site0.05
dim (community_Hymenoptera_singletons_Site0.05)->dims_Site0.05
dims_Site0.05
.rowSums (community_Hymenoptera_singletons_Site0.05,dims_Site0.05[1],dims_Site0.05[2])->sample_richness_singletons_Site0.05 ##summatory by rows
rbind(sample_names_Site0.05,sample_richness_singletons_Site0.05)->richness_singletons_Site0.05
t(richness_singletons_Site0.05)->richness_singletons_Site0.05
colnames(richness_singletons_Site0.05)<-c("sample_names_Site0.05","sample_richness_singletons_Site0.05")
richness_singletons_Site0.05
as.data.frame(richness_singletons_Site0.05)->richness_singletons_Site0.05
#
#**ACCUMULATION CURVES AND EXTRAPOLATED RICHNESS by SITE**
##**Hymenoptera**
##General
specaccum(community_Hymenoptera_Site0.05,"random", permutations=1000)->cum_Site0.05
plot(cum_Site0.05, cex=1.4, cex.lab=1.4, cex.axis=2.3, lwd=3, ylim=c(0,100), main="h_Hymenoptera_Site_0.05")
specpool(community_Hymenoptera_Site0.05)->specpool_Site0.05
specpool_Site0.05$Species/specpool_Site0.05$chao*100
text(x=38, y=5, labels="69.42%", cex=1.5)
mtext(c("Hymenoptera"), side = 4, col = "black", line = 1, cex = 2)
#

#**BETADIVERSITY ORDINATIONS by Sites** 
##**Hymenoptera**
##beta general
beta.multi(community_Hymenoptera_Site0.05, index.family="sorensen")

##**By SITES quitar 122, 124, 36**
##**repeating after removing outlayers from matrix and general habitat table**
community_Hymenoptera_Site0.05[-which(row.names(community_Hymenoptera_Site0.05) %in% c("CON_NTO_AAB_122ACON12", "CON_NTO_AAB_124ACON14", "CON_NTO_TLC_36TCONS6")),]->community_Hymenoptera_sinoutlayer0.05
community_Hymenoptera_sinoutlayer0.05[,which(colSums(community_Hymenoptera_sinoutlayer0.05)!=0)]->community_Hymenoptera_sinoutlayer0.05 ##to remove no data colums
dim(community_Hymenoptera_sinoutlayer0.05)  

general_sample_Mountain1Site0.05[-which(general_sample_Mountain1Site0.05$sample_names %in% c("CON_NTO_AAB_122ACON12", "CON_NTO_AAB_124ACON14", "CON_NTO_TLC_36TCONS6")),]->general_sample_sinoutlayer0.05

beta.pair(community_Hymenoptera_sinoutlayer0.05, index.family="sorensen")->beta.pair  ##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.05

##general_sample_sinoutlayer0.05
##separate(sample_names, c("Conservation","Mountain","Site","ID")
##plot (MDSbetasim0.05, main="Hymenoptera_Site")
##with(general_sample_sinoutlayer0.05,ordispider(MDSbetasim0.05, Mountain, label=T, col="green"))
plot (MDSbetasim0.05, main="Hymenoptera_Site_h")
with(general_sample_sinoutlayer0.05,ordispider(MDSbetasim0.05, Site, label=T, col="blue"))
##plot (MDSbetasim0.05, main="Hymenoptera_Site")
##with(general_sample_sinoutlayer0.05,ordispider(MDSbetasim0.05, ID, label=T, col="red"))
plot (MDSbetasim0.05, main="Hymenoptera_Site_0.05")
x<- MDSbetasim0.05$points[,1]
y<- MDSbetasim0.05$points[,2]
text(x, y, pos = 1, cex=0.7, labels = row.names (community_Hymenoptera_sinoutlayer0.05))

plot (MDSbetasim0.05, main="Hymenoptera_Site_0.05")
with(general_sample_sinoutlayer0.05,ordispider(MDSbetasim0.05, Site, label=T, col="blue"))

plot (MDSbetasim0.05, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), cex.axis=1.4, cex=1.2, cex.lab=1.4, main="Hymenoptera_Site_0.05")
with(general_sample_sinoutlayer0.05,ordispider(MDSbetasim0.05, Site, label=T, cex.lab=0.9, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0")))

plot (MDSbetasim0.05, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5, main="Hymenoptera_Site_0.05")
with(general_sample_sinoutlayer0.05,ordispider(MDSbetasim0.05, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))

##**Anosim**
anosim(beta.pair$beta.sim, general_sample_sinoutlayer0.05$Site, permutations=999)

##I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.16, y=-0.42, labels = mylabel, cex=2)
text(x=0.56, y=-0.42, labels="=0.154 ***", cex=2)
mtext(c("Hymenoptera"), side = 4, col = "black", line = 1, cex = 2)

#**END**
