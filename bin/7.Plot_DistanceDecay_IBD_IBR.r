#**Creating plot**
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
library(ggplot2)
library(patchwork)
library(easypackages)

##**We get the Beta diversity vs Geographic distances, flat, Altitude 3000 and Altitude B from haplotype levels to 7.5% CL for each of the two taxonomic orders studied**

###This script gets all scripts of diversity using 2 arthropods order at multi-hierarchical levels. We used 2 groups: Collembola and Diptera.

#**plot Beta diversity vs conductace flat and Altitude**

png(filename="../figures/Figure5_DistanceDecay_Flat_Altitude.png", width=787, height=467, units="px") # set size of the file to plot 
par(mfrow=c(2,2), mar = c(1.5, 1.5, 1.5, 1.5), omi=c(0.5, 0.5, 0.5, 2)) #the number of rows and columns the figure would have

#**Creating plot**
##Collembola Flat Site

read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_Site/community_Collembola_h_Site.txt")->community_Collembola_h_Site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_Site/community_Collembola0.005_Site.txt")->community_Collembola0.005_Site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_Site/community_Collembola0.015_Site.txt")->community_Collembola0.015_Site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_Site/community_Collembola0.02_Site.txt")->community_Collembola0.02_Site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_Site/community_Collembola0.03_Site.txt")->community_Collembola0.03_Site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_Site/community_Collembola0.05_Site.txt")->community_Collembola0.05_Site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_Site/community_Collembola0.075_Site.txt")->community_Collembola0.075_Site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_Site/community_Collembola0.029_Site.txt")->community_Collembola0.029_Site

#**BETADIVERSITY ORDINATIONS by SITE** 
##**Collembola**

##**beta general_Level_Haplotipos**
beta.pair(community_Collembola_h_Site, index.family="sorensen")->beta.pair_CollembolaSite_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.005_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.005

##**beta general_Level_0.015**
###betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.015_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.015

##**beta general_Level_0.02**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.02_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.03_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.03

##**beta general_Level_0.05**
###betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.05_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.075_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.075

##**beta general_Level_0.029**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.029_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.029

##**FORMAS DE OBTENER LA MATRIX effective resistance*
Resistance_matrix_Site <- read.table("../spatial/IBResistanceFlatMatrix/Result_flat_resistances_CollembolaSite.txt",sep = ",", header=T, row.names = 1)
dim(Resistance_matrix_Site)
class(Resistance_matrix_Site)
Resistance_matrix_Site <- as.matrix(Resistance_matrix_Site)
class(Resistance_matrix_Site)

Resistance_matrix_Site[order(row.names(Resistance_matrix_Site)),order(colnames(Resistance_matrix_Site))]->Resistance_matrix_Site #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_Site[upper.tri(Resistance_matrix_Site)] <- NA 
Resistance_matrix_Site
class(Resistance_matrix_Site)

Resistance_matrix_Site <- as.dist(Resistance_matrix_Site) 
Resistance_matrix_Site

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_CollembolaSite_h$beta.sim +0.001->all_h_betasim
1-beta.pair_CollembolaSite_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_CollembolaSite_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_CollembolaSite_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_CollembolaSite_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_CollembolaSite_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_CollembolaSite_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_CollembolaSite_0.029$beta.sim+0.001->all_0.029_betasim

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.029_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.029_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.029

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5, GMYC**
##**Plot all levels**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(0.4,1.35), pch=20, lwd=4, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
##plot.decay(decay_0.02,add=T,pch=20,lwd=4,col="#d49e57")
plot.decay(decay_0.029,add=T,pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.03,add=T,pch=20,lty=3,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff")
text(x=0.41, y=0.0, labels="(a)", cex=1.8)

#**Creating plot**
##Collembola Altitude 3000 Site

read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_Site/community_Collembola_h_Site.txt")->community_Collembola_h_Site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_Site/community_Collembola0.005_Site.txt")->community_Collembola0.005_Site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_Site/community_Collembola0.015_Site.txt")->community_Collembola0.015_Site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_Site/community_Collembola0.02_Site.txt")->community_Collembola0.02_Site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_Site/community_Collembola0.03_Site.txt")->community_Collembola0.03_Site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_Site/community_Collembola0.05_Site.txt")->community_Collembola0.05_site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_Site/community_Collembola0.075_Site.txt")->community_Collembola0.075_Site
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_Site/community_Collembola0.029_Site.txt")->community_Collembola0.029_Site

#**BETADIVERSITY ORDINATIONS by SITE**
##**Collembola**

##**beta general_Level_Haplotipos**
beta.pair(community_Collembola_h_Site, index.family="sorensen")->beta.pair_CollembolaSite_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.005_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.015_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.015

##**beta general_Level_0.02**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.02_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.03_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.03

##**beta general_Level_0.05**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.05_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.075_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.075

##**beta general_Level_0.029**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.029_Site, index.family="sorensen")->beta.pair_CollembolaSite_0.029

##**FORMAS DE OBTENER LA MATRIX effective resistance** 
Resistance_matrix_Site <- read.table("../spatial/IBResitanceElevationMatrix/Result_3000_resistances.txt",sep = ",", header=T, row.names = 1)
dim(Resistance_matrix_Site)
class(Resistance_matrix_Site)
Resistance_matrix_Site <- as.matrix(Resistance_matrix_Site)
class(Resistance_matrix_Site)

Resistance_matrix_Site[order(row.names(Resistance_matrix_Site)),order(colnames(Resistance_matrix_Site))]->Resistance_matrix_Site #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_Site[upper.tri(Resistance_matrix_Site)] <- NA 
Resistance_matrix_Site
class(Resistance_matrix_Site)

Resistance_matrix_Site <- as.dist(Resistance_matrix_Site) 
Resistance_matrix_Site

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_CollembolaSite_h$beta.sim +0.001->all_h_betasim
1-beta.pair_CollembolaSite_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_CollembolaSite_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_CollembolaSite_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_CollembolaSite_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_CollembolaSite_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_CollembolaSite_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_CollembolaSite_0.029$beta.sim+0.001->all_0.029_betasim

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.029_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.029_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.029

#**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5**
##"H", "NC0.5", "NC1.5", "NC3","NC5", "NC7.5"
##**Plot all levels**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(0.4,1.7), pch=20, lwd=4, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
##plot.decay(decay_0.02,add=T,pch=16,lwd=4,col="#d49e57")
plot.decay(decay_0.029,add=T,pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.03,add=T,pch=20,lty=3,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff") 
text(x=0.41, y=0.0, labels="(b)", cex=1.8)
mtext(c("Collembola"), side = 4, adj = 0, col = "black", line = 2, cex=1.8)

#**Creating plot**
##Diptera Flat Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_Site/community_Diptera_h_Site.txt")->community_Diptera_h_Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_Site/community_Diptera0.005_Site.txt")->community_Diptera0.005_Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_Site/community_Diptera0.015_Site.txt")->community_Diptera0.015_Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_Site/community_Diptera0.02_Site.txt")->community_Diptera0.02_Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_Site/community_Diptera0.03_Site.txt")->community_Diptera0.03_Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_Site/community_Diptera0.05_Site.txt")->community_Diptera0.05_Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_Site/community_Diptera0.075_Site.txt")->community_Diptera0.075_Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_Site/community_Diptera0.0088_Site.txt")->community_Diptera0.0088_Site

#**BETADIVERSITY ORDINATIONS by SITE**
##**Diptera**

##**beta general_Level_Haplotipos**
beta.pair(community_Diptera_h_Site, index.family="sorensen")->beta.pair_DipteraSite_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.005_Site, index.family="sorensen")->beta.pair_DipteraSite_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.015_Site, index.family="sorensen")->beta.pair_DipteraSite_0.015

##**beta general_Level_0.02**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.02_Site, index.family="sorensen")->beta.pair_DipteraSite_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.03_Site, index.family="sorensen")->beta.pair_DipteraSite_0.03

##**beta general_Level_0.05**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.05_Site, index.family="sorensen")->beta.pair_DipteraSite_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.075_Site, index.family="sorensen")->beta.pair_DipteraSite_0.075

##**beta general_Level_0.0088**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.0088_Site, index.family="sorensen")->beta.pair_DipteraSite_0.0088

##**FORMAS DE OBTENER LA MATRIX effective resistance** 
Resistance_matrix_Site <- read.table("../spatial/IBResistanceFlatMatrix/Result_flat_resistances_DipteraSite.txt",sep = ",", header=T, row.names = 1)
dim(Resistance_matrix_Site)
class(Resistance_matrix_Site)
Resistance_matrix_Site <- as.matrix(Resistance_matrix_Site)
class(Resistance_matrix_Site)

Resistance_matrix_Site[order(row.names(Resistance_matrix_Site)),order(colnames(Resistance_matrix_Site))]->Resistance_matrix_Site #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_Site[upper.tri(Resistance_matrix_Site)] <- NA 
Resistance_matrix_Site
class(Resistance_matrix_Site)

Resistance_matrix_Site <- as.dist(Resistance_matrix_Site) 
Resistance_matrix_Site

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_DipteraSite_h$beta.sim +0.001->all_h_betasim
1-beta.pair_DipteraSite_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_DipteraSite_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_DipteraSite_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_DipteraSite_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_DipteraSite_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_DipteraSite_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_DipteraSite_0.0088$beta.sim+0.001->all_0.0088_betasim

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.0088_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.0088_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.0088

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5, GMYC**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(0.4,1.35), pch=20, lwd=4, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.0088,add=T, pch=20, lwd=4, col="#92000A")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
##plot.decay(decay_0.02,add=T,pch=20,lwd=4,col="#d49e57")
plot.decay(decay_0.03,add=T,pch=20,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff")
text(x=0.41, y=0.0, labels="(c)", cex=1.8)
mtext(c("Effective distance flat"), side = 1, col = "black", line = 2.5, cex=1.5)
#
#**Creating plot**
##Diptera Altitude B Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_Site/community_Diptera_h_Site.txt")->community_Diptera_h_Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_Site/community_Diptera0.005_Site.txt")->community_Diptera0.005_Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_Site/community_Diptera0.015_Site.txt")->community_Diptera0.015_Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_Site/community_Diptera0.02_Site.txt")->community_Diptera0.02_Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_Site/community_Diptera0.03_Site.txt")->community_Diptera0.03_Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_Site/community_Diptera0.05_Site.txt")->community_Diptera0.05_Site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_Site/community_Diptera0.075_Site.txt")->community_Diptera0.075_site
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_Site/community_Diptera0.0088_Site.txt")->community_Diptera0.0088_Site

#**BETADIVERSITY ORDINATIONS by SITE**
##**Diptera**

##**beta general_Level_Haplotipos**
beta.pair(community_Diptera_h_Site, index.family="sorensen")->beta.pair_DipteraSite_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.005_Site, index.family="sorensen")->beta.pair_DipteraSite_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.015_Site, index.family="sorensen")->beta.pair_DipteraSite_0.015

##**beta general_Level_0.02**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.02_Site, index.family="sorensen")->beta.pair_DipteraSite_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.03_Site, index.family="sorensen")->beta.pair_DipteraSite_0.03

##**beta general_Level_0.05**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.05_Site, index.family="sorensen")->beta.pair_DipteraSite_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.075_Site, index.family="sorensen")->beta.pair_DipteraSite_0.075

##**beta general_Level_0.0088**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.0088_Site, index.family="sorensen")->beta.pair_DipteraSite_0.0088

##**FORMAS DE OBTENER LA MATRIX effective resistance**
Resistance_matrix_Site <- read.table("../spatial/IBResitanceElevationMatrix/Result_Alt_B_resistances.txt",sep = ",", header=T, row.names = 1)
dim(Resistance_matrix_Site)
class(Resistance_matrix_Site)
Resistance_matrix_Site <- as.matrix(Resistance_matrix_Site)
class(Resistance_matrix_Site)

Resistance_matrix_Site[order(row.names(Resistance_matrix_Site)),order(colnames(Resistance_matrix_Site))]->Resistance_matrix_Site #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_Site[upper.tri(Resistance_matrix_Site)] <- NA 
Resistance_matrix_Site
class(Resistance_matrix_Site)

Resistance_matrix_Site <- as.dist(Resistance_matrix_Site) 
Resistance_matrix_Site

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_DipteraSite_h$beta.sim +0.001->all_h_betasim
1-beta.pair_DipteraSite_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_DipteraSite_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_DipteraSite_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_DipteraSite_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_DipteraSite_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_DipteraSite_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_DipteraSite_0.0088$beta.sim+0.001->all_0.0088_betasim

##**Ddecay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.0088_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.0088_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.0088

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5, GMYC**
plot.decay(decay_h, ylim=c(0,1), xlim=c(1.3,5.6), pch=20, lwd=4, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.0088,add=T,pch=20, lwd=4, col="#92000A")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
##plot.decay(decay_0.02,add=T,pch=20,lwd=4,col="#d49e57")
plot.decay(decay_0.03,add=T,pch=20,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff")
text(x=1.31, y=0.0, labels="(d)", cex=1.8)
mtext(c("Diptera"), side = 4, adj = 0, col = "black", line = 2, cex=1.8)
mtext(c("Effective distance (Altitude rasters)"), side = 1, col = "black", line = 2.5, cex=1.5)
par(xpd=TRUE)
###Legend
legend(6.4, 1, xpd=NA,
       legend=c("Haplotype", "CL 0.5", "CL 1.5", "GMYC", "CL 3","CL 5", "CL 7.5"), 
       col=c("#003695", "#c997a9", "#9d9cc6", "#92000A", "#9cb15b", "#fbd048", "#93dfff"), 
       pch=19,  bty="n", text.font=1.7, lty=1, cex=1.7, lwd=2.2) 
       
###Legend
mtext("Similarity (Simpson's Index)", side=2, outer=TRUE, line=0.7, cex=1.6)

dev.off()
#

#**Plot Beta diversity vs conductace flat in finer scale**

png(filename="../figures/Figure6_DistanceDecay_FinerScaleFlat.png", width=787, height=467, units="px") # set size of the file to plot 
par(mfrow=c(2,2), mar = c(1.5, 1.5, 1.5, 1.5), omi=c(0.5, 0.5, 0.5, 2)) #the number of rows and columns the figure would have

#**Creating plot**
##Finer scale Collembola EAST Flat

read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_East/community_Collembola_h_East.txt")->community_Collembola_h_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_East/community_Collembola0.005_East.txt")->community_Collembola0.005_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_East/community_Collembola0.015_East.txt")->community_Collembola0.015_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_East/community_Collembola0.02_East.txt")->community_Collembola0.02_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_East/community_Collembola0.03_East.txt")->community_Collembola0.03_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_East/community_Collembola0.05_East.txt")->community_Collembola0.05_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_East/community_Collembola0.075_East.txt")->community_Collembola0.075_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_East/community_Collembola0.029_East.txt")->community_Collembola0.029_East

#**BETADIVERSITY ORDINATIONS by SITE**
##**Collembola**

##**beta general_Level_Haplotipos**
beta.pair(community_Collembola_h_East, index.family="sorensen")->beta.pair_CollembolaEast_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.005_East, index.family="sorensen")->beta.pair_CollembolaEast_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.015_East, index.family="sorensen")->beta.pair_CollembolaEast_0.015

##**beta general_Level_0.02**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.02_East, index.family="sorensen")->beta.pair_CollembolaEast_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.03_East, index.family="sorensen")->beta.pair_CollembolaEast_0.03

##**beta general_Level_0.05**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.05_East, index.family="sorensen")->beta.pair_CollembolaEast_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.075_East, index.family="sorensen")->beta.pair_CollembolaEast_0.075

##**beta general_Level_0.029**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.029_East, index.family="sorensen")->beta.pair_CollembolaEast_0.029

##**FORMAS DE OBTENER LA MATRIX effective resistance**
Resistance_matrix_East <- read.table("../spatial/IBResistanceFlatMatrix/Result_flat_resistances_CollembolaEast.txt", sep = ",", header=T, row.names = 1)
dim(Resistance_matrix_East)
class(Resistance_matrix_East)
Resistance_matrix_East <- as.matrix(Resistance_matrix_East)
class(Resistance_matrix_East)

Resistance_matrix_East[order(row.names(Resistance_matrix_East)),order(colnames(Resistance_matrix_East))]->Resistance_matrix_East #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_East[upper.tri(Resistance_matrix_East)] <- NA 
Resistance_matrix_East
class(Resistance_matrix_East)

Resistance_matrix_East <- as.dist(Resistance_matrix_East) 
Resistance_matrix_East

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_CollembolaEast_h$beta.sim +0.001->all_h_betasim
1-beta.pair_CollembolaEast_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_CollembolaEast_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_CollembolaEast_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_CollembolaEast_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_CollembolaEast_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_CollembolaEast_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_CollembolaEast_0.029$beta.sim+0.001->all_0.029_betasim

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.029_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.029_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.029

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5, GMYC**
##**Plot all levels**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(0.4,1.20), pch=20, lwd=4, cex.lab= 1.5, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
#plot.decay(decay_0.02,add=T,pch=20,lwd=4,col="#d49e57")
plot.decay(decay_0.029,add=T,pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.03,add=T,pch=20,lty=3,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff")
text(x=0.41, y=0.0, labels="(a)", cex=1.8)

#**Creating plot**
##Finer scale Collembola WEST Flat

read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_West/community_Collembola_h_West.txt")->community_Collembola_h_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_West/community_Collembola0.005_West.txt")->community_Collembola0.005_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_West/community_Collembola0.015_West.txt")->community_Collembola0.015_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_West/community_Collembola0.02_West.txt")->community_Collembola0.02_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_West/community_Collembola0.03_West.txt")->community_Collembola0.03_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_West/community_Collembola0.05_West.txt")->community_Collembola0.05_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_West/community_Collembola0.075_West.txt")->community_Collembola0.075_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_Flat_West/community_Collembola0.029_West.txt")->community_Collembola0.029_West

#**BETADIVERSITY ORDINATIONS by SITE** 
##**Collembola**

##**beta general_Level_Haplotipos**
beta.pair(community_Collembola_h_West, index.family="sorensen")->beta.pair_CollembolaWest_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.005_West, index.family="sorensen")->beta.pair_CollembolaWest_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.015_West, index.family="sorensen")->beta.pair_CollembolaWest_0.015

##**beta general_Level_0.02**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.02_West, index.family="sorensen")->beta.pair_CollembolaWest_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.03_West, index.family="sorensen")->beta.pair_CollembolaWest_0.03

##**beta general_Level_0.05**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.05_West, index.family="sorensen")->beta.pair_CollembolaWest_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.075_West, index.family="sorensen")->beta.pair_CollembolaWest_0.075

#**beta general_Level_0.029**
###betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.029_West, index.family="sorensen")->beta.pair_CollembolaWest_0.029

##**FORMAS DE OBTENER LA MATRIX effective resistance** 
Resistance_matrix_West <- read.table("../spatial/IBResistanceFlatMatrix/Result_flat_resistances_CollembolaWest.txt", sep = ",", header=T, row.names = 1)
dim(Resistance_matrix_West)
class(Resistance_matrix_West)
Resistance_matrix_West <- as.matrix(Resistance_matrix_West)
class(Resistance_matrix_West)

Resistance_matrix_West[order(row.names(Resistance_matrix_West)),order(colnames(Resistance_matrix_West))]->Resistance_matrix_West #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_West[upper.tri(Resistance_matrix_West)] <- NA 
Resistance_matrix_West
class(Resistance_matrix_West)

Resistance_matrix_West <- as.dist(Resistance_matrix_West) 
Resistance_matrix_West

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_CollembolaWest_h$beta.sim +0.001->all_h_betasim
1-beta.pair_CollembolaWest_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_CollembolaWest_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_CollembolaWest_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_CollembolaWest_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_CollembolaWest_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_CollembolaWest_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_CollembolaWest_0.029$beta.sim+0.001->all_0.029_betasim

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.029_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.029_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.029

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5**
##**Plot all levels**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(0.35,1.10), pch=20, lwd=4, cex.lab= 1.5, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
##plot.decay(decay_0.02,add=T,pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.029,add=T, pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.03,add=T,pch=20,lty=3,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff")
text(x=0.36, y=0.0, labels="(b)", cex=1.8)
mtext(c("Collembola"), side = 4, adj = 0, col = "black", line = 2, cex=1.8)

#**Creating plot**
##Finer scale Diptera EAST Flat
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_East/community_Diptera_h_East.txt")->community_Diptera_h_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_East/community_Diptera0.005_East.txt")->community_Diptera0.005_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_East/community_Diptera0.015_East.txt")->community_Diptera0.015_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_East/community_Diptera0.02_East.txt")->community_Diptera0.02_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_East/community_Diptera0.03_East.txt")->community_Diptera0.03_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_East/community_Diptera0.05_East.txt")->community_Diptera0.05_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_East/community_Diptera0.075_East.txt")->community_Diptera0.075_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_East/community_Diptera0.0088_East.txt")->community_Diptera0.0088_East

#**BETADIVERSITY ORDINATIONS by East** 
##**Diptera**

##**beta general_Level_Haplotipos**
beta.pair(community_Diptera_h_East, index.family="sorensen")->beta.pair_DipteraEast_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.005_East, index.family="sorensen")->beta.pair_DipteraEast_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.015_East, index.family="sorensen")->beta.pair_DipteraEast_0.015

##**beta general_Level_0.02**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.02_East, index.family="sorensen")->beta.pair_DipteraEast_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.03_East, index.family="sorensen")->beta.pair_DipteraEast_0.03

##**beta general_Level_0.05**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.05_East, index.family="sorensen")->beta.pair_DipteraEast_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.075_East, index.family="sorensen")->beta.pair_DipteraEast_0.075

#**beta general_Level_0.0088**
###betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.0088_East, index.family="sorensen")->beta.pair_DipteraEast_0.0088

##**FORMAS DE OBTENER LA MATRIX effective resistance** 
Resistance_matrix_East <- read.table("../spatial/IBResistanceFlatMatrix/Result_flat_resistances_DipteraEast.txt", sep = ",", header=T, row.names = 1)
dim(Resistance_matrix_East)
class(Resistance_matrix_East)
Resistance_matrix_East <- as.matrix(Resistance_matrix_East)
class(Resistance_matrix_East)

Resistance_matrix_East[order(row.names(Resistance_matrix_East)),order(colnames(Resistance_matrix_East))]->Resistance_matrix_East #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_East[upper.tri(Resistance_matrix_East)] <- NA 
Resistance_matrix_East
class(Resistance_matrix_East)

Resistance_matrix_East <- as.dist(Resistance_matrix_East) 
Resistance_matrix_East

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_DipteraEast_h$beta.sim +0.001->all_h_betasim
1-beta.pair_DipteraEast_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_DipteraEast_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_DipteraEast_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_DipteraEast_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_DipteraEast_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_DipteraEast_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_DipteraEast_0.0088$beta.sim+0.001->all_0.0088_betasim

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.0088_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.0088_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.0088

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5, GMYC**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(0.4,1.20), pch=20, lwd=4, cex.lab= 1.5, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.0088,add=T,pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
##plot.decay(decay_0.02,add=T,pch=20,lwd=4,col="#d49e57")
plot.decay(decay_0.03,add=T,pch=20,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff")
text(x=0.41, y=0.0, labels="(c)", cex=1.8)
mtext(c("Effective distance East flat"), side = 1, col = "black", line = 2.5, cex=1.5)

#**Creating plot**
##Finer scale Diptera West Flat

read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_West/community_Diptera_h_West.txt")->community_Diptera_h_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_West/community_Diptera0.005_West.txt")->community_Diptera0.005_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_West/community_Diptera0.015_West.txt")->community_Diptera0.015_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_West/community_Diptera0.02_West.txt")->community_Diptera0.02_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_West/community_Diptera0.03_West.txt")->community_Diptera0.03_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_West/community_Diptera0.05_West.txt")->community_Diptera0.05_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_West/community_Diptera0.075_West.txt")->community_Diptera0.075_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_Flat_West/community_Diptera0.0088_West.txt")->community_Diptera0.0088_West

#**BETADIVERSITY ORDINATIONS by West**
##**Diptera**

##**beta general_Level_Haplotipos**
beta.pair(community_Diptera_h_West, index.family="sorensen")->beta.pair_DipteraWest_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.005_West, index.family="sorensen")->beta.pair_DipteraWest_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.015_West, index.family="sorensen")->beta.pair_DipteraWest_0.015

##**beta general_Level_0.02**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.02_West, index.family="sorensen")->beta.pair_DipteraWest_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.03_West, index.family="sorensen")->beta.pair_DipteraWest_0.03

##**beta general_Level_0.05**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.05_West, index.family="sorensen")->beta.pair_DipteraWest_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.075_West, index.family="sorensen")->beta.pair_DipteraWest_0.075

##**beta general_Level_0.0088**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.0088_West, index.family="sorensen")->beta.pair_DipteraWest_0.0088

##**FORMAS DE OBTENER LA MATRIX effective resistance** 
Resistance_matrix_West <- read.table("../spatial/IBResistanceFlatMatrix/Result_flat_resistances_DipteraWest.txt", sep = ",", header=T, row.names = 1)
dim(Resistance_matrix_West)
class(Resistance_matrix_West)
Resistance_matrix_West <- as.matrix(Resistance_matrix_West)
class(Resistance_matrix_West)

Resistance_matrix_West[order(row.names(Resistance_matrix_West)),order(colnames(Resistance_matrix_West))]->Resistance_matrix_West #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_West[upper.tri(Resistance_matrix_West)] <- NA 
Resistance_matrix_West
class(Resistance_matrix_West)

Resistance_matrix_West <- as.dist(Resistance_matrix_West) 
Resistance_matrix_West

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_DipteraWest_h$beta.sim +0.001->all_h_betasim
1-beta.pair_DipteraWest_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_DipteraWest_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_DipteraWest_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_DipteraWest_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_DipteraWest_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_DipteraWest_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_DipteraWest_0.0088$beta.sim+0.001->all_0.0088_betasim

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.0088_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.0088_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.0088

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(0.35,1.10), pch=20, lwd=4, cex.lab= 1.5, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.0088,add=T,pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
##plot.decay(decay_0.02,add=T,pch=20,lwd=4,col="#d49e57")
plot.decay(decay_0.03,add=T,pch=20,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff")
text(x=0.36, y=0.0, labels="(d)", cex=1.8)
mtext(c("Diptera"), side = 4, adj = 0, col = "black", line = 2, cex=1.8)
mtext(c("Effective distance West flat"), side = 1, col = "black", line = 2.5, cex=1.5)
par(xpd=TRUE)
###Legend
legend(1.24, 1, xpd=NA,
       legend=c("Haplotype", "CL 0.5", "CL 1.5", "GMYC", "CL 3","CL 5", "CL 7.5"), 
       col=c("#003695", "#c997a9", "#9d9cc6", "#92000A", "#9cb15b", "#fbd048", "#93dfff"), 
       pch=19,  bty="n", text.font=1.7, lty=1, cex=1.7, lwd=2.2) 
       
###Legend
mtext("Similarity (Simpson's Index)", side=2, outer=TRUE, line=0.7, cex=1.6)

dev.off()
#

#**Plot Beta diversity vs conductace Altitudes in finer scale**

png(filename="../figures/FigureS6_DistanceDecay_FinerScaleAltitudes.png", width=787, height=467, units="px") # set size of the file to plot 
par(mfrow=c(2,2), mar = c(1.5, 1.5, 1.5, 1.5), omi=c(0.5, 0.5, 0.5, 2)) #the number of rows and columns the figure would have

#**Creating plot**
##Finer scale Collembola EAST Altitude

read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_East/community_Collembola_h_East.txt")->community_Collembola_h_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_East/community_Collembola0.005_East.txt")->community_Collembola0.005_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_East/community_Collembola0.015_East.txt")->community_Collembola0.015_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_East/community_Collembola0.02_East.txt")->community_Collembola0.02_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_East/community_Collembola0.03_East.txt")->community_Collembola0.03_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_East/community_Collembola0.05_East.txt")->community_Collembola0.05_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_East/community_Collembola0.075_East.txt")->community_Collembola0.075_East
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_East/community_Collembola0.029_East.txt")->community_Collembola0.029_East

#**BETADIVERSITY ORDINATIONS by SITE**
##**Collembola**

##**beta general_Level_Haplotipos**
beta.pair(community_Collembola_h_East, index.family="sorensen")->beta.pair_CollembolaEast_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.005_East, index.family="sorensen")->beta.pair_CollembolaEast_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.015_East, index.family="sorensen")->beta.pair_CollembolaEast_0.015

##**beta general_Level_0.02**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.02_East, index.family="sorensen")->beta.pair_CollembolaEast_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.03_East, index.family="sorensen")->beta.pair_CollembolaEast_0.03

##**beta general_Level_0.05**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.05_East, index.family="sorensen")->beta.pair_CollembolaEast_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.075_East, index.family="sorensen")->beta.pair_CollembolaEast_0.075

##**beta general_Level_0.029**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.029_East, index.family="sorensen")->beta.pair_CollembolaEast_0.029

##**FORMAS DE OBTENER LA MATRIX effective resistance** 
Resistance_matrix_East <- read.table("../spatial/IBResitanceElevationMatrix/Result_3000_resist_East.txt", sep = ",", header=T, row.names = 1)
dim(Resistance_matrix_East)
class(Resistance_matrix_East)
Resistance_matrix_East <- as.matrix(Resistance_matrix_East)
class(Resistance_matrix_East)

Resistance_matrix_East[order(row.names(Resistance_matrix_East)),order(colnames(Resistance_matrix_East))]->Resistance_matrix_East #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_East[upper.tri(Resistance_matrix_East)] <- NA 
Resistance_matrix_East
class(Resistance_matrix_East)

Resistance_matrix_East <- as.dist(Resistance_matrix_East) 
Resistance_matrix_East

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_CollembolaEast_h$beta.sim +0.001->all_h_betasim
1-beta.pair_CollembolaEast_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_CollembolaEast_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_CollembolaEast_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_CollembolaEast_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_CollembolaEast_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_CollembolaEast_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_CollembolaEast_0.029$beta.sim+0.001->all_0.029_betasim

##**Ddecay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.029_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.029_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.029

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5, GMYC**
##**Plot all levels**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(0.4,1.2), pch=20, lwd=4, cex.lab= 1.5, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
##plot.decay(decay_0.02,add=T,pch=16,lwd=4,col="#d49e57")
plot.decay(decay_0.029,add=T,pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.03,add=T,pch=20,lty=3,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff") 
text(x=0.41, y=0.0, labels="(a)", cex=1.8)

#**Creating plot**
#**Finer scale Collembola WEST Altitude**

read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_West/community_Collembola_h_West.txt")->community_Collembola_h_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_West/community_Collembola0.005_West.txt")->community_Collembola0.005_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_West/community_Collembola0.015_West.txt")->community_Collembola0.015_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_West/community_Collembola0.02_West.txt")->community_Collembola0.02_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_West/community_Collembola0.03_West.txt")->community_Collembola0.03_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_West/community_Collembola0.05_West.txt")->community_Collembola0.05_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_West/community_Collembola0.075_West.txt")->community_Collembola0.075_West
read.table ("../genetic/Data_out/Collembola/Collembola_IBR_3000_West/community_Collembola0.029_West.txt")->community_Collembola0.029_West

##**BETADIVERSITY ORDINATIONS by SITE**
##**Collembola**

##**beta general_Level_Haplotipos**
beta.pair(community_Collembola_h_West, index.family="sorensen")->beta.pair_CollembolaWest_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.005_West, index.family="sorensen")->beta.pair_CollembolaWest_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.015_West, index.family="sorensen")->beta.pair_CollembolaWest_0.015

##**beta general_Level_0.02**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.02_West, index.family="sorensen")->beta.pair_CollembolaWest_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.03_West, index.family="sorensen")->beta.pair_CollembolaWest_0.03

##**beta general_Level_0.05**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.05_West, index.family="sorensen")->beta.pair_CollembolaWest_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.075_West, index.family="sorensen")->beta.pair_CollembolaWest_0.075

##**beta general_Level_0.029**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Collembola0.029_West, index.family="sorensen")->beta.pair_CollembolaWest_0.029

##**FORMAS DE OBTENER LA MATRIX effective resistance** 
Resistance_matrix_West <- read.table("../spatial/IBResitanceElevationMatrix/Result_3000_resist_West.txt", sep = ",", header=T, row.names = 1)
dim(Resistance_matrix_West)
class(Resistance_matrix_West)
Resistance_matrix_West <- as.matrix(Resistance_matrix_West)
class(Resistance_matrix_West)

Resistance_matrix_West[order(row.names(Resistance_matrix_West)),order(colnames(Resistance_matrix_West))]->Resistance_matrix_West #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_West[upper.tri(Resistance_matrix_West)] <- NA 
Resistance_matrix_West
class(Resistance_matrix_West)

Resistance_matrix_West <- as.dist(Resistance_matrix_West) 
Resistance_matrix_West

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_CollembolaWest_h$beta.sim +0.001->all_h_betasim
1-beta.pair_CollembolaWest_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_CollembolaWest_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_CollembolaWest_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_CollembolaWest_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_CollembolaWest_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_CollembolaWest_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_CollembolaWest_0.029$beta.sim+0.001->all_0.029_betasim

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.029_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.029_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.029

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5**
##"H", "NC0.5", "NC1.5", "NC3","NC5", "NC7.5"
##**Plot all levels**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(0.4,1.3), pch=20, lwd=4, cex.lab= 1.5, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
##plot.decay(decay_0.02,add=T,pch=16,lwd=4,col="#d49e57")
plot.decay(decay_0.029,add=T, pch=20,lwd=4, col="#92000A")
plot.decay(decay_0.03,add=T,pch=20,lwd=4,lty=3,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff") 
text(x=0.41, y=0.0, labels="(b)", cex=1.8)
mtext(c("Collembola"), side = 4, adj = 0, col = "black", line = 2, cex=1.8)

#**Creating plot**
##Finer scale Diptera EAST Altitude B

read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_East/community_Diptera_h_East.txt")->community_Diptera_h_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_East/community_Diptera0.005_East.txt")->community_Diptera0.005_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_East/community_Diptera0.015_East.txt")->community_Diptera0.015_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_East/community_Diptera0.02_East.txt")->community_Diptera0.02_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_East/community_Diptera0.03_East.txt")->community_Diptera0.03_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_East/community_Diptera0.05_East.txt")->community_Diptera0.05_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_East/community_Diptera0.075_East.txt")->community_Diptera0.075_East
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_East/community_Diptera0.0088_East.txt")->community_Diptera0.0088_East

##**BETADIVERSITY ORDINATIONS by East**
##**Diptera**

##**beta general_Level_Haplotipos**
beta.pair(community_Diptera_h_East, index.family="sorensen")->beta.pair_DipteraEast_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.005_East, index.family="sorensen")->beta.pair_DipteraEast_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.015_East, index.family="sorensen")->beta.pair_DipteraEast_0.015

##**beta general_Level_0.02**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.02_East, index.family="sorensen")->beta.pair_DipteraEast_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.03_East, index.family="sorensen")->beta.pair_DipteraEast_0.03

##**beta general_Level_0.05**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.05_East, index.family="sorensen")->beta.pair_DipteraEast_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.075_East, index.family="sorensen")->beta.pair_DipteraEast_0.075

##**beta general_Level_0.0088**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.0088_East, index.family="sorensen")->beta.pair_DipteraEast_0.0088

##**FORMAS DE OBTENER LA MATRIX effective resistance** 
Resistance_matrix_East <- read.table("../spatial/IBResitanceElevationMatrix/Result_Alt_B_resist_East.txt", sep = ",", header=T, row.names = 1)
dim(Resistance_matrix_East)
class(Resistance_matrix_East)
Resistance_matrix_East <- as.matrix(Resistance_matrix_East)
class(Resistance_matrix_East)

Resistance_matrix_East[order(row.names(Resistance_matrix_East)),order(colnames(Resistance_matrix_East))]->Resistance_matrix_East #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_East[upper.tri(Resistance_matrix_East)] <- NA 
Resistance_matrix_East
class(Resistance_matrix_East)

Resistance_matrix_East <- as.dist(Resistance_matrix_East) 
Resistance_matrix_East

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_DipteraEast_h$beta.sim +0.001->all_h_betasim
1-beta.pair_DipteraEast_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_DipteraEast_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_DipteraEast_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_DipteraEast_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_DipteraEast_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_DipteraEast_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_DipteraEast_0.0088$beta.sim+0.001->all_0.0088_betasim

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.0088_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")
decay.model(all_0.0088_betasim,Resistance_matrix_East,model.type = "exp",y.type="sim")->decay_0.0088

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5, GMYC**
plot.decay(decay_h, ylim=c(0,1.0),xlim=c(1.2,4.5), pch=20, lwd=4, cex.lab= 1.5, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.0088,add=T,pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
##plot.decay(decay_0.02,add=T,pch=20,lwd=4,col="#d49e57")
plot.decay(decay_0.03,add=T,pch=20,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff")
text(x=1.21, y=0.0, labels="(c)", cex=1.8)
mtext(c("Effective distance (East altitude rasters)"), side = 1, col = "black", line = 2.5, cex=1.2)

#**Creating plot**
##Finer scale Diptera WEST Altitude

read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_West/community_Diptera_h_West.txt")->community_Diptera_h_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_West/community_Diptera0.005_West.txt")->community_Diptera0.005_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_West/community_Diptera0.015_West.txt")->community_Diptera0.015_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_West/community_Diptera0.02_West.txt")->community_Diptera0.02_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_West/community_Diptera0.03_West.txt")->community_Diptera0.03_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_West/community_Diptera0.05_West.txt")->community_Diptera0.05_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_West/community_Diptera0.075_West.txt")->community_Diptera0.075_West
read.table ("../genetic/Data_out/Diptera/Diptera_IBR_AltB_West/community_Diptera0.0088_West.txt")->community_Diptera0.0088_West

##**BETADIVERSITY ORDINATIONS by West**
##**Diptera**

##**beta general_Level_Haplotipos**
beta.pair(community_Diptera_h_West, index.family="sorensen")->beta.pair_DipteraWest_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.005_West, index.family="sorensen")->beta.pair_DipteraWest_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.015_West, index.family="sorensen")->beta.pair_DipteraWest_0.015

##**beta general_Level_0.02**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.02_West, index.family="sorensen")->beta.pair_DipteraWest_0.02

#**beta general_Level_0.03**
##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.03_West, index.family="sorensen")->beta.pair_DipteraWest_0.03

##**beta general_Level_0.05**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.05_West, index.family="sorensen")->beta.pair_DipteraWest_0.05

##**beta general_Level_0.075**
##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.075_West, index.family="sorensen")->beta.pair_DipteraWest_0.075

##**beta general_Level_0.0088**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Diptera0.0088_West, index.family="sorensen")->beta.pair_DipteraWest_0.0088

##**FORMAS DE OBTENER LA MATRIX effective resistance** 
Resistance_matrix_West <- read.table("../spatial/IBResitanceElevationMatrix/Result_Alt_B_resist_West.txt", sep = ",", header=T, row.names = 1)
dim(Resistance_matrix_West)
class(Resistance_matrix_West)
Resistance_matrix_West <- as.matrix(Resistance_matrix_West)
class(Resistance_matrix_West)

Resistance_matrix_West[order(row.names(Resistance_matrix_West)),order(colnames(Resistance_matrix_West))]->Resistance_matrix_West #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_West[upper.tri(Resistance_matrix_West)] <- NA 
Resistance_matrix_West
class(Resistance_matrix_West)

Resistance_matrix_West <- as.dist(Resistance_matrix_West) 
Resistance_matrix_West

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_DipteraWest_h$beta.sim +0.001->all_h_betasim
1-beta.pair_DipteraWest_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_DipteraWest_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_DipteraWest_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_DipteraWest_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_DipteraWest_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_DipteraWest_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_DipteraWest_0.0088$beta.sim+0.001->all_0.0088_betasim

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.0088_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")
decay.model(all_0.0088_betasim,Resistance_matrix_West,model.type = "exp",y.type="sim")->decay_0.0088

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(1.2,4.2), pch=20, lwd=4, cex.lab= 1.5, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.0088,add=T,pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
##plot.decay(decay_0.02,add=T,pch=20,lwd=4,col="#d49e57")
plot.decay(decay_0.03,add=T,pch=20,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff")
text(x=1.21, y=0.0, labels="(d)", cex=1.8)
mtext(c("Diptera"), side = 4, adj = 0, col = "black", line = 2, cex=1.8)
mtext(c("Effective distance (West altitude rasters)"), side = 1, col = "black", line = 2.5, cex=1.2)
par(xpd=TRUE)
###Legend
legend(4.75, 1, xpd=NA,
       legend=c("Haplotype", "CL 0.5", "CL 1.5", "GMYC", "CL 3","CL 5", "CL 7.5"), 
       col=c("#003695", "#c997a9", "#9d9cc6", "#92000A", "#9cb15b", "#fbd048", "#93dfff"), 
       pch=19,  bty="n", text.font=1.7, lty=1, cex=1.7, lwd=2.2) 

###Legend
mtext("Similarity (Simpson's Index)", side=2, outer=TRUE, line=0.7, cex=1.6)

dev.off()
#


#**Plot Beta diversity vs IBD large scale**

png(filename="../figures/FigureS7_DistanceDecay.png", width=787, height=467, units="px") # set size of the file to plot 
par(mfrow=c(2,2), mar = c(1.5, 1.5, 1.5, 1.5), omi=c(0.5, 0.5, 0.5, 2)) #the number of rows and columns the figure would have

#**Creating plot**
##IBD_Arachnida
read.table ("../genetic/Data_out/Arachnida/Arachnida_IBR_Flat_Site/community_Arachnida_h_Site.txt")->community_Arachnida_h_Site
read.table ("../genetic/Data_out/Arachnida/Arachnida_IBR_Flat_Site/community_Arachnida0.005_Site.txt")->community_Arachnida0.005_Site
read.table ("../genetic/Data_out/Arachnida/Arachnida_IBR_Flat_Site/community_Arachnida0.015_Site.txt")->community_Arachnida0.015_Site
read.table ("../genetic/Data_out/Arachnida/Arachnida_IBR_Flat_Site/community_Arachnida0.015_Site.txt")->community_Arachnida0.015_Site
read.table ("../genetic/Data_out/Arachnida/Arachnida_IBR_Flat_Site/community_Arachnida0.02_Site.txt")->community_Arachnida0.02_Site
read.table ("../genetic/Data_out/Arachnida/Arachnida_IBR_Flat_Site/community_Arachnida0.03_Site.txt")->community_Arachnida0.03_Site
read.table ("../genetic/Data_out/Arachnida/Arachnida_IBR_Flat_Site/community_Arachnida0.05_Site.txt")->community_Arachnida0.05_Site
read.table ("../genetic/Data_out/Arachnida/Arachnida_IBR_Flat_Site/community_Arachnida0.075_Site.txt")->community_Arachnida0.075_Site
read.table ("../genetic/Data_out/Arachnida/Arachnida_IBR_Flat_Site/community_Arachnida0.0126_Site.txt")->community_Arachnida0.0126_Site


#**Arachnida**

##**beta general_Level_Haplotipos**
beta.pair(community_Arachnida_h_Site, index.family="sorensen")->beta.pair_ArachnidaSite_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Arachnida0.005_Site, index.family="sorensen")->beta.pair_ArachnidaSite_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Arachnida0.015_Site, index.family="sorensen")->beta.pair_ArachnidaSite_0.015

##**beta general_Level_0.02**
##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Arachnida0.02_Site, index.family="sorensen")->beta.pair_ArachnidaSite_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Arachnida0.03_Site, index.family="sorensen")->beta.pair_ArachnidaSite_0.03

##**beta general_Level_0.05**
##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Arachnida0.05_Site, index.family="sorensen")->beta.pair_ArachnidaSite_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Arachnida0.075_Site, index.family="sorensen")->beta.pair_ArachnidaSite_0.075

##**beta general_Level_0.0126**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Arachnida0.0126_Site, index.family="sorensen")->beta.pair_ArachnidaSite_0.0126
#

##**FORMAS DE OBTENER LA MATRIX effective resistance** 
Resistance_matrix_Site <- read.table("../spatial/IBResistanceFlatMatrix/Result_flat_resistances_ArachnidaSite.txt", sep= ",", header=T, row.names = 1)
dim(Resistance_matrix_Site)
class(Resistance_matrix_Site)
Resistance_matrix_Site <- as.matrix(Resistance_matrix_Site)
class(Resistance_matrix_Site)

Resistance_matrix_Site[order(row.names(Resistance_matrix_Site)),order(colnames(Resistance_matrix_Site))]->Resistance_matrix_Site #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_Site[upper.tri(Resistance_matrix_Site)] <- NA 
Resistance_matrix_Site
class(Resistance_matrix_Site)

Resistance_matrix_Site <- as.dist(Resistance_matrix_Site) 
Resistance_matrix_Site
#

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_ArachnidaSite_h$beta.sim +0.001->all_h_betasim
1-beta.pair_ArachnidaSite_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_ArachnidaSite_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_ArachnidaSite_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_ArachnidaSite_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_ArachnidaSite_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_ArachnidaSite_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_ArachnidaSite_0.0126$beta.sim+0.001->all_0.0126_betasim

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.0126_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.0126_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.0126
#

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5, GMYC**
##**Plot all levels**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(0.4,1.5), pch=20, lwd=4, cex.lab= 1.5, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.0126,add=T,pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
plot.decay(decay_0.02,add=T,pch=20,lwd=4,col="#d49e57")
plot.decay(decay_0.03,add=T,pch=20,lty=3,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff")
mtext(c("Arachnida"), side = 3, adj = 1, col = "black", line = 0, cex=1.2)
#

#**Creating plot**
#**IBD Coleoptera**

read.table ("../genetic/Data_out/Coleoptera/Coleoptera_IBR_Flat_Site/community_Coleoptera_h_Site.txt")->community_Coleoptera_h_Site
read.table ("../genetic/Data_out/Coleoptera/Coleoptera_IBR_Flat_Site/community_Coleoptera0.005_Site.txt")->community_Coleoptera0.005_Site
read.table ("../genetic/Data_out/Coleoptera/Coleoptera_IBR_Flat_Site/community_Coleoptera0.015_Site.txt")->community_Coleoptera0.015_Site
read.table ("../genetic/Data_out/Coleoptera/Coleoptera_IBR_Flat_Site/community_Coleoptera0.015_Site.txt")->community_Coleoptera0.015_Site
read.table ("../genetic/Data_out/Coleoptera/Coleoptera_IBR_Flat_Site/community_Coleoptera0.02_Site.txt")->community_Coleoptera0.02_Site
read.table ("../genetic/Data_out/Coleoptera/Coleoptera_IBR_Flat_Site/community_Coleoptera0.03_Site.txt")->community_Coleoptera0.03_Site
read.table ("../genetic/Data_out/Coleoptera/Coleoptera_IBR_Flat_Site/community_Coleoptera0.05_Site.txt")->community_Coleoptera0.05_Site
read.table ("../genetic/Data_out/Coleoptera/Coleoptera_IBR_Flat_Site/community_Coleoptera0.075_Site.txt")->community_Coleoptera0.075_Site
read.table ("../genetic/Data_out/Coleoptera/Coleoptera_IBR_Flat_Site/community_Coleoptera0.0073_Site.txt")->community_Coleoptera0.0073_Site

#**BETADIVERSITY ORDINATIONS by SITE**
#**Coleoptera**

##**beta general_Level_Haplotipos**
beta.pair(community_Coleoptera_h_Site, index.family="sorensen")->beta.pair_ColeopteraSite_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Coleoptera0.005_Site, index.family="sorensen")->beta.pair_ColeopteraSite_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Coleoptera0.015_Site, index.family="sorensen")->beta.pair_ColeopteraSite_0.015

##**beta general_Level_0.02**
##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Coleoptera0.02_Site, index.family="sorensen")->beta.pair_ColeopteraSite_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Coleoptera0.03_Site, index.family="sorensen")->beta.pair_ColeopteraSite_0.03

##**beta general_Level_0.05**
##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Coleoptera0.05_Site, index.family="sorensen")->beta.pair_ColeopteraSite_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Coleoptera0.075_Site, index.family="sorensen")->beta.pair_ColeopteraSite_0.075

##**beta general_Level_0.0073**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Coleoptera0.0073_Site, index.family="sorensen")->beta.pair_ColeopteraSite_0.0073
#

##**FORMAS DE OBTENER LA MATRIX effective resistance** 
Resistance_matrix_Site <- read.table("../spatial/IBResistanceFlatMatrix/Result_flat_resistances_ColeopteraSite.txt", sep= ",", header=T, row.names = 1)
dim(Resistance_matrix_Site)
class(Resistance_matrix_Site)
Resistance_matrix_Site <- as.matrix(Resistance_matrix_Site)
class(Resistance_matrix_Site)

Resistance_matrix_Site[order(row.names(Resistance_matrix_Site)),order(colnames(Resistance_matrix_Site))]->Resistance_matrix_Site #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_Site[upper.tri(Resistance_matrix_Site)] <- NA 
Resistance_matrix_Site
class(Resistance_matrix_Site)

Resistance_matrix_Site <- as.dist(Resistance_matrix_Site) 
Resistance_matrix_Site
#

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_ColeopteraSite_h$beta.sim +0.001->all_h_betasim
1-beta.pair_ColeopteraSite_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_ColeopteraSite_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_ColeopteraSite_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_ColeopteraSite_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_ColeopteraSite_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_ColeopteraSite_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_ColeopteraSite_0.0073$beta.sim+0.001->all_0.0073_betasim
#

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.0073_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.0073_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.0073
#

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5, GMYC**
##**Plot all levels**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(0.4,1.5), pch=20, lwd=4, cex.lab= 1.5, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.0073,add=T,pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
plot.decay(decay_0.02,add=T,pch=20,lwd=4,col="#d49e57")
plot.decay(decay_0.03,add=T,pch=20,lty=3,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff")
mtext(c("Coleoptera"), side = 3, adj = 1, col = "black", line = 0, cex=1.2)
#

#**Creating plot**
##IBD Hemiptera

read.table ("../genetic/Data_out/Hemiptera/Hemiptera_IBR_Flat_Site/community_Hemiptera_h_Site.txt")->community_Hemiptera_h_Site
read.table ("../genetic/Data_out/Hemiptera/Hemiptera_IBR_Flat_Site/community_Hemiptera0.005_Site.txt")->community_Hemiptera0.005_Site
read.table ("../genetic/Data_out/Hemiptera/Hemiptera_IBR_Flat_Site/community_Hemiptera0.015_Site.txt")->community_Hemiptera0.015_Site
read.table ("../genetic/Data_out/Hemiptera/Hemiptera_IBR_Flat_Site/community_Hemiptera0.015_Site.txt")->community_Hemiptera0.015_Site
read.table ("../genetic/Data_out/Hemiptera/Hemiptera_IBR_Flat_Site/community_Hemiptera0.02_Site.txt")->community_Hemiptera0.02_Site
read.table ("../genetic/Data_out/Hemiptera/Hemiptera_IBR_Flat_Site/community_Hemiptera0.03_Site.txt")->community_Hemiptera0.03_Site
read.table ("../genetic/Data_out/Hemiptera/Hemiptera_IBR_Flat_Site/community_Hemiptera0.05_Site.txt")->community_Hemiptera0.05_Site
read.table ("../genetic/Data_out/Hemiptera/Hemiptera_IBR_Flat_Site/community_Hemiptera0.075_Site.txt")->community_Hemiptera0.075_Site
read.table ("../genetic/Data_out/Hemiptera/Hemiptera_IBR_Flat_Site/community_Hemiptera0.018_Site.txt")->community_Hemiptera0.018_Site

#**BETADIVERSITY ORDINATIONS by SITE**
#**Hemiptera**

##**beta general_Level_Haplotipos**
beta.pair(community_Hemiptera_h_Site, index.family="sorensen")->beta.pair_HemipteraSite_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hemiptera0.005_Site, index.family="sorensen")->beta.pair_HemipteraSite_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hemiptera0.015_Site, index.family="sorensen")->beta.pair_HemipteraSite_0.015

##**beta general_Level_0.02**
##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hemiptera0.02_Site, index.family="sorensen")->beta.pair_HemipteraSite_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hemiptera0.03_Site, index.family="sorensen")->beta.pair_HemipteraSite_0.03

##**beta general_Level_0.05**
##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hemiptera0.05_Site, index.family="sorensen")->beta.pair_HemipteraSite_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hemiptera0.075_Site, index.family="sorensen")->beta.pair_HemipteraSite_0.075

##**beta general_Level_0.018**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hemiptera0.018_Site, index.family="sorensen")->beta.pair_HemipteraSite_0.018
#

##**FORMAS DE OBTENER LA MATRIX effective resistance** 
Resistance_matrix_Site <- read.table("../spatial/IBResistanceFlatMatrix/Result_flat_resistances_HemipteraSite.txt", sep= ",", header=T, row.names = 1)
dim(Resistance_matrix_Site)
class(Resistance_matrix_Site)
Resistance_matrix_Site <- as.matrix(Resistance_matrix_Site)
class(Resistance_matrix_Site)

Resistance_matrix_Site[order(row.names(Resistance_matrix_Site)),order(colnames(Resistance_matrix_Site))]->Resistance_matrix_Site #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_Site[upper.tri(Resistance_matrix_Site)] <- NA 
Resistance_matrix_Site
class(Resistance_matrix_Site)

Resistance_matrix_Site <- as.dist(Resistance_matrix_Site) 
Resistance_matrix_Site
#

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_HemipteraSite_h$beta.sim +0.001->all_h_betasim
1-beta.pair_HemipteraSite_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_HemipteraSite_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_HemipteraSite_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_HemipteraSite_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_HemipteraSite_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_HemipteraSite_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_HemipteraSite_0.018$beta.sim+0.001->all_0.018_betasim
#

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.018_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.018_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.018
#

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5, GMYC**
##**Plot all levels**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(0.3,1.5), pch=20, lwd=4, cex.lab= 1.5, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
plot.decay(decay_0.018,add=T,pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.02,add=T,pch=20,lwd=4,col="#d49e57")
plot.decay(decay_0.03,add=T,pch=20,lty=3,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff")
mtext(c("Hemiptera"), side = 3, adj = 1, col = "black", line = 0, cex=1.2)
#

#**Creating plot**
##IBD Hymenoptera

read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera_IBR_Flat_Site/community_Hymenoptera_h_Site.txt")->community_Hymenoptera_h_Site
read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera_IBR_Flat_Site/community_Hymenoptera0.005_Site.txt")->community_Hymenoptera0.005_Site
read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera_IBR_Flat_Site/community_Hymenoptera0.015_Site.txt")->community_Hymenoptera0.015_Site
read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera_IBR_Flat_Site/community_Hymenoptera0.015_Site.txt")->community_Hymenoptera0.015_Site
read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera_IBR_Flat_Site/community_Hymenoptera0.02_Site.txt")->community_Hymenoptera0.02_Site
read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera_IBR_Flat_Site/community_Hymenoptera0.03_Site.txt")->community_Hymenoptera0.03_Site
read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera_IBR_Flat_Site/community_Hymenoptera0.05_Site.txt")->community_Hymenoptera0.05_Site
read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera_IBR_Flat_Site/community_Hymenoptera0.075_Site.txt")->community_Hymenoptera0.075_Site
read.table ("../genetic/Data_out/Hymenoptera/Hymenoptera_IBR_Flat_Site/community_Hymenoptera0.0098_Site.txt")->community_Hymenoptera0.0098_Site

#**BETADIVERSITY ORDINATIONS by SITE**
#**Hymenoptera**

##**beta general_Level_Haplotipos**
beta.pair(community_Hymenoptera_h_Site, index.family="sorensen")->beta.pair_HymenopteraSite_h

##**beta general_Level_0.005**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hymenoptera0.005_Site, index.family="sorensen")->beta.pair_HymenopteraSite_0.005

##**beta general_Level_0.015**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hymenoptera0.015_Site, index.family="sorensen")->beta.pair_HymenopteraSite_0.015

##**beta general_Level_0.02**
##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hymenoptera0.02_Site, index.family="sorensen")->beta.pair_HymenopteraSite_0.02

##**beta general_Level_0.03**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hymenoptera0.03_Site, index.family="sorensen")->beta.pair_HymenopteraSite_0.03

##**beta general_Level_0.05**
##betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hymenoptera0.05_Site, index.family="sorensen")->beta.pair_HymenopteraSite_0.05

##**beta general_Level_0.075**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hymenoptera0.075_Site, index.family="sorensen")->beta.pair_HymenopteraSite_0.075

##**beta general_Level_0.0098**
#####betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
beta.pair(community_Hymenoptera0.0098_Site, index.family="sorensen")->beta.pair_HymenopteraSite_0.0098
#

##**FORMAS DE OBTENER LA MATRIX effective resistance** 
Resistance_matrix_Site <- read.table("../spatial/IBResistanceFlatMatrix/Result_flat_resistances_HymenopteraSite.txt", sep= ",", header=T, row.names = 1)
dim(Resistance_matrix_Site)
class(Resistance_matrix_Site)
Resistance_matrix_Site <- as.matrix(Resistance_matrix_Site)
class(Resistance_matrix_Site)

Resistance_matrix_Site[order(row.names(Resistance_matrix_Site)),order(colnames(Resistance_matrix_Site))]->Resistance_matrix_Site #Ordena la Resistance matrix con la de beta. #'important order both matrixes

Resistance_matrix_Site[upper.tri(Resistance_matrix_Site)] <- NA 
Resistance_matrix_Site
class(Resistance_matrix_Site)

Resistance_matrix_Site <- as.dist(Resistance_matrix_Site) 
Resistance_matrix_Site
#

##**Generating similarity values and adding 0.001 to avoid LOG(0)**
1-beta.pair_HymenopteraSite_h$beta.sim +0.001->all_h_betasim
1-beta.pair_HymenopteraSite_0.005$beta.sim+0.001->all_0.005_betasim
1-beta.pair_HymenopteraSite_0.015$beta.sim+0.001->all_0.015_betasim
1-beta.pair_HymenopteraSite_0.02$beta.sim+0.001->all_0.02_betasim
1-beta.pair_HymenopteraSite_0.03$beta.sim+0.001->all_0.03_betasim
1-beta.pair_HymenopteraSite_0.05$beta.sim+0.001->all_0.05_betasim
1-beta.pair_HymenopteraSite_0.075$beta.sim+0.001->all_0.075_betasim
1-beta.pair_HymenopteraSite_0.0098$beta.sim+0.001->all_0.0098_betasim
#

##**Decay using geomatrix**
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_h_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_h

decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.005_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.005

decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.015_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.015

decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.02_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.02

decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.03_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.03

decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.05_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.05

decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.075_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.075

decay.model(all_0.0098_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")
decay.model(all_0.0098_betasim,Resistance_matrix_Site,model.type = "exp",y.type="sim")->decay_0.0098
#

##**Plot with levels: h, 0.5, 1.5, 3, 5, 7.5, GMYC**
##**Plot all levels**
plot.decay(decay_h, ylim=c(0,1.0), xlim=c(0.4,1.5), pch=20, lwd=4, cex.lab= 1.5, cex.axis= 1.5, col="#003695")
plot.decay(decay_0.005,add=T,pch=20,lwd=4,col="#c997a9")
plot.decay(decay_0.0098,add=T,pch=20,lwd=4,col="#92000A")
plot.decay(decay_0.015,add=T,pch=20,lwd=4,col="#9d9cc6")
plot.decay(decay_0.02,add=T,pch=20,lwd=4,col="#d49e57")
plot.decay(decay_0.03,add=T,pch=20,lty=3,lwd=4,col="#9cb15b")
plot.decay(decay_0.05,add=T,pch=20,lwd=4,col="#fbd048")
plot.decay(decay_0.075,add=T,pch=20,lwd=4,col="#93dfff")
#text(x=1.21, y=0.0, labels="(d)", cex=1.8)
mtext(c("Hymenoptera"), side = 3, adj = 1, col = "black", line = 0, cex=1.2)
#mtext(c("Effective distance (West altitude rasters)"), side = 1, col = "black", line = 2.5, cex=1.2)
par(xpd=TRUE)
###Legend
legend(1.6, 1, xpd=NA,
       legend=c("Haplotype", "CL 0.5", "CL 1.5", "GMYC", "CL 3","CL 5", "CL 7.5"), 
       col=c("#003695", "#c997a9", "#9d9cc6", "#92000A", "#9cb15b", "#fbd048", "#93dfff"), 
       pch=19,  bty="n", text.font=1.7, lty=1, cex=1.7, lwd=2.2) 

###Legend
mtext("Similarity (Simpson's Index)", side=2, outer=TRUE, line=0.8, cex=1.6)
mtext(("Effective distance flat"), side=1, outer=TRUE, line=0.9, cex=1.6)

dev.off()



#**END**
