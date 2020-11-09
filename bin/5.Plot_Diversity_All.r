
library(raster)
require(gridExtra)
library(sp)
library(ggplot2)
library(patchwork)
library(easypackages)


#'**Community diversity and composition at the haplotype, 3% and 5% Clustering levels for each of the eight taxonomic orders studied**'#

#This script gets all scripts of diversity using 8 arthropods order at multi-hierarchical levels 
#We used 8 groups: Arachnda, Coleoptera, Collembola, Diptera, Hemiptera, Hymenoptera, Lepidoptera, and Myriapoda. 

################ plots Community diversity and composition #####################


png(filename="../figures/RelativeSpeciesRichness1.png", width=920 , height=1536, units="px") # set size of the file to plot 
par(mfrow=c(8,3), xpd = TRUE) #the number of rows and columns the figure would have


#mar = c(2.5, 3, 3, 3), xpd = TRUE),
#'**Creating plot'**#

#'**Community diversity and composition at the haplotype, 3% and 5% Clustering levels for each of the eight taxonomic orders studied**'#

#This script gets all scripts of diversity using 8 arthropods order at multi-hierarchical levels 
#We used 8 groups: Arachnda, Coleoptera, Collembola, Diptera, Hemiptera, Hymenoptera, Lepidoptera, and Myriapoda. 


################ plot Global Richness by Sites #####################
#'**Diptera'**#
read.table("../genetic/Data_out/Diptera/Diptera_Haplotypes/richness_Site_Diptera.txt",header=TRUE)->richness_Site_h_Diptera 

plot(richness_Site_h_Diptera$Mountain1Site,richness_Site_h_Diptera$sample_richness_Site_h, xaxt="n", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,154))
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h_Diptera)
posthoc.kruskal.nemenyi.test(x=richness_Site_h_Diptera$sample_richness_Site_h, g=richness_Site_h_Diptera$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4.1), y=(152), labels=c("a","a","ab","b"), cex=1.5)
text(x=4.5, y=152, labels="*", cex=2)
mtext(c("Haplotypes"), side = 3, col = "black", line = 1, cex = 1.5)

read.table("../genetic/Data_out/Diptera/Diptera3P/richness_Site0.03_Diptera.txt",header=TRUE)->richness_Site0.03_Diptera

plot(richness_Site0.03_Diptera$Mountain1Site,richness_Site0.03_Diptera$sample_richness_Site0.03, xaxt="n", ylim=c(0,50), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03_Diptera)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03_Diptera$sample_richness_Site0.03, g=richness_Site0.03_Diptera$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4.1), y=(48), labels=c("ab","a","ab","b"), cex=1.5)
text(x=4.5, y=48, labels="*", cex=2)
mtext(c("Clustering 3%"), side = 3, col = "black", line = 1, cex = 1.5)

read.table("../genetic/Data_out/Diptera/Diptera5P/richness_Site0.05_Diptera.txt",header=TRUE)->richness_Site0.05_Diptera

plot(richness_Site0.05_Diptera$Mountain1Site,richness_Site0.05_Diptera$sample_richness_Site0.05, xaxt="n", ylim=c(0,50), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05_Diptera)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05_Diptera$sample_richness_Site0.05, g=richness_Site0.05_Diptera$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4.1), y=(48), labels=c("ab","a","ab","b"), cex=1.5)
text(x=4.5, y=48, labels="*", cex=2)
mtext(c("Clustering 5%"), side = 3, col = "black", line = 1, cex = 1.5)
mtext(c("Diptera"), side = 4, col = "black", line = 1, cex = 1.5)


#'**Collembola'**#

read.table("../genetic/Data_out/Collembola/Collembola_Haplotypes/richness_Site_h_Collembola.txt",header=TRUE)->richness_Site_h_Collembola

plot(richness_Site_h_Collembola$Mountain1Site,richness_Site_h_Collembola$sample_richness_Site_h, xaxt="n", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,60))
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h_Collembola)
posthoc.kruskal.nemenyi.test(x=richness_Site_h_Collembola$sample_richness_Site_h, g=richness_Site_h_Collembola$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4), y=(15.7), labels=c("a","b","ab","b"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=4.5, y=59, labels="ns", cex=2)

read.table("../genetic/Data_out/Collembola/Collembola3P/richness_Site_Collembola0.03.txt",header=TRUE)->richness_Site0.03_Collembola

plot(richness_Site0.03_Collembola$Mountain1Site,richness_Site0.03_Collembola$sample_richness_Site0.03, xaxt="n", ylim=c(0,25), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03_Collembola)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03_Collembola$sample_richness_Site0.03, g=richness_Site0.03_Collembola$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(24), labels=c("a","a","a","b"), cex=1.5)
text(x=4.5, y=24, labels="*", cex=2)

read.table("../genetic/Data_out/Collembola/Collembola5P/richness_Site_Collembola0.05.txt",header=TRUE)->richness_Site0.05_Collembola

plot(richness_Site0.05_Collembola$Mountain1Site,richness_Site0.05_Collembola$sample_richness_Site0.05, xaxt="n", ylim=c(0,25), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05_Collembola)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05_Collembola$sample_richness_Site0.05, g=richness_Site0.05_Collembola$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4), y=(24), labels=c("a","a","a","b"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=4.5, y=24, labels="ns", cex=2)
mtext(c("Collembola"), side = 4, col = "black", line = 1, cex = 1.5)

#'**Arachnida'**#

read.table("../genetic/Data_out/Arachnida/Arachnida_Haplotypes/richness_Site_h_Arachnida_h.txt",header=TRUE)->richness_Site_h_Arachnida

plot(richness_Site_h_Arachnida$Mountain1Site,richness_Site_h_Arachnida$sample_richness_Site_h, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h_Arachnida)
posthoc.kruskal.nemenyi.test(x=richness_Site_h_Arachnida$sample_richness_Site_h, g=richness_Site_h_Arachnida$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(29), labels=c("a","a","b","b"), cex=1.5)
text(x=4.5, y=29, labels="***", cex=2)

read.table("../genetic/Data_out/Arachnida/Arachnida3P/richness_Site_Arachnida_0.03.txt",header=TRUE)->richness_Site0.03_Arachnida

plot(richness_Site0.03_Arachnida$Mountain1Site,richness_Site0.03_Arachnida$sample_richness_Site0.03, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03_Arachnida)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03_Arachnida$sample_richness_Site0.03, g=richness_Site0.03_Arachnida$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(29), labels=c("a","a","b","b"), cex=1.5)
text(x=4.5, y=29, labels="***", cex=2)

read.table("../genetic/Data_out/Arachnida/Arachnida5P/richness_Site_Arachnida0.05.txt",header=TRUE)->richness_Site0.05_Arachnida

plot(richness_Site0.05_Arachnida$Mountain1Site,richness_Site0.05_Arachnida$sample_richness_Site0.05, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05_Arachnida)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05_Arachnida$sample_richness_Site0.05, g=richness_Site0.05_Arachnida$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(29), labels=c("a","a","b","b"), cex=1.5)
text(x=4.5, y=29, labels="***", cex=2)
mtext(c("Arachnida"), side = 4, col = "black", line = 1, cex = 1.5)


#'**Hempitera'**#

read.table("../genetic/Data_out/Hemiptera/Hemiptera_Haplotypes/richness_Site_h_Hemiptera.txt",header=TRUE)->richness_Site_h

plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h, xaxt="n", ylim=c(0,15), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(153), labels=c("a","a","ab","b"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=4.5, y=14, labels="ns", cex=2)

read.table("../genetic/Data_out/Hemiptera/Hemiptera3P/richness_Site0.03_Hemiptera.txt",header=TRUE)->richness_Site0.03

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03, xaxt="n", ylim=c(0,15), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(153), labels=c("a","a","ab","b"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=4.5, y=14, labels="ns", cex=2)

read.table("../genetic/Data_out/Hemiptera/Hemiptera5P/richness_Site0.05_Hemiptera.txt",header=TRUE)->richness_Site0.05

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05, xaxt="n", ylim=c(0,15), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(153), labels=c("a","a","ab","b"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=4.5, y=14, labels="ns", cex=2)
mtext(c("Hemiptera"), side = 4, col = "black", line = 1, cex = 1.5)


#'**Hymenoptera'**#

read.table("../genetic/Data_out/Hymenoptera/Hymenoptera_Haplotypes/richness_Site_h_Hymenoptera.txt",header=TRUE)->richness_Site_h

plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.5)
text(x=4.5, y=29, labels="**", cex=2)

read.table("../genetic/Data_out/Hymenoptera/Hymenoptera3P/richness_Site0.03_Hymenoptera.txt",header=TRUE)->richness_Site0.03

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.5)
text(x=4.5, y=29, labels="*", cex=2)

read.table("../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_Site0.05_Hymenoptera.txt",header=TRUE)->richness_Site0.05

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=4.5, y=29, labels="ns", cex=2)
mtext(c("Hymenoptera"), side = 4, col = "black", line = 1, cex = 1.5)


#'**Coleoptera'**#

read.table("../genetic/Data_out/Coleoptera/Coleoptera_Haplotypes/richness_Site_h_Coleoptera.txt",header=TRUE)->richness_Site_h

plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h, xaxt="n", ylim=c(0,16), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(15), labels=c("a","b","ab","b"), cex=1.5)
text(x=4.5, y=15, labels="**", cex=2)

read.table("../genetic/Data_out/Coleoptera/Coleoptera3P/richness_Site_Coleoptera0.03.txt",header=TRUE)->richness_Site0.03

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03, xaxt="n", ylim=c(0,16), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(15), labels=c("a","b","ab","b"), cex=1.5)
text(x=4.5, y=15, labels="**", cex=2)

read.table("../genetic/Data_out/Coleoptera/Coleoptera5P/richness_Site0.05_Coleoptera.txt",header=TRUE)->richness_Site0.05

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05, xaxt="n", ylim=c(0,16), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(15), labels=c("a","b","ab","b"), cex=1.5)
text(x=4.5, y=15, labels="**", cex=2)
mtext(c("Coleoptera"), side = 4, col = "black", line = 1, cex = 1.5)

#'**Myriapoda'**#

read.table("../genetic/Data_out/Myriapoda/Myriapoda_Haplotypes/richness_Site_h_Myriapoda.txt",header=TRUE)->richness_Site_h

plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h, xaxt="n", ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=4.5, y=7, labels="ns", cex=2)

read.table("../genetic/Data_out/Myriapoda/Myriapoda3P/richness_Site0.03_Myriapoda.txt",header=TRUE)->richness_Site0.03

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03, xaxt="n", ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=4.5, y=7, labels="ns", cex=2)

read.table("../genetic/Data_out/Myriapoda/Myriapoda5P/richness_Site0.05_Myriapoda.txt",header=TRUE)->richness_Site0.05

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05, xaxt="n", ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=4.5, y=7, labels="ns", cex=2)
mtext(c("Myriapoda"), side = 4, col = "black", line = 1, cex = 1.5)

#'**Lepidoptera'**#

read.table("../genetic/Data_out/Lepidoptera/Lepidoptera_Haplotypes/richness_Site_h_Lepidoptera.txt",header=TRUE)->richness_Site_h

plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h, xaxt="n", ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), -1.3, labels=c("AAB","ASB","SJH","TLC"), cex=2, xpd=NA)
text(x=4.5, y=7, labels="ns", cex=2)


read.table("../genetic/Data_out/Lepidoptera/Lepidoptera3P/richness_Site_Lepidoptera.txt",header=TRUE)->richness_Site0.03

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03, xaxt="n", ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), -1.3, labels=c("AAB","ASB","SJH","TLC"), cex=2, xpd=NA)
text(x=4.5, y=7, labels="ns", cex=2)

read.table("../genetic/Data_out/Lepidoptera/Lepidoptera5P/richness_Site_Lepidoptera0.05.txt",header=TRUE)->richness_Site0.05

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05, xaxt="n", tick = TRUE, ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), -1.3, labels=c("AAB","ASB","SJH","TLC"), cex=2, xpd=NA)
text(x=4.5, y=7, labels="ns", cex=2)
mtext(c("Lepidoptera"), side = 4, col = "black", line = 1, cex = 1.5)

mtext("S a m p l e d   s i t e s", side=1, outer=TRUE, line=-1, cex=2)
mtext("R i c h n e s s", side=2, outer=TRUE, line=-2, cex=2) 

dev.off()

####################################END######################################################################



