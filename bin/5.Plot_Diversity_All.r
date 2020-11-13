
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


png(filename="../figures/Figure3.RelativeSpeciesRichness.png", width=920 , height=1536, units="px") # set size of the file to plot 
par(mfrow=c(8,3), mar = c(3, 2.5, 2, 2), omi=c(0.5, 0.5, 0.5, 0.5)) #the number of rows and columns the figure would have

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
text(x=c(1,2,3,4.1), y=(152), labels=c("a","a","ab","b"), cex=2)
text(x=4.5, y=152, labels="*", cex=2)
mtext(c("Haplotypes"), side = 3, col = "black", line = 1, cex=2)

read.table("../genetic/Data_out/Diptera/Diptera3P/richness_Site0.03_Diptera.txt",header=TRUE)->richness_Site0.03_Diptera

plot(richness_Site0.03_Diptera$Mountain1Site,richness_Site0.03_Diptera$sample_richness_Site0.03, xaxt="n", ylim=c(0,50), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03_Diptera)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03_Diptera$sample_richness_Site0.03, g=richness_Site0.03_Diptera$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4.1), y=(48), labels=c("ab","a","ab","b"), cex=2)
text(x=4.5, y=48, labels="*", cex=2)
mtext(c("Clustering 3%"), side = 3, col = "black", line = 1, cex=2)

read.table("../genetic/Data_out/Diptera/Diptera5P/richness_Site0.05_Diptera.txt",header=TRUE)->richness_Site0.05_Diptera

plot(richness_Site0.05_Diptera$Mountain1Site,richness_Site0.05_Diptera$sample_richness_Site0.05, xaxt="n", ylim=c(0,50), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05_Diptera)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05_Diptera$sample_richness_Site0.05, g=richness_Site0.05_Diptera$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4.1), y=(48), labels=c("ab","a","ab","b"), cex=2)
text(x=4.5, y=48, labels="*", cex=2)
mtext(c("Clustering 5%"), side = 3, col = "black", line = 1, cex=2)
mtext(c("Diptera"), side = 4, col = "black", line = 1, cex=2)


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
text(x=c(1,2,3,4), y=(24), labels=c("a","a","a","b"), cex=2)
text(x=4.5, y=24, labels="*", cex=2)

read.table("../genetic/Data_out/Collembola/Collembola5P/richness_Site_Collembola0.05.txt",header=TRUE)->richness_Site0.05_Collembola

plot(richness_Site0.05_Collembola$Mountain1Site,richness_Site0.05_Collembola$sample_richness_Site0.05, xaxt="n", ylim=c(0,25), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05_Collembola)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05_Collembola$sample_richness_Site0.05, g=richness_Site0.05_Collembola$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4), y=(24), labels=c("a","a","a","b"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=4.5, y=24, labels="ns", cex=2)
mtext(c("Collembola"), side = 4, col = "black", line = 1, cex=2)

#'**Arachnida'**#

read.table("../genetic/Data_out/Arachnida/Arachnida_Haplotypes/richness_Site_h_Arachnida_h.txt",header=TRUE)->richness_Site_h_Arachnida

plot(richness_Site_h_Arachnida$Mountain1Site,richness_Site_h_Arachnida$sample_richness_Site_h, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h_Arachnida)
posthoc.kruskal.nemenyi.test(x=richness_Site_h_Arachnida$sample_richness_Site_h, g=richness_Site_h_Arachnida$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(29), labels=c("a","a","b","b"), cex=2)
text(x=4.5, y=29, labels="***", cex=2)

read.table("../genetic/Data_out/Arachnida/Arachnida3P/richness_Site_Arachnida_0.03.txt",header=TRUE)->richness_Site0.03_Arachnida

plot(richness_Site0.03_Arachnida$Mountain1Site,richness_Site0.03_Arachnida$sample_richness_Site0.03, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03_Arachnida)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03_Arachnida$sample_richness_Site0.03, g=richness_Site0.03_Arachnida$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(29), labels=c("a","a","b","b"), cex=2)
text(x=4.5, y=29, labels="***", cex=2)

read.table("../genetic/Data_out/Arachnida/Arachnida5P/richness_Site_Arachnida0.05.txt",header=TRUE)->richness_Site0.05_Arachnida

plot(richness_Site0.05_Arachnida$Mountain1Site,richness_Site0.05_Arachnida$sample_richness_Site0.05, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05_Arachnida)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05_Arachnida$sample_richness_Site0.05, g=richness_Site0.05_Arachnida$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(29), labels=c("a","a","b","b"), cex=2)
text(x=4.5, y=29, labels="***", cex=2)
mtext(c("Arachnida"), side = 4, col = "black", line = 1, cex=2)


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
mtext(c("Hemiptera"), side = 4, col = "black", line = 1, cex=2)


#'**Hymenoptera'**#

read.table("../genetic/Data_out/Hymenoptera/Hymenoptera_Haplotypes/richness_Site_h_Hymenoptera.txt",header=TRUE)->richness_Site_h

plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=2)
text(x=4.5, y=29, labels="**", cex=2)

read.table("../genetic/Data_out/Hymenoptera/Hymenoptera3P/richness_Site0.03_Hymenoptera.txt",header=TRUE)->richness_Site0.03

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=2)
text(x=4.5, y=29, labels="*", cex=2)

read.table("../genetic/Data_out/Hymenoptera/Hymenoptera5P/richness_Site0.05_Hymenoptera.txt",header=TRUE)->richness_Site0.05

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05, xaxt="n", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=4.5, y=29, labels="ns", cex=2)
mtext(c("Hymenoptera"), side = 4, col = "black", line = 1, cex=2)


#'**Coleoptera'**#

read.table("../genetic/Data_out/Coleoptera/Coleoptera_Haplotypes/richness_Site_h_Coleoptera.txt",header=TRUE)->richness_Site_h

plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h, xaxt="n", ylim=c(0,16), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(15), labels=c("a","b","ab","b"), cex=2)
text(x=4.5, y=15, labels="**", cex=2)

read.table("../genetic/Data_out/Coleoptera/Coleoptera3P/richness_Site_Coleoptera0.03.txt",header=TRUE)->richness_Site0.03

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03, xaxt="n", ylim=c(0,16), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against.
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(15), labels=c("a","b","ab","b"), cex=2)
text(x=4.5, y=15, labels="**", cex=2)

read.table("../genetic/Data_out/Coleoptera/Coleoptera5P/richness_Site0.05_Coleoptera.txt",header=TRUE)->richness_Site0.05

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05, xaxt="n", ylim=c(0,16), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
axis(1, at=NULL, label=F, tick=TRUE, lwd=1.3)
text(x=c(1,2,3,4), y=(15), labels=c("a","b","ab","b"), cex=2)
text(x=4.5, y=15, labels="**", cex=2)
mtext(c("Coleoptera"), side = 4, col = "black", line = 1, cex=2)

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
mtext(c("Myriapoda"), side = 4, col = "black", line = 1, cex=2)

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
mtext(c("Lepidoptera"), side = 4, col = "black", line = 1, cex=2)

mtext("S a m p l e d   s i t e s", side=1, outer=TRUE, line=1.5, cex=2.5)
mtext("R i c h n e s s", side=2, outer=TRUE, line=0.7, cex=2.5) 


dev.off()


################ Non-Metric Multidimensional scaling (NMDS) ordinations of community similarity #####################


png(filename="../figures/Figure4.Non-MetricMultidimensionalScaling.png", width=920 , height=1536, units="px") # set size of the file to plot 
par(mfrow=c(8,3), mar = c(2, 2, 2, 2), omi=c(0.5, 0.5, 0.5, 0.5)) #the number of rows and columns the figure would have
#omi=c(0.5, 0.5, 0.5, 0.5)
#mar = c(2, 2, 2, 2)
#'**Creating plot'**#
#'**Diptera'**#
#'
read.table("../genetic/Data_out/Diptera/Diptera_Haplotypes/community_Diptera_h.txt")->community_Diptera_h
read.table("../genetic/Data_out/Diptera/Diptera_Haplotypes/general_sample_Mountain1Site_h.txt",header=TRUE)->general_sample_Mountain1Site_h

beta.multi(community_Diptera_Site_h, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Diptera_Site_h, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_h #'NMDS
plot (MDSbetasim_h, xlim=c(-0.4, 0.4), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain1Site_h,ordispider(MDSbetasim_h, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.595 ***", cex=2)
mtext(c("Haplotypes"), side = 3, col = "black", line = 1, cex = 2)

read.table("../genetic/Data_out/Diptera/Diptera3P/community_Diptera0.03.txt")->community_Diptera0.03
read.table("../genetic/Data_out/Diptera/Diptera3P/general_sample_Mountain1Site0.03.txt",header=TRUE)->general_sample_Mountain1Site0.03

beta.multi(community_Diptera_Site0.03, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Diptera_Site0.03, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.03 #'NMDS
plot (MDSbetasim0.03, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain1Site0.03,ordispider(MDSbetasim0.03, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.372 ***", cex=2)
mtext(c("Clustering 3%"), side = 3, col = "black", line = 1, cex = 2)

read.table("../genetic/Data_out/Diptera/Diptera5P/community_Diptera0.05.txt")->community_Diptera0.05
read.table("../genetic/Data_out/Diptera/Diptera5P/general_sample_Mountain1Site0.05.txt",header=TRUE)->general_sample_Mountain1Site0.05

beta.multi(community_Diptera_Site0.05, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Diptera_Site0.05, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.05 #'NMDS
plot (MDSbetasim0.05, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain1Site0.05,ordispider(MDSbetasim0.05, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.369 ***", cex=2)
mtext(c("Clustering 5%"), side = 3, col = "black", line = 1, cex = 2)
mtext(c("Diptera"), side = 4, col = "black", line = 1, cex = 2)

#'**Collembola'**#
read.table("../genetic/Data_out/Collembola/Collembola_Haplotypes/community_Collembola_h.txt")->community_Collembola_h
read.table("../genetic/Data_out/Collembola/Collembola_Haplotypes/general_sample_Mountain1Site_h.txt",header=TRUE)->general_sample_Mountain1Site_h

beta.multi(community_Collembola_Site_h, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Collembola_Site_h, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_h #'NMDS
plot (MDSbetasim_h, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain1Site_h,ordispider(MDSbetasim_h, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.914 ***", cex=2)

read.table("../genetic/Data_out/Collembola/Collembola3P/community_Collembola0.03.txt")->community_Collembola0.03
read.table("../genetic/Data_out/Collembola/Collembola3P/general_sample_Mountain1Site0.03.txt",header=TRUE)->general_sample_Mountain1Site0.03

beta.multi(community_Collembola_Site0.03, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Collembola_Site0.03, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.03 #'NMDS
plot (MDSbetasim0.03, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain1Site0.03,ordispider(MDSbetasim0.03, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.826 ***", cex=2)

read.table("../genetic/Data_out/Collembola/Collembola5P/community_Collembola0.05.txt")->community_Collembola0.05
read.table("../genetic/Data_out/Collembola/Collembola5P/general_sample_Mountain1Site0.05.txt",header=TRUE)->general_sample_Mountain1Site0.05

beta.multi(community_Collembola_Site0.05, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Collembola_Site0.05, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.05 #'NMDS
plot (MDSbetasim0.05, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain1Site0.05,ordispider(MDSbetasim0.05, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.817 ***", cex=2)
mtext(c("Collembola"), side = 4, col = "black", line = 1, cex = 2)


#'**Arachnida'**#

read.table("../genetic/Data_out/Arachnida/Arachnida_Haplotypes/community_Arachnida_h.txt")->community_Arachnida_h
read.table("../genetic/Data_out/Arachnida/Arachnida_Haplotypes/general_sample_Mountain1Site_h.txt",header=TRUE)->general_sample_Mountain1Site_h

beta.multi(community_Arachnida_Site_h, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Arachnida_Site_h, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_h #'NMDS
plot (MDSbetasim_h, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain1Site_h,ordispider(MDSbetasim_h, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.30, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.569 ***", cex=2)

read.table("../genetic/Data_out/Arachnida/Arachnida3P/community_Arachnida0.03.txt")->community_Arachnida0.03
read.table("../genetic/Data_out/Arachnida/Arachnida3P/general_sample_Mountain1Site0.03.txt",header=TRUE)->general_sample_Mountain1Site0.03

beta.multi(community_Arachnida_Site0.03, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Arachnida_Site0.03, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.03 #'NMDS
plot (MDSbetasim0.03, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain1Site0.03,ordispider(MDSbetasim0.03, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.30, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.452 ***", cex=2)

read.table("../genetic/Data_out/Arachnida/Arachnida5P/community_Arachnida0.05.txt")->community_Arachnida0.05
read.table("../genetic/Data_out/Arachnida/Arachnida5P/general_sample_Mountain1Site0.05.txt",header=TRUE)->general_sample_Mountain1Site0.05

beta.multi(community_Arachnida_Site0.05, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Arachnida_Site0.05, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.05 #'NMDS
plot (MDSbetasim0.05, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain1Site0.05,ordispider(MDSbetasim0.05, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.30, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.461 ***", cex=2)
mtext(c("Arachnida"), side = 4, col = "black", line = 1, cex = 2)

#'**Hempitera'**#
read.table("../genetic/Data_out/Hemiptera/Hemiptera_Haplotypes/community_Hemiptera_h.txt")->community_Hemiptera_h
read.table("../genetic/Data_out/Hemiptera/Hemiptera_Haplotypes/general_sample_Mountain1Site_h.txt",header=TRUE)->general_sample_Mountain1Site_h

beta.multi(community_Hemiptera_Site_h, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Hemiptera_Site_h, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
#'repeating after removing outlayers from matrix and general habitat table
community_Hemiptera_Site_h[-which(row.names(community_Hemiptera_Site_h) %in% c( "CON_NTO_TLC_36TCONS6", "CON_NTO_TLC_31TCONS1", "CON_NTO_SJH_76SHCON7")),]->community_Hemiptera_sinoutlayer_h
community_Hemiptera_sinoutlayer_h[,which(colSums(community_Hemiptera_sinoutlayer_h)!=0)]->community_Hemiptera_sinoutlayer_h #'to remove no data colums
dim(community_Hemiptera_sinoutlayer_h)  

general_sample_Mountain1Site_h[-which(general_sample_Mountain1Site_h$sample_names %in% c("CON_NTO_TLC_36TCONS6", "CON_NTO_TLC_31TCONS1", "CON_NTO_SJH_76SHCON7")),]->general_sample_sinoutlayer_h
beta.pair(community_Hemiptera_sinoutlayer_h, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_h
#with(general_sample_sinoutlayer_h,ordispider(MDSbetasim_h, Mountain, label=T, col="green"))
plot (MDSbetasim_h, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.49), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer_h,ordispider(MDSbetasim_h, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.418 ***", cex=2)


read.table("../genetic/Data_out/Hemiptera/Hemiptera3P/community_Hemiptera0.03.txt")->community_Hemiptera0.03
read.table("../genetic/Data_out/Hemiptera/Hemiptera3P/general_sample_Mountain1Site0.03.txt",header=TRUE)->general_sample_Mountain1Site0.03

beta.multi(community_Hemiptera_Site0.03, index.family="sorensen")
#'repeating after removing outlayers from matrix and general habitat table
community_Hemiptera_Site0.03[-which(row.names(community_Hemiptera_Site0.03) %in% c("CON_NTO_TLC_36TCONS6", "CON_NTO_TLC_31TCONS1")),]->community_Hemiptera_sinoutlayer0.03
community_Hemiptera_sinoutlayer0.03[,which(colSums(community_Hemiptera_sinoutlayer0.03)!=0)]->community_Hemiptera_sinoutlayer0.03 #'to remove no data colums
dim(community_Hemiptera_sinoutlayer0.03)  

general_sample_Mountain1Site0.03[-which(general_sample_Mountain1Site0.03$sample_names %in% c("CON_NTO_TLC_36TCONS6","CON_NTO_TLC_31TCONS1")),]->general_sample_sinoutlayer0.03

beta.pair(community_Hemiptera_sinoutlayer0.03, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.03
plot (MDSbetasim0.03, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer0.03,ordispider(MDSbetasim0.03, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.430 ***", cex=2)


read.table("../genetic/Data_out/Hemiptera/Hemiptera5P/community_Hemiptera0.05.txt")->community_Hemiptera0.05
read.table("../genetic/Data_out/Hemiptera/Hemiptera5P/general_sample_Mountain1Site_0.05.txt",header=TRUE)->general_sample_Mountain1Site_0.05

beta.multi(community_Hemiptera_Site0.05, index.family="sorensen")
#'repeating after removing outlayers from matrix and general habitat table
community_Hemiptera_Site0.05[-which(row.names(community_Hemiptera_Site0.05) %in% c("CON_NTO_TLC_36TCONS6", "CON_NTO_TLC_31TCONS1")),]->community_Hemiptera_sinoutlayer0.05
community_Hemiptera_sinoutlayer0.05[,which(colSums(community_Hemiptera_sinoutlayer0.05)!=0)]->community_Hemiptera_sinoutlayer0.05 #'to remove no data colums
dim(community_Hemiptera_sinoutlayer0.05)  

general_sample_Mountain1Site_0.05[-which(general_sample_Mountain1Site_0.05$sample_names %in% c("CON_NTO_TLC_36TCONS6", "CON_NTO_TLC_31TCONS1")),]->general_sample_sinoutlayer0.05

beta.pair(community_Hemiptera_sinoutlayer0.05, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.05
plot (MDSbetasim0.05, xlim=c(-0.50, 0.5), ylim=c(-0.5, 0.5), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer0.05,ordispider(MDSbetasim0.05, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.418 ***", cex=2)
mtext(c("Hemiptera"), side = 4, col = "black", line = 1, cex = 2)

#'**Hymenoptera'**#
read.table("../genetic/Data_out/Hymenoptera/Hymenoptera_Haplotypes/community_Hymenoptera_h.txt")->community_Hymenoptera_h
read.table("../genetic/Data_out/Hymenoptera/Hymenoptera_Haplotypes/general_sample_Mountain1Site_h.txt",header=TRUE)->general_sample_Mountain1Site_h

beta.multi(community_Hymenoptera_Site_h, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Hymenoptera_Site_h, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_h #'NMDS
#'repeating after removing outlayers from matrix and general habitat table
community_Hymenoptera_Site_h[-which(row.names(community_Hymenoptera_Site_h) %in% c("CON_NTO_AAB_122ACON12")),]->community_Hymenoptera_sinoutlayer_h
community_Hymenoptera_sinoutlayer_h[,which(colSums(community_Hymenoptera_sinoutlayer_h)!=0)]->community_Hymenoptera_sinoutlayer_h #'to remove no data colums
dim(community_Hymenoptera_sinoutlayer_h)  

general_sample_Mountain1Site_h[-which(general_sample_Mountain1Site_h$sample_names %in% c("CON_NTO_AAB_122ACON12")),]->general_sample_sinoutlayer_h

beta.pair(community_Hymenoptera_sinoutlayer_h, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_h
plot (MDSbetasim_h, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer_h,ordispider(MDSbetasim_h, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.297 ***", cex=2)


read.table("../genetic/Data_out/Hymenoptera/Hymenoptera3P/community_Hymenoptera0.03.txt")->community_Hymenoptera0.03
read.table("../genetic/Data_out/Hymenoptera/Hymenoptera3P/general_sample_Mountain1Site0.03.txt",header=TRUE)->general_sample_Mountain1Site0.03

beta.multi(community_Hymenoptera_Site0.03, index.family="sorensen")
#'repeating after removing outlayers from matrix and general habitat table
community_Hymenoptera_Site0.03[-which(row.names(community_Hymenoptera_Site0.03) %in% c("CON_NTO_AAB_122ACON12", "CON_NTO_AAB_124ACON14")),]->community_Hymenoptera_sinoutlayer0.03
community_Hymenoptera_sinoutlayer0.03[,which(colSums(community_Hymenoptera_sinoutlayer0.03)!=0)]->community_Hymenoptera_sinoutlayer0.03 #'to remove no data colums
dim(community_Hymenoptera_sinoutlayer0.03)  

general_sample_Mountain1Site0.03[-which(general_sample_Mountain1Site0.03$sample_names %in% c("CON_NTO_AAB_122ACON12", "CON_NTO_AAB_124ACON14")),]->general_sample_sinoutlayer0.03

beta.pair(community_Hymenoptera_sinoutlayer0.03, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.03
plot (MDSbetasim0.03, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.42), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer0.03,ordispider(MDSbetasim0.03, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.192 ***", cex=2)


read.table("../genetic/Data_out/Hymenoptera/Hymenoptera5P/community_Hymenoptera0.05.txt")->community_Hymenoptera0.05
read.table("../genetic/Data_out/Hymenoptera/Hymenoptera5P/general_sample_Mountain1Site0.05.txt",header=TRUE)->general_sample_Mountain1Site0.05

beta.multi(community_Hymenoptera_Site0.05, index.family="sorensen")
#'repeating after removing outlayers from matrix and general habitat table
community_Hymenoptera_Site0.05[-which(row.names(community_Hymenoptera_Site0.05) %in% c("CON_NTO_AAB_122ACON12", "CON_NTO_AAB_124ACON14", "CON_NTO_TLC_36TCONS6")),]->community_Hymenoptera_sinoutlayer0.05
community_Hymenoptera_sinoutlayer0.05[,which(colSums(community_Hymenoptera_sinoutlayer0.05)!=0)]->community_Hymenoptera_sinoutlayer0.05 #'to remove no data colums
dim(community_Hymenoptera_sinoutlayer0.05)  

general_sample_Mountain1Site0.05[-which(general_sample_Mountain1Site0.05$sample_names %in% c("CON_NTO_AAB_122ACON12", "CON_NTO_AAB_124ACON14", "CON_NTO_TLC_36TCONS6")),]->general_sample_sinoutlayer0.05

beta.pair(community_Hymenoptera_sinoutlayer0.05, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.05
plot (MDSbetasim0.05, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer0.05,ordispider(MDSbetasim0.05, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.16, y=-0.42, labels = mylabel, cex=2)
text(x=0.56, y=-0.42, labels="=0.154 ***", cex=2)
mtext(c("Hymenoptera"), side = 4, col = "black", line = 1, cex = 2)


#'**Coleoptera'**#
read.table("../genetic/Data_out/Coleoptera/Coleoptera_Haplotypes/community_Coleoptera_h.txt")->community_Coleoptera_h
read.table("../genetic/Data_out/Coleoptera/Coleoptera_Haplotypes/general_sample_Mountain1Site_h.txt",header=TRUE)->general_sample_Mountain1Site_h

beta.multi(community_Coleoptera_Site_h, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Coleoptera_Site_h, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_h #'NMDS
#'repeating after removing outlayers from matrix and general habitat table
community_Coleoptera_Site_h[-which(row.names(community_Coleoptera_Site_h) %in% c("CON_NTO_AAB_120ACON10", "CON_NTO_TLC_32TCONS2")),]->community_Coleoptera_sinoutlayer
community_Coleoptera_sinoutlayer[,which(colSums(community_Coleoptera_sinoutlayer)!=0)]->community_Coleoptera_sinoutlayer #'to remove no data colums
dim(community_Coleoptera_sinoutlayer)  

general_sample_Mountain1Site_h[-which(general_sample_Mountain1Site_h$sample_names %in% c("CON_NTO_AAB_120ACON10", "CON_NTO_TLC_32TCONS2")),]->general_sample_sinoutlayer

beta.pair(community_Coleoptera_sinoutlayer, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_h

plot (MDSbetasim_h, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer,ordispider(MDSbetasim_h, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.422 ***", cex=2)

read.table("../genetic/Data_out/Coleoptera/Coleoptera3P/community_Coleoptera0.03.txt")->community_Coleoptera0.03
read.table("../genetic/Data_out/Coleoptera/Coleoptera3P/general_sample_Mountain1Site0.03.txt",header=TRUE)->general_sample_Mountain1Site0.03

beta.multi(community_Coleoptera_Site0.03, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Coleoptera_Site0.03, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.03 #'NMDS
#'repeating after removing outlayers from matrix and general habitat table
community_Coleoptera_Site0.03[-which(row.names(community_Coleoptera_Site0.03) %in% c("CON_NTO_AAB_122ACON12")),]->community_Coleoptera_sinoutlayer0.03
community_Coleoptera_sinoutlayer0.03[,which(colSums(community_Coleoptera_sinoutlayer0.03)!=0)]->community_Coleoptera_sinoutlayer0.03 #'to remove no data colums
dim(community_Coleoptera_sinoutlayer0.03)  

general_sample_Mountain1Site0.03[-which(general_sample_Mountain1Site0.03$sample_names %in% c("CON_NTO_AAB_122ACON12")),]->general_sample_sinoutlayer0.03

beta.pair(community_Coleoptera_sinoutlayer0.03, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.03
plot (MDSbetasim0.03, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.51), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer0.03,ordispider(MDSbetasim0.03, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.45, labels = mylabel, cex=2)
text(x=0.56, y=-0.45, labels="=0.360 ***", cex=2)


read.table("../genetic/Data_out/Coleoptera/Coleoptera5P/community_Coleoptera0.05.txt")->community_Coleoptera0.05
read.table("../genetic/Data_out/Coleoptera/Coleoptera5P/general_sample_Mountain1Site0.05.txt",header=TRUE)->general_sample_Mountain1Site0.05


beta.multi(community_Coleoptera_Site0.05, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Coleoptera_Site0.05, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.05 #'NMDS
#'repeating after removing outlayers from matrix and general habitat table
community_Coleoptera_Site0.05[-which(row.names(community_Coleoptera_Site0.05) %in% c("CON_NTO_AAB_122ACON12", "CON_NTO_AAB_108ACON21")),]->community_Coleoptera_sinoutlayer0.05
community_Coleoptera_sinoutlayer0.05[,which(colSums(community_Coleoptera_sinoutlayer0.05)!=0)]->community_Coleoptera_sinoutlayer0.05 #'to remove no data colums
dim(community_Coleoptera_sinoutlayer0.05)  

general_sample_Mountain1Site0.05[-which(general_sample_Mountain1Site0.05$sample_names %in% c("CON_NTO_AAB_122ACON12", "CON_NTO_AAB_108ACON21")),]->general_sample_sinoutlayer0.05

beta.pair(community_Coleoptera_sinoutlayer0.05, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.05
plot (MDSbetasim0.05, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer0.05,ordispider(MDSbetasim0.05, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.45, labels = mylabel, cex=2)
text(x=0.56, y=-0.45, labels="=0.374 ***", cex=2)
mtext(c("Coleoptera"), side = 4, col = "black", line = 1, cex = 2)

#'**Myriapoda'**#

read.table("../genetic/Data_out/Myriapoda/Myriapoda_Haplotypes/community_Myriapoda_h.txt")->community_Myriapoda_h
read.table("../genetic/Data_out/Myriapoda/Myriapoda_Haplotypes/general_sample_Mountain1Site_h.txt",header=TRUE)->general_sample_Mountain1Site_h

beta.multi(community_Myriapoda_Site_h, index.family="sorensen")
#'repeating after removing outlayers from matrix and general habitat table
community_Myriapoda_Site_h[-which(row.names(community_Myriapoda_Site_h) %in% c("CON_NTO_ASB_86ACON8","CON_NTO_TLC_37TCONS7","CON_NTO_ASB_79ACON1","CON_NTO_AAB_125ACON15")),]->community_Myriapoda_sinoutlayer_h
community_Myriapoda_sinoutlayer_h[,which(colSums(community_Myriapoda_sinoutlayer_h)!=0)]->community_Myriapoda_sinoutlayer_h #'to remove no data colums
dim(community_Myriapoda_sinoutlayer_h)  

general_sample_Mountain1Site_h[-which(general_sample_Mountain1Site_h$sample_names %in% c("CON_NTO_ASB_86ACON8","CON_NTO_TLC_37TCONS7","CON_NTO_ASB_79ACON1","CON_NTO_AAB_125ACON15")),]->general_sample_sinoutlayer_h

beta.pair(community_Myriapoda_sinoutlayer_h, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_h
plot (MDSbetasim_h, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer_h,ordispider(MDSbetasim_h, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.096 ns", cex=2)


read.table("../genetic/Data_out/Myriapoda/Myriapoda3P/community_Myriapoda0.03.txt")->community_Myriapoda0.03
read.table("../genetic/Data_out/Myriapoda/Myriapoda3P/general_sample_Mountain1Site0.03.txt",header=TRUE)->general_sample_Mountain1Site0.0

beta.multi(community_Myriapoda_Site0.03, index.family="sorensen")
#'repeating after removing outlayers from matrix and general habitat table
community_Myriapoda_Site0.03[-which(row.names(community_Myriapoda_Site0.03) %in% c("CON_NTO_ASB_79ACON1")),]->community_Myriapoda_sinoutlayer0.03
community_Myriapoda_sinoutlayer0.03[,which(colSums(community_Myriapoda_sinoutlayer0.03)!=0)]->community_Myriapoda_sinoutlayer0.03 #'to remove no data colums
dim(community_Myriapoda_sinoutlayer0.03)  

general_sample_Mountain1Site0.03[-which(general_sample_Mountain1Site0.03$sample_names %in% c("CON_NTO_ASB_79ACON1")),]->general_sample_sinoutlayer0.03

beta.pair(community_Myriapoda_sinoutlayer0.03, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.03
plot (MDSbetasim0.03, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer0.03,ordispider(MDSbetasim0.03, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.034 ns", cex=2)


read.table("../genetic/Data_out/Myriapoda/Myriapoda5P/community_Myriapoda0.05.txt")->community_Myriapoda0.05
read.table("../genetic/Data_out/Myriapoda/Myriapoda5P/general_sample_Mountain1Site_0.05.txt",header=TRUE)->general_sample_Mountain1Site_0.05

beta.multi(community_Myriapoda_Site0.05, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Myriapoda_Site0.05, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.05 #'NMDS
plot (MDSbetasim0.05, xlim=c(-0.5, 0.5), ylim=c(-0.5, 0.5), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_Mountain1Site_0.05,ordispider(MDSbetasim0.05, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.116 ns", cex=2)
mtext(c("Myriapoda"), side = 4, col = "black", line = 1, cex = 2)


#'**Lepidoptera'**#

read.table("../genetic/Data_out/Lepidoptera/Lepidoptera_Haplotypes/community_Lepidoptera_h.txt")->community_Lepidoptera_h
read.table("../genetic/Data_out/Lepidoptera/Lepidoptera_Haplotypes/general_sample_Mountain1Site_h.txt",header=TRUE)->general_sample_Mountain1Site_h

beta.multi(community_Lepidoptera_Site_h, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Lepidoptera_Site_h, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
#'repeating after removing outlayers from matrix and general habitat table
community_Lepidoptera_Site_h[-which(row.names(community_Lepidoptera_Site_h) %in% c("CON_NTO_SJH_76SHCON7", "CON_NTO_ASB_86ACON8", "CON_NTO_TLC_39TCONS9", "CON_NTO_SJH_74SHCON5", "CON_NTO_AAB_107ACON20",  "CON_NTO_AAB_111ACON24", "CON_NTO_TLC_34TCONS4")),]->community_Lepidoptera_sinoutlayer_h
community_Lepidoptera_sinoutlayer_h[,which(colSums(community_Lepidoptera_sinoutlayer_h)!=0)]->community_Lepidoptera_sinoutlayer_h #'to remove no data colums
dim(community_Lepidoptera_sinoutlayer_h)  

general_sample_Mountain1Site_h[-which(general_sample_Mountain1Site_h$sample_names %in% c("CON_NTO_SJH_76SHCON7", "CON_NTO_ASB_86ACON8", "CON_NTO_TLC_39TCONS9", "CON_NTO_SJH_74SHCON5", "CON_NTO_AAB_107ACON20",  "CON_NTO_AAB_111ACON24", "CON_NTO_TLC_34TCONS4")),]->general_sample_sinoutlayer_h

beta.pair(community_Lepidoptera_sinoutlayer_h, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim_h
plot (MDSbetasim_h, xlim=c(-0.5, 0.5), ylim=c(-0.49, 0.49), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer_h,ordispider(MDSbetasim_h, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.322 *", cex=2)


read.table("../genetic/Data_out/Lepidoptera/Lepidoptera3P/community_Lepidoptera0.03.txt")->community_Lepidoptera0.03
read.table("../genetic/Data_out/Lepidoptera/Lepidoptera3P/sample_names_Mountain1_0.03.txt",header=TRUE)->sample_names_Mountain1_0.03

beta.multi(community_Lepidoptera_Site0.03, index.family="sorensen")
#'turnover by pairs, nmds, anosim
beta.pair(community_Lepidoptera_Site0.03, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
#'repeating after removing outlayers from matrix and general habitat table
community_Lepidoptera_Site0.03[-which(row.names(community_Lepidoptera_Site0.03) %in% c("CON_NTO_SJH_76SHCON7", "CON_NTO_ASB_86ACON8", "CON_NTO_TLC_39TCONS9", "CON_NTO_SJH_74SHCON5", "CON_NTO_AAB_107ACON20",  "CON_NTO_AAB_111ACON24")),]->community_Lepidoptera_sinoutlayer0.03
community_Lepidoptera_sinoutlayer0.03[,which(colSums(community_Lepidoptera_sinoutlayer0.03)!=0)]->community_Lepidoptera_sinoutlayer0.03 #'to remove no data colums
dim(community_Lepidoptera_sinoutlayer0.03)  

sample_names_Mountain1_0.03[-which(sample_names_Mountain1_0.03$sample_names %in% c("CON_NTO_SJH_76SHCON7", "CON_NTO_ASB_86ACON8", "CON_NTO_TLC_39TCONS9", "CON_NTO_SJH_74SHCON5", "CON_NTO_AAB_107ACON20",  "CON_NTO_AAB_111ACON24")),]->general_sample_sinoutlayer0.03

beta.pair(community_Lepidoptera_sinoutlayer0.03, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.03

plot (MDSbetasim0.03, xlim=c(-0.5, 0.5), ylim=c(-0.4, 0.4), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer0.03,ordispider(MDSbetasim0.03, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))

#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.35, labels = mylabel, cex=2)
text(x=0.56, y=-0.35, labels="=0.174 *", cex=2)

read.table ("../genetic/Data_out/Lepidoptera/Lepidoptera5P/community_Lepidoptera0.05.txt")->community_Lepidoptera0.05
read.table("../genetic/Data_out/Lepidoptera/Lepidoptera5P/general_sample_Mountain1Site_0.05.txt",header=TRUE)->general_sample_Mountain1Site_0.05

beta.multi(community_Lepidoptera_Site0.05, index.family="sorensen")
beta.pair(community_Lepidoptera_Site0.05, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
#'repeating after removing outlayers from matrix and general habitat table
community_Lepidoptera_Site0.05[-which(row.names(community_Lepidoptera_Site0.05) %in% c("CON_NTO_SJH_76SHCON7", "CON_NTO_ASB_86ACON8", "CON_NTO_TLC_39TCONS9", "CON_NTO_SJH_74SHCON5", "CON_NTO_AAB_107ACON20",  "CON_NTO_AAB_111ACON24")),]->community_Lepidoptera_sinoutlayer0.05
community_Lepidoptera_sinoutlayer0.05[,which(colSums(community_Lepidoptera_sinoutlayer0.05)!=0)]->community_Lepidoptera_sinoutlayer0.05 #'to remove no data colums
dim(community_Lepidoptera_sinoutlayer0.05)  

general_sample_Mountain1Site_0.05[-which(general_sample_Mountain1Site_0.05$sample_names %in% c("CON_NTO_SJH_76SHCON7", "CON_NTO_ASB_86ACON8", "CON_NTO_TLC_39TCONS9", "CON_NTO_SJH_74SHCON5", "CON_NTO_AAB_107ACON20",  "CON_NTO_AAB_111ACON24")),]->general_sample_sinoutlayer0.05
beta.pair(community_Lepidoptera_sinoutlayer0.05, index.family="sorensen")->beta.pair  #'betadiversity by pair of communities using sorensen on the precense/absence data, with estimation of turnover and nestedness datamatrixes simultaneously
metaMDS (beta.pair$beta.sim)->MDSbetasim0.05
plot (MDSbetasim0.05, xlim=c(-0.5, 0.5), ylim=c(-0.49, 0.49), cex=2, cex.lab=1, cex.axis=2.3, lwd=4.5)
with(general_sample_sinoutlayer0.05,ordispider(MDSbetasim0.05, Site, cex.lab=1, col= c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"), lwd=4.5))
legend(0.21, 0.55,
       legend = c("Agua Bendita", "San Bartolo", "San Juan Huertas", "Tlacotepec"),
       col=c("#153a7b", "#eaa22f", "#97518b", "#81b9d0"),
       bty="n", text.font=2, lty=1, lwd= 7, cex = 1.1)

#I put letter "r" in cursive and r2 value
mylabel = bquote(italic(r)^2)
text(x=0.25, y=-0.42, labels = mylabel, cex=2)
text(x=0.56, y=-0.42, labels="=0.114 ns", cex=2)
mtext(c("Lepidoptera"), side = 4, col = "black", line = 1, cex = 2)

dev.off()



####################################END######################################################################



