
library(raster)
require(gridExtra)
library(sp)
library(ggplot2)
library(patchwork)
library(easypackages)
library(ggplot2)

#'**Community diversity and composition at the haplotype, 3% and 5% Clustering levels for each of the eight taxonomic orders studied**'#

#This script gets all scripts of diversity using 8 arthropods order at multi-hierarchical levels 
#We used 8 groups: Arachnda, Coleoptera, Collembola, Diptera, Hemiptera, Hymenoptera, Lepidoptera, and Myriapoda. 

################ plots Community diversity and composition #####################


png(filename="../figures/RelativeSpeciesRichness.png", width=920 , height=1536, units="px") # set size of the file to plot 
par(mfrow=c(2,3)) #the number of rows and columns the figure would have

#'**Creating plot'**#

#'**Community diversity and composition at the haplotype, 3% and 5% Clustering levels for each of the eight taxonomic orders studied**'#

#This script gets all scripts of diversity using 8 arthropods order at multi-hierarchical levels 
#We used 8 groups: Arachnda, Coleoptera, Collembola, Diptera, Hemiptera, Hymenoptera, Lepidoptera, and Myriapoda. 


################ plot Global Richness by Sites #####################
#'**Diptera'**#
plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h,ylab="richness_Site_h", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,154))
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
text(x=c(1,2,3,4.1), y=(153), labels=c("a","a","ab","b"), cex=1.4)
text(x=4.5, y=153, labels="*", cex=2)
mtext(c("Haplotypes"), side = 3, col = "black", line = 1, cex = 2)

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03,ylab="richness_Site0.03", ylim=c(0,50), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
text(x=4.5, y=49, labels="*", cex=2)
mtext(c("Clustering 3%"), side = 3, col = "black", line = 1, cex = 2)

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05,ylab="richness_Site0.05", ylim=c(0,50), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
text(x=c(1,2,3,4.1), y=(49), labels=c("ab","a","ab","b"), cex=1.4)
text(x=4.5, y=49, labels="*", cex=2)
mtext(c("Clustering 5%"), side = 3, col = "black", line = 1, cex = 2)
mtext(c("Diptera"), side = 4, col = "black", line = 1, cex = 2)


#'**Collemobola'**#
plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h,ylab="richness_Site_h", cex=1.4, cex.axis=2.3, lwd=2.5, ylim=c(0,60))
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4), y=(15.7), labels=c("a","b","ab","b"), cex=1.4)
text(x=4.5, y=59, labels="ns", cex=2)

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03,ylab="richness_Site0.03", ylim=c(0,25), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
text(x=c(1,2,3,4), y=(24), labels=c("a","a","a","b"), cex=1.4)
text(x=4.5, y=24, labels="*", cex=2)

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05,ylab="richness_Site0.05", ylim=c(0,25), cex=1.4, cex.axis=2.3, lwd=2.5)
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4), y=(24), labels=c("a","a","a","b"), cex=1.4)
text(x=4.5, y=24, labels="ns", cex=2)
mtext(c("Collembola"), side = 4, col = "black", line = 1, cex = 2)

dev.off()


#'**Arachnida'**#
plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h,ylab="richness_Site_h", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Arachnida")
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
text(x=c(1,2,3,4), y=(28), labels=c("a","a","b","b"), cex=1.4)
text(x=4.5, y=29.5, labels="***", cex=2)

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03,ylab="richness_Site0.03", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5, main="Richness_Site Arachnida_0.03")
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
text(x=c(1,2,3,4), y=(28), labels=c("a","a","b","b"), cex=1.4)
text(x=4.5, y=29.5, labels="***", cex=2)

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05,ylab="richness_Site0.05", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5, main="Richness_Site Arachnida_0.05")
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
text(x=c(1,2,3,4), y=(28), labels=c("a","a","b","b"), cex=1.4)
text(x=4.5, y=29.5, labels="***", cex=2)


#'**Hempitera'**#
plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h,ylab="richness_Site_h",  ylim=c(0,15), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Hemiptera")
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(153), labels=c("a","a","ab","b"), cex=1.4)
text(x=4.5, y=14, labels="ns", cex=2)

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03,ylab="richness_Site0.03", ylim=c(0,15), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Hemiptera_0.03")
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(153), labels=c("a","a","ab","b"), cex=1.4)
text(x=4.5, y=14, labels="ns", cex=2)

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05,ylab="richness_Site0.05", ylim=c(0,15), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site0.05 Hemiptera_0.05")
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(153), labels=c("a","a","ab","b"), cex=1.4)
text(x=4.5, y=14, labels="ns", cex=2)
mtext(c("Hemiptera"), side = 4, col = "black", line = 1, cex = 2)

#'**Hymenoptera'**#
plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h,ylab="richness_Site_h", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site_h Hymenoptera")
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
text(x=4.5, y=29, labels="**", cex=2)

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03,ylab="richness_Site0.03", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site0.03 Hymenoptera_0.03")
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
text(x=4.5, y=29, labels="*", cex=2)

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05,ylab="richness_Site0.05", ylim=c(0,30), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Hymenoptera_0.05")
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
text(x=4.5, y=29, labels="ns", cex=2)
mtext(c("Hymenoptera"), side = 4, col = "black", line = 1, cex = 2)


#'**Coleoptera'**#
plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h,ylab="richness_Site_h", ylim=c(0,16), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site_h Coleoptera")
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
text(x=c(1,2,3,4), y=(15.7), labels=c("a","b","ab","b"), cex=1.4)
text(x=4.5, y=15.7, labels="**", cex=2)

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03,ylab="richness_Site0.03", ylim=c(0,16), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Coleoptera_0.03")
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
text(x=c(1,2,3,4), y=(15.7), labels=c("a","b","ab","b"), cex=1.4)
text(x=4.5, y=15.7, labels="**", cex=2)

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05,ylab="richness_Site0.05", ylim=c(0,16), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Coleoptera_0.05")
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
text(x=c(1,2,3,4), y=(15.7), labels=c("a","b","ab","b"), cex=1.4)
text(x=4.5, y=15.7, labels="**", cex=2)
mtext(c("Coleoptera"), side = 4, col = "black", line = 1, cex = 2)

#'**Myriapoda'**#
plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h,ylab="richness_Site_h", ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Myriapoda")
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
text(x=4.5, y=7, labels="ns", cex=2)

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03,ylab="richness_Site0.03", ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Myriapoda_0.03")
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
text(x=4.5, y=7, labels="ns", cex=2)

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05,ylab="richness_Site0.05", ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Myriapoda_0.05")
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
text(x=4.5, y=7, labels="ns", cex=2)
mtext(c("Myriapoda"), side = 4, col = "black", line = 1, cex = 2)

#'**Lepidoptera'**#

plot(richness_Site_h$Mountain1Site,richness_Site_h$sample_richness_Site_h,ylab="richness_Site_h", ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Lepidoptera")
kruskal.test(sample_richness_Site_h ~ Mountain1Site, data = richness_Site_h)
posthoc.kruskal.nemenyi.test(x=richness_Site_h$sample_richness_Site_h, g=richness_Site_h$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
text(x=4.5, y=7, labels="ns", cex=2)

plot(richness_Site0.03$Mountain1Site,richness_Site0.03$sample_richness_Site0.03,ylab="richness_Site0.03", ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Lepidoptera_0.03")
kruskal.test(sample_richness_Site0.03 ~ Mountain1Site, data = richness_Site0.03)
posthoc.kruskal.nemenyi.test(x=richness_Site0.03$sample_richness_Site0.03, g=richness_Site0.03$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
text(x=4.5, y=7, labels="ns", cex=2)

plot(richness_Site0.05$Mountain1Site,richness_Site0.05$sample_richness_Site0.05,ylab="richness_Site0.05", ylim=c(0,8), cex=1.4, cex.axis=2.3, lwd=2.5, main="H richness_Site Lepidoptera_0.05")
kruskal.test(sample_richness_Site0.05 ~ Mountain1Site, data = richness_Site0.05)
posthoc.kruskal.nemenyi.test(x=richness_Site0.05$sample_richness_Site0.05, g=richness_Site0.05$Mountain1Site, method="Bonferroni")
# Comparison of each group against. 
#text(x=c(1,2,3,4.1), y=(29), labels=c("a","b","a","a"), cex=1.4)
text(x=4.5, y=7, labels="ns", cex=2)
mtext(c("Lepidoptera"), side = 4, col = "black", line = 1, cex = 2)


dev.off()


####################################END######################################################################



