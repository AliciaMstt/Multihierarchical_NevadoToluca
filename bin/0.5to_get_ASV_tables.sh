#SBATCH --mem=64000
#SBATCH -n 10

#1Diptera2059Reads_realigned.fasta

vsearch --search_exact LENGTH420CatTratamientos137Sopas.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db 1Diptera2059Reads_realigned.fasta --otutabout allele_table_1Diptera2059Reads_realigned.txt


#2Collembola1220Reads_realigned.fasta

vsearch --search_exact LENGTH420CatTratamientos137Sopas.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db 2Collembola1220Reads_realigned.fasta --otutabout allele_table_2Collembola1220Reads_realigned.txt 

#3Arachnida428Reads_realigned.fasta

vsearch --search_exact LENGTH420CatTratamientos137Sopas.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db 3Arachnida428Reads_realigned.fasta --otutabout allele_table_3Arachnida428Reads_realigned.txt

#4Hymenoptera304Reads_realigned.fasta

vsearch --search_exact LENGTH420CatTratamientos137Sopas.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db 4Hymenoptera304Reads_realigned.fasta --otutabout allele_table_4Hymenoptera304Reads_realigned.txt


#5Coleoptera227Reads_realigned.fasta

vsearch --search_exact LENGTH420CatTratamientos137Sopas.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db 5Coleoptera227Reads_realigned.fasta --otutabout allele_table_5Coleoptera227Reads_realigned.txt 


#6Paraneoptera210Reads_realigned.fasta

vsearch --search_exact LENGTH420CatTratamientos137Sopas.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db 6Paraneoptera210Reads_realigned.fasta --otutabout allele_table_6Paraneoptera210Reads_realigned.txt 

#7Lepidoptera59Reads_realigned.fasta

vsearch --search_exact LENGTH420CatTratamientos137Sopas.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db 7Lepidoptera59Reads_realigned.fasta --otutabout allele_table_7Lepidoptera59Reads_realigned.txt 

#8Myriapoda57Reads_realigned.fasta

vsearch --search_exact LENGTH420CatTratamientos137Sopas.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db 8Myriapoda57Reads_realigned.fasta --otutabout allele_table_8Myriapoda57Reads_realigned.txt 

#
