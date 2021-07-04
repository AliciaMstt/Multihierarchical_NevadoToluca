#!/bin/bash

#SBATCH --mem=64000
#SBATCH -n 10

#1Diptera2031Reads_realigned.fasta

vsearch --search_exact ../genetic/Data_in/AlignedSeq/LENGTH420CatConservacion51SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db ../genetic/Data_in/AlignedSeq/1Diptera2031readsRealigned.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_1Diptera2031Reads_realigned.txt


#2Collembola1061Reads_realigned.fasta

vsearch --search_exact ../genetic/Data_in/AlignedSeq/LENGTH420CatConservacion51SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db ../genetic/Data_in/AlignedSeq/2Collembola1061readsRealigned.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_2Collembola1061Reads_realigned.txt 


#3Arachnida336Reads_realigned.fasta

vsearch --search_exact ../genetic/Data_in/AlignedSeq/LENGTH420CatConservacion51SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db ../genetic/Data_in/AlignedSeq/3Arachnida336ReadsRealigned.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_3Arachnida336Reads_realigned.txt


#4Hymenoptera227Reads_realigned.fasta

vsearch --search_exact ../genetic/Data_in/AlignedSeq/LENGTH420CatConservacion51SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db ../genetic/Data_in/AlignedSeq/4Hymenoptera227Realigned.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_4Hymenoptera227Reads_realigned.txt


#5Coleoptera209Reads_realigned.fasta

vsearch --search_exact ../genetic/Data_in/AlignedSeq/LENGTH420CatConservacion51SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db ../genetic/Data_in/AlignedSeq/5Coleoptera209readsRealigned.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_5Coleoptera209Reads_realigned.txt 


#6Paraneoptera246Reads_realigned.fasta

vsearch --search_exact ../genetic/Data_in/AlignedSeq/LENGTH420CatConservacion51SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db ../genetic/Data_in/AlignedSeq/6ParaneopteraHemiptera246Realigned.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_6Paraneoptera246Reads_realigned.txt 


#7Lepidoptera60Reads_realigned.fasta

vsearch --search_exact ../genetic/Data_in/AlignedSeq/LENGTH420CatConservacion51SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db ../genetic/Data_in/AlignedSeq/7Lepidoptera60Realigned.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_7Lepidoptera60Reads_realigned.txt 


#8Myriapoda60Reads_realigned.fasta

vsearch --search_exact ../genetic/Data_in/AlignedSeq/LENGTH420CatConservacion51SOPAS.Ptrim.trimmo.Maxee1.LABEL.sorted420NOPuntoComa.fas --db ../genetic/Data_in/AlignedSeq/8Myriapoda60readsRealigned.fasta --otutabout ../genetic/Data_out/ASVs_Table/allele_table_8Myriapoda60Reads_realigned.txt 

#
