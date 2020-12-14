#!/bin/bash

#map reads to mite genomes and phasable region for ASD project to check if phased regions double cov (paralogs)

#folder needed:
# reads
# mapping
# genomes


#program versions
#bwa 0.7.17
#samtools 1.7


clear
while read line; do
  read -a twoSp <<< $line
  printf "\n****processing ${twoSp[0]}*****\n"
  SP1=${twoSp[0]}
  SP2=${twoSp[1]}

LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH

#index genomes
cd ./genome
#for f in *.fasta; do bwa index $f; done
#for f in *.fasta; do bowtie2-build $f $f; done
cd ..


#mapping (-t is number of CPUs)
printf "**mapping $SP1**\n\n"
#bwa mem -t 12 ./genome/$SP1''_b1v03.fasta ./reads/$SP2''.trim.180.pair1.fastq ./reads/$SP2''.trim.180.pair2.fastq > ./mapping/$SP1''_aln.sam;
bowtie2 --local -x ./genome/$SP1''_b1v03.fasta -q --local --threads 24 \
-1 ~/Data/mites/kraken_filtered/$SP2''.trim.180.krak.pair1.fastq.gz \
-2 ~/Data/mites/kraken_filtered/$SP2''.trim.180.krak.pair2.fastq.gz \
-S ./mapping/$SP1''_aln.sam;


#bowtie2 -x ./genome/On_contigs_maskedphases_all_b1v03.fasta -q --local --threads 12 \
#-1 ~/Data/mites/kraken_filtered/On6.trim.180.krak.pair1.fastq \
#-2 ~/Data/mites/kraken_filtered/On6.trim.180.krak.pair2.fastq \
#-S ./mapping/On_contigs_maskedphases_all_aln.sam


#exclude multi-mapped reads
cd ./mapping
#cat $SP1''_aln.sam | grep -v -e 'XA:Z:' -e 'SA:Z:' > $SP1''_aln_unq.sam
cd ..

#samtools
printf "**samtools $SP1**\n\n"
cd ./mapping
samtools view -@ 24 -q 10 -bS -F 4 $SP1''_aln.sam -o $SP1''_aln.bam #with qual MAPQ score filter filtering mostly unique reads
samtools sort -@ 24 $SP1''_aln.bam -o $SP1''_aln.sorted.bam
samtools index $SP1''_aln.sorted.bam
cd ..

#remove unnecessary files
  rm ./mapping/$SP1''_aln.sam;
  rm ./mapping/$SP1''_aln.bam;



#picard remove duplicates
printf "**remove duplicates $SP1**\n\n"
cd ./mapping
java -jar ~/Software/tools/picard.jar MarkDuplicates \
I=$SP1''_aln.sorted.bam \
O=$SP1''_aln.sorted.marked_duplicates.bam \
M=$SP1''.marked_dup_metrics.txt \
REMOVE_DUPLICATES=true \
REMOVE_SEQUENCING_DUPLICATES=true
#java -jar ~/Software/tools/picard.jar MarkDuplicates -I $SP1''_aln.sorted.bam -O $SP1''_aln.sorted.marked_duplicates.bam -M $SP1''.marked_dup_metrics.txt \
#-REMOVE_DUPLICATES true -REMOVE_SEQUENCING_DUPLICATES true
cd ..


#coverage
cd ./mapping
bedtools genomecov -d -ibam $SP1''_aln.sorted.marked_duplicates.bam -g ../genome/$SP1''_b1v03.fasta > $SP1''.cov
cat $SP1''.cov | awk 'NF==3{sum[$1]+=$3;sum2[$1]+=$3*$3;N[$1]++;}END{for(key in sum){print(key,sum[key]/N[key],sqrt((N[key]*sum2[key]-sum[key]*sum[key])/(N[key]*N[key]-N[key])));}}' | sort > $SP1''.cov.mn
cd ..



done < ./samplenames_phases_new




printf "****DONE****\n"



#afterwards a bit of editing:
cd ./mapping
#cat $SP1''.cov.mn | awk '{print $0, $SP1}' > $SP1''.cov.mn2
#cat $SP1''.cov.mn2 | awk '{printf "%s %.0f %.0f %s\n",$1,$2,$3,$4}' > $SP1''.cov.rd
cd ..
#cat On_contigs_maskedphases_all.cov.mn | awk '{print $0, "Contigs_phases_masked"}' > On_contigs_maskedphases_all.cov.mn2
#cat phaseable_regions.cov.mn | awk '{print $0, "PhasedRegion"}' > phaseable_regions.cov.mn2
#cat 1_On.cov.mn2 phaseable_regions.cov.mn2 > On_phreg.cov
#cat On_phreg.cov | awk '{printf "%s %.0f %.0f %s\n",$1,$2,$3,$4}' > On_phreg.cov.rd

