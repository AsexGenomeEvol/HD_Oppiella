## Table of contents

1. [QC pre-processing on raw data](#1_preprocessing)
	* [Paired end](#11_pe)
	* [Mate pair](#12_mp)
	
2. [Assembly](#2_assembly)
	* [Contig](#21_contig)
	* [Scaffolding and Gap-closing](#22_scaffold)

3. [Filtration](#3_filtration)
	* [Contamination removal](#31_contamination)		
	* [Size filtration](#32_size)	
	
4. [Evaluation](#4_evaluation)
					
## <a name="1_preprocessing"></a>1) QC pre-processing on raw data

#### <a name="11_pe"></a>1.1) Paired end

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) was used for paired end filtering.

```
java -jar trimmomatic.jar PE -threads 1 -phred33 *.raw.180.pair1.fastq *.raw.180.pair2.fastq \
trim/180/*.trim.180.pair1.fastq trim/180/*.trim.180.single1.fastq trim/180/*.trim.180.pair2.fastq trim/180/*.trim.180.single2.fastq \
ILLUMINACLIP:trimmomatic/0.36/adapters/all-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:3:20 MINLEN:100
```

#### <a name="12_mp"></a>1.2) Mate pair

[Nxtrim](https://github.com/sequencing/NxTrim) and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) were used for mate pair filtering.

```
module add UHTS/Quality_control/NxTrim/0.4.1

nxtrim -1 *.pair1.fastq -2 *.pair2.fastq -O name --separate --preserve-mp --minlength 40
```

Categories MP and UNKNOWN were concatenated (as suggested by the authors) then [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) was used for further fitlering.

```
java -jar trimmomatic.jar PE -threads 1 -phred33 *.nxtrim.3000.pair1.fastq *.nxtrim.3000.pair2.fastq \
*.nxtrim.trim.3000.pair1.fastq *.nxtrim.trim.3000.single1.fastq *.nxtrim.trim.3000.pair2.fastq *.nxtrim.trim.3000.single2.fastq \
ILLUMINACLIP:trimmomatic/0.36/adapters/all-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:60
```

## <a name="2_assembly"></a>2) Assembly

#### <a name="21_contig"></a>2.1) Contig

PCR whole genome amplification produce unveven coverage. [BBnorm](https://sourceforge.net/projects/bbmap/) was used to reduce that coverage bias as most of assemblers assume uniform coverage.

1) Reformat the overlapped library to interleaved format and merge the overlapped reads:

```
module add UHTS/Aligner/BBMap/36.59

reformat.sh in1=*.trim.180.pair1.fastq in2=*.trim.180.pair2.fastq out=*.trim.180.pair12.fastq 2> *.reformat.log

bbmerge.sh in=*.trim.180.pair12.fastq out=*.trim.merged.180.pair12.fastq minoverlap=15 mismatches=0 ecct strict -Xmx90g threads=35 2> *.bbmerge.log
```

2) Reverse complement the S2:

```
module add UHTS/Analysis/fastx_toolkit/0.0.13.2;
 
fastx_reverse_complement -i *.trim.180.single2.fastq -o *.trim.rc.180.single2.fastq -Q33
```

3) Concatenate the merged PE and the single reads:

```
cat *.trim.merged.180.pair12.fastq *.trim.180.single1.fastq *.trim.rc.180.single2.fastq > *.trim.180.all.fastq
```

4) Normalized the concatenate file (composed of single reads):

```
bbnorm.sh in=*.trim.180.all.fastq out=*.trim.norm.180.all.fastq target=65 min=3 -Xmx150g threads=20 2> bbnorm.log
```

5) Assembly with [SPAdes](http://bioinf.spbau.ru/spades):

```
export PATH=$PATH:/scratch/beegfs/monthly/ptranvan/Software/SPAdes-3.10.1-Linux/bin
 
spades.py -m 400 -t 27 --careful -k 21,33,55,77,99,111,127 -o * --s1 *.trim.norm.180.all.fastq

```

#### <a name="22_scaffold"></a>2.2) Scaffolding and Gap-closing

The other libraries (where insert_size>2*len(reads)) were used for scaffolding and gap-closing.

1) Scaffolding with [SSPACE](https://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE):

```
export PERL5LIB=/home/ptranvan/perl5/lib/perl5
module add UHTS/Aligner/bwa/0.7.13
 
perl /scratch/beegfs/monthly/ptranvan/Software/sspace/3.0/SSPACE_Standard_v3.0.pl -l *.sspace.txt -s */contigs.fasta -b * -T 25
```

2) Gap-closing with [Gapcloser](http://soap.genomics.org.cn/soapdenovo.html):

```
/scratch/beegfs/monthly/ptranvan/Software/GapCloser/1.12-r6/GapCloser -a ../*.final.scaffolds.fasta -b *.gapclose.txt -o *v01.fasta -l 125 -t 25 >> sm2.gapclose.log
```

## <a name="3_filtration"></a>3) Filtration

#### <a name="31_contamination"></a>3.1) Contamination removal

[BlobTools](https://github.com/DRL/blobtools) was used for contamination checking.

1) Map reads back to the genome:

```
module add UHTS/Aligner/bwa/0.7.15

bwa mem -M db_bwa/* *.trim.350.pair1.fastq.gz *.trim.350.pair2.fastq.gz -t 15 2> bwa/*.trim.350.pair12.bwa.output.log | samtools view -bS - > bwa/*.trim.350.pair12.bam
```

2) Compute coverage for each scaffold:

```
cat ../*.sam > merge.sam

module add UHTS/Analysis/BBMap/37.82
module add UHTS/Analysis/samtools/1.4

pileup.sh in=merge.sam out=merge.stats.txt hist=histogram.txt

# Reformat the output to <Name_scaffold>\t<Coverage>

awk {'printf ("%s\t%s\n", $1, $2)'} merge.stats.txt | awk '{if(NR>1)print}' > merge.stats_parse.txt
```

3) Blastn against nt:

```
module add Blast/ncbi-blast/2.7.1+;
 
export BLASTDB=/archive/dee/schwander/ptranvan/Database/taxdb:$BLASTDB

blastn -query *v01.fasta -db /archive/dee/schwander/ptranvan/Database/nt/nt \
-outfmt '6 qseqid staxids bitscore evalue std sscinames sskingdoms stitle' -max_target_seqs 10 \
-max_hsps 1 -evalue 1e-25 -num_threads 15 -out *.vs.nt.max10.1e25.blastn.out
```

4) Run BlobTools:

```
source /scratch/beegfs/monthly/ptranvan/Software/blobtools/1.0.sh

blobtools create -i *v01.fasta -t *.vs.nt.max10.1e25.blastn.out \
--nodes /scratch/beegfs/monthly/ptranvan/Software/blobtools/1.0/data/nodes.dmp \
--names /scratch/beegfs/monthly/ptranvan/Software/blobtools/1.0/data/names.dmp \
-c merge.stats_parse.txt -x bestsumorder -o *

blobtools blobplot -i *.blobDB.json --sort count --hist count -x bestsumorder

blobtools view -i *.blobDB.json --hits --rank all -x bestsumorder
```

5) Scaffolds without hits to metazoans were filtered out.

```
python contamination_filtration.py -s contamination_identification -i1 *.blobDB.table.txt

module add UHTS/Analysis/BBMap/37.82

filterbyname.sh in=*v01.fasta names=contaminant_scaffolds.txt out=*.filtered.fasta include=f -Xmx20g

```

6) Sort by descreasing size and rename the scaffolds header.

```
seqkit sort --by-length --reverse *.filtered.fasta -o - | rename_headers.awk -v sp=* version=b1v02 - > *v02.fasta
```

#### <a name="32_size"></a>3.2) Size filtration

1) Scaffolds > 500 bp were kept:

```
prinseq-lite/0.20.4/prinseq-lite.pl -fasta *.fasta -min_len 501 -out_good *.501 -out_bad *.bad 2> prinseq.stat       
```

2) Sort by descreasing size and rename the scaffolds header.

```
seqkit sort --by-length --reverse *.501 -o - | rename_headers.awk -v sp=* version=b1v03 - > *v03.fasta
```                

:arrow_right: v03 is the final assembly.

## <a name="4_evaluation"></a>4) Evaluation

[Cegma](http://korflab.ucdavis.edu/datasets/cegma/):

```
module add SequenceAnalysis/CEGMA/2.5

cegma -g *v03.fasta -T 10 -o *v03 2> cegma.output_log.txt
```

[Busco](http://busco.ezlab.org/):

```
source /scratch/beegfs/monthly/ptranvan/Software/busco/3.0.2b.sh

run_BUSCO.py --long -i *v03.fasta -o output -l /scratch/beegfs/monthly/ptranvan/Software/busco/3.0/arthropoda_odb9 -m geno -c 10

```
