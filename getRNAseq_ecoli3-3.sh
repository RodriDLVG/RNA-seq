#####################################################################################
#####################################################################################
##### Eukaryotic gene expression analysis ######
##### Install bowtie2
 sudo apt install bowtie2

##### Install EBseq
git clone https://github.com/deweylab/RSEM.git
sudo chmod -R 777 RSEM
# add path to 
#export PATH="/u01/RSEM:$PATH"
#source ~/.bashrc
## Follow the instructions here but starting from "Build References"

# https://github.com/bli25broad/RSEM_tutorial
#####################################################################################
##### Prokaryotic gene expression analysis ######

##### RNAseq using UCSD's E coli experimental dataset
# create a folder for this analysis

mkdir RNAseq && cd RNAseq && ls

mkdir reads refgenome
cd refgenome
wget ftp://ftp.ensemblgenomes.org/pub/release-49/bacteria//gtf/bacteria_79_collection/escherichia_coli_str_k_12_substr_w3110_gca_000010245/Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.49.gtf.gz
wget ftp://ftp.ensemblgenomes.org/pub/release-49/bacteria//fasta/bacteria_79_collection/escherichia_coli_str_k_12_substr_w3110_gca_000010245/dna/Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.dna.chromosome.Chromosome.fa.gz
gunzip Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.dna.chromosome.Chromosome.fa.gz 
mv Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.dna.chromosome.Chromosome.fa ecW3110.fna
gunzip Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.49.gtf.gz 
mv Escherichia_coli_str_k_12_substr_w3110_gca_000010245.ASM1024v1.49.gtf ecW3110.gtf


cd ../reads
## data sets are here: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=SRP270849&o=acc_s%3Aa
## we will use three of wildtype and endpoint strains
# get to know the what each of the files mean in biology

fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files --accession SRR12170033 &
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files --accession SRR12170034 &
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files --accession SRR12170038 &
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files --accession SRR12170037 &
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files --accession SRR12170040 &
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files --accession SRR12170042 &

## quality control reads using fastqc

for i in *.fastq; do fastqc -t 8 -f fastq $i;done
# fastqc -t 8 -f fastq file_name optional_output_file_name
#for i in *.fastq 
#do
#    fastqc -t 8 -f fastq $i
#done

## Based on fastqc output, fileter and trim reads using Trimmomatic
# HEADCROP: Cut the specified number of bases from the start of the read
# TRAILING: Cut bases off the end of a read, if below a threshold quality
# SLIDINGWINDOW: Performs a sliding window trimming approach. It starts
#   scanning at the 5‟ end and clips the read once the average quality within the window
#   falls below a threshold. 
# MINLEN: Drop the read if it is below a specified length
# Trimmomatic is run wit 
   ACTGGGGGGGGGGNTTTCCCAATTTT
# path, like this: ~/miniconda3/bin/trimmomatic
# we are still in the ./reads folder
#xxxxx_1.fastq
#xxxxx_2.fastq
readR=${i/_1/_2}
new=${i%%fastq}q26tr.fq;
xxxxx_1.q26tr.fq
for i in *_1.fastq; do  readR=${i%%_1.fastq}_2.fastq; new=${i%%fastq}q26tr.fq; java -jar /home/test/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
-threads 8 $i $readR -baseout $new HEADCROP:14 TRAILING:26 SLIDINGWINDOW:4:26 MINLEN:75; done

######## Run RSEM #############
### create index
# download genome files
# cd ref 
# scp vcm@dku-vcm-914.vm.duke.edu:/u01/RNAseq/ref/ecoli* .
# scp vcm@dku-vcm-914.vm.duke.edu:/u01/RNAseq/reads/*.fastq .
# create index, needing genome fna and gtf files and bowtie2 path
# rsem will call bowtie-build of the bowtie package, so it is the bowtie-build path is required
cd ..
mkdir ref
sudo chmod -R 777 ref
rsem-prepare-reference --gtf refgenome/ecW3110.gtf --bowtie2 --bowtie2-path /usr/bin refgenome/ecW3110.fna ref/ecoliref
# assemble reads using RSEM with Bowtie2
# rsem-calculate-expression -p 1 --paired-end --bowtie2 --bowtie2-path /usr/bin/ --append-names --output-genome-bam reads/SRR12170033_1.q26tr_1P.fq reads/SRR12170033_1.q26tr_2P.fq ref/ecoliref exp/WT1 &
# rsem-calculate-expression -p 1 --paired-end --bowtie2 --bowtie2-path /usr/bin/ --append-names --output-genome-bam reads/reads/SRR12170034_1.q26tr_1P.fq reads/SRR12170034_1.q26tr_2P.fq ref/ecoliref exp/WT2 &
mkdir exp
sudo chmod -R 777 exp
cd 
for i in reads/*_1P.fq; do rev2=$(basename $i); r2=reads/${rev2%%1P.fq}2P.fq; out=${rev2%%_1.q26tr_1P.fq}; time rsem-calculate-expression -p 8 --paired-end --bowtie2 --bowtie2-path /usr/bin --append-names --output-genome-bam $i $r2 ref/ecoliref exp/$out ;done &
