#!/bin/bash


# Make a new directory for the reference genomes.
mkdir -p ~/refs

# Assign the file name to a variable
ACC=AF086833

# The reference genome stored in the local variable "REF", using its absolute location.
REF=~/refs/${ACC}.fa

# Fetch the reference ebola sequence, using the efetch tool of the Entrez-Direct tool suite.
# The ebola sequence is retrieved from the nucleotide Entrez database, and stored in the FASTA format in the local variable "REF"
efetch -db=nuccore -format=fasta -id=$ACC > $REF

#Construct the BWT index of the Ebola reference genome.
bwa index $REF

# Obtain all runs from the Ebola project using the esearch tool of the Entrez Direct tool suite.
# Interrogate the SRA database for all runs matching bioproject accession " PRJNA257197"
# Pipe this output to efetch and download all metadata for all runs in CSV format.
esearch -db sra -query PRJNA257197  | efetch -format runinfo > runinfo.csv

# Download first 10,000 reads from run SRR1972739 and convert to FASTQ format.
# Split paired-end data to yield SRR1972739_1.fastq and SRR1972739_2.fastq
fastq-dump -X 10000 --split-files SRR1972739

# Create local shortcuts to files.
R1=SRR1972739_1.fastq
R2=SRR1972739_2.fastq

# Align the data to the reference genome using BWA-MEM, which seeds alignments with maximal exact matches (MEMs) and then extends seeds with the affine-gap Smith-Wtaerman algorithm (SW).
# Store output in bwa.sam
bwa mem $REF $R1 $R2 > bwa.sam
