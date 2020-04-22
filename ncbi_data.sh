#!/bin/bash

#Returns the the data for accession number AF086833.2 in the FASTA format:
curl -s https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?id=AF086833.2\&db=nuccore\&rettype=fasta | head

efetch -db=nuccore -format=gb -id=AF086833 | head

# Accession number AF086833 in Genbank format.
efetch -db=nuccore -format=gb -id=AF086833 > AF086833.gb

# Accession number AF086833 in Fasta format.
efetch -db=nuccore -format=fasta -id=AF086833 > AF086833.fa

# efetch can take additional parameters and select a section of the sequence.
efetch -db=nuccore -format=fasta -id=AF086833 -seq_start=1 -seq_stop=3

#efetch can produce the sequence from reverse strands:
efetch -db=nuccore -format=fasta -id=AF086833 -seq_start=1 -seq_stop=5 -strand=1
efetch -db=nuccore -format=fasta -id=AF086833 -seq_start=1 -seq_stop=5 -strand=2

#entrez can search for the data associated with an accession number.
#These can be located in published papers.
esearch -help

esearch -db nucleotide -query PRJNA257197 | efetch -format=fasta > genomes.fa

esearch -db protein -query PRJNA257197 | efetch -format=fasta > proteins.fa

# the xtract tool in Entrez Direct allows navigating and selecting parts of an XML file.
efetch -db taxonomy -id 9606,7227,10090 -format xml | xtract -Pattern Taxon -first TaxId ScientificName GenbankCommonName Division


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