# CBD_PROG_1: Ebola Virus Genomic Data Processing Pipeline

## Overview
The "CBD_PROG_1" project facilitates the automated fetching, processing, and analysis of genomic data related to the Ebola virus from the National Center for Biotechnology Information (NCBI) databases. The pipeline is encapsulated in a shell script (`ncbi_data.sh`) designed to streamline several bioinformatics tasks including the downloading of reference genomes, indexing for alignment, fetching sequencing run metadata, downloading and processing sequencing reads, and aligning these reads to the reference genome.

## Features
- Download the Ebola virus reference genome.
- Index the reference genome for efficient alignment.
- Fetch metadata for sequencing runs associated with specific bioprojects.
- Download and split sequencing reads for analysis.
- Align sequencing reads to the reference genome using BWA-MEM.

## Prerequisites
To use this pipeline, you will need:
- Unix-like operating system (Linux, macOS)
- Installed software:
  - [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/) for fetching data from NCBI.
  - [BWA](http://bio-bwa.sourceforge.net/) for indexing the reference genome and aligning reads.
  - [SRA Toolkit](https://github.com/ncbi/sra-tools) for downloading sequencing reads.
- Internet connection for data download.

## Installation
1. **Install Required Software**:
   - **Entrez Direct**: Follow the [installation guide](https://www.ncbi.nlm.nih.gov/books/NBK179288/).
   - **BWA**: Installation instructions are available on the [BWA homepage](http://bio-bwa.sourceforge.net/).
   - **SRA Toolkit**: Installation can be done following the instructions on [GitHub](https://github.com/ncbi/sra-tools).

2. **Download the Project**:
   ```bash
   git clone https://github.com/StewartONeill24/CBD_PROG_1.git
   cd CBD_PROG_1
   ```

## Usage
Before running the script, ensure all prerequisites are installed and properly configured. Execute the script from the terminal:

```bash
bash ncbi_data.sh
```

The script performs the following operations:
1. Creates a directory for storing reference genomes.
2. Downloads the Ebola virus reference genome using its accession number.
3. Indexes the downloaded reference genome for alignment.
4. Fetches metadata for sequencing runs from a specific bioproject.
5. Downloads and processes a subset of sequencing reads.
6. Aligns the reads to the reference genome and outputs the alignment in SAM format.

## Output
The pipeline generates several files throughout its execution:
- A FASTA file of the Ebola virus reference genome.
- Indexed reference genome files created by BWA.
- A CSV file containing metadata of sequencing runs.
- FASTQ files of the downloaded reads.
- A SAM file containing the alignment of reads to the reference genome.

## Contributing
Contributions to the "CBD_PROG_1" project are welcome. Please open an issue to discuss your ideas or submit a pull request with your improvements.

## License
This project is released under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments
- This pipeline utilizes data and tools provided by NCBI, demonstrating the power of public genomic databases and bioinformatics tools in viral research.
