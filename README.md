# **STEC_KMA**

[![CircleCI](https://dl.circleci.com/status-badge/img/gh/OLC-Bioinformatics/STEC_KMA/tree/main.svg?style=shield)](https://dl.circleci.com/status-badge/redirect/gh/OLC-Bioinformatics/STEC_KMA/tree/main)
[![codecov](https://codecov.io/gh/OLC-Bioinformatics/STEC_KMA/graph/badge.svg?token=EZWA724UQC)](https://codecov.io/gh/OLC-Bioinformatics/STEC_KMA)
[![Anaconda-Server Badge](https://img.shields.io/badge/install%20with-conda-brightgreen)](https://anaconda.org/olcbioinformatics/stec_kma)
[![GitHub Release](https://img.shields.io/github/v/release/OLC-Bioinformatics/STEC_KMA?display_name=release&label=version&color=%20dark%20green
)](https://github.com/OLC-Bioinformatics/STEC_KMA/releases)
[![GitHub issues](https://img.shields.io/github/issues/OLC-Bioinformatics/STEC_KMA)](https://github.com/OLC-LOC-Bioinformatics/STEC_KMA/issues)
[![license](https://img.shields.io/badge/license-MIT-brightgreen)](https://github.com/OLC-Bioinformatics/STEC_KMA/blob/main/LICENSE)

`STEC_KMA` is a bioinformatics tool designed to process sequencing data for Shiga toxin-producing *Escherichia coli* (STEC). It uses the **KMA** (K-mer Alignment) algorithm to map sequencing reads to a combined allele database, identify the best hits, extract relevant reads, and confirm results through alignment using **BWA** (Burrows-Wheeler Aligner). The tool also identifies alleles with insertions and generates comprehensive reports.

---

## **Features**
- **Read Baiting**: Uses BBMap to bait reads from sequencing data based on a reference allele database.
- **KMA Mapping**: Maps reads to allele databases using KMA to identify the best hits.
- **BWA Alignment**: Aligns reads to reference alleles to confirm results and identify insertions.
- **Consensus Sequence Extraction**: Extracts consensus sequences and identifies insertions from mapped reads.
- **Stx Profile Calculation**: Calculates Shiga toxin (stx) profiles for each sample.
- **Comprehensive Reporting**: Generates detailed reports, including allele matches, percent identity, and insertion details.

---

## **Installation**

### **Step 1: Install Conda**
Ensure you have Conda installed. If not, download and install Miniconda or Anaconda from [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html).

### **Step 2: Create the Conda Environment**
Create a Conda environment named `stec_kma` and install the required dependencies:
```bash
conda create -n stec_kma olcbioinformatics::stec_kma
conda activate stec_kma
```

---

## **Usage**
Run the `stec_kma.py` script with the required arguments:

```bash
python src/stec_kma.py -s <sequence_path> -d <database_path> -r <report_path> -t <threads> -ID <identity> -c <min_coverage>
```

### **Arguments**
- `-s, --sequence_path`: Path to the folder containing sequencing reads (required).
- `-d, --database_path`: **Name and path of the indexed KMA database** (required). This is the value provided to the `-o` argument when running the `kma index` command. **Ensure that the database has been processed with `kma index` before running this tool.**
- `-r, --report_path`: Path to the folder where reports will be written (optional, defaults to `sequence_path/reports`).
- `-t, --threads`: Number of threads to use (default: number of CPU cores).
- `-ID, --identity`: Minimum identity percentage for KMA hits (default: 90%).
- `-c, --min_coverage`: Minimum fraction of reads required to call a consensus insertion (default: 0.7).
- `--verbosity`: Logging verbosity level (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`; default: `INFO`).

---

## **Pipeline Overview**
The `stec_kma.py` script orchestrates the following steps:

1. **Locate Samples**:
   - Organizes FASTQ files into a dictionary and creates subdirectories for each sample.

2. **Bait Reads**:
   - Uses BBMap to bait reads from the input FASTQ files based on the reference allele database.

3. **Reverse Bait Targets**:
   - Baits targets from the allele database using the previously baited reads.

4. **Index Baited Sequences**:
   - Indexes the baited allele sequences using KMA.

5. **Map Reads with KMA**:
   - Maps reads to the indexed allele database using KMA.

6. **Extract Allele Sequences**:
   - Extracts allele sequences corresponding to the best hits from the KMA reports.

7. **Map Reads with BWA**:
   - Aligns reads to reference alleles using BWA to confirm results and identify insertions.

8. **Extract Consensus Sequences**:
   - Extracts consensus sequences and identifies insertions from the BWA output.

9. **Calculate Stx Profiles**:
   - Calculates Shiga toxin (stx) profiles for each sample.

10. **Generate Reports**:
    - Writes detailed reports summarizing the results.

---

## **Output**
The tool generates the following outputs:
1. **Reports**:
   - A tab-delimited report (`stec_kma_report.tsv`) summarizing:
     - Sample name
     - Best allele match
     - Percent identity
     - Stx profiles
     - Notes (e.g., insertions, truncations, or internal stop codons)

2. **FASTA Files**:
   - Nucleotide and amino acid sequences for alleles with insertions.

3. **Intermediate Files**:
   - Baited reads, mapped reads, and alignment files for further analysis.

---

## **Example**
Before running `stec_kma.py`, ensure that `kma index` has been run on your database. For example:
```bash
kma index -i /path/to/allele_db.fasta -o /path/to/indexed_db
```

Then, run the tool:
```bash
python src/stec_kma.py \
    -s /path/to/sequence_data \
    -d /path/to/indexed_db \
    -r /path/to/output_reports \
    -t 8 \
    -ID 95 \
    -c 0.8 \
    --verbosity DEBUG
```

---

## **Logging**
The tool supports configurable logging levels (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). Logs include detailed information about each step of the pipeline, including system commands and outputs.

---

## **License**
This project is licensed under the MIT License. See the `LICENSE` file for details.

---

## **Contact**
For questions or issues, please contact:
- **Author**: Adam Koziol
- **Email**: adam.koziol@inspection.gc.ca