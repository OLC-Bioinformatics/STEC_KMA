#!/usr/bin/env python

"""
Methods for the STEC_KMA package
"""

# Standard imports
from collections import defaultdict
from glob import glob
import gzip
import logging
import os
import re
import shlex
import statistics
import subprocess
from typing import (
    Dict,
    List,
    Union
)

# Third-party imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam

__author__ = 'adamkoziol'


def _setup_logging(
    # *,  Enforce keyword arguments
    verbosity: str
) -> None:
    """
    Set up and initialise logging

    Args:
        verbosity: The logging level
    """
    # Create a logging format
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    # Set the logging level
    log_level = getattr(logging, verbosity)
    # Set up the logging
    logging.basicConfig(level=log_level, format=log_format)


def _find_samples(
    sequence_path: str
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Organizes fastq files into a dictionary and creates sub-folders for each
    sample.

    Args:
        fastq_dir (str): Path to the directory containing fastq files.

    Returns:
        dict: A dictionary where keys are sample names and values are
        dictionaries containing lists of fastq file paths and the output path
        for the sample.
    """
    # Create a dictionary to store the samples
    sample_dict = {}

    # Glob all fastq files
    fastq_files = glob(os.path.join(sequence_path, "*.fastq.gz"))

    for fastq_file in sorted(fastq_files):
        # Extract sample name from the fastq file name
        file_name = os.path.basename(fastq_file)
        print(file_name)

        # Find the sample name using regex
        match = re.search(r"(\d{4}-[A-Za-z]{3,5}-\d{4})", file_name)

        if match:
            sample_name = match.group(1)
        else:
            logging.warning(
                "Warning: Could not extract sample name from: %s", file_name
            )
            continue

        # Create output directory for the sample
        output_path = os.path.join(sequence_path, sample_name)
        os.makedirs(output_path, exist_ok=True)

        # Create the symlinks in the output directory
        file_name = os.path.basename(fastq_file)
        link_path = os.path.join(output_path, file_name)
        if not os.path.exists(link_path):
            os.symlink(fastq_file, link_path)

        # Update the sample dictionary
        if sample_name not in sample_dict:
            sample_dict[sample_name] = {
                "files": [],
                "output_path": output_path
            }

        # Add the file to the sample dictionary
        sample_dict[sample_name]["files"].append(link_path)

    return sample_dict


def _bait_reads_bbmap(
    *,  # Enforce keyword arguments
    db_fasta: str,
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]]
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Use BBMap to bait reads with the combined allele database
    """
    for sample, sequence_dict in sample_dict.items():

        # Output path for the sample
        baited_reads = os.path.join(
            sequence_dict["output_path"],
            f'{sample}_baited.fastq.gz'
        )

        # Update the sample dictionary
        sequence_dict['baited_reads'] = baited_reads

        # Create the system call
        bbduk_call = (
            f'bbduk.sh in={sequence_dict['files'][0]} '
            f'in2={sequence_dict['files'][1]} '
            f'outm={baited_reads} '
            f'ref={db_fasta} '
            '--fixjunk'
        )

        logging.info('System call: %s', bbduk_call)

        # Run the command if the outputs do not already exist
        if not os.path.exists(baited_reads):
            output = run_command(
                command=bbduk_call, split=False
            )
            logging.debug('BBMap output: %s', output.stdout)

    return sample_dict


def _reverse_bait_targets_bbmap(
    *,  # Enforce keyword arguments
    db_fasta: str,
    identity: int,
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]]
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Use BBMap to bait targets from the combined allele database with the
    previously baited reads

    Args:
        db_fasta: Path to the gene-specific FASTA file
        identity: Minimum identity percentage
        sample_dict: Dictionary containing the sample information

    Returns:
        dict: The updated sample dictionary
    """
    for sample, sequence_dict in sample_dict.items():

        # Output path for the baited targets
        baited_targets = os.path.join(
            sequence_dict['output_path'],
            f'{sample}_baited_targets.fasta'
        )

        # Update the sample dictionary
        sequence_dict['baited_targets'] = baited_targets

        # Convert the identity to a fraction
        fraction = identity / 100

        # Create the system call
        bbduk_call = (
            f'bbduk.sh ref={sequence_dict['baited_reads']} '
            f'outm={baited_targets} '
            f'in={db_fasta} '
            f'mincovfraction={fraction} '
            '--maskmiddle=f '
            '--fixjunk'
        )

        logging.info('System call: %s', bbduk_call)

        # Run the command if the outputs do not already exist
        if not os.path.exists(baited_targets):
            output = run_command(
                command=bbduk_call, split=False
            )
            logging.debug('BBMap output: %s', output.stdout)

    return sample_dict


def _index_baited_sequences(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]],
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Use KMA to index the allele sequences

    Args:
        sample_dict: Dictionary containing the sample information

    Returns:
        dict: The updated sample dictionary
    """

    for _, sequence_dict in sample_dict.items():

        # Extract the baited targets
        baited_targets = sequence_dict['baited_targets']

        # Check if the baited targets file exists and is not empty
        if (
            not os.path.exists(baited_targets)
            or os.path.getsize(baited_targets) == 0
        ):
            logging.warning(
                'Warning: Baited targets file either does not exist or is '
                'empty: %s',
                baited_targets
            )

            # Set the KMA sample database to None
            sequence_dict['kma_sample_db'] = None

            # Skip to the next sample
            continue

        # Remove the extension from the allele name
        allele_db = baited_targets.split('.')[0] + '_db'

        # Update the sample dictionary
        sequence_dict['kma_sample_db'] = allele_db

        # Construct the KMA index command
        kma_index_cmd = (
            f'kma index -i {baited_targets} -o {allele_db}'
        )

        logging.debug('KMA index command: %s', kma_index_cmd)

        # Construct the report file name
        report_file = allele_db + '.name'

        # Run the analyses
        if not os.path.exists(report_file):
            output = run_command(
                command=kma_index_cmd,
                split=False
            )
            logging.debug('KMA index output: %s', output.stdout)

    return sample_dict


def _map_allele_sequences_kma(
    *,  # Enforce keyword arguments
    identity: int,
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]],
    threads: int
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Run KMA in parallel to map the reads to the sample-specific mapped allele
    databases

    Args:
        identity: Minimum identity percentage for KMA hits
        sample_dict: Dictionary containing the sample information
        threads: Number of threads to use

    Returns:
        dict: The updated sample dictionary
    """
    for _, sequence_dict in sample_dict.items():

        # Ensure that the sample has a KMA database
        if sequence_dict['kma_sample_db'] is None:
            continue

        # Extract the indexed allele database
        allele_db = sequence_dict['kma_sample_db']

        # Extract the output path
        output_path = os.path.splitext(
            sequence_dict['baited_targets']
        )[0]

        # Construct the KMA command
        kma_cmd = (
            f'kma -int {sequence_dict['baited_reads']} '
            f'-o {output_path} '
            f'-t_db {allele_db} '
            f'-ID {identity} '  # Identity threshold
            f'-t {threads} '     # Number of threads
            '-mem_mode '        # Memory mode
            '-a '               # Output all matches
        )

        # Construct the report file name
        report_file = output_path + '.res'

        # Update the sample dictionary
        sequence_dict['kma_report_file'] = report_file

        # Run the analyses
        if not os.path.exists(report_file):
            logging.debug('KMA command: %s', kma_cmd)
            output = run_command(
                command=kma_cmd,
                split=False
            )

            logging.debug('KMA output: %s', output.stdout)

    return sample_dict


def _parse_allele_reports(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]]
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Parse the KMA reports to find the best hits

    Args:
        sample_dict: Dictionary containing the sample information

    Returns:
        dict: The updated sample dictionary
    """
    for sample, sequence_dict in sample_dict.items():

        # Ensure that the sample has a KMA database
        if sequence_dict['kma_sample_db'] is None:
            continue

        # Create a dictionaries to store the best hits, best scores, and
        # best hit files
        best_hits = {}
        best_scores = {}
        best_hit_files = {}

        # Extract the report file
        report_file = sequence_dict['kma_report_file']

        # Read the report file
        with open(report_file, 'r', encoding='utf-8') as report_fh:
            for line in report_fh:

                # Skip the header
                if line.startswith('#'):
                    continue

                try:
                    # Extract the reference allele and percentage identity
                    # Template	Score	Expected	Template_length
                    # Template_Identity	Template_Coverage	Query_Identity
                    # Query_Coverage	Depth	q_value	p_value
                    reference, _, _, _, t_identity, \
                        _, _, _, _, _, _ = line.strip().split('\t')

                    # Extract the gene name from the report file
                    gene_name = _find_gene_name(allele=reference)

                    # Initialise the dictionaries for the gene
                    if gene_name not in best_hits:
                        best_hits[gene_name] = {}
                    if gene_name not in best_scores:
                        best_scores[gene_name] = 0
                    if gene_name not in best_hit_files:
                        best_hit_files[gene_name] = ''

                    # Update the best hit if the score is higher
                    if float(t_identity) > float(best_scores[gene_name]):
                        best_hits[gene_name] = reference
                        best_scores[gene_name] = float(t_identity)
                        best_hit_files[gene_name] = report_file
                except ValueError:
                    pass

        # Update the sample dictionary
        sequence_dict['best_hits'] = best_hits
        sequence_dict['best_scores'] = best_scores
        sequence_dict['best_hit_files'] = best_hit_files

        # Log the best hits
        logging.debug('Best hits for sample %s:', sample)
        for gene_name, best_hit in best_hits.items():
            logging.debug(
                '\t%s: %s with score %s',
                gene_name, best_hit, best_scores[gene_name]
            )

    return sample_dict


def _extract_allele_sequence(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]]
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Extract the allele sequences from the baited targets FASTA file
    that correspond to the best hits from the KMA report files

    Args:
        sample_dict: Dictionary containing the sample information

    Returns:
        dict: The updated sample dictionary
    """

    for _, sequence_dict in sample_dict.items():

        # Ensure that the sample has a KMA database
        if sequence_dict['kma_sample_db'] is None:
            continue

        # Extract the baited targets
        baited_targets = sequence_dict['baited_targets']

        # Create a dictionary to store the allele sequences
        allele_sequences = {}

        # Read the baited targets FASTA file
        for record in SeqIO.parse(baited_targets, 'fasta'):
            # Extract the allele name from the record
            allele_name = record.id

            # Check if the allele name is in the best hits
            if allele_name in sequence_dict['best_hits'].values():

                # Write the allele sequence to a FASTA file
                output_file = _write_allele_sequences(
                    record=record,
                    sequence_dict=sequence_dict
                )
                allele_sequences[allele_name] = output_file

        # Update the sample dictionary
        sequence_dict['allele_sequences'] = allele_sequences

    return sample_dict


def _write_allele_sequences(
    *,  # Enforce keyword arguments
    record: SeqIO.SeqRecord,
    sequence_dict: Dict[str, Dict[str, Union[List[str], str]]]
) -> str:
    """
    Write the allele sequences to a FASTA file
    Args:
        record: The SeqIO record to write
        sequence_dict: The sample dictionary

    Returns:
        str: The path to the output file
    """
    # Extract the gene name from the allele name
    gene_name = _find_gene_name(allele=record.id)

    # Set the name of the output directory
    output_dir = os.path.join(
        sequence_dict['output_path'],
        gene_name
    )

    # Create the output directory if it does not exist
    os.makedirs(output_dir, exist_ok=True)

    # Ensure that the allele name does not contain any special characters or
    # punctuation
    allele_name = re.sub(r'[^a-zA-Z0-9_]', '_', record.id)

    # Create the output file path
    output_file = os.path.join(
        output_dir,
        f'{allele_name}.fasta'
    )

    # Write the record to the output file
    SeqIO.write(record, output_file, 'fasta')

    return output_file


def _extract_kma_mapped_reads(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]]
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Extract the mapped reads from the KMA report files

    Args:
        sample_dict: Dictionary containing the sample information

    Returns:
        dict: The updated sample dictionary
    """
    for _, sequence_dict in sample_dict.items():

        # Ensure that the sample has a KMA database
        if sequence_dict['kma_sample_db'] is None:
            continue

        # Extract the KMA report file
        kma_report_file = sequence_dict['kma_report_file']

        # Set the name of the .frag.gz file
        frag_file = os.path.splitext(kma_report_file)[0] + '.frag.gz'

        # Initialise a dictionary to store the sample reads
        sequence_dict['mapped_reads'] = {}

        # Read in the frag file
        with gzip.open(frag_file, 'rt') as frag_fh:
            for line in frag_fh:

                # Split the line into fields
                info = line.strip().split('\t')

                # Extract the allele name and read name
                allele = info[-2]
                read_name = info[-1]

                # Check if the allele is in the best hits
                if allele in sequence_dict['best_hits'].values():
                    gene = _find_gene_name(allele=allele)

                    # Update the sample dictionary
                    if gene not in sequence_dict['mapped_reads']:
                        sequence_dict['mapped_reads'][gene] = {}

                    if allele not in sequence_dict['mapped_reads'][gene]:
                        sequence_dict['mapped_reads'][gene][allele] = []

                    # Append the read name to the list
                    sequence_dict['mapped_reads'][gene][allele].append(
                        read_name
                    )

    return sample_dict


def _write_read_names_file(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]],
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Write the KMA-mapped read names to a file

    Args:
        sample_dict: Dictionary containing the sample information

    Returns:
        dict: The updated sample dictionary
    """
    for _, sequence_dict in sample_dict.items():

        # Ensure that the sample has a KMA database
        if sequence_dict['kma_sample_db'] is None:
            continue

        # Initialise a dictionary to store the mapped read files
        sequence_dict['mapped_read_files'] = {}

        # Create a dictionary to store the output directories
        sequence_dict['gene_output_dir'] = {}

        for gene, allele_dict in sequence_dict['mapped_reads'].items():

            # Set the name of the output directory
            sequence_dict['gene_output_dir'][gene] = os.path.join(
                sequence_dict['output_path'],
                gene
            )

            # Ensure the output directory exists
            os.makedirs(sequence_dict['gene_output_dir'][gene], exist_ok=True)

            # Update the sample dictionary
            if gene not in sequence_dict['mapped_read_files']:
                sequence_dict['mapped_read_files'][gene] = {}

            for allele, read_list in allele_dict.items():

                # Set the name of the output file
                output_file = os.path.join(
                    sequence_dict['gene_output_dir'][gene],
                    f'{gene}_{allele}.names'
                )

                # Write the read names to the output file
                with open(output_file, 'w', encoding='utf-8') as out_fh:
                    for read_name in read_list:
                        out_fh.write(f'{read_name}\n')

                # Update the sample dictionary
                sequence_dict['mapped_read_files'][gene][allele] = output_file

    return sample_dict


def _extract_reads_by_name(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]],
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Use filterbyname to extract the KMA-mapped reads from the FASTQ files

    Args:
        sample_dict: Dictionary containing the sample information

    Returns:
        dict: The updated sample dictionary
    """
    for _, sequence_dict in sample_dict.items():

        # Ensure that the sample has a KMA database
        if sequence_dict['kma_sample_db'] is None:
            continue

        # Extract the baited reads
        baited_reads = sequence_dict['baited_reads']

        # Create a dictionary to store the extracted reads
        extracted_reads = {}

        for gene, allele_dict in sequence_dict['mapped_reads'].items():

            # Update the extracted reads dictionary
            if gene not in extracted_reads:
                extracted_reads[gene] = {}

            for allele in allele_dict:

                # Extract the read names file
                read_file = sequence_dict['mapped_read_files'][gene][allele]

                # Set the name of the output file
                output_file = os.path.join(
                    sequence_dict['gene_output_dir'][gene],
                    f'{gene}_{allele}.fastq.gz'
                )

                # Update the sample dictionary
                extracted_reads[gene][allele] = output_file

                # Construct the filterbyname command
                filterbyname_cmd = (
                    f'filterbyname.sh in={baited_reads} '
                    f'out={output_file} '
                    f'include=t names={read_file}'
                )

                logging.debug('Filterbyname command: %s', filterbyname_cmd)

                # Run the command if the output does not already exist
                if not os.path.exists(output_file):
                    output = run_command(
                        command=filterbyname_cmd,
                        split=False
                    )
                    logging.debug('Filterbyname output: %s', output.stdout)

        # Update the sample dictionary
        sequence_dict['extracted_reads'] = extracted_reads

    return sample_dict


def _find_gene_name(
    *,  # Enforce keyword arguments
    allele: str
) -> str:
    """
    Extract the gene name from the allele name

    Args:
        alleles: string of an allele name

    Returns:
        string: gene name
    """
    return os.path.basename(allele).split('_')[0]


def _map_reads_bwa(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]],
    threads: int
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Use BWA to map the reads to the reference allele

    Args:
        sample_dict: Dictionary containing the sample information
        threads: Number of threads to use

    Returns:
        dict: The updated sample dictionary
    """
    for _, sequence_dict in sample_dict.items():

        # Ensure that the sample has a KMA database
        if sequence_dict['kma_sample_db'] is None:
            continue

        # Initialise a dictionary to store the BWA output
        sequence_dict['bwa_output'] = {}

        for __, allele_dict in sequence_dict['extracted_reads'].items():

            for allele in allele_dict:

                # Extract the extracted reads database
                allele_db = sequence_dict['allele_sequences'][allele]

                # Extract the gene name from the allele name
                gene = _find_gene_name(allele=allele)

                # Extract the parent directory of the reference allele
                reference_allele_dir = sequence_dict['gene_output_dir'][gene]

                # Index the reference allele if it was not already done
                if not os.path.exists(allele_db + '.bwt'):
                    _index_reference_allele(reference_allele=allele_db)

                # Set the name of the output file
                output_bam = os.path.join(
                    reference_allele_dir,
                    f'{allele}.bam'
                )

                # Update the sample dictionaries
                sequence_dict['bwa_output'][allele] = output_bam

                # Construct BWA command
                bwa_cmd = [
                    'bwa', 'mem', '-t', str(threads), allele_db,
                    sequence_dict['baited_reads']
                ]

                # Construct Samtools view command
                samtools_view_cmd = ['samtools', 'view', '-Sb', '-']

                # Construct Samtools sort command
                samtools_sort_cmd = ['samtools', 'sort', '-o', output_bam]

                logging.debug('BWA command: %s', ' '.join(bwa_cmd))
                logging.debug(
                    'Samtools view command: %s', ' '.join(samtools_view_cmd)
                )
                logging.debug(
                    'Samtools sort command: %s', ' '.join(samtools_sort_cmd)
                )

                # Run the command if the output does not already exist
                if not os.path.exists(output_bam):
                    try:
                        # Execute BWA
                        bwa_process = subprocess.Popen(
                            bwa_cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE
                        )

                        # Execute Samtools view
                        samtools_view_process = subprocess.Popen(
                            samtools_view_cmd,
                            stdin=bwa_process.stdout,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE
                        )
                        bwa_process.stdout.close()  # Allow BWA to exit

                        # Execute Samtools sort
                        samtools_sort_process = subprocess.Popen(
                            samtools_sort_cmd,
                            stdin=samtools_view_process.stdout,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE
                        )
                        # Allow samtools view to exit
                        samtools_view_process.stdout.close()

                        # Get the output and error messages
                        _, stderr = samtools_sort_process.communicate()

                        if stderr:
                            logging.error(
                                "Error: %s", stderr.decode()
                            )
                        else:
                            logging.debug(
                                'BWA/Samtools pipeline completed successfully'
                            )

                        # Index the BAM file
                        _index_bam_file(bam_file=output_bam)

                    except FileNotFoundError as exc:
                        logging.error(
                            "Error: One of the required files not found: %s",
                            exc
                        )
                        raise SystemExit from exc
                    except Exception as exc:
                        logging.error("An error occurred: %s", exc)
                        raise SystemExit from exc

    return sample_dict


def _index_reference_allele(
    *,  # Enforce keyword arguments
    reference_allele: str
) -> None:
    """
    Run BWA index on the reference allele

    Args:
        reference_allele: Path to the reference allele
    """
    # Create the system call
    index_cmd = f'bwa index {reference_allele}'

    logging.debug('System call: %s', index_cmd)

    try:
        # Run the command
        output = run_command(
            command=index_cmd,
            split=False
        )
        logging.debug('BWA index output: %s', output.stdout)
    except subprocess.CalledProcessError as exc:
        logging.error("Command failed: %s", exc)
        raise


def _index_bam_file(
    *,  # Enforce keyword arguments
    bam_file: str
) -> None:
    """
    Index the sorted BAM file using samtools index

    Args:
        bam_file: Path to the sorted BAM file
    """
    # Create the samtools index command
    index_cmd = ['samtools', 'index', bam_file]

    logging.info('Samtools index command: %s', ' '.join(index_cmd))

    # Only run the command if the index file does not already exist
    if os.path.exists(bam_file + '.bai'):
        return

    try:
        # Run the samtools index command
        output = subprocess.run(
            index_cmd, capture_output=True, text=True, check=True
        )
        logging.debug('Samtools index output: %s', output.stdout)
    except subprocess.CalledProcessError as exc:
        logging.error("Samtools index command failed: %s", exc)
        raise


def _extract_consensus_insertions_and_sequence(
    *,  # Enforce keyword arguments,
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]],
    min_coverage: float = 0.5,
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Iterate through the BWA output files to extract the consensus insertions
    and sequences

    Args:
        sample_dict: Dictionary containing the sample information
         min_coverage: Minimum *fraction* of reads required to call a consensus
        insertion. For example, 0.5 means the insertion must be present in at
        least 50% of reads covering that position.

    Returns:
        dict: The updated sample dictionary
    """
    for _, sequence_dict in sample_dict.items():

        # Ensure that the sample has a KMA database
        if sequence_dict['kma_sample_db'] is None:
            continue

        # Initialise a dictionary to store outputs
        sequence_dict['insertions'] = {}
        sequence_dict['consensus'] = {}
        sequence_dict['aa_sequence'] = {}
        sequence_dict['aa_sequence_qa'] = {}

        for allele, bam_file in sequence_dict['bwa_output'].items():

            # Calculate the consensus internal insertions
            insertions, consensus_sequence = \
                _find_consensus_internal_insertions_and_sequence(
                    bam_file=bam_file,
                    min_coverage=min_coverage,
                )

            # Calculate the amino acid sequence
            aa_sequence = _translate_allele_sequence(
                nt_sequence=consensus_sequence
            )

            # Check if the amino acid sequence is valid
            aa_sequence_qa = _detect_truncated_sequences(
                protein_sequence=aa_sequence,
            )

            # Update the sample dictionaries
            sequence_dict['insertions'][allele] = insertions
            sequence_dict['consensus'][allele] = consensus_sequence
            sequence_dict['aa_sequence'][allele] = aa_sequence
            sequence_dict['aa_sequence_qa'][allele] = aa_sequence_qa

    return sample_dict


def _find_consensus_internal_insertions_and_sequence(
    bam_file: str,
    min_coverage: float = 0.5,
) -> tuple[list[tuple[int, int]], str]:
    """
    Identifies consensus internal insertions (from CINS and CSCLIP) from a
    sorted BAM file and calculates the consensus sequence.

    Args:
        bam_file: Path to the sorted BAM file.
        min_coverage: Minimum *fraction* of reads required to call a consensus
        insertion. For example, 0.5 means the insertion must be present in at
        least 50% of reads covering that position.

    Returns:
        A tuple containing:
        - A list of tuples, where each tuple contains the reference position
        of the insertion and the median length of the insertions at that
        position. **Positions are reported in 1-based coordinates.**
        - A string representing the consensus sequence.
    """
    # Open the BAM file
    sam_file = pysam.AlignmentFile(bam_file, "rb")

    # Get the reference length and name
    ref_len = sam_file.lengths[0]
    ref_name = sam_file.references[0]

    # Calculate the depth at each position
    total_reads_at_position = _calculate_depth_per_position(
        bam_file=bam_file,
        ref_name=ref_name,
        ref_len=ref_len
    )

    # Store lengths of insertions at each position
    insertion_positions = defaultdict(list)

    # Store base counts for consensus sequence calculation
    base_counts = defaultdict(lambda: defaultdict(int))

    for read in sam_file.fetch(ref_name, 0, ref_len):

        query_sequence = read.query_sequence

        # Extract CIGAR operations
        cigar = read.cigartuples

        # Initialize positions
        query_pos = 0
        ref_pos = read.reference_start

        # Iterate through CIGAR operations
        for operation, length in cigar:

            # Match or equal
            if operation in [pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF]:
                for _ in range(length):
                    if ref_pos < ref_len:
                        base_counts[ref_pos][query_sequence[query_pos]] += 1
                    ref_pos += 1
                    query_pos += 1

            # Deletion
            elif operation == pysam.CDEL:
                ref_pos += length

            # Insertion or soft clip
            elif operation in [pysam.CINS, pysam.CSOFT_CLIP]:

                # Determine if the insertion is at the 5' or 3' end. Allow a
                # 10bp buffer for the 5' and 3' ends
                five_prime_clip = ref_pos <= 10
                three_prime_clip = ref_pos >= ref_len - 10

                # Extract internal soft clips
                if not five_prime_clip and not three_prime_clip:
                    # Store the insertion length at the next reference position
                    # (1-based coordinates)
                    insertion_positions[ref_pos].append(length)

                query_pos += length

            # Hard clip
            elif operation == pysam.CREF_SKIP:
                ref_pos += length

            # Pad
            elif operation in [pysam.CHARD_CLIP, pysam.CPAD]:
                pass
            else:
                raise ValueError(f"Unexpected CIGAR operation: {operation}")

    # Close the BAM file
    sam_file.close()

    # Filter insertions based on minimum coverage
    consensus_insertions = []

    for ref_pos, lengths in sorted(insertion_positions.items()):

        # Calculate the total number of reads at this position
        num_reads = len(lengths)

        # Extract the total number of reads at this position
        total_reads = total_reads_at_position[ref_pos]

        # Check if the position is covered by enough reads to pass the
        # minimum coverage threshold
        if (
            ref_pos in total_reads_at_position
            and (num_reads / total_reads) >= min_coverage
        ):

            # Calculate median length of insertions
            median_length = int(statistics.median(lengths))

            # Append the reference position and median length to the list
            consensus_insertions.append((ref_pos, median_length))

    # Build the consensus sequence
    consensus_sequence = []
    for pos in range(ref_len):
        if pos in base_counts:
            # Get the most frequent base at this position
            most_frequent_base = max(
                base_counts[pos], key=base_counts[pos].get
            )
            consensus_sequence.append(most_frequent_base)
        else:
            # If no base is present, use a gap
            consensus_sequence.append('-')

    # Join the consensus sequence into a string
    consensus_sequence = ''.join(consensus_sequence)

    return consensus_insertions, consensus_sequence


def _calculate_depth_per_position(
    bam_file: str,
    ref_name: str,
    ref_len: int
) -> Dict[int, int]:
    """
    Calculate the depth (number of reads) at each position in the reference.

    Args:
        bam_file: Path to the sorted BAM file.
        ref_name: Name of the reference sequence.
        ref_len: Length of the reference sequence.

    Returns:
        A dictionary where keys are reference positions (0-based) and values
        are the depth (number of reads covering that position).
    """
    # Open the BAM file
    sam_file = pysam.AlignmentFile(bam_file, "rb")

    # Initialize a dictionary to store the depth at each position
    depth_per_position = defaultdict(int)

    # Iterate through the pileup columns for the reference
    for pileup_column in sam_file.pileup(ref_name, 0, ref_len):
        # Get the reference position (0-based) and the depth
        ref_pos = pileup_column.reference_pos + 1
        depth = pileup_column.nsegments

        # Update the depth dictionary
        depth_per_position[ref_pos] = depth

    # Close the BAM file
    sam_file.close()

    return depth_per_position


def _translate_allele_sequence(
    *,  # Enforce keyword arguments
    nt_sequence: str
) -> str:
    """
    Translate the allele sequence to protein.

    Args:
        nt_sequence: Nucleotide sequence of the allele.

    Returns:
        str: The translated protein sequence.
    """
    # Ensure the sequence length is a multiple of three
    remainder = len(nt_sequence) % 3
    if remainder != 0:
        # Pad with trailing 'N' to make the length a multiple of three
        nt_sequence += 'N' * (3 - remainder)

    # Create a Seq object from the nucleotide sequence
    seq = Seq(nt_sequence)

    # Translate the sequence to protein
    protein_sequence = seq.translate()

    return str(protein_sequence)


def _detect_truncated_sequences(
    *,  # Enforce keyword arguments
    protein_sequence: str,
    min_length: int = 313
) -> str:
    """
    Detect truncated sequences based on the length of the protein sequence

    Args:
        protein_sequence: Protein sequence of the allele
        min_length: Minimum length of the protein sequence

    Returns:
        string describing the truncation/internal stop codon
    """
    # Check if the protein sequence is truncated
    if len(protein_sequence) < min_length:
        return f'Truncated: {len(protein_sequence)}'

    # Check for stop codons prior to the min length
    for i in range(0, min(len(protein_sequence), min_length)):
        if protein_sequence[i] == '*':
            return f'Internal stop codon at position: {i + 1}'

    return False


def _write_novel_allele_sequences(
    *,  # Enforce keyword arguments
    report_path: str,
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]],
) -> None:
    """
    Write the nucleotide and amino acid sequence of any alleles that are not
    a 100% match to the reference allele

    Args:
        report_path: Directory into which the report file is to be written
        sample_dict: Dictionary containing the sample information
    """
    for sample, sequence_dict in sample_dict.items():

        # Ensure that the sample has a KMA database
        if sequence_dict['kma_sample_db'] is None:
            continue

        # Extract the allele identities
        for gene, identity in sequence_dict['best_scores'].items():

            # Extract the best hit allele
            allele = sequence_dict['best_hits'][gene]

            # Check if the identity is less than 100%
            if identity == 100 or identity < 90:
                continue

            # Extract the allele sequence
            allele_sequence = sequence_dict['consensus'][allele]
            aa_allele_seq = sequence_dict['aa_sequence'][allele]

            # Extract the gene name from the allele name
            gene = _find_gene_name(allele=allele)

            # Define the output file names
            nucleotide_file = os.path.join(
                report_path,
                f'{sample}_{gene}_nucleotide.fasta'
            )
            amino_acid_file = os.path.join(
                report_path,
                f'{sample}_{allele}_amino_acid.fasta'
            )

            # Create SeqRecords of the allele sequences
            nucleotide_record = SeqRecord(
                Seq(allele_sequence),
                id=f'{sample}_{gene}',
                description=''
            )
            amino_acid_record = SeqRecord(
                Seq(aa_allele_seq),
                id=f'{sample}_{gene}',
                description=''
            )

            # Write the nucleotide sequence to a FASTA file
            with open(nucleotide_file, 'w', encoding='utf-8') as out_fh:
                SeqIO.write(nucleotide_record, out_fh, 'fasta')

            # Write the amino acid sequence to a FASTA file
            with open(amino_acid_file, 'w', encoding='utf-8') as out_fh:
                SeqIO.write(amino_acid_record, out_fh, 'fasta')


def _calculate_stx_profile(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]]
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Calculate the stx profile for each sample

    Args:
        sample_dict: Dictionary containing the sample information

    Returns:
        dict: The updated sample dictionary
    """

    for _, sequence_dict in sample_dict.items():

        # Ensure that the sample has a KMA database
        if sequence_dict['kma_sample_db'] is None:
            continue

        # Create a dictionary to store the stx profiles
        stx_profiles = {}

        # Extract the best hits
        best_hits = sequence_dict['best_hits']

        # Iterate through the best hits
        for gene, allele in best_hits.items():

            # Extract the best score
            best_score = sequence_dict['best_scores'][gene]

            # Check if the best score is above the threshold
            if best_score < 90:
                continue

            # Update the stx profiles dictionary
            if gene not in stx_profiles:
                stx_profiles[gene] = []

            # Append the allele to the list
            stx_profiles[gene].append(allele.lower())

        # Update the sample dictionary
        sequence_dict['stx_profiles'] = sorted(stx_profiles)

    return sample_dict


def _write_report(
    *,  # Enforce keyword arguments
    report_path: str,
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]]
):
    """
    Write the report to file

    Args:
        report_path: Directory into which the report file is to be written
        sample_dict: Dictionary containing the sample information
    """

    # Set the name of the combined report file
    combined_report_file = os.path.join(
        report_path,
        'stec_kma_report.tsv'
    )

    # Set the headers for the report file
    headers = [
        'Sample',
        'Allele',
        'PercentIdentity',
        'ExpectedProfile',
        'Notes'
    ]

    # Write the report file
    with open(combined_report_file, 'w', encoding='utf-8') as out_fh:

        # Write the header
        out_fh.write('\t'.join(headers) + '\n')
        print('\t'.join(headers))

        for sample, sequence_dict in sample_dict.items():

            # Ensure that the sample has a KMA database
            if sequence_dict['kma_sample_db'] is None:
                out_fh.write(f'{sample}\tND\n')
                print(f'{sample}\tND')
                continue

            # Extract the best hits
            best_hits = sequence_dict['best_hits']

            # Extract the stx profiles
            stx_profiles = ';'.join(sequence_dict['stx_profiles'])

            # Iterate through the best hits
            for gene, allele in best_hits.items():

                # Extract the best score
                best_score = sequence_dict['best_scores'][gene]

                # Check if the best score is above the threshold
                if best_score < 90:
                    continue

                # Extract the insertions
                insertions = sequence_dict['insertions'][allele]

                # Format the insertions
                insertions = ';'.join(
                    [
                        f'Insertion:Position({pos}):MedianLength({length})'
                        for pos, length in insertions
                    ]
                )

                # Extract the consensus sequence QA
                qa = sequence_dict['aa_sequence_qa'][allele]

                # Check if the sequence is truncated or has an internal stop
                # codon
                if qa:
                    if insertions:
                        insertions += ';'
                    insertions += qa

                # Write the report line
                out_fh.write(
                    f'{sample}\t{allele}\t{best_score}\t'
                    f'{stx_profiles}\t{insertions}\n'
                )
                print(
                    f'{sample}\t{allele}\t{best_score}\t'
                    f'{stx_profiles}\t{insertions}')


def run_command(
    *,  # Enforce keyword arguments
    command: str,
    capture_output: bool = True,
    split: bool = True,
    text: bool = True
) -> subprocess.CompletedProcess:
    """
    Runs a shell command using subprocess.

    Args:
        command: The command to run as a string.
        capture_output: Whether to capture stdout and stderr.
        split: Whether to split the command into a list.
        text: Whether to decode stdout and stderr as text.

    Returns:
        A subprocess.CompletedProcess instance.

    Raises:
        subprocess.CalledProcessError: If the command returns a non-zero exit
        code.
    """
    try:
        if split:
            command = shlex.split(command)
            result = subprocess.run(
                command,
                capture_output=capture_output,
                text=text,
                check=True
            )
        else:
            result = subprocess.run(
                command,
                capture_output=capture_output,
                text=text,
                check=True,
                shell=True
            )
        return result
    except subprocess.CalledProcessError as exc:
        logging.error("Command failed: %s",  exc)
        raise


def _tilde_expand(*, path: str) -> str:
    """
    Expand the tilde (~) in the supplied path to the full user directory path.

    :param path: The path that may contain a tilde
    :return: The expanded absolute path

    Example usage:
    >>> expanded_path = tilde_expand('~/my_dir')
    >>> print(expanded_path)  # Output: '/home/user/my_dir'
    """
    # Check if the path starts with a tilde (~)
    if path.startswith('~'):
        # Expand the tilde to the full user directory path
        return_path = os.path.abspath(
            os.path.expanduser(
                os.path.join(path)
            )
        )
    else:
        # Return the absolute path if no tilde is present
        return_path = os.path.abspath(
            os.path.join(path)
        )

    return return_path