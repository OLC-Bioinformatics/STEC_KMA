#!/usr/bin/env python3

"""
This is the main module of the STEC_KMA package. It runs KMA with the combined
allele database, finds the best hits, extracts FASTQ reads, and runs BWA to
align the reads to the reference allele to confirm the outputs, and also return
any alleles with insertions
"""

# Standard imports
from argparse import ArgumentParser
from collections import defaultdict
from concurrent.futures import as_completed, ThreadPoolExecutor
from glob import glob
import logging
import os
import re
import shlex
import subprocess
from typing import (
    Dict,
    List,
    Union
)

# Third-party imports
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pysam

# Local imports
from src.version import __version__

__author__ = 'adamkoziol'


def setup_logging(
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


def main(
    *,  # Enforce keyword arguments
    sequence_path: str,
    database_path: str,
    threads: int,
    identity: int
) -> None:
    """
    Run KMA with the combined allele database, find the best hits, extract
    FASTQ reads, and run BWA to align the reads to the reference allele to
    confirm the outputs, and also return the full sequence of any alleles
    with insertions

    Args:
        sequence_path: Path of the folder containing sequencing reads
        database_path: Path of the folder containing the database
        threads: Number of threads to use
        identity: Minimum identity percentage for KMA hits
    """
    # Extract the database path
    db_path = os.path.dirname(database_path)

    # Create a list to store the reference files
    ref_files = []

    # Run glob on the database path to find the reference file
    for ext in ['fa', 'fna', 'fsa', 'fasta']:
        ref_files.extend(glob(os.path.join(db_path, f'*.{ext}')))
    db_fasta = ref_files[0]

    logging.info('Locating samples in the sequence path, %s', sequence_path)
    sample_dict = _find_samples(sequence_path=sequence_path)

    logging.info('Baiting reads with BBMap')
    sample_dict = _bait_reads_bbmap(
        db_fasta=db_fasta,
        sample_dict=sample_dict
    )

    logging.info('Baiting targets with BBMap')
    sample_dict = _reverse_bait_targets_bbmap(
        db_fasta=db_fasta,
        identity=identity,
        sample_dict=sample_dict
    )

    logging.info(
        'Splitting baited target database into individual alleles for '
        'KMA analyses'
    )
    sample_dict = _extract_sequences(
        sample_dict=sample_dict
        )
    sample_dict = _write_allele_sequences(sample_dict=sample_dict)

    logging.info('Indexing the allele sequences with KMA')
    sample_dict = _index_allele_sequences(
        sample_dict=sample_dict,
        threads=threads
    )

    logging.info('Mapping reads to the reference alleles with KMA')
    sample_dict = _map_allele_sequences_kma(
        identity=identity,
        sample_dict=sample_dict,
        threads=threads
    )
    sample_dict = _parse_allele_reports(sample_dict=sample_dict)

    logging.info('Mapping reads to the reference alleles with BWA')
    sample_dict = _map_reads_bwa(
        sample_dict=sample_dict,
        threads=threads
    )

    logging.info('Extracting consensus sequences from the BWA output')
    sample_dict = _consensus_sequence(sample_dict=sample_dict)


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


def _extract_sequences(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]]
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Extract the sample-specific matching alleles from the database file

    Args:
        sample_dict (dict): Dictionary containing the sample information

    Returns:
        dict: The updated sample dictionary
    """

    for _, sequence_dict in sample_dict.items():

        # Extract the baited targets
        baited_targets = sequence_dict['baited_targets']

        with open(baited_targets, 'r', encoding='utf-8') as fasta_file:
            # Create a list to store the sequences
            sequences = []

            # Iterate through the FASTA records
            records = SeqIO.parse(fasta_file, 'fasta')
            for record in records:
                sequences.append(record)

            # Update the sample dictionary
            sequence_dict['allele_sequences'] = sequences

    return sample_dict


def _write_allele_sequences(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[SeqRecord], str]]]
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Write the sample-specific matching alleles to individual FASTA files,
    named based on the record.id.

    Args:
        sample_dict: Dictionary containing the sample information, where
        'allele_sequences' contains a list of SeqRecord objects.

    Returns:
        dict: The updated sample dictionary
    """
    for _, sequence_dict in sample_dict.items():

        # Initialize lists to store the output file paths and names
        sequence_dict['allele_output_paths'] = []
        sequence_dict['allele_sequence_files'] = []

        # Extract the output path
        output_path = sequence_dict['output_path']

        # Iterate through each sequence record
        for record in sequence_dict['allele_sequences']:

            # Replace the '.' in the record ID with an underscore
            record_id = record.id.replace('.', '_')

            # Set the name and create the allele-specific output directory
            allele_output_path = os.path.join(
                output_path,
                record_id,
                record_id
            )

            # Create the directory if it does not already exist
            os.makedirs(allele_output_path, exist_ok=True)

            # Set the name of the output file based on the record ID
            output_fasta = os.path.join(
                allele_output_path,
                f'{record_id}.fasta'
            )

            # Update the sample dictionaries
            sequence_dict['allele_output_paths'].append(allele_output_path)
            sequence_dict['allele_sequence_files'].append(output_fasta)

            # Only write the sequence if the file does not already exist
            if not os.path.exists(output_fasta):
                # Write the sequence to the output file
                with open(output_fasta, 'w', encoding='utf-8') as fasta_file:
                    # Write a list containing the single record
                    SeqIO.write([record], fasta_file, 'fasta')

    return sample_dict


def _index_allele_sequences(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]],
    threads: int
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Use KMA to index the allele sequences

    Args:
        sample_dict: Dictionary containing the sample information
        threads: Number of threads to use

    Returns:
        dict: The updated sample dictionary
    """

    # Create a list to store the KMA commands
    kma_commands = []

    # Create a list to store the report files
    report_files = []

    for _, sequence_dict in sample_dict.items():
        # Initialize a list to store the allele databases
        sequence_dict['kma_sample_db'] = []

        # Iterate through the allele sequences
        for allele_fasta in sequence_dict['allele_sequence_files']:

            # Remove the extension from the allele name
            allele_db = allele_fasta.split('.')[0] + '_db'

            # Update the sample dictionary
            sequence_dict['kma_sample_db'].append(allele_db)

            # Construct the KMA index command
            kma_index_cmd = (
                f'kma index -i {allele_fasta} -o {allele_db}'
            )

            logging.debug('KMA index command: %s', kma_index_cmd)

            # Add the KMA command to the list
            kma_commands.append(kma_index_cmd)

            # Construct the report file name
            report_file = allele_db + '.name'

            # Add the report file to the list
            report_files.append(report_file)

    # Run the KMA commands in parallel
    _run_kma_commands_in_parallel(
        kma_commands=kma_commands,
        report_files=report_files,
        threads=threads
    )

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
    # Create a list to store the KMA commands
    kma_commands = []

    # Create a list to store the report files
    report_files = []

    for _, sequence_dict in sample_dict.items():

        sequence_dict['kma_report_files'] = []

        # Iterate through the allele databases
        for iterator, allele_db in enumerate(sequence_dict['kma_sample_db']):

            # Extract the output path
            output_path = sequence_dict['allele_output_paths'][iterator]

            # Construct the KMA command
            kma_cmd = (
                f'kma -int {sequence_dict['baited_reads']} '
                f'-o {output_path} '
                f'-t_db {allele_db} '
                f'-ID {identity} '  # Identity threshold
            )
            kma_commands.append(kma_cmd)

            # Construct the report file name
            report_file = output_path + '.res'
            report_files.append(report_file)

            # Update the sample dictionary
            sequence_dict['kma_report_files'].append(report_file)

    # Run the KMA commands in parallel
    _run_kma_commands_in_parallel(
        kma_commands=kma_commands,
        report_files=report_files,
        threads=threads
    )

    return sample_dict


def _run_kma_commands_in_parallel(
    *,  # Enforce keyword arguments
    kma_commands: List[str],
    report_files: List[str],
    threads: int
) -> None:
    """
    Run KMA commands in parallel

    Args:
        kma_commands (List[str]): List of KMA commands
        report_files (List[str]): List of report files
        threads: int
    """
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                run_command,
                command=cmd,
                split=False
            )
            for cmd, report in zip(kma_commands, report_files)
            if not os.path.isfile(report)
        ]
        for future in as_completed(futures):
            future.result()


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
    for _, sequence_dict in sample_dict.items():
        print(_)
        # Create a dictionaries to store the best hits, best scores, and
        # best hit files
        best_hits = {}
        best_scores = {}
        best_hit_files = {}
        allele_db = {}

        # Iterate through the report files
        for report_file in sequence_dict['kma_report_files']:

            # Extract the gene name from the report file
            gene_name = _find_gene_name(allele=report_file)

            # Initialise the dictionaries for the gene
            if gene_name not in best_hits:
                best_hits[gene_name] = {}
            if gene_name not in best_scores:
                best_scores[gene_name] = 0
            if gene_name not in best_hit_files:
                best_hit_files[gene_name] = ''

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

                        # Update the best hit if the score is higher
                        if float(t_identity) > float(best_scores[gene_name]):
                            best_hits[gene_name] = reference
                            best_scores[gene_name] = float(t_identity)
                            best_hit_files[gene_name] = report_file
                    except ValueError:
                        pass

        # Create a dictionary of allele name: allele db path
        for _, allele in best_hits.items():
            for allele_db_file in sequence_dict['allele_output_paths']:
                # Ensure that the allele exists
                if not allele:
                    continue
                # Check if the allele is in the allele database
                if allele.replace('.', '_') in allele_db_file:
                    allele_db[allele] = allele_db_file

        # Update the sample dictionary
        sequence_dict['best_hits'] = best_hits
        sequence_dict['best_scores'] = best_scores
        sequence_dict['best_hit_files'] = best_hit_files
        sequence_dict['allele_db'] = allele_db

        for allele, identity in sequence_dict['best_scores'].items():
            print('\t', allele, identity, sequence_dict['best_hits'][allele])

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

        # Initialise a dictionary to store the BWA output
        sequence_dict['bwa_output'] = {}

        # Initialise a dictionary to store the allele database
        sequence_dict['reference_allele_file'] = {}

        for allele, allele_db in sequence_dict['allele_db'].items():

            # Set the name of the .fsa output file from KMA to be used as the
            # reference alleles
            reference_allele = allele_db + '.fsa'

            # Extract the parent directory of the reference allele
            reference_allele_dir = os.path.dirname(allele_db)

            # Index the reference allele if it was not already done
            if not os.path.exists(reference_allele + '.bwt'):
                _index_reference_allele(reference_allele=reference_allele)

            # Set the name of the output file
            output_bam = os.path.join(
                reference_allele_dir,
                f'{allele}.bam'
            )

            # Update the sample dictionaries
            sequence_dict['bwa_output'][allele] = output_bam
            sequence_dict['reference_allele_file'][allele] = reference_allele

            # Construct BWA command
            bwa_cmd = [
                'bwa', 'mem', '-t', str(threads), reference_allele,
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
                        bwa_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                    )

                    # Execute Samtools view
                    samtools_view_process = subprocess.Popen(
                        samtools_view_cmd, stdin=bwa_process.stdout,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE
                    )
                    bwa_process.stdout.close()  # Allow BWA to exit

                    # Execute Samtools sort
                    samtools_sort_process = subprocess.Popen(
                        samtools_sort_cmd, stdin=samtools_view_process.stdout,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE
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
                        "Error: One of the required executables not found: %s",
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


def _consensus_sequence(
    *,  # Enforce keyword arguments,
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]],
    min_coverage: int = 1,
    min_quality: int = 20
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Iterate through the BWA output files to extract the consensus sequences

    Args:
        sample_dict: Dictionary containing the sample information
        min_coverage: Minimum coverage required to call a base
        min_quality: Minimum base quality required to consider a base

    Returns:
        dict: The updated sample dictionary
    """
    for _, sequence_dict in sample_dict.items():
        print(_)
        # Initialize a dictionary to store the consensus sequences
        sequence_dict['consensus_sequences'] = {}

        for allele, bam_file in sequence_dict['bwa_output'].items():
            print('\t', allele)
            # Extract the consensus sequence
            consensus_sequence = _extract_consensus_sequence(
                bam_file=bam_file,
                min_coverage=min_coverage,
                min_quality=min_quality
            )
            # print('\t', allele, '\n\t', consensus_sequence)

            # Update the sample dictionary
            sequence_dict['consensus_sequences'][allele] = consensus_sequence

    return sample_dict


def reverse_complement(
    *,  # Enforce keyword arguments
    seq: str
):
    """
    Returns the reverse complement of a DNA sequence.

    Args:
        seq: A string representing the DNA sequence.

    Returns:
        A string representing the reverse complement of the input sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))


def create_pseudo_alignment(soft_clips_by_position):
    """
    Creates a pseudo-alignment of soft-clipped sequences, grouped by reference
    position.

    Args:
        soft_clips_by_position: A dictionary where keys are reference
        positions and values are lists of soft clip data dictionaries.

    Returns:
        A tuple containing:
            - A string representing the column-aligned sequences.
            - A list of the padded sequences.
    """

    pseudo_alignments = {}

    for ref_pos, soft_clip_data_list in soft_clips_by_position.items():
        # Determine the minimum relative alignment start
        min_start = min(data['relative_alignment_start'] for data in soft_clip_data_list)

        # Determine the maximum length of the aligned sequences
        max_length = 0
        for data in soft_clip_data_list:
            length = len(data['sequence']) - data['relative_alignment_start']
            if length > max_length:
                max_length = length

        # Create padded sequences
        padded_sequences = []
        for data in soft_clip_data_list:
            sequence = data['sequence']
            start = data['relative_alignment_start']
            padding_left = '-' * abs(min_start + start)
            padding_right = '-' * (max_length - len(sequence) + start)
            padded_sequence = padding_left + sequence + padding_right
            padded_sequences.append(padded_sequence)

        # Create column-aligned string
        column_alignment = ""
        for seq in padded_sequences:
            column_alignment += seq + "\n"  # Add the whole sequence and a newline

        # Store the pseudo-alignment
        pseudo_alignments[ref_pos] = (column_alignment, padded_sequences)

    return pseudo_alignments

def _extract_consensus_sequence(
    *,  # Enforce keyword arguments
    bam_file: str,
    min_coverage: int = 1,
    min_quality: int = 20,
    debug_allele: str = "Stx2a_1_X07865"  # Added debug_allele parameter
) -> str:
    """
    Extracts the consensus sequence from a sorted BAM file, including
    internal soft clips, but excluding soft clips that extend beyond the
    reference sequence boundaries.

    Args:
        bam_file: Path to the sorted BAM file.
        min_coverage: Minimum coverage required to call a base.
        min_quality: Minimum base quality for a base to be considered.
        debug_allele: The allele to print debug information for.

    Returns:
        The consensus sequence as a string.
    """

    # Open the BAM file
    sam_file = pysam.AlignmentFile(bam_file, "rb")

    # Get the reference sequence length. Assuming only one reference
    # sequence in the BAM file
    ref_len = sam_file.lengths[0]
    ref_name = sam_file.references[0]

    # Initialize a dictionary to store base counts at each position
    base_counts = defaultdict(
        lambda: {
            'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0
        }
    )

    # Initialize a dictionary to store soft-clipped sequences, grouped by reference position
    soft_clips_by_position = defaultdict(list)

    # Iterate through the aligned reads
    for read in sam_file.fetch(ref_name, 0, ref_len):

        # Get the aligned positions
        aligned_positions = read.get_reference_positions(full_length=True)

        # Get the query sequence
        query_sequence = read.query_sequence

        # Get the CIGAR string
        cigar = read.cigartuples

        # Iterate through the CIGAR operations
        query_pos = 0
        ref_pos = read.reference_start  # Start position on the reference
        ref_end = read.reference_end

        for operation, length in cigar:
            if operation == pysam.CMATCH or operation == pysam.CEQUAL or operation == pysam.CDIFF:
                # Match or mismatch
                for i in range(length):
                    if ref_pos < ref_len:  # Only consider positions within the reference
                        base = query_sequence[query_pos]
                        base_counts[ref_pos][base] += 1
                    query_pos += 1
                    ref_pos += 1
            elif operation == pysam.CINS:
                # Insertion (soft clip within the alignment)
                insertion_seq = query_sequence[query_pos:query_pos + length]
                # Consider the insertion as contributing to the *next* reference position
                if ref_pos < ref_len:
                    for base in insertion_seq:
                        base_counts[ref_pos]['I:' + base] = base_counts[ref_pos].get('I:' + base, 0) + 1
                query_pos += length
            elif operation == pysam.CDEL:
                # Deletion
                ref_pos += length
            elif operation == pysam.CSOFT_CLIP:
                # Soft clip - ignore if at the start or end of the alignment
                five_prime_clip = (ref_pos <= 10)  # Within 10 bases of the start
                three_prime_clip = (ref_pos + length >= ref_len - 10)  # Within 10 bases of the end
                cigar_left_clip = (query_pos == 0)
                cigar_right_clip = (query_pos + length == len(query_sequence))

                # Calculate the soft clip position
                soft_clip_pos = ref_pos

                if debug_allele in ref_name and not five_prime_clip and not three_prime_clip:
                    print(f"Debugging read: {read.query_name}")
                    print(f"  Cigar string: {read.cigartuples}")
                    print(f"  Reference start: {read.reference_start}")
                    print(f"  Reference end: {read.reference_end}")
                    print(f"  Query sequence: {read.query_sequence}")
                    print(f"  Soft clip at ref_pos {soft_clip_pos}, length {length}")
                    print(f"  Ref pos: {ref_pos}, ref length: {ref_len}")
                    print(f"    Is left clip: {five_prime_clip}, Is right clip: {three_prime_clip}")
                    print(f"    Is cigar left clip: {cigar_left_clip}, is cigar right clip: {cigar_right_clip}")
                    clip_seq = query_sequence[query_pos:query_pos + length]
                    # Adjust clip_seq for reverse reads
                    if read.is_reverse:
                        clip_seq = reverse_complement(clip_seq)
                    # Calculate relative alignment start
                    if cigar_left_clip:
                        relative_alignment_start = -length
                    else:
                        relative_alignment_start = 0
                    print(f"    Clip seq : {clip_seq}")
                    print(f"    Orientation: {'+' if not read.is_reverse else '-'}")
                    print(f"    Relative alignment start: {relative_alignment_start}")

                # Extract internal soft clips
                if not five_prime_clip and not three_prime_clip:
                    clip_seq = query_sequence[query_pos:query_pos + length]
                    # Adjust clip_seq for reverse reads
                    if read.is_reverse:
                        clip_seq = reverse_complement(clip_seq)

                    # Calculate relative alignment start
                    if cigar_left_clip:
                        relative_alignment_start = -length
                    else:
                        relative_alignment_start = 0

                    soft_clip_data = {
                        "sequence": clip_seq,
                        "ref_pos": ref_pos,
                        "query_pos": query_pos,
                        "length": length,
                        "five_prime_clip": five_prime_clip,
                        "three_prime_clip": three_prime_clip,
                        "cigar_left_clip": cigar_left_clip,
                        "cigar_right_clip": cigar_right_clip,
                        "orientation": "+" if not read.is_reverse else "-",  # Add orientation
                        "read_name": read.query_name,
                        "relative_alignment_start": relative_alignment_start
                    }
                    soft_clips_by_position[ref_pos].append(soft_clip_data)
                    # Again, contribute to the *next* reference position
                    if ref_pos < ref_len:
                        for base in clip_seq:
                            base_counts[ref_pos]['S:' + base] = base_counts[ref_pos].get('S:' + base, 0) + 1
                query_pos += length
            elif operation == pysam.CREF_SKIP:
                # Reference skip
                ref_pos += length
            elif operation == pysam.CHARD_CLIP or operation == pysam.CPAD:
                # Hard clip or padding - ignore
                pass
            else:
                raise ValueError(f"Unexpected CIGAR operation: {operation}")

    # Determine the consensus base at each position
    consensus_sequence = ""
    for pos in range(ref_len):
        # Find the base with the highest count
        max_base = 'N'
        max_count = 0
        for base, count in base_counts[pos].items():
            if count > max_count:
                max_base = base
                max_count = count
        consensus_sequence += max_base

    # Close the BAM file
    sam_file.close()

    # Create pseudo-alignments
    pseudo_alignments = create_pseudo_alignment(soft_clips_by_position)
    if debug_allele in ref_name:
        for pos, (column_alignment, padded_sequences) in pseudo_alignments.items():
            if pos == 906:
                print(f"Pseudo-alignments at position {pos}:")
                print(column_alignment)
                # print("Padded Sequences:")
                # for seq in padded_sequences:
                    # print(seq)

    return consensus_sequence


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


def cli():
    """
    Collect arguments from the command line
    """
    parser = ArgumentParser(
        description='Run KMA with the combined allele database, find the '
        'best hits, extract FASTQ reads, and run BWA to align the reads to '
        'the reference allele to confirm the outputs, and also return the '
        'full sequence of any alleles with insertions'
    )
    parser.add_argument(
        '-v', '--version', action='version',
        version=f'%(prog)s commit {__version__}'
    )
    parser.add_argument(
        '-s', '--sequence_path',
        required=True,
        help='Path of the folder containing sequencing reads'
    )
    parser.add_argument(
        '-r', '--reference_path',
        required=True,
        help='Name and path of the index KMA database to use'
    )
    parser.add_argument(
        '-t', '--threads',
        default=os.cpu_count(),
        help='Number of threads to use. Default is the number of cores in the '
        'system'
    )
    parser.add_argument(
        '-ID', '--identity',
        default=99,
        help='Minimum identity percentage for KMA hits. Default is 99%'
    )
    parser.add_argument(
        '--verbosity',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO',
        help='Set the logging verbosity level. Default is INFO'
    )

    # Parse the arguments into an object
    arguments = parser.parse_args()

    # Set up the logging
    setup_logging(
        verbosity=arguments.verbosity
    )

    return arguments


def tilde_expand(*, path: str) -> str:
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


if __name__ == '__main__':
    args = cli()

    # Expand the tilde in the paths
    args.sequence_path = tilde_expand(path=args.sequence_path)
    args.database_path = tilde_expand(path=args.reference_path)

    main(
        sequence_path=args.sequence_path,
        database_path=args.database_path,
        threads=args.threads,
        identity=args.identity
    )
