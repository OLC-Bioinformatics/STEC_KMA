#!/usr/bin/env python3

"""
This is the main module of the STEC_KMA package. It runs KMA with the combined
allele database, finds the best hits, extracts FASTQ reads, and runs BWA to
align the reads to the reference allele to confirm the outputs, and also return
any alleles with insertions
"""

# Standard imports
from argparse import ArgumentParser
from concurrent.futures import as_completed, ThreadPoolExecutor
from glob import glob
import gzip
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

    logging.info('Running KMA with the combined allele database')
    sample_dict = _run_combined_db_kma(
        database_path=database_path,
        sample_dict=sample_dict,
        threads=threads
    )

    logging.info('Extracting mapped reads')
    sample_dict = _read_frag_files(sample_dict=sample_dict)
    sample_dict = _ensure_forward_reverse(sample_dict=sample_dict)
    sample_dict = extract_reads_bbmap(sample_dict=sample_dict)

    logging.info(
        'Splitting reference database into individual alleles for KMA analyses'
    )
    sample_dict = _extract_sequences(
        db_fasta=db_fasta,
        sample_dict=sample_dict
        )
    sample_dict = _write_allele_sequences(sample_dict=sample_dict)
    sample_dict = _index_allele_sequences(sample_dict=sample_dict)

    logging.info('Mapping reads to the reference alleles with KMA')    
    sample_dict = _map_allele_sequences_kma(
        identity=identity,
        sample_dict=sample_dict,
        threads=threads
    )
    sample_dict = _parse_allele_reports(sample_dict=sample_dict)
    quit()
    logging.info('Mapping reads to the reference alleles with BWA')
    # sample_dict = map_reads_bwa(
    #     sample_dict=sample_dict,
    #     threads=threads
    # )


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


def _run_combined_db_kma(
    database_path: str,
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]],
    threads: int,
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Run KMA with the combined allele database

    Args:
        database_path: Path of the folder containing the database
        sample_dict: Dictionary containing the sample information
        threads: Number of threads to use

    Returns:
        dict: The updated sample dictionary
    """
    for sample, sequence_dict in sample_dict.items():

        # Output path for the sample
        output_path = os.path.join(sequence_dict["output_path"], sample)
        sequence_dict['kma_outputs'] = output_path

        # Create the system call
        combined_kma_call = (
            f'kma -ipe {" ".join(sequence_dict["files"])} '
            f'-o {output_path} '
            f'-t_db {database_path} '
            f'-t {threads} '    # Threading
            f'-ID 50 '  # Identity threshold
            '-mem_mode '        # Memory mode
            '-a '               # Output all matches
        )

        logging.info('System call: %s', combined_kma_call)

        # Run the command if the outputs do not already exist
        if not os.path.exists(output_path + '.res'):
            output = run_command(command=combined_kma_call)
            logging.debug('KMA output: %s', output.stdout)

    return sample_dict


def _read_frag_files(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]]
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Extract the mapped read names from the KMA .frag.gz files
    """
    # Iterate through the samples
    for _, sequence_dict in sample_dict.items():
        # Initialize sets to store the read names and allele hits
        read_names = set()
        alleles = set()

        # Create the path to the .frag.gz file
        frag_file = sequence_dict['kma_outputs'] + '.frag.gz'

        # Read in the .frag.gz file
        with gzip.open(frag_file, 'rt') as frag_fh:
            for line in frag_fh:
                # Extract the read name and add it to the list
                read_name = line.strip().split('\t')[-1]
                read_names.add(read_name)

                # Extract the allele name and add it to the list
                allele = line.strip().split('\t')[-2]
                alleles.add(allele)

        # Update the sample dictionary
        sequence_dict['read_names'] = sorted(list(read_names))
        sequence_dict['alleles'] = sorted(list(alleles))

    return sample_dict


def _ensure_forward_reverse(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]]
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Ensure that both forward and reverse reads are present in the read_names
    list. If one direction is missing, add it to the list.

    Args:
        sample_dict: Dictionary containing the sample information

    Returns:
        dict: The updated sample dictionary
    """
    # Iterate through the samples
    for _, sequence_dict in sample_dict.items():

        # Use a set for efficient membership checking
        read_names = set(sequence_dict['read_names'])

        # Create a copy
        updated_read_names = set(read_names)

        # Iterate through a copy of the set
        for read_name in list(read_names):
            if ' 1' in read_name:
                reverse_read_name = read_name.replace(' 1', ' 2')
                updated_read_names.add(reverse_read_name)
            elif ' 2' in read_name:
                forward_read_name = read_name.replace(' 2', ' 1')
                updated_read_names.add(forward_read_name)

        # Convert back to a sorted list
        sequence_dict['read_names'] = sorted(list(updated_read_names))

    return sample_dict


def extract_reads_bbmap(
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]],
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Extracts reads from paired FASTQ files using BBMap's filterbyname.sh.

    Args:
        sample_dict: Dict with sample information.

    Returns:
        Dict with updated sample information.
    """
    for sample, sequence_dict in sample_dict.items():

        # Ensure there are two files (paired-end reads)
        if len(sequence_dict['files']) != 2:
            logging.error(
                "Error: Expected two FASTQ files for paired-end reads."
            )
            continue  # Skip to the next sample

        fastq_file1 = sequence_dict['files'][0]
        fastq_file2 = sequence_dict['files'][1]

        # Set the name of the output files
        output_fastq1 = os.path.join(
            sequence_dict['output_path'],
            f'{sample}_R1.fastq.gz'
        )
        output_fastq2 = os.path.join(
            sequence_dict['output_path'],
            f'{sample}_R2.fastq.gz'
        )

        # Add the output files to the list
        sequence_dict['mapped_reads'] = [output_fastq1, output_fastq2]

        # Only extract reads if the output files do not already exist
        if os.path.exists(output_fastq1) and os.path.exists(output_fastq2):
            if (
                os.path.getsize(output_fastq1) > 0
                and os.path.getsize(output_fastq2) > 0
            ):
                logging.info(
                    "Output files already exist for sample: %s", sample
                )
                continue

        # Extract the read names
        read_names = sequence_dict['read_names']

        # Create a file to store the read names
        names_file = os.path.join(
            sequence_dict['output_path'],
            f'{sample}_names.txt'
        )

        try:
            # Write the read names to file
            with open(names_file, 'w', encoding='utf-8') as reads_file:
                for read_name in read_names:
                    reads_file.write(read_name + '\n')

            # Construct the BBMap command
            bbmap_cmd = [
                'filterbyname.sh',
                f'in={fastq_file1}',
                f'in2={fastq_file2}',
                f'out={output_fastq1}',
                f'out2={output_fastq2}',
                f'names={names_file}',
                'overwrite=t'  # Add overwrite flag
            ]

            # Execute the BBMap command
            with subprocess.Popen(
                bbmap_cmd, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE
            ) as process:
                _, stderr = process.communicate()

                if stderr:
                    logging.error(
                        "BBMap error: %s", stderr.decode()
                    )

        except FileNotFoundError:
            logging.error(
                "Error: Input FASTQ not found."
            )
        except Exception as exc:
            logging.error("An error occurred: %s", exc)
            raise SystemExit from exc

        # Sort the list of mapped reads
        sequence_dict['mapped_reads'] = sorted(
            sequence_dict['mapped_reads']
        )

    return sample_dict


def _extract_sequences(
    *,  # Enforce keyword arguments
    db_fasta: str,
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]]
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Extract the sample-specific matching alleles from the database file

    Args:
        db_fasta (str): Path to the gene-specific FASTA file
        sample_dict (dict): Dictionary containing the sample information

    Returns:
        dict: The updated sample dictionary
    """

    for _, sequence_dict in sample_dict.items():

        with open(db_fasta, 'r', encoding='utf-8') as fasta_file:
            # Create a list to store the sequences
            sequences = []

            # Extract the allele names
            alleles = set(sequence_dict['alleles'])

            # Iterate through the FASTA records
            records = SeqIO.parse(fasta_file, 'fasta')
            for record in records:
                if record.id in alleles:
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

            # Write the sequence to the output file
            with open(output_fasta, 'w', encoding='utf-8') as fasta_file:
                # Write a list containing the single record
                SeqIO.write([record], fasta_file, 'fasta')

    return sample_dict


def _index_allele_sequences(
    *,  # Enforce keyword arguments
    sample_dict: Dict[str, Dict[str, Union[List[str], str]]]
) -> Dict[str, Dict[str, Union[List[str], str]]]:
    """
    Use KMA to index the allele sequences

    Args:
        sample_dict: Dictionary containing the sample information

    Returns:
        dict: The updated sample dictionary
    """
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

            logging.info('KMA index command: %s', kma_index_cmd)

            # Run the command if the index does not already exist
            if not os.path.exists(allele_db + '.name'):
                try:
                    # Execute the KMA index command
                    output = run_command(command=kma_index_cmd)
                    logging.debug('KMA index output: %s', output.stdout)
                except subprocess.CalledProcessError as exc:
                    logging.error("Command failed: %s", exc)
                    raise

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
        # Create a list to store the KMA commands
        kma_commands = []

        # Create a list to store the report files
        report_files = []
        sequence_dict['kma_report_files'] = []

        # Iterate through the allele databases
        for iterator, allele_db in enumerate(sequence_dict['kma_sample_db']):

            # Extract the output path
            output_path = sequence_dict['allele_output_paths'][iterator]

            # Construct the KMA command
            kma_cmd = (
                f'kma -ipe {" ".join(sequence_dict["files"])} '
                f'-o {output_path} '
                f'-t_db {allele_db} '
                f'-ID {identity} '  # Identity threshold
            )
            kma_commands.append(kma_cmd)

            # Construct the report file name
            report_file = output_path + '.res'
            report_files.append(report_file)

            # Update the sample dictionary
            sequence_dict['kma_report_files'] = report_files

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
        # Create a dictionaries to store the best hits, best scores, and
        # best hit files
        best_hits = {}
        best_scores = {}
        best_hit_files = {}
        
        # Iterate through the report files
        for report_file in sequence_dict['kma_report_files']:
            
            # Extract the gene name from the report file
            gene_name = _find_gene_name(allele=report_file)

            # Initialise the dictionaries for the gene
            if gene_name not in best_hits:
                best_hits[gene_name] = {}
                best_scores[gene_name] = ''
                best_hit_files[gene_name] = 0
            
            # Read the report file
            with open(report_file, 'r', encoding='utf-8') as report_fh:
                for line in report_fh:
                    print('line:', line)
                    # Skip the header
                    if line.startswith('#'):
                        continue
                    
                    try:
                        # Extract the reference allele and percentage identity
                        # Template	Score	Expected	Template_length
                        # Template_Identity	Template_Coverage	Query_Identity
                        # Query_Coverage	Depth	q_value	p_value
                        print(line.rstrip())
                        quit()
                        reference, _, _, _, t_identity, \
                            _, _, _, _, _, = line.strip().split('\t')
                        print(gene_name, reference, t_identity)
                        quit()
                        # Update the best hit if the score is higher
                        if int(t_identity) > int(best_scores[gene_name]):
                            best_hits[gene_name] = reference
                            best_scores[gene_name] = t_identity
                            best_hit_files[gene_name] = report_file
                    except ValueError:
                        pass
            
        # Update the sample dictionary
        sequence_dict['best_hits'] = best_hits
        sequence_dict['best_scores'] = best_scores
        sequence_dict['best_hit_files'] = best_hit_files
        
        print(_, 'best_hits:', sequence_dict['best_hits'])
    quit()
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


def map_reads_bwa(
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
    for sample, sequence_dict in sample_dict.items():

        # Set the name of the .fsa output file from KMA to be used as the
        # reference alleles
        reference_allele = sequence_dict['kma_sample_db'] + '.fsa'

        # Index the reference allele if it was not already done
        if not os.path.exists(reference_allele + '.bwt'):
            _index_reference_allele(reference_allele=reference_allele)

        # Set the name of the output file
        output_bam = os.path.join(
            sequence_dict['output_path'],
            f'{sample}.bam'
        )

        # Update the sample dictionary
        sequence_dict['bwa_output'] = output_bam

        # Construct BWA command
        bwa_cmd = [
            'bwa', 'mem', '-t', str(threads), reference_allele,
            sequence_dict["mapped_reads"][0],
            sequence_dict["mapped_reads"][1]
        ]

        # Construct Samtools view command
        samtools_view_cmd = ['samtools', 'view', '-Sb', '-']

        # Construct Samtools sort command
        samtools_sort_cmd = ['samtools', 'sort', '-o', output_bam]

        logging.info('BWA command: %s', ' '.join(bwa_cmd))
        logging.info('Samtools view command: %s', ' '.join(samtools_view_cmd))
        logging.info('Samtools sort command: %s', ' '.join(samtools_sort_cmd))

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
                    "Error: One of the required executables not found: %s", exc
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

    logging.info('System call: %s', ' '.join(index_cmd))

    try:
        # Run the command
        output = run_command(
            command=index_cmd,
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

    try:
        # Run the samtools index command
        output = subprocess.run(
            index_cmd, capture_output=True, text=True, check=True
        )
        logging.debug('Samtools index output: %s', output.stdout)
    except subprocess.CalledProcessError as exc:
        logging.error("Samtools index command failed: %s", exc)
        raise


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
        default=95,
        help='Minimum identity percentage for KMA hits. Default is 95%'
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
