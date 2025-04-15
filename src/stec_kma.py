#!/usr/bin/env python3

"""
This is the main module of the STEC_KMA package. It runs KMA with the combined
allele database, finds the best hits, extracts FASTQ reads, and runs BWA to
align the reads to the reference allele to confirm the outputs, and also return
any alleles with insertions
"""

# Standard imports
from argparse import ArgumentParser
from glob import glob
import logging
import os
import time

# Local imports
from src.methods import (
    _bait_reads_bbmap,
    _calculate_stx_profile,
    _extract_allele_sequence,
    _extract_consensus_insertions_and_sequence,
    _extract_kma_mapped_reads,
    _extract_reads_by_name,
    _find_samples,
    _index_baited_sequences,
    _map_allele_sequences_kma,
    _map_reads_bwa,
    _parse_allele_reports,
    _reverse_bait_targets_bbmap,
    _setup_logging,
    _write_novel_allele_sequences,
    _write_read_names_file,
    _write_report,
)
from src.version import __version__

__author__ = 'adamkoziol'


def main(
    *,  # Enforce keyword arguments
    sequence_path: str,
    database_path: str,
    report_path: str,
    min_coverage: float,
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
        report_path: Path of the folder to write the reports
        min_coverage: Minimum fraction of reads required to call a consensus
        insertion.
        threads: Number of threads to use
        identity: Minimum identity percentage for KMA hits
    """
    # Determine the starting time (to be used when calculating the elapsed
    # time)
    start_time = time.time()

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

    logging.info('Indexing the baited sequences with KMA')
    sample_dict = _index_baited_sequences(
        sample_dict=sample_dict
    )

    logging.info('Mapping reads to the reference alleles with KMA')
    sample_dict = _map_allele_sequences_kma(
        identity=50,
        sample_dict=sample_dict,
        threads=threads
    )
    sample_dict = _parse_allele_reports(sample_dict=sample_dict)

    logging.info('Extracting the allele sequences from the baited targets')
    sample_dict = _extract_allele_sequence(sample_dict=sample_dict)

    logging.info('Extracting the mapped reads from the baited FASTQ files ')
    sample_dict = _extract_kma_mapped_reads(sample_dict=sample_dict)
    sample_dict = _write_read_names_file(sample_dict=sample_dict)
    sample_dict = _extract_reads_by_name(sample_dict=sample_dict)

    logging.info('Mapping reads to the reference alleles with BWA')
    sample_dict = _map_reads_bwa(
        sample_dict=sample_dict,
        threads=threads
    )

    logging.info(
        'Extracting consensus insertions and sequence from the BWA output'
    )
    sample_dict = _extract_consensus_insertions_and_sequence(
        min_coverage=min_coverage,
        sample_dict=sample_dict
        )

    logging.info('Writing novel allele sequences to FASTA files')
    _write_novel_allele_sequences(
        sample_dict=sample_dict,
        report_path=report_path
    )

    logging.info('Calculating stx profiles')
    sample_dict = _calculate_stx_profile(
        sample_dict=sample_dict
    )

    logging.info('Writing reports')
    _write_report(
        report_path=report_path,
        sample_dict=sample_dict
    )

    # Calculate and log the elapsed time
    elapsed_time = time.time() - start_time
    logging.info(
        'STEC calculation completed in %.2f seconds (%.2f minutes)',
        elapsed_time, elapsed_time / 60
    )


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
        '-d', '--database_path',
        required=True,
        help='Name and path of the index KMA database to use. This should be '
        'the path to the folder as well as the index name. For example, '
        '/path/to/database/allele_db'
    )
    parser.add_argument(
        '-r', '--report_path',
        help='Path to folder into which the report file is to be written. '
        'Default is the same folder as the input sequence files'
    )
    parser.add_argument(
        '-t', '--threads',
        default=os.cpu_count(),
        help='Number of threads to use. Default is the number of cores in the '
        'system'
    )
    parser.add_argument(
        '-ID', '--identity',
        default=90,
        help='Minimum identity percentage for KMA hits. Default is 90%%'
    )
    parser.add_argument(
        '-c', '--min_coverage',
        default=0.7,
        help='Minimum fraction of reads required to call a consensus '
        'insertion. For example, 0.7 means the insertion must be present in '
        'at least 70%%  of the reads covering that position. '
        'Default is 0.7'
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
    _setup_logging(
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
    args.database_path = tilde_expand(path=args.database_path)

    # Set the report path to the sequence path if not provided
    if args.report_path is None:
        args.report_path = os.path.join(
            args.sequence_path,
            'reports'
        )
    else:
        args.report_path = tilde_expand(path=args.report_path)

    # Create the report path if it does not exist
    os.makedirs(args.report_path, exist_ok=True)

    main(
        sequence_path=args.sequence_path,
        database_path=args.database_path,
        report_path=args.report_path,
        min_coverage=args.min_coverage,
        threads=args.threads,
        identity=args.identity
    )
