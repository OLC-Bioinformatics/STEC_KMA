#!/usr/bin/env python

"""
Integration test for the STEC_KMA pipeline using randomreads.sh from BBMap.
This test generates synthetic FASTQ files from actual FASTA genome files,
runs the STEC_KMA pipeline, and validates the outputs against expected results.
"""

# Standard imports
import csv
import os
import re
from pathlib import Path
import shutil
from subprocess import run

# Local imports
from src.methods import _tilde_expand
from src.stec_kma import main


def generate_fastq_with_bbmap(
    *,  # Enforce keyword arguments
    fasta_file: str,
    output_dir: str,
    sample_name: str,
    sample_number: int = 0,
    coverage: int = 20,
    read_length: int = 150
) -> tuple[str, str]:
    """
    Generate synthetic FASTQ files from a FASTA file using randomreads.sh from
    BBMap.

    Args:
        fasta_file: Path to the input FASTA file.
        output_dir: Directory where the output FASTQ files will be saved.
        sample_name: Name of the sample (used in output filenames).
        sample_number: Sample number (default: 0).
        coverage: Desired coverage for the generated reads (default: 20).
        read_length: Length of the generated reads (default: 150).

    Returns:
        tuple[str, str]: Paths to the generated R1 and R2 FASTQ files.
    """
    # Set the name and path of the output FASTQ files
    r1_file = os.path.join(
        output_dir,
        f"{sample_name}_S{sample_number}_L001_R1_001.fastq.gz"
    )
    r2_file = os.path.join(
        output_dir,
        f"{sample_name}_S{sample_number}_L001_R2_001.fastq.gz"
    )

    # Create the system command to run randomreads.sh
    call = [
        "randomreads.sh",
        f"ref={fasta_file}",
        f"out1={r1_file}",
        f"out2={r2_file}",
        f"length={read_length}",
        f"coverage={coverage}",
        "paired=t",
        "gzip=t",
        "q=30",  # Good quality reads
        "adderrors=t",  # Add realistic errors
        "maxsnps=0",  # No SNPs
        "maxinss=0",  # No insertions
        "maxdels=0",  # No deletions
        "maxsubs=0",  # No substitutions
    ]

    # Run randomreads.sh to generate paired-end FASTQ files
    run(call, check=True)

    return r1_file, r2_file


def get_expected_results():
    """
    Creates a dictionary of expected results for each test sample.
    """
    expected = {
        "2000-SEQ-0000": {
            "profile": "ND",
            "alleles": []
        },
        "2000-SEQ-0001": {
            "profile": "Stx1c;Stx2a;Stx2b",
            "alleles": [
                {"name": "Stx1c_12_Z36901", "identity": 100.0, "notes": ""},
                {"name": "Stx2b_116_AF043627", "identity": 100.0, "notes": ""},
                {"name": "Stx2a_1_AB030484", "identity": 99.6, "notes": ""}
            ]
        },
        "2000-SEQ-0002": {
            "profile": "Stx1a;Stx2a",
            "alleles": [
                {
                    "name": "Stx2a_1_X07865",
                    "identity": 99.92,
                    "notes": "Insertion:Position(901):MedianLength(134)"},
                {"name": "Stx1a_9_AB035142", "identity": 100.0, "notes": ""}
            ]
        },
        "2000-SEQ-0003": {
            "profile": "Stx2a",
            "alleles": [
                {"name": "Stx2a_1_AB030484", "identity": 100.0, "notes": ""}
            ]
        }
    }
    return expected


def validate_results(
    *,  # Enforce keyword arguments
    report_file: str,
    sample_name: str,
    expected_results: dict
) -> bool:
    """
    Validates STEC_KMA results against expected results using standard library.

    Args:
        report_file: Path to the generated report file
        sample_name: Name of the sample being tested
        expected_results: Dictionary with expected results

    Returns:
        bool: True if results match expectations, False otherwise
    """
    if not os.path.exists(report_file):
        print(f"ERROR: Report file {report_file} does not exist")
        return False

    # Read the TSV report file and filter for the current sample
    sample_results = []

    try:
        with open(report_file, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row['Sample'] == sample_name:
                    sample_results.append(row)
    except FileNotFoundError as exc:
        print(
            f"ERROR: Failed to locate report file: {exc}"
        )
        return False

    # Get expected results for this sample
    expected = expected_results.get(sample_name, {})

    # Check for ND (Not Detected) case
    if expected.get("profile") == "ND":
        if len(sample_results) == 1 and sample_results[0]['Allele'] == "ND":
            print(
                f"✓ {sample_name}: Correctly identified as negative (ND)")
            return True
        print(
            f"✗ {sample_name}: Expected ND but got "
            f"{len(sample_results)} results"
        )
        return False

    # For positive samples, check alleles
    if len(sample_results) != len(expected.get("alleles", [])):
        print(
            f"✗ {sample_name}: Expected {len(expected.get('alleles', []))} "
            f"alleles but got {len(sample_results)}"
        )
        return False

    # Check each allele
    all_match = True
    for exp_allele in expected.get("alleles", []):
        # Find matching allele in results
        matched = False
        for result in sample_results:
            if result['Allele'] == exp_allele['name']:
                # Check identity (with some tolerance)
                result_identity = float(result['PercentIdentity'])
                if abs(result_identity - exp_allele['identity']) > 0.5:
                    print(
                        f"✗ {sample_name}: Allele {exp_allele['name']} "
                        f"identity mismatch. Expected "
                        f"{exp_allele['identity']}, got {result_identity}"
                    )
                    all_match = False

                # Check notes if applicable
                notes = result.get('Notes', '')
                if (
                    exp_allele['notes'] and
                    exp_allele['notes'] not in str(notes)
                ):
                    print(
                        f"✗ {sample_name}: Notes mismatch for "
                        f"{exp_allele['name']}. Expected "
                        f"'{exp_allele['notes']}', got '{notes}'"
                    )
                    all_match = False

                matched = True
                break

        if not matched:
            print(
                f"✗ {sample_name}: Expected allele {exp_allele['name']} "
                f"not found in results"
            )
            all_match = False

    return all_match


def test_stec_kma_pipeline_with_real_genomes():
    """
    Integration test for the STEC_KMA pipeline using real test genomes.
    """
    # Set up local temporary directories
    temp_dir = os.path.join(Path.cwd(), "tests", "temp")
    os.makedirs(temp_dir, exist_ok=True)
    sequence_dir = os.path.join(temp_dir, "sequences")
    database_dir = os.path.join(temp_dir, "database")
    report_dir = os.path.join(temp_dir, "reports")
    os.makedirs(sequence_dir, exist_ok=True)
    os.makedirs(database_dir, exist_ok=True)
    os.makedirs(report_dir, exist_ok=True)

    # Path to test files directory
    test_files_dir = os.path.join(Path.cwd(), "tests", "files")
    database_path = os.path.join(
        test_files_dir, 'ref', 'StxDB'
    )

    # Get expected results
    expected_results = get_expected_results()

    all_tests_passed = True

    # Process each test genome
    for i, fasta_file in enumerate(
        sorted(Path(test_files_dir).glob("*.fasta"))
    ):
        # Extract sample name (without the part after underscore)
        sample_name = re.match(r"(\d+-SEQ-\d+)_", fasta_file.name).group(1)

        # Generate synthetic FASTQ files
        generate_fastq_with_bbmap(
            fasta_file=str(fasta_file),
            output_dir=sequence_dir,
            sample_name=sample_name,
            sample_number=i + 1
        )

    # Run the pipeline
    main(
        sequence_path=_tilde_expand(path=sequence_dir),
        database_path=_tilde_expand(path=database_path),
        report_path=_tilde_expand(path=report_dir),
        min_coverage=0.7,
        threads=4,
        identity=90
    )

    # Validate the results
    report_file = os.path.join(report_dir, "stec_kma_report.tsv")
    result = validate_results(
        report_file=report_file,
        sample_name=sample_name,
        expected_results=expected_results
    )
    all_tests_passed = all_tests_passed and result

    assert all_tests_passed, \
        "One or more tests failed. See the output above for details."

    # Clean up temporary directories
    shutil.rmtree(temp_dir, ignore_errors=True)
