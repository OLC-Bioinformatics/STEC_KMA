#!/usr/bin/env python

"""
Test suite for the methods in the STEC-KMA package.
"""

# Standard imports
import os
import re
import subprocess
from unittest.mock import (
    call,
    mock_open,
    patch,
    MagicMock
)

# Third-party imports
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam
import pytest

# Local imports
from src.methods import (
    _bait_reads_bbmap,
    _calculate_depth_per_position,
    _calculate_stx_profile,
    _detect_truncated_sequences,
    _extract_allele_sequence,
    _extract_consensus_insertions_and_sequence,
    _extract_kma_mapped_reads,
    _extract_reads_by_name,
    _find_consensus_internal_insertions_and_sequence,
    _find_samples,
    _index_baited_sequences,
    _index_bam_file,
    _index_reference_allele,
    _map_allele_sequences_kma,
    _map_reads_bwa,
    _parse_allele_reports,
    _reverse_bait_targets_bbmap,
    _setup_logging,
    _tilde_expand,
    _translate_allele_sequence,
    _write_allele_sequences,
    _write_novel_allele_sequences,
    _write_read_names_file,
    _write_report,
    run_command,
)
from src.version import __version__


def test_setup_logging():
    """
    Test _setup_logging function to ensure logging is set up correctly.
    """
    with patch("logging.basicConfig") as mock_logging:
        _setup_logging(verbosity="DEBUG")
        mock_logging.assert_called_once_with(
            level=10,  # DEBUG level
            format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        )


def test_find_samples_no_fastq_files(setup_temp_dir):
    """
    Test _find_samples when no .fastq.gz files are present in the directory.
    """
    sample_dict = _find_samples(sequence_path=setup_temp_dir)
    assert sample_dict == {}, \
        "Expected an empty dictionary when no files exist."


def test_find_samples_misnamed_fastq_files(setup_temp_dir):
    """
    Test _find_samples when .fastq.gz files are misnamed.
    """
    # Cre7ate mock .fastq.gz files
    fastq_file_1 = os.path.join(setup_temp_dir, "sample1_R1.fastq.gz")
    fastq_file_2 = os.path.join(setup_temp_dir, "sample1_R2.fastq.gz")

    # Create empty files to simulate fastq files
    with (
        open(fastq_file_1, "w", encoding='utf-8'),
        open(fastq_file_2, "w", encoding='utf-8')
    ):
        pass

    sample_dict = _find_samples(sequence_path=setup_temp_dir)
    assert sample_dict == {}, \
        f"Warning: Could not extract sample name from: {fastq_file_1}"


def test_find_samples_with_fastq_files(setup_temp_dir):
    """
    Test _find_samples when .fastq.gz files are present in the directory.
    """
    # Create mock .fastq.gz files
    fastq_file_1 = os.path.join(setup_temp_dir, "2222-FAKE-0000_R1.fastq.gz")
    fastq_file_2 = os.path.join(setup_temp_dir, "2222-FAKE-0000_R2.fastq.gz")

    # Create empty files to simulate fastq files
    with (
        open(fastq_file_1, "w", encoding='utf-8'),
        open(fastq_file_2, "w", encoding='utf-8')
    ):
        pass

    sample_dict = _find_samples(sequence_path=setup_temp_dir)
    assert "2222-FAKE-0000" in sample_dict, \
        "Expected 2222-FAKE-0000 to be in the sample dictionary."
    assert len(sample_dict["2222-FAKE-0000"]["files"]) == 2, \
        "Expected two fastq files."


@patch("subprocess.run")
def test_run_command(mock_subprocess_run):
    """
    Test run_command to ensure shell commands are executed correctly.
    """
    mock_subprocess_run.return_value.returncode = 0
    result = run_command(command="echo 'Hello, World!'", split=False)
    assert result.returncode == 0, \
        "Expected the command to execute successfully."
    mock_subprocess_run.assert_called_once_with(
        "echo 'Hello, World!'",
        capture_output=True,
        text=True,
        check=True,
        shell=True
    )


@patch("subprocess.run")
def test_run_command_split_true(mock_subprocess_run):
    """
    Test run_command with split=True to ensure the command is split correctly.
    """
    # Mock subprocess.run to simulate a successful command execution
    mock_subprocess_run.return_value.returncode = 0
    mock_subprocess_run.return_value.stdout = "Command executed successfully"

    # Call the function with split=True
    result = run_command(command="echo Hello, World!", split=True)

    # Verify the command was split and executed
    mock_subprocess_run.assert_called_once_with(
        ["echo", "Hello,", "World!"],
        capture_output=True,
        text=True,
        check=True
    )
    assert result.returncode == 0, \
        "Expected the command to execute successfully."
    assert result.stdout == "Command executed successfully", \
        "Expected the correct stdout output."


@patch("subprocess.run")
def test_run_command_split_false(mock_subprocess_run):
    """
    Test run_command with split=False to ensure the command is executed as a
    single string.
    """
    # Mock subprocess.run to simulate a successful command execution
    mock_subprocess_run.return_value.returncode = 0
    mock_subprocess_run.return_value.stdout = "Command executed successfully"

    # Call the function with split=False
    result = run_command(command="echo 'Hello, World!'", split=False)

    # Verify the command was executed as a single string
    mock_subprocess_run.assert_called_once_with(
        "echo 'Hello, World!'",
        capture_output=True,
        text=True,
        check=True,
        shell=True
    )
    assert result.returncode == 0, \
        "Expected the command to execute successfully."
    assert result.stdout == "Command executed successfully", \
        "Expected the correct stdout output."


@patch("subprocess.run")
def test_run_command_calledprocesserror(mock_subprocess_run, caplog):
    """
    Test run_command handles subprocess.CalledProcessError correctly.
    """
    # Mock subprocess.run to raise a CalledProcessError
    mock_subprocess_run.side_effect = subprocess.CalledProcessError(
        returncode=1,
        cmd="failing_command",
        output="Error occurred",
        stderr="Command failed"
    )

    # Call the function and expect it to raise the exception
    with caplog.at_level("ERROR"):
        with pytest.raises(subprocess.CalledProcessError):
            run_command(command="failing_command", split=True)


@patch("src.methods.run_command")
def test_bait_reads_bbmap(mock_run_command, setup_temp_dir):
    """
    Test _bait_reads_bbmap to ensure BBMap is called correctly.
    """
    # Mock sample dictionary
    sample_dict = {
        "2222-FAKE-0000": {
            "files": [
                os.path.join(setup_temp_dir, "2222-FAKE-0000_R1.fastq.gz"),
                os.path.join(setup_temp_dir, "2222-FAKE-0000_R2.fastq.gz"),
            ],
            "output_path": setup_temp_dir,
        }
    }

    # Call the function
    updated_sample_dict = _bait_reads_bbmap(
        db_fasta="mock_db.fasta", sample_dict=sample_dict
    )

    # Assert that run_command was called
    mock_run_command.assert_called_once()
    assert updated_sample_dict == sample_dict, \
        "Expected the 2222-FAKE-0000 dictionary to remain unchanged."


# Add this test to the existing test_methods.py file

@patch("src.methods.run_command")
def test_reverse_bait_targets_bbmap(mock_run_command, setup_temp_dir):
    """
    Test _reverse_bait_targets_bbmap to ensure BBMap is called correctly.
    """
    # Mock sample dictionary
    sample_dict = {
        "2222-FAKE-0000": {
            "baited_reads": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited.fastq.gz"
            ),
            "output_path": setup_temp_dir,
        }
    }

    # Mock the database FASTA file
    db_fasta = os.path.join(setup_temp_dir, "mock_db.fasta")

    # Call the function
    updated_sample_dict = _reverse_bait_targets_bbmap(
        db_fasta=db_fasta,
        identity=90,
        sample_dict=sample_dict
    )

    # Assert that the baited targets file path is correctly updated
    baited_targets = os.path.join(
        setup_temp_dir,
        "2222-FAKE-0000_baited_targets.fasta"
    )
    assert updated_sample_dict["2222-FAKE-0000"]["baited_targets"] == \
        baited_targets, \
        "Expected baited_targets path to be updated in the sample dictionary."

    # Assert that run_command was called with the correct BBMap command
    expected_command = (
        f"bbduk.sh ref={sample_dict['2222-FAKE-0000']['baited_reads']} "
        f"outm={baited_targets} "
        f"in={db_fasta} "
        f"mincovfraction=0.9 "
        "--maskmiddle=f "
        "--fixjunk"
    )
    mock_run_command.assert_called_once_with(
        command=expected_command, split=False
    )


@patch("src.methods.run_command")
@patch("os.path.exists")
@patch("os.path.getsize")
def test_index_baited_sequences(
    mock_getsize,
    mock_exists,
    mock_run_command,
    setup_temp_dir
):
    """
    Test _index_baited_sequences to ensure KMA indexing is called correctly.
    """
    # Mock sample dictionary
    sample_dict = {
        "2222-FAKE-0000": {
            "baited_targets": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets.fasta"
            ),
            "output_path": setup_temp_dir,
        }
    }

    # Mock the baited targets file to exist and have a non-zero size
    def mock_exists_side_effect(path):
        if path == sample_dict["2222-FAKE-0000"]["baited_targets"]:
            return True  # Baited targets file exists
        elif path == os.path.join(
            setup_temp_dir,
            "2222-FAKE-0000_baited_targets_db.name"
        ):
            return False  # Report file does not exist
        return False

    mock_exists.side_effect = mock_exists_side_effect
    mock_getsize.return_value = 100  # Non-zero size for baited targets file

    # Call the function
    updated_sample_dict = _index_baited_sequences(sample_dict=sample_dict)

    # Assert that the KMA database path is correctly updated
    expected_db_path = os.path.join(
        setup_temp_dir,
        "2222-FAKE-0000_baited_targets_db"
    )
    assert updated_sample_dict["2222-FAKE-0000"]["kma_sample_db"] == \
        expected_db_path, \
        "Expected kma_sample_db path to be updated in the sample dictionary."

    # Assert that run_command was called with the correct KMA index command
    expected_command = (
        f"kma index -i {sample_dict['2222-FAKE-0000']['baited_targets']} "
        f"-o {expected_db_path}"
    )
    mock_run_command.assert_called_once_with(
        command=expected_command,
        split=False
    )

    # Assert that the report file is checked for existence
    expected_report_file = expected_db_path + ".name"
    mock_exists.assert_any_call(expected_report_file)


@patch("os.path.exists")
@patch("os.path.getsize")
def test_index_baited_sequences_missing_file(
    mock_getsize,
    mock_exists,
    setup_temp_dir
):
    """
    Test _index_baited_sequences when the baited targets file is missing
    or empty.
    """
    # Mock sample dictionary
    sample_dict = {
        "2222-FAKE-0000": {
            "baited_targets": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets.fasta"
            ),
            "output_path": setup_temp_dir,
        }
    }

    # Mock the baited targets file to not exist or have zero size
    mock_exists.return_value = False
    mock_getsize.return_value = 0  # Zero size

    # Call the function
    updated_sample_dict = _index_baited_sequences(sample_dict=sample_dict)

    # Assert that the KMA database path is set to None
    assert updated_sample_dict["2222-FAKE-0000"]["kma_sample_db"] is None, \
        "Expected kma_sample_db to be None for missing or empty baited " \
        "targets file."


@patch("src.methods.run_command")
@patch("os.path.exists")
def test_map_allele_sequences_kma(
    mock_exists,
    mock_run_command,
    setup_temp_dir
):
    """
    Test _map_allele_sequences_kma to ensure KMA mapping is called correctly.
    """
    # Mock sample dictionary
    sample_dict = {
        "2222-FAKE-0000": {
            "kma_sample_db": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets_db"
            ),
            "baited_reads": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited.fastq.gz"
            ),
            "baited_targets": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets.fasta"
            ),
            "output_path": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_output"
            ),
        }
    }

    # Mock the report file to not exist
    def mock_exists_side_effect(path):
        if path.endswith(".res"):
            return False  # Report file does not exist
        return True

    mock_exists.side_effect = mock_exists_side_effect

    # Call the function
    updated_sample_dict = _map_allele_sequences_kma(
        identity=90,
        sample_dict=sample_dict,
        threads=4
    )

    # Assert that the KMA report file path is correctly updated
    expected_report_file = os.path.join(
        setup_temp_dir,
        "2222-FAKE-0000_baited_targets.res"
    )
    assert updated_sample_dict["2222-FAKE-0000"]["kma_report_file"] == \
        expected_report_file, \
        "Expected kma_report_file path to be updated in the sample dictionary."

    # Assert that run_command was called with the correct KMA mapping command
    expected_command = (
        f"kma -int {sample_dict['2222-FAKE-0000']['baited_reads']} "
        f"-o {os.path.splitext(
            sample_dict['2222-FAKE-0000']['baited_targets']
        )[0]} "
        f"-t_db {sample_dict['2222-FAKE-0000']['kma_sample_db']} "
        f"-ID 90 "
        f"-t 4 "
        "-mem_mode "
        "-a "
    )
    mock_run_command.assert_called_once_with(
        command=expected_command,
        split=False
    )


@patch("src.methods.run_command")
@patch("os.path.exists")
def test_map_allele_sequences_kma_no_kma_db(
    mock_exists,
    mock_run_command,
    setup_temp_dir
):
    """
    Test _map_allele_sequences_kma when sequence_dict['kma_sample_db'] is None.
    """
    # Mock sample dictionary with kma_sample_db set to None
    sample_dict = {
        "2222-FAKE-0000": {
            "kma_sample_db": None,  # No KMA database
            "baited_reads": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited.fastq.gz"
            ),
            "baited_targets": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets.fasta"
            ),
            "output_path": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_output"
            ),
        }
    }

    # Mock the report file to not exist
    mock_exists.return_value = False

    # Call the function
    updated_sample_dict = _map_allele_sequences_kma(
        identity=90,
        sample_dict=sample_dict,
        threads=4
    )

    # Assert that the KMA report file path is not updated
    assert "kma_report_file" not in updated_sample_dict["2222-FAKE-0000"], \
        "Expected kma_report_file to not be set when kma_sample_db is None."

    # Assert that run_command was not called
    mock_run_command.assert_not_called()


@patch("builtins.open", create=True)
def test_parse_allele_reports(mock_open, setup_temp_dir):
    """
    Test _parse_allele_reports to ensure KMA reports are parsed correctly.
    """
    # Mock sample dictionary
    sample_dict = {
        "2222-FAKE-0000": {
            "kma_sample_db": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets_db"
            ),
            "kma_report_file": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets.res"
            ),
        }
    }

    # Mock the content of the KMA report file
    mock_report_content = (
        "# Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\t"
        "Template_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\t"
        "p_value\n"
        "gene1_geneA_1\t100\t0.01\t1000\t99.5\t95.0\t99.0\t90.0\t30\t0.001\t"
        "0.05\n"
        "gene2_geneB_2\t200\t0.02\t1500\t98.0\t90.0\t97.0\t85.0\t25\t0.002\t"
        "0.04\n"
    )
    mock_open.return_value.__enter__.return_value = \
        mock_report_content.splitlines()

    # Call the function
    updated_sample_dict = _parse_allele_reports(sample_dict=sample_dict)

    # Assert that the best hits, scores, and hit files are updated correctly
    assert "gene1" in updated_sample_dict["2222-FAKE-0000"]["best_hits"], \
        "Expected gene1 to be in the best hits."
    assert updated_sample_dict["2222-FAKE-0000"]["best_hits"]["gene1"] == \
        "gene1_geneA_1", \
        "Expected best hit for gene1 to be gene1_geneA_1."
    assert updated_sample_dict["2222-FAKE-0000"]["best_scores"]["gene1"] == \
        99.5, \
        "Expected best score for gene1 to be 99.5."
    assert updated_sample_dict[
        "2222-FAKE-0000"
        ]["best_hit_files"]["gene1"] == \
        sample_dict["2222-FAKE-0000"]["kma_report_file"], \
        "Expected best hit file for gene1 to match the report file."


def test_parse_allele_reports_no_kma_db(setup_temp_dir):
    """
    Test _parse_allele_reports when sequence_dict['kma_sample_db'] is None.
    """
    # Mock sample dictionary with kma_sample_db set to None
    sample_dict = {
        "2222-FAKE-0000": {
            "kma_sample_db": None,  # No KMA database
            "kma_report_file": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets.res"
            ),
        }
    }

    # Call the function
    updated_sample_dict = _parse_allele_reports(sample_dict=sample_dict)

    # Assert that no changes are made to the sample dictionary
    assert "best_hits" not in updated_sample_dict["2222-FAKE-0000"], \
        "Expected best_hits to not be set when kma_sample_db is None."
    assert "best_scores" not in updated_sample_dict["2222-FAKE-0000"], \
        "Expected best_scores to not be set when kma_sample_db is None."
    assert "best_hit_files" not in updated_sample_dict["2222-FAKE-0000"], \
        "Expected best_hit_files to not be set when kma_sample_db is None."


@patch("builtins.open", create=True)
def test_parse_allele_reports_value_error(mock_open, setup_temp_dir):
    """
    Test _parse_allele_reports to ensure it handles ValueError gracefully
    when the KMA report file contains malformed lines.
    """
    # Mock sample dictionary
    sample_dict = {
        "2222-FAKE-0000": {
            "kma_sample_db": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets_db"
            ),
            "kma_report_file": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets.res"
            ),
        }
    }

    # Mock the content of the KMA report file with a malformed line
    mock_report_content = (
        "# Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\t"
        "Template_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\t"
        "p_value\n"
        "gene1_geneA_1\t100\t0.01\t1000\t99.5\t95.0\t99.0\t90.0\t30\t0.001\t"
        "0.05\n"
        "malformed_line_without_enough_columns\n"  # Malformed line
    )
    mock_open.return_value.__enter__.return_value = \
        mock_report_content.splitlines()

    # Call the function
    updated_sample_dict = _parse_allele_reports(sample_dict=sample_dict)

    # Assert that the function handled the malformed line gracefully
    assert "gene1" in updated_sample_dict["2222-FAKE-0000"]["best_hits"], \
        "Expected gene1 to be in the best hits despite the malformed line."
    assert updated_sample_dict["2222-FAKE-0000"]["best_hits"]["gene1"] == \
        "gene1_geneA_1", \
        "Expected best hit for gene1 to be gene1_geneA_1."
    assert updated_sample_dict["2222-FAKE-0000"]["best_scores"]["gene1"] == \
        99.5, \
        "Expected best score for gene1 to be 99.5."
    assert updated_sample_dict[
        "2222-FAKE-0000"
        ]["best_hit_files"]["gene1"] == \
        sample_dict["2222-FAKE-0000"]["kma_report_file"], \
        "Expected best hit file for gene1 to match the report file."


@patch("src.methods.SeqIO.parse")
@patch("src.methods._write_allele_sequences")
def test_extract_allele_sequence(
    mock_write_allele_sequences,
    mock_seqio_parse,
    setup_temp_dir
):
    """
    Test _extract_allele_sequence to ensure allele sequences are extracted
    correctly.
    """
    # prepare sample_dict
    sample_dict = {
        "2222-FAKE-0000": {
            "kma_sample_db": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets_db"
            ),
            "baited_targets": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets.fasta"
            ),
            "best_hits": {"gene1": "geneA_1", "gene2": "geneB_2"},
            "output_path": setup_temp_dir,
        }
    }

    # simulate FASTA records
    mock_seqio_parse.return_value = [
        SeqRecord(Seq("ATGC"), id="geneA_1"),
        SeqRecord(Seq("CGTA"), id="geneB_2"),
        SeqRecord(Seq("GCTA"), id="allele3"),
    ]

    # side effect writes and returns a path
    def write_side_effect(record, sequence_dict):
        assert sequence_dict
        return os.path.join(setup_temp_dir, f"{record.id}.fasta")

    mock_write_allele_sequences.side_effect = write_side_effect

    updated = _extract_allele_sequence(sample_dict=sample_dict)

    # verify allele_sequences entries
    expected = {
        "geneA_1": os.path.join(setup_temp_dir, "geneA_1.fasta"),
        "geneB_2": os.path.join(setup_temp_dir, "geneB_2.fasta"),
    }
    assert updated["2222-FAKE-0000"]["allele_sequences"] == expected

    # Check that _write_allele_sequences was called for the correct allele ids
    calls = mock_write_allele_sequences.call_args_list
    called_ids = {call.kwargs['record'].id for call in calls}
    assert called_ids == {"geneA_1", "geneB_2"}
    assert mock_write_allele_sequences.call_count == 2


def test_extract_allele_sequence_no_kma_db(setup_temp_dir):
    """
    Test _extract_allele_sequence skips samples lacking a KMA DB.
    """
    sample_dict = {
        "2222-FAKE-0000": {
            "kma_sample_db": None,
            "baited_targets": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets.fasta"
            ),
            "best_hits": {"gene1": "geneA_1"},
            "output_path": setup_temp_dir,
        }
    }

    updated = _extract_allele_sequence(sample_dict=sample_dict)

    # no allele_sequences key added
    assert "allele_sequences" not in updated["2222-FAKE-0000"]


def test_write_allele_sequences(tmp_path):
    """
    Test _write_allele_sequences writes a FASTA file in the correct directory
    and returns the correct path.
    """
    # Prepare a mock SeqRecord
    record = SeqRecord(Seq("ATGCATGC"), id="geneX_alleleY", description="desc")

    # Prepare a mock sequence_dict
    sequence_dict = {
        "output_path": str(tmp_path)
    }

    # Call the function
    output_file = _write_allele_sequences(
        record=record, sequence_dict=sequence_dict
    )

    # The output directory should be output_path/geneX
    expected_dir = os.path.join(str(tmp_path), "geneX")
    expected_file = os.path.join(expected_dir, "geneX_alleleY.fasta")

    # Check that the returned path is correct
    assert output_file == expected_file

    # Check that the file exists
    assert os.path.exists(output_file)

    # Check that the file contains the correct FASTA record
    with open(output_file, "r", encoding="utf-8") as f:
        contents = f.read()
        assert ">geneX_alleleY" in contents
        assert "ATGCATGC" in contents


@patch("gzip.open")
def test_extract_kma_mapped_reads(mock_gzip_open, setup_temp_dir):
    """
    Test _extract_kma_mapped_reads to ensure mapped reads are extracted
    correctly.
    """
    # Prepare a mock frag file content (tab-separated fields, last two are
    # allele and read name)
    frag_lines = [
        "field1\tfield2\tgeneA_1\treadA\n",
        "field1\tfield2\tgeneB_2\treadB\n",
        "field1\tfield2\tgeneA_1\treadC\n",
        "field1\tfield2\tallele3\treadD\n",  # Not in best_hits
    ]
    # Mock gzip.open to return our lines
    mock_gzip_open.return_value.__enter__.return_value = frag_lines

    # Prepare sample_dict with best_hits
    sample_dict = {
        "2222-FAKE-0000": {
            "kma_sample_db": "mock_db",
            "kma_report_file": os.path.join(
                setup_temp_dir, "2222-FAKE-0000_baited_targets.res"
            ),
            "best_hits": {"geneA": "geneA_1", "geneB": "geneB_2"},
        }
    }

    # Patch os.path.splitext to return the correct frag file path
    with patch(
        "os.path.splitext",
        side_effect=lambda x: (x.replace(".res", ""), ".res")
    ):
        updated = _extract_kma_mapped_reads(sample_dict=sample_dict)

    mapped_reads = updated["2222-FAKE-0000"]["mapped_reads"]
    # Should only contain alleles in best_hits
    assert "geneA" in mapped_reads
    assert "geneB" in mapped_reads
    assert "geneA_1" in mapped_reads["geneA"]
    assert "geneB_2" in mapped_reads["geneB"]
    assert mapped_reads["geneA"]["geneA_1"] == ["readA", "readC"]
    assert mapped_reads["geneB"]["geneB_2"] == ["readB"]
    # allele3 should not be present
    for gene_dict in mapped_reads.values():
        assert "allele3" not in gene_dict


def test_extract_kma_mapped_reads_no_kma_db(setup_temp_dir):
    """
    Test _extract_kma_mapped_reads skips samples lacking a KMA DB.
    """
    sample_dict = {
        "2222-FAKE-0000": {
            "kma_sample_db": None,
            "kma_report_file": os.path.join(
                setup_temp_dir,
                "2222-FAKE-0000_baited_targets.res"
            ),
            "best_hits": {"geneA": "geneA_1"},
        }
    }
    updated = _extract_kma_mapped_reads(sample_dict=sample_dict)
    # mapped_reads key should not be added
    assert "mapped_reads" not in updated["2222-FAKE-0000"]


def test_write_read_names_file(tmp_path):
    """
    Test _write_read_names_file writes read names to the correct files.
    """
    # Prepare a sample_dict with mapped_reads
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "output_path": str(tmp_path),
            "mapped_reads": {
                "geneA": {
                    "geneA_1": ["read1", "read2"],
                    "geneB_2": ["read3"]
                },
                "geneB": {
                    "allele3": ["read4"]
                }
            }
        }
    }

    updated = _write_read_names_file(sample_dict=sample_dict)

    # Check that mapped_read_files and gene_output_dir are set
    mapped_files = updated["sample1"]["mapped_read_files"]
    gene_dirs = updated["sample1"]["gene_output_dir"]

    # Check files exist and contain correct read names
    for gene, allele_dict in sample_dict["sample1"]["mapped_reads"].items():
        gene_dir = os.path.join(tmp_path, gene)
        assert gene_dirs[gene] == gene_dir
        for allele, reads in allele_dict.items():
            out_file = os.path.join(gene_dir, f"{gene}_{allele}.names")
            assert mapped_files[gene][allele] == out_file
            assert os.path.exists(out_file)
            with open(out_file, "r", encoding="utf-8") as f:
                lines = [line.strip() for line in f.readlines()]
                assert lines == reads


def test_write_read_names_file_no_kma_db(tmp_path):
    """
    Test _write_read_names_file skips samples lacking a KMA DB.
    """
    sample_dict = {
        "sample1": {
            "kma_sample_db": None,
            "output_path": str(tmp_path),
            "mapped_reads": {
                "geneA": {
                    "geneA_1": ["read1"]
                }
            }
        }
    }
    updated = _write_read_names_file(sample_dict=sample_dict)
    # mapped_read_files and gene_output_dir should not be set
    assert "mapped_read_files" not in updated["sample1"]
    assert "gene_output_dir" not in updated["sample1"]


@patch("src.methods.run_command")
@patch("os.path.exists")
def test_extract_reads_by_name(mock_exists, mock_run_command, tmp_path):
    """
    Test _extract_reads_by_name runs filterbyname and updates extracted_reads.
    """
    # Prepare a sample_dict with mapped_read_files and gene_output_dir
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "baited_reads": str(tmp_path / "sample1_baited.fastq.gz"),
            "mapped_reads": {
                "geneA": {"geneA_1": ["read1", "read2"]},
                "geneB": {"geneB_2": ["read3"]}
            },
            "mapped_read_files": {
                "geneA": {
                    "geneA_1": str(tmp_path / "geneA" / "geneA_geneA_1.names")
                },
                "geneB": {
                    "geneB_2": str(tmp_path / "geneB" / "geneB_geneB_2.names")
                }
            },
            "gene_output_dir": {
                "geneA": str(tmp_path / "geneA"),
                "geneB": str(tmp_path / "geneB")
            }
        }
    }

    # Simulate that output files do not exist
    mock_exists.return_value = False

    updated = _extract_reads_by_name(sample_dict=sample_dict)

    # Check that extracted_reads is set correctly
    extracted = updated["sample1"]["extracted_reads"]
    assert "geneA" in extracted and "geneA_1" in extracted["geneA"]
    assert "geneB" in extracted and "geneB_2" in extracted["geneB"]
    assert extracted["geneA"]["geneA_1"].endswith("geneA_geneA_1.fastq.gz")
    assert extracted["geneB"]["geneB_2"].endswith("geneB_geneB_2.fastq.gz")

    # Check that run_command was called for each allele
    assert mock_run_command.call_count == 2
    for call in mock_run_command.call_args_list:
        assert "filterbyname.sh" in call.kwargs["command"]


@patch("src.methods.run_command")
def test_extract_reads_by_name_no_kma_db(mock_run_command, tmp_path):
    """
    Test _extract_reads_by_name skips samples lacking a KMA DB.
    """
    sample_dict = {
        "sample1": {
            "kma_sample_db": None,
            "baited_reads": str(tmp_path / "sample1_baited.fastq.gz"),
            "mapped_reads": {
                "geneA": {"geneA_1": ["read1"]}
            },
            "mapped_read_files": {
                "geneA": {
                    "geneA_1": str(tmp_path / "geneA" / "geneA_geneA_1.names")
                }
            },
            "gene_output_dir": {
                "geneA": str(tmp_path / "geneA")
            }
        }
    }
    updated = _extract_reads_by_name(sample_dict=sample_dict)
    # extracted_reads key should not be added
    assert "extracted_reads" not in updated["sample1"]
    mock_run_command.assert_not_called()


@patch("src.methods._index_reference_allele")
@patch("src.methods._index_bam_file")
@patch("os.path.exists")
@patch("subprocess.Popen")
def test_map_reads_bwa(
    mock_popen,
    mock_exists,
    mock_index_bam_file,
    mock_index_reference_allele,
    tmp_path
):
    """
    Test _map_reads_bwa runs BWA/Samtools pipeline and updates bwa_output.
    """
    # Prepare a sample_dict with required fields
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "baited_reads": str(tmp_path / "sample1_baited.fastq.gz"),
            "allele_sequences": {
                "geneA_1": str(tmp_path / "geneA" / "geneA_1.fasta"),
                "geneB_2": str(tmp_path / "geneB" / "geneB_2.fasta"),
            },
            "gene_output_dir": {
                "geneA": str(tmp_path / "geneA"),
                "geneB": str(tmp_path / "geneB"),
            },
            "extracted_reads": {
                "geneA": {"geneA_1": "foo"},
                "geneB": {"geneB_2": "bar"},
            }
        }
    }

    # Simulate that output BAM files do not exist
    def exists_side_effect(path):
        # Only .bam and .bwt files do not exist
        if path.endswith(".bam") or path.endswith(".bwt"):
            return False
        return True
    mock_exists.side_effect = exists_side_effect

    # Mock subprocess.Popen to simulate BWA/Samtools pipeline
    mock_proc = MagicMock()
    mock_proc.stdout = MagicMock()
    mock_proc.communicate.return_value = (b"", b"")
    mock_proc.stderr = b""
    mock_popen.return_value = mock_proc

    updated = _map_reads_bwa(sample_dict=sample_dict, threads=2)

    bwa_output = updated["sample1"]["bwa_output"]
    # Check that output BAM paths are set correctly
    assert bwa_output["geneA_1"].endswith("geneA/geneA_1.bam")
    assert bwa_output["geneB_2"].endswith("geneB/geneB_2.bam")
    # Check that _index_reference_allele and _index_bam_file were called
    assert mock_index_reference_allele.called
    assert mock_index_bam_file.called
    # Check that subprocess.Popen was called for BWA, samtools view, and sort
    assert mock_popen.call_count == 3 * 2  # 3 per allele, 2 alleles

@patch("src.methods._index_reference_allele")
@patch("src.methods._index_bam_file")
@patch("subprocess.Popen")
def test_map_reads_bwa_no_kma_db(
    mock_popen,
    mock_index_bam_file,
    mock_index_reference_allele,
    tmp_path
):
    """
    Test _map_reads_bwa skips samples lacking a KMA DB.
    """
    sample_dict = {
        "sample1": {
            "kma_sample_db": None,
            "baited_reads": str(tmp_path / "sample1_baited.fastq.gz"),
            "allele_sequences": {
                "geneA_1": str(tmp_path / "geneA" / "geneA_1.fasta"),
            },
            "gene_output_dir": {
                "geneA": str(tmp_path / "geneA"),
            },
            "extracted_reads": {
                "geneA": {"geneA_1": "foo"},
            }
        }
    }
    updated = _map_reads_bwa(sample_dict=sample_dict, threads=2)
    # bwa_output key should not be added
    assert "bwa_output" not in updated["sample1"]
    mock_popen.assert_not_called()
    mock_index_reference_allele.assert_not_called()
    mock_index_bam_file.assert_not_called()


@patch("src.methods._index_reference_allele")
@patch("src.methods._index_bam_file")
@patch("os.path.exists")
@patch("subprocess.Popen")
def test_map_reads_bwa_logs_stderr(
    mock_popen,
    mock_exists,
    _,
    __,
    tmp_path,
    caplog
):
    """
    Test _map_reads_bwa logs an error if stderr is present from samtools sort.
    """
    # Prepare a sample_dict with required fields
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "baited_reads": str(tmp_path / "sample1_baited.fastq.gz"),
            "allele_sequences": {
                "geneA_1": str(tmp_path / "geneA" / "geneA_1.fasta"),
            },
            "gene_output_dir": {
                "geneA": str(tmp_path / "geneA"),
            },
            "extracted_reads": {
                "geneA": {"geneA_1": "foo"},
            }
        }
    }

    # Simulate that output BAM and .bwt files do not exist
    def exists_side_effect(path):
        if path.endswith(".bam") or path.endswith(".bwt"):
            return False
        return True
    mock_exists.side_effect = exists_side_effect

    # Mock subprocess.Popen to simulate BWA/Samtools pipeline with stderr
    mock_proc = MagicMock()
    mock_proc.stdout = MagicMock()
    mock_proc.communicate.return_value = (b"", b"simulated error")
    mock_proc.stderr = b"simulated error"
    mock_popen.return_value = mock_proc

    with caplog.at_level("ERROR"):
        _map_reads_bwa(sample_dict=sample_dict, threads=2)
        # Check that the error was logged
        assert any(
            "Error: simulated error" in msg for msg in caplog.text.splitlines()
        )


@patch("src.methods._index_reference_allele")
@patch("src.methods._index_bam_file")
@patch("os.path.exists")
@patch("subprocess.Popen")
def test_map_reads_bwa_filenotfounderror(
    mock_popen,
    mock_exists,
    mock_index_bam_file,
    mock_index_reference_allele,
    tmp_path,
    caplog
):
    """
    Test _map_reads_bwa handles FileNotFoundError and logs the error.
    """
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "baited_reads": str(tmp_path / "sample1_baited.fastq.gz"),
            "allele_sequences": {
                "geneA_1": str(tmp_path / "geneA" / "geneA_1.fasta"),
            },
            "gene_output_dir": {
                "geneA": str(tmp_path / "geneA"),
            },
            "extracted_reads": {
                "geneA": {"geneA_1": "foo"},
            }
        }
    }

    # Simulate that output BAM and .bwt files do not exist
    def exists_side_effect(path):
        if path.endswith(".bam") or path.endswith(".bwt"):
            return False
        return True
    mock_exists.side_effect = exists_side_effect

    # Mock subprocess.Popen to raise FileNotFoundError
    mock_popen.side_effect = FileNotFoundError("bwa not found")

    with caplog.at_level("ERROR"):
        try:
            _map_reads_bwa(sample_dict=sample_dict, threads=2)
        except SystemExit:
            pass
        else:
            assert False, "SystemExit not raised on FileNotFoundError"
        assert any(
            "Error: One of the required files not found: bwa not found" in msg
            for msg in caplog.text.splitlines()
        )


@patch("src.methods._index_reference_allele")
@patch("src.methods._index_bam_file")
@patch("os.path.exists")
@patch("subprocess.Popen")
def test_map_reads_bwa_generic_exception(
    mock_popen,
    mock_exists,
    mock_index_bam_file,
    mock_index_reference_allele,
    tmp_path,
    caplog
):
    """
    Test _map_reads_bwa handles generic Exception and logs the error.
    """
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "baited_reads": str(tmp_path / "sample1_baited.fastq.gz"),
            "allele_sequences": {
                "geneA_1": str(tmp_path / "geneA" / "geneA_1.fasta"),
            },
            "gene_output_dir": {
                "geneA": str(tmp_path / "geneA"),
            },
            "extracted_reads": {
                "geneA": {"geneA_1": "foo"},
            }
        }
    }

    # Simulate that output BAM and .bwt files do not exist
    def exists_side_effect(path):
        if path.endswith(".bam") or path.endswith(".bwt"):
            return False
        return True
    mock_exists.side_effect = exists_side_effect

    # Mock subprocess.Popen to raise a generic Exception
    mock_popen.side_effect = Exception("unexpected error")

    with caplog.at_level("ERROR"):
        try:
            _map_reads_bwa(sample_dict=sample_dict, threads=2)
        except SystemExit:
            pass
        else:
            assert False, "SystemExit not raised on Exception"
        assert any(
            "An error occurred: unexpected error" in msg
            for msg in caplog.text.splitlines()
        )


@patch("src.methods.run_command")
def test_index_reference_allele_calls_run_command(mock_run_command):
    """
    Test _index_reference_allele runs the correct BWA index command.
    """
    reference_allele = "/tmp/mock_allele.fasta"
    _index_reference_allele(reference_allele=reference_allele)

    expected_command = f"bwa index {reference_allele}"
    mock_run_command.assert_called_once_with(
        command=expected_command, split=False
    )


@patch("src.methods.run_command")
def test_index_reference_allele_command_failure(mock_run_command, caplog):
    """
    Test _index_reference_allele logs and raises on command failure.
    """
    reference_allele = "/tmp/mock_allele.fasta"
    mock_run_command.side_effect = subprocess.CalledProcessError(
        1, "bwa index"
    )

    with caplog.at_level("ERROR"):
        try:
            _index_reference_allele(reference_allele=reference_allele)
        except subprocess.CalledProcessError:
            pass
        else:
            assert False, "Expected CalledProcessError to be raised"

        assert any(
            "Command failed:" in msg for msg in caplog.text.splitlines()
        )


@patch("os.path.exists")
@patch("subprocess.run")
def test_index_bam_file_runs_command(mock_run, mock_exists, tmp_path):
    """
    Test _index_bam_file runs samtools index if .bai does not exist.
    """

    bam_file = str(tmp_path / "test.bam")
    bai_file = bam_file + ".bai"

    # Simulate that the .bai file does not exist
    mock_exists.return_value = False

    # Simulate successful subprocess.run
    mock_run.return_value.returncode = 0
    mock_run.return_value.stdout = "index ok"

    _index_bam_file(bam_file=bam_file)

    mock_exists.assert_called_with(bai_file)
    mock_run.assert_called_once_with(
        ['samtools', 'index', bam_file],
        capture_output=True,
        text=True,
        check=True
    )


@patch("os.path.exists")
@patch("subprocess.run")
def test_index_bam_file_skips_if_bai_exists(mock_run, mock_exists, tmp_path):
    """
    Test _index_bam_file does nothing if .bai already exists.
    """

    bam_file = str(tmp_path / "test.bam")
    bai_file = bam_file + ".bai"

    # Simulate that the .bai file exists
    mock_exists.return_value = True

    _index_bam_file(bam_file=bam_file)

    mock_exists.assert_called_with(bai_file)
    mock_run.assert_not_called()


@patch("os.path.exists")
@patch("subprocess.run")
def test_index_bam_file_command_failure(
    mock_run,
    mock_exists,
    tmp_path,
    caplog
):
    """
    Test _index_bam_file logs and raises on samtools index failure.
    """

    bam_file = str(tmp_path / "test.bam")

    # Simulate that the .bai file does not exist
    mock_exists.return_value = False

    # Simulate subprocess.run raising CalledProcessError
    mock_run.side_effect = subprocess.CalledProcessError(1, "samtools index")

    with caplog.at_level("ERROR"):
        try:
            _index_bam_file(bam_file=bam_file)
        except subprocess.CalledProcessError:
            pass
        else:
            assert False, "Expected CalledProcessError to be raised"

        assert any(
            "Samtools index command failed:" in msg for msg in
            caplog.text.splitlines()
        )


@patch("src.methods._find_consensus_internal_insertions_and_sequence")
@patch("src.methods._translate_allele_sequence")
@patch("src.methods._detect_truncated_sequences")
def test_extract_consensus_insertions_and_sequence(
    mock_detect_trunc,
    mock_translate,
    mock_find_consensus,
    tmp_path
):
    """
    Test _extract_consensus_insertions_and_sequence populates insertions,
    consensus, aa_sequence, and aa_sequence_qa for each allele.
    """

    # Mock return values
    mock_find_consensus.return_value = ([(10, 2), (50, 1)], "ATGCATGC")
    mock_translate.return_value = "MRT*"
    mock_detect_trunc.return_value = False

    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "bwa_output": {
                "allele1": str(tmp_path / "geneA" / "allele1.bam"),
                "allele2": str(tmp_path / "geneB" / "allele2.bam"),
            }
        }
    }

    updated = _extract_consensus_insertions_and_sequence(
        sample_dict=sample_dict,
        min_coverage=0.5
    )

    for allele in ["allele1", "allele2"]:
        assert updated["sample1"]["insertions"][allele] == [(10, 2), (50, 1)]
        assert updated["sample1"]["consensus"][allele] == "ATGCATGC"
        assert updated["sample1"]["aa_sequence"][allele] == "MRT*"
        assert updated["sample1"]["aa_sequence_qa"][allele] is False

    # Check that the mocks were called for each allele
    assert mock_find_consensus.call_count == 2
    assert mock_translate.call_count == 2
    assert mock_detect_trunc.call_count == 2


def test_extract_consensus_insertions_and_sequence_no_kma_db(tmp_path):
    """
    Test _extract_consensus_insertions_and_sequence skips samples lacking a
    KMA DB.
    """
    sample_dict = {
        "sample1": {
            "kma_sample_db": None,
            "bwa_output": {
                "allele1": str(tmp_path / "geneA" / "allele1.bam"),
            }
        }
    }

    updated = _extract_consensus_insertions_and_sequence(
        sample_dict=sample_dict,
        min_coverage=0.5
    )

    # Should not add insertions, consensus, aa_sequence, or aa_sequence_qa
    assert "insertions" not in updated["sample1"]
    assert "consensus" not in updated["sample1"]
    assert "aa_sequence" not in updated["sample1"]
    assert "aa_sequence_qa" not in updated["sample1"]


@patch("pysam.AlignmentFile")
@patch("src.methods._calculate_depth_per_position")
def test_find_consensus_internal_insertions_and_sequence_basic(
    mock_calc_depth,
    mock_alignmentfile,
):
    """
    Test _find_consensus_internal_insertions_and_sequence returns expected
    insertions and consensus.
    """

    # Mock reference info
    mock_samfile = MagicMock()
    mock_samfile.lengths = [30]
    mock_samfile.references = ["ref1"]

    class FakeRead:
        query_sequence = "A" * 5
        cigartuples = [
            (pysam.CMATCH, 2),   # 2M
            (pysam.CINS, 1),     # 1I
            (pysam.CMATCH, 2),   # 2M
        ]
        reference_start = 15
        query_name = "read1"

    # Make sure fetch returns the fake read regardless of arguments
    mock_samfile.fetch.return_value = [FakeRead()]
    mock_alignmentfile.return_value = mock_samfile

    mock_calc_depth.return_value = {i: 2 for i in range(1, 31)}

    insertions, consensus = _find_consensus_internal_insertions_and_sequence(
        bam_file="dummy.bam",
        min_coverage=0.5
    )

    assert insertions == [(17, 1)]
    assert isinstance(consensus, str)
    assert len(consensus) == 30


@patch("pysam.AlignmentFile")
@patch("src.methods._calculate_depth_per_position")
def test_find_consensus_internal_insertions_and_sequence_no_insertions(
    mock_calc_depth,
    mock_alignmentfile,
    tmpdir
):
    """
    Test _find_consensus_internal_insertions_and_sequence with no insertions.
    """

    # Mock reference info
    mock_samfile = MagicMock()
    mock_samfile.lengths = [3]
    mock_samfile.references = ["ref1"]
    mock_alignmentfile.return_value = mock_samfile

    # Mock depth per position
    mock_calc_depth.return_value = {1: 2, 2: 2, 3: 2}

    # Create a fake read with a CIGAR string: 3M (all match)
    class FakeRead:
        query_sequence = "ATG"
        cigartuples = [
            (pysam.CMATCH, 3),   # 3M
        ]
        reference_start = 0
        query_name = "read1"

    mock_samfile.fetch.return_value = [FakeRead()]

    insertions, consensus = _find_consensus_internal_insertions_and_sequence(
        bam_file="dummy.bam",
        min_coverage=0.5
    )

    assert insertions == []
    assert isinstance(consensus, str)
    assert len(consensus) == 3


@patch("src.methods.pysam.AlignmentFile")
@patch("src.methods._calculate_depth_per_position")
def test_find_consensus_internal_insertions_and_sequence_min_coverage(
    mock_calc_depth,
    mock_alignmentfile,
    tmpdir
):
    mock_samfile = MagicMock()
    mock_samfile.lengths = [30]
    mock_samfile.references = ["ref1"]

    class FakeRead1:
        query_sequence = "ATGC"
        cigartuples = [
            (pysam.CMATCH, 1),
            (pysam.CINS, 1),
            (pysam.CMATCH, 2),
        ]
        reference_start = 15
        query_name = "read1"

    class FakeRead2:
        query_sequence = "ATGC"
        cigartuples = [
            (pysam.CMATCH, 4),
        ]
        reference_start = 15
        query_name = "read2"

    mock_samfile.fetch.return_value = [FakeRead1(), FakeRead2()]
    mock_alignmentfile.return_value = mock_samfile
    mock_calc_depth.return_value = {i: 2 for i in range(0, 30)}

    insertions, _ = _find_consensus_internal_insertions_and_sequence(
        bam_file="dummy.bam",
        min_coverage=0.5
    )
    assert insertions == [(16, 1)]

    insertions, _ = _find_consensus_internal_insertions_and_sequence(
        bam_file="dummy.bam",
        min_coverage=0.6
    )
    assert insertions == []


@patch("pysam.AlignmentFile")
@patch("src.methods._calculate_depth_per_position")
def test_find_consensus_internal_insertions_and_sequence_deletion(
    mock_calc_depth,
    mock_alignmentfile,
):
    """
    Test _find_consensus_internal_insertions_and_sequence handles deletions
    (CDEL).
    """

    # Mock reference info
    mock_samfile = MagicMock()
    mock_samfile.lengths = [10]
    mock_samfile.references = ["ref1"]
    mock_alignmentfile.return_value = mock_samfile

    # Mock depth per position
    mock_calc_depth.return_value = {i: 2 for i in range(10)}

    # Create a fake read with a CIGAR string: 2M 1D 2M (match, deletion, match)
    class FakeRead:
        query_sequence = "ATGC"
        cigartuples = [
            (pysam.CMATCH, 2),   # 2M
            (pysam.CDEL, 1),     # 1D
            (pysam.CMATCH, 2),   # 2M
        ]
        reference_start = 0
        query_name = "read1"

    mock_samfile.fetch.return_value = [FakeRead()]

    insertions, consensus = _find_consensus_internal_insertions_and_sequence(
        bam_file="dummy.bam",
        min_coverage=0.5
    )

    # No insertions expected, and consensus should reflect the deletion
    assert insertions == []
    assert isinstance(consensus, str)
    assert len(consensus) == 10  # Reference length


@patch("pysam.AlignmentFile")
@patch("src.methods._calculate_depth_per_position")
def test_find_consensus_internal_insertions_and_sequence_hard_clip_and_pad(
    mock_calc_depth,
    mock_alignmentfile,
):
    """
    Test _find_consensus_internal_insertions_and_sequence handles hard clips
    (CHARD_CLIP), reference skips (CREF_SKIP), and pads (CPAD).
    """

    # Mock reference info
    mock_samfile = MagicMock()
    mock_samfile.lengths = [10]
    mock_samfile.references = ["ref1"]
    mock_alignmentfile.return_value = mock_samfile

    # Mock depth per position
    mock_calc_depth.return_value = {i: 2 for i in range(10)}

    # Create a fake read with a CIGAR string: 2M 1H 1N 1P 2M (match, hard clip, skip, pad, match)
    class FakeRead:
        query_sequence = "ATGC"
        cigartuples = [
            (pysam.CMATCH, 2),       # 2M
            (pysam.CHARD_CLIP, 1),   # 1H
            (pysam.CREF_SKIP, 1),    # 1N
            (pysam.CPAD, 1),         # 1P
            (pysam.CMATCH, 2),       # 2M
        ]
        reference_start = 0
        query_name = "read1"

    mock_samfile.fetch.return_value = [FakeRead()]

    insertions, consensus = _find_consensus_internal_insertions_and_sequence(
        bam_file="dummy.bam",
        min_coverage=0.5
    )

    # No insertions expected, and consensus should reflect the reference length
    assert insertions == []
    assert isinstance(consensus, str)
    assert len(consensus) == 10  # Reference length


@patch("pysam.AlignmentFile")
@patch("src.methods._calculate_depth_per_position")
def test_find_consensus_internal_insertions_and_sequence_unexpected_cigar(
    mock_calc_depth,
    mock_alignmentfile,
):
    """
    Test _find_consensus_internal_insertions_and_sequence raises ValueError
    for unexpected CIGAR operations.
    """

    # Mock reference info
    mock_samfile = MagicMock()
    mock_samfile.lengths = [10]
    mock_samfile.references = ["ref1"]
    mock_alignmentfile.return_value = mock_samfile

    # Mock depth per position
    mock_calc_depth.return_value = {i: 2 for i in range(10)}

    # Create a fake read with an unexpected CIGAR operation
    class FakeRead:
        query_sequence = "ATGC"
        cigartuples = [
            (99, 1),  # Unexpected operation
        ]
        reference_start = 0
        query_name = "read1"

    mock_samfile.fetch.return_value = [FakeRead()]

    with pytest.raises(ValueError, match="Unexpected CIGAR operation: 99"):
        _find_consensus_internal_insertions_and_sequence(
            bam_file="dummy.bam",
            min_coverage=0.5
        )


@patch("pysam.AlignmentFile")
@patch("src.methods._calculate_depth_per_position")
def test_find_consensus_internal_insertions_and_sequence_no_reads(
    mock_calc_depth,
    mock_alignmentfile,
):
    """
    Test _find_consensus_internal_insertions_and_sequence handles positions
    with no reads.
    """

    # Mock reference info
    mock_samfile = MagicMock()
    mock_samfile.lengths = [10]
    mock_samfile.references = ["ref1"]
    mock_alignmentfile.return_value = mock_samfile

    # Mock depth per position with no reads at position 5
    mock_calc_depth.return_value = {i: 2 for i in range(10)}
    mock_calc_depth.return_value[5] = 0  # No reads at position 5

    # Create a fake read with a CIGAR string: 10M (all match)
    class FakeRead:
        query_sequence = "ATGCATGCAT"
        cigartuples = [
            (pysam.CMATCH, 10),  # 10M
        ]
        reference_start = 0
        query_name = "read1"

    mock_samfile.fetch.return_value = [FakeRead()]

    insertions, consensus = _find_consensus_internal_insertions_and_sequence(
        bam_file="dummy.bam",
        min_coverage=0.5
    )

    # No insertions expected, and consensus should reflect the reference length
    assert insertions == []
    assert isinstance(consensus, str)
    assert len(consensus) == 10  # Reference length


@patch("pysam.AlignmentFile")
def test_calculate_depth_per_position_basic(mock_alignmentfile):
    """
    Test _calculate_depth_per_position calculates depth correctly for a simple
    case.
    """

    # Mock the BAM file and pileup columns
    mock_samfile = MagicMock()
    mock_alignmentfile.return_value = mock_samfile

    # Mock pileup columns
    mock_pileup_column_1 = MagicMock()
    mock_pileup_column_1.reference_pos = 0  # 0-based position
    mock_pileup_column_1.nsegments = 5  # Depth at position 1

    mock_pileup_column_2 = MagicMock()
    mock_pileup_column_2.reference_pos = 1  # 0-based position
    mock_pileup_column_2.nsegments = 3  # Depth at position 2

    mock_samfile.pileup.return_value = [
        mock_pileup_column_1, mock_pileup_column_2
    ]

    # Call the function
    depth = _calculate_depth_per_position(
        bam_file="dummy.bam",
        ref_name="ref1",
        ref_len=10
    )

    # Assert the depth is calculated correctly
    assert depth == {1: 5, 2: 3}, \
        "Expected depth at positions 1 and 2 to be 5 and 3."


@patch("pysam.AlignmentFile")
def test_calculate_depth_per_position_no_coverage(mock_alignmentfile):
    """
    Test _calculate_depth_per_position handles positions with no coverage.
    """

    # Mock the BAM file and pileup columns
    mock_samfile = MagicMock()
    mock_alignmentfile.return_value = mock_samfile

    # Mock pileup columns (no coverage at any position)
    mock_samfile.pileup.return_value = []

    # Call the function
    depth = _calculate_depth_per_position(
        bam_file="dummy.bam",
        ref_name="ref1",
        ref_len=10
    )

    # Assert the depth is empty
    assert depth == {}, "Expected an empty depth dictionary for no coverage."


@patch("pysam.AlignmentFile")
def test_calculate_depth_per_position_partial_coverage(mock_alignmentfile):
    """
    Test _calculate_depth_per_position handles partial coverage across the
    reference.
    """

    # Mock the BAM file and pileup columns
    mock_samfile = MagicMock()
    mock_alignmentfile.return_value = mock_samfile

    # Mock pileup columns
    mock_pileup_column_1 = MagicMock()
    mock_pileup_column_1.reference_pos = 0  # 0-based position
    mock_pileup_column_1.nsegments = 5  # Depth at position 1

    mock_pileup_column_3 = MagicMock()
    mock_pileup_column_3.reference_pos = 2  # 0-based position
    mock_pileup_column_3.nsegments = 2  # Depth at position 3

    mock_samfile.pileup.return_value = [
        mock_pileup_column_1, mock_pileup_column_3
    ]

    # Call the function
    depth = _calculate_depth_per_position(
        bam_file="dummy.bam",
        ref_name="ref1",
        ref_len=10
    )

    # Assert the depth is calculated correctly
    assert depth == {1: 5, 3: 2}, \
        "Expected depth at positions 1 and 3 to be and 2."


@patch("pysam.AlignmentFile")
def test_calculate_depth_per_position_large_reference(mock_alignmentfile):
    """
    Test _calculate_depth_per_position handles a large reference sequence.
    """

    # Mock the BAM file and pileup columns
    mock_samfile = MagicMock()
    mock_alignmentfile.return_value = mock_samfile

    # Mock pileup columns for a large reference
    mock_pileup_columns = []
    for i in range(1000):  # Simulate 1000 positions
        mock_pileup_column = MagicMock()
        mock_pileup_column.reference_pos = i
        mock_pileup_column.nsegments = i % 10  # Depth varies from 0 to 9
        mock_pileup_columns.append(mock_pileup_column)

    mock_samfile.pileup.return_value = mock_pileup_columns

    # Call the function
    depth = _calculate_depth_per_position(
        bam_file="dummy.bam",
        ref_name="ref1",
        ref_len=1000
    )

    # Assert the depth is calculated correctly for a few positions
    assert depth[1] == 0, "Expected depth at position 1 to be 1."
    assert depth[10] == 9, "Expected depth at position 10 to be 9."
    assert depth[999] == 8, "Expected depth at position 999 to be 8."


def test_translate_allele_sequence_basic():
    """
    Test _translate_allele_sequence with a basic nucleotide sequence.
    """
    nt_sequence = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    expected_protein = "MAIVMGR*KGAR*"
    protein_sequence = _translate_allele_sequence(nt_sequence=nt_sequence)
    assert protein_sequence == expected_protein, \
        f"Expected {expected_protein}, but got {protein_sequence}"


def test_translate_allele_sequence_with_padding():
    """
    Test _translate_allele_sequence with a sequence length not divisible by 3.
    """
    nt_sequence = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATA"  # Length 41
    expected_protein = "MAIVMGR*KGARX"  # 'X' represents padding with 'N'
    protein_sequence = _translate_allele_sequence(nt_sequence=nt_sequence)
    assert protein_sequence == expected_protein, \
        f"Expected {expected_protein}, but got {protein_sequence}"


def test_translate_allele_sequence_empty():
    """
    Test _translate_allele_sequence with an empty sequence.
    """
    nt_sequence = ""
    expected_protein = ""
    protein_sequence = _translate_allele_sequence(nt_sequence=nt_sequence)
    assert protein_sequence == expected_protein, \
        f"Expected an empty string, but got {protein_sequence}"


def test_translate_allele_sequence_with_stop_codon():
    """
    Test _translate_allele_sequence with a sequence containing a stop codon.
    """
    nt_sequence = "ATGTAA"  # 'TAA' is a stop codon
    expected_protein = "M*"
    protein_sequence = _translate_allele_sequence(nt_sequence=nt_sequence)
    assert protein_sequence == expected_protein, \
        f"Expected {expected_protein}, but got {protein_sequence}"


def test_translate_allele_sequence_with_non_standard_bases():
    """
    Test _translate_allele_sequence with a sequence containing non-standard
    bases.
    """
    nt_sequence = "ATGNNN"  # 'N' represents any base
    expected_protein = "MX"
    protein_sequence = _translate_allele_sequence(nt_sequence=nt_sequence)
    assert protein_sequence == expected_protein, \
        f"Expected {expected_protein}, but got {protein_sequence}"


def test_detect_truncated_sequences_no_truncation():
    """
    Test _detect_truncated_sequences with a valid protein sequence of
    sufficient length.
    """
    protein_sequence = "M" * 313  # Minimum length
    result = _detect_truncated_sequences(protein_sequence=protein_sequence)
    assert result is False, "Expected no truncation for a valid sequence."


def test_detect_truncated_sequences_truncated():
    """
    Test _detect_truncated_sequences with a truncated protein sequence.
    """
    protein_sequence = "M" * 200  # Shorter than the minimum length
    result = _detect_truncated_sequences(protein_sequence=protein_sequence)
    assert result == \
        "Truncated: 200", f"Expected 'Truncated: 200', but got {result}"


def test_detect_truncated_sequences_internal_stop_codon():
    """
    Test _detect_truncated_sequences with an internal stop codon.
    """
    # Stop codon at position 101
    protein_sequence = "M" * 100 + "*" + "M" * 213
    result = _detect_truncated_sequences(protein_sequence=protein_sequence)
    assert result == "Internal stop codon at position: 101", \
        f"Expected 'Internal stop codon at position: 101', but got {result}"


def test_detect_truncated_sequences_stop_codon_after_min_length():
    """
    Test _detect_truncated_sequences with a stop codon after the minimum length
    """
    protein_sequence = "M" * 313 + "*"  # Stop codon after the minimum length
    result = _detect_truncated_sequences(protein_sequence=protein_sequence)
    assert result is False, \
        "Expected no truncation for a stop codon after the minimum length."


def test_detect_truncated_sequences_empty_sequence():
    """
    Test _detect_truncated_sequences with an empty protein sequence.
    """
    protein_sequence = ""
    result = _detect_truncated_sequences(protein_sequence=protein_sequence)
    assert result == \
        "Truncated: 0", f"Expected 'Truncated: 0', but got {result}"


def test_detect_truncated_sequences_custom_min_length():
    """
    Test _detect_truncated_sequences with a custom minimum length.
    """
    protein_sequence = "M" * 150  # Shorter than the custom minimum length
    result = _detect_truncated_sequences(
        protein_sequence=protein_sequence,
        min_length=200
    )
    assert result == \
        "Truncated: 150", f"Expected 'Truncated: 150', but got {result}"


@patch("src.methods.SeqIO.write")
@patch("os.makedirs")
def test_write_novel_allele_sequences_basic(
    mock_makedirs,
    mock_seqio_write,
    tmp_path
):
    """
    Test _write_novel_allele_sequences writes nucleotide and amino acid
    sequences for alleles with identity less than 100% but greater than or
    equal to 90%.
    """
    # Prepare a mock sample_dict
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "best_scores": {"geneA": 95.0, "geneB": 100.0, "geneC": 89.0},
            "best_hits": {
                "geneA": "geneA_1", "geneB": "geneB_2", "geneC": "geneC_3"
            },
            "consensus": {
                "geneA_1": "ATGC", "geneB_2": "CGTA", "geneC_3": "GCTA"
            },
            "aa_sequence": {
                "geneA_1": "MRT", "geneB_2": "KLM", "geneC_3": "XYZ"
            },
        }
    }

    # Call the function
    _write_novel_allele_sequences(
        report_path=str(tmp_path),
        sample_dict=sample_dict,
    )

    # Check that SeqIO.write is called for the correct sequences
    calls = mock_seqio_write.call_args_list
    assert len(calls) == 2, \
        "Expected SeqIO.write to be called twice (nucleotide and amino acid)."

    # Verify the nucleotide sequence
    nucleotide_record = calls[0].args[0]
    assert str(nucleotide_record.seq) == "ATGC", \
        "Expected nucleotide sequence to be 'ATGC'."
    assert nucleotide_record.id == "sample1_geneA", \
        "Expected nucleotide record ID to be 'sample1_gene1'."

    # Verify the amino acid sequence
    amino_acid_record = calls[1].args[0]
    assert str(amino_acid_record.seq) == "MRT", \
        "Expected amino acid sequence to be 'MRT'."
    assert amino_acid_record.id == "sample1_geneA", \
        "Expected amino acid record ID to be 'sample1_gene1'."


@patch("src.methods.SeqIO.write")
@patch("os.makedirs")
def test_write_novel_allele_sequences_no_kma_db(
    mock_makedirs,
    mock_seqio_write,
    tmp_path
):
    """
    Test _write_novel_allele_sequences skips samples where kma_sample_db is
    None.
    """
    # Prepare a mock sample_dict
    sample_dict = {
        "sample1": {
            "kma_sample_db": None,  # No KMA database
            "best_scores": {"gene1": 95.0},
            "best_hits": {"gene1": "allele1"},
            "consensus": {"allele1": "ATGC"},
            "aa_sequence": {"allele1": "MRT"},
        }
    }

    # Call the function
    _write_novel_allele_sequences(
        report_path=str(tmp_path),
        sample_dict=sample_dict,
    )

    # Ensure no directories are created and SeqIO.write is not called
    mock_makedirs.assert_not_called()
    mock_seqio_write.assert_not_called()


@patch("src.methods._calculate_stx_profile")
def test_calculate_stx_profile_basic(mock_calculate_stx_profile):
    """
    Test _calculate_stx_profile calculates stx profiles correctly for valid data.
    """

    # Prepare a mock sample_dict
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "best_hits": {"geneA": "geneA_1", "geneB": "geneB_2"},
            "best_scores": {"geneA": 95.0, "geneB": 100.0},
        },
        "sample2": {
            "kma_sample_db": "mock_db",
            "best_hits": {"geneC": "geneC_1"},
            "best_scores": {"geneC": 92.0},
        },
    }

    # Call the function
    updated_sample_dict = _calculate_stx_profile(sample_dict=sample_dict)

    # Verify the stx profiles
    assert updated_sample_dict["sample1"]["stx_profiles"] == ["geneA", "geneB"], \
        "Expected stx profiles for sample1 to include geneA and geneB."
    assert updated_sample_dict["sample2"]["stx_profiles"] == ["geneC"], \
        "Expected stx profiles for sample2 to include geneC."


def test_calculate_stx_profile_no_kma_db():
    """
    Test _calculate_stx_profile skips samples where kma_sample_db is None.
    """

    # Prepare a mock sample_dict
    sample_dict = {
        "sample1": {
            "kma_sample_db": None,  # No KMA database
            "best_hits": {"geneA": "geneA_1"},
            "best_scores": {"geneA": 95.0},
        }
    }

    # Call the function
    updated_sample_dict = _calculate_stx_profile(sample_dict=sample_dict)

    # Verify that stx_profiles is not added
    assert "stx_profiles" not in updated_sample_dict["sample1"], \
        "Expected stx_profiles to not be set when kma_sample_db is None."


def test_calculate_stx_profile_scores_below_threshold():
    """
    Test _calculate_stx_profile skips genes with scores below the threshold.
    """

    # Prepare a mock sample_dict
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "best_hits": {"geneA": "geneA_1", "geneB": "geneB_2"},
            "best_scores": {"geneA": 85.0, "geneB": 89.0},  # Below threshold
        }
    }

    # Call the function
    updated_sample_dict = _calculate_stx_profile(sample_dict=sample_dict)

    # Verify that stx_profiles is empty
    assert updated_sample_dict["sample1"]["stx_profiles"] == [], \
        "Expected stx_profiles to be empty for scores below the threshold."


def test_calculate_stx_profile_mixed_scores():
    """
    Test _calculate_stx_profile handles a mix of valid and invalid scores.
    """

    # Prepare a mock sample_dict
    # geneB below threshold
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "best_hits": {"geneA": "geneA_1", "geneB": "geneB_2"},
            "best_scores": {"geneA": 95.0, "geneB": 85.0},
        }
    }

    # Call the function
    updated_sample_dict = _calculate_stx_profile(sample_dict=sample_dict)

    # Verify that only valid scores are included in stx_profiles
    assert updated_sample_dict["sample1"]["stx_profiles"] == ["geneA"], \
        "Expected stx_profiles to include only geneA."


def test_calculate_stx_profile_empty_best_hits():
    """
    Test _calculate_stx_profile handles samples with no best hits.
    """

    # Prepare a mock sample_dict
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "best_hits": {},  # No best hits
            "best_scores": {},
        }
    }

    # Call the function
    updated_sample_dict = _calculate_stx_profile(sample_dict=sample_dict)

    # Verify that stx_profiles is empty
    assert updated_sample_dict["sample1"]["stx_profiles"] == [], \
        "Expected stx_profiles to be empty for samples with no best hits."


def test_calculate_stx_profile_multiple_samples():
    """
    Test _calculate_stx_profile handles multiple samples correctly.
    """

    # Prepare a mock sample_dict
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "best_hits": {"geneA": "geneA_1"},
            "best_scores": {"geneA": 95.0},
        },
        "sample2": {
            "kma_sample_db": "mock_db",
            "best_hits": {"geneB": "geneB_1"},
            "best_scores": {"geneB": 85.0},  # Below threshold
        },
    }

    # Call the function
    updated_sample_dict = _calculate_stx_profile(sample_dict=sample_dict)

    # Verify the stx profiles for each sample
    assert updated_sample_dict["sample1"]["stx_profiles"] == ["geneA"], \
        "Expected stx_profiles for sample1 to include geneA."
    assert updated_sample_dict["sample2"]["stx_profiles"] == [], \
        "Expected stx_profiles for sample2 to be empty."


def test_calculate_stx_profile_case_insensitivity():
    """
    Test _calculate_stx_profile ensures case insensitivity in stx profiles.
    """

    # Prepare a mock sample_dict
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "best_hits": {"GeneA": "geneA_1", "GENEB": "geneB_2"},
            "best_scores": {"GeneA": 95.0, "GENEB": 100.0},
        }
    }

    # Call the function
    updated_sample_dict = _calculate_stx_profile(sample_dict=sample_dict)

    # Verify that stx_profiles are case-insensitive
    assert updated_sample_dict["sample1"]["stx_profiles"] == [
        "GENEB", "GeneA"], \
        "Expected stx_profiles to include GeneA and GENEB."


@patch("src.methods.open", new_callable=mock_open)
def test_write_report_basic(mock_open, tmp_path):
    """
    Test _write_report writes the report file correctly for a basic case.
    """

    # Prepare a mock sample_dict
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "best_hits": {"geneA": "geneA_1", "geneB": "geneB_2"},
            "best_scores": {"geneA": 95.0, "geneB": 100.0},
            "stx_profiles": ["stx1", "stx2"],
            "insertions": {
                "geneA_1": [(10, 2), (50, 1)],
                "geneB_2": []
            },
            "aa_sequence_qa": {
                "geneA_1": "Truncated: 200",
                "geneB_2": False
            }
        }
    }

    # Call the function
    _write_report(
        report_path=str(tmp_path),
        sample_dict=sample_dict
    )

    # Check that the report file is written correctly
    combined_report_file = os.path.join(tmp_path, "stec_kma_report.tsv")
    mock_open.assert_called_once_with(
        combined_report_file, "w", encoding="utf-8")

    # Verify the content written to the file
    handle = mock_open()
    handle.write.assert_has_calls([
        call("Sample\tAllele\tPercentIdentity\tExpectedProfile\tNotes\n"),
        call(
            "sample1\tgeneA_1\t95.0\tstx1;stx2\t"
            "Insertion:Position(10):MedianLength(2);Insertion:Position(50):"
            "MedianLength(1);Truncated: 200\n"
        ),
        call("sample1\tgeneB_2\t100.0\tstx1;stx2\t\n"),
    ])


@patch("src.methods.open", new_callable=mock_open)
def test_write_report_no_kma_db(mock_open, tmp_path):
    """
    Test _write_report skips samples where kma_sample_db is None.
    """

    # Prepare a mock sample_dict
    sample_dict = {
        "sample1": {
            "kma_sample_db": None,  # No KMA database
        }
    }

    # Call the function
    _write_report(
        report_path=str(tmp_path),
        sample_dict=sample_dict
    )

    # Check that the report file is written correctly
    combined_report_file = os.path.join(tmp_path, "stec_kma_report.tsv")
    mock_open.assert_called_once_with(
        combined_report_file, "w", encoding="utf-8"
    )

    # Verify the content written to the file
    handle = mock_open()
    handle.write.assert_has_calls([
        call("Sample\tAllele\tPercentIdentity\tExpectedProfile\tNotes\n"),
        call("sample1\tND\n"),
    ])


@patch("src.methods.open", new_callable=mock_open)
def test_write_report_no_hits_above_threshold(mock_open, tmp_path):
    """
    Test _write_report skips alleles with scores below the threshold.
    """

    # Prepare a mock sample_dict
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "best_hits": {"geneA": "allele1", "geneB": "allele2"},
            "best_scores": {"geneA": 85.0, "geneB": 89.0},  # Below threshold
            "stx_profiles": ["stx1", "stx2"],
            "insertions": {},
            "aa_sequence_qa": {}
        }
    }

    # Call the function
    _write_report(
        report_path=str(tmp_path),
        sample_dict=sample_dict
    )

    # Check that the report file is written correctly
    combined_report_file = os.path.join(tmp_path, "stec_kma_report.tsv")
    mock_open.assert_called_once_with(
        combined_report_file, "w", encoding="utf-8"
    )

    # Verify the content written to the file
    handle = mock_open()
    handle.write.assert_has_calls([
        call("Sample\tAllele\tPercentIdentity\tExpectedProfile\tNotes\n"),
    ])
    # No alleles should be written since all scores are below the threshold


@patch("src.methods.open", new_callable=mock_open)
def test_write_report_with_insertions_and_qa(mock_open, tmp_path):
    """
    Test _write_report includes insertions and QA notes in the report.
    """

    # Prepare a mock sample_dict
    sample_dict = {
        "sample1": {
            "kma_sample_db": "mock_db",
            "best_hits": {"geneA": "allele1"},
            "best_scores": {"geneA": 95.0},
            "stx_profiles": ["stx1"],
            "insertions": {
                "allele1": [(10, 2), (20, 3)]
            },
            "aa_sequence_qa": {
                "allele1": "Internal stop codon at position: 50"
            }
        }
    }

    # Call the function
    _write_report(
        report_path=str(tmp_path),
        sample_dict=sample_dict
    )

    # Check that the report file is written correctly
    combined_report_file = os.path.join(tmp_path, "stec_kma_report.tsv")
    mock_open.assert_called_once_with(
        combined_report_file, "w", encoding="utf-8"
    )

    # Verify the content written to the file
    handle = mock_open()
    handle.write.assert_has_calls([
        call("Sample\tAllele\tPercentIdentity\tExpectedProfile\tNotes\n"),
        call(
            "sample1\tallele1\t95.0\tstx1\t"
            "Insertion:Position(10):MedianLength(2);Insertion:Position(20):"
            "MedianLength(3);Internal stop codon at position: 50\n"
        ),
    ])


def test_tilde_expand_with_tilde(monkeypatch):
    """
    Test _tilde_expand when the path starts with a tilde (~).
    """
    # Mock the user's home directory
    monkeypatch.setattr(
        os.path, "expanduser", lambda path: "/mock/home" + path[1:]
    )
    monkeypatch.setattr(os.path, "abspath", lambda path: path)

    # Input path with a tilde
    input_path = "~/my_dir"
    expected_path = "/mock/home/my_dir"

    # Call the function
    result = _tilde_expand(path=input_path)

    # Assert the result matches the expected expanded path
    assert result == \
        expected_path, f"Expected {expected_path}, but got {result}"


def test_tilde_expand_without_tilde(monkeypatch):
    """
    Test _tilde_expand when the path does not start with a tilde (~).
    """
    # Mock os.path.abspath to return the absolute path
    monkeypatch.setattr(os.path, "abspath", lambda path: "/absolute" + path)

    # Input path without a tilde
    input_path = "/some/dir"
    expected_path = "/absolute/some/dir"

    # Call the function
    result = _tilde_expand(path=input_path)

    # Assert the result matches the expected absolute path
    assert result == \
        expected_path, f"Expected {expected_path}, but got {result}"


def test_version():
    """
    Test that the __version__ variable is defined and matches the expected
    format.
    """
    # Ensure the version is defined
    assert __version__, "The __version__ variable should be defined."

    # Ensure the version matches the expected format (e.g., YYYY.MM.DD.X)
    version_pattern = r"^\d{4}\.\d{2}\.\d{2}\.\d+$"
    assert re.match(version_pattern, __version__), \
        f"The version '{__version__}' does not match the expected format."
