#!/usr/bin/env python

"""
Unit tests for the CLI of the STEC_KMA package.
"""

# Standard imports
import os
import sys
from unittest.mock import patch

# Third-party imports
import pytest

# Local imports
from src.stec_kma import (
    cli,
    restricted_float
)


def test_cli_returns_default_arguments():
    """
    Test that cli() correctly parses and returns arguments with default values.
    """
    test_args = [
        "stec_kma.py",
        "--sequence_path", "/path/to/sequences",
        "--database_path", "/path/to/database",
        "--report_path", "/path/to/reports"
    ]

    with patch.object(sys, 'argv', test_args):
        args = cli()

        # Verify arguments are parsed correctly
        assert args.sequence_path == "/path/to/sequences"
        assert args.database_path == "/path/to/database"
        assert args.report_path == "/path/to/reports"
        assert args.min_coverage == 0.7  # default value
        assert args.threads == os.cpu_count()  # default value
        assert args.identity == 90.0  # default value


def test_cli_returns_custom_arguments():
    """
    Test that cli() correctly parses and returns custom arguments.
    """
    test_args = [
        "stec_kma.py",
        "--sequence_path", "/path/to/sequences",
        "--database_path", "/path/to/database",
        "--report_path", "/path/to/reports",
        "--min_coverage", "0.75",
        "--identity", "95"
    ]

    with patch.object(sys, 'argv', test_args):
        args = cli()

        # Verify arguments are parsed correctly
        assert args.sequence_path == "/path/to/sequences"
        assert args.database_path == "/path/to/database"
        assert args.report_path == "/path/to/reports"
        assert args.min_coverage == 0.75
        assert args.identity == "95"


def test_cli_handles_missing_report_path():
    """
    Test that cli() correctly handles when report_path is not provided.
    """
    test_args = [
        "stec_kma.py",
        "--sequence_path", "/path/to/sequences",
        "--database_path", "/path/to/database"
    ]

    with patch.object(sys, 'argv', test_args):
        args = cli()

        # Report path should be None when not provided
        assert args.sequence_path == "/path/to/sequences"
        assert args.database_path == "/path/to/database"
        assert args.report_path is None


def test_cli_with_invalid_arguments():
    """
    Test that cli() handles invalid arguments appropriately.
    """
    test_args = [
        "stec_kma.py",
        "--sequence_path", "/path/to/sequences",
        # Missing required database_path
    ]

    with patch.object(sys, 'argv', test_args):
        with pytest.raises(SystemExit):
            cli()


def test_cli_with_invalid_min_coverage():
    """
    Test that cli() validates the min_coverage value.
    """
    test_args = [
        "stec_kma.py",
        "--sequence_path", "/path/to/sequences",
        "--database_path", "/path/to/database",
        "--bin_coverage", "0.7"  # Invalid
    ]

    with patch.object(sys, 'argv', test_args):

        # Check the argparse validation
        with pytest.raises(SystemExit):
            cli()


def test_restricted_float_valid():
    """
    Test that restricted_float() correctly converts valid string inputs
    """
    assert restricted_float("0.0") == 0.0
    assert restricted_float("1.0") == 1.0
    assert restricted_float("0.5") == 0.5


def test_restricted_float_invalid_valueerror():
    """
    Test that restricted_float() raises ValueError for invalid string inputs
    """
    with pytest.raises(ValueError, match="not_a_float is not a valid float"):
        restricted_float("not_a_float")


def test_restricted_float_out_of_range_low():
    """
    Test that restricted_float() raises ValueError for out-of-range low values
    """
    with pytest.raises(ValueError, match="-0.1 is not between 0 and 1"):
        restricted_float("-0.1")


def test_restricted_float_out_of_range_high():
    """
    Test that restricted_float() raises ValueError for out-of-range high values
    """
    with pytest.raises(ValueError, match="1.1 is not between 0 and 1"):
        restricted_float("1.1")
