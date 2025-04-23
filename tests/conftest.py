#!/usr/bin/env python

"""
Shared fixtures for testing the STEC-KMA package.
"""

# Standard imports
import os
import tempfile
from typing import Generator

# Third-party imports
import pytest

__author__ = "adamkoziol"


@pytest.fixture(name="variables", scope="module")
def setup_variables():
    """
    Sets up the necessary variables for testing.

    Returns:
        Variables: An instance of the Variables class with all the necessary
        attributes set.
    """
    class Variables:
        """
        Class to hold test variables and paths.
        This includes paths to test files, mock arguments, and other
        configurations needed for testing.
        """
        def __init__(self):
            # Set up test paths
            self.test_path = os.path.abspath(os.path.dirname(__file__))
            self.file_path = os.path.join(self.test_path, "test_files")

            # Example paths for testing
            self.sequence_path = os.path.join(self.file_path, "sequences")
            self.database_path = os.path.join(self.file_path, "database")
            self.report_path = os.path.join(self.file_path, "reports")

            # Mock arguments
            self.args = {
                "sequence_path": self.sequence_path,
                "database_path": self.database_path,
                "report_path": self.report_path,
                "min_coverage": 0.7,
                "threads": 4,
                "identity": 90,
            }

    return Variables()


@pytest.fixture
def setup_temp_dir() -> Generator[str, None, None]:
    """
    Pytest fixture to create a temporary directory for each test.
    The directory is automatically cleaned up after each test.
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        yield temp_dir
