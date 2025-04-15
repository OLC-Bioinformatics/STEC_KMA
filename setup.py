#!/usr/bin/env python3

"""
Setup script for the STEC_KMA package.
"""

# Standard imports
import os
import re
from setuptools import setup, find_packages

__author__ = 'adamkoziol'
__email__ = 'adam.koziol@inspection.gc.ca'

# Find the version
version = {}
version_script = os.path.join('src', 'version.py')

# Open the version script
with open(version_script, 'r', encoding='utf-8') as version_file:
    version_content = version_file.read()
    version_match = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]",
        version_content,
        re.M
    )
    if version_match:
        version['__version__'] = version_match.group(1)
    else:
        raise RuntimeError(
            f"Unable to find version string in {version_script}.")

setup(
    name="STEC_KMA",
    version=version['__version__'],
    packages=find_packages(),
    include_package_data=True,
    scripts=[
        os.path.join('src', 'stec_kma.py'),
    ],
    license='MIT',
    author=__author__,
    author_email=__email__,
    description='Description of your STEC_KMA tool',
    url='https://github.com/OLC-LOC-Bioinformatics/STEC_KMA',
    long_description=open('README.md', encoding='utf-8').read(),
    long_description_content_type='text/markdown',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
)
