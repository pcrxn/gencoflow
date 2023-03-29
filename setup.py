#!/usr/bin/env python
"""
Set up the package
"""

# Standard inputs
from distutils.util import convert_path
import os

# Third party inputs
from setuptools import setup, find_packages

# Find the version
version = {}
with open(convert_path(os.path.join('src', 'version.py')), 'r') as version_file:
    exec(version_file.read(), version)

setup(
    name="GenCoFlow",
    version=version['__version__'],
    scripts=[
        os.path.join('src', 'create_windows.py'),
        os.path.join('src', 'parse_gbk.py'),
        os.path.join('src', 'gencoflow.py')
    ],
    packages=find_packages(),
    include_package_data=True,
    author="Liam Brown",
    author_email="liam.brown@inspection.gc.ca",
    url="https://github.com/pcrxn/gencoflow",
)