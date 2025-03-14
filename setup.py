from setuptools import setup, find_packages
import os
import re

# Read the version from __init__.py
with open(os.path.join('hapli', '__init__.py'), 'r') as f:
    version_file = f.read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        version = version_match.group(1)
    else:
        version = "0.1.0"  # Default if not found

# Read the long description from README.md
with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name="hapli",
    version=version,
    packages=find_packages(),
    install_requires=[
        "biopython>=1.79",
        "numpy>=1.20.0",
    ],
    extras_require={
        'rdf': ['rdflib>=6.0.0'],
        'dev': [
            'pytest>=7.0.0',
            'pytest-cov>=2.12.0',
        ],
    },
    entry_points={
        'console_scripts': [
            'hapli=hapli.cli:main',
        ],
    },
    author="Claude",
    author_email="claude@example.com",
    description="Haplotype and Genotype Functional Annotation Tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="bioinformatics, genomics, haplotype, annotation, GFA, GFF3",
    url="https://github.com/yourusername/hapli",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires=">=3.8",
    scripts=['hapli.py'],
)
