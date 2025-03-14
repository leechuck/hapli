#!/usr/bin/env python3
"""
Setup script for variant-effect-report package.
"""

from setuptools import setup, find_packages

setup(
    name="variant-effect-report",
    version="0.1.0",
    description="A tool for analyzing variants in GFA files and reporting their effects on genomic features",
    author="Variant Effect Report Team",
    author_email="example@example.com",
    url="https://github.com/example/variant-effect-report",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.79",
        "rdflib>=6.0.0",
    ],
    entry_points={
        "console_scripts": [
            "variant-effect-report=variant_effect_report.cli:main",
            "generate-test-data=variant_effect_report.tools.generate_test_data:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.7",
)
