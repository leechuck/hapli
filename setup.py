from setuptools import setup, find_packages

setup(
    name="hapli",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.79",
        "rdflib>=6.0.0",
    ],
    entry_points={
        'console_scripts': [
            'hapli=hapli.cli:main',
        ],
    },
    author="Claude",
    author_email="claude@example.com",
    description="Haplotype and Genotype Functional Annotation Tool",
    keywords="bioinformatics, genomics, haplotype, annotation, GFA, GFF3",
    url="https://github.com/yourusername/hapli",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires=">=3.6",
)
