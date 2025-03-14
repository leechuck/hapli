#!/usr/bin/env python3
"""
hapli - Haplotype and Genotype Functional Annotation Tool

Main entry point for the hapli tool.
"""

import sys
import argparse
import logging
from hapli.cli import main
from hapli.tools import generate_test_data, vcf_to_gfa, gfa_to_gff3

def setup_subcommands(parser):
    """Set up subcommands for the hapli CLI."""
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Main hapli command (default)
    main_parser = subparsers.add_parser('analyze', help='Run hapli analysis')
    
    # Generate test data command
    test_data_parser = subparsers.add_parser('generate-test-data', help='Generate test data for hapli')
    test_data_parser.add_argument('--output-dir', default='testdata', help='Output directory for test files')
    test_data_parser.add_argument('--prefix', default='test', help='Prefix for output files')
    test_data_parser.add_argument('--length', type=int, default=10000, help='Length of reference sequence')
    
    # File-specific options
    test_data_parser.add_argument('--gfa', action='store_true', help='Generate GFA file')
    test_data_parser.add_argument('--gff3', action='store_true', help='Generate GFF3 file')
    test_data_parser.add_argument('--vcf', action='store_true', help='Generate VCF file')
    test_data_parser.add_argument('--fasta', action='store_true', help='Generate FASTA file')
    
    # GFA options
    test_data_parser.add_argument('--segments', type=int, default=5, help='Number of segments in GFA')
    test_data_parser.add_argument('--variants', type=int, default=10, help='Number of variants')
    test_data_parser.add_argument('--paths', type=int, default=2, help='Number of paths in GFA')
    
    # GFF3 options
    test_data_parser.add_argument('--genes', type=int, default=5, help='Number of genes in GFF3')
    test_data_parser.add_argument('--exons', type=int, default=3, help='Number of exons per gene')
    
    # VCF options
    test_data_parser.add_argument('--samples', type=int, default=2, help='Number of samples in VCF')
    test_data_parser.add_argument('--unphased', action='store_true', help='Generate unphased genotypes')
    
    # Debug and logging options
    test_data_parser.add_argument('--debug', action='store_true', help='Enable debug output')
    test_data_parser.add_argument('--verbose', action='store_true', help='Enable verbose output without full debug')
    test_data_parser.add_argument('--log-file', help='Write log to this file')
    test_data_parser.add_argument('--seed', type=int, help='Random seed for reproducible data generation')
    
    # VCF to GFA command
    vcf_to_gfa_parser = subparsers.add_parser('vcf-to-gfa', help='Convert VCF variants to GFA paths')
    vcf_to_gfa_parser.add_argument('gfa_file', help='Input GFA file')
    vcf_to_gfa_parser.add_argument('vcf_file', help='Input VCF file with variants')
    vcf_to_gfa_parser.add_argument('output_file', help='Output GFA file')
    vcf_to_gfa_parser.add_argument('--ref-path', default='REF', help='Name of the reference path in GFA')
    vcf_to_gfa_parser.add_argument('--debug', action='store_true', help='Enable debug output')
    
    # GFA to GFF3 command
    gfa_to_gff3_parser = subparsers.add_parser('gfa-to-gff3', help='Generate GFF3 file from GFA')
    gfa_to_gff3_parser.add_argument('gfa_file', help='Input GFA file')
    gfa_to_gff3_parser.add_argument('output_file', help='Output GFF3 file')
    gfa_to_gff3_parser.add_argument('--ref-path', default='REF', help='Name of the reference path in GFA')
    gfa_to_gff3_parser.add_argument('--debug', action='store_true', help='Enable debug output')
    
    return subparsers

def cli():
    """Command-line interface for hapli."""
    parser = argparse.ArgumentParser(description='hapli - Haplotype and Genotype Functional Annotation Tool')
    subparsers = setup_subcommands(parser)
    
    # If no arguments are provided, show help
    if len(sys.argv) == 1:
        parser.print_help()
        return 0
    
    args = parser.parse_args()
    
    if args.command == 'generate-test-data':
        # Pass the arguments directly to generate_test_data.main
        sys.argv[0] = 'generate_test_data'
        return generate_test_data.main()
    elif args.command == 'vcf-to-gfa':
        # Pass the arguments directly to vcf_to_gfa.main
        sys.argv[0] = 'vcf_to_gfa'
        return vcf_to_gfa.main()
    elif args.command == 'gfa-to-gff3':
        # Pass the arguments directly to gfa_to_gff3.main
        sys.argv[0] = 'gfa_to_gff3'
        return gfa_to_gff3.main()
    elif args.command == 'analyze':
        return main()
    else:
        parser.print_help()
        return 1

if __name__ == "__main__":
    sys.exit(cli())
