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
    test_data_parser.add_argument('--debug', action='store_true', help='Enable debug output')
    
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
    setup_subcommands(parser)
    
    args = parser.parse_args()
    
    if args.command == 'generate-test-data':
        return generate_test_data.main()
    elif args.command == 'vcf-to-gfa':
        return vcf_to_gfa.main()
    elif args.command == 'gfa-to-gff3':
        return gfa_to_gff3.main()
    else:
        # Default to main hapli functionality
        return main()

if __name__ == "__main__":
    sys.exit(cli())
