#!/usr/bin/env python3
"""
Command-line interface for variant effect reporting.
"""

import argparse
import sys
import logging

from variant_effect_report.core.utils import setup_logging
from variant_effect_report.core.analyzer import VariantEffectAnalyzer
from variant_effect_report.rdf.schemas import get_shex_schema


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Generate a variant effect report based on GFA and GFF3 files.')
    parser.add_argument('gfa_file', help='Input GFA file with variants')
    parser.add_argument('gff_file', help='Input GFF3 file with feature annotations')
    
    # Output options
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument('--output', '-o', help='Output file (default: stdout)')
    output_group.add_argument('--format', '-f', choices=['text', 'rdf'], default='text',
                         help='Output format: text (default) or RDF')
    output_group.add_argument('--rdf-format', choices=['turtle', 'n3', 'xml', 'json-ld', 'ntriples'], default='n3',
                         help='RDF serialization format (default: n3)')
    output_group.add_argument('--base-uri', default='http://example.org/genomics/',
                         help='Base URI for RDF output (default: http://example.org/genomics/)')
    output_group.add_argument('--consolidated', action='store_true',
                         help='Output a single consolidated RDF file for all samples')
    
    # Sample options
    sample_group = parser.add_argument_group('Sample Options')
    sample_group.add_argument('--sample-reports', action='store_true', 
                         help='Generate individual reports for each sample')
    sample_group.add_argument('--output-prefix', help='Prefix for sample report files')
    sample_group.add_argument('--samples', help='Comma-separated list of sample names to process (default: all)')
    
    # Analysis options
    analysis_group = parser.add_argument_group('Analysis Options')
    analysis_group.add_argument('--use-alignment', action='store_true',
                         help='Use sequence alignment for variant effect analysis')
    
    # Debug and logging options
    debug_group = parser.add_argument_group('Debug Options')
    debug_group.add_argument('--debug', action='store_true', help='Enable debug output')
    debug_group.add_argument('--verbose', action='store_true', help='Enable verbose output without full debug')
    debug_group.add_argument('--log-file', help='Write log to this file')
    debug_group.add_argument('--save-shex', help='Save ShEx validation schema to the specified file')
    
    return parser.parse_args()


def main():
    """Main function to run the variant effect report generation."""
    args = parse_args()
    
    # Setup logging
    logger = setup_logging(debug=args.debug, log_file=args.log_file, verbose=args.verbose)
    
    # If --save-shex is specified, save the ShEx schema and exit
    if args.save_shex:
        with open(args.save_shex, 'w') as f:
            f.write(get_shex_schema())
        logging.info(f"ShEx schema saved to {args.save_shex}")
        return 0
    
    try:
        # Create analyzer
        analyzer = VariantEffectAnalyzer(use_alignment=args.use_alignment)
        
        # Load data
        analyzer.load_data(args.gfa_file, args.gff_file)
        
        # Process samples
        if args.samples:
            sample_list = [s.strip() for s in args.samples.split(',')]
            for sample_name in sample_list:
                if args.output_prefix:
                    output = f"{args.output_prefix}_{sample_name}"
                else:
                    output = args.output
                
                analyzer.generate_report(
                    sample_name=sample_name,
                    output=output,
                    format=args.format,
                    rdf_format=args.rdf_format,
                    base_uri=args.base_uri
                )
        elif args.sample_reports:
            # Generate reports for all samples
            analyzer.generate_report(
                sample_name=None,  # Process all samples
                output=args.output_prefix,
                format=args.format,
                rdf_format=args.rdf_format,
                base_uri=args.base_uri
            )
        else:
            # Generate a single report
            analyzer.generate_report(
                sample_name=None,  # Use first available sample
                output=args.output,
                format=args.format,
                rdf_format=args.rdf_format,
                base_uri=args.base_uri
            )
        
        return 0
    
    except Exception as e:
        logging.error(f"Error during report generation: {e}")
        if args.debug:
            import traceback
            logging.error(traceback.format_exc())
        return 1


if __name__ == "__main__":
    sys.exit(main())
