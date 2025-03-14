#!/usr/bin/env python3
"""
Generate Test Data for Variant Effect Analysis Testing
Creates a simple GFA graph and VCF file with structural variants
"""

import argparse
import os
import sys
import logging
from typing import List, Dict, Tuple, Optional

from variant_effect_report.core.test_data import (
    SimpleGFAGenerator, 
    load_reference_sequence,
    parse_var_types,
    parse_size_range
)
from variant_effect_report.core.utils import setup_logging


def main():
    """Main function to generate test data files"""
    parser = argparse.ArgumentParser(description='Generate test data for variant effect analysis')
    
    # Input/output options
    io_group = parser.add_argument_group('Input/Output Options')
    io_group.add_argument('--output-dir', '-o', type=str, default='.',
                        help='Output directory for generated files (default: current directory)')
    io_group.add_argument('--prefix', '-p', type=str, default='test',
                        help='Prefix for output filenames (default: test)')
    io_group.add_argument('--gfa', action='store_true',
                        help='Generate GFA file')
    io_group.add_argument('--vcf', action='store_true',
                        help='Generate VCF file')
    io_group.add_argument('--fasta', action='store_true',
                        help='Generate FASTA file with reference sequence')
    io_group.add_argument('--reference', '-r', type=str,
                        help='Input reference FASTA file (.fa or .fa.gz) to use instead of generating random sequence')
    
    # Sequence/graph options
    seq_group = parser.add_argument_group('Sequence and Graph Options')
    seq_group.add_argument('--seq-length', '-l', type=int, default=1000,
                        help='Length of reference sequence (ignored if --reference is provided) (default: 1000)')
    seq_group.add_argument('--node-size', '-n', type=int, default=100,
                        help='Average node size (default: 100)')
    seq_group.add_argument('--contig-id', type=str, default='1',
                        help='Contig ID to use in VCF (default: 1)')
    seq_group.add_argument('--hgvs-prefix', type=str, default='NC_000001.11',
                        help='Prefix to use for HGVS notation (default: NC_000001.11)')
    
    # Variant options
    var_group = parser.add_argument_group('Variant Options')
    var_group.add_argument('--num-variants', '-v', type=int, default=5,
                        help='Number of variants to generate (default: 5)')
    var_group.add_argument('--var-types', type=str, default='SNP,DEL,INS,DUP,INV',
                        help='Comma-separated list of variant types to generate (default: SNP,DEL,INS,DUP,INV)')
    var_group.add_argument('--size-ranges', type=str, default='',
                        help='Size ranges for structural variants in format "TYPE:MIN-MAX,TYPE:MIN-MAX" (e.g., DEL:1-50,INS:1-20)')
    var_group.add_argument('--seed', '-s', type=int, default=42,
                        help='Random seed (default: 42)')
    
    # Sample options
    sample_group = parser.add_argument_group('Sample Options')
    sample_group.add_argument('--samples', type=str, default='',
                        help='Comma-separated list of sample names for multi-sample VCF (default: no samples)')
    sample_group.add_argument('--num-samples', type=int, default=0,
                        help='Number of samples to generate with auto-named IDs (e.g., Sample1, Sample2)')
    sample_group.add_argument('--phased', action='store_true',
                        help='Generate phased genotypes (|) instead of unphased (/) in the VCF')
    
    # Debug options
    debug_group = parser.add_argument_group('Debug Options')
    debug_group.add_argument('--debug', action='store_true', help='Enable debug output')
    debug_group.add_argument('--verbose', action='store_true', help='Enable verbose output without full debug')
    debug_group.add_argument('--log-file', help='Write log to this file')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(debug=args.debug, log_file=args.log_file, verbose=args.verbose)
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Determine which files to generate
    generate_all = not (args.gfa or args.vcf or args.fasta)
    
    # Parse variant types and size ranges
    try:
        var_types = parse_var_types(args.var_types)
    except ValueError as e:
        parser.error(str(e))
        return 1
    
    size_ranges = parse_size_range(args.size_ranges)
    
    # Create file paths
    gfa_path = os.path.join(args.output_dir, f"{args.prefix}.gfa")
    vcf_path = os.path.join(args.output_dir, f"{args.prefix}.vcf")
    fasta_path = os.path.join(args.output_dir, f"{args.prefix}.fa")
    
    # Load or generate reference sequence
    ref_seq = None
    sequence_id = "REF"
    
    if args.reference:
        try:
            logging.info(f"Loading reference sequence from {args.reference}")
            sequence_id, ref_seq = load_reference_sequence(args.reference)
            logging.info(f"Loaded reference sequence '{sequence_id}' of length {len(ref_seq)}")
        except Exception as e:
            parser.error(f"Error loading reference sequence: {str(e)}")
            return 1
    else:
        logging.info(f"Generating random reference sequence of length {args.seq_length}")
    
    # Prepare sample list for VCF
    samples = []
    if args.samples:
        samples = [s.strip() for s in args.samples.split(',') if s.strip()]
    
    if args.num_samples > 0:
        # Add auto-generated sample names
        existing_count = len(samples)
        for i in range(1, args.num_samples + 1):
            samples.append(f"Sample{existing_count + i}")
    
    if samples:
        logging.info(f"VCF will include {len(samples)} samples: {', '.join(samples)}")
        if args.phased:
            logging.info("Genotypes will be phased")
        else:
            logging.info("Genotypes will be unphased")
    
    # Create generator with reference or random sequence
    generator = SimpleGFAGenerator(
        ref_seq=ref_seq,
        seq_length=args.seq_length,
        node_size=args.node_size,
        seed=args.seed,
        sequence_id=sequence_id
    )
    
    # Generate GFA if requested or if generating all
    if args.gfa or generate_all:
        logging.info(f"Exporting GFA to {gfa_path}")
        generator.export_gfa(gfa_path)
    
    # Generate FASTA if requested or if generating all
    if args.fasta or generate_all:
        logging.info(f"Exporting reference FASTA to {fasta_path}")
        generator.export_fasta(fasta_path)
    
    # Generate VCF if requested or if generating all
    if args.vcf or generate_all:
        logging.info(f"Generating {args.num_variants} structural variants of types: {', '.join(var_types)}")
        
        if size_ranges:
            size_info = ", ".join([f"{k}:{v[0]}-{v[1]}" for k, v in size_ranges.items()])
            logging.info(f"Using custom size ranges: {size_info}")
        
        variants = generator.generate_vcf_variants(
            num_variants=args.num_variants,
            var_types=var_types,
            size_ranges=size_ranges,
            hgvs_prefix=args.hgvs_prefix
        )
        
        logging.info(f"Exporting VCF to {vcf_path}")
        generator.export_vcf(variants, vcf_path, args.contig_id, samples, args.phased)
        
        # Report details of variants
        logging.info("\nGenerated variants:")
        for i, var in enumerate(variants):
            logging.info(f"{i+1}. {var['type']} at position {var['pos']}: {var['hgvs']}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
