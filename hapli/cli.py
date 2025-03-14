#!/usr/bin/env python3
"""
hapli - Haplotype and Genotype Functional Annotation Tool

Main command-line interface for the hapli tool.
"""

import argparse
import sys
import os
import logging
from hapli.parsers.gfa_parser import parse_gfa
from hapli.parsers.gff_parser import parse_gff3
from hapli.analysis.sequence_analysis import (
    build_path_sequence, 
    analyze_haplotype_differences,
    analyze_haplotype_differences_with_alignment
)
from hapli.analysis.variant_analysis import (
    calculate_variant_length_changes,
    infer_variant_types,
    identify_repeat_sequences
)
from hapli.analysis.feature_analysis import (
    identify_first_cds_in_genes,
    analyze_sample_haplotypes
)
from hapli.reporting.text_report import generate_variant_effect_report
from hapli.reporting.rdf_report import (
    create_rdf_variant_report,
    create_consolidated_rdf_report,
    output_rdf_report,
    get_shex_schema
)
from hapli.utils.logging import setup_logging


def main(args=None):
    """Main function to run the variant effect report generation."""
    if args is None:
        # Parse arguments if not provided
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
    
    args = parser.parse_args()
    
    # If --save-shex is specified, save the ShEx schema and exit
    if args.save_shex:
        with open(args.save_shex, 'w') as f:
            f.write(get_shex_schema())
        logging.info(f"ShEx schema saved to {args.save_shex}")
        return 0
    
    # Setup logging
    logger = setup_logging(debug=args.debug, log_file=args.log_file, verbose=args.verbose)
    
    try:
        # Parse input files
        segments, links, paths, variants, samples, haplotypes = parse_gfa(args.gfa_file)
        features, feature_by_id, children_by_parent = parse_gff3(args.gff_file)
        
        # Mark first CDS in each gene for start codon analysis
        features = identify_first_cds_in_genes(features, feature_by_id, children_by_parent)
        
        # Calculate accurate length changes for variants
        variants = calculate_variant_length_changes(variants, segments)
        
        # Infer variant types based on segment information and length changes
        variants = infer_variant_types(variants, segments)
        
        # Analyze repetitive sequences in variants and ensure variant info is complete
        for variant in variants:
            # Analyze repetitive sequences for INS and DUP variants
            if variant.get('type') in ['DUP', 'INS']:
                repeat_info = identify_repeat_sequences(variant, segments)
                if repeat_info:
                    variant['repeat_info'] = repeat_info
            
            # Ensure we have segment associations
            if 'segments' not in variant:
                variant['segments'] = []
                
            # Make sure all variants have position information
            if 'pos' not in variant or not variant['pos']:
                # Try to infer position from segment info
                for seg_id in variant.get('segments', []):
                    if seg_id in segments and segments[seg_id].get('original_segment'):
                        orig_seg = segments[seg_id]['original_segment']
                        if orig_seg and orig_seg.startswith('S'):
                            try:
                                pos = int(orig_seg[1:])
                                variant['pos'] = pos
                                variant['end'] = pos + len(segments[seg_id]['sequence']) - 1
                                break
                            except ValueError:
                                pass
        
        # Log variant information
        for var in variants:
            logging.info(f"Variant {var.get('id', 'unknown')}: {var.get('type', 'UNKNOWN')}, length change: {var.get('length_change', 0)} bp")
        
        # Filter samples if requested
        if args.samples and args.sample_reports:
            requested_samples = [s.strip() for s in args.samples.split(',')]
            filtered_samples = {name: paths for name, paths in samples.items() if name in requested_samples}
            if not filtered_samples:
                logging.warning(f"None of the requested samples found in GFA")
            samples = filtered_samples
        
        # Select the appropriate analysis function based on the --use-alignment flag
        if args.use_alignment:
            logging.info("Using sequence alignment for variant effect analysis")
            analyze_function = analyze_haplotype_differences_with_alignment
        else:
            logging.info("Using position-based approach for variant effect analysis")
            analyze_function = analyze_haplotype_differences
        
        # Handle consolidated output if requested (only for RDF format)
        if args.consolidated and args.format == 'rdf':
            logging.info("Generating consolidated RDF report for all samples")
            
            # Create consolidated graph
            consolidated_graph = create_consolidated_rdf_report(
                variants,
                features,
                samples,
                haplotypes,
                paths,
                segments,
                feature_by_id,
                children_by_parent,
                base_uri=args.base_uri
            )
            
            # Output consolidated graph
            output_rdf_report(consolidated_graph, args.output, args.rdf_format)
            
            # Skip other report generation
            return 0
            
        # Generate sample-specific reports if requested
        if args.sample_reports:
            if samples:
                sample_reports = {}
                
                for sample_name, sample_paths in samples.items():
                    logging.info(f"Processing sample: {sample_name}")
                    
                    # Get reference path
                    ref_path_name = 'REF'
                    if ref_path_name not in paths:
                        logging.error("REF path not found in GFA")
                        continue
                    
                    # Build reference sequence
                    ref_path_segments = paths[ref_path_name]['segments']
                    ref_seq, _ = build_path_sequence(segments, ref_path_segments)
                    
                    # Determine if we have phased haplotypes for this sample
                    sample_haplotypes = haplotypes.get(sample_name, {})
                    has_phased_haplotypes = len(sample_haplotypes) >= 2
                    
                    # Collect variant segments
                    variant_segments = {}
                    for seg_id, seg_data in segments.items():
                        if seg_data.get('variant_id'):
                            variant_segments[seg_id] = seg_data
                    
                    # Prepare output file path
                    if args.output_prefix:
                        if args.format == 'text':
                            outfile = f"{args.output_prefix}_{sample_name}.txt"
                        else:  # RDF format
                            outfile = f"{args.output_prefix}_{sample_name}.{args.rdf_format}"
                    else:
                        outfile = None  # Use stdout
                    
                    if has_phased_haplotypes:
                        # Process each haplotype
                        feature_effects_by_haplotype = {}
                        
                        for hap_name, path_name in sample_haplotypes.items():
                            # Build path sequence
                            path_segments = paths[path_name]['segments']
                            hap_seq, segment_offsets = build_path_sequence(segments, path_segments)
                            
                            # Analyze differences using the selected analysis method
                            hap_effects = analyze_function(
                                ref_seq, 
                                hap_seq, 
                                features, 
                                segment_offsets, 
                                variant_segments,
                                segments
                            )
                            
                            feature_effects_by_haplotype[hap_name] = hap_effects
                        
                        # Analyze homozygous vs heterozygous effects
                        zygosity_effects = analyze_sample_haplotypes(
                            feature_effects_by_haplotype, 
                            features,
                            feature_by_id,
                            children_by_parent
                        )
                        
                        # Generate output
                        if args.format == 'text':
                            # Generate text report
                            generate_variant_effect_report([], variants, outfile, sample_name, zygosity_effects)
                            sample_reports[sample_name] = outfile
                        else:  # RDF format
                            # Create RDF graph
                            graph = create_rdf_variant_report(
                                variants,
                                features,
                                feature_effects=None,
                                sample_name=sample_name,
                                sample_effects=zygosity_effects,
                                base_uri=args.base_uri
                            )
                            # Output RDF graph
                            output_rdf_report(graph, outfile, args.rdf_format)
                            sample_reports[sample_name] = outfile
                    else:
                        # Process single haplotype
                        path_name = sample_paths[0]
                        path_segments = paths[path_name]['segments']
                        
                        # Build path sequence
                        alt_seq, segment_offsets = build_path_sequence(segments, path_segments)
                        
                        # Analyze differences using the selected analysis method
                        feature_effects = analyze_function(
                            ref_seq, 
                            alt_seq, 
                            features, 
                            segment_offsets, 
                            variant_segments,
                            segments
                        )
                        
                        # Generate output
                        if args.format == 'text':
                            # Generate text report
                            generate_variant_effect_report(feature_effects, variants, outfile, sample_name)
                            sample_reports[sample_name] = outfile
                        else:  # RDF format
                            # Create RDF graph
                            graph = create_rdf_variant_report(
                                variants,
                                features,
                                feature_effects=feature_effects,
                                sample_name=sample_name,
                                sample_effects=None,
                                base_uri=args.base_uri
                            )
                            # Output RDF graph
                            output_rdf_report(graph, outfile, args.rdf_format)
                            sample_reports[sample_name] = outfile
                
                logging.info(f"Generated reports for {len(sample_reports)} samples")
            else:
                logging.warning("No samples found in GFA, cannot generate sample reports")
        
        # Generate main report if requested
        if not args.sample_reports and not args.consolidated:
            # Get reference path
            ref_path_name = 'REF'
            if ref_path_name not in paths:
                logging.error("REF path not found in GFA")
                return 1
            
            # Build reference sequence
            ref_path_segments = paths[ref_path_name]['segments']
            ref_seq, ref_offsets = build_path_sequence(segments, ref_path_segments)
            
            # Collect variant segments
            variant_segments = {}
            for seg_id, seg_data in segments.items():
                if seg_data.get('variant_id'):
                    variant_segments[seg_id] = seg_data
            
            # Process each sample or first available sample if none specified
            processed_sample = False
            
            for sample_name, sample_haplotypes in haplotypes.items():
                if not args.samples or sample_name in args.samples.split(','):
                    logging.info(f"Processing sample: {sample_name}")
                    
                    # Process all haplotypes for compound heterozygous detection
                    if len(sample_haplotypes) >= 2:
                        feature_effects_by_haplotype = {}
                        
                        for hap_name, path_name in sample_haplotypes.items():
                            # Build path sequence
                            path_segments = paths[path_name]['segments']
                            hap_seq, segment_offsets = build_path_sequence(segments, path_segments)
                            
                            # Analyze differences using the selected analysis method
                            hap_effects = analyze_function(
                                ref_seq, 
                                hap_seq, 
                                features, 
                                segment_offsets, 
                                variant_segments,
                                segments
                            )
                            
                            feature_effects_by_haplotype[hap_name] = hap_effects
                        
                        # Analyze homozygous vs heterozygous effects
                        zygosity_effects = analyze_sample_haplotypes(
                            feature_effects_by_haplotype, 
                            features,
                            feature_by_id,
                            children_by_parent
                        )
                        
                        # Generate output
                        if args.format == 'text':
                            # Generate text report
                            generate_variant_effect_report([], variants, args.output, sample_name, zygosity_effects)
                        else:  # RDF format
                            # Create RDF graph
                            graph = create_rdf_variant_report(
                                variants,
                                features,
                                feature_effects=None,
                                sample_name=sample_name,
                                sample_effects=zygosity_effects,
                                base_uri=args.base_uri
                            )
                            # Output RDF graph
                            output_rdf_report(graph, args.output, args.rdf_format)
                    else:
                        # Process single haplotype
                        hap_name, path_name = next(iter(sample_haplotypes.items()))
                        
                        # Build path sequence
                        path_segments = paths[path_name]['segments']
                        alt_seq, segment_offsets = build_path_sequence(segments, path_segments)
                        
                        # Analyze differences using the selected analysis method
                        feature_effects = analyze_function(
                            ref_seq, 
                            alt_seq, 
                            features, 
                            segment_offsets, 
                            variant_segments,
                            segments
                        )
                        
                        # Generate output
                        if args.format == 'text':
                            # Generate text report
                            generate_variant_effect_report(feature_effects, variants, args.output, sample_name)
                        else:  # RDF format
                            # Create RDF graph
                            graph = create_rdf_variant_report(
                                variants,
                                features,
                                feature_effects=feature_effects,
                                sample_name=sample_name,
                                sample_effects=None,
                                base_uri=args.base_uri
                            )
                            # Output RDF graph
                            output_rdf_report(graph, args.output, args.rdf_format)
                    
                    processed_sample = True
                    break
            
            # If no samples processed, use the first available path
            if not processed_sample:
                for path_name, path_data in paths.items():
                    if path_name != 'REF' and 'segments' in path_data:
                        # Build path sequence
                        path_segments = path_data['segments']
                        alt_seq, segment_offsets = build_path_sequence(segments, path_segments)
                        
                        # Analyze differences using the selected analysis method
                        feature_effects = analyze_function(
                            ref_seq, 
                            alt_seq, 
                            features, 
                            segment_offsets, 
                            variant_segments,
                            segments
                        )
                        
                        # Generate output
                        if args.format == 'text':
                            # Generate text report
                            generate_variant_effect_report(feature_effects, variants, args.output)
                        else:  # RDF format
                            # Create RDF graph
                            graph = create_rdf_variant_report(
                                variants,
                                features,
                                feature_effects=feature_effects,
                                sample_name=None,
                                sample_effects=None,
                                base_uri=args.base_uri
                            )
                            # Output RDF graph
                            output_rdf_report(graph, args.output, args.rdf_format)
                        break
        
        return 0
        
    except Exception as e:
        logging.error(f"Error during report generation: {e}")
        if args.debug:
            import traceback
            logging.error(traceback.format_exc())
        return 1


if __name__ == "__main__":
    sys.exit(main())
