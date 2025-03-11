#!/usr/bin/env python
"""
Generate a GFF3 file based on GFA paths and VCF variants.
This script takes a GFA file and optionally a VCF file, builds
the reference sequence, and generates synthetic gene features
for demonstration purposes.

Author: Claude
Date: 2025-03-11
"""

import argparse
import re
import sys
import os
import logging
import random
from collections import defaultdict
import time
from Bio import SeqIO
from Bio.Seq import Seq

def setup_logging(debug=False, log_file=None, verbose=False):
    """Configure logging based on debug flag and optional log file."""
    if debug:
        log_level = logging.DEBUG
        log_format = '%(asctime)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s'
    elif verbose:
        log_level = logging.INFO
        log_format = '%(asctime)s - %(levelname)s - %(message)s'
    else:
        log_level = logging.INFO
        log_format = '%(message)s'
    
    # Configure root logger
    logging.basicConfig(level=log_level, format=log_format)
    logger = logging.getLogger()
    
    # Clear any existing handlers
    logger.handlers = []
    
    # Add console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_formatter = logging.Formatter(log_format)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    # Add file handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s')
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    
    return logger

def parse_gfa(gfa_file, validate=True):
    """Parse a GFA file into segments, links, and paths."""
    segments = {}
    links = []
    paths = {}
    
    start_time = time.time()
    logging.info(f"Parsing GFA file: {gfa_file}")
    
    with open(gfa_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            record_type = fields[0]
            
            if record_type == 'S':  # Segment
                if len(fields) < 3:
                    logging.warning(f"Line {line_num}: Invalid segment record, missing fields")
                    continue
                    
                seg_id = fields[1]
                sequence = fields[2]
                
                # Extract optional tags if present
                tags = {}
                for field in fields[3:]:
                    if ':' in field:
                        tag_parts = field.split(':', 2)
                        if len(tag_parts) >= 3:
                            tag_name, tag_type, tag_value = tag_parts
                            tags[tag_name] = (tag_type, tag_value)
                
                segments[seg_id] = {
                    'sequence': sequence,
                    'length': len(sequence),
                    'tags': tags
                }
                
                logging.debug(f"Parsed segment: {seg_id} (length: {len(sequence)}, tags: {len(tags)})")
                
            elif record_type == 'L':  # Link
                if len(fields) < 6:
                    logging.warning(f"Line {line_num}: Invalid link record, missing fields")
                    continue
                    
                from_id = fields[1]
                from_dir = fields[2]
                to_id = fields[3]
                to_dir = fields[4]
                overlap = fields[5]
                
                links.append((from_id, from_dir, to_id, to_dir, overlap))
                logging.debug(f"Parsed link: {from_id}{from_dir} -> {to_id}{to_dir}")
                
            elif record_type == 'P':  # Path
                if len(fields) < 3:
                    logging.warning(f"Line {line_num}: Invalid path record, missing fields")
                    continue
                    
                path_name = fields[1]
                path_segments = []
                for seg in fields[2].split(','):
                    # Split into segment ID and orientation
                    if not seg or len(seg) < 2:
                        logging.warning(f"Line {line_num}: Invalid segment in path: {seg}")
                        continue
                        
                    seg_id = seg[:-1]
                    orientation = seg[-1]
                    
                    # Validate segment exists
                    if validate and seg_id not in segments:
                        logging.warning(f"Line {line_num}: Path {path_name} references unknown segment: {seg_id}")
                    
                    path_segments.append((seg_id, orientation))
                paths[path_name] = path_segments
                logging.debug(f"Parsed path: {path_name} with {len(path_segments)} segments")
    
    elapsed = time.time() - start_time
    logging.info(f"Finished parsing GFA in {elapsed:.2f}s: {len(segments)} segments, {len(links)} links, {len(paths)} paths")
    
    return segments, links, paths

def parse_vcf(vcf_file):
    """Parse a VCF file to extract variant information."""
    variants = []
    
    start_time = time.time()
    logging.info(f"Parsing VCF file: {vcf_file}")
    
    with open(vcf_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 8:
                logging.warning(f"Line {line_num}: Invalid VCF record, missing fields")
                continue
                
            chrom = fields[0]
            try:
                pos = int(fields[1])
            except ValueError:
                logging.warning(f"Line {line_num}: Invalid position: {fields[1]}")
                continue
                
            var_id = fields[2]
            ref = fields[3]
            alt = fields[4]
            
            # Parse INFO field
            info_dict = {}
            for item in fields[7].split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
                else:
                    info_dict[item] = True
            
            # Determine variant type
            variant_type = info_dict.get('SVTYPE', 'UNKNOWN')
            
            # Extract END position
            try:
                end_pos = int(info_dict.get('END', pos + len(ref) - 1))
            except ValueError:
                logging.warning(f"Line {line_num}: Invalid END position: {info_dict.get('END')}")
                end_pos = pos + len(ref) - 1
            
            # Extract SVLEN if available
            try:
                svlen = int(info_dict.get('SVLEN', len(alt) - len(ref)))
            except ValueError:
                logging.warning(f"Line {line_num}: Invalid SVLEN: {info_dict.get('SVLEN')}")
                svlen = len(alt) - len(ref)
            
            variant = {
                'id': var_id,
                'chrom': chrom,
                'pos': pos,
                'end': end_pos,
                'ref': ref,
                'alt': alt,
                'type': variant_type,
                'svlen': svlen,
                'info': info_dict
            }
            
            variants.append(variant)
            logging.debug(f"Parsed variant: {var_id} {variant_type} at position {pos}-{end_pos}")
    
    elapsed = time.time() - start_time
    logging.info(f"Finished parsing VCF in {elapsed:.2f}s: {len(variants)} variants")
    
    return variants

def build_reference_sequence(segments, path_segments):
    """Build a reference sequence from segments and calculate offsets."""
    reference_seq = ""
    segment_offsets = {}
    current_offset = 0
    
    for seg_id, orientation in path_segments:
        if seg_id not in segments:
            logging.error(f"Missing segment: {seg_id}")
            continue
            
        segment_seq = segments[seg_id]['sequence']
        
        # Handle reverse complement if needed
        if orientation == '-':
            seq_obj = Seq(segment_seq)
            segment_seq = str(seq_obj.reverse_complement())
            
        segment_offsets[seg_id] = {
            'start': current_offset,
            'end': current_offset + len(segment_seq) - 1,
            'orientation': orientation
        }
        reference_seq += segment_seq
        current_offset += len(segment_seq)
    
    logging.info(f"Built reference sequence of length {len(reference_seq)}bp from {len(path_segments)} segments")
    return reference_seq, segment_offsets

def create_synthetic_features(reference_seq, num_genes=5, num_exons_per_gene=3, 
                             min_gene_length=200, max_gene_length=500):
    """Create synthetic gene features for demonstration purposes."""
    
    seq_length = len(reference_seq)
    if seq_length < min_gene_length:
        logging.error(f"Reference sequence too short ({seq_length}bp) for generating features")
        return []
    
    logging.info(f"Generating {num_genes} synthetic genes with {num_exons_per_gene} exons each")
    
    # Calculate boundaries for placing genes
    usable_length = seq_length - max_gene_length
    if usable_length <= 0:
        logging.warning(f"Sequence too short for max gene length, reducing gene length")
        max_gene_length = int(seq_length * 0.8)
        min_gene_length = min(min_gene_length, int(max_gene_length * 0.5))
        usable_length = seq_length - max_gene_length
    
    # Calculate gene start positions to distribute across sequence
    gene_starts = []
    for i in range(num_genes):
        # Calculate target position based on even distribution
        target_pos = int(usable_length * i / (num_genes - 1 if num_genes > 1 else 1))
        
        # Add some randomness (+/- 10% of usable length)
        if num_genes > 1:  # Add randomness only if multiple genes
            jitter = int(usable_length * 0.1)
            jitter_amount = random.randint(-jitter, jitter)
            target_pos = max(0, min(usable_length, target_pos + jitter_amount))
            
        gene_starts.append(target_pos)
    
    # Sort gene starts to ensure they're in order
    gene_starts.sort()
    
    features = []
    
    # Generate genes and their sub-features
    gene_id = 1
    for start_pos in gene_starts:
        # Determine gene length
        gene_length = random.randint(min_gene_length, max_gene_length)
        end_pos = start_pos + gene_length
        
        # Ensure gene doesn't exceed sequence length
        if end_pos >= seq_length:
            end_pos = seq_length - 1
            gene_length = end_pos - start_pos + 1
        
        # Skip if gene is too short after adjustment
        if gene_length < min_gene_length:
            continue
        
        gene_name = f"gene_{gene_id}"
        gene_id_str = f"gene{gene_id}"
        strand = "+" if random.random() > 0.3 else "-"  # 70% on + strand
        
        # Create gene feature
        gene_feature = {
            'type': 'gene',
            'start': start_pos + 1,  # Convert to 1-based
            'end': end_pos + 1,      # Convert to 1-based
            'strand': strand,
            'attributes': {
                'ID': gene_id_str,
                'Name': gene_name,
                'locus_tag': f"LOCUS{gene_id:04d}"
            }
        }
        
        features.append(gene_feature)
        
        # Create mRNA feature
        mrna_id = f"mRNA{gene_id}"
        mrna_feature = {
            'type': 'mRNA',
            'start': gene_feature['start'],
            'end': gene_feature['end'],
            'strand': strand,
            'attributes': {
                'ID': mrna_id,
                'Parent': gene_id_str,
                'Name': f"{gene_name}_transcript"
            }
        }
        
        features.append(mrna_feature)
        
        # Create exons and CDS features
        # First, decide how to divide the gene into exons
        exon_boundaries = []
        gene_size = gene_feature['end'] - gene_feature['start'] + 1
        
        if num_exons_per_gene == 1:
            # Single exon gene
            exon_boundaries = [(gene_feature['start'], gene_feature['end'])]
        else:
            # Multi-exon gene - generate random cut points
            cuts = sorted([random.randint(gene_feature['start'], gene_feature['end']) 
                          for _ in range(num_exons_per_gene - 1)])
            
            # Create exon boundaries from cuts
            exon_start = gene_feature['start']
            for cut in cuts:
                exon_boundaries.append((exon_start, cut))
                exon_start = cut + 1
            
            # Add the final exon
            exon_boundaries.append((exon_start, gene_feature['end']))
            
            # Ensure no exons are too small (minimum 15 bp)
            min_exon_size = 15
            valid_boundaries = []
            for start, end in exon_boundaries:
                if end - start + 1 >= min_exon_size:
                    valid_boundaries.append((start, end))
            
            exon_boundaries = valid_boundaries
        
        # Create exon and CDS features
        for i, (exon_start, exon_end) in enumerate(exon_boundaries, 1):
            exon_id = f"exon{gene_id}.{i}"
            
            # Exon feature
            exon_feature = {
                'type': 'exon',
                'start': exon_start,
                'end': exon_end,
                'strand': strand,
                'attributes': {
                    'ID': exon_id,
                    'Parent': mrna_id,
                    'Name': f"{gene_name}_exon{i}"
                }
            }
            
            features.append(exon_feature)
            
            # CDS feature (coding sequence)
            # For simplicity, we'll make the CDS the same as the exon
            cds_feature = {
                'type': 'CDS',
                'start': exon_start,
                'end': exon_end,
                'strand': strand,
                'attributes': {
                    'ID': f"cds{gene_id}.{i}",
                    'Parent': mrna_id,
                    'Name': f"{gene_name}_cds{i}"
                }
            }
            
            features.append(cds_feature)
        
        gene_id += 1
    
    # Add some regulatory features (promoters, terminators)
    for gene_feature in [f for f in features if f['type'] == 'gene']:
        gene_id = gene_feature['attributes']['ID']
        strand = gene_feature['strand']
        
        # Promoter (upstream of gene)
        promoter_length = random.randint(50, 150)
        if strand == '+':
            promoter_start = max(1, gene_feature['start'] - promoter_length)
            promoter_end = gene_feature['start'] - 1
        else:
            promoter_start = gene_feature['end'] + 1
            promoter_end = min(len(reference_seq), gene_feature['end'] + promoter_length)
        
        # Only add if we have space
        if promoter_end >= promoter_start:
            promoter_feature = {
                'type': 'promoter',
                'start': promoter_start,
                'end': promoter_end,
                'strand': strand,
                'attributes': {
                    'ID': f"promoter_{gene_id}",
                    'Name': f"promoter_{gene_id}",
                    'Associated_gene': gene_id
                }
            }
            features.append(promoter_feature)
        
        # Terminator (downstream of gene)
        terminator_length = random.randint(30, 80)
        if strand == '+':
            terminator_start = gene_feature['end'] + 1
            terminator_end = min(len(reference_seq), gene_feature['end'] + terminator_length)
        else:
            terminator_start = max(1, gene_feature['start'] - terminator_length)
            terminator_end = gene_feature['start'] - 1
        
        # Only add if we have space
        if terminator_end >= terminator_start:
            terminator_feature = {
                'type': 'terminator',
                'start': terminator_start,
                'end': terminator_end,
                'strand': strand,
                'attributes': {
                    'ID': f"terminator_{gene_id}",
                    'Name': f"terminator_{gene_id}",
                    'Associated_gene': gene_id
                }
            }
            features.append(terminator_feature)
    
    logging.info(f"Generated {len(features)} synthetic features")
    return features

def create_variant_features(variants, reference_length):
    """Create features from variant information."""
    features = []
    
    for variant in variants:
        var_type = 'variation'
        
        # Convert between VCF variant types and GFF3 feature types
        if variant['type'] == 'DEL':
            var_type = 'deletion'
        elif variant['type'] == 'INS':
            var_type = 'insertion'
        elif variant['type'] == 'INV':
            var_type = 'inversion'
        elif variant['type'] == 'SNP':
            var_type = 'SNV'
        
        # Create feature
        feature = {
            'type': var_type,
            'start': variant['pos'],  # VCF positions are 1-based
            'end': variant['end'],
            'strand': '+',  # Variants in VCF don't have strand information
            'attributes': {
                'ID': variant['id'],
                'Name': variant['id'],
                'Reference_seq': variant['ref'],
                'Alternate_seq': variant['alt'],
                'Variant_type': variant['type']
            }
        }
        
        # Add additional info attributes
        for key, value in variant['info'].items():
            if key not in ['SVTYPE', 'END'] and isinstance(value, str):
                # Clean attribute name and value
                attr_name = key.replace(' ', '_')
                attr_value = value.replace(' ', '_').replace(';', ',')
                feature['attributes'][attr_name] = attr_value
        
        features.append(feature)
    
    logging.info(f"Created {len(features)} variant features")
    return features

def apply_variant_effects(features, variants):
    """Add variant effect annotations to features that overlap with variants."""
    # Create a map of positions to variants
    variant_positions = defaultdict(list)
    for variant in variants:
        for pos in range(variant['pos'], variant['end'] + 1):
            variant_positions[pos].append(variant)
    
    # Check each feature for overlaps with variants
    for feature in features:
        overlapping_variants = set()
        
        # Check if any position in the feature overlaps with a variant
        for pos in range(feature['start'], feature['end'] + 1):
            if pos in variant_positions:
                for variant in variant_positions[pos]:
                    overlapping_variants.add(variant['id'])
        
        # Add annotation if there are overlapping variants
        if overlapping_variants:
            feature['attributes']['Overlapping_variants'] = ','.join(sorted(overlapping_variants))
    
    logging.info(f"Annotated features with variant effects")
    return features

def write_gff3(features, output_file, seqid="1", source="gfa_to_gff3"):
    """Write features to a GFF3 file."""
    # Sort features by start position and then by end position
    features.sort(key=lambda f: (f['start'], f['end']))
    
    with open(output_file, 'w') as f:
        # Write GFF3 header
        f.write("##gff-version 3\n")
        f.write(f"##sequence-region {seqid} 1 1000000\n")  # Placeholder for sequence length
        
        # Write features
        for feature in features:
            # Format attributes
            attributes = []
            for key, value in feature['attributes'].items():
                attributes.append(f"{key}={value}")
            
            # Join attributes with semicolons
            attributes_str = ';'.join(attributes)
            
            # Write GFF3 line
            line = [
                seqid,                  # seqid
                source,                 # source
                feature['type'],        # type
                str(feature['start']),  # start
                str(feature['end']),    # end
                '.',                    # score
                feature['strand'],      # strand
                '.',                    # phase
                attributes_str          # attributes
            ]
            
            f.write('\t'.join(line) + '\n')
    
    logging.info(f"Wrote {len(features)} features to GFF3 file: {output_file}")
    return True

def export_fasta(sequence, output_file, seqid="1", description="Reference sequence"):
    """Export sequence to a FASTA file using BioPython."""
    from Bio.SeqRecord import SeqRecord
    
    # Create a SeqRecord
    record = SeqRecord(
        Seq(sequence),
        id=seqid,
        name=seqid,
        description=description
    )
    
    # Write to file
    with open(output_file, 'w') as f:
        SeqIO.write(record, f, "fasta")
    
    logging.info(f"Wrote sequence of length {len(sequence)}bp to FASTA file: {output_file}")
    return True

def main():
    """Main function to run the GFA to GFF3 conversion."""
    parser = argparse.ArgumentParser(description='Generate a GFF3 file based on GFA paths and VCF variants.')
    parser.add_argument('gfa_file', help='Input GFA file')
    parser.add_argument('output_file', help='Output GFF3 file')
    
    # Input options
    parser.add_argument('--vcf', help='Input VCF file with variants')
    parser.add_argument('--ref-path', default='REF', help='Name of the reference path in GFA (default: %(default)s)')
    
    # Feature generation options
    parser.add_argument('--genes', type=int, default=5, help='Number of synthetic genes to generate (default: %(default)s)')
    parser.add_argument('--exons', type=int, default=3, help='Number of exons per gene (default: %(default)s)')
    parser.add_argument('--min-gene-length', type=int, default=200, help='Minimum gene length (default: %(default)s bp)')
    parser.add_argument('--max-gene-length', type=int, default=500, help='Maximum gene length (default: %(default)s bp)')
    parser.add_argument('--include-variants', action='store_true', help='Include variant features in GFF3')
    parser.add_argument('--annotate-effects', action='store_true', help='Annotate features with variant effects')
    
    # Output options
    parser.add_argument('--seq-id', default='1', help='Sequence ID to use in GFF3 (default: %(default)s)')
    parser.add_argument('--source', default='gfa_to_gff3', help='Source field for GFF3 features (default: %(default)s)')
    parser.add_argument('--export-fasta', help='Export reference sequence to this FASTA file')
    
    # Debug and logging options
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output without full debug')
    parser.add_argument('--log-file', help='Write log to this file')
    parser.add_argument('--seed', type=int, help='Random seed for reproducible feature generation')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(debug=args.debug, log_file=args.log_file, verbose=args.verbose)
    
    # Set random seed if provided
    if args.seed:
        random.seed(args.seed)
        logging.info(f"Using random seed: {args.seed}")
    
    try:
        # Parse GFA file
        segments, links, paths = parse_gfa(args.gfa_file)
        
        # Get reference path
        if args.ref_path not in paths:
            available_paths = ", ".join(paths.keys())
            if available_paths:
                logging.error(f"Reference path '{args.ref_path}' not found. Available paths: {available_paths}")
            else:
                logging.error(f"No paths found in GFA file")
            return 1
        
        ref_path = paths[args.ref_path]
        
        # Build reference sequence and segment offsets
        reference_seq, segment_offsets = build_reference_sequence(segments, ref_path)
        
        # Parse variants if provided
        variants = []
        if args.vcf:
            variants = parse_vcf(args.vcf)
        
        # Create synthetic features
        features = create_synthetic_features(
            reference_seq,
            num_genes=args.genes,
            num_exons_per_gene=args.exons,
            min_gene_length=args.min_gene_length,
            max_gene_length=args.max_gene_length
        )
        
        # Add variant features if requested
        if args.include_variants and variants:
            variant_features = create_variant_features(variants, len(reference_seq))
            features.extend(variant_features)
        
        # Annotate variant effects if requested
        if args.annotate_effects and variants:
            features = apply_variant_effects(features, variants)
        
        # Write GFF3 file
        write_gff3(features, args.output_file, seqid=args.seq_id, source=args.source)
        
        # Export FASTA if requested
        if args.export_fasta:
            export_fasta(reference_seq, args.export_fasta, seqid=args.seq_id)
        
        logging.info("GFF3 generation completed successfully")
        return 0
        
    except Exception as e:
        logging.error(f"Error during GFF3 generation: {e}")
        if args.debug:
            import traceback
            logging.error(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())
