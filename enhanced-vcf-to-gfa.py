#!/usr/bin/env python
"""
Convert variants from VCF to a new path in GFA format.
This script parses VCF variants using the hgvs library and BioPython,
then applies them to create an alternate path in GFA.

Author: Claude
Date: 2025-03-11
"""

import argparse
import re
import sys
import os
import logging
from collections import defaultdict
import hgvs.parser
import hgvs.validator
import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.Align import substitution_matrices
import multiprocessing as mp

def setup_logging(debug=False, log_file=None, verbose=False):
    """Configure logging based on debug flag and optional log file.
    
    Args:
        debug (bool): Enable debug mode logging
        log_file (str): Optional file path for logging
        verbose (bool): Enable verbose output without full debug
    """
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
    """Parse a GFA file into segments, links, and paths with optional validation.
    
    Args:
        gfa_file (str): Path to GFA file
        validate (bool): Whether to validate GFA structure
        
    Returns:
        tuple: (segments, links, paths) dictionaries and lists
    """
    segments = {}
    links = []
    paths = {}
    segment_counts = {}  # For validation
    
    start_time = time.time()
    logging.info(f"Parsing GFA file: {gfa_file}")
    
    with open(gfa_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
                
            fields = line.split('\t')
            record_type = fields[0]
            
            if record_type == 'S':  # Segment
                if len(fields) < 3:
                    logging.warning(f"Line {line_num}: Invalid segment record, missing fields")
                    continue
                    
                seg_id = fields[1]
                sequence = fields[2]
                
                # Store segment sequence and track counts for validation
                segments[seg_id] = sequence
                segment_counts[seg_id] = segment_counts.get(seg_id, 0) + 1
                
                # Extract optional tags if present
                tags = {}
                for field in fields[3:]:
                    if ':' in field:
                        tag_parts = field.split(':', 2)
                        if len(tag_parts) >= 3:
                            tag_name, tag_type, tag_value = tag_parts
                            tags[tag_name] = (tag_type, tag_value)
                
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
                
                # Extract optional tags
                tags = {}
                for field in fields[6:]:
                    if ':' in field:
                        tag_parts = field.split(':', 2)
                        if len(tag_parts) >= 3:
                            tag_name, tag_type, tag_value = tag_parts
                            tags[tag_name] = (tag_type, tag_value)
                
                links.append((from_id, from_dir, to_id, to_dir, overlap, tags))
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
                
                # Extract path overlaps if present
                overlaps = fields[3] if len(fields) > 3 else "*"
                
                # Store with overlaps
                paths[path_name] = {'segments': path_segments, 'overlaps': overlaps}
                logging.debug(f"Parsed path: {path_name} with {len(path_segments)} segments")
    
    # Validate GFA structure if requested
    if validate:
        validation_errors = []
        
        # Check for duplicate segments
        for seg_id, count in segment_counts.items():
            if count > 1:
                validation_errors.append(f"Duplicate segment ID: {seg_id} (appears {count} times)")
        
        # Check for links referencing non-existent segments
        for from_id, from_dir, to_id, to_dir, overlap, tags in links:
            if from_id not in segments:
                validation_errors.append(f"Link references non-existent segment: {from_id}")
            if to_id not in segments:
                validation_errors.append(f"Link references non-existent segment: {to_id}")
        
        if validation_errors:
            for error in validation_errors:
                logging.warning(f"GFA validation error: {error}")
    
    elapsed = time.time() - start_time
    logging.info(f"Finished parsing GFA in {elapsed:.2f}s: {len(segments)} segments, {len(links)} links, {len(paths)} paths")
    
    # Simplify paths structure for compatibility with rest of code
    simplified_paths = {}
    for path_name, path_data in paths.items():
        simplified_paths[path_name] = path_data['segments']
    
    return segments, links, simplified_paths

def build_reference_sequence(segments, path_segments):
    """Build a reference sequence from segments and calculate offsets.
    
    Args:
        segments (dict): Dictionary of segment ID to sequence
        path_segments (list): List of (segment_id, orientation) tuples
        
    Returns:
        tuple: (reference_seq, segment_offsets)
    """
    reference_seq = ""
    segment_offsets = {}
    current_offset = 0
    
    for seg_id, orientation in path_segments:
        if seg_id not in segments:
            logging.error(f"Missing segment: {seg_id}")
            continue
            
        segment_seq = segments[seg_id]
        # Handle reverse complement if needed
        if orientation == '-':
            segment_seq = reverse_complement(segment_seq)
            
        segment_offsets[seg_id] = (current_offset, current_offset + len(segment_seq) - 1, orientation)
        reference_seq += segment_seq
        current_offset += len(segment_seq)
    
    logging.info(f"Built reference sequence of length {len(reference_seq)}bp from {len(path_segments)} segments")
    return reference_seq, segment_offsets

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence using BioPython.
    
    Args:
        sequence (str): DNA sequence
        
    Returns:
        str: Reverse complemented sequence
    """
    # Use BioPython for reverse complement
    seq_obj = Seq(sequence)
    return str(seq_obj.reverse_complement())

def parse_vcf(vcf_file, strict_hgvs=False, max_variants=None, chrom_filter=None):
    """Parse VCF file and extract variants using hgvs library.
    
    Args:
        vcf_file (str): Path to the VCF file
        strict_hgvs (bool): If True, fail on HGVS parse errors
        max_variants (int): Maximum number of variants to process
        chrom_filter (str): Only process variants from this chromosome
        
    Returns:
        list: List of variant dictionaries
    """
    variants = []
    parser = hgvs.parser.Parser()
    
    start_time = time.time()
    logging.info(f"Parsing VCF file: {vcf_file}")
    
    vcf_stats = {
        'total': 0,
        'passed': 0,
        'filtered': 0,
        'by_type': defaultdict(int),
        'by_chrom': defaultdict(int),
        'invalid': 0
    }
    
    with open(vcf_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith('#'):
                continue
                
            vcf_stats['total'] += 1
            
            # Stop if we've reached max variants
            if max_variants and vcf_stats['passed'] >= max_variants:
                logging.info(f"Reached maximum of {max_variants} variants, stopping")
                break
                
            fields = line.strip().split('\t')
            if len(fields) < 8:
                logging.warning(f"Line {line_num}: Invalid VCF record, missing fields")
                vcf_stats['invalid'] += 1
                continue
                
            chrom = fields[0]
            
            # Apply chromosome filter if provided
            if chrom_filter and chrom != chrom_filter:
                vcf_stats['filtered'] += 1
                continue
                
            try:
                pos = int(fields[1])
            except ValueError:
                logging.warning(f"Line {line_num}: Invalid position: {fields[1]}")
                vcf_stats['invalid'] += 1
                continue
                
            var_id = fields[2] if fields[2] != '.' else f"var_{line_num}"
            ref = fields[3]
            alt = fields[4]
            
            if not ref or not alt or alt == '.':
                logging.warning(f"Line {line_num}: Missing or invalid REF or ALT")
                vcf_stats['invalid'] += 1
                continue
            
            # Handle multiple ALT alleles
            alt_alleles = alt.split(',')
            if len(alt_alleles) > 1:
                logging.debug(f"Line {line_num}: Processing {len(alt_alleles)} alternate alleles")
            
            # Parse INFO field
            info_dict = {}
            for item in fields[7].split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    info_dict[key] = value
                else:
                    info_dict[item] = True
            
            # Process each ALT allele
            for idx, alt_allele in enumerate(alt_alleles):
                current_id = var_id if len(alt_alleles) == 1 else f"{var_id}_{idx+1}"
                
                # Determine variant type based on length and SVTYPE field
                if 'SVTYPE' in info_dict:
                    variant_type = info_dict['SVTYPE']
                else:
                    if len(ref) == 1 and len(alt_allele) == 1:
                        variant_type = 'SNP'
                    elif len(ref) > len(alt_allele):
                        variant_type = 'DEL'
                    elif len(ref) < len(alt_allele):
                        variant_type = 'INS'
                    elif len(ref) == len(alt_allele) and len(ref) > 1:
                        variant_type = 'MNP'
                    else:
                        variant_type = 'OTHER'
                
                try:
                    end_pos = int(info_dict.get('END', pos + len(ref) - 1))
                except ValueError:
                    logging.warning(f"Line {line_num}: Invalid END position: {info_dict.get('END')}")
                    end_pos = pos + len(ref) - 1
                    
                hgvs_notation = info_dict.get('HGVS', '')
                hgvs_obj = None
                
                # Parse the HGVS notation using the hgvs library
                if hgvs_notation:
                    try:
                        hgvs_obj = parser.parse_hgvs_variant(hgvs_notation)
                        logging.debug(f"Parsed HGVS: {hgvs_notation} -> {hgvs_obj}")
                    except Exception as e:
                        error_msg = f"Could not parse HGVS notation for {current_id}: {e}"
                        if strict_hgvs:
                            logging.error(error_msg)
                            raise ValueError(error_msg)
                        else:
                            logging.warning(error_msg)
                
                variant = {
                    'id': current_id,
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt_allele,
                    'type': variant_type,
                    'end': end_pos,
                    'hgvs': hgvs_notation,
                    'hgvs_obj': hgvs_obj,
                    'info': info_dict,
                    'line_num': line_num
                }
                
                variants.append(variant)
                vcf_stats['passed'] += 1
                vcf_stats['by_type'][variant_type] += 1
                vcf_stats['by_chrom'][chrom] += 1
                
                logging.debug(f"Parsed variant: {current_id} {variant_type} at position {pos}")
    
    # Sort variants by position (descending) to avoid position shifts
    variants = sorted(variants, key=lambda v: v['pos'], reverse=True)
    
    elapsed = time.time() - start_time
    logging.info(f"Finished parsing VCF in {elapsed:.2f}s: {len(variants)} variants")
    
    # Log statistics
    logging.info(f"VCF Statistics:")
    logging.info(f"  Total records: {vcf_stats['total']}")
    logging.info(f"  Passed variants: {vcf_stats['passed']}")
    logging.info(f"  Filtered variants: {vcf_stats['filtered']}")
    logging.info(f"  Invalid records: {vcf_stats['invalid']}")
    logging.info(f"  Variant types:")
    for var_type, count in sorted(vcf_stats['by_type'].items()):
        logging.info(f"    {var_type}: {count}")
    
    return variants

def find_segment_for_position(position, segment_offsets, fuzzy_search=False, search_window=10):
    """Find which segment contains a given position and calculate relative position.
    
    Args:
        position (int): Genome position (1-based)
        segment_offsets (dict): Dictionary of segment offsets
        fuzzy_search (bool): If True, allow searching in neighboring segments
        search_window (int): Position window to search if fuzzy_search is True
        
    Returns:
        tuple: (segment_id, relative_position, orientation)
    """
    # First, try exact match
    for seg_id, (start, end, orientation) in segment_offsets.items():
        if start <= position - 1 <= end:  # Convert to 0-based for internal calculations
            rel_pos = position - 1 - start
            # For reverse orientation, we need to recalculate relative position
            if orientation == '-':
                # Length of the segment
                seg_len = end - start + 1
                # Calculate position from end of segment
                rel_pos = seg_len - rel_pos - 1
            return seg_id, rel_pos, orientation
    
    # If fuzzy search is enabled and no exact match found
    if fuzzy_search:
        logging.debug(f"Position {position} not found exactly, trying fuzzy search")
        
        # Find segments close to position
        candidates = []
        pos_0based = position - 1
        
        for seg_id, (start, end, orientation) in segment_offsets.items():
            if abs(start - pos_0based) <= search_window or abs(end - pos_0based) <= search_window:
                # Calculate distance to segment boundary
                if pos_0based < start:
                    dist = start - pos_0based
                    # This would be position 0 in the segment
                    rel_pos = 0
                elif pos_0based > end:
                    dist = pos_0based - end
                    # This would be the last position in the segment
                    rel_pos = end - start
                else:
                    # Should not happen as we already checked for exact matches
                    dist = 0
                    rel_pos = pos_0based - start
                
                # For reverse orientation, recalculate relative position
                if orientation == '-':
                    seg_len = end - start + 1
                    rel_pos = seg_len - rel_pos - 1
                
                candidates.append((seg_id, rel_pos, orientation, dist))
        
        # Sort by distance (ascending)
        candidates.sort(key=lambda x: x[3])
        
        if candidates:
            best_match = candidates[0]
            logging.debug(f"Found closest segment {best_match[0]} at distance {best_match[3]}")
            return best_match[0], best_match[1], best_match[2]
    
    return None, None, None

def apply_variants_to_segments(segments, variants, segment_offsets, path_segments, 
                             allow_mismatches=False, match_threshold=0.8, parallel=False):
    """Apply variants to segments and create modified segments.
    
    Args:
        segments (dict): Dictionary of segment ID to sequence
        variants (list): List of variant dictionaries
        segment_offsets (dict): Dictionary of segment ID to (start, end, orientation) tuple
        path_segments (list): List of (segment_id, orientation) tuples
        allow_mismatches (bool): Whether to allow reference mismatches
        match_threshold (float): Minimum sequence similarity for fuzzy matching
        parallel (bool): Whether to use parallel processing
        
    Returns:
        tuple: (modified_segments, new_path_segments, variant_application_log)
    """
    start_time = time.time()
    logging.info("Applying variants to segments...")
    
    # Group variants by segment
    variants_by_segment = defaultdict(list)
    unmapped_variants = []
    
    # Try to map variants to segments
    for variant in variants:
        # First try exact mapping
        seg_id, rel_pos, orientation = find_segment_for_position(variant['pos'], segment_offsets)
        
        # If unsuccessful, try fuzzy mapping if allowed
        if not seg_id and allow_mismatches:
            logging.debug(f"Trying fuzzy segment mapping for variant {variant['id']}")
            seg_id, rel_pos, orientation = find_segment_for_position(
                variant['pos'], segment_offsets, fuzzy_search=True, search_window=10
            )
            if seg_id:
                logging.warning(f"Variant {variant['id']} mapped to segment {seg_id} using fuzzy positioning")
        
        if seg_id:
            variants_by_segment[seg_id].append((variant, rel_pos, orientation))
            logging.debug(f"Mapped variant {variant['id']} to segment {seg_id} at relative position {rel_pos}")
        else:
            unmapped_variants.append(variant)
            logging.warning(f"Could not map variant {variant['id']} (pos {variant['pos']}) to any segment")
    
    if unmapped_variants:
        logging.warning(f"Failed to map {len(unmapped_variants)} variants to segments")
    
    # Create modified segments
    modified_segments = {}
    variant_application_log = []
    
    # Function to process a single segment
    def process_segment(seg_id, var_list):
        # Sort by relative position (descending)
        var_list.sort(key=lambda x: x[1], reverse=True)
        
        # Get the original sequence
        segment_seq = segments[seg_id]
        segment_seq_original = segment_seq  # Keep original for comparison
        
        logging.info(f"Modifying segment {seg_id} with {len(var_list)} variants")
        
        segment_log = []
        
        # Apply variants
        for variant, rel_pos, orientation in var_list:
            expected_ref = segment_seq[rel_pos:rel_pos + len(variant['ref'])]
            
            # Check for reference match
            if expected_ref != variant['ref']:
                # Try fuzzy matching for mismatches
                alignment_score = 0
                if len(expected_ref) > 0 and len(variant['ref']) > 0:
                    # Use BioPython for sequence alignment to find best match
                    alignments = pairwise2.align.localms(
                        expected_ref, variant['ref'], 
                        2, -1, -2, -0.5,  # Match, mismatch, gap penalties
                        one_alignment_only=True
                    )
                    if alignments:
                        alignment_score = alignments[0].score / (2 * max(len(expected_ref), len(variant['ref'])))
                
                if alignment_score < match_threshold:
                    msg = f"Reference mismatch for {variant['id']} in segment {seg_id} at position {rel_pos}."
                    msg += f"\n  Expected: {variant['ref']}, Found: {expected_ref}"
                    msg += f"\n  Similarity score: {alignment_score:.2f}"
                    
                    if not allow_mismatches:
                        logging.error(msg)
                        raise ValueError(f"Reference mismatch for {variant['id']}. Use --allow-mismatches to force application.")
                    else:
                        logging.warning(msg)
                else:
                    logging.info(f"Fuzzy match for {variant['id']}: Score {alignment_score:.2f}")
            
            # Apply the variant based on type
            old_seq = segment_seq
            
            if variant['type'] == 'SNP':
                logging.debug(f"Applying SNP {variant['id']}: {variant['ref']} -> {variant['alt']} at position {rel_pos}")
                segment_seq = segment_seq[:rel_pos] + variant['alt'] + segment_seq[rel_pos + len(variant['ref']):]
            elif variant['type'] in ['INV', 'MNP']:
                logging.debug(f"Applying {variant['type']} {variant['id']}: {variant['ref']} -> {variant['alt']} at position {rel_pos}")
                segment_seq = segment_seq[:rel_pos] + variant['alt'] + segment_seq[rel_pos + len(variant['ref']):]
            elif variant['type'] == 'DEL':
                logging.debug(f"Applying DEL {variant['id']}: Removing {len(variant['ref'])} bp at position {rel_pos}")
                segment_seq = segment_seq[:rel_pos] + (variant['alt'] if variant['alt'] != '.' else '') + segment_seq[rel_pos + len(variant['ref']):]
            elif variant['type'] == 'INS':
                logging.debug(f"Applying INS {variant['id']}: Inserting {len(variant['alt'])-1} bp at position {rel_pos}")
                # Handle INS where ref is the base before insertion
                if len(variant['ref']) == 1:
                    segment_seq = segment_seq[:rel_pos + 1] + variant['alt'][1:] + segment_seq[rel_pos + 1:]
                else:
                    segment_seq = segment_seq[:rel_pos] + variant['alt'] + segment_seq[rel_pos + len(variant['ref']):]
            else:
                logging.debug(f"Applying generic variant {variant['id']} ({variant['type']}): {variant['ref']} -> {variant['alt']} at position {rel_pos}")
                segment_seq = segment_seq[:rel_pos] + variant['alt'] + segment_seq[rel_pos + len(variant['ref']):]
            
            # Log the change
            segment_log.append({
                'variant_id': variant['id'],
                'segment_id': seg_id,
                'position': variant['pos'],
                'rel_position': rel_pos,
                'orientation': orientation,
                'type': variant['type'],
                'ref': variant['ref'],
                'alt': variant['alt'],
                'seq_before': old_seq[max(0, rel_pos-5):min(len(old_seq), rel_pos+len(variant['ref'])+5)],
                'seq_after': segment_seq[max(0, rel_pos-5):min(len(segment_seq), rel_pos+len(variant['alt'])+5)]
            })
        
        # Create a new segment ID
        var_ids = '_'.join(v[0]['id'] for v in var_list)
        new_seg_id = f"{seg_id}_{var_ids}"
        
        # Store the modified segment
        modified_segment = {
            'new_id': new_seg_id,
            'sequence': segment_seq,
            'original_sequence': segment_seq_original,
            'variants': [v[0] for v in var_list],
            'length_diff': len(segment_seq) - len(segment_seq_original)
        }
        
        logging.info(f"Created modified segment {new_seg_id}, length diff: {modified_segment['length_diff']} bp")
        
        return seg_id, modified_segment, segment_log
    
    # Process segments (either in parallel or serially)
    if parallel and len(variants_by_segment) > 1:
        logging.info(f"Processing {len(variants_by_segment)} segments in parallel")
        with mp.Pool(processes=min(mp.cpu_count(), len(variants_by_segment))) as pool:
            results = pool.starmap(process_segment, [(seg_id, var_list) for seg_id, var_list in variants_by_segment.items()])
            
            for seg_id, mod_seg, seg_log in results:
                modified_segments[seg_id] = mod_seg
                variant_application_log.extend(seg_log)
    else:
        for seg_id, var_list in variants_by_segment.items():
            seg_id, mod_seg, seg_log = process_segment(seg_id, var_list)
            modified_segments[seg_id] = mod_seg
            variant_application_log.extend(seg_log)
    
    # Create a new path with modified segments
    new_path_segments = []
    for seg_id, orientation in path_segments:
        if seg_id in modified_segments:
            new_path_segments.append((modified_segments[seg_id]['new_id'], orientation))
        else:
            new_path_segments.append((seg_id, orientation))
    
    elapsed = time.time() - start_time
    logging.info(f"Created new path with {len(modified_segments)} modified segments in {elapsed:.2f}s")
    
    return modified_segments, new_path_segments, variant_application_log

def generate_new_gfa(original_segments, modified_segments, links, paths, new_path_segments, 
                output_file, alt_path_name="ALT", keep_original=True, add_metadata=True, 
                add_variant_annotations=True, variants=None):
    """Generate a new GFA file with the modified segments and paths.
    
    Args:
        original_segments (dict): Dictionary of original segment ID to sequence
        modified_segments (dict): Dictionary of modified segment information
        links (list): List of link tuples
        paths (dict): Dictionary of path name to segments
        new_path_segments (list): List of (segment_id, orientation) tuples for new path
        output_file (str): Path to output GFA file
        alt_path_name (str): Name for the alternate path
        keep_original (bool): Whether to keep original segments and paths
        add_metadata (bool): Whether to add metadata tags to GFA
    """
    start_time = time.time()
    logging.info(f"Writing new GFA to {output_file}")
    
    with open(output_file, 'w') as f:
        # Write header with metadata
        if add_metadata:
            timestamp = time.strftime("%Y-%m-%dT%H:%M:%S")
            f.write(f"H\tVN:Z:1.0\tTS:Z:{timestamp}\tPG:Z:vcf_to_gfa.py\n")
        else:
            f.write("H\tVN:Z:1.0\n")
            
        # Add variant annotations if requested and variants are provided
        if add_variant_annotations and variants:
            f.write("# Variant Annotations\n")
            f.write("# Format: #VAR <variant_id> <type> <position> <reference> <alternate> [additional_info]\n")
            for var in variants:
                # Format variant information
                var_id = var['id']
                var_type = var['type']
                var_pos = var['pos']
                var_ref = var['ref']
                var_alt = var['alt']
                
                # Add basic variant line
                f.write(f"# VAR\t{var_id}\t{var_type}\t{var_pos}\t{var_ref}\t{var_alt}")
                
                # Add additional information if available
                if 'info' in var and var['info']:
                    info_str = "\t" + ";".join(f"{k}={v}" for k, v in var['info'].items() 
                                              if k not in ['SVTYPE', 'END', 'HGVS'] and isinstance(v, str))
                    if info_str != "\t":
                        f.write(info_str)
                
                f.write("\n")
            
            f.write("# End of Variant Annotations\n")
        
        # Write original segments if requested
        segments_written = 0
        if keep_original:
            for seg_id, sequence in original_segments.items():
                f.write(f"S\t{seg_id}\t{sequence}\tLN:i:{len(sequence)}\n")
                segments_written += 1
                logging.debug(f"Wrote original segment {seg_id}")
        
        # Write modified segments
        for seg_id, mod_seg in modified_segments.items():
            new_id = mod_seg['new_id']
            sequence = mod_seg['sequence']
            
            # Add segment metadata
            tags = [f"LN:i:{len(sequence)}"]
            
            if add_metadata:
                # Add variant metadata
                variant_ids = ';'.join(v['id'] for v in mod_seg['variants'])
                tags.append(f"VA:Z:{variant_ids}")
                
                # Add length difference
                if mod_seg['length_diff'] != 0:
                    tags.append(f"LD:i:{mod_seg['length_diff']}")
                
                # Add original segment reference
                tags.append(f"OR:Z:{seg_id}")
                
                # Add detailed variant information for preservation
                for i, variant in enumerate(mod_seg['variants']):
                    # Add position, ref and alt for each variant
                    var_prefix = f"V{i+1}"
                    tags.append(f"{var_prefix}ID:Z:{variant['id']}")
                    tags.append(f"{var_prefix}POS:i:{variant['pos']}")
                    tags.append(f"{var_prefix}TYPE:Z:{variant['type']}")
                    tags.append(f"{var_prefix}REF:Z:{variant['ref']}")
                    tags.append(f"{var_prefix}ALT:Z:{variant['alt']}")
                    
                    # Add additional metadata if available
                    if 'info' in variant and variant['info']:
                        for key, value in variant['info'].items():
                            if key not in ['SVTYPE', 'END', 'HGVS'] and isinstance(value, str):
                                # Clean and truncate value if needed
                                clean_value = value.replace(' ', '_').replace('\t', '_')
                                if len(clean_value) > 50:  # Limit tag length
                                    clean_value = clean_value[:47] + "..."
                                tags.append(f"{var_prefix}{key}:Z:{clean_value}")
            
            # Write the segment with tags
            f.write(f"S\t{new_id}\t{sequence}\t{' '.join(tags)}\n")
            segments_written += 1
            logging.debug(f"Wrote modified segment {new_id}")
        
        # Write original links if requested
        links_written = 0
        if keep_original:
            for link_data in links:
                # Handle both old and new link formats
                if len(link_data) == 5:
                    from_id, from_dir, to_id, to_dir, overlap = link_data
                    tags = {}
                else:
                    from_id, from_dir, to_id, to_dir, overlap, tags = link_data
                
                # Convert tags to string format
                tag_str = ""
                if tags:
                    tag_str = " " + " ".join(f"{name}:{type_val[0]}:{type_val[1]}" for name, type_val in tags.items())
                
                f.write(f"L\t{from_id}\t{from_dir}\t{to_id}\t{to_dir}\t{overlap}{tag_str}\n")
                links_written += 1
                logging.debug(f"Wrote original link {from_id}{from_dir} -> {to_id}{to_dir}")
        
        # Write new links for modified segments
        ref_path_name = args.ref_path if 'args' in locals() else next(iter(paths.keys()))
        if "REF" in paths and ref_path_name != "REF":
            ref_path_name = "REF"
            
        old_path = paths[ref_path_name]
        
        new_links_added = 0
        for i in range(len(old_path) - 1):
            from_seg, from_dir = old_path[i]
            to_seg, to_dir = old_path[i + 1]
            
            # If either segment is modified, create a new link
            if from_seg in modified_segments or to_seg in modified_segments:
                new_from = modified_segments[from_seg]['new_id'] if from_seg in modified_segments else from_seg
                new_to = modified_segments[to_seg]['new_id'] if to_seg in modified_segments else to_seg
                
                # Add metadata to link
                tags = []
                if add_metadata:
                    tags.append("SV:Z:variant")
                
                tag_str = " " + " ".join(tags) if tags else ""
                
                # Avoid duplicating links
                f.write(f"L\t{new_from}\t{from_dir}\t{new_to}\t{to_dir}\t0M{tag_str}\n")
                new_links_added += 1
                links_written += 1
                logging.debug(f"Wrote new link {new_from}{from_dir} -> {new_to}{to_dir}")
        
        logging.info(f"Added {new_links_added} new links")
        
        # Write original paths if requested
        paths_written = 0
        if keep_original:
            for path_name, path_segs in paths.items():
                segments_str = ','.join(f"{seg}{dir}" for seg, dir in path_segs)
                f.write(f"P\t{path_name}\t{segments_str}\t*\n")
                paths_written += 1
                logging.debug(f"Wrote original path {path_name}")
        
        # Write new path with variants
        new_path_str = ','.join(f"{seg}{dir}" for seg, dir in new_path_segments)
        
        # Add path metadata tags
        tags = []
        if add_metadata:
            tags.append(f"VN:Z:variant_path")
            tags.append(f"SR:Z:{ref_path_name}")
            
        tag_str = "\t" + "\t".join(tags) if tags else ""
        
        f.write(f"P\t{alt_path_name}\t{new_path_str}\t*{tag_str}\n")
        paths_written += 1
        logging.info(f"Wrote new path {alt_path_name} with {len(new_path_segments)} segments")
    
    elapsed = time.time() - start_time
    logging.info(f"GFA file written successfully to {output_file} in {elapsed:.2f}s")
    logging.info(f"  Segments: {segments_written}")
    logging.info(f"  Links: {links_written}")
    logging.info(f"  Paths: {paths_written}")
    
    # Validate file was written
    if not os.path.exists(output_file):
        logging.error(f"Failed to write output file: {output_file}")
        return False
    
    return True

def write_variant_report(variant_log, report_file, include_sequences=True, format='markdown'):
    """Write a detailed report of variant applications.
    
    Args:
        variant_log (list): List of variant application log entries
        report_file (str): Output report file path
        include_sequences (bool): Whether to include sequence context
        format (str): Output format ('markdown', 'tsv', or 'json')
    """
    start_time = time.time()
    
    if format == 'markdown':
        with open(report_file, 'w') as f:
            f.write("# Variant Application Report\n\n")
            
            # Group variants by segment
            variants_by_segment = defaultdict(list)
            for entry in variant_log:
                variants_by_segment[entry['segment_id']].append(entry)
            
            # Write summary table
            f.write("## Summary\n\n")
            f.write("| Segment | Variants Applied |\n")
            f.write("|---------|------------------|\n")
            for seg_id, variants in variants_by_segment.items():
                f.write(f"| {seg_id} | {len(variants)} |\n")
            f.write("\n")
            
            # Write detailed variant information
            f.write("## Variant Details\n\n")
            for entry in variant_log:
                f.write(f"### Variant: {entry['variant_id']} ({entry['type']})\n")
                f.write(f"Segment: {entry['segment_id']}\n")
                f.write(f"Position: {entry['position']} (relative: {entry['rel_position']})\n")
                f.write(f"Orientation: {entry['orientation']}\n")
                f.write(f"Reference: {entry['ref']}\n")
                f.write(f"Alternate: {entry['alt']}\n")
                
                if include_sequences:
                    f.write("Sequence context before:\n")
                    f.write(f"```\n{entry['seq_before']}\n```\n")
                    f.write("Sequence context after:\n")
                    f.write(f"```\n{entry['seq_after']}\n```\n")
                f.write("\n")
    
    elif format == 'tsv':
        with open(report_file, 'w') as f:
            # Write header
            header = ["variant_id", "segment_id", "position", "rel_position", "orientation", 
                     "type", "ref", "alt"]
            if include_sequences:
                header.extend(["seq_before", "seq_after"])
            f.write('\t'.join(header) + '\n')
            
            # Write data
            for entry in variant_log:
                row = [
                    entry['variant_id'],
                    entry['segment_id'],
                    str(entry['position']),
                    str(entry['rel_position']),
                    entry['orientation'],
                    entry['type'],
                    entry['ref'],
                    entry['alt']
                ]
                if include_sequences:
                    row.extend([entry['seq_before'], entry['seq_after']])
                f.write('\t'.join(row) + '\n')
    
    elif format == 'json':
        import json
        with open(report_file, 'w') as f:
            if not include_sequences:
                # Remove sequence fields if not requested
                for entry in variant_log:
                    entry.pop('seq_before', None)
                    entry.pop('seq_after', None)
            json.dump(variant_log, f, indent=2)
    
    else:
        logging.error(f"Unknown report format: {format}")
        return False
    
    elapsed = time.time() - start_time
    logging.info(f"Variant report written to {report_file} in {elapsed:.2f}s")
    return True

def export_sequences(segments, original_path, modified_path, output_prefix, format='fasta', annotate_variants=True, variant_log=None):
    """Export reference and modified sequences to FASTA/GenBank files using BioPython.
    
    Args:
        segments (dict): Dictionary of segment ID to sequence
        original_path (list): List of (segment_id, orientation) tuples for reference path
        modified_path (list): List of (segment_id, orientation) tuples for alternate path
        output_prefix (str): Prefix for output files
        format (str): Output format ('fasta' or 'genbank')
    
    Returns:
        tuple: (ref_file, alt_file) paths
    """
    # Build reference sequence
    ref_seq = ""
    for seg_id, orientation in original_path:
        if seg_id in segments:
            seq = segments[seg_id]
            if orientation == '-':
                seq = reverse_complement(seq)
            ref_seq += seq
    
    # Build alternate sequence
    alt_seq = ""
    for seg_id, orientation in modified_path:
        if seg_id in segments:
            seq = segments[seg_id]
            if orientation == '-':
                seq = reverse_complement(seq)
            alt_seq += seq
    
    # Create BioPython SeqRecord objects
    ref_record = SeqRecord(
        Seq(ref_seq),
        id="REF",
        name="reference",
        description="Reference sequence"
    )
    
    # Create description with variant information if available
    alt_description = "Alternate sequence with variants"
    if annotate_variants and variant_log:
        # Group variants by ID
        var_groups = {}
        for entry in variant_log:
            var_id = entry['variant_id']
            if var_id not in var_groups:
                var_groups[var_id] = entry
        
        # Add variant info to description (limited to keep it reasonable)
        if var_groups:
            var_info = ", ".join([f"{v_id}:{v['type']}@{v['position']}" 
                               for v_id, v in list(var_groups.items())[:5]])
            if len(var_groups) > 5:
                var_info += f" and {len(var_groups)-5} more"
            alt_description = f"Alternate sequence with variants: {var_info}"
    
    alt_record = SeqRecord(
        Seq(alt_seq),
        id="ALT",
        name="alternate",
        description=alt_description
    )
    
    # Set output file extensions
    if format.lower() == 'fasta':
        ext = 'fa'
        output_format = 'fasta'
    elif format.lower() in ['genbank', 'gb']:
        ext = 'gb'
        output_format = 'genbank'
    else:
        logging.error(f"Unknown sequence format: {format}")
        ext = 'fa'
        output_format = 'fasta'
    
    # Write sequences to files
    ref_file = f"{output_prefix}_ref.{ext}"
    alt_file = f"{output_prefix}_alt.{ext}"
    
    with open(ref_file, 'w') as f:
        SeqIO.write(ref_record, f, output_format)
    
    with open(alt_file, 'w') as f:
        SeqIO.write(alt_record, f, output_format)
    
    logging.info(f"Reference sequence ({len(ref_seq)} bp) exported to {ref_file}")
    logging.info(f"Modified sequence ({len(alt_seq)} bp) exported to {alt_file}")
    
    return ref_file, alt_file

def compare_sequences(ref_file, alt_file, output_file=None):
    """Compare reference and alternate sequences using BioPython.
    
    Args:
        ref_file (str): Path to reference sequence file
        alt_file (str): Path to alternate sequence file
        output_file (str): Optional path to write alignment results
        
    Returns:
        dict: Comparison statistics
    """
    # Read sequences
    ref_record = SeqIO.read(ref_file, "fasta")
    alt_record = SeqIO.read(alt_file, "fasta")
    
    ref_seq = str(ref_record.seq)
    alt_seq = str(alt_record.seq)
    
    # Calculate basic statistics
    ref_len = len(ref_seq)
    alt_len = len(alt_seq)
    length_diff = alt_len - ref_len
    
    logging.info(f"Comparing sequences:")
    logging.info(f"  Reference length: {ref_len} bp")
    logging.info(f"  Alternate length: {alt_len} bp")
    logging.info(f"  Length difference: {length_diff} bp ({(length_diff/ref_len*100):.2f}%)")
    
    # Count nucleotide composition
    ref_composition = {base: ref_seq.count(base) for base in 'ACGTN'}
    alt_composition = {base: alt_seq.count(base) for base in 'ACGTN'}
    
    # Perform sequence alignment to find differences
    alignments = pairwise2.align.globalms(
        ref_seq[:10000] if len(ref_seq) > 10000 else ref_seq,  # Limit size for performance
        alt_seq[:10000] if len(alt_seq) > 10000 else alt_seq,
        2, -1, -2, -0.5,  # Match, mismatch, gap penalties
        one_alignment_only=True
    )
    
    # Calculate identity and similarity
    if alignments:
        alignment = alignments[0]
        matches = sum(a == b for a, b in zip(alignment.seqA, alignment.seqB) 
                      if a != '-' and b != '-')
        aligned_length = sum(1 for a, b in zip(alignment.seqA, alignment.seqB) 
                           if a != '-' or b != '-')
        identity = matches / aligned_length if aligned_length > 0 else 0
    else:
        identity = 0
    
    logging.info(f"  Sequence identity: {identity:.2%}")
    
    # Write output if requested
    if output_file:
        with open(output_file, 'w') as f:
            f.write("# Sequence Comparison Report\n\n")
            f.write(f"Reference length: {ref_len} bp\n")
            f.write(f"Alternate length: {alt_len} bp\n")
            f.write(f"Length difference: {length_diff} bp ({(length_diff/ref_len*100):.2f}%)\n")
            f.write(f"Sequence identity: {identity:.2%}\n\n")
            
            # Write composition table
            f.write("## Nucleotide Composition\n\n")
            f.write("| Base | Reference | Alternate | Difference |\n")
            f.write("|------|-----------|-----------|------------|\n")
            for base in 'ACGTN':
                ref_count = ref_composition[base]
                alt_count = alt_composition[base]
                diff = alt_count - ref_count
                diff_pct = (diff / ref_count * 100) if ref_count > 0 else float('inf')
                f.write(f"| {base} | {ref_count} | {alt_count} | {diff:+} ({diff_pct:+.2f}%) |\n")
    
    return {
        'ref_length': ref_len,
        'alt_length': alt_len,
        'length_diff': length_diff,
        'identity': identity,
        'ref_composition': ref_composition,
        'alt_composition': alt_composition
    }

def main():
    """Main function to run the VCF to GFA conversion."""
    parser = argparse.ArgumentParser(description='Convert VCF variants to a new path in GFA.')
    parser.add_argument('gfa_file', help='Input GFA file')
    parser.add_argument('vcf_file', help='Input VCF file with variants')
    parser.add_argument('output_file', help='Output GFA file')
    
    # Input/output options
    parser.add_argument('--ref-path', default='REF', help='Name of the reference path in GFA (default: %(default)s)')
    parser.add_argument('--alt-path', default='ALT', help='Name for the alternate path in output GFA (default: %(default)s)')
    parser.add_argument('--report', help='Write a detailed variant application report to this file')
    parser.add_argument('--report-format', choices=['markdown', 'tsv', 'json'], default='markdown',
                        help='Format for the variant report (default: %(default)s)')
    parser.add_argument('--export-sequences', help='Export reference and alternate sequences to FASTA files with this prefix')
    parser.add_argument('--seq-format', choices=['fasta', 'genbank'], default='fasta',
                        help='Format for exported sequences (default: %(default)s)')
    parser.add_argument('--compare', action='store_true', help='Compare reference and alternate sequences')
    
    # Variant processing options
    parser.add_argument('--allow-mismatches', action='store_true', help='Allow reference mismatches when applying variants')
    parser.add_argument('--match-threshold', type=float, default=0.8, 
                        help='Minimum sequence similarity for fuzzy matching (default: %(default)s)')
    parser.add_argument('--strict-hgvs', action='store_true', help='Fail on HGVS parse errors')
    parser.add_argument('--max-variants', type=int, help='Maximum number of variants to process')
    parser.add_argument('--chrom', help='Only process variants from this chromosome')
    
    # GFA output options
    parser.add_argument('--no-original', action='store_true', help='Do not include original segments and paths in output')
    parser.add_argument('--no-metadata', action='store_true', help='Do not add metadata tags to GFA')
    parser.add_argument('--no-variant-annotations', action='store_true', help='Do not add variant annotations to GFA as comments')
    
    # Performance options
    parser.add_argument('--parallel', action='store_true', help='Use parallel processing for segment modification')
    
    # Debug and logging options
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output without full debug')
    parser.add_argument('--log-file', help='Write log to this file')
    parser.add_argument('--validate', action='store_true', help='Validate GFA structure during parsing')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(debug=args.debug, log_file=args.log_file, verbose=args.verbose)
    
    try:
        # Parse input files
        segments, links, paths = parse_gfa(args.gfa_file, validate=args.validate)
        
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
        
        # Parse variants
        variants = parse_vcf(
            args.vcf_file, 
            strict_hgvs=args.strict_hgvs,
            max_variants=args.max_variants,
            chrom_filter=args.chrom
        )
        
        if not variants:
            logging.warning("No variants found or passed filters")
            return 0
        
        # Apply variants to segments
        modified_segments, new_path_segments, variant_log = apply_variants_to_segments(
            segments, 
            variants, 
            segment_offsets, 
            ref_path, 
            allow_mismatches=args.allow_mismatches,
            match_threshold=args.match_threshold,
            parallel=args.parallel
        )
        
        # Add new segments to segments dictionary for output
        for seg_id, mod_seg in modified_segments.items():
            segments[mod_seg['new_id']] = mod_seg['sequence']
        
        # Generate new GFA
        generate_new_gfa(
            segments, 
            modified_segments, 
            links, 
            paths, 
            new_path_segments,
            args.output_file,
            alt_path_name=args.alt_path,
            keep_original=not args.no_original,
            add_metadata=not args.no_metadata,
            add_variant_annotations=not args.no_variant_annotations,
            variants=variants
        )
        
        # Write variant report if requested
        if args.report:
            write_variant_report(
                variant_log, 
                args.report, 
                include_sequences=True,
                format=args.report_format
            )
        
        # Export sequences if requested
        if args.export_sequences:
            ref_file, alt_file = export_sequences(
                segments,
                ref_path,
                new_path_segments,
                args.export_sequences,
                format=args.seq_format,
                annotate_variants=True,
                variant_log=variant_log
            )
            
            # Compare sequences if requested
            if args.compare:
                compare_file = f"{args.export_sequences}_comparison.md"
                compare_sequences(ref_file, alt_file, compare_file)
        
        logging.info("Conversion completed successfully")
        return 0
        
    except Exception as e:
        logging.error(f"Error during conversion: {e}")
        if args.debug:
            import traceback
            logging.error(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())
