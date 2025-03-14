#!/usr/bin/env python3
"""
Convert variants from VCF to new paths in GFA format.

This script parses VCF variants and applies them to create alternate paths in GFA.
It supports multi-sample VCFs and can create paths for phased/unphased haplotypes.

Usage:
  python -m hapli.tools.vcf_to_gfa [options] <gfa_file> <vcf_file> <output_file>
"""

import argparse
import logging
import os
import sys
import time
from collections import defaultdict
import multiprocessing as mp

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2

# Required imports
try:
    import cyvcf2
    CYVCF2_AVAILABLE = True
except ImportError:
    CYVCF2_AVAILABLE = False
    logging.warning("cyvcf2 not available. Please install with: pip install cyvcf2")

# Optional import for HGVS support
try:
    import hgvs.parser
    import hgvs.validator
    import hgvs.exceptions
    HGVS_AVAILABLE = True
except ImportError:
    HGVS_AVAILABLE = False
    logging.warning("hgvs package not available. HGVS validation will be skipped.")

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
    """Parse VCF file and extract variants using cyvcf2.
    
    Args:
        vcf_file (str): Path to the VCF file
        strict_hgvs (bool): If True, fail on HGVS parse errors
        max_variants (int): Maximum number of variants to process
        chrom_filter (str): Only process variants from this chromosome
        
    Returns:
        tuple: (variants, samples, format_fields)
           - variants: List of variant dictionaries
           - samples: List of sample names
           - format_fields: List of format fields
    """
    variants = []
    samples = []
    format_fields = []
    
    # Check if cyvcf2 is available
    if not CYVCF2_AVAILABLE:
        logging.error("cyvcf2 is required for VCF parsing. Please install with: pip install cyvcf2")
        raise ImportError("cyvcf2 is required for VCF parsing")
    
    # Initialize HGVS parser if available
    hgvs_parser = None
    if HGVS_AVAILABLE:
        try:
            hgvs_parser = hgvs.parser.Parser()
            logging.debug("HGVS parser initialized successfully")
        except Exception as e:
            logging.warning(f"Failed to initialize HGVS parser: {e}")
            hgvs_parser = None
    
    start_time = time.time()
    logging.info(f"Parsing VCF file: {vcf_file}")
    
    vcf_stats = {
        'total': 0,
        'passed': 0,
        'filtered': 0,
        'by_type': defaultdict(int),
        'by_chrom': defaultdict(int),
        'invalid': 0,
        'phased': 0,
        'unphased': 0,
        'multiallelic': 0
    }
    
    try:
        # Open VCF file with cyvcf2
        vcf = cyvcf2.VCF(vcf_file)
        
        # Get samples
        samples = vcf.samples
        if samples:
            logging.info(f"Found {len(samples)} samples in VCF: {', '.join(samples)}")
        
        # Get format fields
        format_fields = vcf.FORMAT
        
        # Process variants
        for variant in vcf:
            vcf_stats['total'] += 1
            
            # Stop if we've reached max variants
            if max_variants and vcf_stats['passed'] >= max_variants:
                logging.info(f"Reached maximum of {max_variants} variants, stopping")
                break
            
            chrom = variant.CHROM
            
            # Apply chromosome filter if provided
            if chrom_filter and chrom != chrom_filter:
                vcf_stats['filtered'] += 1
                continue
            
            pos = variant.POS
            var_id = variant.ID if variant.ID else f"var_{vcf_stats['total']}"
            ref = variant.REF
            alts = variant.ALT  # List of alternate alleles
            
            # Check for valid REF/ALT
            if not ref or not alts or (len(alts) == 1 and alts[0] == '.'):
                logging.warning(f"Variant {var_id}: Missing or invalid REF or ALT")
                vcf_stats['invalid'] += 1
                continue
            
            # Handle multiple ALT alleles
            if len(alts) > 1:
                logging.debug(f"Variant {var_id}: Processing {len(alts)} alternate alleles")
                vcf_stats['multiallelic'] += 1
            
            # Get INFO fields
            info_dict = {}
            for field in variant.INFO:
                if field:
                    info_dict[field] = variant.INFO.get(field)
            
            # Check if any genotypes are phased
            is_phased = any(variant.gt_phases)
            if is_phased:
                vcf_stats['phased'] += 1
            else:
                vcf_stats['unphased'] += 1
            
            # Process each ALT allele
            for idx, alt_allele in enumerate(alts):
                current_id = var_id if len(alts) == 1 else f"{var_id}_{idx+1}"
                
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
                
                # Get END position
                end_pos = variant.INFO.get('END', pos + len(ref) - 1)
                
                # Handle HGVS notation if present
                hgvs_notation = variant.INFO.get('HGVS', '')
                hgvs_obj = None
                
                # Parse HGVS notation if available and parser is initialized
                if hgvs_notation and hgvs_parser and HGVS_AVAILABLE:
                    try:
                        # Add proper prefix for HGVS parsing if needed
                        if not hgvs_notation.startswith('chr'):
                            prefixed_hgvs = f"chr{chrom}:{hgvs_notation}"
                        else:
                            prefixed_hgvs = hgvs_notation
                            
                        hgvs_obj = hgvs_parser.parse_hgvs_variant(prefixed_hgvs)
                        logging.debug(f"Parsed HGVS: {hgvs_notation} -> {hgvs_obj}")
                    except hgvs.exceptions.HGVSParseError as e:
                        error_msg = f"Could not parse HGVS notation for {current_id}: {hgvs_notation}: {e}"
                        if strict_hgvs:
                            logging.error(error_msg)
                            raise ValueError(error_msg)
                        else:
                            logging.warning(error_msg)
                
                # Process genotype data for this allele
                genotypes = []
                for i, sample_name in enumerate(samples):
                    # Get genotype for this sample
                    gt_type = variant.gt_types[i]  # 0=HOM_REF, 1=HET, 2=HOM_ALT, 3=UNKNOWN
                    gt_bases = variant.gt_bases[i]  # e.g., "A/C"
                    
                    # Parse allele indices
                    if '|' in gt_bases:  # Phased
                        phase_type = 'phased'
                        allele_indices = []
                        for allele in gt_bases.split('|'):
                            if allele == ref:
                                allele_indices.append(0)
                            elif allele == alt_allele:
                                allele_indices.append(idx + 1)
                            else:
                                # Check if it's another alt allele
                                try:
                                    alt_idx = alts.index(allele) + 1
                                    allele_indices.append(alt_idx)
                                except ValueError:
                                    allele_indices.append(None)
                    elif '/' in gt_bases:  # Unphased
                        phase_type = 'unphased'
                        allele_indices = []
                        for allele in gt_bases.split('/'):
                            if allele == ref:
                                allele_indices.append(0)
                            elif allele == alt_allele:
                                allele_indices.append(idx + 1)
                            else:
                                # Check if it's another alt allele
                                try:
                                    alt_idx = alts.index(allele) + 1
                                    allele_indices.append(alt_idx)
                                except ValueError:
                                    allele_indices.append(None)
                    else:
                        # Single allele (haploid)
                        phase_type = 'unknown'
                        if gt_bases == ref:
                            allele_indices = [0]
                        elif gt_bases == alt_allele:
                            allele_indices = [idx + 1]
                        else:
                            # Check if it's another alt allele
                            try:
                                alt_idx = alts.index(gt_bases) + 1
                                allele_indices = [alt_idx]
                            except ValueError:
                                allele_indices = [None]
                    
                    # Check if this specific ALT allele is present in genotype
                    has_alt = (idx + 1) in allele_indices
                    
                    # Get FORMAT fields for this sample
                    format_data = {}
                    for field in format_fields:
                        try:
                            format_data[field] = variant.format(field)[i]
                        except:
                            pass
                    
                    genotypes.append({
                        'sample': sample_name,
                        'gt': '|'.join(map(str, allele_indices)) if phase_type == 'phased' else '/'.join(map(str, allele_indices)),
                        'has_alt': has_alt,
                        'allele_indices': allele_indices,
                        'phase_type': phase_type,
                        'format_data': format_data
                    })
                
                # Create variant dictionary
                variant_dict = {
                    'id': current_id,
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt_allele,
                    'type': variant_type,
                    'end': end_pos,
                    'genotypes': genotypes,
                    'is_phased': is_phased,
                    'allele_index': idx + 1,
                    'info': info_dict,
                    'line_num': vcf_stats['total']
                }
                
                variants.append(variant_dict)
                vcf_stats['passed'] += 1
                vcf_stats['by_type'][variant_type] += 1
                vcf_stats['by_chrom'][chrom] += 1
                
                logging.debug(f"Parsed variant: {current_id} {variant_type} at position {pos}")
        
        # Sort variants by position (ascending)
        variants = sorted(variants, key=lambda v: v['pos'])
        
    except Exception as e:
        logging.error(f"Error parsing VCF file: {e}")
        if strict_hgvs:
            raise
    
    elapsed = time.time() - start_time
    logging.info(f"Finished parsing VCF in {elapsed:.2f}s: {len(variants)} variants")
    
    # Log statistics
    logging.info(f"VCF Statistics:")
    logging.info(f"  Total records: {vcf_stats['total']}")
    logging.info(f"  Passed variants: {vcf_stats['passed']}")
    logging.info(f"  Filtered variants: {vcf_stats['filtered']}")
    logging.info(f"  Invalid records: {vcf_stats['invalid']}")
    logging.info(f"  Multiallelic sites: {vcf_stats['multiallelic']}")
    logging.info(f"  Phased genotypes: {vcf_stats['phased']}")
    logging.info(f"  Unphased genotypes: {vcf_stats['unphased']}")
    logging.info(f"  Variant types:")
    for var_type, count in sorted(vcf_stats['by_type'].items()):
        logging.info(f"    {var_type}: {count}")
    
    return variants, samples, format_fields

def find_segment_for_position(position, segment_offsets):
    """Find which segment contains a given position and calculate relative position.
    
    Args:
        position (int): Genome position (1-based)
        segment_offsets (dict): Dictionary of segment offsets
        
    Returns:
        tuple: (segment_id, relative_position, orientation)
    """
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
    return None, None, None

def phase_variants(variants, samples):
    """Group variants by their phase blocks.
    
    Args:
        variants (list): List of variant dictionaries
        samples (list): List of sample names
        
    Returns:
        dict: Sample name -> Phase block ID -> List of variants
    """
    phased_variants = {sample: defaultdict(list) for sample in samples}
    unphased_variants = {sample: [] for sample in samples}
    
    for variant in variants:
        # Process each sample's genotype
        for genotype in variant.get('genotypes', []):
            sample = genotype['sample']
            phase_type = genotype.get('phase_type', 'unknown')
            
            # Only process variants that are present in this sample
            if not genotype.get('has_alt', False):
                continue
            
            if phase_type == 'phased':
                # Use chromosome as the block ID for phased variants (simplified model)
                # In a real implementation, we would use phase block information if available
                block_id = variant['chrom']
                phased_variants[sample][block_id].append(variant)
            else:
                # For unphased variants, keep them separate
                unphased_variants[sample].append(variant)
                
    return phased_variants, unphased_variants

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
    
    for variant in variants:
        seg_id, rel_pos, orientation = find_segment_for_position(variant['pos'], segment_offsets)
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
        # Sort by relative position (ascending)
        var_list.sort(key=lambda x: x[1])
        
        # Get the original sequence
        segment_seq = segments[seg_id]
        segment_seq_original = segment_seq  # Keep original for comparison
        
        logging.info(f"Modifying segment {seg_id} with {len(var_list)} variants")
        
        segment_log = []
        
        # Track applied variants and their positions to detect conflicts
        applied_variants = []
        
        # Apply variants from left to right (ascending position)
        for variant, rel_pos, orientation in var_list:
            # Adjust position based on previous modifications
            adjusted_pos = rel_pos
            for applied_var, applied_pos, applied_len_diff in applied_variants:
                if applied_pos < rel_pos:
                    adjusted_pos += applied_len_diff
            
            # Get the expected reference at the adjusted position
            expected_ref = segment_seq[adjusted_pos:adjusted_pos + len(variant['ref'])]
            
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
                    msg = f"Reference mismatch for {variant['id']} in segment {seg_id} at position {adjusted_pos}."
                    msg += f"\n  Expected: {variant['ref']}, Found: {expected_ref}"
                    msg += f"\n  Similarity score: {alignment_score:.2f}"
                    msg += f"\n  Original position: {rel_pos}, Adjusted position: {adjusted_pos}"
                    
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
                logging.debug(f"Applying SNP {variant['id']}: {variant['ref']} -> {variant['alt']} at position {adjusted_pos}")
                segment_seq = segment_seq[:adjusted_pos] + variant['alt'] + segment_seq[adjusted_pos + len(variant['ref']):]
                len_diff = len(variant['alt']) - len(variant['ref'])
            elif variant['type'] in ['INV', 'MNP']:
                logging.debug(f"Applying {variant['type']} {variant['id']}: {variant['ref']} -> {variant['alt']} at position {adjusted_pos}")
                segment_seq = segment_seq[:adjusted_pos] + variant['alt'] + segment_seq[adjusted_pos + len(variant['ref']):]
                len_diff = len(variant['alt']) - len(variant['ref'])
            elif variant['type'] == 'DEL':
                logging.debug(f"Applying DEL {variant['id']}: Removing {len(variant['ref'])} bp at position {adjusted_pos}")
                segment_seq = segment_seq[:adjusted_pos] + (variant['alt'] if variant['alt'] != '.' else '') + segment_seq[adjusted_pos + len(variant['ref']):]
                len_diff = len(variant['alt']) - len(variant['ref'])
            elif variant['type'] == 'INS':
                logging.debug(f"Applying INS {variant['id']}: Inserting {len(variant['alt'])-1} bp at position {adjusted_pos}")
                # Handle INS where ref is the base before insertion
                if len(variant['ref']) == 1:
                    segment_seq = segment_seq[:adjusted_pos + 1] + variant['alt'][1:] + segment_seq[adjusted_pos + 1:]
                    len_diff = len(variant['alt']) - 1
                else:
                    segment_seq = segment_seq[:adjusted_pos] + variant['alt'] + segment_seq[adjusted_pos + len(variant['ref']):]
                    len_diff = len(variant['alt']) - len(variant['ref'])
            else:
                logging.debug(f"Applying generic variant {variant['id']} ({variant['type']}): {variant['ref']} -> {variant['alt']} at position {adjusted_pos}")
                segment_seq = segment_seq[:adjusted_pos] + variant['alt'] + segment_seq[adjusted_pos + len(variant['ref']):]
                len_diff = len(variant['alt']) - len(variant['ref'])
            
            # Track this variant for future position adjustments
            applied_variants.append((variant, adjusted_pos, len_diff))
            
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

def apply_single_variant(segments, variant, segment_offsets, path_segments, allow_mismatches=False, match_threshold=0.8):
    """Apply a single variant to create a modified path.
    
    This function is used for unphased variants to create individual paths for each variant.
    
    Args:
        segments (dict): Dictionary of segment ID to sequence
        variant (dict): Single variant dictionary
        segment_offsets (dict): Dictionary of segment offsets
        path_segments (list): Reference path segments
        allow_mismatches (bool): Whether to allow reference mismatches
        match_threshold (float): Threshold for sequence similarity
        
    Returns:
        tuple: (modified_segment, new_path_segments, log_entry)
    """
    # Find which segment contains this variant
    seg_id, rel_pos, orientation = find_segment_for_position(variant['pos'], segment_offsets)
    if not seg_id:
        logging.warning(f"Could not map variant {variant['id']} to any segment")
        return None, None, None
    
    # Get the original sequence
    segment_seq = segments[seg_id]
    segment_seq_original = segment_seq
    
    # Check reference sequence
    expected_ref = segment_seq[rel_pos:rel_pos + len(variant['ref'])]
    if expected_ref != variant['ref']:
        # Try fuzzy matching
        alignment_score = 0
        if len(expected_ref) > 0 and len(variant['ref']) > 0:
            alignments = pairwise2.align.localms(
                expected_ref, variant['ref'], 
                2, -1, -2, -0.5,
                one_alignment_only=True
            )
            if alignments:
                alignment_score = alignments[0].score / (2 * max(len(expected_ref), len(variant['ref'])))
        
        if alignment_score < match_threshold:
            msg = f"Reference mismatch for {variant['id']} in segment {seg_id}."
            msg += f"\n  Expected: {variant['ref']}, Found: {expected_ref}"
            
            if not allow_mismatches:
                logging.error(msg)
                return None, None, None
            else:
                logging.warning(msg)
    
    # Apply the variant based on type
    if variant['type'] == 'SNP':
        segment_seq = segment_seq[:rel_pos] + variant['alt'] + segment_seq[rel_pos + len(variant['ref']):]
    elif variant['type'] in ['INV', 'MNP']:
        segment_seq = segment_seq[:rel_pos] + variant['alt'] + segment_seq[rel_pos + len(variant['ref']):]
    elif variant['type'] == 'DEL':
        segment_seq = segment_seq[:rel_pos] + (variant['alt'] if variant['alt'] != '.' else '') + segment_seq[rel_pos + len(variant['ref']):]
    elif variant['type'] == 'INS':
        if len(variant['ref']) == 1:
            segment_seq = segment_seq[:rel_pos + 1] + variant['alt'][1:] + segment_seq[rel_pos + 1:]
        else:
            segment_seq = segment_seq[:rel_pos] + variant['alt'] + segment_seq[rel_pos + len(variant['ref']):]
    else:
        segment_seq = segment_seq[:rel_pos] + variant['alt'] + segment_seq[rel_pos + len(variant['ref']):]
    
    # Create a new segment ID
    new_seg_id = f"{seg_id}_{variant['id']}"
    
    # Store the modified segment
    modified_segment = {
        'seg_id': seg_id,
        'new_id': new_seg_id,
        'sequence': segment_seq,
        'original_sequence': segment_seq_original,
        'variants': [variant],
        'length_diff': len(segment_seq) - len(segment_seq_original)
    }
    
    # Create a new path
    new_path_segments = []
    for curr_seg_id, orient in path_segments:
        if curr_seg_id == seg_id:
            new_path_segments.append((new_seg_id, orient))
        else:
            new_path_segments.append((curr_seg_id, orient))
    
    # Create log entry
    log_entry = {
        'variant_id': variant['id'],
        'segment_id': seg_id,
        'position': variant['pos'],
        'rel_position': rel_pos,
        'orientation': orientation,
        'type': variant['type'],
        'ref': variant['ref'],
        'alt': variant['alt'],
        'seq_before': segment_seq_original[max(0, rel_pos-5):min(len(segment_seq_original), rel_pos+len(variant['ref'])+5)],
        'seq_after': segment_seq[max(0, rel_pos-5):min(len(segment_seq), rel_pos+len(variant['alt'])+5)]
    }
    
    return modified_segment, new_path_segments, log_entry

def generate_new_gfa(original_segments, modified_segments, links, paths, new_path_segments, 
                output_file, alt_path_name="ALT", keep_original=True, add_metadata=True,
                sample_paths=None):
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
        sample_paths (dict): Optional dictionary of sample -> path name -> path segments
                            for multi-sample VCF support
    """
    start_time = time.time()
    logging.info(f"Writing new GFA to {output_file}")
    
    with open(output_file, 'w') as f:
        # Write header with metadata
        if add_metadata:
            timestamp = time.strftime("%Y-%m-%dT%H:%M:%S")
            f.write(f"H\tVN:Z:1.0\tTS:Z:{timestamp}\tPG:Z:hapli.tools.vcf_to_gfa\n")
        else:
            f.write("H\tVN:Z:1.0\n")
        
        # Collect all segments to write (original + modified)
        all_segments = {}
        if keep_original:
            all_segments.update(original_segments)
        
        # Collect segments from sample paths if provided
        if sample_paths:
            for sample, path_dict in sample_paths.items():
                for path_name, path_info in path_dict.items():
                    if 'segments' in path_info:
                        for seg in path_info['segments']:
                            if isinstance(seg, dict) and 'new_id' in seg:
                                all_segments[seg['new_id']] = seg['sequence']
        
        # Add segments from the main modified segments dict
        for seg_id, mod_seg in modified_segments.items():
            new_id = mod_seg['new_id']
            all_segments[new_id] = mod_seg['sequence']
        
        # Write all segments
        segments_written = 0
        for seg_id, sequence in all_segments.items():
            # Check if this is a modified segment
            is_modified = False
            mod_data = None
            
            # Look in main modified segments
            for orig_id, mod_seg in modified_segments.items():
                if mod_seg['new_id'] == seg_id:
                    is_modified = True
                    mod_data = mod_seg
                    break
            
            # Look in sample paths if not found
            if not is_modified and sample_paths:
                for sample, path_dict in sample_paths.items():
                    for path_name, path_info in path_dict.items():
                        if 'segments' in path_info:
                            for seg in path_info['segments']:
                                if isinstance(seg, dict) and seg.get('new_id') == seg_id:
                                    is_modified = True
                                    mod_data = seg
                                    break
            
            # Write segment with appropriate tags
            if is_modified and mod_data:
                # Add segment metadata
                tags = [f"LN:i:{len(sequence)}"]
                
                if add_metadata:
                    # Add variant metadata
                    variant_ids = ';'.join(v['id'] for v in mod_data['variants'])
                    tags.append(f"VA:Z:{variant_ids}")
                    
                    # Add length difference
                    if mod_data.get('length_diff', 0) != 0:
                        tags.append(f"LD:i:{mod_data['length_diff']}")
                    
                    # Add original segment reference if available
                    if 'seg_id' in mod_data:
                        tags.append(f"OR:Z:{mod_data['seg_id']}")
                    
                f.write(f"S\t{seg_id}\t{sequence}\t{' '.join(tags)}\n")
            else:
                # Write regular segment
                f.write(f"S\t{seg_id}\t{sequence}\tLN:i:{len(sequence)}\n")
            
            segments_written += 1
            logging.debug(f"Wrote segment {seg_id}")
        
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

        # Helper function to create links for a path
        def create_links_for_path(path):
            links_created = 0
            for i in range(len(path) - 1):
                # Handle both tuple and dict formats for segments
                if isinstance(path[i], tuple):
                    from_seg, from_dir = path[i]
                else:
                    from_seg, from_dir = path[i]['new_id'], path[i]['orientation']
                
                if isinstance(path[i+1], tuple):
                    to_seg, to_dir = path[i+1]
                else:
                    to_seg, to_dir = path[i+1]['new_id'], path[i+1]['orientation']
                
                # Create link
                tag_str = " SV:Z:variant" if add_metadata else ""
                f.write(f"L\t{from_seg}\t{from_dir}\t{to_seg}\t{to_dir}\t0M{tag_str}\n")
                links_created += 1
                logging.debug(f"Created link {from_seg}{from_dir} -> {to_seg}{to_dir}")
            
            return links_created
        
        # Create links for the main path
        ref_path_name = next(iter(paths.keys()))  # Use first path if REF not found
        if "REF" in paths:
            ref_path_name = "REF"
        
        # Create links for the main new path
        new_links_added = 0
        if new_path_segments:
            new_links_added = create_links_for_path(new_path_segments)
            links_written += new_links_added
            logging.info(f"Added {new_links_added} links for main path")
        
        # Create links for sample paths
        if sample_paths:
            for sample, path_dict in sample_paths.items():
                for path_name, path_info in path_dict.items():
                    if 'segments' in path_info:
                        sample_links = create_links_for_path(path_info['segments'])
                        links_written += sample_links
                        logging.info(f"Added {sample_links} links for {sample} path {path_name}")
        
        # Write original paths if requested
        paths_written = 0
        if keep_original:
            for path_name, path_segs in paths.items():
                segments_str = ','.join(f"{seg}{dir}" for seg, dir in path_segs)
                f.write(f"P\t{path_name}\t{segments_str}\t*\n")
                paths_written += 1
                logging.debug(f"Wrote original path {path_name}")
        
        # Write new main path with variants if it exists
        if new_path_segments:
            new_path_str = ','.join(f"{seg}{dir}" for seg, dir in new_path_segments)
            
            # Add path metadata tags
            tags = []
            if add_metadata:
                tags.append(f"VN:Z:variant_path")
                tags.append(f"SR:Z:{ref_path_name}")
                
            tag_str = "\t" + "\t".join(tags) if tags else ""
            
            f.write(f"P\t{alt_path_name}\t{new_path_str}\t*{tag_str}\n")
            paths_written += 1
            logging.info(f"Wrote main path {alt_path_name} with {len(new_path_segments)} segments")
        
        # Write sample paths if provided
        if sample_paths:
            for sample, path_dict in sample_paths.items():
                for path_name, path_info in path_dict.items():
                    if 'segments' in path_info:
                        # Convert path segments to string
                        if isinstance(path_info['segments'][0], tuple):
                            path_str = ','.join(f"{seg}{dir}" for seg, dir in path_info['segments'])
                        else:
                            path_str = ','.join(f"{seg['new_id']}{seg['orientation']}" for seg in path_info['segments'])
                        
                        # Add sample metadata
                        meta_tags = []
                        if add_metadata:
                            meta_tags.append(f"SM:Z:{sample}")
                            if 'phase' in path_info:
                                meta_tags.append(f"PH:Z:{path_info['phase']}")
                            if 'variant_ids' in path_info:
                                meta_tags.append(f"VA:Z:{';'.join(path_info['variant_ids'])}")
                            meta_tags.append(f"SR:Z:{ref_path_name}")
                        
                        meta_str = "\t" + "\t".join(meta_tags) if meta_tags else ""
                        
                        f.write(f"P\t{path_name}\t{path_str}\t*{meta_str}\n")
                        paths_written += 1
                        logging.debug(f"Wrote sample path {path_name} for {sample}")
    
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

def export_sequences(segments, original_path, modified_path, output_prefix, format='fasta'):
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
    
    alt_record = SeqRecord(
        Seq(alt_seq),
        id="ALT",
        name="alternate",
        description="Alternate sequence with variants"
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

def main():
    """Main function to run the VCF to GFA conversion."""
    # Check for required dependencies
    if not CYVCF2_AVAILABLE:
        logging.error("cyvcf2 is required for VCF parsing. Please install with: pip install cyvcf2")
        return 1
    
    parser = argparse.ArgumentParser(description='Convert VCF variants to new paths in GFA.')
    parser.add_argument('gfa_file', help='Input GFA file')
    parser.add_argument('vcf_file', help='Input VCF file with variants')
    parser.add_argument('output_file', help='Output GFA file')
    
    # Input/output options
    parser.add_argument('--ref-path', default='REF', help='Name of the reference path in GFA (default: %(default)s)')
    parser.add_argument('--alt-path', default='ALT', help='Name for the alternate path in output GFA (default: %(default)s)')
    parser.add_argument('--export-sequences', help='Export reference and alternate sequences to FASTA files with this prefix')
    parser.add_argument('--seq-format', choices=['fasta', 'genbank'], default='fasta',
                        help='Format for exported sequences (default: %(default)s)')
    
    # Variant processing options
    parser.add_argument('--allow-mismatches', action='store_true', help='Allow reference mismatches when applying variants')
    parser.add_argument('--match-threshold', type=float, default=0.8, 
                        help='Minimum sequence similarity for fuzzy matching (default: %(default)s)')
    parser.add_argument('--max-variants', type=int, help='Maximum number of variants to process')
    parser.add_argument('--chrom', help='Only process variants from this chromosome')
    
    # Sample and genotype options
    parser.add_argument('--samples', help='Comma-separated list of samples to process (default: all samples)')
    parser.add_argument('--respect-phasing', action='store_true', 
                       help='Respect phasing information in VCF (default is to create separate paths for each variant)')
    parser.add_argument('--haplotype-paths', action='store_true',
                       help='Create separate paths for each haplotype in each sample')
    parser.add_argument('--sample-prefix', default='SAMPLE_',
                       help='Prefix for sample-specific path names (default: %(default)s)')
    
    # GFA output options
    parser.add_argument('--no-original', action='store_true', help='Do not include original segments and paths in output')
    parser.add_argument('--no-metadata', action='store_true', help='Do not add metadata tags to GFA')
    
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
        variants, samples, format_fields = parse_vcf(
            args.vcf_file, 
            max_variants=args.max_variants,
            chrom_filter=args.chrom
        )
        
        if not variants:
            logging.warning("No variants found or passed filters")
            return 0
        
        # Filter samples if requested
        if args.samples:
            requested_samples = [s.strip() for s in args.samples.split(',')]
            samples = [s for s in samples if s in requested_samples]
            if not samples:
                logging.warning(f"None of the requested samples ({args.samples}) found in VCF")
                samples = []  # Process variants without sample information
        
        # Create combined path with all variants (traditional approach)
        modified_segments = {}
        variant_log = []
        new_path_segments = []
        
        if not samples or (not args.respect_phasing and not args.haplotype_paths):
            # Apply all variants to create a single ALT path (traditional approach)
            logging.info("Creating unified ALT path with all variants")
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
            
            # Generate new GFA with single ALT path
            generate_new_gfa(
                segments, 
                modified_segments, 
                links, 
                paths, 
                new_path_segments,
                args.output_file,
                alt_path_name=args.alt_path,
                keep_original=not args.no_original,
                add_metadata=not args.no_metadata
            )
        else:
            # Handle sample-specific paths
            sample_paths = {}
            all_modified_segments = {}
            all_variant_logs = []
            
            if args.respect_phasing:
                # Group variants by phase blocks
                logging.info("Grouping variants by phase blocks")
                phased_variants, unphased_variants = phase_variants(variants, samples)
                
                # Process phased variants
                for sample in samples:
                    sample_paths[sample] = {}
                    
                    # Process phased variants for this sample
                    for block_id, block_variants in phased_variants[sample].items():
                        if not block_variants:
                            continue
                        
                        logging.info(f"Processing phase block {block_id} with {len(block_variants)} variants for sample {sample}")
                        
                        # Create haplotype paths if requested
                        if args.haplotype_paths:
                            # Create separate paths for haplotype 1 and 2
                            for hap_idx in range(2):  # Assuming diploid
                                # Filter variants for this haplotype
                                hap_variants = []
                                for var in block_variants:
                                    for gt in var.get('genotypes', []):
                                        if gt['sample'] == sample and hap_idx < len(gt['allele_indices']):
                                            if gt['allele_indices'][hap_idx] == var['allele_index']:
                                                hap_variants.append(var)
                                                break
                                
                                if not hap_variants:
                                    continue
                                
                                path_name = f"{args.sample_prefix}{sample}_HAP{hap_idx+1}_{block_id}"
                                logging.info(f"Creating haplotype path {path_name} with {len(hap_variants)} variants")
                                
                                # Apply variants for this haplotype
                                hap_modified_segments, hap_path_segments, hap_log = apply_variants_to_segments(
                                    segments,
                                    hap_variants,
                                    segment_offsets,
                                    ref_path,
                                    allow_mismatches=args.allow_mismatches,
                                    match_threshold=args.match_threshold,
                                    parallel=args.parallel
                                )
                                
                                # Store results
                                all_modified_segments.update(hap_modified_segments)
                                all_variant_logs.extend(hap_log)
                                
                                # Add to sample paths
                                sample_paths[sample][path_name] = {
                                    'segments': hap_path_segments,
                                    'phase': f"haplotype_{hap_idx+1}",
                                    'block_id': block_id,
                                    'variant_ids': [v['id'] for v in hap_variants]
                                }
                                
                                # Add segments to main dictionary
                                for seg_id, mod_seg in hap_modified_segments.items():
                                    segments[mod_seg['new_id']] = mod_seg['sequence']
                        else:
                            # Create a single path for all variants in this phase block
                            path_name = f"{args.sample_prefix}{sample}_{block_id}"
                            logging.info(f"Creating path {path_name} with {len(block_variants)} phased variants")
                            
                            # Apply variants
                            block_modified_segments, block_path_segments, block_log = apply_variants_to_segments(
                                segments,
                                block_variants,
                                segment_offsets,
                                ref_path,
                                allow_mismatches=args.allow_mismatches,
                                match_threshold=args.match_threshold,
                                parallel=args.parallel
                            )
                            
                            # Store results
                            all_modified_segments.update(block_modified_segments)
                            all_variant_logs.extend(block_log)
                            
                            # Add to sample paths
                            sample_paths[sample][path_name] = {
                                'segments': block_path_segments,
                                'phase': 'phased',
                                'block_id': block_id,
                                'variant_ids': [v['id'] for v in block_variants]
                            }
                            
                            # Add segments to main dictionary
                            for seg_id, mod_seg in block_modified_segments.items():
                                segments[mod_seg['new_id']] = mod_seg['sequence']
                    
                    # Process unphased variants for this sample
                    for var in unphased_variants[sample]:
                        # Create a separate path for each unphased variant
                        path_name = f"{args.sample_prefix}{sample}_VAR_{var['id']}"
                        logging.info(f"Creating individual path {path_name} for unphased variant {var['id']}")
                        
                        # Apply single variant
                        mod_segment, var_path_segments, var_log = apply_single_variant(
                            segments,
                            var,
                            segment_offsets,
                            ref_path,
                            allow_mismatches=args.allow_mismatches,
                            match_threshold=args.match_threshold
                        )
                        
                        if mod_segment:
                            # Store the modified segment
                            segments[mod_segment['new_id']] = mod_segment['sequence']
                            
                            # Add to sample paths
                            sample_paths[sample][path_name] = {
                                'segments': var_path_segments,
                                'phase': 'unphased',
                                'variant_ids': [var['id']]
                            }
                            
                            all_variant_logs.append(var_log)
            else:
                # Create separate paths for each variant regardless of phasing
                logging.info("Creating individual paths for each variant (ignoring phasing)")
                
                for sample in samples:
                    sample_paths[sample] = {}
                    
                    # Find variants for this sample
                    sample_variants = []
                    for var in variants:
                        for gt in var.get('genotypes', []):
                            if gt['sample'] == sample and gt.get('has_alt', False):
                                sample_variants.append(var)
                                break
                    
                    logging.info(f"Processing {len(sample_variants)} variants for sample {sample}")
                    
                    # Create a separate path for each variant
                    for var in sample_variants:
                        path_name = f"{args.sample_prefix}{sample}_VAR_{var['id']}"
                        logging.info(f"Creating individual path {path_name} for variant {var['id']}")
                        
                        # Apply single variant
                        mod_segment, var_path_segments, var_log = apply_single_variant(
                            segments,
                            var,
                            segment_offsets,
                            ref_path,
                            allow_mismatches=args.allow_mismatches,
                            match_threshold=args.match_threshold
                        )
                        
                        if mod_segment:
                            # Store the modified segment
                            segments[mod_segment['new_id']] = mod_segment['sequence']
                            
                            # Add to sample paths
                            sample_paths[sample][path_name] = {
                                'segments': var_path_segments,
                                'phase': 'individual',
                                'variant_ids': [var['id']]
                            }
                            
                            all_variant_logs.append(var_log)
            
            # Generate new GFA with all paths
            generate_new_gfa(
                segments,
                all_modified_segments,  # This might be empty if we only created individual paths
                links,
                paths,
                [],  # No main new path
                args.output_file,
                alt_path_name=args.alt_path,
                keep_original=not args.no_original,
                add_metadata=not args.no_metadata,
                sample_paths=sample_paths
            )
            
            # Update for reporting
            variant_log = all_variant_logs
            modified_segments = all_modified_segments
        
        # Export sequences if requested
        if args.export_sequences:
            if not samples or (not args.respect_phasing and not args.haplotype_paths):
                # Export main paths
                ref_file, alt_file = export_sequences(
                    segments,
                    ref_path,
                    new_path_segments,
                    args.export_sequences,
                    format=args.seq_format
                )
            else:
                # Export sample-specific paths
                for sample in samples:
                    if sample in sample_paths and sample_paths[sample]:
                        sample_prefix = f"{args.export_sequences}_{sample}"
                        for path_name, path_info in sample_paths[sample].items():
                            if 'segments' in path_info:
                                path_segments = path_info['segments']
                                path_id = path_name.replace(args.sample_prefix, '')
                                out_prefix = f"{sample_prefix}_{path_id}"
                                
                                export_sequences(
                                    segments,
                                    ref_path,
                                    path_segments,
                                    out_prefix,
                                    format=args.seq_format
                                )
        
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
