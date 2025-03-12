#!/usr/bin/env python3
"""
GFA Variant Effect Analyzer

Analyzes variants in GFA files with sample haplotype information and reports their effects 
on genomic features defined in GFF3 files.

Usage:
    python gfa_variant_analyzer.py out.gfa test.gff3 --output variant_report.txt
"""

import argparse
import sys
import os
import logging
import time
import re
from collections import defaultdict, Counter
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

def parse_gfa(gfa_file):
    """Parse a GFA file and extract segments, paths, and variant information."""
    segments = {}
    links = []
    paths = {}
    variants = []
    samples = defaultdict(list)
    haplotypes = defaultdict(dict)
    
    start_time = time.time()
    logging.info(f"Parsing GFA file: {gfa_file}")
    
    # First pass: collect segment information
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
                
                # Extract variant info from segment name
                var_id = None
                var_type = None
                length_change = 0
                
                # Extract optional tags
                tags = {}
                for field in fields[3:]:
                    if ':' in field:
                        tag_parts = field.split(':', 2)
                        if len(tag_parts) >= 3:
                            tag_name, tag_type, tag_value = tag_parts
                            tags[tag_name] = (tag_type, tag_value)
                
                # Check if segment name contains variant info (e.g., S25_VAR16)
                if '_VAR' in seg_id:
                    parts = seg_id.split('_VAR')
                    if len(parts) > 1:
                        # The original segment ID is the part before _VAR
                        orig_seg_id = parts[0]
                        # Extract variant ID from segment name
                        var_id = 'VAR' + parts[1].split('_')[0]  # Handle multiple VAR tags
                    
                    # Try to infer variant type from segment name or tags
                    if "SNP" in seg_id or "SNV" in seg_id:
                        var_type = "SNP"
                    elif "DEL" in seg_id:
                        var_type = "DEL"
                    elif "INS" in seg_id:
                        var_type = "INS"
                    elif "DUP" in seg_id:
                        var_type = "DUP"
                    elif "INV" in seg_id:
                        var_type = "INV"
                
                # Extract variant info from tags
                if 'VA' in tags and tags['VA'][0] == 'Z':
                    variant_ids = tags['VA'][1].split(';')
                    if var_id is None and variant_ids:
                        var_id = variant_ids[0]
                
                # Extract length difference if available
                if 'LD' in tags and tags['LD'][0] == 'i':
                    try:
                        length_change = int(tags['LD'][1])
                    except ValueError:
                        pass
                
                # Store segment information
                segments[seg_id] = {
                    'sequence': sequence,
                    'length': len(sequence),
                    'tags': tags,
                    'variant_id': var_id,
                    'variant_type': var_type,
                    'length_change': length_change
                }
                
                # Add to variants list if this is a variant segment
                if var_id:
                    # Check if variant already exists
                    existing_var = next((v for v in variants if v['id'] == var_id), None)
                    if existing_var:
                        # Update existing variant entry
                        existing_var['segments'].append(seg_id)
                    else:
                        # Create new variant entry
                        variants.append({
                            'id': var_id,
                            'type': var_type,
                            'segments': [seg_id],
                            'length_change': length_change
                        })
            
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
    
    # Second pass: collect path information
    with open(gfa_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            record_type = fields[0]
            
            if record_type == 'P':  # Path
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
                    
                    path_segments.append((seg_id, orientation))
                
                # Extract tags if present
                path_tags = {}
                for field in fields[3:]:
                    if ':' in field:
                        tag_parts = field.split(':', 2)
                        if len(tag_parts) >= 3:
                            tag_name, tag_type, tag_value = tag_parts
                            path_tags[tag_name] = (tag_type, tag_value)
                
                # Extract sample and haplotype information
                sample_name = None
                haplotype = None
                variant_ids = []
                
                # First try to get info from tags
                if 'SM' in path_tags and path_tags['SM'][0] == 'Z':
                    sample_name = path_tags['SM'][1]
                
                if 'PH' in path_tags and path_tags['PH'][0] == 'Z':
                    haplotype = path_tags['PH'][1]
                
                if 'VA' in path_tags and path_tags['VA'][0] == 'Z':
                    variant_ids = path_tags['VA'][1].split(';')
                
                # If no sample info from tags, try to extract from path name
                if not sample_name:
                    # Try format like SAMPLE_Sample1_HAP1_1
                    match = re.match(r'SAMPLE_(.+)_HAP(\d+)_\d+', path_name)
                    if match:
                        sample_name = match.group(1)
                        hap_num = match.group(2)
                        haplotype = f"haplotype_{hap_num}"
                
                # If still no sample info but path contains variant segments, 
                # create a generic sample name
                if not sample_name:
                    has_variant_segments = any(segments.get(s[0], {}).get('variant_id') for s in path_segments)
                    if has_variant_segments and path_name != 'REF':
                        sample_name = f"Sample_{path_name}"
                        haplotype = "haplotype_1"
                
                # Register sample and haplotype
                if sample_name:
                    samples[sample_name].append(path_name)
                    if haplotype:
                        haplotypes[sample_name][haplotype] = path_name
                
                # Extract variant IDs from path segments if not already defined
                if not variant_ids:
                    for seg_id, _ in path_segments:
                        var_id = segments.get(seg_id, {}).get('variant_id')
                        if var_id and var_id not in variant_ids:
                            variant_ids.append(var_id)
                
                paths[path_name] = {
                    'segments': path_segments,
                    'tags': path_tags,
                    'sample': sample_name,
                    'haplotype': haplotype,
                    'variant_ids': variant_ids
                }
    
    # Find or create REF path
    if 'REF' not in paths:
        logging.info("REF path not found, looking for reference path")
        # Try to find a path without variants
        for name, path_data in paths.items():
            if not path_data['variant_ids'] and not path_data['sample']:
                paths['REF'] = path_data
                logging.info(f"Using {name} as reference path")
                break
    
    # If still no REF path, create one from non-variant segments
    if 'REF' not in paths:
        logging.info("Creating synthetic REF path from non-variant segments")
        ref_segments = []
        non_variant_segments = [seg_id for seg_id, data in segments.items() 
                               if not data.get('variant_id')]
        
        # Sort segments by any available numbering
        def extract_number(seg_id):
            match = re.search(r'S(\d+)', seg_id)
            return int(match.group(1)) if match else 0
        
        sorted_segments = sorted(non_variant_segments, key=extract_number)
        
        if sorted_segments:
            ref_segments = [(seg_id, '+') for seg_id in sorted_segments]
            paths['REF'] = {
                'segments': ref_segments,
                'tags': {},
                'sample': None,
                'haplotype': None,
                'variant_ids': []
            }
        else:
            logging.warning("Could not create REF path, no non-variant segments found")
    
    elapsed = time.time() - start_time
    logging.info(f"Finished parsing GFA in {elapsed:.2f}s: {len(segments)} segments, {len(paths)} paths, {len(variants)} variants")
    logging.info(f"Found {len(samples)} samples with paths")
    
    return segments, links, paths, variants, samples, haplotypes

def parse_gff3(gff_file):
    """Parse a GFF3 file into a list of features."""
    features = []
    
    start_time = time.time()
    logging.info(f"Parsing GFF3 file: {gff_file}")
    
    with open(gff_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            
            if len(fields) < 9:
                logging.warning(f"Line {line_num}: Invalid GFF3 record, missing fields")
                continue
                
            seqid = fields[0]
            source = fields[1]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            score = fields[5]
            strand = fields[6]
            phase = fields[7]
            
            # Parse attributes
            attributes = {}
            for attr in fields[8].split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attributes[key] = value
            
            feature = {
                'seqid': seqid,
                'source': source,
                'type': feature_type,
                'start': start,
                'end': end,
                'score': score,
                'strand': strand,
                'phase': phase,
                'attributes': attributes,
                'line_num': line_num
            }
            
            features.append(feature)
    
    elapsed = time.time() - start_time
    logging.info(f"Finished parsing GFF3 in {elapsed:.2f}s: {len(features)} features")
    
    # Build feature hierarchy
    feature_by_id = {}
    children_by_parent = defaultdict(list)
    
    for feature in features:
        if 'ID' in feature['attributes']:
            feature_id = feature['attributes']['ID']
            feature_by_id[feature_id] = feature
            
        if 'Parent' in feature['attributes']:
            parent_id = feature['attributes']['Parent']
            children_by_parent[parent_id].append(feature)
    
    # Attach children to their parents
    for feature in features:
        if 'ID' in feature['attributes']:
            feature_id = feature['attributes']['ID']
            feature['children'] = children_by_parent.get(feature_id, [])
    
    return features, feature_by_id, children_by_parent

def build_path_sequence(segments, path_segments):
    """Build a sequence from path segments."""
    sequence = ""
    segment_offsets = []
    current_offset = 0
    
    for seg_id, orientation in path_segments:
        if seg_id not in segments:
            logging.warning(f"Missing segment: {seg_id}")
            continue
            
        segment_seq = segments[seg_id]['sequence']
        
        # Handle reverse complement if needed
        if orientation == '-':
            seq_obj = Seq(segment_seq)
            segment_seq = str(seq_obj.reverse_complement())
            
        # Store segment offset information
        segment_offsets.append({
            'seg_id': seg_id,
            'orientation': orientation,
            'start': current_offset,
            'end': current_offset + len(segment_seq) - 1,
            'length': len(segment_seq),
            'variant_id': segments[seg_id].get('variant_id')
        })
        
        sequence += segment_seq
        current_offset += len(segment_seq)
    
    return sequence, segment_offsets

def identify_first_cds_in_genes(features, feature_by_id, children_by_parent):
    """Identify the first CDS feature in each gene to check for start codon changes."""
    # For each mRNA, find all its CDS children and mark the first one
    for feature in features:
        if feature['type'] == 'mRNA':
            mrna_id = feature['attributes'].get('ID')
            if not mrna_id:
                continue
                
            # Get all CDS features for this mRNA
            cds_features = [f for f in children_by_parent.get(mrna_id, []) if f['type'] == 'CDS']
            
            # Sort by position, considering strand
            if feature['strand'] == '+':
                cds_features.sort(key=lambda x: x['start'])
            else:
                cds_features.sort(key=lambda x: x['start'], reverse=True)
            
            # Mark the first CDS
            if cds_features:
                cds_features[0]['is_first_cds'] = True
    
    return features

def translate_sequence(nucleotide_seq, strand='+'):
    """Translate a nucleotide sequence to amino acids, considering strand."""
    if not nucleotide_seq:
        return ""
    
    # Handle reverse strand
    if strand == '-':
        seq_obj = Seq(nucleotide_seq)
        nucleotide_seq = str(seq_obj.reverse_complement())
    
    # Ensure sequence length is a multiple of 3
    remainder = len(nucleotide_seq) % 3
    if remainder > 0:
        nucleotide_seq = nucleotide_seq[:-remainder]
    
    if not nucleotide_seq:
        return ""
    
    # Translate to amino acids
    seq_obj = Seq(nucleotide_seq)
    return str(seq_obj.translate())

def compare_amino_acid_sequences(ref_aa, alt_aa):
    """Compare two amino acid sequences and identify changes."""
    changes = []
    substitutions = []
    inserted = []
    deleted = []
    
    # Handle case where sequences are identical
    if ref_aa == alt_aa:
        return {
            'changes': 0,
            'substitutions': 0,
            'insertions': 0,
            'deletions': 0,
            'details': []
        }
    
    # Find the common prefix and suffix
    prefix_len = 0
    while prefix_len < min(len(ref_aa), len(alt_aa)) and ref_aa[prefix_len] == alt_aa[prefix_len]:
        prefix_len += 1
    
    suffix_len = 0
    while (suffix_len < min(len(ref_aa), len(alt_aa)) - prefix_len and 
           ref_aa[len(ref_aa) - 1 - suffix_len] == alt_aa[len(alt_aa) - 1 - suffix_len]):
        suffix_len += 1
    
    # Extract the changed regions
    ref_changed = ref_aa[prefix_len:len(ref_aa) - suffix_len]
    alt_changed = alt_aa[prefix_len:len(alt_aa) - suffix_len]
    
    # Determine the type of changes
    if len(ref_changed) == 0 and len(alt_changed) > 0:
        # Pure insertion
        changes.append(f"Insertion of {len(alt_changed)} AA at position {prefix_len+1}")
        inserted = list(alt_changed)
    elif len(alt_changed) == 0 and len(ref_changed) > 0:
        # Pure deletion
        changes.append(f"Deletion of {len(ref_changed)} AA at position {prefix_len+1}")
        deleted = list(ref_changed)
    elif len(ref_changed) == len(alt_changed):
        # Pure substitution(s)
        for i in range(len(ref_changed)):
            if ref_changed[i] != alt_changed[i]:
                pos = prefix_len + i + 1
                changes.append(f"{ref_changed[i]}{pos}{alt_changed[i]}")
                substitutions.append((pos, ref_changed[i], alt_changed[i]))
    else:
        # Complex change (indel + substitution)
        changes.append(f"Complex change: {ref_changed} -> {alt_changed}")
        
        if len(ref_changed) < len(alt_changed):
            # Net insertion
            inserted = list(alt_changed)
            for aa in ref_changed:
                if aa in inserted:
                    inserted.remove(aa)
        elif len(ref_changed) > len(alt_changed):
            # Net deletion
            deleted = list(ref_changed)
            for aa in alt_changed:
                if aa in deleted:
                    deleted.remove(aa)
    
    return {
        'changes': len(changes),
        'substitutions': len(substitutions),
        'insertions': len(inserted),
        'deletions': len(deleted),
        'details': changes
    }

def identify_variants_by_position(segment_offsets, ref_segments, variant_segments):
    """Identify variants in an alternate sequence based on segment differences."""
    variants_by_position = []
    
    # If no reference segments provided, return empty list
    if not ref_segments:
        return variants_by_position
        
    # Create a map of reference positions
    ref_pos_map = {}
    current_pos = 0
    
    for seg_info in ref_segments:
        seg_id = seg_info['seg_id']
        length = seg_info.get('length', 0)
        ref_pos_map[seg_id] = (current_pos, current_pos + length - 1)
        current_pos += length
    
    # Find alternate segments that differ from reference
    for seg_info in segment_offsets:
        seg_id = seg_info['seg_id']
        start_pos = seg_info['start']
        length = seg_info['length']
        variant_id = seg_info.get('variant_id')
        
        # Only process segments with variant IDs
        if variant_id and seg_id in variant_segments:
            var_data = variant_segments[seg_id]
            ref_id = var_data.get('original_segment', '')
            
            # Find reference position for this variant
            if ref_id in ref_pos_map:
                ref_start, ref_end = ref_pos_map[ref_id]
                # Create a variant entry
                variant = {
                    'id': variant_id,
                    'type': var_data.get('variant_type', 'UNKNOWN'),
                    'pos': ref_start + 1,  # Convert to 1-based
                    'end': ref_end + 1,    # Convert to 1-based
                    'length_change': var_data.get('length_change', 0),
                    'segments': [seg_id]
                }
                variants_by_position.append(variant)
    
    return variants_by_position

def analyze_haplotype_differences(ref_seq, alt_seq, features, segment_offsets, variant_segments, segments):
    """Analyze differences between reference and alternate haplotype sequences."""
    feature_effects = []
    
    # Identify variants by comparing segment positions
    ref_segments = []
    for seg_info in segment_offsets:
        seg_id = seg_info['seg_id']
        if seg_id in segments and not segments[seg_id].get('variant_id'):
            ref_segments.append({
                'seg_id': seg_id,
                'orientation': seg_info['orientation'],
                'length': segments[seg_id]['length']
            })
    
    variants = identify_variants_by_position(segment_offsets, ref_segments, variant_segments)
    
    for feature in features:
        # Extract feature boundaries
        start = feature['start'] - 1  # Convert to 0-based
        end = feature['end'] - 1      # Convert to 0-based
        feature_type = feature['type']
        strand = feature['strand']
        
        # Skip if feature is outside sequence bounds
        if start >= len(ref_seq) or end >= len(ref_seq):
            logging.warning(f"Feature {feature['attributes'].get('ID', 'unknown')} is outside sequence bounds, skipping")
            continue
        
        # Extract feature sequence from reference
        ref_feature_seq = ref_seq[start:end+1]
        
        # Extract feature sequence from alternate
        alt_feature_seq = ""
        if end < len(alt_seq):
            alt_feature_seq = alt_seq[start:end+1]
        else:
            # Handle case where alternate sequence is shorter than reference
            if start < len(alt_seq):
                alt_feature_seq = alt_seq[start:]
            # Feature is beyond the end of the alternate sequence
            else:
                alt_feature_seq = ""
        
        # Find overlapping variants
        overlapping_variants = []
        for variant in variants:
            var_start = variant.get('pos', 0) - 1  # Convert to 0-based
            var_end = variant.get('end', var_start)
            
            # Check for overlap
            if not (var_end < start or var_start > end):
                overlapping_variants.append(variant)
        
        if ref_feature_seq == alt_feature_seq and not overlapping_variants:
            # No changes to this feature
            effect = {
                'feature': feature,
                'feature_type': feature_type,
                'ref_feature_seq': ref_feature_seq,
                'alt_feature_seq': alt_feature_seq,
                'variants': [],
                'effects': ['no_change'],
                'details': {}
            }
        else:
            # Analyze effects of changes
            effects = []
            effect_details = {}
            
            # Calculate basic effects
            length_change = len(alt_feature_seq) - len(ref_feature_seq)
            
            if length_change != 0:
                effects.append('length_change')
                effect_details['length_change'] = length_change
            
            # Check for more specific effects
            if feature_type == 'CDS':
                # Check if length change is a multiple of 3 (in-frame)
                if length_change % 3 == 0:
                    effects.append('in_frame_change')
                else:
                    effects.append('frame_shift')
                    effect_details['frame_shift'] = length_change % 3
                
                # Translate sequences to check for amino acid changes
                ref_aa = translate_sequence(ref_feature_seq, strand)
                alt_aa = translate_sequence(alt_feature_seq, strand)
                
                # Check for premature stop codons
                if '*' in alt_aa and (not '*' in ref_aa or alt_aa.index('*') < ref_aa.index('*') if '*' in ref_aa else True):
                    effects.append('premature_stop_codon')
                    effect_details['premature_stop_position'] = alt_aa.index('*') * 3
                
                # Check for amino acid changes
                aa_changes = compare_amino_acid_sequences(ref_aa, alt_aa)
                if aa_changes['changes'] > 0:
                    effects.append('amino_acid_change')
                    effect_details['amino_acid_changes'] = aa_changes
                
                # Check start codon disruption
                if feature.get('is_first_cds', False) and alt_aa and (not alt_aa.startswith('M')):
                    effects.append('start_codon_disruption')
            
            # Check for regulatory region effects
            if feature_type in ['promoter', 'terminator']:
                effects.append(f"{feature_type}_affected")
            
            # Check for splicing region effects
            if feature_type == 'exon':
                effects.append('splicing_affected')
            
            # Remove duplicates
            effects = list(set(effects))
            
            effect = {
                'feature': feature,
                'feature_type': feature_type,
                'ref_feature_seq': ref_feature_seq,
                'alt_feature_seq': alt_feature_seq,
                'variants': overlapping_variants,
                'effects': effects,
                'details': effect_details
            }
        
        feature_effects.append(effect)
    
    return feature_effects

def analyze_sample_haplotypes(feature_effects_by_haplotype, features):
    """Analyze differences between sample haplotypes to identify homo/heterozygous effects."""
    if len(feature_effects_by_haplotype) < 2:
        # Cannot determine zygosity with less than 2 haplotypes
        return {'homozygous': [], 'heterozygous': [], 'incomplete': True}
    
    # Get haplotype names
    haplotype_names = list(feature_effects_by_haplotype.keys())
    
    # Map feature IDs to effects for each haplotype
    feature_map = {}
    for hap_name, effects in feature_effects_by_haplotype.items():
        for effect in effects:
            feature_id = effect['feature']['attributes'].get('ID', 'unknown')
            if feature_id not in feature_map:
                feature_map[feature_id] = {}
            
            feature_map[feature_id][hap_name] = effect
    
    # Identify homozygous vs heterozygous effects
    homozygous = []
    heterozygous = []
    
    for feature_id, hap_effects in feature_map.items():
        # Skip if feature is not affected in all haplotypes
        if len(hap_effects) < len(haplotype_names):
            continue
        
        # Check if effects are identical across haplotypes
        effect_signatures = {}
        for hap_name, effect in hap_effects.items():
            # Create a signature of the effect for comparison
            sig = (
                tuple(sorted(effect['effects'])),
                effect['alt_feature_seq']
            )
            effect_signatures[hap_name] = sig
        
        # Check if all signatures are the same
        is_homozygous = len(set(effect_signatures.values())) == 1
        
        if is_homozygous:
            # Use the effect from the first haplotype
            first_hap = haplotype_names[0]
            effect = hap_effects[first_hap]
            effect['zygosity'] = 'homozygous'
            homozygous.append(effect)
        else:
            # Create a heterozygous effect entry
            # We'll use the effect from the first haplotype but note the differences
            first_hap = haplotype_names[0]
            effect = hap_effects[first_hap].copy()
            effect['zygosity'] = 'heterozygous'
            effect['haplotype_effects'] = hap_effects
            
            # Analyze differences
            diff_effects = set()
            for hap_name, hap_effect in hap_effects.items():
                for e in hap_effect['effects']:
                    diff_effects.add(e)
            
            effect['combined_effects'] = list(diff_effects)
            heterozygous.append(effect)
    
    return {
        'homozygous': homozygous,
        'heterozygous': heterozygous,
        'incomplete': False
    }

def generate_variant_effect_report(feature_effects, variants, outfile=None, sample_name=None, sample_effects=None):
    """Generate a comprehensive report of variant effects on features."""
    # Open output file if specified
    out = open(outfile, 'w') if outfile else sys.stdout
    
    # Write report header
    if sample_name:
        out.write(f"# Variant Effect Report for Sample: {sample_name}\n")
    else:
        out.write("# Variant Effect Report\n")
    out.write("#" + "=" * 79 + "\n\n")
    
    # Write variant summary
    out.write("## Variant Summary\n")
    out.write("-" * 80 + "\n")
    
    # Filter variants if sample-specific
    if sample_name and sample_effects:
        # Extract variant IDs from sample effects
        variant_ids = set()
        for effect_type in ['homozygous', 'heterozygous']:
            for effect in sample_effects.get(effect_type, []):
                for var in effect.get('variants', []):
                    variant_ids.add(var['id'])
        
        # Filter variants
        filtered_variants = [v for v in variants if v['id'] in variant_ids]
        if filtered_variants:
            for variant in filtered_variants:
                zygosity = ""
                # Determine zygosity for this variant
                for effect_type in ['homozygous', 'heterozygous']:
                    for effect in sample_effects.get(effect_type, []):
                        if any(v['id'] == variant['id'] for v in effect.get('variants', [])):
                            zygosity = effect_type.upper()
                            break
                    if zygosity:
                        break
                
                out.write(f"Variant: {variant['id']} ({variant.get('type', 'UNKNOWN')}){' - ' + zygosity if zygosity else ''}\n")
                if 'pos' in variant and variant['pos'] > 0:
                    out.write(f"Position: {variant.get('pos', 'unknown')}-{variant.get('end', 'unknown')}\n")
                out.write(f"Length Change: {variant.get('length_change', 'unknown')} bp\n")
                out.write(f"Affected Segments: {', '.join(variant.get('segments', []))}\n")
                out.write("\n")
        else:
            out.write("No variants for this sample.\n\n")
    else:
        # Write all variants
        for variant in variants:
            out.write(f"Variant: {variant['id']} ({variant.get('type', 'UNKNOWN')})\n")
            if 'pos' in variant and variant['pos'] > 0:
                out.write(f"Position: {variant.get('pos', 'unknown')}-{variant.get('end', 'unknown')}\n")
            out.write(f"Length Change: {variant.get('length_change', 'unknown')} bp\n")
            out.write(f"Affected Segments: {', '.join(variant.get('segments', []))}\n")
            out.write("\n")
    
    # Group features by type
    features_by_type = defaultdict(list)
    
    # Use sample_effects if provided, otherwise use feature_effects
    if sample_name and sample_effects and not sample_effects.get('incomplete'):
        effects_to_process = []
        effects_to_process.extend(sample_effects.get('homozygous', []))
        effects_to_process.extend(sample_effects.get('heterozygous', []))
    else:
        effects_to_process = feature_effects
    
    for effect in effects_to_process:
        features_by_type[effect['feature_type']].append(effect)
    
    # Write feature effect summary
    if sample_name:
        out.write(f"\n## Feature Effect Summary for Sample: {sample_name}\n")
    else:
        out.write("\n## Feature Effect Summary\n")
    out.write("-" * 80 + "\n")
    
    # Calculate statistics
    effect_counts = Counter()
    affected_feature_counts = Counter()
    zygosity_counts = Counter()
    
    for effect in effects_to_process:
        if effect['effects'] and effect['effects'] != ['no_change']:  # If feature is affected
            affected_feature_counts[effect['feature_type']] += 1
            for e in effect['effects']:
                effect_counts[e] += 1
            
            # Count zygosity if available
            if 'zygosity' in effect:
                zygosity_counts[effect['zygosity']] += 1
    
    # Write statistics
    out.write(f"Total Features Analyzed: {len(effects_to_process)}\n")
    out.write(f"Features Affected by Variants: {sum(affected_feature_counts.values())}\n\n")
    
    # Add zygosity statistics for sample-specific reports
    if sample_name and zygosity_counts:
        out.write("Zygosity Summary:\n")
        out.write(f"  Homozygous Effects: {zygosity_counts.get('homozygous', 0)}\n")
        out.write(f"  Heterozygous Effects: {zygosity_counts.get('heterozygous', 0)}\n\n")
    
    out.write("Affected Features by Type:\n")
    for feature_type, count in sorted(affected_feature_counts.items()):
        total = len(features_by_type[feature_type])
        percentage = (count / total * 100) if total > 0 else 0
        out.write(f"  {feature_type}: {count}/{total} ({percentage:.1f}%)\n")
    
    out.write("\nVariant Effects by Type:\n")
    for effect_type, count in sorted(effect_counts.items(), key=lambda x: x[1], reverse=True):
        out.write(f"  {effect_type}: {count}\n")
    
    # Write detailed feature effects
    if sample_name:
        out.write(f"\n\n## Detailed Feature Effects for Sample: {sample_name}\n")
    else:
        out.write("\n\n## Detailed Feature Effects\n")
    
    # Process each feature type
    for feature_type, effects in sorted(features_by_type.items()):
        affected_effects = [e for e in effects if e['effects'] and e['effects'] != ['no_change']]
        
        if not affected_effects:
            continue  # Skip feature types not affected by variants
        
        out.write(f"\n### {feature_type.upper()} Features\n")
        out.write("-" * 80 + "\n")
        
        for effect in affected_effects:
            feature = effect['feature']
            feature_id = feature['attributes'].get('ID', 'unknown')
            feature_name = feature['attributes'].get('Name', feature_id)
            
            out.write(f"\nFeature: {feature_name} ({feature_id})\n")
            out.write(f"Location: {feature['start']}-{feature['end']} ({feature['strand']})\n")
            
            # Add zygosity information if available
            if 'zygosity' in effect:
                out.write(f"Zygosity: {effect['zygosity'].upper()}\n")
            
            # List variants affecting this feature
            if effect.get('variants'):
                out.write("Affected by variants:\n")
                for variant in effect['variants']:
                    out.write(f"  - {variant['id']} ({variant.get('type', 'UNKNOWN')})\n")
            
            # List effects
            out.write("Effects:\n")
            
            # Use combined_effects for heterozygous variants if available
            if 'zygosity' in effect and effect['zygosity'] == 'heterozygous' and 'combined_effects' in effect:
                effect_list = effect['combined_effects']
            else:
                effect_list = effect['effects']
                
            for effect_type in sorted(effect_list):
                if effect_type == 'no_change':
                    out.write("  - No effect on this feature\n")
                    continue
                
                # Format effect with details if available
                detail_str = ""
                if effect_type in effect['details']:
                    detail = effect['details'][effect_type]
                    if isinstance(detail, dict):
                        detail_str = ": " + ", ".join(f"{k}={v}" for k, v in detail.items())
                    elif isinstance(detail, list):
                        detail_str = ": " + ", ".join(str(d) for d in detail)
                    else:
                        detail_str = f": {detail}"
                
                out.write(f"  - {effect_type}{detail_str}\n")
            
            # For heterozygous effects, show haplotype-specific differences
            if 'zygosity' in effect and effect['zygosity'] == 'heterozygous' and 'haplotype_effects' in effect:
                out.write("\nHaplotype-specific Effects:\n")
                for hap_name, hap_effect in effect['haplotype_effects'].items():
                    out.write(f"  {hap_name}:\n")
                    for e in sorted(hap_effect['effects']):
                        if e == 'no_change':
                            continue
                        out.write(f"    - {e}\n")
            
            # Add sequence details for CDS features
            if feature_type == 'CDS' and 'amino_acid_change' in effect['effects']:
                aa_changes = effect['details'].get('amino_acid_changes', {})
                if aa_changes:
                    out.write("\nAmino Acid Changes:\n")
                    for detail in aa_changes.get('details', []):
                        out.write(f"  - {detail}\n")
            
            # Write sequence changes
            if len(effect['ref_feature_seq']) <= 50 and len(effect['alt_feature_seq']) <= 50:
                out.write("\nSequence Changes:\n")
                out.write(f"  Reference: {effect['ref_feature_seq']}\n")
                out.write(f"  Alternate: {effect['alt_feature_seq']}\n")
            else:
                out.write("\nSequence Changes (truncated):\n")
                out.write(f"  Reference: {effect['ref_feature_seq'][:25]}...{effect['ref_feature_seq'][-25:]}\n")
                out.write(f"  Alternate: {effect['alt_feature_seq'][:25]}...{effect['alt_feature_seq'][-25:]}\n")
            
            # Add frame information for CDS features
            if feature_type == 'CDS':
                if 'frame_shift' in effect['effects']:
                    frame_shift = effect['details'].get('frame_shift', 'unknown')
                    out.write(f"\nFrame Shift: {frame_shift} bases\n")
                elif 'in_frame_change' in effect['effects']:
                    out.write("\nIn-frame change (multiple of 3 bases)\n")
                
                if 'premature_stop_codon' in effect['effects']:
                    stop_pos = effect['details'].get('premature_stop_position', 'unknown')
                    out.write(f"Premature stop codon introduced at position: {stop_pos}\n")
    
    # Close output file if opened
    if outfile:
        out.close()
        logging.info(f"Report written to {outfile}")

def generate_variant_effect_report(feature_effects, variants, outfile=None, sample_name=None, sample_effects=None):
    """Generate a comprehensive report of variant effects on features."""
    # Open output file if specified
    out = open(outfile, 'w') if outfile else sys.stdout
    
    # Write report header
    if sample_name:
        out.write(f"# Variant Effect Report for Sample: {sample_name}\n")
    else:
        out.write("# Variant Effect Report\n")
    out.write("#" + "=" * 79 + "\n\n")
    
    # Write variant summary
    out.write("## Variant Summary\n")
    out.write("-" * 80 + "\n")
    
    # Filter variants if sample-specific
    if sample_name and sample_effects:
        # Extract variant IDs from sample effects
        variant_ids = set()
        for effect_type in ['homozygous', 'heterozygous']:
            for effect in sample_effects.get(effect_type, []):
                for var in effect.get('variants', []):
                    variant_ids.add(var['id'])
        
        # Filter variants
        filtered_variants = [v for v in variants if v['id'] in variant_ids]
        if filtered_variants:
            for variant in filtered_variants:
                zygosity = ""
                # Determine zygosity for this variant
                for effect_type in ['homozygous', 'heterozygous']:
                    for effect in sample_effects.get(effect_type, []):
                        if any(v['id'] == variant['id'] for v in effect.get('variants', [])):
                            zygosity = effect_type.upper()
                            break
                    if zygosity:
                        break
                
                out.write(f"Variant: {variant['id']} ({variant.get('type', 'UNKNOWN')}){' - ' + zygosity if zygosity else ''}\n")
                if 'pos' in variant and variant['pos'] > 0:
                    out.write(f"Position: {variant.get('pos', 'unknown')}-{variant.get('end', 'unknown')}\n")
                out.write(f"Length Change: {variant.get('length_change', 'unknown')} bp\n")
                out.write(f"Affected Segments: {', '.join(variant.get('segments', []))}\n")
                out.write("\n")
        else:
            out.write("No variants for this sample.\n\n")
    else:
        # Write all variants
        for variant in variants:
            out.write(f"Variant: {variant['id']} ({variant.get('type', 'UNKNOWN')})\n")
            if 'pos' in variant and variant['pos'] > 0:
                out.write(f"Position: {variant.get('pos', 'unknown')}-{variant.get('end', 'unknown')}\n")
            out.write(f"Length Change: {variant.get('length_change', 'unknown')} bp\n")
            out.write(f"Affected Segments: {', '.join(variant.get('segments', []))}\n")
            out.write("\n")
    
    # Group features by type
    features_by_type = defaultdict(list)
    
    # Use sample_effects if provided, otherwise use feature_effects
    if sample_name and sample_effects and not sample_effects.get('incomplete'):
        effects_to_process = []
        effects_to_process.extend(sample_effects.get('homozygous', []))
        effects_to_process.extend(sample_effects.get('heterozygous', []))
    else:
        effects_to_process = feature_effects
    
    for effect in effects_to_process:
        features_by_type[effect['feature_type']].append(effect)
    
    # Write feature effect summary
    if sample_name:
        out.write(f"\n## Feature Effect Summary for Sample: {sample_name}\n")
    else:
        out.write("\n## Feature Effect Summary\n")
    out.write("-" * 80 + "\n")
    
    # Calculate statistics
    effect_counts = Counter()
    affected_feature_counts = Counter()
    zygosity_counts = Counter()
    
    for effect in effects_to_process:
        if effect['effects'] and effect['effects'] != ['no_change']:  # If feature is affected
            affected_feature_counts[effect['feature_type']] += 1
            for e in effect['effects']:
                effect_counts[e] += 1
            
            # Count zygosity if available
            if 'zygosity' in effect:
                zygosity_counts[effect['zygosity']] += 1
    
    # Write statistics
    out.write(f"Total Features Analyzed: {len(effects_to_process)}\n")
    out.write(f"Features Affected by Variants: {sum(affected_feature_counts.values())}\n\n")
    
    # Add zygosity statistics for sample-specific reports
    if sample_name and zygosity_counts:
        out.write("Zygosity Summary:\n")
        out.write(f"  Homozygous Effects: {zygosity_counts.get('homozygous', 0)}\n")
        out.write(f"  Heterozygous Effects: {zygosity_counts.get('heterozygous', 0)}\n\n")
    
    out.write("Affected Features by Type:\n")
    for feature_type, count in sorted(affected_feature_counts.items()):
        total = len(features_by_type[feature_type])
        percentage = (count / total * 100) if total > 0 else 0
        out.write(f"  {feature_type}: {count}/{total} ({percentage:.1f}%)\n")
    
    out.write("\nVariant Effects by Type:\n")
    for effect_type, count in sorted(effect_counts.items(), key=lambda x: x[1], reverse=True):
        out.write(f"  {effect_type}: {count}\n")
    
    # Write detailed feature effects
    if sample_name:
        out.write(f"\n\n## Detailed Feature Effects for Sample: {sample_name}\n")
    else:
        out.write("\n\n## Detailed Feature Effects\n")
    
    # Process each feature type
    for feature_type, effects in sorted(features_by_type.items()):
        affected_effects = [e for e in effects if e['effects'] and e['effects'] != ['no_change']]
        
        if not affected_effects:
            continue  # Skip feature types not affected by variants
        
        out.write(f"\n### {feature_type.upper()} Features\n")
        out.write("-" * 80 + "\n")
        
        for effect in affected_effects:
            feature = effect['feature']
            feature_id = feature['attributes'].get('ID', 'unknown')
            feature_name = feature['attributes'].get('Name', feature_id)
            
            out.write(f"\nFeature: {feature_name} ({feature_id})\n")
            out.write(f"Location: {feature['start']}-{feature['end']} ({feature['strand']})\n")
            
            # Add zygosity information if available
            if 'zygosity' in effect:
                out.write(f"Zygosity: {effect['zygosity'].upper()}\n")
            
            # List variants affecting this feature
            if effect.get('variants'):
                out.write("Affected by variants:\n")
                for variant in effect['variants']:
                    out.write(f"  - {variant['id']} ({variant.get('type', 'UNKNOWN')})\n")
            
            # List effects
            out.write("Effects:\n")
            
            # Use combined_effects for heterozygous variants if available
            if 'zygosity' in effect and effect['zygosity'] == 'heterozygous' and 'combined_effects' in effect:
                effect_list = effect['combined_effects']
            else:
                effect_list = effect['effects']
                
            for effect_type in sorted(effect_list):
                if effect_type == 'no_change':
                    out.write("  - No effect on this feature\n")
                    continue
                
                # Format effect with details if available
                detail_str = ""
                if effect_type in effect['details']:
                    detail = effect['details'][effect_type]
                    if isinstance(detail, dict):
                        detail_str = ": " + ", ".join(f"{k}={v}" for k, v in detail.items())
                    elif isinstance(detail, list):
                        detail_str = ": " + ", ".join(str(d) for d in detail)
                    else:
                        detail_str = f": {detail}"
                
                out.write(f"  - {effect_type}{detail_str}\n")
            
            # For heterozygous effects, show haplotype-specific differences
            if 'zygosity' in effect and effect['zygosity'] == 'heterozygous' and 'haplotype_effects' in effect:
                out.write("\nHaplotype-specific Effects:\n")
                for hap_name, hap_effect in effect['haplotype_effects'].items():
                    out.write(f"  {hap_name}:\n")
                    for e in sorted(hap_effect['effects']):
                        if e == 'no_change':
                            continue
                        out.write(f"    - {e}\n")
            
            # Add sequence details for CDS features
            if feature_type == 'CDS' and 'amino_acid_change' in effect['effects']:
                aa_changes = effect['details'].get('amino_acid_changes', {})
                if aa_changes:
                    out.write("\nAmino Acid Changes:\n")
                    for detail in aa_changes.get('details', []):
                        out.write(f"  - {detail}\n")
            
            # Write sequence changes
            if len(effect['ref_feature_seq']) <= 50 and len(effect['alt_feature_seq']) <= 50:
                out.write("\nSequence Changes:\n")
                out.write(f"  Reference: {effect['ref_feature_seq']}\n")
                out.write(f"  Alternate: {effect['alt_feature_seq']}\n")
            else:
                out.write("\nSequence Changes (truncated):\n")
                out.write(f"  Reference: {effect['ref_feature_seq'][:25]}...{effect['ref_feature_seq'][-25:]}\n")
                out.write(f"  Alternate: {effect['alt_feature_seq'][:25]}...{effect['alt_feature_seq'][-25:]}\n")
            
            # Add frame information for CDS features
            if feature_type == 'CDS':
                if 'frame_shift' in effect['effects']:
                    frame_shift = effect['details'].get('frame_shift', 'unknown')
                    out.write(f"\nFrame Shift: {frame_shift} bases\n")
                elif 'in_frame_change' in effect['effects']:
                    out.write("\nIn-frame change (multiple of 3 bases)\n")
                
                if 'premature_stop_codon' in effect['effects']:
                    stop_pos = effect['details'].get('premature_stop_position', 'unknown')
                    out.write(f"Premature stop codon introduced at position: {stop_pos}\n")
    
    # Close output file if opened
    if outfile:
        out.close()
        logging.info(f"Report written to {outfile}")

def generate_sample_reports(segments, paths, variants, features, outfile_prefix=None):
    """Generate sample-specific variant effect reports."""
    # Identify samples with paths
    samples = defaultdict(list)
    haplotypes = defaultdict(dict)
    
    for path_name, path_data in paths.items():
        sample_name = path_data.get('sample')
        if sample_name:
            samples[sample_name].append(path_name)
            
            # Track haplotype paths
            haplotype = path_data.get('haplotype')
            if haplotype and haplotype.startswith('haplotype_'):
                haplotypes[sample_name][haplotype] = path_name
    
    if not samples:
        logging.warning("No sample-specific paths found in GFA")
        return {}
    
    logging.info(f"Found {len(samples)} samples with paths")
    sample_reports = {}
    
    # Get reference path
    ref_path_name = 'REF'
    if ref_path_name not in paths:
        logging.error("REF path not found in GFA")
        return {}
    
    ref_path_segments = paths[ref_path_name]['segments']
    ref_seq, _ = build_path_sequence(segments, ref_path_segments)
    
    # Process each sample
    for sample_name, path_list in samples.items():
        logging.info(f"Generating report for sample: {sample_name}")
        
        # Determine if we have phased haplotypes for this sample
        sample_haplotypes = haplotypes.get(sample_name, {})
        has_phased_haplotypes = len(sample_haplotypes) >= 2
        
        if has_phased_haplotypes:
            logging.info(f"Found {len(sample_haplotypes)} phased haplotypes for sample {sample_name}")
            
            # Process each haplotype
            feature_effects_by_haplotype = {}
            variant_segments = {}
            
            # First collect variant segments
            for seg_id, seg_data in segments.items():
                if seg_data.get('variant_id'):
                    variant_segments[seg_id] = seg_data
            
            for hap_name, path_name in sample_haplotypes.items():
                logging.info(f"Processing haplotype: {hap_name} (path: {path_name})")
                
                # Build path sequence
                path_segments = paths[path_name]['segments']
                hap_seq, segment_offsets = build_path_sequence(segments, path_segments)
                
                # Analyze haplotype differences
                hap_effects = analyze_haplotype_differences(
                    ref_seq, 
                    hap_seq, 
                    features, 
                    segment_offsets, 
                    variant_segments,
                    segments
                )
                
                feature_effects_by_haplotype[hap_name] = hap_effects
            
            # Analyze homozygous vs heterozygous effects
            zygosity_effects = analyze_sample_haplotypes(feature_effects_by_haplotype, features)
            
            # Generate report for this sample
            outfile = f"{outfile_prefix}_{sample_name}.txt" if outfile_prefix else None
            generate_variant_effect_report([], variants, outfile, sample_name, zygosity_effects)
            
            sample_reports[sample_name] = outfile
        else:
            # Process single haplotype
            logging.info(f"Only one haplotype found for sample {sample_name}, processing as single path")
            
            # Get the first path for this sample
            path_name = path_list[0]
            path_segments = paths[path_name]['segments']
            
            # Collect variant segments
            variant_segments = {}
            for seg_id, seg_data in segments.items():
                if seg_data.get('variant_id'):
                    variant_segments[seg_id] = seg_data
            
            # Build path sequence
            alt_seq, segment_offsets = build_path_sequence(segments, path_segments)
            
            # Analyze differences
            feature_effects = analyze_haplotype_differences(
                ref_seq, 
                alt_seq, 
                features, 
                segment_offsets, 
                variant_segments,
                segments
            )
            
            # Generate report
            outfile = f"{outfile_prefix}_{sample_name}.txt" if outfile_prefix else None
            generate_variant_effect_report(feature_effects, variants, outfile, sample_name)
            
            sample_reports[sample_name] = outfile
    
    return sample_reports

def main():
    """Main function to run the variant effect report generation."""
    parser = argparse.ArgumentParser(description='Generate a variant effect report based on GFA and GFF3 files.')
    parser.add_argument('gfa_file', help='Input GFA file with variants')
    parser.add_argument('gff_file', help='Input GFF3 file with feature annotations')
    parser.add_argument('--output', help='Output report file (default: stdout)')
    
    # Sample options
    parser.add_argument('--sample-reports', action='store_true', 
                       help='Generate individual reports for each sample')
    parser.add_argument('--output-prefix', help='Prefix for sample report files')
    parser.add_argument('--samples', help='Comma-separated list of sample names to process (default: all)')
    
    # Debug and logging options
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output without full debug')
    parser.add_argument('--log-file', help='Write log to this file')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(debug=args.debug, log_file=args.log_file, verbose=args.verbose)
    
    try:
        # Parse input files
        segments, links, paths, variants, samples, haplotypes = parse_gfa(args.gfa_file)
        features, feature_by_id, children_by_parent = parse_gff3(args.gff_file)
        
        # Mark first CDS in each gene for start codon analysis
        features = identify_first_cds_in_genes(features, feature_by_id, children_by_parent)
        
        # Filter samples if requested
        if args.samples and args.sample_reports:
            requested_samples = [s.strip() for s in args.samples.split(',')]
            filtered_samples = {name: paths for name, paths in samples.items() if name in requested_samples}
            if not filtered_samples:
                logging.warning(f"None of the requested samples found in GFA")
            samples = filtered_samples
        
        # Generate sample-specific reports if requested
        if args.sample_reports:
            if samples:
                sample_reports = generate_sample_reports(
                    segments, 
                    paths, 
                    variants, 
                    features, 
                    outfile_prefix=args.output_prefix
                )
                logging.info(f"Generated {len(sample_reports)} sample reports")
            else:
                logging.warning("No samples found in GFA, cannot generate sample reports")
        
        # Generate main report if requested
        if args.output or not args.sample_reports:
            # Get reference path
            ref_path_name = 'REF'
            if ref_path_name not in paths:
                logging.error("REF path not found in GFA")
                return 1
            
            # Collect variant segments
            variant_segments = {}
            for seg_id, seg_data in segments.items():
                if seg_data.get('variant_id'):
                    variant_segments[seg_id] = seg_data
            
            # Build reference sequence
            ref_path_segments = paths[ref_path_name]['segments']
            ref_seq, ref_offsets = build_path_sequence(segments, ref_path_segments)
            
            # Process each sample or first available sample if none specified
            processed_sample = False
            
            for sample_name, sample_haplotypes in haplotypes.items():
                if not args.samples or sample_name in args.samples.split(','):
                    logging.info(f"Processing sample: {sample_name}")
                    
                    # Process first haplotype
                    hap_name, path_name = next(iter(sample_haplotypes.items()))
                    
                    # Build path sequence
                    path_segments = paths[path_name]['segments']
                    alt_seq, segment_offsets = build_path_sequence(segments, path_segments)
                    
                    # Analyze differences
                    feature_effects = analyze_haplotype_differences(
                        ref_seq, 
                        alt_seq, 
                        features, 
                        segment_offsets, 
                        variant_segments,
                        segments
                    )
                    
                    # Generate report
                    outfile = f"{args.output}_{sample_name}.txt" if args.output else None
                    generate_variant_effect_report(feature_effects, variants, outfile, sample_name)
                    
                    processed_sample = True
                    break
            
            # If no samples processed, use the first available path
            if not processed_sample:
                for path_name, path_data in paths.items():
                    if path_name != 'REF' and 'segments' in path_data:
                        # Build path sequence
                        path_segments = path_data['segments']
                        alt_seq, segment_offsets = build_path_sequence(segments, path_segments)
                        
                        # Analyze differences
                        feature_effects = analyze_haplotype_differences(
                            ref_seq, 
                            alt_seq, 
                            features, 
                            segment_offsets, 
                            variant_segments,
                            segments
                        )
                        
                        # Generate report
                        generate_variant_effect_report(feature_effects, variants, args.output)
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
