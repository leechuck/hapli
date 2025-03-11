#!/usr/bin/env python
"""
Generate a comprehensive variant effect report by analyzing GFA and GFF3 files.
This script identifies and reports all functional effects of variants on genomic features.

Author: Claude
Date: 2025-03-11
"""

import argparse
import re
import sys
import os
import logging
import time
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
    """Parse a GFA file into segments, links, paths, and extract variant information."""
    segments = {}
    links = []
    paths = {}
    variants = []
    
    start_time = time.time()
    logging.info(f"Parsing GFA file: {gfa_file}")
    
    with open(gfa_file, 'r') as f:
        in_variant_annotations = False
        
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
                
            # Parse variant annotations from comments
            if line.startswith('# Variant Annotations'):
                in_variant_annotations = True
                continue
            elif line.startswith('# End of Variant Annotations'):
                in_variant_annotations = False
                continue
            elif in_variant_annotations and line.startswith('# VAR'):
                # Extract variant information from comment line
                parts = line.split('\t')
                if len(parts) >= 6:
                    var_id = parts[1]
                    var_type = parts[2]
                    pos = int(parts[3])
                    ref = parts[4]
                    alt = parts[5]
                    
                    # Parse additional info if available
                    info = {}
                    if len(parts) > 6:
                        for item in parts[6:]:
                            if '=' in item:
                                key, value = item.split('=', 1)
                                info[key] = value
                    
                    variant = {
                        'id': var_id,
                        'type': var_type,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'info': info,
                        'end': pos + len(ref) - 1,
                        'length_change': len(alt) - len(ref)
                    }
                    
                    variants.append(variant)
                continue
            elif line.startswith('#'):
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
                
                # Extract variant information from segment tags
                var_info = []
                if 'VA' in tags and tags['VA'][0] == 'Z':
                    var_ids = tags['VA'][1].split(';')
                    
                    # Look for variant details in tags
                    for i in range(1, 10):  # Assume maximum 9 variants per segment
                        var_prefix = f"V{i}"
                        
                        if f"{var_prefix}ID" in tags and tags[f"{var_prefix}ID"][0] == 'Z':
                            var_id = tags[f"{var_prefix}ID"][1]
                            var_type = tags.get(f"{var_prefix}TYPE", ('Z', 'UNKNOWN'))[1]
                            var_pos = int(tags.get(f"{var_prefix}POS", ('i', '0'))[1])
                            var_ref = tags.get(f"{var_prefix}REF", ('Z', ''))[1]
                            var_alt = tags.get(f"{var_prefix}ALT", ('Z', ''))[1]
                            
                            var_detail = {
                                'id': var_id,
                                'type': var_type,
                                'pos': var_pos,
                                'ref': var_ref,
                                'alt': var_alt,
                                'end': var_pos + len(var_ref) - 1,
                                'length_change': len(var_alt) - len(var_ref)
                            }
                            
                            var_info.append(var_detail)
                            
                            # Add to global variants if not already present
                            if not any(v['id'] == var_id for v in variants):
                                variants.append(var_detail)
                        else:
                            break
                
                segments[seg_id] = {
                    'sequence': sequence,
                    'length': len(sequence),
                    'tags': tags,
                    'variants': var_info
                }
                
                # If this is a modified segment, note the original segment
                if 'OR' in tags and tags['OR'][0] == 'Z':
                    original_seg = tags['OR'][1]
                    segments[seg_id]['original_segment'] = original_seg
                
                logging.debug(f"Parsed segment: {seg_id} (length: {len(sequence)}, variants: {len(var_info)})")
                
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
                    
                    path_segments.append((seg_id, orientation))
                
                # Extract tags if present
                path_tags = {}
                for field in fields[3:]:
                    if ':' in field:
                        tag_parts = field.split(':', 2)
                        if len(tag_parts) >= 3:
                            tag_name, tag_type, tag_value = tag_parts
                            path_tags[tag_name] = (tag_type, tag_value)
                
                paths[path_name] = {
                    'segments': path_segments,
                    'tags': path_tags
                }
                logging.debug(f"Parsed path: {path_name} with {len(path_segments)} segments")
    
    elapsed = time.time() - start_time
    logging.info(f"Finished parsing GFA in {elapsed:.2f}s: {len(segments)} segments, {len(links)} links, {len(paths)} paths, {len(variants)} variants")
    
    return segments, links, paths, variants

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
            logging.error(f"Missing segment: {seg_id}")
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
            'length': len(segment_seq)
        })
        
        sequence += segment_seq
        current_offset += len(segment_seq)
    
    return sequence, segment_offsets

def map_features_to_paths(features, ref_seq, alt_seq, variants):
    """Map features from reference to alternate path and identify effects."""
    feature_effects = []
    
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
        
        # Find overlapping variants
        overlapping_variants = []
        for variant in variants:
            var_start = variant['pos'] - 1  # Convert to 0-based
            var_end = var_start + len(variant['ref']) - 1
            
            # Check for overlap
            if not (var_end < start or var_start > end):
                overlapping_variants.append(variant)
        
        if not overlapping_variants:
            # No variants overlap this feature
            effect = {
                'feature': feature,
                'feature_type': feature_type,
                'ref_feature_seq': ref_feature_seq,
                'alt_feature_seq': ref_feature_seq,  # No change
                'variants': [],
                'effects': ['no_change'],
                'details': {}
            }
        else:
            # Identify impact of overlapping variants
            effects = []
            effect_details = {}
            
            # Find and apply all variants to the feature sequence
            alt_feature_seq = apply_variants_to_sequence(ref_feature_seq, overlapping_variants, start)
            
            # Calculate basic effects
            length_change = len(alt_feature_seq) - len(ref_feature_seq)
            
            if length_change != 0:
                effects.append('length_change')
                effect_details['length_change'] = length_change
            
            # Analyze more specific effects
            for variant in overlapping_variants:
                var_effects = analyze_variant_effect(variant, feature, start, end, ref_feature_seq, alt_feature_seq)
                effects.extend(var_effects['effects'])
                
                # Merge effect details
                for k, v in var_effects['details'].items():
                    if k in effect_details:
                        # If already exists, convert to list or append to list
                        if not isinstance(effect_details[k], list):
                            effect_details[k] = [effect_details[k]]
                        if isinstance(v, list):
                            effect_details[k].extend(v)
                        else:
                            effect_details[k].append(v)
                    else:
                        effect_details[k] = v
            
            # Check for frame effects in coding features
            if feature_type == 'CDS':
                # Get the original phase
                phase = int(feature['phase']) if feature['phase'] != '.' else 0
                
                # Check if length change is a multiple of 3 for in-frame changes
                if length_change % 3 == 0:
                    effects.append('in_frame_change')
                else:
                    effects.append('frame_shift')
                    effect_details['frame_shift'] = length_change % 3
                
                # Check for premature stop codons
                ref_cds = ref_feature_seq
                alt_cds = alt_feature_seq
                
                # Adjust for phase
                if phase > 0:
                    ref_cds = ref_cds[phase:]
                    # For alt_cds, we need to be careful about variant positions
                    alt_cds = alt_cds[phase:]
                
                # Translate and check for stop codons
                ref_aa = translate_sequence(ref_cds, strand)
                alt_aa = translate_sequence(alt_cds, strand)
                
                # Check for premature stop codons in alt sequence
                if '*' in alt_aa and (not '*' in ref_aa or alt_aa.index('*') < ref_aa.index('*')):
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
            
            # Remove duplicates from effects list
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

def apply_variants_to_sequence(sequence, variants, offset):
    """Apply variants to a sequence, accounting for position shifts."""
    # Sort variants by position (ascending)
    sorted_variants = sorted(variants, key=lambda v: v['pos'])
    
    # Apply variants in order of position
    result = sequence
    position_shift = 0
    
    for variant in sorted_variants:
        var_start = variant['pos'] - 1 - offset  # Convert to 0-based and adjust for feature offset
        var_end = var_start + len(variant['ref'])
        
        # Adjust for previous variants
        var_start += position_shift
        var_end += position_shift
        
        # Ensure variant is within bounds
        if var_start < 0:
            var_start = 0
        if var_end > len(result):
            var_end = len(result)
        
        # Apply the variant
        if var_start <= len(result) and var_start >= 0:
            if variant['type'] == 'SNP':
                result = result[:var_start] + variant['alt'] + result[var_end:]
            elif variant['type'] == 'DEL':
                result = result[:var_start] + variant['alt'] + result[var_end:]
            elif variant['type'] == 'INS':
                result = result[:var_start] + variant['alt'] + result[var_start:]
            elif variant['type'] == 'INV':
                # For inversion, we reverse complement the reference sequence
                result = result[:var_start] + variant['alt'] + result[var_end:]
            else:
                # Default handling for other variant types
                result = result[:var_start] + variant['alt'] + result[var_end:]
            
            # Update position shift
            position_shift += len(variant['alt']) - (var_end - var_start)
    
    return result

def analyze_variant_effect(variant, feature, feature_start, feature_end, ref_feature_seq, alt_feature_seq):
    """Analyze the specific effect of a variant on a feature."""
    effects = []
    details = {}
    
    variant_type = variant['type']
    var_start = variant['pos'] - 1  # Convert to 0-based
    var_end = var_start + len(variant['ref']) - 1
    
    # Basic effect based on variant type
    if variant_type == 'SNP':
        effects.append('substitution')
        details['substitution'] = f"{variant['ref']} to {variant['alt']}"
        
        # Check if SNP changes amino acids
        if feature['type'] == 'CDS':
            codon_pos = (var_start - feature_start) % 3
            details['codon_position'] = codon_pos
    
    elif variant_type == 'DEL':
        effects.append('deletion')
        details['deletion_length'] = len(variant['ref']) - len(variant['alt'])
        
        # Check if deletion removes a full codon
        if feature['type'] == 'CDS':
            is_codon_aligned = (var_start - feature_start) % 3 == 0 and len(variant['ref']) % 3 == 0
            details['codon_aligned_deletion'] = is_codon_aligned
    
    elif variant_type == 'INS':
        effects.append('insertion')
        details['insertion_length'] = len(variant['alt']) - len(variant['ref'])
        
        # Check if insertion is codon-aligned
        if feature['type'] == 'CDS':
            is_codon_aligned = (var_start - feature_start) % 3 == 0 and len(variant['alt']) % 3 == 0
            details['codon_aligned_insertion'] = is_codon_aligned
    
    elif variant_type == 'INV':
        effects.append('inversion')
        details['inversion_length'] = len(variant['ref'])
    
    elif variant_type == 'DUP':
        effects.append('duplication')
        details['duplication_length'] = len(variant['alt']) - len(variant['ref'])
        
        # Detect tandem duplications
        if len(variant['alt']) > len(variant['ref']):
            # Check if the extra sequence is a repeat of part of the reference
            extra_seq = variant['alt'][len(variant['ref']):]
            if extra_seq in ref_feature_seq:
                effects.append('tandem_duplication')
                details['tandem_repeat_count'] = len(variant['alt']) // len(variant['ref'])
    
    # Analyze location of variant relative to feature
    if var_start <= feature_start and var_end >= feature_end:
        effects.append('feature_span')  # Variant spans entire feature
    elif var_start <= feature_start < var_end < feature_end:
        effects.append('feature_start_disruption')  # Variant disrupts start
    elif feature_start < var_start < feature_end and var_end >= feature_end:
        effects.append('feature_end_disruption')  # Variant disrupts end
    elif feature_start < var_start and var_end < feature_end:
        effects.append('feature_internal')  # Variant is internal to feature
    
    # Special analysis for regulatory regions
    if feature['type'] in ['promoter', 'terminator']:
        effects.append(f"{feature['type']}_affected")
    
    # Special analysis for splicing regions (if near exon boundaries)
    if feature['type'] == 'exon':
        splice_region_size = 3  # nucleotides at exon boundaries considered splice regions
        
        if abs(var_start - feature_start) < splice_region_size or abs(var_end - feature_start) < splice_region_size:
            effects.append('splice_acceptor_affected')
        
        if abs(var_start - feature_end) < splice_region_size or abs(var_end - feature_end) < splice_region_size:
            effects.append('splice_donor_affected')
    
    return {'effects': effects, 'details': details}

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

def generate_variant_effect_report(feature_effects, variants, outfile=None):
    """Generate a comprehensive report of variant effects on features."""
    # Open output file if specified
    out = open(outfile, 'w') if outfile else sys.stdout
    
    # Write report header
    out.write("# Variant Effect Report\n")
    out.write("#" + "=" * 79 + "\n\n")
    
    # Write variant summary
    out.write("## Variant Summary\n")
    out.write("-" * 80 + "\n")
    
    for variant in variants:
        out.write(f"Variant: {variant['id']} ({variant['type']})\n")
        out.write(f"Position: {variant['pos']}-{variant['end']}\n")
        out.write(f"Reference: {variant['ref']}\n")
        out.write(f"Alternate: {variant['alt']}\n")
        out.write(f"Length Change: {variant['length_change']} bp\n")
        out.write("\n")
    
    # Group features by type
    features_by_type = defaultdict(list)
    for effect in feature_effects:
        features_by_type[effect['feature_type']].append(effect)
    
    # Write feature effect summary
    out.write("\n## Feature Effect Summary\n")
    out.write("-" * 80 + "\n")
    
    # Calculate statistics
    effect_counts = Counter()
    affected_feature_counts = Counter()
    
    for effect in feature_effects:
        if effect['variants']:  # If feature is affected by variants
            affected_feature_counts[effect['feature_type']] += 1
            for e in effect['effects']:
                effect_counts[e] += 1
    
    # Write statistics
    out.write(f"Total Features Analyzed: {len(feature_effects)}\n")
    out.write(f"Features Affected by Variants: {sum(affected_feature_counts.values())}\n\n")
    
    out.write("Affected Features by Type:\n")
    for feature_type, count in sorted(affected_feature_counts.items()):
        total = len(features_by_type[feature_type])
        percentage = (count / total * 100) if total > 0 else 0
        out.write(f"  {feature_type}: {count}/{total} ({percentage:.1f}%)\n")
    
    out.write("\nVariant Effects by Type:\n")
    for effect_type, count in sorted(effect_counts.items(), key=lambda x: x[1], reverse=True):
        out.write(f"  {effect_type}: {count}\n")
    
    # Write detailed feature effects
    out.write("\n\n## Detailed Feature Effects\n")
    
    # Process each feature type
    for feature_type, effects in sorted(features_by_type.items()):
        affected_effects = [e for e in effects if e['variants']]
        
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
            
            # List variants affecting this feature
            out.write("Affected by variants:\n")
            for variant in effect['variants']:
                out.write(f"  - {variant['id']} ({variant['type']})\n")
            
            # List effects
            out.write("Effects:\n")
            for effect_type in sorted(effect['effects']):
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
    
    # Write summary of the most severe effects
    out.write("\n\n## Most Severe Effects Summary\n")
    out.write("-" * 80 + "\n")
    
    # Define severity order
    severity_order = [
        'premature_stop_codon',
        'start_codon_disruption',
        'frame_shift',
        'splice_acceptor_affected',
        'splice_donor_affected',
        'feature_span',
        'feature_start_disruption',
        'feature_end_disruption',
        'amino_acid_change',
        'in_frame_change',
        'substitution',
        'deletion',
        'insertion',
        'inversion',
        'duplication',
        'promoter_affected',
        'terminator_affected',
        'feature_internal',
        'length_change',
        'no_change'
    ]
    
    # Group effects by variant
    effects_by_variant = defaultdict(list)
    
    for effect in feature_effects:
        for variant in effect['variants']:
            effects_by_variant[variant['id']].append(effect)
    
    # Write severe effects by variant
    for var_id, effects in sorted(effects_by_variant.items()):
        # Find most severe effect for each feature
        severe_effects = []
        
        for effect in effects:
            # Find the most severe effect according to severity_order
            most_severe = None
            lowest_severity = float('inf')
            
            for e in effect['effects']:
                try:
                    severity = severity_order.index(e)
                    if severity < lowest_severity:
                        lowest_severity = severity
                        most_severe = e
                except ValueError:
                    # Effect not in severity order, assign low priority
                    pass
            
            if most_severe:
                feature_id = effect['feature']['attributes'].get('ID', 'unknown')
                feature_type = effect['feature_type']
                severe_effects.append((most_severe, feature_id, feature_type))
        
        # Write summary for this variant
        variant = next(v for v in variants if v['id'] == var_id)
        out.write(f"Variant {var_id} ({variant['type']}) - Position: {variant['pos']}-{variant['end']}\n")
        
        # Count effect types
        effect_type_counts = Counter([e[0] for e in severe_effects])
        out.write("Effect types:\n")
        for effect_type, count in sorted(effect_type_counts.items(), key=lambda x: severity_order.index(x[0]) if x[0] in severity_order else 999):
            out.write(f"  - {effect_type}: {count} features\n")
        
        # List affected features by effect type
        if severe_effects:
            out.write("Affected features by severity:\n")
            for effect_type in severity_order:
                features_with_effect = [(f_id, f_type) for e, f_id, f_type in severe_effects if e == effect_type]
                if features_with_effect:
                    out.write(f"  {effect_type}:\n")
                    for f_id, f_type in features_with_effect:
                        out.write(f"    - {f_id} ({f_type})\n")
        
        out.write("\n")
    
    # Close output file if opened
    if outfile:
        out.close()

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

def export_fasta_sequences(ref_seq, alt_seq, prefix):
    """Export reference and alternate sequences to FASTA files."""
    ref_file = f"{prefix}_ref.fa"
    alt_file = f"{prefix}_alt.fa"
    
    # Export reference sequence
    with open(ref_file, 'w') as f:
        f.write(">REF Reference sequence\n")
        # Write in wrapped format (80 chars per line)
        for i in range(0, len(ref_seq), 80):
            f.write(f"{ref_seq[i:i+80]}\n")
    
    # Export alternate sequence
    with open(alt_file, 'w') as f:
        f.write(">ALT Alternate sequence with variants\n")
        # Write in wrapped format (80 chars per line)
        for i in range(0, len(alt_seq), 80):
            f.write(f"{alt_seq[i:i+80]}\n")
    
    logging.info(f"Reference sequence exported to {ref_file}")
    logging.info(f"Alternate sequence exported to {alt_file}")
    
    return ref_file, alt_file

def main():
    """Main function to run the variant effect report generation."""
    parser = argparse.ArgumentParser(description='Generate a variant effect report based on GFA and GFF3 files.')
    parser.add_argument('gfa_file', help='Input GFA file with variants')
    parser.add_argument('gff_file', help='Input GFF3 file with feature annotations')
    parser.add_argument('--output', help='Output report file (default: stdout)')
    
    # Path options
    parser.add_argument('--ref-path', default='REF', help='Name of the reference path in GFA (default: %(default)s)')
    parser.add_argument('--alt-path', default='ALT', help='Name of the alternate path in GFA (default: %(default)s)')
    
    # Report options
    parser.add_argument('--include-sequences', action='store_true', help='Include full sequences in the report')
    parser.add_argument('--simple-report', action='store_true', help='Generate a simplified summary report')
    parser.add_argument('--export-fasta', help='Export reference and alternate sequences to FASTA files with this prefix')
    
    # Output format options
    parser.add_argument('--format', choices=['text', 'json', 'jsonld', 'rdf'], default='text',
                       help='Output format (default: %(default)s)')
    parser.add_argument('--rdf-format', choices=['turtle', 'xml', 'n3', 'nt', 'json-ld'], default='turtle',
                       help='RDF serialization format when using --format=rdf (default: %(default)s)')
    
    # Debug and logging options
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output without full debug')
    parser.add_argument('--log-file', help='Write log to this file')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(debug=args.debug, log_file=args.log_file, verbose=args.verbose)
    
    try:
        # Parse input files
        segments, links, paths, variants = parse_gfa(args.gfa_file)
        features, feature_by_id, children_by_parent = parse_gff3(args.gff_file)
        
        # Mark first CDS in each gene for start codon analysis
        features = identify_first_cds_in_genes(features, feature_by_id, children_by_parent)
        
        # Get reference and alternate paths
        if args.ref_path not in paths:
            available_paths = ", ".join(paths.keys())
            if available_paths:
                logging.error(f"Reference path '{args.ref_path}' not found. Available paths: {available_paths}")
            else:
                logging.error(f"No paths found in GFA file")
            return 1
            
        if args.alt_path not in paths:
            available_paths = ", ".join(paths.keys())
            if available_paths:
                logging.error(f"Alternate path '{args.alt_path}' not found. Available paths: {available_paths}")
            else:
                logging.error(f"No paths found in GFA file")
            return 1
        
        ref_path = paths[args.ref_path]['segments']
        alt_path = paths[args.alt_path]['segments']
        
        # Build reference and alternate sequences
        ref_seq, ref_offsets = build_path_sequence(segments, ref_path)
        alt_seq, alt_offsets = build_path_sequence(segments, alt_path)
        
        logging.info(f"Reference sequence length: {len(ref_seq)} bp")
        logging.info(f"Alternate sequence length: {len(alt_seq)} bp")
        logging.info(f"Length difference: {len(alt_seq) - len(ref_seq)} bp")
        
        # Map features to paths and identify effects
        feature_effects = map_features_to_paths(features, ref_seq, alt_seq, variants)
        
        # Export FASTA sequences if requested
        if args.export_fasta:
            export_fasta_sequences(ref_seq, alt_seq, args.export_fasta)
        
        # Generate output in the requested format
        if args.format == 'text':
            generate_variant_effect_report(feature_effects, variants, args.output)
        elif args.format == 'json':
            generate_json_output(feature_effects, variants, args.output)
        elif args.format == 'jsonld':
            generate_jsonld_output(feature_effects, variants, args.output)
        elif args.format == 'rdf':
            generate_rdf_output(feature_effects, variants, args.output, args.rdf_format)
        else:
            logging.error(f"Unknown output format: {args.format}")
            return 1
        
        logging.info(f"Variant effect report generated successfully in {args.format} format")
        return 0
        
    except Exception as e:
        logging.error(f"Error during report generation: {e}")
        if args.debug:
            import traceback
            logging.error(traceback.format_exc())
        return 1

def generate_json_output(feature_effects, variants, outfile=None):
    """Generate a JSON format report of variant effects on features."""
    import json
    
    # Build the output structure
    output = {
        "meta": {
            "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "tool": "variant-effect-report.py",
            "version": "1.0"
        },
        "variants": [],
        "feature_effects": []
    }
    
    # Add variants
    for variant in variants:
        output["variants"].append({
            "id": variant['id'],
            "type": variant['type'],
            "position": variant['pos'],
            "end": variant['end'],
            "reference": variant['ref'],
            "alternate": variant['alt'],
            "length_change": variant['length_change']
        })
    
    # Add effects
    for effect in feature_effects:
        if not effect['variants']:
            continue
            
        feature = effect['feature']
        feature_id = feature['attributes'].get('ID', 'unknown')
        feature_name = feature['attributes'].get('Name', feature_id)
        
        effect_entry = {
            "feature": {
                "id": feature_id,
                "name": feature_name,
                "type": effect['feature_type'],
                "location": {
                    "start": feature['start'],
                    "end": feature['end'],
                    "strand": feature['strand']
                },
                "attributes": feature['attributes']
            },
            "variants": [v['id'] for v in effect['variants']],
            "effects": effect['effects'],
            "details": effect['details'],
            "sequences": {
                "reference": effect['ref_feature_seq'],
                "alternate": effect['alt_feature_seq']
            }
        }
        
        # Add amino acid changes for CDS features
        if effect['feature_type'] == 'CDS' and 'amino_acid_changes' in effect['details']:
            effect_entry["amino_acid_changes"] = effect['details']['amino_acid_changes']
        
        output["feature_effects"].append(effect_entry)
    
    # Write to file or stdout
    if outfile:
        with open(outfile, 'w') as f:
            json.dump(output, f, indent=2)
        logging.info(f"JSON output written to {outfile}")
    else:
        print(json.dumps(output, indent=2))
    
    return output

def generate_jsonld_output(feature_effects, variants, outfile=None):
    """Generate a JSON-LD format report of variant effects on features."""
    import json
    
    # Build the JSON-LD context
    context = {
        "@context": {
            "@vocab": "http://example.org/vep/",
            "xsd": "http://www.w3.org/2001/XMLSchema#",
            "so": "http://purl.obolibrary.org/obo/SO_",
            "faldo": "http://biohackathon.org/resource/faldo#",
            "variant": {"@id": "hasVariant", "@type": "@id"},
            "feature": {"@id": "affectsFeature", "@type": "@id"},
            "effect": {"@id": "hasEffect", "@type": "@id"},
            "position": {"@id": "faldo:position", "@type": "xsd:integer"},
            "start": {"@id": "faldo:begin", "@type": "xsd:integer"},
            "end": {"@id": "faldo:end", "@type": "xsd:integer"},
            "reference": "hasReferenceSequence",
            "alternate": "hasAlternateSequence",
            "strand": "onStrand",
            "variants": {"@id": "hasVariant", "@type": "@id", "@container": "@set"},
            "effects": {"@id": "hasEffect", "@type": "@id", "@container": "@set"},
            "feature_effects": {"@id": "hasFeatureEffect", "@container": "@set"}
        }
    }
    
    # Map variant types to Sequence Ontology terms
    so_term_map = {
        "SNP": "so:0001483",         # SNV
        "DEL": "so:0000159",         # deletion
        "INS": "so:0000667",         # insertion
        "INV": "so:1000036",         # inversion
        "DUP": "so:1000035"          # duplication
    }
    
    # Map effect types to Sequence Ontology terms
    effect_term_map = {
        "no_change": "so:0001878",               # feature_variant
        "substitution": "so:1000002",            # substitution
        "deletion": "so:0000159",                # deletion
        "insertion": "so:0000667",               # insertion
        "inversion": "so:1000036",               # inversion
        "duplication": "so:1000035",             # duplication
        "frame_shift": "so:0001589",             # frameshift_variant
        "in_frame_change": "so:0001650",         # inframe_variant
        "premature_stop_codon": "so:0001587",    # stop_gained
        "start_codon_disruption": "so:0002012",  # start_lost
        "amino_acid_change": "so:0001992",       # nonsynonymous_variant
        "splice_acceptor_affected": "so:0001574", # splice_acceptor_variant
        "splice_donor_affected": "so:0001575",    # splice_donor_variant
        "feature_span": "so:0001537",            # structural_variant
        "feature_start_disruption": "so:0001580", # coding_sequence_variant
        "feature_end_disruption": "so:0001580",  # coding_sequence_variant
        "feature_internal": "so:0001580",        # coding_sequence_variant
        "promoter_affected": "so:0001631",       # upstream_gene_variant
        "terminator_affected": "so:0001632",     # downstream_gene_variant
        "length_change": "so:0001059"            # sequence_alteration
    }
    
    # Build the output structure
    output = {
        "@context": context["@context"],
        "@id": "urn:x-local:variant-effect-report",
        "@type": "VariantEffectReport",
        "generatedAt": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "variants": [],
        "feature_effects": []
    }
    
    # Add variants
    for variant in variants:
        var_type = variant['type']
        var_entry = {
            "@id": f"variant:{variant['id']}",
            "@type": so_term_map.get(var_type, "SequenceVariant"),
            "id": variant['id'],
            "type": var_type,
            "position": variant['pos'],
            "end": variant['end'],
            "reference": variant['ref'],
            "alternate": variant['alt'],
            "length_change": variant['length_change']
        }
        output["variants"].append(var_entry)
    
    # Add effects
    for effect_index, effect in enumerate(feature_effects):
        if not effect['variants']:
            continue
            
        feature = effect['feature']
        feature_id = feature['attributes'].get('ID', 'unknown')
        feature_name = feature['attributes'].get('Name', feature_id)
        
        # Create feature entry
        feature_entry = {
            "@id": f"feature:{feature_id}",
            "@type": "GenomicFeature",
            "id": feature_id,
            "name": feature_name,
            "featureType": effect['feature_type'],
            "location": {
                "@type": "Region",
                "start": feature['start'],
                "end": feature['end'],
                "strand": feature['strand']
            }
        }
        
        # Create effect entry
        effect_entry = {
            "@id": f"effect:{effect_index}",
            "@type": "VariantEffect",
            "feature": feature_entry["@id"],
            "variants": [f"variant:{v['id']}" for v in effect['variants']],
            "effects": [effect_term_map.get(e, e) for e in effect['effects']],
            "sequences": {
                "reference": effect['ref_feature_seq'],
                "alternate": effect['alt_feature_seq']
            }
        }
        
        # Add details as separate properties
        for k, v in effect['details'].items():
            if k == 'amino_acid_changes' and isinstance(v, dict):
                if v.get('details'):
                    effect_entry["aminoAcidChanges"] = v['details']
            elif isinstance(v, (int, str, bool, float)):
                effect_entry[k] = v
        
        output["feature_effects"].append(effect_entry)
    
    # Write to file or stdout
    if outfile:
        with open(outfile, 'w') as f:
            json.dump(output, f, indent=2)
        logging.info(f"JSON-LD output written to {outfile}")
    else:
        print(json.dumps(output, indent=2))
    
    return output

def generate_rdf_output(feature_effects, variants, outfile=None, format='turtle'):
    """
    Generate RDF output of variant effects on features.
    
    Args:
        feature_effects: List of feature effect dictionaries
        variants: List of variant dictionaries
        outfile: Output file path (default: stdout)
        format: RDF serialization format ('turtle', 'xml', 'n3', 'nt', 'json-ld')
    """
    try:
        from rdflib import Graph, Namespace, URIRef, Literal, BNode
        from rdflib.namespace import RDF, RDFS, XSD
    except ImportError:
        logging.error("RDF output requires rdflib package. Please install with: pip install rdflib")
        return None
    
    # Create the RDF graph
    g = Graph()
    
    # Define namespaces
    VEP = Namespace("http://example.org/vep/")
    SO = Namespace("http://purl.obolibrary.org/obo/SO_")
    FALDO = Namespace("http://biohackathon.org/resource/faldo#")
    
    # Bind namespaces to prefixes
    g.bind("vep", VEP)
    g.bind("so", SO)
    g.bind("faldo", FALDO)
    g.bind("rdf", RDF)
    g.bind("rdfs", RDFS)
    
    # Map variant types to Sequence Ontology terms
    so_term_map = {
        "SNP": SO["0001483"],         # SNV
        "DEL": SO["0000159"],         # deletion
        "INS": SO["0000667"],         # insertion
        "INV": SO["1000036"],         # inversion
        "DUP": SO["1000035"]          # duplication
    }
    
    # Map effect types to Sequence Ontology terms
    effect_term_map = {
        "no_change": SO["0001878"],               # feature_variant
        "substitution": SO["1000002"],            # substitution
        "deletion": SO["0000159"],                # deletion
        "insertion": SO["0000667"],               # insertion
        "inversion": SO["1000036"],               # inversion
        "duplication": SO["1000035"],             # duplication
        "frame_shift": SO["0001589"],             # frameshift_variant
        "in_frame_change": SO["0001650"],         # inframe_variant
        "premature_stop_codon": SO["0001587"],    # stop_gained
        "start_codon_disruption": SO["0002012"],  # start_lost
        "amino_acid_change": SO["0001992"],       # nonsynonymous_variant
        "splice_acceptor_affected": SO["0001574"], # splice_acceptor_variant
        "splice_donor_affected": SO["0001575"],    # splice_donor_variant
        "feature_span": SO["0001537"],            # structural_variant
        "feature_start_disruption": SO["0001580"], # coding_sequence_variant
        "feature_end_disruption": SO["0001580"],  # coding_sequence_variant
        "feature_internal": SO["0001580"],        # coding_sequence_variant
        "promoter_affected": SO["0001631"],       # upstream_gene_variant
        "terminator_affected": SO["0001632"],     # downstream_gene_variant
        "length_change": SO["0001059"]            # sequence_alteration
    }
    
    # Create report node
    report_uri = URIRef("urn:x-local:variant-effect-report")
    g.add((report_uri, RDF.type, VEP.VariantEffectReport))
    g.add((report_uri, VEP.generatedAt, Literal(time.strftime("%Y-%m-%dT%H:%M:%S"), datatype=XSD.dateTime)))
    
    # Add variants
    for variant in variants:
        var_uri = URIRef(f"urn:variant:{variant['id']}")
        var_type = variant['type']
        
        # Add type
        g.add((var_uri, RDF.type, so_term_map.get(var_type, SO["0001060"])))  # Default to sequence_variant
        
        # Add variant properties
        g.add((var_uri, VEP.id, Literal(variant['id'])))
        g.add((var_uri, VEP.type, Literal(var_type)))
        g.add((var_uri, VEP.position, Literal(variant['pos'], datatype=XSD.integer)))
        g.add((var_uri, VEP.end, Literal(variant['end'], datatype=XSD.integer)))
        g.add((var_uri, VEP.reference, Literal(variant['ref'])))
        g.add((var_uri, VEP.alternate, Literal(variant['alt'])))
        g.add((var_uri, VEP.lengthChange, Literal(variant['length_change'], datatype=XSD.integer)))
        
        # Link to report
        g.add((report_uri, VEP.hasVariant, var_uri))
    
    # Add effects
    for effect_index, effect in enumerate(feature_effects):
        if not effect['variants']:
            continue
            
        feature = effect['feature']
        feature_id = feature['attributes'].get('ID', 'unknown')
        feature_name = feature['attributes'].get('Name', feature_id)
        
        # Create feature node
        feature_uri = URIRef(f"urn:feature:{feature_id}")
        g.add((feature_uri, RDF.type, VEP.GenomicFeature))
        g.add((feature_uri, VEP.id, Literal(feature_id)))
        g.add((feature_uri, VEP.name, Literal(feature_name)))
        g.add((feature_uri, VEP.featureType, Literal(effect['feature_type'])))
        
        # Create location node
        location_node = BNode()
        g.add((location_node, RDF.type, FALDO.Region))
        g.add((location_node, FALDO.begin, Literal(feature['start'], datatype=XSD.integer)))
        g.add((location_node, FALDO.end, Literal(feature['end'], datatype=XSD.integer)))
        g.add((location_node, FALDO.strand, Literal(feature['strand'])))
        g.add((feature_uri, FALDO.location, location_node))
        
        # Create effect node
        effect_uri = URIRef(f"urn:effect:{effect_index}")
        g.add((effect_uri, RDF.type, VEP.VariantEffect))
        g.add((effect_uri, VEP.affectsFeature, feature_uri))
        
        # Link to variants
        for variant in effect['variants']:
            var_uri = URIRef(f"urn:variant:{variant['id']}")
            g.add((effect_uri, VEP.hasVariant, var_uri))
        
        # Add effects
        for e in effect['effects']:
            effect_term = effect_term_map.get(e, VEP[e])
            g.add((effect_uri, VEP.hasEffect, effect_term))
        
        # Add sequences
        g.add((effect_uri, VEP.referenceSequence, Literal(effect['ref_feature_seq'])))
        g.add((effect_uri, VEP.alternateSequence, Literal(effect['alt_feature_seq'])))
        
        # Add details
        for k, v in effect['details'].items():
            if k == 'amino_acid_changes' and isinstance(v, dict):
                if v.get('details'):
                    for change in v['details']:
                        g.add((effect_uri, VEP.aminoAcidChange, Literal(change)))
            elif isinstance(v, (int, str, bool)):
                g.add((effect_uri, VEP[k], Literal(v)))
        
        # Link to report
        g.add((report_uri, VEP.hasFeatureEffect, effect_uri))
    
    # Serialize the graph
    if outfile:
        g.serialize(destination=outfile, format=format)
        logging.info(f"RDF output ({format}) written to {outfile}")
    else:
        print(g.serialize(format=format).decode('utf-8'))
    
    return g
    
if __name__ == "__main__":
    sys.exit(main())
