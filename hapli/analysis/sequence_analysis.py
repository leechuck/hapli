"""
Sequence analysis functions for comparing reference and alternate sequences.
"""

import logging
import re
from Bio import Seq, Align
from Bio.Seq import Seq

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

def reverse_complement(sequence):
    """Return the reverse complement of a DNA sequence."""
    seq_obj = Seq(sequence)
    return str(seq_obj.reverse_complement())

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

def identify_variants_by_position(segment_offsets, ref_segments, variant_segments, segments=None):
    """Identify variants in an alternate sequence based on segment differences."""
    variants_by_position = []
    
    # If no reference segments provided, return empty list
    if not ref_segments:
        return variants_by_position
        
    # Create a map of reference positions
    ref_pos_map = {}
    current_pos = 0
    
    for seg_info in ref_segments:
        seg_id = seg_info.get('seg_id')
        length = seg_info.get('length', 0)
        ref_pos_map[seg_id] = (current_pos, current_pos + length - 1)
        current_pos += length
    
    # Find alternate segments that differ from reference
    for seg_info in segment_offsets:
        seg_id = seg_info.get('seg_id')
        start_pos = seg_info.get('start', 0)
        length = seg_info.get('length', 0)
        variant_id = seg_info.get('variant_id')
        
        # Skip segments without variant IDs
        if not variant_id:
            variant_id = segments.get(seg_id, {}).get('variant_id') if segments else None
            if not variant_id:
                continue
        
        # Get variant information
        var_data = variant_segments.get(seg_id, {})
        if not var_data and segments:
            var_data = segments.get(seg_id, {})
        
        # Try to find the original segment
        ref_id = var_data.get('original_segment')
        if not ref_id and '_VAR' in seg_id:
            # Attempt to extract original segment from ID
            parts = seg_id.split('_VAR')
            if len(parts) > 1:
                ref_id = parts[0]
        
        # Find reference position for this variant
        if ref_id and ref_id in ref_pos_map:
            ref_start, ref_end = ref_pos_map[ref_id]
            # Create a variant entry
            variant = {
                'id': variant_id,
                'type': var_data.get('variant_type', 'UNKNOWN'),
                'pos': ref_start + 1,  # Convert to 1-based
                'end': ref_end + 1,    # Convert to 1-based
                'length_change': var_data.get('length_change', 0),
                'segments': [seg_id],
                'original_segment': ref_id
            }
            variants_by_position.append(variant)
    
    return variants_by_position

def analyze_sequence_alignment(ref_seq, alt_seq, gap_open=-10, gap_extend=-0.5, match=2, mismatch=-1):
    """
    Perform global sequence alignment to identify differences between reference and alternate sequences.
    
    Uses Bio.Align.PairwiseAligner instead of the deprecated Bio.pairwise2
    
    Args:
        ref_seq: Reference sequence string
        alt_seq: Alternate sequence string
        gap_open: Gap opening penalty (default: -10)
        gap_extend: Gap extension penalty (default: -0.5)
        match: Match score (default: 2)
        mismatch: Mismatch score (default: -1)
        
    Returns:
        Dictionary with alignment information and identified changes
    """
    # Create aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    
    # Perform alignment
    alignments = aligner.align(ref_seq, alt_seq)
    
    # Get the best alignment
    if not alignments:
        return {
            'alignment_score': 0,
            'changes': 'complete_replacement',
            'ref_aligned': ref_seq,
            'alt_aligned': alt_seq,
            'insertions': [],
            'deletions': [],
            'substitutions': [],
            'inversions': [],
            'length_change': len(alt_seq) - len(ref_seq)
        }
    
    best_alignment = alignments[0]
    score = best_alignment.score
    
    # Convert alignment to string representation
    alignment_strings = str(best_alignment).split('\n')
    
    # Extract aligned sequences
    ref_aligned = ""
    alt_aligned = ""
    
    # Process the alignment path to get aligned sequences with gaps
    ref_idx = 0
    alt_idx = 0
    
    for op in best_alignment.path:
        if op[0] == op[1]:  # Match or mismatch
            ref_aligned += ref_seq[ref_idx]
            alt_aligned += alt_seq[alt_idx]
            ref_idx += 1
            alt_idx += 1
        elif op[0] > op[1]:  # Gap in alt sequence
            ref_aligned += ref_seq[ref_idx]
            alt_aligned += "-"
            ref_idx += 1
        else:  # Gap in ref sequence
            ref_aligned += "-"
            alt_aligned += alt_seq[alt_idx]
            alt_idx += 1
    
    # Fallback to the string representation if the above method fails
    if not ref_aligned or not alt_aligned:
        if len(alignment_strings) >= 3:
            ref_aligned = alignment_strings[0]
            alt_aligned = alignment_strings[2]
        else:
            # Last resort fallback
            ref_aligned = ref_seq
            alt_aligned = alt_seq
    
    # Identify changes from the alignment
    changes = analyze_alignment_changes(ref_aligned, alt_aligned)
    
    return {
        'alignment_score': score,
        'ref_aligned': ref_aligned,
        'alt_aligned': alt_aligned,
        'insertions': changes['insertions'],
        'deletions': changes['deletions'],
        'substitutions': changes['substitutions'],
        'inversions': changes['inversions'],
        'complex_regions': changes['complex_regions'],
        'length_change': len(alt_seq) - len(ref_seq)
    }

def analyze_alignment_changes(ref_aligned, alt_aligned):
    """
    Analyze an alignment to identify mutations, insertions, deletions, and potential inversions.
    
    Args:
        ref_aligned: Reference sequence with alignment gaps
        alt_aligned: Alternate sequence with alignment gaps
        
    Returns:
        Dictionary with lists of detected changes
    """
    insertions = []
    deletions = []
    substitutions = []
    complex_regions = []
    
    # Find simple insertions, deletions, and substitutions
    current_ins = None
    current_del = None
    
    for i in range(len(ref_aligned)):
        ref_base = ref_aligned[i]
        alt_base = alt_aligned[i]
        
        if ref_base == '-' and alt_base != '-':
            # Insertion in alternate sequence
            if current_ins is None:
                current_ins = {'start': i, 'sequence': alt_base}
            else:
                current_ins['sequence'] += alt_base
        elif ref_base != '-' and alt_base == '-':
            # Deletion in alternate sequence
            if current_del is None:
                current_del = {'start': i, 'sequence': ref_base}
            else:
                current_del['sequence'] += ref_base
        elif ref_base != alt_base:
            # Substitution
            substitutions.append({'position': i, 'ref': ref_base, 'alt': alt_base})
            
            # Close any open indels
            if current_ins:
                insertions.append(current_ins)
                current_ins = None
            if current_del:
                deletions.append(current_del)
                current_del = None
        else:
            # Matching position
            # Close any open indels
            if current_ins:
                insertions.append(current_ins)
                current_ins = None
            if current_del:
                deletions.append(current_del)
                current_del = None
    
    # Close any open indels at the end
    if current_ins:
        insertions.append(current_ins)
    if current_del:
        deletions.append(current_del)
    
    # Look for potential inversions
    inversions = detect_inversions(ref_aligned, alt_aligned)
    
    # Identify complex regions (potential rearrangements)
    complex_regions = identify_complex_regions(ref_aligned, alt_aligned)
    
    return {
        'insertions': insertions,
        'deletions': deletions,
        'substitutions': substitutions,
        'inversions': inversions,
        'complex_regions': complex_regions
    }

def detect_inversions(ref_aligned, alt_aligned):
    """
    Detect potential inversions by looking for regions where the alternate sequence
    matches the reverse complement of the reference.
    
    Args:
        ref_aligned: Reference sequence with alignment gaps
        alt_aligned: Alternate sequence with alignment gaps
        
    Returns:
        List of detected inversions
    """
    inversions = []
    
    # Strip gaps for sequence comparison
    ref_seq = ref_aligned.replace('-', '')
    alt_seq = alt_aligned.replace('-', '')
    
    # Minimum inversion size to consider (to avoid random matches)
    min_inversion_size = 10
    
    # Look for potential inversions of various sizes
    for i in range(len(ref_seq) - min_inversion_size):
        for j in range(i + min_inversion_size, min(i + 200, len(ref_seq) + 1)):
            ref_segment = ref_seq[i:j]
            ref_revcomp = str(Seq(ref_segment).reverse_complement())
            
            # Search for reverse complement in alternate sequence
            for match in re.finditer(re.escape(ref_revcomp), alt_seq):
                start_alt = match.start()
                end_alt = match.end()
                
                inversions.append({
                    'ref_start': i,
                    'ref_end': j,
                    'alt_start': start_alt,
                    'alt_end': end_alt,
                    'length': j - i,
                    'ref_segment': ref_segment,
                    'alt_segment': ref_revcomp
                })
    
    # Filter overlapping inversions to keep only the most significant ones
    if inversions:
        inversions.sort(key=lambda x: x['length'], reverse=True)
        filtered_inversions = [inversions[0]]
        
        for inv in inversions[1:]:
            # Check if this inversion overlaps with any already selected
            overlaps = False
            for selected in filtered_inversions:
                # Check ref overlap
                ref_overlap = not (inv['ref_end'] <= selected['ref_start'] or inv['ref_start'] >= selected['ref_end'])
                # Check alt overlap
                alt_overlap = not (inv['alt_end'] <= selected['alt_start'] or inv['alt_start'] >= selected['alt_end'])
                
                if ref_overlap or alt_overlap:
                    overlaps = True
                    break
            
            if not overlaps:
                filtered_inversions.append(inv)
        
        return filtered_inversions
    
    return inversions

def identify_complex_regions(ref_aligned, alt_aligned):
    """
    Identify complex regions that might represent rearrangements or complex changes.
    
    Args:
        ref_aligned: Reference sequence with alignment gaps
        alt_aligned: Alternate sequence with alignment gaps
        
    Returns:
        List of complex regions
    """
    complex_regions = []
    
    # Initialize variables to track complex regions
    in_complex_region = False
    current_region = None
    mismatch_count = 0
    window_size = 20
    complexity_threshold = 8  # Number of mismatches in window_size to consider complex
    
    for i in range(len(ref_aligned)):
        # Check for region with high density of mismatches
        start_window = max(0, i - window_size)
        mismatches_in_window = sum(1 for j in range(start_window, i + 1) 
                                   if j < len(ref_aligned) and ref_aligned[j] != alt_aligned[j])
        
        if mismatches_in_window >= complexity_threshold:
            if not in_complex_region:
                in_complex_region = True
                current_region = {
                    'start': max(0, i - window_size),
                    'ref_sequence': ref_aligned[max(0, i - window_size):i + 1],
                    'alt_sequence': alt_aligned[max(0, i - window_size):i + 1]
                }
            else:
                # Extend current complex region
                current_region['ref_sequence'] += ref_aligned[i:i+1]
                current_region['alt_sequence'] += alt_aligned[i:i+1]
        else:
            if in_complex_region:
                in_complex_region = False
                current_region['end'] = i
                # Remove gaps for a cleaner representation
                current_region['ref_sequence_nogap'] = current_region['ref_sequence'].replace('-', '')
                current_region['alt_sequence_nogap'] = current_region['alt_sequence'].replace('-', '')
                complex_regions.append(current_region)
                current_region = None
    
    # Close any open complex region at the end
    if in_complex_region:
        current_region['end'] = len(ref_aligned)
        current_region['ref_sequence_nogap'] = current_region['ref_sequence'].replace('-', '')
        current_region['alt_sequence_nogap'] = current_region['alt_sequence'].replace('-', '')
        complex_regions.append(current_region)
    
    return complex_regions

def analyze_feature_with_alignment(feature, ref_seq, alt_seq):
    """
    Analyze a feature using sequence alignment to determine the effects of variants.
    
    Args:
        feature: Feature dictionary
        ref_seq: Reference sequence
        alt_seq: Alternate sequence
        
    Returns:
        Dictionary with effect analysis results
    """
    # Extract feature boundaries
    start = feature['start'] - 1  # Convert to 0-based
    end = feature['end'] - 1      # Convert to 0-based
    feature_type = feature['type']
    strand = feature['strand']
    
    # Skip if feature is outside sequence bounds
    if start >= len(ref_seq) or end >= len(ref_seq):
        logging.warning(f"Feature {feature['attributes'].get('ID', 'unknown')} is outside sequence bounds, skipping")
        return {
            'feature': feature,
            'feature_type': feature_type,
            'effects': ['no_change'],
            'details': {}
        }
    
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
    
    # If the sequences are identical, no change
    if ref_feature_seq == alt_feature_seq:
        return {
            'feature': feature,
            'feature_type': feature_type,
            'ref_feature_seq': ref_feature_seq,
            'alt_feature_seq': alt_feature_seq,
            'effects': ['no_change'],
            'details': {}
        }
    
    # Perform sequence alignment for detailed analysis
    alignment_result = analyze_sequence_alignment(ref_feature_seq, alt_feature_seq)
    
    # Determine effects based on alignment
    effects = []
    effect_details = {}
    
    # Calculate length change
    length_change = len(alt_feature_seq) - len(ref_feature_seq)
    if length_change != 0:
        effects.append('length_change')
        effect_details['length_change'] = length_change
    
    # Check for inversions
    if alignment_result['inversions']:
        effects.append('inversion')
        effect_details['inversions'] = alignment_result['inversions']
    
    # Check for complex regions (potential rearrangements)
    if alignment_result['complex_regions']:
        effects.append('complex_rearrangement')
        effect_details['complex_regions'] = alignment_result['complex_regions']
    
    # Check for insertions
    if alignment_result['insertions']:
        effects.append('insertion')
        effect_details['insertions'] = alignment_result['insertions']
    
    # Check for deletions
    if alignment_result['deletions']:
        effects.append('deletion')
        effect_details['deletions'] = alignment_result['deletions']
    
    # Check for substitutions
    if alignment_result['substitutions']:
        effects.append('substitution')
        effect_details['substitutions'] = alignment_result['substitutions']
    
    # Add alignment information
    effect_details['alignment'] = {
        'score': alignment_result['alignment_score'],
        'ref_aligned': alignment_result['ref_aligned'],
        'alt_aligned': alignment_result['alt_aligned']
    }
    
    # Check for more specific effects based on feature type
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
    
    return {
        'feature': feature,
        'feature_type': feature_type,
        'ref_feature_seq': ref_feature_seq,
        'alt_feature_seq': alt_feature_seq,
        'effects': effects,
        'details': effect_details
    }

def analyze_haplotype_differences_with_alignment(ref_seq, alt_seq, features, segment_offsets=None, variant_segments=None, segments=None):
    """
    Analyze differences between reference and alternate haplotype sequences using sequence alignment.
    
    Args:
        ref_seq: Reference sequence string
        alt_seq: Alternate sequence string
        features: List of feature dictionaries
        segment_offsets: Optional segment offset information
        variant_segments: Optional variant segment information
        segments: Optional segment data
        
    Returns:
        List of feature effect dictionaries
    """
    feature_effects = []
    
    # Find variants from segment information if available
    variants = []
    if segment_offsets and variant_segments and segments:
        # Find segments in the alternate sequence that correspond to variants
        ref_segments = []
        for seg_id in sorted(segments.keys()):
            seg_data = segments[seg_id]
            if not seg_data.get('variant_id') and not '_VAR' in seg_id:
                ref_segments.append({
                    'seg_id': seg_id,
                    'orientation': '+',
                    'length': seg_data['length']
                })
        
        variants = identify_variants_by_position(segment_offsets, ref_segments, variant_segments, segments)
    
    # Analyze each feature using alignment
    for feature in features:
        effect = analyze_feature_with_alignment(feature, ref_seq, alt_seq)
        
        # Add variants if available
        if variants:
            # Find variants that overlap with this feature
            feature_start = feature['start'] - 1  # Convert to 0-based
            feature_end = feature['end'] - 1      # Convert to 0-based
            
            overlapping_variants = []
            for variant in variants:
                var_start = variant.get('pos', 0) - 1  # Convert to 0-based
                var_end = variant.get('end', var_start)
                
                # If positions are unknown (0), check if the variant's segments are in this path
                if var_start == -1 and var_end == -1:
                    # Include the variant if it appears in our segment offsets
                    var_segments = variant.get('segments', [])
                    if segment_offsets and any(seg['seg_id'] in var_segments for seg in segment_offsets):
                        overlapping_variants.append(variant)
                    continue
                
                # Check for overlap with normal positions
                if not (var_end < feature_start or var_start > feature_end):
                    overlapping_variants.append(variant)
            
            effect['variants'] = overlapping_variants
        
        feature_effects.append(effect)
    
    return feature_effects

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
