"""
Variant analysis functions for identifying and characterizing variants.
"""

import logging
import re
from Bio.Seq import Seq

def calculate_variant_length_changes(variants, segments):
    """
    Calculate length changes for each variant by comparing variant segments to their original segments.
    
    Args:
        variants: List of variant dictionaries
        segments: Dictionary of segment data
        
    Returns:
        Updated list of variants with accurate length change values
    """
    for variant in variants:
        # Skip if length change is already set and non-zero
        if variant.get('length_change', 0) != 0:
            continue
            
        var_segments = variant.get('segments', [])
        if not var_segments:
            continue
            
        total_length_change = 0
        
        for seg_id in var_segments:
            if seg_id not in segments:
                continue
                
            # Get variant segment length
            var_seg_length = segments[seg_id].get('length', 0)
            
            # Try to get original segment ID
            orig_seg_id = segments[seg_id].get('original_segment')
            
            # If not directly specified, try to extract from segment name
            if not orig_seg_id and '_VAR' in seg_id:
                parts = seg_id.split('_VAR')
                if len(parts) > 1:
                    orig_seg_id = parts[0]
            
            # Calculate length change if we found the original segment
            if orig_seg_id and orig_seg_id in segments:
                orig_seg_length = segments[orig_seg_id].get('length', 0)
                segment_length_change = var_seg_length - orig_seg_length
                total_length_change += segment_length_change
            # Get length change from the segment's tag if available
            elif 'LD' in segments[seg_id].get('tags', {}):
                tag_data = segments[seg_id]['tags']['LD']
                if tag_data[0] == 'i':  # Integer type
                    try:
                        segment_length_change = int(tag_data[1])
                        total_length_change += segment_length_change
                    except ValueError:
                        pass
        
        # Update the variant with the calculated length change
        variant['length_change'] = total_length_change
    
    return variants

def infer_variant_types(variants, segments):
    """
    Infer variant types based on segment information and length changes.
    
    Args:
        variants: List of variant dictionaries
        segments: Dictionary of segment data
        
    Returns:
        Updated list of variants with inferred types
    """
    for variant in variants:
        # Skip if type is already defined and not UNKNOWN
        if variant.get('type') and variant.get('type') != 'UNKNOWN':
            continue
            
        length_change = variant.get('length_change', 0)
        var_segments = variant.get('segments', [])
        
        # Look for variant type hints in segment names
        for seg_id in var_segments:
            if seg_id not in segments:
                continue
                
            # Check segment name for type hints
            if "SNP" in seg_id or "SNV" in seg_id:
                variant['type'] = "SNP"
                break
            elif "DEL" in seg_id:
                variant['type'] = "DEL"
                break
            elif "INS" in seg_id:
                variant['type'] = "INS"
                break
            elif "DUP" in seg_id:
                variant['type'] = "DUP"
                break
            elif "INV" in seg_id:
                variant['type'] = "INV"
                break
        
        # If no type found from names, infer from length change
        if (not variant.get('type') or variant.get('type') == 'UNKNOWN') and length_change != 0:
            if length_change > 0:
                # Check for duplication - if the variant segment is much longer
                for seg_id in var_segments:
                    if seg_id in segments:
                        seg_data = segments[seg_id]
                        orig_seg_id = seg_data.get('original_segment')
                        
                        if orig_seg_id and orig_seg_id in segments:
                            var_seq = seg_data.get('sequence', '')
                            orig_seq = segments[orig_seg_id].get('sequence', '')
                            
                            # If the variant contains the original sequence multiple times
                            if orig_seq and len(orig_seq) > 0 and var_seq.count(orig_seq) > 1:
                                variant['type'] = "DUP"
                                break
                
                # Default to insertion if not a duplication
                if not variant.get('type') or variant.get('type') == 'UNKNOWN':
                    variant['type'] = "INS"
            elif length_change < 0:
                variant['type'] = "DEL"
        
        # Default to SNP for small zero-length-change variants without a type
        if (not variant.get('type') or variant.get('type') == 'UNKNOWN') and length_change == 0:
            variant['type'] = "SNP"
    
    return variants

def identify_repeat_sequences(variant, segments):
    """Identify repeat patterns in insertion or duplication variants."""
    if variant['type'] not in ['DUP', 'INS']:
        return None
    
    # Try to find the sequence for this variant
    variant_sequence = ""
    for seg_id in variant.get('segments', []):
        if seg_id in segments:
            variant_sequence = segments[seg_id]['sequence']
            break
    
    if not variant_sequence:
        return None
    
    # Look for repeating patterns
    repeat_info = {}
    
    # Try different pattern lengths, starting from shorter patterns
    for pattern_len in range(1, min(len(variant_sequence) // 2 + 1, 20)):
        pattern = variant_sequence[:pattern_len]
        
        # Count how many times the pattern repeats
        count = 0
        pos = 0
        while pos + pattern_len <= len(variant_sequence) and variant_sequence[pos:pos+pattern_len] == pattern:
            count += 1
            pos += pattern_len
        
        # If we found a significant repeat (at least 3 times)
        if count >= 3 and count * pattern_len >= len(variant_sequence) * 0.6:
            repeat_info['sequence'] = pattern
            repeat_info['count'] = count
            break
    
    return repeat_info or None
