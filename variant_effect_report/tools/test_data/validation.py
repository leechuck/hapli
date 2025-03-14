"""
Biological consistency validation for test data.
"""

import logging
from typing import Dict, List, Set, Tuple, Optional
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


class BiologicalValidator:
    """Validate biological consistency of variants"""
    
    @staticmethod
    def check_variant_overlaps(variants: List[Dict]) -> List[Tuple[int, int]]:
        """
        Check for overlapping variants
        
        Args:
            variants: List of variant dictionaries
            
        Returns:
            List of (index1, index2) tuples for overlapping variants
        """
        overlaps = []
        
        # Sort variants by position
        sorted_variants = sorted(enumerate(variants), key=lambda x: x[1]['pos'])
        
        # Check for overlaps
        for i in range(len(sorted_variants) - 1):
            idx1, var1 = sorted_variants[i]
            pos1 = var1['pos']
            end1 = pos1 + len(var1['ref']) - 1
            
            for j in range(i + 1, len(sorted_variants)):
                idx2, var2 = sorted_variants[j]
                pos2 = var2['pos']
                
                # If the next variant starts after the current one ends, no overlap
                if pos2 > end1:
                    break
                    
                # Otherwise, we have an overlap
                overlaps.append((idx1, idx2))
        
        return overlaps
    
    @staticmethod
    def validate_variants(variants: List[Dict], reference_sequence: Optional[str] = None) -> Tuple[bool, List[str]]:
        """
        Validate a list of variants for biological consistency
        
        Args:
            variants: List of variant dictionaries
            reference_sequence: Optional reference sequence to validate against
            
        Returns:
            Tuple of (is_valid, error_messages)
        """
        errors = []
        
        # Check for overlapping variants
        overlaps = BiologicalValidator.check_variant_overlaps(variants)
        if overlaps:
            for idx1, idx2 in overlaps:
                var1 = variants[idx1]
                var2 = variants[idx2]
                errors.append(
                    f"Overlapping variants: {var1['type']} at position {var1['pos']} "
                    f"overlaps with {var2['type']} at position {var2['pos']}"
                )
        
        # Check for valid reference bases
        for i, var in enumerate(variants):
            ref = var['ref']
            if not all(base in 'ACGTNacgtn' for base in ref):
                errors.append(f"Invalid reference bases in variant {i+1}: {ref}")
        
        # Check for valid alternate bases
        for i, var in enumerate(variants):
            alt = var['alt']
            if not all(base in 'ACGTNacgtn' for base in alt):
                errors.append(f"Invalid alternate bases in variant {i+1}: {alt}")
        
        # If reference sequence is provided, validate reference bases
        if reference_sequence:
            for i, var in enumerate(variants):
                pos = var['pos'] - 1  # Convert to 0-based
                ref = var['ref']
                
                # Check if position is within reference
                if pos >= len(reference_sequence):
                    errors.append(f"Variant {i+1} position {var['pos']} exceeds reference length {len(reference_sequence)}")
                    continue
                
                # Check if reference matches
                if pos + len(ref) <= len(reference_sequence):
                    ref_in_seq = reference_sequence[pos:pos+len(ref)]
                    if ref != ref_in_seq:
                        errors.append(
                            f"Reference mismatch for variant {i+1} at position {var['pos']}: "
                            f"variant has '{ref}', reference has '{ref_in_seq}'"
                        )
            
            # Check GC content of insertions (should be similar to reference)
            ref_gc = gc_fraction(reference_sequence) * 100  # Convert to percentage
            for i, var in enumerate(variants):
                if var['type'] == 'INS':
                    ins_seq = var['alt'][1:]  # Remove first base which is reference
                    if len(ins_seq) >= 10:  # Only check longer insertions
                        ins_gc = gc_fraction(ins_seq) * 100  # Convert to percentage
                        if abs(ins_gc - ref_gc) > 20:  # More than 20% difference
                            errors.append(
                                f"Unusual GC content in insertion at position {var['pos']}: "
                                f"insertion has {ins_gc:.1f}%, reference has {ref_gc:.1f}%"
                            )
        
        return len(errors) == 0, errors
