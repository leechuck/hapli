"""
Utility functions for variant effect analysis.
"""

import logging
import sys
from typing import Dict, List, Optional, Tuple
from Bio.Seq import Seq


def setup_logging(debug: bool = False, log_file: Optional[str] = None, verbose: bool = False) -> logging.Logger:
    """
    Configure logging based on debug flag and optional log file.
    
    Args:
        debug: Enable debug output
        log_file: Write log to this file
        verbose: Enable verbose output without full debug
        
    Returns:
        Configured logger
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


def reverse_complement(sequence: str) -> str:
    """
    Get the reverse complement of a DNA sequence.
    
    Args:
        sequence: DNA sequence
        
    Returns:
        Reverse complement of the sequence
    """
    seq_obj = Seq(sequence)
    return str(seq_obj.reverse_complement())


def translate_sequence(nucleotide_seq: str, strand: str = '+') -> str:
    """
    Translate a nucleotide sequence to amino acids, considering strand.
    
    Args:
        nucleotide_seq: Nucleotide sequence
        strand: Strand ('+' or '-')
        
    Returns:
        Amino acid sequence
    """
    if not nucleotide_seq:
        return ""
    
    # Handle reverse strand
    if strand == '-':
        nucleotide_seq = reverse_complement(nucleotide_seq)
    
    # Ensure sequence length is a multiple of 3
    remainder = len(nucleotide_seq) % 3
    if remainder > 0:
        nucleotide_seq = nucleotide_seq[:-remainder]
    
    if not nucleotide_seq:
        return ""
    
    # Translate to amino acids
    seq_obj = Seq(nucleotide_seq)
    return str(seq_obj.translate())


def compare_amino_acid_sequences(ref_aa: str, alt_aa: str) -> Dict:
    """
    Compare two amino acid sequences and identify changes.
    
    Args:
        ref_aa: Reference amino acid sequence
        alt_aa: Alternate amino acid sequence
        
    Returns:
        Dictionary with comparison results
    """
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


def build_path_sequence(segments: Dict, path_segments: List[Tuple[str, str]]) -> Tuple[str, List[Dict]]:
    """
    Build a sequence from path segments.
    
    Args:
        segments: Dictionary of segment data
        path_segments: List of (segment_id, orientation) tuples
        
    Returns:
        Tuple of (sequence, segment_offsets)
    """
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
            segment_seq = reverse_complement(segment_seq)
            
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
