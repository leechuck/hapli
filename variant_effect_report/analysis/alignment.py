"""
Sequence alignment functionality for variant effect analysis.
"""

import re
from typing import Dict, List, Tuple
from Bio import Align
from Bio.Seq import Seq


class SequenceAligner:
    """Perform sequence alignment and analyze differences."""
    
    def __init__(self, gap_open=-10, gap_extend=-0.5, match=2, mismatch=-1):
        """
        Initialize the sequence aligner.
        
        Args:
            gap_open: Gap opening penalty
            gap_extend: Gap extension penalty
            match: Match score
            mismatch: Mismatch score
        """
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.match = match
        self.mismatch = mismatch
    
    def align_sequences(self, ref_seq: str, alt_seq: str) -> Dict:
        """
        Perform global sequence alignment to identify differences between reference and alternate sequences.
        
        Args:
            ref_seq: Reference sequence string
            alt_seq: Alternate sequence string
            
        Returns:
            Dictionary with alignment information and identified changes
        """
        # Create aligner
        aligner = Align.PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = self.match
        aligner.mismatch_score = self.mismatch
        aligner.open_gap_score = self.gap_open
        aligner.extend_gap_score = self.gap_extend
        
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
        
        # The alignment output format is:
        # Line 0: reference sequence with gaps
        # Line 1: alignment characters (| for match, space for mismatch)
        # Line 2: query sequence with gaps
        if len(alignment_strings) >= 3:
            ref_aligned = alignment_strings[0]
            alt_aligned = alignment_strings[2]
        else:
            # Fallback if alignment string format is different
            ref_aligned = ref_seq
            alt_aligned = alt_seq
        
        # Identify changes from the alignment
        changes = self.analyze_alignment_changes(ref_aligned, alt_aligned)
        
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
    
    def analyze_alignment_changes(self, ref_aligned: str, alt_aligned: str) -> Dict:
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
        inversions = self.detect_inversions(ref_aligned, alt_aligned)
        
        # Identify complex regions (potential rearrangements)
        complex_regions = self.identify_complex_regions(ref_aligned, alt_aligned)
        
        return {
            'insertions': insertions,
            'deletions': deletions,
            'substitutions': substitutions,
            'inversions': inversions,
            'complex_regions': complex_regions
        }
    
    def detect_inversions(self, ref_aligned: str, alt_aligned: str) -> List[Dict]:
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
    
    def identify_complex_regions(self, ref_aligned: str, alt_aligned: str) -> List[Dict]:
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
