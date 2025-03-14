"""
Feature effect analysis for variant effect reporting.
"""

import logging
from typing import Dict, List, Optional, Tuple

from variant_effect_report.core.models import Feature, Variant, FeatureEffect
from variant_effect_report.core.utils import translate_sequence, compare_amino_acid_sequences
from variant_effect_report.analysis.alignment import SequenceAligner


class FeatureEffectAnalyzer:
    """Analyze effects of variants on genomic features."""
    
    def __init__(self, use_alignment: bool = False):
        """
        Initialize the feature effect analyzer.
        
        Args:
            use_alignment: Whether to use sequence alignment for analysis
        """
        self.use_alignment = use_alignment
        self.sequence_aligner = SequenceAligner() if use_alignment else None
    
    def analyze_feature(self, feature: Feature, ref_seq: str, alt_seq: str, 
                       overlapping_variants: Optional[List[Variant]] = None) -> FeatureEffect:
        """
        Analyze a feature to determine the effects of variants.
        
        Args:
            feature: Feature to analyze
            ref_seq: Reference sequence
            alt_seq: Alternate sequence
            overlapping_variants: List of variants overlapping with this feature
            
        Returns:
            FeatureEffect object with analysis results
        """
        # Extract feature boundaries
        start = feature.start - 1  # Convert to 0-based
        end = feature.end - 1      # Convert to 0-based
        feature_type = feature.type
        strand = feature.strand
        
        # Skip if feature is outside sequence bounds
        if start >= len(ref_seq) or end >= len(ref_seq):
            logging.warning(f"Feature {feature.attributes.get('ID', 'unknown')} is outside sequence bounds, skipping")
            return FeatureEffect(
                feature=feature,
                feature_type=feature_type,
                effects=['no_change'],
                details={}
            )
        
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
        if ref_feature_seq == alt_feature_seq and not overlapping_variants:
            return FeatureEffect(
                feature=feature,
                feature_type=feature_type,
                ref_feature_seq=ref_feature_seq,
                alt_feature_seq=alt_feature_seq,
                effects=['no_change'],
                details={},
                variants=overlapping_variants or []
            )
        
        # Determine effects based on sequence differences
        if self.use_alignment and self.sequence_aligner:
            return self._analyze_with_alignment(feature, ref_feature_seq, alt_feature_seq, overlapping_variants)
        else:
            return self._analyze_basic_effects(feature, ref_feature_seq, alt_feature_seq, overlapping_variants)
    
    def _analyze_with_alignment(self, feature: Feature, ref_feature_seq: str, alt_feature_seq: str,
                              overlapping_variants: Optional[List[Variant]] = None) -> FeatureEffect:
        """
        Analyze a feature using sequence alignment.
        
        Args:
            feature: Feature to analyze
            ref_feature_seq: Reference feature sequence
            alt_feature_seq: Alternate feature sequence
            overlapping_variants: List of variants overlapping with this feature
            
        Returns:
            FeatureEffect object with analysis results
        """
        # Perform sequence alignment for detailed analysis
        alignment_result = self.sequence_aligner.align_sequences(ref_feature_seq, alt_feature_seq)
        
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
        self._analyze_feature_type_specific_effects(feature, ref_feature_seq, alt_feature_seq, 
                                                  length_change, effects, effect_details)
        
        # Remove duplicates
        effects = list(set(effects))
        
        return FeatureEffect(
            feature=feature,
            feature_type=feature.type,
            ref_feature_seq=ref_feature_seq,
            alt_feature_seq=alt_feature_seq,
            effects=effects,
            details=effect_details,
            variants=overlapping_variants or []
        )
    
    def _analyze_basic_effects(self, feature: Feature, ref_feature_seq: str, alt_feature_seq: str,
                             overlapping_variants: Optional[List[Variant]] = None) -> FeatureEffect:
        """
        Analyze basic effects without using sequence alignment.
        
        Args:
            feature: Feature to analyze
            ref_feature_seq: Reference feature sequence
            alt_feature_seq: Alternate feature sequence
            overlapping_variants: List of variants overlapping with this feature
            
        Returns:
            FeatureEffect object with analysis results
        """
        effects = []
        effect_details = {}
        
        # Calculate basic effects
        length_change = len(alt_feature_seq) - len(ref_feature_seq)
        
        if length_change != 0:
            effects.append('length_change')
            effect_details['length_change'] = length_change
        
        # Check for more specific effects based on feature type
        self._analyze_feature_type_specific_effects(feature, ref_feature_seq, alt_feature_seq, 
                                                  length_change, effects, effect_details)
        
        # Remove duplicates
        effects = list(set(effects))
        
        return FeatureEffect(
            feature=feature,
            feature_type=feature.type,
            ref_feature_seq=ref_feature_seq,
            alt_feature_seq=alt_feature_seq,
            effects=effects,
            details=effect_details,
            variants=overlapping_variants or []
        )
    
    def _analyze_feature_type_specific_effects(self, feature: Feature, ref_feature_seq: str, 
                                             alt_feature_seq: str, length_change: int,
                                             effects: List[str], effect_details: Dict) -> None:
        """
        Analyze feature-type specific effects.
        
        Args:
            feature: Feature to analyze
            ref_feature_seq: Reference feature sequence
            alt_feature_seq: Alternate feature sequence
            length_change: Length change between reference and alternate sequences
            effects: List of effects to update
            effect_details: Dictionary of effect details to update
        """
        feature_type = feature.type
        strand = feature.strand
        
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
            if feature.is_first_cds and alt_aa and (not alt_aa.startswith('M')):
                effects.append('start_codon_disruption')
        
        # Check for regulatory region effects
        if feature_type in ['promoter', 'terminator']:
            effects.append(f"{feature_type}_affected")
        
        # Check for splicing region effects
        if feature_type == 'exon':
            effects.append('splicing_affected')
