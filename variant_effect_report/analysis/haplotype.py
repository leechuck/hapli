"""
Haplotype analysis for variant effect reporting.
"""

import logging
from typing import Dict, List, Optional, Tuple, Any

from variant_effect_report.core.models import Feature, Variant, FeatureEffect, SegmentOffset
from variant_effect_report.core.utils import build_path_sequence
from variant_effect_report.analysis.feature_effects import FeatureEffectAnalyzer


class HaplotypeAnalyzer:
    """Analyze haplotype differences."""
    
    def __init__(self, use_alignment: bool = False):
        """
        Initialize the haplotype analyzer.
        
        Args:
            use_alignment: Whether to use sequence alignment for analysis
        """
        self.use_alignment = use_alignment
        self.feature_effect_analyzer = FeatureEffectAnalyzer(use_alignment=use_alignment)
    
    def analyze_haplotype_differences(self, ref_seq: str, alt_seq: str, features: List[Feature], 
                                     segment_offsets: Optional[List[Dict]] = None, 
                                     variant_segments: Optional[Dict] = None,
                                     segments: Optional[Dict] = None) -> List[FeatureEffect]:
        """
        Analyze differences between reference and alternate haplotype sequences.
        
        Args:
            ref_seq: Reference sequence string
            alt_seq: Alternate sequence string
            features: List of feature objects
            segment_offsets: Optional segment offset information
            variant_segments: Optional variant segment information
            segments: Optional segment data
            
        Returns:
            List of feature effect objects
        """
        feature_effects = []
        
        # Identify variants from segment information if available
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
            
            variants = self._identify_variants_by_position(segment_offsets, ref_segments, variant_segments, segments)
        
        # Analyze each feature
        for feature in features:
            # Extract feature boundaries
            start = feature.start - 1  # Convert to 0-based
            end = feature.end - 1      # Convert to 0-based
            
            # Skip if feature is outside sequence bounds
            if start >= len(ref_seq) or end >= len(ref_seq):
                logging.warning(f"Feature {feature.attributes.get('ID', 'unknown')} is outside sequence bounds, skipping")
                continue
            
            # Find variants that overlap with this feature
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
                if not (var_end < start or var_start > end):
                    overlapping_variants.append(variant)
            
            # Analyze feature effects
            effect = self.feature_effect_analyzer.analyze_feature(
                feature, ref_seq, alt_seq, overlapping_variants
            )
            
            feature_effects.append(effect)
        
        return feature_effects
    
    def _identify_variants_by_position(self, segment_offsets: List[Dict], ref_segments: List[Dict], 
                                     variant_segments: Dict, segments: Optional[Dict] = None) -> List[Dict]:
        """
        Identify variants in an alternate sequence based on segment differences.
        
        Args:
            segment_offsets: List of segment offset dictionaries
            ref_segments: List of reference segment dictionaries
            variant_segments: Dictionary of variant segment data
            segments: Optional dictionary of all segment data
            
        Returns:
            List of variant dictionaries
        """
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
