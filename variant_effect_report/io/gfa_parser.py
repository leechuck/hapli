"""
Parser for GFA (Graphical Fragment Assembly) files.
"""

import logging
import re
import time
from typing import Dict, List, Tuple, Set, Optional
from collections import defaultdict

from variant_effect_report.core.models import Segment, Link, Path, Variant, VariantType


class GFAParser:
    """Parser for GFA files."""
    
    def __init__(self):
        """Initialize the GFA parser."""
        self.segments = {}
        self.links = []
        self.paths = {}
        self.variants = []
        self.samples = defaultdict(list)
        self.haplotypes = defaultdict(dict)
    
    def parse(self, gfa_file: str) -> Tuple[Dict, List, Dict, List, Dict, Dict]:
        """
        Parse a GFA file and extract segments, paths, and variant information.
        
        Args:
            gfa_file: Path to the GFA file
            
        Returns:
            Tuple of (segments, links, paths, variants, samples, haplotypes)
        """
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
                    self._parse_segment(fields, line_num)
                elif record_type == 'L':  # Link
                    self._parse_link(fields, line_num)
        
        # Second pass: collect path information
        with open(gfa_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                    
                fields = line.split('\t')
                record_type = fields[0]
                
                if record_type == 'P':  # Path
                    self._parse_path(fields, line_num)
        
        # Find or create REF path
        self._ensure_ref_path()
        
        # Calculate accurate length changes for variants
        self._calculate_variant_length_changes()
        
        # Infer variant types
        self._infer_variant_types()
        
        elapsed = time.time() - start_time
        logging.info(f"Finished parsing GFA in {elapsed:.2f}s: {len(self.segments)} segments, "
                     f"{len(self.paths)} paths, {len(self.variants)} variants")
        logging.info(f"Found {len(self.samples)} samples with paths")
        
        return (self.segments, self.links, self.paths, self.variants, 
                self.samples, self.haplotypes)
    
    def _parse_segment(self, fields: List[str], line_num: int) -> None:
        """Parse a segment record."""
        if len(fields) < 3:
            logging.warning(f"Line {line_num}: Invalid segment record, missing fields")
            return
            
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
        orig_seg_id = None
        if '_VAR' in seg_id:
            parts = seg_id.split('_VAR')
            if len(parts) > 1:
                # The original segment ID is the part before _VAR
                orig_seg_id = parts[0]
                # Extract variant ID from segment name
                var_id = 'VAR' + parts[1].split('_')[0]  # Handle multiple VAR tags
            
            # Try to infer variant type from segment name or tags
            if "SNP" in seg_id or "SNV" in seg_id:
                var_type = VariantType.SNP
            elif "DEL" in seg_id:
                var_type = VariantType.DEL
            elif "INS" in seg_id:
                var_type = VariantType.INS
            elif "DUP" in seg_id:
                var_type = VariantType.DUP
            elif "INV" in seg_id:
                var_type = VariantType.INV
        
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
        self.segments[seg_id] = {
            'sequence': sequence,
            'length': len(sequence),
            'tags': tags,
            'variant_id': var_id,
            'variant_type': var_type,
            'length_change': length_change,
            'original_segment': orig_seg_id
        }
        
        # Add to variants list if this is a variant segment
        if var_id:
            # Check if variant already exists
            existing_var = next((v for v in self.variants if v['id'] == var_id), None)
            if existing_var:
                # Update existing variant entry
                if seg_id not in existing_var['segments']:
                    existing_var['segments'].append(seg_id)
            else:
                # Create new variant entry
                self.variants.append({
                    'id': var_id,
                    'type': var_type.value if var_type else VariantType.UNKNOWN.value,
                    'segments': [seg_id],
                    'length_change': length_change
                })
    
    def _parse_link(self, fields: List[str], line_num: int) -> None:
        """Parse a link record."""
        if len(fields) < 6:
            logging.warning(f"Line {line_num}: Invalid link record, missing fields")
            return
            
        from_id = fields[1]
        from_dir = fields[2]
        to_id = fields[3]
        to_dir = fields[4]
        overlap = fields[5]
        
        self.links.append((from_id, from_dir, to_id, to_dir, overlap))
    
    def _parse_path(self, fields: List[str], line_num: int) -> None:
        """Parse a path record."""
        if len(fields) < 3:
            logging.warning(f"Line {line_num}: Invalid path record, missing fields")
            return
            
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
            has_variant_segments = any(seg_id in self.segments and self.segments[seg_id].get('variant_id') 
                                       for seg_id, _ in path_segments)
            if has_variant_segments and path_name != 'REF':
                sample_name = f"Sample_{path_name}"
                haplotype = "haplotype_1"
        
        # Register sample and haplotype
        if sample_name:
            self.samples[sample_name].append(path_name)
            if haplotype:
                self.haplotypes[sample_name][haplotype] = path_name
        
        # Extract variant IDs from path segments if not already defined
        if not variant_ids:
            for seg_id, _ in path_segments:
                if seg_id in self.segments:
                    var_id = self.segments[seg_id].get('variant_id')
                    if var_id and var_id not in variant_ids:
                        variant_ids.append(var_id)
        
        self.paths[path_name] = {
            'segments': path_segments,
            'tags': path_tags,
            'sample': sample_name,
            'haplotype': haplotype,
            'variant_ids': variant_ids
        }
        
        # Associate variants with paths and samples
        for var_id in variant_ids:
            # Find the corresponding variant
            variant = next((v for v in self.variants if v['id'] == var_id), None)
            if variant:
                # Add path information to variant
                if 'paths' not in variant:
                    variant['paths'] = []
                if path_name not in variant['paths']:
                    variant['paths'].append(path_name)
                
                # Add sample information to variant
                if sample_name:
                    if 'samples' not in variant:
                        variant['samples'] = []
                    if sample_name not in variant['samples']:
                        variant['samples'].append(sample_name)
    
    def _ensure_ref_path(self) -> None:
        """Find or create a reference path."""
        if 'REF' not in self.paths:
            logging.info("REF path not found, looking for reference path")
            # Try to find a path without variants
            for name, path_data in self.paths.items():
                if not path_data['variant_ids'] and not path_data['sample']:
                    self.paths['REF'] = path_data
                    logging.info(f"Using {name} as reference path")
                    break
        
        # If still no REF path, create one from non-variant segments
        if 'REF' not in self.paths:
            logging.info("Creating synthetic REF path from non-variant segments")
            ref_segments = []
            non_variant_segments = [seg_id for seg_id, data in self.segments.items() 
                                   if not data.get('variant_id')]
            
            # Sort segments by any available numbering
            def extract_number(seg_id):
                match = re.search(r'S(\d+)', seg_id)
                return int(match.group(1)) if match else 0
            
            sorted_segments = sorted(non_variant_segments, key=extract_number)
            
            if sorted_segments:
                ref_segments = [(seg_id, '+') for seg_id in sorted_segments]
                self.paths['REF'] = {
                    'segments': ref_segments,
                    'tags': {},
                    'sample': None,
                    'haplotype': None,
                    'variant_ids': []
                }
            else:
                logging.warning("Could not create REF path, no non-variant segments found")
    
    def _calculate_variant_length_changes(self) -> None:
        """Calculate length changes for variants."""
        for variant in self.variants:
            # Skip if length change is already set and non-zero
            if variant.get('length_change', 0) != 0:
                continue
                
            var_segments = variant.get('segments', [])
            if not var_segments:
                continue
                
            total_length_change = 0
            
            for seg_id in var_segments:
                if seg_id not in self.segments:
                    continue
                    
                # Get variant segment length
                var_seg_length = self.segments[seg_id].get('length', 0)
                
                # Try to get original segment ID
                orig_seg_id = self.segments[seg_id].get('original_segment')
                
                # If not directly specified, try to extract from segment name
                if not orig_seg_id and '_VAR' in seg_id:
                    parts = seg_id.split('_VAR')
                    if len(parts) > 1:
                        orig_seg_id = parts[0]
                
                # Calculate length change if we found the original segment
                if orig_seg_id and orig_seg_id in self.segments:
                    orig_seg_length = self.segments[orig_seg_id].get('length', 0)
                    segment_length_change = var_seg_length - orig_seg_length
                    total_length_change += segment_length_change
                # Get length change from the segment's tag if available
                elif 'LD' in self.segments[seg_id].get('tags', {}):
                    tag_data = self.segments[seg_id]['tags']['LD']
                    if tag_data[0] == 'i':  # Integer type
                        try:
                            segment_length_change = int(tag_data[1])
                            total_length_change += segment_length_change
                        except ValueError:
                            pass
            
            # Update the variant with the calculated length change
            variant['length_change'] = total_length_change
    
    def _infer_variant_types(self) -> None:
        """Infer variant types based on segment information and length changes."""
        for variant in self.variants:
            # Skip if type is already defined and not UNKNOWN
            if variant.get('type') and variant.get('type') != VariantType.UNKNOWN.value:
                continue
                
            length_change = variant.get('length_change', 0)
            var_segments = variant.get('segments', [])
            
            # Look for variant type hints in segment names
            for seg_id in var_segments:
                if seg_id not in self.segments:
                    continue
                    
                # Check segment name for type hints
                if "SNP" in seg_id or "SNV" in seg_id:
                    variant['type'] = VariantType.SNP.value
                    break
                elif "DEL" in seg_id:
                    variant['type'] = VariantType.DEL.value
                    break
                elif "INS" in seg_id:
                    variant['type'] = VariantType.INS.value
                    break
                elif "DUP" in seg_id:
                    variant['type'] = VariantType.DUP.value
                    break
                elif "INV" in seg_id:
                    variant['type'] = VariantType.INV.value
                    break
            
            # If no type found from names, infer from length change
            if (not variant.get('type') or variant.get('type') == VariantType.UNKNOWN.value) and length_change != 0:
                if length_change > 0:
                    # Check for duplication - if the variant segment is much longer
                    for seg_id in var_segments:
                        if seg_id in self.segments:
                            seg_data = self.segments[seg_id]
                            orig_seg_id = seg_data.get('original_segment')
                            
                            if orig_seg_id and orig_seg_id in self.segments:
                                var_seq = seg_data.get('sequence', '')
                                orig_seq = self.segments[orig_seg_id].get('sequence', '')
                                
                                # If the variant contains the original sequence multiple times
                                if orig_seq and len(orig_seq) > 0 and var_seq.count(orig_seq) > 1:
                                    variant['type'] = VariantType.DUP.value
                                    break
                    
                    # Default to insertion if not a duplication
                    if not variant.get('type') or variant.get('type') == VariantType.UNKNOWN.value:
                        variant['type'] = VariantType.INS.value
                elif length_change < 0:
                    variant['type'] = VariantType.DEL.value
            
            # Default to SNP for small zero-length-change variants without a type
            if (not variant.get('type') or variant.get('type') == VariantType.UNKNOWN.value) and length_change == 0:
                variant['type'] = VariantType.SNP.value
