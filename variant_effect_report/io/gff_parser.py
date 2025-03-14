"""
Parser for GFF3 (General Feature Format version 3) files.
"""

import logging
import time
from typing import Dict, List, Tuple
from collections import defaultdict

from variant_effect_report.core.models import Feature


class GFFParser:
    """Parser for GFF3 files."""
    
    def __init__(self):
        """Initialize the GFF parser."""
        self.features = []
        self.feature_by_id = {}
        self.children_by_parent = defaultdict(list)
    
    def parse(self, gff_file: str) -> Tuple[List[Feature], Dict[str, Feature], Dict[str, List[Feature]]]:
        """
        Parse a GFF3 file into a list of features.
        
        Args:
            gff_file: Path to the GFF3 file
            
        Returns:
            Tuple of (features, feature_by_id, children_by_parent)
        """
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
                
                feature = Feature(
                    seqid=seqid,
                    source=source,
                    type=feature_type,
                    start=start,
                    end=end,
                    score=score,
                    strand=strand,
                    phase=phase,
                    attributes=attributes,
                    line_num=line_num
                )
                
                self.features.append(feature)
        
        elapsed = time.time() - start_time
        logging.info(f"Finished parsing GFF3 in {elapsed:.2f}s: {len(self.features)} features")
        
        # Build feature hierarchy
        self._build_feature_hierarchy()
        
        # Identify first CDS in genes
        self._identify_first_cds_in_genes()
        
        return self.features, self.feature_by_id, self.children_by_parent
    
    def _build_feature_hierarchy(self) -> None:
        """Build the feature hierarchy based on Parent attributes."""
        # First pass: index features by ID
        for feature in self.features:
            if 'ID' in feature.attributes:
                feature_id = feature.attributes['ID']
                self.feature_by_id[feature_id] = feature
                
            if 'Parent' in feature.attributes:
                parent_id = feature.attributes['Parent']
                self.children_by_parent[parent_id].append(feature)
        
        # Second pass: attach children to their parents
        for feature in self.features:
            if 'ID' in feature.attributes:
                feature_id = feature.attributes['ID']
                feature.children = self.children_by_parent.get(feature_id, [])
    
    def _identify_first_cds_in_genes(self) -> None:
        """Identify the first CDS feature in each gene to check for start codon changes."""
        # For each mRNA, find all its CDS children and mark the first one
        for feature in self.features:
            if feature.type == 'mRNA':
                mrna_id = feature.attributes.get('ID')
                if not mrna_id:
                    continue
                    
                # Get all CDS features for this mRNA
                cds_features = [f for f in self.children_by_parent.get(mrna_id, []) if f.type == 'CDS']
                
                # Sort by position, considering strand
                if feature.strand == '+':
                    cds_features.sort(key=lambda x: x.start)
                else:
                    cds_features.sort(key=lambda x: x.start, reverse=True)
                
                # Mark the first CDS
                if cds_features:
                    cds_features[0].is_first_cds = True
