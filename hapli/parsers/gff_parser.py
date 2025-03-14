"""
Parser for GFF3 (General Feature Format version 3) files.
"""

import time
import logging
from collections import defaultdict

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
