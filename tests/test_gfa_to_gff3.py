#!/usr/bin/env python3
"""
Tests for the GFA to GFF3 conversion functionality.
"""

import os
import tempfile
import unittest
from hapli.tools import gfa_to_gff3

class GfaToGff3Tests(unittest.TestCase):
    """Test cases for GFA to GFF3 conversion."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.temp_dir.name
        
        # Create test GFA file
        self.gfa_file = os.path.join(self.output_dir, "test.gfa")
        self.output_file = os.path.join(self.output_dir, "test.gff3")
        
        with open(self.gfa_file, 'w') as f:
            f.write("H\tVN:Z:1.0\n")
            f.write("S\ts1\tACGTACGTACGT\tLN:i:12\n")
            f.write("S\ts2\tGCTAGCTAGCTA\tLN:i:12\n")
            f.write("L\ts1\t+\ts2\t+\t0M\n")
            f.write("P\tREF\ts1+,s2+\t*\n")
    
    def tearDown(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()
    
    def test_parse_gfa(self):
        """Test parsing a GFA file."""
        segments, links, paths = gfa_to_gff3.parse_gfa(self.gfa_file)
        
        # Check segments
        self.assertEqual(len(segments), 2)
        self.assertIn('s1', segments)
        self.assertIn('s2', segments)
        self.assertEqual(segments['s1'], 'ACGTACGTACGT')
        self.assertEqual(segments['s2'], 'GCTAGCTAGCTA')
        
        # Check links
        self.assertEqual(len(links), 1)
        self.assertEqual(links[0], ('s1', '+', 's2', '+', '0M'))
        
        # Check paths
        self.assertEqual(len(paths), 1)
        self.assertIn('REF', paths)
        self.assertEqual(len(paths['REF']), 2)
        self.assertEqual(paths['REF'][0], ('s1', '+'))
        self.assertEqual(paths['REF'][1], ('s2', '+'))
    
    def test_build_reference_sequence(self):
        """Test building a reference sequence from GFA segments."""
        segments, links, paths = gfa_to_gff3.parse_gfa(self.gfa_file)
        ref_path = paths['REF']
        
        reference_seq = gfa_to_gff3.build_reference_sequence(segments, ref_path)
        
        # Expected sequence: s1 + s2
        expected = 'ACGTACGTACGTGCTAGCTAGCTA'
        
        self.assertEqual(reference_seq, expected)
    
    def test_create_synthetic_features(self):
        """Test creating synthetic features."""
        reference_seq = 'ACGTACGTACGTGCTAGCTAGCTA'
        
        features = gfa_to_gff3.create_synthetic_features(
            reference_seq, num_genes=2, num_exons_per_gene=2
        )
        
        # Check features
        self.assertGreaterEqual(len(features), 2 * (1 + 1 + 2 * 2))  # 2 genes, 2 mRNAs, 4 exons, 4 CDSs
        
        # Check feature types
        feature_types = [f['type'] for f in features]
        self.assertIn('gene', feature_types)
        self.assertIn('mRNA', feature_types)
        self.assertIn('exon', feature_types)
        self.assertIn('CDS', feature_types)
    
    def test_write_gff3(self):
        """Test writing features to a GFF3 file."""
        # Create some test features
        features = [
            {
                'type': 'gene',
                'start': 1,
                'end': 100,
                'strand': '+',
                'attributes': {'ID': 'gene1', 'Name': 'test_gene'}
            },
            {
                'type': 'mRNA',
                'start': 1,
                'end': 100,
                'strand': '+',
                'attributes': {'ID': 'mRNA1', 'Parent': 'gene1'}
            },
            {
                'type': 'exon',
                'start': 1,
                'end': 50,
                'strand': '+',
                'attributes': {'ID': 'exon1', 'Parent': 'mRNA1'}
            },
            {
                'type': 'CDS',
                'start': 1,
                'end': 50,
                'strand': '+',
                'phase': 0,
                'attributes': {'ID': 'cds1', 'Parent': 'mRNA1'}
            }
        ]
        
        # Write to GFF3
        gfa_to_gff3.write_gff3(features, self.output_file)
        
        # Check that the file was created
        self.assertTrue(os.path.exists(self.output_file))
        
        # Check file content
        with open(self.output_file, 'r') as f:
            content = f.read()
            self.assertIn('##gff-version 3', content)
            self.assertIn('gene', content)
            self.assertIn('mRNA', content)
            self.assertIn('exon', content)
            self.assertIn('CDS', content)
            self.assertIn('ID=gene1', content)
            self.assertIn('Parent=gene1', content)
            self.assertIn('Parent=mRNA1', content)

if __name__ == '__main__':
    unittest.main()
