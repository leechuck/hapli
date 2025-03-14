#!/usr/bin/env python3
"""
Tests for the GFF parser functionality.
"""

import os
import tempfile
import unittest
from hapli.parsers import gff_parser

class GFFParserTests(unittest.TestCase):
    """Test cases for GFF parser functionality."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.temp_dir.name
        
        # Create test GFF file
        self.gff_file = os.path.join(self.output_dir, "test.gff3")
        
        with open(self.gff_file, 'w') as f:
            f.write("##gff-version 3\n")
            f.write("##sequence-region 1 1 1000\n")
            f.write("1\ttest\tgene\t1\t500\t.\t+\t.\tID=gene1;Name=test_gene\n")
            f.write("1\ttest\tmRNA\t1\t500\t.\t+\t.\tID=mRNA1;Parent=gene1\n")
            f.write("1\ttest\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=mRNA1\n")
            f.write("1\ttest\texon\t201\t300\t.\t+\t.\tID=exon2;Parent=mRNA1\n")
            f.write("1\ttest\texon\t401\t500\t.\t+\t.\tID=exon3;Parent=mRNA1\n")
            f.write("1\ttest\tCDS\t1\t100\t.\t+\t0\tID=cds1;Parent=mRNA1\n")
            f.write("1\ttest\tCDS\t201\t300\t.\t+\t0\tID=cds2;Parent=mRNA1\n")
            f.write("1\ttest\tCDS\t401\t500\t.\t+\t0\tID=cds3;Parent=mRNA1\n")
    
    def tearDown(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()
    
    def test_parse_gff3(self):
        """Test parsing a GFF3 file."""
        features, feature_by_id, children_by_parent = gff_parser.parse_gff3(self.gff_file)
        
        # Check that features were parsed correctly
        self.assertEqual(len(features), 8)  # 1 gene, 1 mRNA, 3 exons, 3 CDS
        
        # Check feature_by_id
        self.assertIn('gene1', feature_by_id)
        self.assertIn('mRNA1', feature_by_id)
        self.assertIn('exon1', feature_by_id)
        self.assertIn('cds1', feature_by_id)
        
        # Check children_by_parent
        self.assertIn('gene1', children_by_parent)
        self.assertEqual(len(children_by_parent['gene1']), 1)  # 1 mRNA
        
        self.assertIn('mRNA1', children_by_parent)
        self.assertEqual(len(children_by_parent['mRNA1']), 6)  # 3 exons + 3 CDS
        
        # Check feature attributes
        gene = feature_by_id['gene1']
        self.assertEqual(gene['type'], 'gene')
        self.assertEqual(gene['start'], 1)
        self.assertEqual(gene['end'], 500)
        self.assertEqual(gene['strand'], '+')
        self.assertEqual(gene['attributes']['ID'], 'gene1')
        self.assertEqual(gene['attributes']['Name'], 'test_gene')
    
    def test_parse_gff3_with_invalid_file(self):
        """Test parsing an invalid GFF3 file."""
        # Create an invalid GFF file
        invalid_gff = os.path.join(self.output_dir, "invalid.gff3")
        
        with open(invalid_gff, 'w') as f:
            f.write("This is not a valid GFF3 file\n")
        
        # This should not raise an exception, but should return empty structures
        features, feature_by_id, children_by_parent = gff_parser.parse_gff3(invalid_gff)
        
        self.assertEqual(len(features), 0)
        self.assertEqual(len(feature_by_id), 0)
        self.assertEqual(len(children_by_parent), 0)
    
    def test_parse_gff3_with_missing_parent(self):
        """Test parsing a GFF3 file with missing parent references."""
        # Create a GFF file with missing parent
        missing_parent_gff = os.path.join(self.output_dir, "missing_parent.gff3")
        
        with open(missing_parent_gff, 'w') as f:
            f.write("##gff-version 3\n")
            f.write("1\ttest\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=missing_mRNA\n")
        
        # This should not raise an exception
        features, feature_by_id, children_by_parent = gff_parser.parse_gff3(missing_parent_gff)
        
        # The feature should be parsed, but not linked to a parent
        self.assertEqual(len(features), 1)
        self.assertIn('exon1', feature_by_id)
        self.assertNotIn('missing_mRNA', children_by_parent)

if __name__ == '__main__':
    unittest.main()
