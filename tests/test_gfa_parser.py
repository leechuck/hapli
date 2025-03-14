#!/usr/bin/env python3
"""
Tests for the GFA parser functionality.
"""

import os
import tempfile
import unittest
from hapli.parsers import gfa_parser

class GFAParserTests(unittest.TestCase):
    """Test cases for GFA parser functionality."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.temp_dir.name
        
        # Create test GFA file
        self.gfa_file = os.path.join(self.output_dir, "test.gfa")
        
        with open(self.gfa_file, 'w') as f:
            f.write("H\tVN:Z:1.0\n")
            f.write("S\ts1\tACGT\tLN:i:4\n")
            f.write("S\ts2\tGCTA\tLN:i:4\n")
            f.write("L\ts1\t+\ts2\t+\t0M\n")
            f.write("P\tREF\ts1+,s2+\t*\n")
            f.write("P\tALT\ts1+,s2+\t*\tVN:Z:test\n")
    
    def tearDown(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()
    
    def test_parse_gfa(self):
        """Test parsing a GFA file."""
        result = gfa_parser.parse_gfa(self.gfa_file)
        
        # Check that the result contains the expected components
        self.assertIn('segments', result)
        self.assertIn('links', result)
        self.assertIn('paths', result)
        
        # Check segments
        self.assertEqual(len(result['segments']), 2)
        self.assertIn('s1', result['segments'])
        self.assertIn('s2', result['segments'])
        self.assertEqual(result['segments']['s1'], 'ACGT')
        self.assertEqual(result['segments']['s2'], 'GCTA')
        
        # Check links
        self.assertEqual(len(result['links']), 1)
        self.assertEqual(result['links'][0], ('s1', '+', 's2', '+', '0M'))
        
        # Check paths
        self.assertEqual(len(result['paths']), 2)
        self.assertIn('REF', result['paths'])
        self.assertIn('ALT', result['paths'])
        
        # Check path segments
        ref_path = result['paths']['REF']
        self.assertEqual(len(ref_path), 2)
        self.assertEqual(ref_path[0], ('s1', '+'))
        self.assertEqual(ref_path[1], ('s2', '+'))
    
    def test_parse_gfa_with_missing_ref_path(self):
        """Test parsing a GFA file with missing REF path."""
        # Create a GFA file without a REF path
        no_ref_gfa = os.path.join(self.output_dir, "no_ref.gfa")
        
        with open(no_ref_gfa, 'w') as f:
            f.write("H\tVN:Z:1.0\n")
            f.write("S\ts1\tACGT\tLN:i:4\n")
            f.write("S\ts2\tGCTA\tLN:i:4\n")
            f.write("L\ts1\t+\ts2\t+\t0M\n")
            f.write("P\tALT\ts1+,s2+\t*\n")
        
        result = gfa_parser.parse_gfa(no_ref_gfa)
        
        # Check that a synthetic REF path was created
        self.assertIn('REF', result['paths'])
        ref_path = result['paths']['REF']
        self.assertEqual(len(ref_path), 2)
    
    def test_parse_gfa_with_invalid_file(self):
        """Test parsing an invalid GFA file."""
        # Create an invalid GFA file
        invalid_gfa = os.path.join(self.output_dir, "invalid.gfa")
        
        with open(invalid_gfa, 'w') as f:
            f.write("This is not a valid GFA file\n")
        
        # This should not raise an exception, but should return empty structures
        result = gfa_parser.parse_gfa(invalid_gfa)
        
        self.assertEqual(len(result['segments']), 0)
        self.assertEqual(len(result['links']), 0)
        self.assertEqual(len(result['paths']), 0)

if __name__ == '__main__':
    unittest.main()
