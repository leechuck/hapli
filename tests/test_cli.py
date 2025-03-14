#!/usr/bin/env python3
"""
Tests for the CLI functionality.
"""

import os
import sys
import tempfile
import unittest
from unittest.mock import patch, MagicMock
from io import StringIO
from hapli.cli import main

class CLITests(unittest.TestCase):
    """Test cases for CLI functionality."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.temp_dir.name
        
        # Create test files
        self.gfa_file = os.path.join(self.output_dir, "test.gfa")
        self.gff_file = os.path.join(self.output_dir, "test.gff3")
        
        # Create minimal test files
        with open(self.gfa_file, 'w') as f:
            f.write("H\tVN:Z:1.0\n")
            f.write("S\ts1\tACGT\n")
            f.write("P\tREF\ts1+\t*\n")
        
        with open(self.gff_file, 'w') as f:
            f.write("##gff-version 3\n")
            f.write("##sequence-region 1 1 100\n")
            f.write("1\ttest\tgene\t1\t100\t.\t+\t.\tID=gene1;Name=test_gene\n")
    
    def tearDown(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()
    
    @patch('sys.stdout', new_callable=StringIO)
    def test_main_function_basic(self, mock_stdout):
        """Test that the main function runs without errors."""
        # Create a mock args object
        args = MagicMock()
        args.gfa_file = self.gfa_file
        args.gff_file = self.gff_file
        args.output = None
        args.output_dir = self.output_dir
        args.output_format = 'text'
        args.format = 'text'
        args.rdf_format = 'turtle'
        args.base_uri = 'http://example.org/genomics/'
        args.consolidated = False
        args.sample_reports = False
        args.output_prefix = None
        args.samples = None
        args.use_alignment = False
        args.debug = False
        args.verbose = False
        args.log_file = None
        args.save_shex = None
        
        # Run the main function
        result = main(args)
        
        # Check that the function returned successfully
        self.assertEqual(result, 0)
    
    @patch('sys.stdout', new_callable=StringIO)
    def test_main_function_with_output_file(self, mock_stdout):
        """Test main function with output to file."""
        # Create a mock args object
        args = MagicMock()
        args.gfa_file = self.gfa_file
        args.gff_file = self.gff_file
        args.output = os.path.join(self.output_dir, "output.txt")
        args.output_dir = self.output_dir
        args.output_format = 'text'
        args.format = 'text'
        args.rdf_format = 'turtle'
        args.base_uri = 'http://example.org/genomics/'
        args.consolidated = False
        args.sample_reports = False
        args.output_prefix = None
        args.samples = None
        args.use_alignment = False
        args.debug = False
        args.verbose = False
        args.log_file = None
        args.save_shex = None
        
        # Run the main function
        result = main(args)
        
        # Check that the function returned successfully
        self.assertEqual(result, 0)
        
        # Check that the output file was created
        self.assertTrue(os.path.exists(args.output))
    
    @patch('sys.stdout', new_callable=StringIO)
    def test_main_function_with_rdf_output(self, mock_stdout):
        """Test main function with RDF output."""
        # Create a mock args object
        args = MagicMock()
        args.gfa_file = self.gfa_file
        args.gff_file = self.gff_file
        args.output = os.path.join(self.output_dir, "output.ttl")
        args.output_dir = self.output_dir
        args.output_format = 'rdf'
        args.format = 'rdf'
        args.rdf_format = 'turtle'
        args.base_uri = 'http://example.org/genomics/'
        args.consolidated = False
        args.sample_reports = False
        args.output_prefix = None
        args.samples = None
        args.use_alignment = False
        args.debug = False
        args.verbose = False
        args.log_file = None
        args.save_shex = None
        
        # Run the main function
        try:
            result = main(args)
            # Check that the function returned successfully
            self.assertEqual(result, 0)
            # Check that the output file was created
            self.assertTrue(os.path.exists(args.output))
        except ImportError:
            # Skip if RDFLib is not available
            self.skipTest("RDFLib not available")

if __name__ == '__main__':
    unittest.main()
