#!/usr/bin/env python3
"""
Tests for the text reporting functionality.
"""

import os
import tempfile
import unittest
from io import StringIO
from unittest.mock import patch
from hapli.reporting import text_report

class TextReportTests(unittest.TestCase):
    """Test cases for text reporting functionality."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.temp_dir.name
    
    def tearDown(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()
    
    @patch('sys.stdout', new_callable=StringIO)
    def test_generate_variant_effect_report(self, mock_stdout):
        """Test generating a variant effect report."""
        # Define test data
        feature_effects = {
            'feature1': {
                'effects': ['effect1', 'effect2'],
                'zygosity': 'HOMOZYGOUS',
                'feature': {
                    'id': 'feature1',
                    'type': 'gene',
                    'start': 100,
                    'end': 200,
                    'strand': '+'
                },
                'sequence_changes': {
                    'reference': 'ACGT',
                    'alternate': 'ACTT'
                },
                'amino_acid_changes': [
                    ('A', 'T', 1)
                ]
            },
            'feature2': {
                'effects': ['effect3'],
                'zygosity': 'HETEROZYGOUS',
                'feature': {
                    'id': 'feature2',
                    'type': 'exon',
                    'start': 150,
                    'end': 180,
                    'strand': '+'
                },
                'sequence_changes': {
                    'reference': 'ACGT',
                    'alternate': 'ACGTA'
                }
            }
        }
        
        variants = [
            {
                'id': 'var1',
                'type': 'SNP',
                'pos': 100,
                'ref': 'A',
                'alt': 'T'
            },
            {
                'id': 'var2',
                'type': 'INS',
                'pos': 180,
                'ref': 'G',
                'alt': 'GA'
            }
        ]
        
        # Generate report to stdout
        text_report.generate_variant_effect_report(feature_effects, variants)
        
        # Check output
        output = mock_stdout.getvalue()
        self.assertIn('Variant Effect Report', output)
        self.assertIn('feature1', output)
        self.assertIn('feature2', output)
        self.assertIn('HOMOZYGOUS', output)
        self.assertIn('HETEROZYGOUS', output)
        self.assertIn('effect1', output)
        self.assertIn('effect2', output)
        self.assertIn('effect3', output)
    
    def test_generate_variant_effect_report_to_file(self):
        """Test generating a variant effect report to a file."""
        # Define test data
        feature_effects = {
            'feature1': {
                'effects': ['effect1', 'effect2'],
                'zygosity': 'HOMOZYGOUS',
                'feature': {
                    'id': 'feature1',
                    'type': 'gene',
                    'start': 100,
                    'end': 200,
                    'strand': '+'
                },
                'sequence_changes': {
                    'reference': 'ACGT',
                    'alternate': 'ACTT'
                },
                'amino_acid_changes': [
                    ('A', 'T', 1)
                ]
            }
        }
        
        variants = [
            {
                'id': 'var1',
                'type': 'SNP',
                'pos': 100,
                'ref': 'A',
                'alt': 'T'
            }
        ]
        
        # Output file
        output_file = os.path.join(self.output_dir, 'report.txt')
        
        # Generate report to file
        text_report.generate_variant_effect_report(feature_effects, variants, outfile=output_file)
        
        # Check that the file was created
        self.assertTrue(os.path.exists(output_file))
        
        # Check file content
        with open(output_file, 'r') as f:
            content = f.read()
            self.assertIn('Variant Effect Report', content)
            self.assertIn('feature1', content)
            self.assertIn('HOMOZYGOUS', content)
            self.assertIn('effect1', content)
            self.assertIn('effect2', content)

if __name__ == '__main__':
    unittest.main()
