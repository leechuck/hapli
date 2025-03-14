#!/usr/bin/env python3
"""
Tests for the variant analysis functionality.
"""

import unittest
from hapli.analysis import variant_analysis

class VariantAnalysisTests(unittest.TestCase):
    """Test cases for variant analysis functionality."""
    
    def test_calculate_variant_length_changes(self):
        """Test calculating length changes for variants."""
        # Define test variants and segments
        variants = [
            # SNP (no length change)
            {'id': 'var1', 'type': 'SNP', 'ref': 'A', 'alt': 'G'},
            
            # Insertion (length increase)
            {'id': 'var2', 'type': 'INS', 'ref': 'A', 'alt': 'ACGT'},
            
            # Deletion (length decrease)
            {'id': 'var3', 'type': 'DEL', 'ref': 'ACGT', 'alt': 'A'},
            
            # Complex variant
            {'id': 'var4', 'type': 'OTHER', 'ref': 'ACGT', 'alt': 'TGCA'}
        ]
        
        segments = {}  # Not used in this function
        
        # Calculate length changes
        result = variant_analysis.calculate_variant_length_changes(variants, segments)
        
        # Check results
        self.assertEqual(result['var1'], 0)  # SNP: no change
        self.assertEqual(result['var2'], 3)  # INS: +3 bases
        self.assertEqual(result['var3'], -3)  # DEL: -3 bases
        self.assertEqual(result['var4'], 0)  # Complex: same length
    
    def test_infer_variant_types(self):
        """Test inferring variant types."""
        # Define test variants
        variants = [
            # Explicit types
            {'id': 'var1', 'type': 'SNP', 'ref': 'A', 'alt': 'G'},
            {'id': 'var2', 'type': 'INS', 'ref': 'A', 'alt': 'ACGT'},
            {'id': 'var3', 'type': 'DEL', 'ref': 'ACGT', 'alt': 'A'},
            
            # Inferred types
            {'id': 'var4', 'ref': 'A', 'alt': 'G'},  # SNP
            {'id': 'var5', 'ref': 'A', 'alt': 'ACGT'},  # INS
            {'id': 'var6', 'ref': 'ACGT', 'alt': 'A'},  # DEL
            {'id': 'var7', 'ref': 'ACGT', 'alt': 'TGCA'}  # MNP
        ]
        
        segments = {}  # Not used in this function
        
        # Infer variant types
        result = variant_analysis.infer_variant_types(variants, segments)
        
        # Check results
        self.assertEqual(result[0]['type'], 'SNP')  # Explicit SNP
        self.assertEqual(result[1]['type'], 'INS')  # Explicit INS
        self.assertEqual(result[2]['type'], 'DEL')  # Explicit DEL
        self.assertEqual(result[3]['type'], 'SNP')  # Inferred SNP
        self.assertEqual(result[4]['type'], 'INS')  # Inferred INS
        self.assertEqual(result[5]['type'], 'DEL')  # Inferred DEL
        self.assertEqual(result[6]['type'], 'MNP')  # Inferred MNP
    
    def test_identify_repeat_sequences(self):
        """Test identifying repeat sequences around variants."""
        # Define test variants and segments
        segments = {
            's1': 'ACGTACGTACGT'  # 12 bases with repeating ACGT pattern
        }
        
        variants = [
            # Variant in a repetitive region
            {
                'id': 'var1',
                'segment': 's1',
                'position': 5,  # 0-based position in the segment
                'ref': 'A',
                'alt': 'G'
            },
            
            # Variant not in a repetitive region
            {
                'id': 'var2',
                'segment': 's1',
                'position': 1,
                'ref': 'C',
                'alt': 'T'
            }
        ]
        
        # Identify repeat sequences
        for variant in variants:
            result = variant_analysis.identify_repeat_sequences(variant, segments)
            
            # Check that the function returns a dictionary with expected keys
            self.assertIsInstance(result, dict)
            self.assertIn('left_repeat', result)
            self.assertIn('right_repeat', result)
            self.assertIn('is_in_repeat', result)
            
            # For var1, we expect to find repeats
            if variant['id'] == 'var1':
                self.assertTrue(result['is_in_repeat'])
                self.assertGreater(len(result['left_repeat']), 0)
                self.assertGreater(len(result['right_repeat']), 0)

if __name__ == '__main__':
    unittest.main()
