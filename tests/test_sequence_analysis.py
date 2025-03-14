#!/usr/bin/env python3
"""
Tests for the sequence analysis functionality.
"""

import unittest
from hapli.analysis import sequence_analysis

class SequenceAnalysisTests(unittest.TestCase):
    """Test cases for sequence analysis functionality."""
    
    def test_build_path_sequence(self):
        """Test building a sequence from path segments."""
        # Define test segments and path
        segments = {
            's1': 'ACGT',
            's2': 'GCTA',
            's3': 'TGCA'
        }
        
        path_segments = [
            ('s1', '+'),  # Forward orientation
            ('s2', '+'),  # Forward orientation
            ('s3', '-')   # Reverse orientation (should be reverse complemented)
        ]
        
        # Expected result: s1 + s2 + reverse_complement(s3)
        expected = 'ACGTGCTATGCA'  # TGCA reverse complemented is TGCA
        
        # Test the function
        result = sequence_analysis.build_path_sequence(segments, path_segments)
        
        # Check the result
        self.assertEqual(result, expected)
    
    def test_reverse_complement(self):
        """Test reverse complementing a DNA sequence."""
        # Test cases
        test_cases = [
            ('ACGT', 'ACGT'),  # Palindromic sequence
            ('AAAAAA', 'TTTTTT'),
            ('GCTA', 'TAGC'),
            ('ATGC', 'GCAT'),
            ('N', 'N'),  # Ambiguous base
            ('', '')  # Empty string
        ]
        
        for seq, expected in test_cases:
            result = sequence_analysis.reverse_complement(seq)
            self.assertEqual(result, expected)
    
    def test_translate_sequence(self):
        """Test translating a DNA sequence to amino acids."""
        # Test cases for forward strand
        test_cases_forward = [
            ('ATGGCCTAA', 'MA*'),  # Start codon, alanine, stop codon
            ('ATGCCCGGG', 'MPG'),  # No stop codon
            ('CCCTTT', 'PF'),      # No start/stop codons
            ('', '')               # Empty string
        ]
        
        for seq, expected in test_cases_forward:
            result = sequence_analysis.translate_sequence(seq, strand='+')
            self.assertEqual(result, expected)
        
        # Test cases for reverse strand
        test_cases_reverse = [
            ('TTACCGGAT', 'MA*'),  # Reverse complement of ATGGCCTAA
            ('CCCGGGCAT', 'MPG'),  # Reverse complement of ATGCCCGGG
            ('AAAGGG', 'PF'),      # Reverse complement of CCCTTT
            ('', '')               # Empty string
        ]
        
        for seq, expected in test_cases_reverse:
            result = sequence_analysis.translate_sequence(seq, strand='-')
            self.assertEqual(result, expected)
    
    def test_compare_amino_acid_sequences(self):
        """Test comparing amino acid sequences."""
        # Test cases
        test_cases = [
            # Identical sequences
            ('MAPLG*', 'MAPLG*', {'identical': True, 'changes': []}),
            
            # Single substitution
            ('MAPLG*', 'MASLG*', {'identical': False, 'changes': [('P', 'S', 3)]}),
            
            # Multiple substitutions
            ('MAPLG*', 'MASLR*', {'identical': False, 'changes': [('P', 'S', 3), ('G', 'R', 5)]}),
            
            # Different lengths
            ('MAPLG*', 'MAPLGK*', {'identical': False, 'changes': [(None, 'K', 6)]}),
            ('MAPLGK*', 'MAPLG*', {'identical': False, 'changes': [('K', None, 6)]})
        ]
        
        for ref_aa, alt_aa, expected in test_cases:
            result = sequence_analysis.compare_amino_acid_sequences(ref_aa, alt_aa)
            self.assertEqual(result['identical'], expected['identical'])
            self.assertEqual(len(result['changes']), len(expected['changes']))
            
            # Check each change
            for i, (ref, alt, pos) in enumerate(expected['changes']):
                self.assertEqual(result['changes'][i][0], ref)
                self.assertEqual(result['changes'][i][1], alt)
                self.assertEqual(result['changes'][i][2], pos)
    
    def test_analyze_sequence_alignment(self):
        """Test analyzing sequence alignment."""
        # Test cases
        test_cases = [
            # Identical sequences
            ('ACGT', 'ACGT'),
            
            # Single substitution
            ('ACGT', 'ACTT'),
            
            # Insertion
            ('ACGT', 'ACGGT'),
            
            # Deletion
            ('ACGT', 'ACT'),
            
            # Complex changes
            ('ACGTACGT', 'ACTTACGGT')
        ]
        
        for ref_seq, alt_seq in test_cases:
            result = sequence_analysis.analyze_sequence_alignment(ref_seq, alt_seq)
            
            # Basic checks
            self.assertIsNotNone(result)
            self.assertIn('ref_aligned', result)
            self.assertIn('alt_aligned', result)
            self.assertIn('alignment_score', result)
            
            # Check that aligned sequences have the same length
            self.assertEqual(len(result['ref_aligned']), len(result['alt_aligned']))
            
            # Check that the alignment preserves the original sequences (ignoring gaps)
            self.assertEqual(ref_seq, result['ref_aligned'].replace('-', ''))
            self.assertEqual(alt_seq, result['alt_aligned'].replace('-', ''))

if __name__ == '__main__':
    unittest.main()
