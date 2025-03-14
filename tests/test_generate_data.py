#!/usr/bin/env python3
"""
Tests for the test data generation functionality.
"""

import os
import tempfile
import unittest
from hapli.tools import generate_test_data

class TestDataGenerationTests(unittest.TestCase):
    """Test cases for test data generation."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.temp_dir.name
    
    def tearDown(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()
    
    def test_reference_sequence_generation(self):
        """Test that reference sequences are generated correctly."""
        seq_length = 1000
        sequence = generate_test_data.generate_random_sequence(seq_length)
        self.assertEqual(len(sequence), seq_length)
        # Check that sequence only contains valid DNA characters
        self.assertTrue(all(c in 'ACGT' for c in sequence))
    
    def test_variant_reference_matching(self):
        """Test that variants match the reference sequence."""
        # Generate a reference sequence
        seq_length = 1000
        reference_seq = generate_test_data.generate_random_sequence(seq_length)
        
        # Generate variants
        variants = []
        for i in range(1, 11):
            pos = i * 50  # Evenly space variants
            
            # SNP variant
            if i % 3 == 0:
                ref_base = reference_seq[pos-1]
                alt_bases = [b for b in 'ACGT' if b != ref_base]
                alt_base = alt_bases[0]
                variants.append({
                    'id': f"var{i}",
                    'type': 'SNP',
                    'pos': pos,
                    'ref': ref_base,
                    'alt': alt_base
                })
            # INS variant
            elif i % 3 == 1:
                ref_base = reference_seq[pos-1]
                ins_seq = "ACGT"
                variants.append({
                    'id': f"var{i}",
                    'type': 'INS',
                    'pos': pos,
                    'ref': ref_base,
                    'alt': ref_base + ins_seq
                })
            # DEL variant
            else:
                del_length = 3
                ref_seq = reference_seq[pos-1:pos+del_length]
                alt_base = reference_seq[pos-1]
                variants.append({
                    'id': f"var{i}",
                    'type': 'DEL',
                    'pos': pos,
                    'ref': ref_seq,
                    'alt': alt_base
                })
        
        # Validate variants
        is_valid, errors = generate_test_data.validate_variants_against_reference(variants, reference_seq)
        self.assertTrue(is_valid, f"Variants do not match reference: {errors}")
        
        # Test with deliberately mismatched variants
        bad_variants = variants.copy()
        bad_variants[0]['ref'] = 'X' + bad_variants[0]['ref'][1:] if bad_variants[0]['ref'] else 'X'
        
        is_valid, errors = generate_test_data.validate_variants_against_reference(bad_variants, reference_seq)
        self.assertFalse(is_valid, "Validation failed to detect mismatched variants")
    
    def test_gfa_generation(self):
        """Test GFA file generation."""
        gfa_file = os.path.join(self.output_dir, "test.gfa")
        segments, variants, ref_path, _ = generate_test_data.generate_gfa(
            gfa_file, seq_length=1000, num_segments=3, num_variants=5
        )
        
        # Check that the file was created
        self.assertTrue(os.path.exists(gfa_file))
        
        # Check that segments were created
        self.assertGreaterEqual(len(segments), 3)
        
        # Check that variants were created
        self.assertGreaterEqual(len(variants), 1)
        
        # Check that reference path was created
        self.assertEqual(len(ref_path), len(segments))
    
    def test_vcf_generation(self):
        """Test VCF file generation."""
        # Generate a reference sequence
        seq_length = 1000
        reference_seq = generate_test_data.generate_random_sequence(seq_length)
        
        # Generate VCF file
        vcf_file = os.path.join(self.output_dir, "test.vcf")
        variants, samples = generate_test_data.generate_vcf(
            vcf_file, sequence_length=seq_length, num_variants=5, 
            num_samples=2, reference_seq=reference_seq
        )
        
        # Check that the file was created
        self.assertTrue(os.path.exists(vcf_file))
        
        # Check that variants were created
        self.assertGreaterEqual(len(variants), 1)
        
        # Check that samples were created
        self.assertEqual(len(samples), 2)
        
        # Validate variants against reference
        is_valid, errors = generate_test_data.validate_variants_against_reference(variants, reference_seq)
        self.assertTrue(is_valid, f"Variants do not match reference: {errors}")
    
    def test_gff3_generation(self):
        """Test GFF3 file generation."""
        gff_file = os.path.join(self.output_dir, "test.gff3")
        num_genes, num_exons = generate_test_data.generate_gff3(
            gff_file, sequence_length=1000, num_genes=3, num_exons_per_gene=2
        )
        
        # Check that the file was created
        self.assertTrue(os.path.exists(gff_file))
        
        # Check that genes and exons were created
        self.assertEqual(num_genes, 3)
        self.assertEqual(num_exons, 2)
    
    def test_complete_dataset_generation(self):
        """Test generation of a complete dataset."""
        files = generate_test_data.generate_test_dataset(
            self.output_dir, prefix="test", sequence_length=1000
        )
        
        # Check that all files were created
        self.assertTrue(os.path.exists(files['fasta']))
        self.assertTrue(os.path.exists(files['gfa']))
        self.assertTrue(os.path.exists(files['gff3']))
        self.assertTrue(os.path.exists(files['vcf']))

if __name__ == '__main__':
    unittest.main()
