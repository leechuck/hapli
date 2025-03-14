#!/usr/bin/env python3
"""
Tests for the feature analysis functionality.
"""

import unittest
from hapli.analysis import feature_analysis

class FeatureAnalysisTests(unittest.TestCase):
    """Test cases for feature analysis functionality."""
    
    def test_identify_first_cds_in_genes(self):
        """Test identifying the first CDS in each gene."""
        # Define test features
        features = [
            {
                'id': 'gene1',
                'type': 'gene',
                'start': 1,
                'end': 1000,
                'strand': '+',
                'attributes': {'ID': 'gene1'}
            },
            {
                'id': 'mRNA1',
                'type': 'mRNA',
                'start': 1,
                'end': 1000,
                'strand': '+',
                'attributes': {'ID': 'mRNA1', 'Parent': 'gene1'}
            },
            {
                'id': 'cds1.1',
                'type': 'CDS',
                'start': 100,
                'end': 200,
                'strand': '+',
                'attributes': {'ID': 'cds1.1', 'Parent': 'mRNA1'}
            },
            {
                'id': 'cds1.2',
                'type': 'CDS',
                'start': 300,
                'end': 400,
                'strand': '+',
                'attributes': {'ID': 'cds1.2', 'Parent': 'mRNA1'}
            },
            {
                'id': 'gene2',
                'type': 'gene',
                'start': 2000,
                'end': 3000,
                'strand': '-',  # Reverse strand
                'attributes': {'ID': 'gene2'}
            },
            {
                'id': 'mRNA2',
                'type': 'mRNA',
                'start': 2000,
                'end': 3000,
                'strand': '-',
                'attributes': {'ID': 'mRNA2', 'Parent': 'gene2'}
            },
            {
                'id': 'cds2.1',
                'type': 'CDS',
                'start': 2100,
                'end': 2200,
                'strand': '-',
                'attributes': {'ID': 'cds2.1', 'Parent': 'mRNA2'}
            },
            {
                'id': 'cds2.2',
                'type': 'CDS',
                'start': 2300,
                'end': 2400,
                'strand': '-',
                'attributes': {'ID': 'cds2.2', 'Parent': 'mRNA2'}
            }
        ]
        
        # Create feature_by_id and children_by_parent dictionaries
        feature_by_id = {f['id']: f for f in features}
        
        children_by_parent = {}
        for f in features:
            if 'Parent' in f['attributes']:
                parent_id = f['attributes']['Parent']
                if parent_id not in children_by_parent:
                    children_by_parent[parent_id] = []
                children_by_parent[parent_id].append(f)
        
        # Identify first CDS in genes
        first_cds = feature_analysis.identify_first_cds_in_genes(
            features, feature_by_id, children_by_parent
        )
        
        # Check results
        self.assertEqual(len(first_cds), 2)  # Two genes
        
        # For gene1 (forward strand), the first CDS should be the one with the lowest start position
        self.assertEqual(first_cds['gene1']['id'], 'cds1.1')
        
        # For gene2 (reverse strand), the first CDS should be the one with the highest end position
        self.assertEqual(first_cds['gene2']['id'], 'cds2.2')
    
    def test_analyze_sample_haplotypes(self):
        """Test analyzing sample haplotypes."""
        # Define test data
        feature_effects_by_haplotype = {
            'sample1': {
                'haplotype1': {
                    'feature1': {'effects': ['effect1', 'effect2']},
                    'feature2': {'effects': ['effect3']}
                },
                'haplotype2': {
                    'feature1': {'effects': ['effect1']},
                    'feature3': {'effects': ['effect4']}
                }
            },
            'sample2': {
                'haplotype1': {
                    'feature2': {'effects': ['effect3']},
                    'feature3': {'effects': ['effect4']}
                },
                'haplotype2': {
                    'feature1': {'effects': ['effect1', 'effect2']},
                    'feature2': {'effects': ['effect3']}
                }
            }
        }
        
        features = [
            {'id': 'feature1', 'type': 'gene'},
            {'id': 'feature2', 'type': 'exon'},
            {'id': 'feature3', 'type': 'CDS'}
        ]
        
        # Analyze sample haplotypes
        result = feature_analysis.analyze_sample_haplotypes(
            feature_effects_by_haplotype, features
        )
        
        # Check results
        self.assertEqual(len(result), 2)  # Two samples
        
        # Check sample1
        sample1 = result['sample1']
        self.assertEqual(sample1['total_features'], 3)
        self.assertEqual(sample1['affected_features'], 3)
        self.assertEqual(sample1['homozygous_effects'], 1)  # feature2 affected in both haplotypes
        self.assertEqual(sample1['heterozygous_effects'], 2)  # feature1 and feature3 affected in one haplotype
        
        # Check sample2
        sample2 = result['sample2']
        self.assertEqual(sample2['total_features'], 3)
        self.assertEqual(sample2['affected_features'], 3)
        self.assertEqual(sample2['homozygous_effects'], 1)  # feature2 affected in both haplotypes
        self.assertEqual(sample2['heterozygous_effects'], 2)  # feature1 and feature3 affected in one haplotype
    
    def test_identify_compound_heterozygous_effects(self):
        """Test identifying compound heterozygous effects."""
        # Define test data
        feature_effects_by_haplotype = {
            'sample1': {
                'haplotype1': {
                    'gene1': {'effects': ['effect1']},
                    'gene2': {'effects': ['effect2']}
                },
                'haplotype2': {
                    'gene1': {'effects': ['effect3']},
                    'gene3': {'effects': ['effect4']}
                }
            }
        }
        
        features = [
            {'id': 'gene1', 'type': 'gene'},
            {'id': 'gene2', 'type': 'gene'},
            {'id': 'gene3', 'type': 'gene'}
        ]
        
        feature_by_id = {f['id']: f for f in features}
        children_by_parent = {}
        
        # Identify compound heterozygous effects
        result = feature_analysis.identify_compound_heterozygous_effects(
            feature_effects_by_haplotype, features, feature_by_id, children_by_parent
        )
        
        # Check results
        self.assertEqual(len(result), 1)  # One sample
        
        # Check sample1
        sample1 = result['sample1']
        self.assertEqual(len(sample1), 1)  # One gene with compound heterozygous effects
        
        # gene1 should have compound heterozygous effects
        self.assertIn('gene1', sample1)
        self.assertEqual(sample1['gene1']['haplotype1']['effects'], ['effect1'])
        self.assertEqual(sample1['gene1']['haplotype2']['effects'], ['effect3'])

if __name__ == '__main__':
    unittest.main()
