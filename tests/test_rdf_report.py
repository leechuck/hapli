#!/usr/bin/env python3
"""
Tests for the RDF reporting functionality.
"""

import os
import tempfile
import unittest
from unittest.mock import patch
from hapli.reporting import rdf_report

class RDFReportTests(unittest.TestCase):
    """Test cases for RDF reporting functionality."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.temp_dir.name
        
        # Skip tests if RDFLib is not available
        try:
            import rdflib
        except ImportError:
            self.skipTest("RDFLib not available")
    
    def tearDown(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()
    
    def test_create_rdf_variant_report(self):
        """Test creating an RDF variant report."""
        # Define test data
        variants = [
            {
                'id': 'var1',
                'type': 'SNP',
                'pos': 100,
                'ref': 'A',
                'alt': 'G',
                'chrom': '1'
            },
            {
                'id': 'var2',
                'type': 'INS',
                'pos': 200,
                'ref': 'C',
                'alt': 'CGTA',
                'chrom': '1'
            }
        ]
        
        features = [
            {
                'id': 'gene1',
                'type': 'gene',
                'start': 50,
                'end': 150,
                'strand': '+',
                'attributes': {'ID': 'gene1', 'Name': 'test_gene'}
            },
            {
                'id': 'exon1',
                'type': 'exon',
                'start': 100,
                'end': 150,
                'strand': '+',
                'attributes': {'ID': 'exon1', 'Parent': 'gene1'}
            }
        ]
        
        # Create RDF report
        graph = rdf_report.create_rdf_variant_report(variants, features)
        
        # Check that the graph was created
        self.assertIsNotNone(graph)
        self.assertGreater(len(graph), 0)
    
    def test_output_rdf_report(self):
        """Test outputting an RDF report to a file."""
        # Create a simple RDF graph
        import rdflib
        graph = rdflib.Graph()
        graph.add((rdflib.URIRef('http://example.org/subject'),
                  rdflib.URIRef('http://example.org/predicate'),
                  rdflib.URIRef('http://example.org/object')))
        
        # Output to file
        output_file = os.path.join(self.output_dir, 'test.ttl')
        rdf_report.output_rdf_report(graph, output_file, format='turtle')
        
        # Check that the file was created
        self.assertTrue(os.path.exists(output_file))
        
        # Check file content
        with open(output_file, 'r') as f:
            content = f.read()
            self.assertIn('http://example.org/subject', content)
            self.assertIn('http://example.org/predicate', content)
            self.assertIn('http://example.org/object', content)
    
    def test_get_shex_schema(self):
        """Test getting the ShEx schema."""
        schema = rdf_report.get_shex_schema()
        
        # Check that the schema was returned
        self.assertIsNotNone(schema)
        self.assertIsInstance(schema, str)
        self.assertGreater(len(schema), 0)
        
        # Check for expected content
        self.assertIn('PREFIX', schema)
        self.assertIn('Variant', schema)
        self.assertIn('Feature', schema)

if __name__ == '__main__':
    unittest.main()
