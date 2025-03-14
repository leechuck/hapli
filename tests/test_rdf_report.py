#!/usr/bin/env python3
"""
Tests for the RDF reporting functionality.
"""

import os
import tempfile
import unittest
from unittest.mock import patch, MagicMock
from hapli.reporting import rdf_report

# Check if RDFLib is available
try:
    import rdflib
    RDFLIB_AVAILABLE = True
except ImportError:
    RDFLIB_AVAILABLE = False

@unittest.skipIf(not RDFLIB_AVAILABLE, "RDFLib not available")
class RDFReportTests(unittest.TestCase):
    """Test cases for RDF reporting functionality."""
    
    def setUp(self):
        """Set up test environment."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = self.temp_dir.name
    
    def tearDown(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()
    
    @patch('hapli.reporting.rdf_report.rdflib')
    def test_create_rdf_variant_report(self, mock_rdflib):
        """Test creating an RDF variant report."""
        # Setup mock graph
        mock_graph = MagicMock()
        mock_rdflib.Graph.return_value = mock_graph
        mock_graph.__len__.return_value = 10  # Simulate some triples in the graph
        
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
        self.assertEqual(len(graph), 10)
        
        # Verify that the appropriate methods were called
        mock_rdflib.Graph.assert_called_once()
    
    @patch('hapli.reporting.rdf_report.rdflib')
    def test_output_rdf_report(self, mock_rdflib):
        """Test outputting an RDF report to a file."""
        # Create a mock graph
        mock_graph = MagicMock()
        
        # Setup mock serialization
        def mock_serialize(destination=None, format='turtle'):
            if destination:
                with open(destination, 'w') as f:
                    f.write('<http://example.org/subject> <http://example.org/predicate> <http://example.org/object> .')
            return '<http://example.org/subject> <http://example.org/predicate> <http://example.org/object> .'
        
        mock_graph.serialize.side_effect = mock_serialize
        
        # Output to file
        output_file = os.path.join(self.output_dir, 'test.ttl')
        rdf_report.output_rdf_report(mock_graph, output_file, format='turtle')
        
        # Check that the file was created
        self.assertTrue(os.path.exists(output_file))
        
        # Check file content
        with open(output_file, 'r') as f:
            content = f.read()
            self.assertIn('http://example.org/subject', content)
            self.assertIn('http://example.org/predicate', content)
            self.assertIn('http://example.org/object', content)
        
        # Verify that serialize was called with the right parameters
        mock_graph.serialize.assert_called_once_with(destination=output_file, format='turtle')
    
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
        
    @patch('hapli.reporting.rdf_report.rdflib')
    def test_consolidated_rdf_report(self, mock_rdflib):
        """Test creating a consolidated RDF report."""
        # Setup mock graph
        mock_graph = MagicMock()
        mock_rdflib.Graph.return_value = mock_graph
        mock_graph.__len__.return_value = 15  # Simulate some triples in the graph
        
        # Test data
        variants = [{'id': 'var1', 'type': 'SNP', 'pos': 100, 'ref': 'A', 'alt': 'G'}]
        features = [{'id': 'gene1', 'type': 'gene', 'start': 50, 'end': 150}]
        samples = {'sample1': {'haplotypes': ['hap1', 'hap2']}}
        haplotypes = {'hap1': {'variants': ['var1']}}
        paths = {'REF': [], 'ALT': []}
        segments = {'seg1': 'ACGT'}
        
        # Create consolidated report
        graph = rdf_report.create_consolidated_rdf_report(
            variants, features, samples, haplotypes, paths, segments
        )
        
        # Check that the graph was created
        self.assertIsNotNone(graph)
        self.assertEqual(len(graph), 15)
        
        # Verify that the appropriate methods were called
        mock_rdflib.Graph.assert_called_once()

if __name__ == '__main__':
    unittest.main()
