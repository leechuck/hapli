"""
Core analysis functionality for variant effect reporting.
"""

import logging
import time
from typing import Dict, List, Optional, Tuple, Any

from variant_effect_report.core.models import Feature, Variant, FeatureEffect
from variant_effect_report.core.utils import build_path_sequence
from variant_effect_report.io.gfa_parser import GFAParser
from variant_effect_report.io.gff_parser import GFFParser
from variant_effect_report.io.report_writer import TextReportWriter, RDFReportWriter
from variant_effect_report.analysis.haplotype import HaplotypeAnalyzer
from variant_effect_report.analysis.zygosity import ZygosityAnalyzer


class VariantEffectAnalyzer:
    """Main analysis class for variant effect reporting."""
    
    def __init__(self, use_alignment: bool = False):
        """
        Initialize the variant effect analyzer.
        
        Args:
            use_alignment: Whether to use sequence alignment for analysis
        """
        self.use_alignment = use_alignment
        self.gfa_parser = GFAParser()
        self.gff_parser = GFFParser()
        self.haplotype_analyzer = HaplotypeAnalyzer(use_alignment=use_alignment)
        self.zygosity_analyzer = ZygosityAnalyzer()
        
        # Data storage
        self.segments = {}
        self.links = []
        self.paths = {}
        self.variants = []
        self.samples = {}
        self.haplotypes = {}
        self.features = []
        self.feature_by_id = {}
        self.children_by_parent = {}
    
    def load_data(self, gfa_file: str, gff_file: str) -> None:
        """
        Load data from GFA and GFF files.
        
        Args:
            gfa_file: Path to the GFA file
            gff_file: Path to the GFF file
        """
        # Parse GFA file
        start_time = time.time()
        logging.info(f"Loading data from {gfa_file} and {gff_file}")
        
        self.segments, self.links, self.paths, self.variants, self.samples, self.haplotypes = self.gfa_parser.parse(gfa_file)
        
        # Parse GFF file
        self.features, self.feature_by_id, self.children_by_parent = self.gff_parser.parse(gff_file)
        
        elapsed = time.time() - start_time
        logging.info(f"Data loaded in {elapsed:.2f}s")
    
    def analyze_sample(self, sample_name: str) -> Dict:
        """
        Analyze a sample to determine variant effects.
        
        Args:
            sample_name: Name of the sample to analyze
            
        Returns:
            Dictionary with analysis results
        """
        logging.info(f"Analyzing sample: {sample_name}")
        
        # Check if sample exists
        if sample_name not in self.samples:
            logging.error(f"Sample {sample_name} not found")
            return {'error': f"Sample {sample_name} not found"}
        
        # Get reference path
        ref_path_name = 'REF'
        if ref_path_name not in self.paths:
            logging.error("REF path not found in GFA")
            return {'error': "REF path not found in GFA"}
        
        # Build reference sequence
        ref_path_segments = self.paths[ref_path_name]['segments']
        ref_seq, _ = build_path_sequence(self.segments, ref_path_segments)
        
        # Collect variant segments
        variant_segments = {}
        for seg_id, seg_data in self.segments.items():
            if seg_data.get('variant_id'):
                variant_segments[seg_id] = seg_data
        
        # Determine if we have phased haplotypes for this sample
        sample_haplotypes = self.haplotypes.get(sample_name, {})
        has_phased_haplotypes = len(sample_haplotypes) >= 2
        
        if has_phased_haplotypes:
            logging.info(f"Found {len(sample_haplotypes)} phased haplotypes for sample {sample_name}")
            
            # Process each haplotype
            feature_effects_by_haplotype = {}
            
            for hap_name, path_name in sample_haplotypes.items():
                logging.info(f"Processing haplotype: {hap_name} (path: {path_name})")
                
                # Build path sequence
                path_segments = self.paths[path_name]['segments']
                hap_seq, segment_offsets = build_path_sequence(self.segments, path_segments)
                
                # Analyze haplotype differences
                hap_effects = self.haplotype_analyzer.analyze_haplotype_differences(
                    ref_seq, 
                    hap_seq, 
                    self.features, 
                    segment_offsets, 
                    variant_segments,
                    self.segments
                )
                
                feature_effects_by_haplotype[hap_name] = hap_effects
            
            # Analyze homozygous vs heterozygous effects
            zygosity_effects = self.zygosity_analyzer.analyze_sample_haplotypes(
                feature_effects_by_haplotype, 
                self.features,
                self.feature_by_id,
                self.children_by_parent
            )
            
            return {
                'sample_name': sample_name,
                'has_phased_haplotypes': True,
                'haplotypes': list(sample_haplotypes.keys()),
                'feature_effects_by_haplotype': feature_effects_by_haplotype,
                'zygosity_effects': zygosity_effects
            }
        else:
            # Process single haplotype
            logging.info(f"Only one haplotype found for sample {sample_name}, processing as single path")
            
            # Get the first path for this sample
            path_name = self.samples[sample_name][0]
            path_segments = self.paths[path_name]['segments']
            
            # Build path sequence
            alt_seq, segment_offsets = build_path_sequence(self.segments, path_segments)
            
            # Analyze differences
            feature_effects = self.haplotype_analyzer.analyze_haplotype_differences(
                ref_seq, 
                alt_seq, 
                self.features, 
                segment_offsets, 
                variant_segments,
                self.segments
            )
            
            return {
                'sample_name': sample_name,
                'has_phased_haplotypes': False,
                'feature_effects': feature_effects
            }
    
    def generate_report(self, sample_name: Optional[str] = None, output: Optional[str] = None, 
                       format: str = 'text', rdf_format: str = 'turtle', 
                       base_uri: str = "http://example.org/genomics/") -> None:
        """
        Generate a report for a sample or all samples.
        
        Args:
            sample_name: Name of the sample to analyze (None for all samples)
            output: Output file path (None for stdout)
            format: Output format ('text' or 'rdf')
            rdf_format: RDF serialization format
            base_uri: Base URI for RDF output
        """
        if sample_name:
            # Analyze a single sample
            result = self.analyze_sample(sample_name)
            
            if 'error' in result:
                logging.error(f"Error analyzing sample {sample_name}: {result['error']}")
                return
            
            # Create appropriate report writer
            if format.lower() == 'text':
                writer = TextReportWriter(output=output)
            else:  # RDF format
                writer = RDFReportWriter(output=output, format=rdf_format, base_uri=base_uri)
            
            # Generate report
            if result['has_phased_haplotypes']:
                writer.write_report(
                    self.variants,
                    self.features,
                    feature_effects=None,
                    sample_name=sample_name,
                    sample_effects=result['zygosity_effects']
                )
            else:
                writer.write_report(
                    self.variants,
                    self.features,
                    feature_effects=result['feature_effects'],
                    sample_name=sample_name
                )
        else:
            # Process all samples
            for sample_name in self.samples:
                output_file = f"{output}_{sample_name}" if output else None
                self.generate_report(sample_name, output_file, format, rdf_format, base_uri)
