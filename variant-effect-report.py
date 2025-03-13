#!/usr/bin/env python3
"""
GFA Variant Effect Analyzer

Analyzes variants in GFA files with sample haplotype information and reports their effects 
on genomic features defined in GFF3 files.

Usage:
    python gfa_variant_analyzer.py out.gfa test.gff3 --output variant_report.txt
"""

import argparse
import sys
import os
import logging
import time
import re
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2
from Bio import Align

import json
import rdflib
from rdflib import Graph, Literal, BNode, Namespace, RDF, URIRef, XSD
from rdflib.namespace import FOAF, RDFS, SKOS
import uuid
from datetime import datetime


# Define namespaces for genomic data and sequence features
class GenomicNamespaces:
    """Provides standard namespaces for genomic data in RDF."""
    
    def __init__(self, base_uri="http://example.org/genomics/"):
        # Base namespace for our data
        self.base = Namespace(base_uri)
        
        # Standard ontologies
        self.rdf = RDF
        self.rdfs = RDFS
        self.xsd = XSD
        self.skos = SKOS
        
        # Genomic ontologies
        self.faldo = Namespace("http://biohackathon.org/resource/faldo#")
        self.so = Namespace("http://purl.obolibrary.org/obo/SO_")  # Sequence Ontology
        self.obo = Namespace("http://purl.obolibrary.org/obo/")
        self.dc = Namespace("http://purl.org/dc/terms/")
        self.sio = Namespace("http://semanticscience.org/resource/")  # Semantic Science Integrated Ontology
        self.gff = Namespace("http://www.sequenceontology.org/gff3/1.0/")  # GFF3 terms
        
        # VCF/Variant specific terms
        self.variant = Namespace(f"{base_uri}variant/")
        self.sample = Namespace(f"{base_uri}sample/")
        self.feature = Namespace(f"{base_uri}feature/")
        self.effect = Namespace(f"{base_uri}effect/")
        self.position = Namespace(f"{base_uri}position/")
        
    def bind_to_graph(self, g):
        """Bind all namespaces to a graph."""
        g.bind("rdf", self.rdf)
        g.bind("rdfs", self.rdfs)
        g.bind("xsd", self.xsd)
        g.bind("skos", self.skos)
        g.bind("faldo", self.faldo)
        g.bind("so", self.so)
        g.bind("obo", self.obo)
        g.bind("dc", self.dc)
        g.bind("sio", self.sio)
        g.bind("gff", self.gff)
        g.bind("variant", self.variant)
        g.bind("sample", self.sample)
        g.bind("feature", self.feature)
        g.bind("effect", self.effect)
        g.bind("position", self.position)
        g.bind("base", self.base)
        
    def get_so_term(self, variant_type):
        """Map variant types to Sequence Ontology terms."""
        # Map common variant types to SO terms
        so_mapping = {
            "SNP": "0001483",      # SNV - single nucleotide variant
            "SNV": "0001483",      # SNV - single nucleotide variant
            "DEL": "0000159",      # deletion
            "INS": "0000667",      # insertion
            "INV": "1000036",      # inversion
            "DUP": "1000035",      # duplication
            "UNKNOWN": "0001059"   # sequence_alteration
        }
        
        # Return the mapped SO term or the default term for unknown types
        return self.so[so_mapping.get(variant_type, "0001059")]

    def get_effect_term(self, effect_type):
        """Map effect types to Sequence Ontology terms."""
        # Map effect types to SO terms
        effect_mapping = {
            "no_change": "0000992",            # silent_mutation
            "amino_acid_change": "0001583",    # missense_variant
            "frame_shift": "0001589",          # frameshift_variant
            "in_frame_change": "0001650",      # inframe_variant
            "premature_stop_codon": "0001587", # stop_gained
            "start_codon_disruption": "0002012", # start_lost
            "promoter_affected": "0001566",    # regulatory_region_variant
            "terminator_affected": "0001566",  # regulatory_region_variant
            "splicing_affected": "0001568",    # splicing_variant
            "length_change": "0001059",        # sequence_alteration
            "sequence_change": "0001059",      # sequence_alteration
            "compound_heterozygous": "0001034" # compound_heterozygous
        }
        
        # Return the mapped SO term or a generic term for unknown effects
        return self.so[effect_mapping.get(effect_type, "0001564")]  # gene_variant as default


def create_rdf_variant_report(variants, features, feature_effects=None, sample_name=None, sample_effects=None, base_uri="http://example.org/genomics/"):
    """
    Create an RDF graph of variant effects.
    
    Args:
        variants: List of variant dictionaries
        features: List of feature dictionaries
        feature_effects: List of feature effect dictionaries (optional)
        sample_name: Name of the sample being analyzed (optional)
        sample_effects: Dictionary of sample-specific effects (optional)
        base_uri: Base URI for the RDF graph
        
    Returns:
        RDF graph object
    """
    # Initialize graph and namespaces
    g = Graph()
    ns = GenomicNamespaces(base_uri)
    ns.bind_to_graph(g)
    
    # Create report node
    report_id = str(uuid.uuid4())
    report_uri = URIRef(f"{base_uri}report/{report_id}")
    g.add((report_uri, RDF.type, ns.sio["SequenceVariantAnalysisReport"]))
    g.add((report_uri, ns.dc.created, Literal(datetime.now().isoformat(), datatype=XSD.dateTime)))
    
    # Add sample information if available
    if sample_name:
        sample_uri = URIRef(f"{ns.sample}{sample_name.replace(' ', '_')}")
        g.add((sample_uri, RDF.type, ns.sio.Sample))
        g.add((sample_uri, RDFS.label, Literal(sample_name)))
        g.add((report_uri, ns.sio.refersTo, sample_uri))
    
    # Process variants
    for variant in variants:
        var_id = variant.get('id', f"unknown_{str(uuid.uuid4())}")
        var_type = variant.get('type', 'UNKNOWN')
        
        # Create variant node
        var_uri = URIRef(f"{ns.variant}{var_id}")
        g.add((var_uri, RDF.type, ns.so.sequence_variant))
        g.add((var_uri, RDFS.label, Literal(var_id)))
        g.add((var_uri, ns.variant.variantType, ns.get_so_term(var_type)))
        g.add((var_uri, ns.variant.lengthChange, Literal(variant.get('length_change', 0), datatype=XSD.integer)))
        
        # Add variant to report
        g.add((report_uri, ns.variant.hasVariant, var_uri))
        
        # Add position information if available
        if 'pos' in variant and variant.get('pos', 0) > 0:
            pos_uri = URIRef(f"{ns.position}{var_id}")
            g.add((pos_uri, RDF.type, ns.faldo.Region))
            g.add((var_uri, ns.faldo.location, pos_uri))
            
            # Add start position
            start_uri = URIRef(f"{ns.position}{var_id}_start")
            g.add((start_uri, RDF.type, ns.faldo.ExactPosition))
            g.add((start_uri, ns.faldo.position, Literal(variant.get('pos', 0), datatype=XSD.integer)))
            g.add((pos_uri, ns.faldo.begin, start_uri))
            
            # Add end position
            end_uri = URIRef(f"{ns.position}{var_id}_end")
            g.add((end_uri, RDF.type, ns.faldo.ExactPosition))
            g.add((end_uri, ns.faldo.position, Literal(variant.get('end', variant.get('pos', 0)), datatype=XSD.integer)))
            g.add((pos_uri, ns.faldo.end, end_uri))
        
        # Add affected segments
        if 'segments' in variant and variant['segments']:
            for seg_id in variant['segments']:
                seg_uri = URIRef(f"{ns.base}segment/{seg_id}")
                g.add((seg_uri, RDF.type, ns.base.Segment))
                g.add((seg_uri, RDFS.label, Literal(seg_id)))
                g.add((var_uri, ns.variant.affectsSegment, seg_uri))
    
    # Process features
    for feature in features:
        feature_id = feature['attributes'].get('ID', f"unknown_{str(uuid.uuid4())}")
        feature_name = feature['attributes'].get('Name', feature_id)
        feature_type = feature['type']
        
        # Create feature node
        feature_uri = URIRef(f"{ns.feature}{feature_id}")
        g.add((feature_uri, RDF.type, URIRef(f"{ns.gff}{feature_type}")))
        g.add((feature_uri, RDFS.label, Literal(feature_name)))
        
        # Add feature to report
        g.add((report_uri, ns.feature.hasFeature, feature_uri))
        
        # Add location information
        loc_uri = URIRef(f"{ns.position}feature_{feature_id}")
        g.add((loc_uri, RDF.type, ns.faldo.Region))
        g.add((feature_uri, ns.faldo.location, loc_uri))
        
        # Add start position
        start_uri = URIRef(f"{ns.position}feature_{feature_id}_start")
        g.add((start_uri, RDF.type, ns.faldo.ExactPosition))
        g.add((start_uri, ns.faldo.position, Literal(feature['start'], datatype=XSD.integer)))
        g.add((loc_uri, ns.faldo.begin, start_uri))
        
        # Add end position
        end_uri = URIRef(f"{ns.position}feature_{feature_id}_end")
        g.add((end_uri, RDF.type, ns.faldo.ExactPosition))
        g.add((end_uri, ns.faldo.position, Literal(feature['end'], datatype=XSD.integer)))
        g.add((loc_uri, ns.faldo.end, end_uri))
        
        # Add strand information
        if feature['strand'] == '+':
            g.add((feature_uri, ns.faldo.strand, ns.faldo.ForwardStrand))
        elif feature['strand'] == '-':
            g.add((feature_uri, ns.faldo.strand, ns.faldo.ReverseStrand))
        else:
            g.add((feature_uri, ns.faldo.strand, ns.faldo.BothStrand))
    
    # Process effects
    # Use sample_effects if provided, otherwise use feature_effects
    effects_to_process = []
    
    if sample_name and sample_effects and not sample_effects.get('incomplete'):
        effects_to_process.extend(sample_effects.get('homozygous', []))
        effects_to_process.extend(sample_effects.get('heterozygous', []))
        
        # Add compound heterozygous effects
        if 'gene_compound_heterozygous' in sample_effects:
            for comp_het in sample_effects['gene_compound_heterozygous']:
                # Create a synthetic effect entry
                comp_effect = {
                    'feature': comp_het['feature'],
                    'feature_type': 'gene',
                    'effects': ['compound_heterozygous'],
                    'variants': comp_het['variants'],
                    'zygosity': 'compound_heterozygous',
                    'details': comp_het['details']
                }
                effects_to_process.append(comp_effect)
                
        # Add feature-level compound heterozygous effects
        if 'feature_compound_heterozygous' in sample_effects:
            for comp_het in sample_effects['feature_compound_heterozygous']:
                # Create a synthetic effect entry
                comp_effect = {
                    'feature': comp_het['feature'],
                    'feature_type': comp_het['feature_type'],
                    'effects': ['compound_heterozygous'],
                    'variants': comp_het['variants'],
                    'zygosity': 'compound_heterozygous',
                    'details': comp_het['details']
                }
                effects_to_process.append(comp_effect)
    else:
        effects_to_process = feature_effects if feature_effects else []
    
    # Process each effect
    for effect in effects_to_process:
        if not effect.get('effects') or effect['effects'] == ['no_change']:
            continue  # Skip effects with no change
            
        feature = effect['feature']
        feature_id = feature['attributes'].get('ID', 'unknown')
        feature_type = feature['type']
        feature_uri = URIRef(f"{ns.feature}{feature_id}")
        
        # For each effect type
        for effect_type in effect['effects']:
            if effect_type == 'no_change':
                continue
                
            effect_id = f"{feature_id}_{effect_type}_{str(uuid.uuid4())}"
            effect_uri = URIRef(f"{ns.effect}{effect_id}")
            
            # Create effect node
            g.add((effect_uri, RDF.type, ns.effect.VariantEffect))
            g.add((effect_uri, ns.effect.effectType, ns.get_effect_term(effect_type)))
            g.add((effect_uri, ns.effect.affectsFeature, feature_uri))
            
            # Add effect to report
            g.add((report_uri, ns.effect.hasEffect, effect_uri))
            
            # Add detailed description if available
            if effect_type in effect.get('details', {}):
                detail = effect['details'][effect_type]
                if isinstance(detail, dict):
                    desc = ", ".join(f"{k}={v}" for k, v in detail.items())
                elif isinstance(detail, list):
                    desc = ", ".join(str(d) for d in detail)
                else:
                    desc = str(detail)
                g.add((effect_uri, ns.dc.description, Literal(desc)))
            
            # Add zygosity if available
            if 'zygosity' in effect:
                g.add((effect_uri, ns.variant.zygosity, Literal(effect['zygosity'])))
            
            # Add sequence information
            if 'ref_feature_seq' in effect and 'alt_feature_seq' in effect:
                g.add((effect_uri, ns.variant.referenceSequence, Literal(effect['ref_feature_seq'])))
                g.add((effect_uri, ns.variant.alternateSequence, Literal(effect['alt_feature_seq'])))
            
            # Add variants causing this effect
            if 'variants' in effect and effect['variants']:
                for var in effect['variants']:
                    var_id = var.get('id', 'unknown')
                    var_uri = URIRef(f"{ns.variant}{var_id}")
                    g.add((effect_uri, ns.effect.causedBy, var_uri))
    
    return g


def output_rdf_report(graph, output=None, format='turtle'):
    """
    Output an RDF graph in the specified format.
    
    Args:
        graph: RDF graph object
        output: Output file (default: stdout)
        format: RDF serialization format (default: turtle)
        
    Returns:
        None
    """
    # Map format names to rdflib serialization format names
    format_map = {
        'turtle': 'turtle',
        'ttl': 'turtle',
        'n3': 'n3',
        'xml': 'xml',
        'rdf': 'xml',
        'rdfxml': 'xml',
        'jsonld': 'json-ld',
        'json-ld': 'json-ld',
        'nt': 'nt',
        'ntriples': 'nt'
    }
    
    # Get the serialization format
    rdf_format = format_map.get(format.lower(), 'turtle')
    
    # Output to file or stdout
    if output:
        graph.serialize(destination=output, format=rdf_format)
    else:
        # Serialize to string and print to stdout
        output_str = graph.serialize(format=rdf_format)
        sys.stdout.write(output_str.decode('utf-8') if isinstance(output_str, bytes) else output_str)


def create_consolidated_rdf_report(variants, features, samples, haplotypes, paths, segments, 
                                  feature_by_id, children_by_parent, base_uri="http://example.org/genomics/"):
    """
    Create a consolidated RDF graph with variant effects for all samples.
    
    Args:
        variants: List of variant dictionaries
        features: List of feature dictionaries
        samples: Dictionary of sample names and their paths
        haplotypes: Dictionary of samples and their haplotypes
        paths: Dictionary of path names and their data
        segments: Dictionary of segment names and their data
        feature_by_id: Dictionary of features indexed by their ID
        children_by_parent: Dictionary of child features indexed by parent ID
        base_uri: Base URI for the RDF graph
        
    Returns:
        RDF graph object with all sample data
    """
    # Initialize graph and namespaces
    g = Graph()
    ns = GenomicNamespaces(base_uri)
    ns.bind_to_graph(g)
    
    # Create report node
    report_id = str(uuid.uuid4())
    report_uri = URIRef(f"{base_uri}report/{report_id}")
    g.add((report_uri, RDF.type, ns.sio["SequenceVariantAnalysisReport"]))
    g.add((report_uri, ns.dc.created, Literal(datetime.now().isoformat(), datatype=XSD.dateTime)))
    g.add((report_uri, RDFS.label, Literal("Consolidated Variant Effect Report")))
    
    # Get reference path
    ref_path_name = 'REF'
    if ref_path_name in paths:
        ref_path_segments = paths[ref_path_name]['segments']
        ref_seq, _ = build_path_sequence(segments, ref_path_segments)
    else:
        logging.error("REF path not found in GFA")
        return g
    
    # Collect variant segments
    variant_segments = {}
    for seg_id, seg_data in segments.items():
        if seg_data.get('variant_id'):
            variant_segments[seg_id] = seg_data
    
    # Process variants (just once for all samples)
    processed_variants = {}
    for variant in variants:
        var_id = variant.get('id', f"unknown_{str(uuid.uuid4())}")
        var_type = variant.get('type', 'UNKNOWN')
        
        # Create variant node
        var_uri = URIRef(f"{ns.variant}{var_id}")
        g.add((var_uri, RDF.type, ns.so.sequence_variant))
        g.add((var_uri, RDFS.label, Literal(var_id)))
        g.add((var_uri, ns.variant.variantType, ns.get_so_term(var_type)))
        g.add((var_uri, ns.variant.lengthChange, Literal(variant.get('length_change', 0), datatype=XSD.integer)))
        
        # Add variant to report
        g.add((report_uri, ns.variant.hasVariant, var_uri))
        
        # Add position information if available
        if 'pos' in variant and variant.get('pos', 0) > 0:
            pos_uri = URIRef(f"{ns.position}{var_id}")
            g.add((pos_uri, RDF.type, ns.faldo.Region))
            g.add((var_uri, ns.faldo.location, pos_uri))
            
            # Add start position
            start_uri = URIRef(f"{ns.position}{var_id}_start")
            g.add((start_uri, RDF.type, ns.faldo.ExactPosition))
            g.add((start_uri, ns.faldo.position, Literal(variant.get('pos', 0), datatype=XSD.integer)))
            g.add((pos_uri, ns.faldo.begin, start_uri))
            
            # Add end position
            end_uri = URIRef(f"{ns.position}{var_id}_end")
            g.add((end_uri, RDF.type, ns.faldo.ExactPosition))
            g.add((end_uri, ns.faldo.position, Literal(variant.get('end', variant.get('pos', 0)), datatype=XSD.integer)))
            g.add((pos_uri, ns.faldo.end, end_uri))
        
        # Add affected segments
        if 'segments' in variant and variant['segments']:
            for seg_id in variant['segments']:
                seg_uri = URIRef(f"{ns.base}segment/{seg_id}")
                g.add((seg_uri, RDF.type, ns.base.Segment))
                g.add((seg_uri, RDFS.label, Literal(seg_id)))
                g.add((var_uri, ns.variant.affectsSegment, seg_uri))
        
        # Store the variant URI for later reference
        processed_variants[var_id] = var_uri
        
    # Process features (just once for all samples)
    processed_features = {}
    for feature in features:
        feature_id = feature['attributes'].get('ID', f"unknown_{str(uuid.uuid4())}")
        feature_name = feature['attributes'].get('Name', feature_id)
        feature_type = feature['type']
        
        # Create feature node
        feature_uri = URIRef(f"{ns.feature}{feature_id}")
        g.add((feature_uri, RDF.type, URIRef(f"{ns.gff}{feature_type}")))
        g.add((feature_uri, RDFS.label, Literal(feature_name)))
        
        # Add feature to report
        g.add((report_uri, ns.feature.hasFeature, feature_uri))
        
        # Add location information
        loc_uri = URIRef(f"{ns.position}feature_{feature_id}")
        g.add((loc_uri, RDF.type, ns.faldo.Region))
        g.add((feature_uri, ns.faldo.location, loc_uri))
        
        # Add start position
        start_uri = URIRef(f"{ns.position}feature_{feature_id}_start")
        g.add((start_uri, RDF.type, ns.faldo.ExactPosition))
        g.add((start_uri, ns.faldo.position, Literal(feature['start'], datatype=XSD.integer)))
        g.add((loc_uri, ns.faldo.begin, start_uri))
        
        # Add end position
        end_uri = URIRef(f"{ns.position}feature_{feature_id}_end")
        g.add((end_uri, RDF.type, ns.faldo.ExactPosition))
        g.add((end_uri, ns.faldo.position, Literal(feature['end'], datatype=XSD.integer)))
        g.add((loc_uri, ns.faldo.end, end_uri))
        
        # Add strand information
        if feature['strand'] == '+':
            g.add((feature_uri, ns.faldo.strand, ns.faldo.ForwardStrand))
        elif feature['strand'] == '-':
            g.add((feature_uri, ns.faldo.strand, ns.faldo.ReverseStrand))
        else:
            g.add((feature_uri, ns.faldo.strand, ns.faldo.BothStrand))
            
        # Store the feature URI for later reference
        processed_features[feature_id] = feature_uri
    
    # Process each sample
    for sample_name, sample_paths in samples.items():
        logging.info(f"Adding sample {sample_name} to consolidated report")
        
        # Create sample node
        sample_uri = URIRef(f"{ns.sample}{sample_name.replace(' ', '_')}")
        g.add((sample_uri, RDF.type, ns.sio.Sample))
        g.add((sample_uri, RDFS.label, Literal(sample_name)))
        g.add((report_uri, ns.sio.refersTo, sample_uri))
        
        # Determine if we have phased haplotypes for this sample
        sample_haplotypes = haplotypes.get(sample_name, {})
        has_phased_haplotypes = len(sample_haplotypes) >= 2
        
        if has_phased_haplotypes:
            # Process each haplotype
            feature_effects_by_haplotype = {}
            
            for hap_name, path_name in sample_haplotypes.items():
                # Create haplotype node
                hap_uri = URIRef(f"{ns.sample}{sample_name.replace(' ', '_')}/haplotype/{hap_name.replace(' ', '_')}")
                g.add((hap_uri, RDF.type, ns.base.Haplotype))
                g.add((hap_uri, RDFS.label, Literal(f"{sample_name} - {hap_name}")))
                g.add((sample_uri, ns.base.hasHaplotype, hap_uri))
                
                # Build path sequence
                path_segments = paths[path_name]['segments']
                hap_seq, segment_offsets = build_path_sequence(segments, path_segments)
                
                # Analyze differences
                hap_effects = analyze_haplotype_differences(
                    ref_seq, 
                    hap_seq, 
                    features, 
                    segment_offsets, 
                    variant_segments,
                    segments
                )
                
                feature_effects_by_haplotype[hap_name] = hap_effects
            
            # Analyze homozygous vs heterozygous effects
            zygosity_effects = analyze_sample_haplotypes(
                feature_effects_by_haplotype, 
                features,
                feature_by_id,
                children_by_parent
            )
            
            # Add effects to graph
            effect_counter = 0
            
            # Process homozygous effects
            for effect in zygosity_effects.get('homozygous', []):
                if not effect.get('effects') or effect['effects'] == ['no_change']:
                    continue  # Skip effects with no change
                    
                feature = effect['feature']
                feature_id = feature['attributes'].get('ID', 'unknown')
                feature_uri = processed_features.get(feature_id)
                
                if not feature_uri:
                    continue
                
                # For each effect type
                for effect_type in effect['effects']:
                    if effect_type == 'no_change':
                        continue
                    
                    effect_counter += 1
                    effect_id = f"{sample_name}_{feature_id}_{effect_type}_{effect_counter}"
                    effect_uri = URIRef(f"{ns.effect}{effect_id}")
                    
                    # Create effect node
                    g.add((effect_uri, RDF.type, ns.effect.VariantEffect))
                    g.add((effect_uri, ns.effect.effectType, ns.get_effect_term(effect_type)))
                    g.add((effect_uri, ns.effect.affectsFeature, feature_uri))
                    g.add((effect_uri, ns.variant.zygosity, Literal("homozygous")))
                    g.add((effect_uri, ns.effect.inSample, sample_uri))
                    
                    # Add effect to report
                    g.add((report_uri, ns.effect.hasEffect, effect_uri))
                    
                    # Add detailed description if available
                    if effect_type in effect.get('details', {}):
                        detail = effect['details'][effect_type]
                        if isinstance(detail, dict):
                            desc = ", ".join(f"{k}={v}" for k, v in detail.items())
                        elif isinstance(detail, list):
                            desc = ", ".join(str(d) for d in detail)
                        else:
                            desc = str(detail)
                        g.add((effect_uri, ns.dc.description, Literal(desc)))
                    
                    # Add sequence information
                    if 'ref_feature_seq' in effect and 'alt_feature_seq' in effect:
                        g.add((effect_uri, ns.variant.referenceSequence, Literal(effect['ref_feature_seq'])))
                        g.add((effect_uri, ns.variant.alternateSequence, Literal(effect['alt_feature_seq'])))
                    
                    # Add variants causing this effect
                    if 'variants' in effect and effect['variants']:
                        for var in effect['variants']:
                            var_id = var.get('id', 'unknown')
                            var_uri = processed_variants.get(var_id)
                            if var_uri:
                                g.add((effect_uri, ns.effect.causedBy, var_uri))
            
            # Process heterozygous effects
            for effect in zygosity_effects.get('heterozygous', []):
                if not effect.get('effects') or effect['effects'] == ['no_change']:
                    continue  # Skip effects with no change
                    
                feature = effect['feature']
                feature_id = feature['attributes'].get('ID', 'unknown')
                feature_uri = processed_features.get(feature_id)
                
                if not feature_uri:
                    continue
                
                # For each effect type
                effect_list = effect.get('combined_effects', effect['effects'])
                for effect_type in effect_list:
                    if effect_type == 'no_change':
                        continue
                    
                    effect_counter += 1
                    effect_id = f"{sample_name}_{feature_id}_{effect_type}_{effect_counter}"
                    effect_uri = URIRef(f"{ns.effect}{effect_id}")
                    
                    # Create effect node
                    g.add((effect_uri, RDF.type, ns.effect.VariantEffect))
                    g.add((effect_uri, ns.effect.effectType, ns.get_effect_term(effect_type)))
                    g.add((effect_uri, ns.effect.affectsFeature, feature_uri))
                    g.add((effect_uri, ns.variant.zygosity, Literal("heterozygous")))
                    g.add((effect_uri, ns.effect.inSample, sample_uri))
                    
                    # Add effect to report
                    g.add((report_uri, ns.effect.hasEffect, effect_uri))
                    
                    # Add detailed description if available
                    if effect_type in effect.get('details', {}):
                        detail = effect['details'][effect_type]
                        if isinstance(detail, dict):
                            desc = ", ".join(f"{k}={v}" for k, v in detail.items())
                        elif isinstance(detail, list):
                            desc = ", ".join(str(d) for d in detail)
                        else:
                            desc = str(detail)
                        g.add((effect_uri, ns.dc.description, Literal(desc)))
                    
                    # Add sequence information
                    if 'ref_feature_seq' in effect and 'alt_feature_seq' in effect:
                        g.add((effect_uri, ns.variant.referenceSequence, Literal(effect['ref_feature_seq'])))
                        g.add((effect_uri, ns.variant.alternateSequence, Literal(effect['alt_feature_seq'])))
                    
                    # Add variants causing this effect
                    if 'variants' in effect and effect['variants']:
                        for var in effect['variants']:
                            var_id = var.get('id', 'unknown')
                            var_uri = processed_variants.get(var_id)
                            if var_uri:
                                g.add((effect_uri, ns.effect.causedBy, var_uri))
                    
                    # For heterozygous effects, add haplotype-specific information
                    if 'haplotype_effects' in effect:
                        for hap_name, hap_effect in effect['haplotype_effects'].items():
                            hap_uri = URIRef(f"{ns.sample}{sample_name.replace(' ', '_')}/haplotype/{hap_name.replace(' ', '_')}")
                            g.add((effect_uri, ns.effect.inHaplotype, hap_uri))
            
            # Process compound heterozygous effects
            if 'gene_compound_heterozygous' in zygosity_effects:
                for comp_het in zygosity_effects['gene_compound_heterozygous']:
                    gene_id = comp_het.get('gene_id', 'unknown')
                    gene_uri = processed_features.get(gene_id)
                    
                    if not gene_uri:
                        continue
                    
                    effect_counter += 1
                    effect_id = f"{sample_name}_{gene_id}_compound_heterozygous_{effect_counter}"
                    effect_uri = URIRef(f"{ns.effect}{effect_id}")
                    
                    # Create effect node
                    g.add((effect_uri, RDF.type, ns.effect.VariantEffect))
                    g.add((effect_uri, ns.effect.effectType, ns.get_effect_term('compound_heterozygous')))
                    g.add((effect_uri, ns.effect.affectsFeature, gene_uri))
                    g.add((effect_uri, ns.variant.zygosity, Literal("compound_heterozygous")))
                    g.add((effect_uri, ns.effect.inSample, sample_uri))
                    
                    # Add effect to report
                    g.add((report_uri, ns.effect.hasEffect, effect_uri))
                    
                    # Add description
                    g.add((effect_uri, ns.dc.description, Literal("Different variants affecting the same gene on different haplotypes")))
                    
                    # Add variants by haplotype
                    for hap_name, hap_variants in comp_het.get('haplotype_variants', {}).items():
                        hap_uri = URIRef(f"{ns.sample}{sample_name.replace(' ', '_')}/haplotype/{hap_name.replace(' ', '_')}")
                        g.add((effect_uri, ns.effect.inHaplotype, hap_uri))
                        
                        for var in hap_variants:
                            var_id = var.get('id', 'unknown')
                            var_uri = processed_variants.get(var_id)
                            if var_uri:
                                # Create a blank node for the haplotype-variant association
                                bnode = BNode()
                                g.add((bnode, RDF.type, ns.effect.HaplotypeVariantAssociation))
                                g.add((bnode, ns.effect.haplotype, hap_uri))
                                g.add((bnode, ns.effect.variant, var_uri))
                                g.add((effect_uri, ns.effect.variantInHaplotype, bnode))
            
            # Process feature-level compound heterozygous effects
            if 'feature_compound_heterozygous' in zygosity_effects:
                for comp_het in zygosity_effects['feature_compound_heterozygous']:
                    feature_id = comp_het.get('feature_id', 'unknown')
                    feature_uri = processed_features.get(feature_id)
                    
                    if not feature_uri:
                        continue
                    
                    effect_counter += 1
                    effect_id = f"{sample_name}_{feature_id}_compound_heterozygous_{effect_counter}"
                    effect_uri = URIRef(f"{ns.effect}{effect_id}")
                    
                    # Create effect node
                    g.add((effect_uri, RDF.type, ns.effect.VariantEffect))
                    g.add((effect_uri, ns.effect.effectType, ns.get_effect_term('compound_heterozygous')))
                    g.add((effect_uri, ns.effect.affectsFeature, feature_uri))
                    g.add((effect_uri, ns.variant.zygosity, Literal("compound_heterozygous")))
                    g.add((effect_uri, ns.effect.inSample, sample_uri))
                    
                    # Add effect to report
                    g.add((report_uri, ns.effect.hasEffect, effect_uri))
                    
                    # Add description
                    feature_type = comp_het.get('feature_type', 'feature')
                    g.add((effect_uri, ns.dc.description, Literal(f"Different variants affecting the same {feature_type} on different haplotypes")))
                    
                    # Add variants by haplotype
                    for hap_name, hap_variants in comp_het.get('haplotype_variants', {}).items():
                        hap_uri = URIRef(f"{ns.sample}{sample_name.replace(' ', '_')}/haplotype/{hap_name.replace(' ', '_')}")
                        g.add((effect_uri, ns.effect.inHaplotype, hap_uri))
                        
                        for var in hap_variants:
                            var_id = var.get('id', 'unknown')
                            var_uri = processed_variants.get(var_id)
                            if var_uri:
                                # Create a blank node for the haplotype-variant association
                                bnode = BNode()
                                g.add((bnode, RDF.type, ns.effect.HaplotypeVariantAssociation))
                                g.add((bnode, ns.effect.haplotype, hap_uri))
                                g.add((bnode, ns.effect.variant, var_uri))
                                g.add((effect_uri, ns.effect.variantInHaplotype, bnode))
                                
        else:
            # Process single haplotype
            path_name = sample_paths[0]
            path_segments = paths[path_name]['segments']
            
            # Build path sequence
            alt_seq, segment_offsets = build_path_sequence(segments, path_segments)
            
            # Analyze differences
            feature_effects = analyze_haplotype_differences(
                ref_seq, 
                alt_seq, 
                features, 
                segment_offsets, 
                variant_segments,
                segments
            )
            
            # Add effects to graph
            effect_counter = 0
            
            for effect in feature_effects:
                if not effect.get('effects') or effect['effects'] == ['no_change']:
                    continue  # Skip effects with no change
                    
                feature = effect['feature']
                feature_id = feature['attributes'].get('ID', 'unknown')
                feature_uri = processed_features.get(feature_id)
                
                if not feature_uri:
                    continue
                
                # For each effect type
                for effect_type in effect['effects']:
                    if effect_type == 'no_change':
                        continue
                    
                    effect_counter += 1
                    effect_id = f"{sample_name}_{feature_id}_{effect_type}_{effect_counter}"
                    effect_uri = URIRef(f"{ns.effect}{effect_id}")
                    
                    # Create effect node
                    g.add((effect_uri, RDF.type, ns.effect.VariantEffect))
                    g.add((effect_uri, ns.effect.effectType, ns.get_effect_term(effect_type)))
                    g.add((effect_uri, ns.effect.affectsFeature, feature_uri))
                    g.add((effect_uri, ns.effect.inSample, sample_uri))
                    
                    # Add effect to report
                    g.add((report_uri, ns.effect.hasEffect, effect_uri))
                    
                    # Add detailed description if available
                    if effect_type in effect.get('details', {}):
                        detail = effect['details'][effect_type]
                        if isinstance(detail, dict):
                            desc = ", ".join(f"{k}={v}" for k, v in detail.items())
                        elif isinstance(detail, list):
                            desc = ", ".join(str(d) for d in detail)
                        else:
                            desc = str(detail)
                        g.add((effect_uri, ns.dc.description, Literal(desc)))
                    
                    # Add sequence information
                    if 'ref_feature_seq' in effect and 'alt_feature_seq' in effect:
                        g.add((effect_uri, ns.variant.referenceSequence, Literal(effect['ref_feature_seq'])))
                        g.add((effect_uri, ns.variant.alternateSequence, Literal(effect['alt_feature_seq'])))
                    
                    # Add variants causing this effect
                    if 'variants' in effect and effect['variants']:
                        for var in effect['variants']:
                            var_id = var.get('id', 'unknown')
                            var_uri = processed_variants.get(var_id)
                            if var_uri:
                                g.add((effect_uri, ns.effect.causedBy, var_uri))
    
    return g

def analyze_sequence_alignment(ref_seq, alt_seq, gap_open=-10, gap_extend=-0.5, match=2, mismatch=-1):
    """
    Perform global sequence alignment to identify differences between reference and alternate sequences.
    
    Uses Bio.Align.PairwiseAligner instead of the deprecated Bio.pairwise2
    
    Args:
        ref_seq: Reference sequence string
        alt_seq: Alternate sequence string
        gap_open: Gap opening penalty (default: -10)
        gap_extend: Gap extension penalty (default: -0.5)
        match: Match score (default: 2)
        mismatch: Mismatch score (default: -1)
        
    Returns:
        Dictionary with alignment information and identified changes
    """
    # Create aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.open_gap_score = gap_open
    aligner.extend_gap_score = gap_extend
    
    # Perform alignment
    alignments = aligner.align(ref_seq, alt_seq)
    
    # Get the best alignment
    if not alignments:
        return {
            'alignment_score': 0,
            'changes': 'complete_replacement',
            'ref_aligned': ref_seq,
            'alt_aligned': alt_seq,
            'insertions': [],
            'deletions': [],
            'substitutions': [],
            'inversions': [],
            'length_change': len(alt_seq) - len(ref_seq)
        }
    
    best_alignment = alignments[0]
    score = best_alignment.score
    
    # Convert alignment to string representation
    alignment_strings = str(best_alignment).split('\n')
    
    # The alignment output format is:
    # Line 0: reference sequence with gaps
    # Line 1: alignment characters (| for match, space for mismatch)
    # Line 2: query sequence with gaps
    if len(alignment_strings) >= 3:
        ref_aligned = alignment_strings[0]
        alt_aligned = alignment_strings[2]
    else:
        # Fallback if alignment string format is different
        ref_aligned = ref_seq
        alt_aligned = alt_seq
    
    # Identify changes from the alignment
    changes = analyze_alignment_changes(ref_aligned, alt_aligned)
    
    return {
        'alignment_score': score,
        'ref_aligned': ref_aligned,
        'alt_aligned': alt_aligned,
        'insertions': changes['insertions'],
        'deletions': changes['deletions'],
        'substitutions': changes['substitutions'],
        'inversions': changes['inversions'],
        'complex_regions': changes['complex_regions'],
        'length_change': len(alt_seq) - len(ref_seq)
    }

def analyze_alignment_changes(ref_aligned, alt_aligned):
    """
    Analyze an alignment to identify mutations, insertions, deletions, and potential inversions.
    
    Args:
        ref_aligned: Reference sequence with alignment gaps
        alt_aligned: Alternate sequence with alignment gaps
        
    Returns:
        Dictionary with lists of detected changes
    """
    insertions = []
    deletions = []
    substitutions = []
    complex_regions = []
    
    # Find simple insertions, deletions, and substitutions
    current_ins = None
    current_del = None
    
    for i in range(len(ref_aligned)):
        ref_base = ref_aligned[i]
        alt_base = alt_aligned[i]
        
        if ref_base == '-' and alt_base != '-':
            # Insertion in alternate sequence
            if current_ins is None:
                current_ins = {'start': i, 'sequence': alt_base}
            else:
                current_ins['sequence'] += alt_base
        elif ref_base != '-' and alt_base == '-':
            # Deletion in alternate sequence
            if current_del is None:
                current_del = {'start': i, 'sequence': ref_base}
            else:
                current_del['sequence'] += ref_base
        elif ref_base != alt_base:
            # Substitution
            substitutions.append({'position': i, 'ref': ref_base, 'alt': alt_base})
            
            # Close any open indels
            if current_ins:
                insertions.append(current_ins)
                current_ins = None
            if current_del:
                deletions.append(current_del)
                current_del = None
        else:
            # Matching position
            # Close any open indels
            if current_ins:
                insertions.append(current_ins)
                current_ins = None
            if current_del:
                deletions.append(current_del)
                current_del = None
    
    # Close any open indels at the end
    if current_ins:
        insertions.append(current_ins)
    if current_del:
        deletions.append(current_del)
    
    # Look for potential inversions
    inversions = detect_inversions(ref_aligned, alt_aligned)
    
    # Identify complex regions (potential rearrangements)
    complex_regions = identify_complex_regions(ref_aligned, alt_aligned)
    
    return {
        'insertions': insertions,
        'deletions': deletions,
        'substitutions': substitutions,
        'inversions': inversions,
        'complex_regions': complex_regions
    }

def detect_inversions(ref_aligned, alt_aligned):
    """
    Detect potential inversions by looking for regions where the alternate sequence
    matches the reverse complement of the reference.
    
    Args:
        ref_aligned: Reference sequence with alignment gaps
        alt_aligned: Alternate sequence with alignment gaps
        
    Returns:
        List of detected inversions
    """
    inversions = []
    
    # Strip gaps for sequence comparison
    ref_seq = ref_aligned.replace('-', '')
    alt_seq = alt_aligned.replace('-', '')
    
    # Minimum inversion size to consider (to avoid random matches)
    min_inversion_size = 10
    
    # Look for potential inversions of various sizes
    for i in range(len(ref_seq) - min_inversion_size):
        for j in range(i + min_inversion_size, min(i + 200, len(ref_seq) + 1)):
            ref_segment = ref_seq[i:j]
            ref_revcomp = str(Seq(ref_segment).reverse_complement())
            
            # Search for reverse complement in alternate sequence
            for match in re.finditer(re.escape(ref_revcomp), alt_seq):
                start_alt = match.start()
                end_alt = match.end()
                
                inversions.append({
                    'ref_start': i,
                    'ref_end': j,
                    'alt_start': start_alt,
                    'alt_end': end_alt,
                    'length': j - i,
                    'ref_segment': ref_segment,
                    'alt_segment': ref_revcomp
                })
    
    # Filter overlapping inversions to keep only the most significant ones
    if inversions:
        inversions.sort(key=lambda x: x['length'], reverse=True)
        filtered_inversions = [inversions[0]]
        
        for inv in inversions[1:]:
            # Check if this inversion overlaps with any already selected
            overlaps = False
            for selected in filtered_inversions:
                # Check ref overlap
                ref_overlap = not (inv['ref_end'] <= selected['ref_start'] or inv['ref_start'] >= selected['ref_end'])
                # Check alt overlap
                alt_overlap = not (inv['alt_end'] <= selected['alt_start'] or inv['alt_start'] >= selected['alt_end'])
                
                if ref_overlap or alt_overlap:
                    overlaps = True
                    break
            
            if not overlaps:
                filtered_inversions.append(inv)
        
        return filtered_inversions
    
    return inversions

def identify_complex_regions(ref_aligned, alt_aligned):
    """
    Identify complex regions that might represent rearrangements or complex changes.
    
    Args:
        ref_aligned: Reference sequence with alignment gaps
        alt_aligned: Alternate sequence with alignment gaps
        
    Returns:
        List of complex regions
    """
    complex_regions = []
    
    # Initialize variables to track complex regions
    in_complex_region = False
    current_region = None
    mismatch_count = 0
    window_size = 20
    complexity_threshold = 8  # Number of mismatches in window_size to consider complex
    
    for i in range(len(ref_aligned)):
        # Check for region with high density of mismatches
        start_window = max(0, i - window_size)
        mismatches_in_window = sum(1 for j in range(start_window, i + 1) 
                                   if j < len(ref_aligned) and ref_aligned[j] != alt_aligned[j])
        
        if mismatches_in_window >= complexity_threshold:
            if not in_complex_region:
                in_complex_region = True
                current_region = {
                    'start': max(0, i - window_size),
                    'ref_sequence': ref_aligned[max(0, i - window_size):i + 1],
                    'alt_sequence': alt_aligned[max(0, i - window_size):i + 1]
                }
            else:
                # Extend current complex region
                current_region['ref_sequence'] += ref_aligned[i:i+1]
                current_region['alt_sequence'] += alt_aligned[i:i+1]
        else:
            if in_complex_region:
                in_complex_region = False
                current_region['end'] = i
                # Remove gaps for a cleaner representation
                current_region['ref_sequence_nogap'] = current_region['ref_sequence'].replace('-', '')
                current_region['alt_sequence_nogap'] = current_region['alt_sequence'].replace('-', '')
                complex_regions.append(current_region)
                current_region = None
    
    # Close any open complex region at the end
    if in_complex_region:
        current_region['end'] = len(ref_aligned)
        current_region['ref_sequence_nogap'] = current_region['ref_sequence'].replace('-', '')
        current_region['alt_sequence_nogap'] = current_region['alt_sequence'].replace('-', '')
        complex_regions.append(current_region)
    
    return complex_regions

def analyze_feature_with_alignment(feature, ref_seq, alt_seq):
    """
    Analyze a feature using sequence alignment to determine the effects of variants.
    
    Args:
        feature: Feature dictionary
        ref_seq: Reference sequence
        alt_seq: Alternate sequence
        
    Returns:
        Dictionary with effect analysis results
    """
    # Extract feature boundaries
    start = feature['start'] - 1  # Convert to 0-based
    end = feature['end'] - 1      # Convert to 0-based
    feature_type = feature['type']
    strand = feature['strand']
    
    # Skip if feature is outside sequence bounds
    if start >= len(ref_seq) or end >= len(ref_seq):
        logging.warning(f"Feature {feature['attributes'].get('ID', 'unknown')} is outside sequence bounds, skipping")
        return {
            'feature': feature,
            'feature_type': feature_type,
            'effects': ['no_change'],
            'details': {}
        }
    
    # Extract feature sequence from reference
    ref_feature_seq = ref_seq[start:end+1]
    
    # Extract feature sequence from alternate
    alt_feature_seq = ""
    if end < len(alt_seq):
        alt_feature_seq = alt_seq[start:end+1]
    else:
        # Handle case where alternate sequence is shorter than reference
        if start < len(alt_seq):
            alt_feature_seq = alt_seq[start:]
        # Feature is beyond the end of the alternate sequence
        else:
            alt_feature_seq = ""
    
    # If the sequences are identical, no change
    if ref_feature_seq == alt_feature_seq:
        return {
            'feature': feature,
            'feature_type': feature_type,
            'ref_feature_seq': ref_feature_seq,
            'alt_feature_seq': alt_feature_seq,
            'effects': ['no_change'],
            'details': {}
        }
    
    # Perform sequence alignment for detailed analysis
    alignment_result = analyze_sequence_alignment(ref_feature_seq, alt_feature_seq)
    
    # Determine effects based on alignment
    effects = []
    effect_details = {}
    
    # Calculate length change
    length_change = len(alt_feature_seq) - len(ref_feature_seq)
    if length_change != 0:
        effects.append('length_change')
        effect_details['length_change'] = length_change
    
    # Check for inversions
    if alignment_result['inversions']:
        effects.append('inversion')
        effect_details['inversions'] = alignment_result['inversions']
    
    # Check for complex regions (potential rearrangements)
    if alignment_result['complex_regions']:
        effects.append('complex_rearrangement')
        effect_details['complex_regions'] = alignment_result['complex_regions']
    
    # Check for insertions
    if alignment_result['insertions']:
        effects.append('insertion')
        effect_details['insertions'] = alignment_result['insertions']
    
    # Check for deletions
    if alignment_result['deletions']:
        effects.append('deletion')
        effect_details['deletions'] = alignment_result['deletions']
    
    # Check for substitutions
    if alignment_result['substitutions']:
        effects.append('substitution')
        effect_details['substitutions'] = alignment_result['substitutions']
    
    # Add alignment information
    effect_details['alignment'] = {
        'score': alignment_result['alignment_score'],
        'ref_aligned': alignment_result['ref_aligned'],
        'alt_aligned': alignment_result['alt_aligned']
    }
    
    # Check for more specific effects based on feature type
    if feature_type == 'CDS':
        # Check if length change is a multiple of 3 (in-frame)
        if length_change % 3 == 0:
            effects.append('in_frame_change')
        else:
            effects.append('frame_shift')
            effect_details['frame_shift'] = length_change % 3
        
        # Translate sequences to check for amino acid changes
        ref_aa = translate_sequence(ref_feature_seq, strand)
        alt_aa = translate_sequence(alt_feature_seq, strand)
        
        # Check for premature stop codons
        if '*' in alt_aa and (not '*' in ref_aa or alt_aa.index('*') < ref_aa.index('*') if '*' in ref_aa else True):
            effects.append('premature_stop_codon')
            effect_details['premature_stop_position'] = alt_aa.index('*') * 3
        
        # Check for amino acid changes
        aa_changes = compare_amino_acid_sequences(ref_aa, alt_aa)
        if aa_changes['changes'] > 0:
            effects.append('amino_acid_change')
            effect_details['amino_acid_changes'] = aa_changes
        
        # Check start codon disruption
        if feature.get('is_first_cds', False) and alt_aa and (not alt_aa.startswith('M')):
            effects.append('start_codon_disruption')
    
    # Check for regulatory region effects
    if feature_type in ['promoter', 'terminator']:
        effects.append(f"{feature_type}_affected")
    
    # Check for splicing region effects
    if feature_type == 'exon':
        effects.append('splicing_affected')
    
    # Remove duplicates
    effects = list(set(effects))
    
    return {
        'feature': feature,
        'feature_type': feature_type,
        'ref_feature_seq': ref_feature_seq,
        'alt_feature_seq': alt_feature_seq,
        'effects': effects,
        'details': effect_details
    }

def analyze_haplotype_differences_with_alignment(ref_seq, alt_seq, features, segment_offsets=None, variant_segments=None, segments=None):
    """
    Analyze differences between reference and alternate haplotype sequences using sequence alignment.
    
    Args:
        ref_seq: Reference sequence string
        alt_seq: Alternate sequence string
        features: List of feature dictionaries
        segment_offsets: Optional segment offset information
        variant_segments: Optional variant segment information
        segments: Optional segment data
        
    Returns:
        List of feature effect dictionaries
    """
    feature_effects = []
    
    # Find variants from segment information if available
    variants = []
    if segment_offsets and variant_segments and segments:
        # Find segments in the alternate sequence that correspond to variants
        ref_segments = []
        for seg_id in sorted(segments.keys()):
            seg_data = segments[seg_id]
            if not seg_data.get('variant_id') and not '_VAR' in seg_id:
                ref_segments.append({
                    'seg_id': seg_id,
                    'orientation': '+',
                    'length': seg_data['length']
                })
        
        variants = identify_variants_by_position(segment_offsets, ref_segments, variant_segments, segments)
    
    # Analyze each feature using alignment
    for feature in features:
        effect = analyze_feature_with_alignment(feature, ref_seq, alt_seq)
        
        # Add variants if available
        if variants:
            # Find variants that overlap with this feature
            feature_start = feature['start'] - 1  # Convert to 0-based
            feature_end = feature['end'] - 1      # Convert to 0-based
            
            overlapping_variants = []
            for variant in variants:
                var_start = variant.get('pos', 0) - 1  # Convert to 0-based
                var_end = variant.get('end', var_start)
                
                # If positions are unknown (0), check if the variant's segments are in this path
                if var_start == -1 and var_end == -1:
                    # Include the variant if it appears in our segment offsets
                    var_segments = variant.get('segments', [])
                    if segment_offsets and any(seg['seg_id'] in var_segments for seg in segment_offsets):
                        overlapping_variants.append(variant)
                    continue
                
                # Check for overlap with normal positions
                if not (var_end < feature_start or var_start > feature_end):
                    overlapping_variants.append(variant)
            
            effect['variants'] = overlapping_variants
        
        feature_effects.append(effect)
    
    return feature_effects
        
def setup_logging(debug=False, log_file=None, verbose=False):
    """Configure logging based on debug flag and optional log file."""
    if debug:
        log_level = logging.DEBUG
        log_format = '%(asctime)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s'
    elif verbose:
        log_level = logging.INFO
        log_format = '%(asctime)s - %(levelname)s - %(message)s'
    else:
        log_level = logging.INFO
        log_format = '%(message)s'
    
    # Configure root logger
    logging.basicConfig(level=log_level, format=log_format)
    logger = logging.getLogger()
    
    # Clear any existing handlers
    logger.handlers = []
    
    # Add console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_formatter = logging.Formatter(log_format)
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)
    
    # Add file handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_level)
        file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s')
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    
    return logger

def parse_gfa(gfa_file):
    """Parse a GFA file and extract segments, paths, and variant information."""
    segments = {}
    links = []
    paths = {}
    variants = []
    samples = defaultdict(list)
    haplotypes = defaultdict(dict)
    
    start_time = time.time()
    logging.info(f"Parsing GFA file: {gfa_file}")
    
    # First pass: collect segment information
    with open(gfa_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            record_type = fields[0]
            
            if record_type == 'S':  # Segment
                if len(fields) < 3:
                    logging.warning(f"Line {line_num}: Invalid segment record, missing fields")
                    continue
                    
                seg_id = fields[1]
                sequence = fields[2]
                
                # Extract variant info from segment name
                var_id = None
                var_type = None
                length_change = 0
                
                # Extract optional tags
                tags = {}
                for field in fields[3:]:
                    if ':' in field:
                        tag_parts = field.split(':', 2)
                        if len(tag_parts) >= 3:
                            tag_name, tag_type, tag_value = tag_parts
                            tags[tag_name] = (tag_type, tag_value)
                
                # Check if segment name contains variant info (e.g., S25_VAR16)
                if '_VAR' in seg_id:
                    parts = seg_id.split('_VAR')
                    if len(parts) > 1:
                        # The original segment ID is the part before _VAR
                        orig_seg_id = parts[0]
                        # Extract variant ID from segment name
                        var_id = 'VAR' + parts[1].split('_')[0]  # Handle multiple VAR tags
                    
                    # Try to infer variant type from segment name or tags
                    if "SNP" in seg_id or "SNV" in seg_id:
                        var_type = "SNP"
                    elif "DEL" in seg_id:
                        var_type = "DEL"
                    elif "INS" in seg_id:
                        var_type = "INS"
                    elif "DUP" in seg_id:
                        var_type = "DUP"
                    elif "INV" in seg_id:
                        var_type = "INV"
                
                # Extract variant info from tags
                if 'VA' in tags and tags['VA'][0] == 'Z':
                    variant_ids = tags['VA'][1].split(';')
                    if var_id is None and variant_ids:
                        var_id = variant_ids[0]
                
                # Extract length difference if available
                if 'LD' in tags and tags['LD'][0] == 'i':
                    try:
                        length_change = int(tags['LD'][1])
                    except ValueError:
                        pass
                
                # Store segment information
                segments[seg_id] = {
                    'sequence': sequence,
                    'length': len(sequence),
                    'tags': tags,
                    'variant_id': var_id,
                    'variant_type': var_type,
                    'length_change': length_change,
                    'original_segment': orig_seg_id if '_VAR' in seg_id else None
                }
                
                # Add to variants list if this is a variant segment
                if var_id:
                    # Check if variant already exists
                    existing_var = next((v for v in variants if v['id'] == var_id), None)
                    if existing_var:
                        # Update existing variant entry
                        if seg_id not in existing_var['segments']:
                            existing_var['segments'].append(seg_id)
                    else:
                        # Create new variant entry
                        variants.append({
                            'id': var_id,
                            'type': var_type if var_type else 'UNKNOWN',
                            'segments': [seg_id],
                            'length_change': length_change
                        })
            
            elif record_type == 'L':  # Link
                if len(fields) < 6:
                    logging.warning(f"Line {line_num}: Invalid link record, missing fields")
                    continue
                    
                from_id = fields[1]
                from_dir = fields[2]
                to_id = fields[3]
                to_dir = fields[4]
                overlap = fields[5]
                
                links.append((from_id, from_dir, to_id, to_dir, overlap))
    
    # Second pass: collect path information
    with open(gfa_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            record_type = fields[0]
            
            if record_type == 'P':  # Path
                if len(fields) < 3:
                    logging.warning(f"Line {line_num}: Invalid path record, missing fields")
                    continue
                    
                path_name = fields[1]
                path_segments = []
                
                for seg in fields[2].split(','):
                    # Split into segment ID and orientation
                    if not seg or len(seg) < 2:
                        logging.warning(f"Line {line_num}: Invalid segment in path: {seg}")
                        continue
                        
                    seg_id = seg[:-1]
                    orientation = seg[-1]
                    
                    path_segments.append((seg_id, orientation))
                
                # Extract tags if present
                path_tags = {}
                for field in fields[3:]:
                    if ':' in field:
                        tag_parts = field.split(':', 2)
                        if len(tag_parts) >= 3:
                            tag_name, tag_type, tag_value = tag_parts
                            path_tags[tag_name] = (tag_type, tag_value)
                
                # Extract sample and haplotype information
                sample_name = None
                haplotype = None
                variant_ids = []
                
                # First try to get info from tags
                if 'SM' in path_tags and path_tags['SM'][0] == 'Z':
                    sample_name = path_tags['SM'][1]
                
                if 'PH' in path_tags and path_tags['PH'][0] == 'Z':
                    haplotype = path_tags['PH'][1]
                
                if 'VA' in path_tags and path_tags['VA'][0] == 'Z':
                    variant_ids = path_tags['VA'][1].split(';')
                
                # If no sample info from tags, try to extract from path name
                if not sample_name:
                    # Try format like SAMPLE_Sample1_HAP1_1
                    match = re.match(r'SAMPLE_(.+)_HAP(\d+)_\d+', path_name)
                    if match:
                        sample_name = match.group(1)
                        hap_num = match.group(2)
                        haplotype = f"haplotype_{hap_num}"
                
                # If still no sample info but path contains variant segments, 
                # create a generic sample name
                if not sample_name:
                    has_variant_segments = any(seg_id in segments and segments[seg_id].get('variant_id') 
                                               for seg_id, _ in path_segments)
                    if has_variant_segments and path_name != 'REF':
                        sample_name = f"Sample_{path_name}"
                        haplotype = "haplotype_1"
                
                # Register sample and haplotype
                if sample_name:
                    samples[sample_name].append(path_name)
                    if haplotype:
                        haplotypes[sample_name][haplotype] = path_name
                
                # Extract variant IDs from path segments if not already defined
                if not variant_ids:
                    for seg_id, _ in path_segments:
                        if seg_id in segments:
                            var_id = segments[seg_id].get('variant_id')
                            if var_id and var_id not in variant_ids:
                                variant_ids.append(var_id)
                
                paths[path_name] = {
                    'segments': path_segments,
                    'tags': path_tags,
                    'sample': sample_name,
                    'haplotype': haplotype,
                    'variant_ids': variant_ids
                }
                
                # Associate variants with paths and samples
                for var_id in variant_ids:
                    # Find the corresponding variant
                    variant = next((v for v in variants if v['id'] == var_id), None)
                    if variant:
                        # Add path information to variant
                        if 'paths' not in variant:
                            variant['paths'] = []
                        if path_name not in variant['paths']:
                            variant['paths'].append(path_name)
                        
                        # Add sample information to variant
                        if sample_name:
                            if 'samples' not in variant:
                                variant['samples'] = []
                            if sample_name not in variant['samples']:
                                variant['samples'].append(sample_name)
    
    # Find or create REF path
    if 'REF' not in paths:
        logging.info("REF path not found, looking for reference path")
        # Try to find a path without variants
        for name, path_data in paths.items():
            if not path_data['variant_ids'] and not path_data['sample']:
                paths['REF'] = path_data
                logging.info(f"Using {name} as reference path")
                break
    
    # If still no REF path, create one from non-variant segments
    if 'REF' not in paths:
        logging.info("Creating synthetic REF path from non-variant segments")
        ref_segments = []
        non_variant_segments = [seg_id for seg_id, data in segments.items() 
                               if not data.get('variant_id')]
        
        # Sort segments by any available numbering
        def extract_number(seg_id):
            match = re.search(r'S(\d+)', seg_id)
            return int(match.group(1)) if match else 0
        
        sorted_segments = sorted(non_variant_segments, key=extract_number)
        
        if sorted_segments:
            ref_segments = [(seg_id, '+') for seg_id in sorted_segments]
            paths['REF'] = {
                'segments': ref_segments,
                'tags': {},
                'sample': None,
                'haplotype': None,
                'variant_ids': []
            }
        else:
            logging.warning("Could not create REF path, no non-variant segments found")
    
    elapsed = time.time() - start_time
    logging.info(f"Finished parsing GFA in {elapsed:.2f}s: {len(segments)} segments, {len(paths)} paths, {len(variants)} variants")
    logging.info(f"Found {len(samples)} samples with paths")
    
    return segments, links, paths, variants, samples, haplotypes

def parse_gff3(gff_file):
    """Parse a GFF3 file into a list of features."""
    features = []
    
    start_time = time.time()
    logging.info(f"Parsing GFF3 file: {gff_file}")
    
    with open(gff_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            fields = line.split('\t')
            
            if len(fields) < 9:
                logging.warning(f"Line {line_num}: Invalid GFF3 record, missing fields")
                continue
                
            seqid = fields[0]
            source = fields[1]
            feature_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            score = fields[5]
            strand = fields[6]
            phase = fields[7]
            
            # Parse attributes
            attributes = {}
            for attr in fields[8].split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attributes[key] = value
            
            feature = {
                'seqid': seqid,
                'source': source,
                'type': feature_type,
                'start': start,
                'end': end,
                'score': score,
                'strand': strand,
                'phase': phase,
                'attributes': attributes,
                'line_num': line_num
            }
            
            features.append(feature)
    
    elapsed = time.time() - start_time
    logging.info(f"Finished parsing GFF3 in {elapsed:.2f}s: {len(features)} features")
    
    # Build feature hierarchy
    feature_by_id = {}
    children_by_parent = defaultdict(list)
    
    for feature in features:
        if 'ID' in feature['attributes']:
            feature_id = feature['attributes']['ID']
            feature_by_id[feature_id] = feature
            
        if 'Parent' in feature['attributes']:
            parent_id = feature['attributes']['Parent']
            children_by_parent[parent_id].append(feature)
    
    # Attach children to their parents
    for feature in features:
        if 'ID' in feature['attributes']:
            feature_id = feature['attributes']['ID']
            feature['children'] = children_by_parent.get(feature_id, [])
    
    return features, feature_by_id, children_by_parent

def build_path_sequence(segments, path_segments):
    """Build a sequence from path segments."""
    sequence = ""
    segment_offsets = []
    current_offset = 0
    
    for seg_id, orientation in path_segments:
        if seg_id not in segments:
            logging.warning(f"Missing segment: {seg_id}")
            continue
            
        segment_seq = segments[seg_id]['sequence']
        
        # Handle reverse complement if needed
        if orientation == '-':
            seq_obj = Seq(segment_seq)
            segment_seq = str(seq_obj.reverse_complement())
            
        # Store segment offset information
        segment_offsets.append({
            'seg_id': seg_id,
            'orientation': orientation,
            'start': current_offset,
            'end': current_offset + len(segment_seq) - 1,
            'length': len(segment_seq),
            'variant_id': segments[seg_id].get('variant_id')
        })
        
        sequence += segment_seq
        current_offset += len(segment_seq)
    
    return sequence, segment_offsets

def identify_first_cds_in_genes(features, feature_by_id, children_by_parent):
    """Identify the first CDS feature in each gene to check for start codon changes."""
    # For each mRNA, find all its CDS children and mark the first one
    for feature in features:
        if feature['type'] == 'mRNA':
            mrna_id = feature['attributes'].get('ID')
            if not mrna_id:
                continue
                
            # Get all CDS features for this mRNA
            cds_features = [f for f in children_by_parent.get(mrna_id, []) if f['type'] == 'CDS']
            
            # Sort by position, considering strand
            if feature['strand'] == '+':
                cds_features.sort(key=lambda x: x['start'])
            else:
                cds_features.sort(key=lambda x: x['start'], reverse=True)
            
            # Mark the first CDS
            if cds_features:
                cds_features[0]['is_first_cds'] = True
    
    return features

def translate_sequence(nucleotide_seq, strand='+'):
    """Translate a nucleotide sequence to amino acids, considering strand."""
    if not nucleotide_seq:
        return ""
    
    # Handle reverse strand
    if strand == '-':
        seq_obj = Seq(nucleotide_seq)
        nucleotide_seq = str(seq_obj.reverse_complement())
    
    # Ensure sequence length is a multiple of 3
    remainder = len(nucleotide_seq) % 3
    if remainder > 0:
        nucleotide_seq = nucleotide_seq[:-remainder]
    
    if not nucleotide_seq:
        return ""
    
    # Translate to amino acids
    seq_obj = Seq(nucleotide_seq)
    return str(seq_obj.translate())

def compare_amino_acid_sequences(ref_aa, alt_aa):
    """Compare two amino acid sequences and identify changes."""
    changes = []
    substitutions = []
    inserted = []
    deleted = []
    
    # Handle case where sequences are identical
    if ref_aa == alt_aa:
        return {
            'changes': 0,
            'substitutions': 0,
            'insertions': 0,
            'deletions': 0,
            'details': []
        }
    
    # Find the common prefix and suffix
    prefix_len = 0
    while prefix_len < min(len(ref_aa), len(alt_aa)) and ref_aa[prefix_len] == alt_aa[prefix_len]:
        prefix_len += 1
    
    suffix_len = 0
    while (suffix_len < min(len(ref_aa), len(alt_aa)) - prefix_len and 
           ref_aa[len(ref_aa) - 1 - suffix_len] == alt_aa[len(alt_aa) - 1 - suffix_len]):
        suffix_len += 1
    
    # Extract the changed regions
    ref_changed = ref_aa[prefix_len:len(ref_aa) - suffix_len]
    alt_changed = alt_aa[prefix_len:len(alt_aa) - suffix_len]
    
    # Determine the type of changes
    if len(ref_changed) == 0 and len(alt_changed) > 0:
        # Pure insertion
        changes.append(f"Insertion of {len(alt_changed)} AA at position {prefix_len+1}")
        inserted = list(alt_changed)
    elif len(alt_changed) == 0 and len(ref_changed) > 0:
        # Pure deletion
        changes.append(f"Deletion of {len(ref_changed)} AA at position {prefix_len+1}")
        deleted = list(ref_changed)
    elif len(ref_changed) == len(alt_changed):
        # Pure substitution(s)
        for i in range(len(ref_changed)):
            if ref_changed[i] != alt_changed[i]:
                pos = prefix_len + i + 1
                changes.append(f"{ref_changed[i]}{pos}{alt_changed[i]}")
                substitutions.append((pos, ref_changed[i], alt_changed[i]))
    else:
        # Complex change (indel + substitution)
        changes.append(f"Complex change: {ref_changed} -> {alt_changed}")
        
        if len(ref_changed) < len(alt_changed):
            # Net insertion
            inserted = list(alt_changed)
            for aa in ref_changed:
                if aa in inserted:
                    inserted.remove(aa)
        elif len(ref_changed) > len(alt_changed):
            # Net deletion
            deleted = list(ref_changed)
            for aa in alt_changed:
                if aa in deleted:
                    deleted.remove(aa)
    
    return {
        'changes': len(changes),
        'substitutions': len(substitutions),
        'insertions': len(inserted),
        'deletions': len(deleted),
        'details': changes
    }

def identify_variants_by_position(segment_offsets, ref_segments, variant_segments, segments=None):
    """Identify variants in an alternate sequence based on segment differences."""
    variants_by_position = []
    
    # If no reference segments provided, return empty list
    if not ref_segments:
        return variants_by_position
        
    # Create a map of reference positions
    ref_pos_map = {}
    current_pos = 0
    
    for seg_info in ref_segments:
        seg_id = seg_info.get('seg_id')
        length = seg_info.get('length', 0)
        ref_pos_map[seg_id] = (current_pos, current_pos + length - 1)
        current_pos += length
    
    # Find alternate segments that differ from reference
    for seg_info in segment_offsets:
        seg_id = seg_info.get('seg_id')
        start_pos = seg_info.get('start', 0)
        length = seg_info.get('length', 0)
        variant_id = seg_info.get('variant_id')
        
        # Skip segments without variant IDs
        if not variant_id:
            variant_id = segments.get(seg_id, {}).get('variant_id') if segments else None
            if not variant_id:
                continue
        
        # Get variant information
        var_data = variant_segments.get(seg_id, {})
        if not var_data and segments:
            var_data = segments.get(seg_id, {})
        
        # Try to find the original segment
        ref_id = var_data.get('original_segment')
        if not ref_id and '_VAR' in seg_id:
            # Attempt to extract original segment from ID
            parts = seg_id.split('_VAR')
            if len(parts) > 1:
                ref_id = parts[0]
        
        # Find reference position for this variant
        if ref_id and ref_id in ref_pos_map:
            ref_start, ref_end = ref_pos_map[ref_id]
            # Create a variant entry
            variant = {
                'id': variant_id,
                'type': var_data.get('variant_type', 'UNKNOWN'),
                'pos': ref_start + 1,  # Convert to 1-based
                'end': ref_end + 1,    # Convert to 1-based
                'length_change': var_data.get('length_change', 0),
                'segments': [seg_id],
                'original_segment': ref_id
            }
            variants_by_position.append(variant)
    
    return variants_by_position


def analyze_haplotype_differences(ref_seq, alt_seq, features, segment_offsets, variant_segments, segments):
    """Analyze differences between reference and alternate haplotype sequences."""
    feature_effects = []
    
    # Identify variants by comparing segment positions
    ref_segments = []
    for seg_info in segment_offsets:
        seg_id = seg_info['seg_id']
        if seg_id in segments and not segments[seg_id].get('variant_id'):
            ref_segments.append({
                'seg_id': seg_id,
                'orientation': seg_info['orientation'],
                'length': segments[seg_id]['length']
            })
    
    variants = identify_variants_by_position(segment_offsets, ref_segments, variant_segments)
    
    for feature in features:
        # Extract feature boundaries
        start = feature['start'] - 1  # Convert to 0-based
        end = feature['end'] - 1      # Convert to 0-based
        feature_type = feature['type']
        strand = feature['strand']
        
        # Skip if feature is outside sequence bounds
        if start >= len(ref_seq) or end >= len(ref_seq):
            logging.warning(f"Feature {feature['attributes'].get('ID', 'unknown')} is outside sequence bounds, skipping")
            continue
        
        # Extract feature sequence from reference
        ref_feature_seq = ref_seq[start:end+1]
        
        # Extract feature sequence from alternate
        alt_feature_seq = ""
        if end < len(alt_seq):
            alt_feature_seq = alt_seq[start:end+1]
        else:
            # Handle case where alternate sequence is shorter than reference
            if start < len(alt_seq):
                alt_feature_seq = alt_seq[start:]
            # Feature is beyond the end of the alternate sequence
            else:
                alt_feature_seq = ""
        
        # Find overlapping variants
        overlapping_variants = []
        for variant in variants:
            var_start = variant.get('pos', 0) - 1  # Convert to 0-based
            var_end = variant.get('end', var_start)
            
            # Check for overlap
            if not (var_end < start or var_start > end):
                overlapping_variants.append(variant)
        
        if ref_feature_seq == alt_feature_seq and not overlapping_variants:
            # No changes to this feature
            effect = {
                'feature': feature,
                'feature_type': feature_type,
                'ref_feature_seq': ref_feature_seq,
                'alt_feature_seq': alt_feature_seq,
                'variants': [],
                'effects': ['no_change'],
                'details': {}
            }
        else:
            # Analyze effects of changes
            effects = []
            effect_details = {}
            
            # Calculate basic effects
            length_change = len(alt_feature_seq) - len(ref_feature_seq)
            
            if length_change != 0:
                effects.append('length_change')
                effect_details['length_change'] = length_change
            
            # Check for more specific effects
            if feature_type == 'CDS':
                # Check if length change is a multiple of 3 (in-frame)
                if length_change % 3 == 0:
                    effects.append('in_frame_change')
                else:
                    effects.append('frame_shift')
                    effect_details['frame_shift'] = length_change % 3
                
                # Translate sequences to check for amino acid changes
                ref_aa = translate_sequence(ref_feature_seq, strand)
                alt_aa = translate_sequence(alt_feature_seq, strand)
                
                # Check for premature stop codons
                if '*' in alt_aa and (not '*' in ref_aa or alt_aa.index('*') < ref_aa.index('*') if '*' in ref_aa else True):
                    effects.append('premature_stop_codon')
                    effect_details['premature_stop_position'] = alt_aa.index('*') * 3
                
                # Check for amino acid changes
                aa_changes = compare_amino_acid_sequences(ref_aa, alt_aa)
                if aa_changes['changes'] > 0:
                    effects.append('amino_acid_change')
                    effect_details['amino_acid_changes'] = aa_changes
                
                # Check start codon disruption
                if feature.get('is_first_cds', False) and alt_aa and (not alt_aa.startswith('M')):
                    effects.append('start_codon_disruption')
            
            # Check for regulatory region effects
            if feature_type in ['promoter', 'terminator']:
                effects.append(f"{feature_type}_affected")
            
            # Check for splicing region effects
            if feature_type == 'exon':
                effects.append('splicing_affected')
            
            # Remove duplicates
            effects = list(set(effects))
            
            effect = {
                'feature': feature,
                'feature_type': feature_type,
                'ref_feature_seq': ref_feature_seq,
                'alt_feature_seq': alt_feature_seq,
                'variants': overlapping_variants,
                'effects': effects,
                'details': effect_details
            }
        
        feature_effects.append(effect)
    
    return feature_effects

def analyze_sample_haplotypes(feature_effects_by_haplotype, features, feature_by_id=None, children_by_parent=None):
    """
    Analyze differences between sample haplotypes to identify homo/heterozygous effects
    and compound heterozygous variants.
    """
    if len(feature_effects_by_haplotype) < 2:
        # Cannot determine zygosity with less than 2 haplotypes
        return {'homozygous': [], 'heterozygous': [], 'incomplete': True}
    
    # Get haplotype names
    haplotype_names = list(feature_effects_by_haplotype.keys())
    
    # Map feature IDs to effects for each haplotype
    feature_map = {}
    for hap_name, effects in feature_effects_by_haplotype.items():
        for effect in effects:
            feature_id = effect['feature']['attributes'].get('ID', 'unknown')
            if feature_id not in feature_map:
                feature_map[feature_id] = {}
            
            feature_map[feature_id][hap_name] = effect
    
    # Identify homozygous vs heterozygous effects
    homozygous = []
    heterozygous = []
    
    for feature_id, hap_effects in feature_map.items():
        # Check if effects are identical across haplotypes
        effect_signatures = {}
        for hap_name, effect in hap_effects.items():
            # Create a signature of the effect for comparison
            # Include sequence and effect types to determine if truly identical
            sig = (
                tuple(sorted(effect['effects'])),
                effect['alt_feature_seq'],
                tuple(v.get('id', '') for v in effect.get('variants', []))
            )
            effect_signatures[hap_name] = sig
        
        # Skip features without effects on all haplotypes
        if len(hap_effects) < len(haplotype_names):
            continue
        
        # Check if all signatures are the same
        is_homozygous = len(set(effect_signatures.values())) == 1
        
        if is_homozygous:
            # Use the effect from the first haplotype
            first_hap = haplotype_names[0]
            effect = hap_effects[first_hap].copy()
            effect['zygosity'] = 'homozygous'
            homozygous.append(effect)
        else:
            # Create a heterozygous effect entry
            # We'll use the effect from the first haplotype but note the differences
            first_hap = haplotype_names[0]
            effect = hap_effects[first_hap].copy()
            effect['zygosity'] = 'heterozygous'
            effect['haplotype_effects'] = hap_effects
            
            # Analyze differences
            diff_effects = set()
            for hap_name, hap_effect in hap_effects.items():
                for e in hap_effect['effects']:
                    diff_effects.add(e)
            
            effect['combined_effects'] = list(diff_effects)
            heterozygous.append(effect)
    
    # Identify compound heterozygous effects if we have feature hierarchy info
    compound_heterozygous = []
    if feature_by_id and children_by_parent:
        compound_heterozygous = identify_compound_heterozygous_effects(
            feature_effects_by_haplotype, 
            features, 
            feature_by_id, 
            children_by_parent
        )
    
    return {
        'homozygous': homozygous,
        'heterozygous': heterozygous,
        'compound_heterozygous': compound_heterozygous,
        'incomplete': False
    }

def identify_compound_heterozygous_effects(feature_effects_by_haplotype, features, feature_by_id, children_by_parent):
    """
    Identify compound heterozygous effects where different variants affect the same gene in different haplotypes.
    
    Compound heterozygosity occurs when a gene has different variants on each chromosome copy.
    This is important for recessive conditions where two different variants can result in disease.
    """
    compound_het_effects = []
    
    # Need at least 2 haplotypes
    if len(feature_effects_by_haplotype) < 2:
        return compound_het_effects
    
    # Get list of haplotypes
    haplotype_names = list(feature_effects_by_haplotype.keys())
    
    # Map genes to their CDS features
    gene_to_cds = {}
    for feature in features:
        if feature['type'] == 'gene':
            gene_id = feature['attributes'].get('ID')
            if not gene_id:
                continue
                
            # Find all CDS features associated with this gene
            cds_features = []
            # First find mRNAs for this gene
            mrnas = children_by_parent.get(gene_id, [])
            for mrna in mrnas:
                mrna_id = mrna['attributes'].get('ID')
                if mrna_id:
                    # Find CDS features for this mRNA
                    cds_list = [f for f in children_by_parent.get(mrna_id, []) if f['type'] == 'CDS']
                    cds_features.extend(cds_list)
            
            if cds_features:
                gene_to_cds[gene_id] = cds_features
    
    # Create a map of gene_id -> haplotype -> variants
    gene_variants_by_haplotype = {}
    
    # For each haplotype's effects
    for hap_name, effects in feature_effects_by_haplotype.items():
        # For each effect
        for effect in effects:
            feature = effect['feature']
            feature_id = feature['attributes'].get('ID', 'unknown')
            feature_type = feature['type']
            
            # Skip if no variants or feature is not affected
            if not effect.get('variants') or effect.get('effects') == ['no_change']:
                continue
                
            # Only concerned with CDS features
            if feature_type != 'CDS':
                continue
                
            # Find the gene this CDS belongs to
            gene_id = None
            for gene, cds_list in gene_to_cds.items():
                if any(cds['attributes'].get('ID') == feature_id for cds in cds_list):
                    gene_id = gene
                    break
            
            if not gene_id:
                continue
                
            # Add variants to the map
            if gene_id not in gene_variants_by_haplotype:
                gene_variants_by_haplotype[gene_id] = {}
            
            if hap_name not in gene_variants_by_haplotype[gene_id]:
                gene_variants_by_haplotype[gene_id][hap_name] = []
            
            # Add only variants we haven't seen for this haplotype and gene
            for variant in effect.get('variants', []):
                # Skip if we've already added this variant
                if any(v.get('id') == variant.get('id') for v in gene_variants_by_haplotype[gene_id][hap_name]):
                    continue
                gene_variants_by_haplotype[gene_id][hap_name].append(variant)
    
    # Identify compound heterozygous effects
    for gene_id, haplotype_variants in gene_variants_by_haplotype.items():
        # Need variants in at least 2 haplotypes
        if len(haplotype_variants) < 2:
            continue
            
        # Check that all haplotypes have at least one variant
        if any(len(variants) == 0 for variants in haplotype_variants.values()):
            continue
        
        # Get variant IDs for each haplotype
        variant_ids_by_haplotype = {}
        for hap, variants in haplotype_variants.items():
            variant_ids_by_haplotype[hap] = [v.get('id') for v in variants if v.get('id')]
        
        # Check if different variants are affecting different haplotypes
        is_compound_het = False
        
        # Check if the variants are different between haplotypes
        variant_sets = [set(ids) for ids in variant_ids_by_haplotype.values()]
        
        # If any haplotype has a unique variant not shared with other haplotypes
        # This is a potential compound heterozygous situation
        if len(variant_sets) >= 2:
            for i, hap_variants in enumerate(variant_sets):
                other_variants = set().union(*[s for j, s in enumerate(variant_sets) if j != i])
                unique_variants = hap_variants - other_variants
                if unique_variants:
                    is_compound_het = True
                    break
        
        if is_compound_het:
            # Create compound heterozygous effect entry
            gene_feature = next((f for f in features if f['type'] == 'gene' and f['attributes'].get('ID') == gene_id), None)
            if not gene_feature:
                continue
                
            # Get gene name
            gene_name = gene_feature['attributes'].get('Name', gene_id)
            
            # Get all variants affecting this gene across haplotypes
            all_variants = []
            for hap, variants in haplotype_variants.items():
                for var in variants:
                    if not any(v.get('id') == var.get('id') for v in all_variants):
                        all_variants.append(var)
            
            # Create effect entry with detailed information
            compound_effect = {
                'gene_id': gene_id,
                'gene_name': gene_name,
                'haplotype_variants': haplotype_variants,
                'variants': all_variants,
                'compound_heterozygous': True,
                'feature': gene_feature,
                # Add more detailed analysis of the combined effect
                'details': {
                    'variants_by_haplotype': {hap: [v.get('id') for v in vars] 
                                             for hap, vars in haplotype_variants.items()},
                }
            }
            
            compound_het_effects.append(compound_effect)
    
    return compound_het_effects




def identify_repeat_sequences(variant, segments):
    """Identify repeat patterns in insertion or duplication variants."""
    if variant['type'] not in ['DUP', 'INS']:
        return None
    
    # Try to find the sequence for this variant
    variant_sequence = ""
    for seg_id in variant.get('segments', []):
        if seg_id in segments:
            variant_sequence = segments[seg_id]['sequence']
            break
    
    if not variant_sequence:
        return None
    
    # Look for repeating patterns
    repeat_info = {}
    
    # Try different pattern lengths, starting from shorter patterns
    for pattern_len in range(1, min(len(variant_sequence) // 2 + 1, 20)):
        pattern = variant_sequence[:pattern_len]
        
        # Count how many times the pattern repeats
        count = 0
        pos = 0
        while pos + pattern_len <= len(variant_sequence) and variant_sequence[pos:pos+pattern_len] == pattern:
            count += 1
            pos += pattern_len
        
        # If we found a significant repeat (at least 3 times)
        if count >= 3 and count * pattern_len >= len(variant_sequence) * 0.6:
            repeat_info['sequence'] = pattern
            repeat_info['count'] = count
            break
    
    return repeat_info or None

def generate_variant_effect_report(feature_effects, variants, outfile=None, sample_name=None, sample_effects=None):
    """Generate a comprehensive report of variant effects on features."""
    # Open output file if specified
    out = open(outfile, 'w') if outfile else sys.stdout
    
    # Write report header
    if sample_name:
        out.write(f"# Variant Effect Report for Sample: {sample_name}\n")
    else:
        out.write("# Variant Effect Report\n")
    out.write("#" + "=" * 79 + "\n\n")
    
    # Write variant summary
    out.write("## Variant Summary\n")
    out.write("-" * 80 + "\n")
    
    # Get all variants affecting this sample
    sample_variants = []
    if sample_name:
        # First approach: use direct sample association from GFA parsing
        sample_variants = [v for v in variants if 'samples' in v and sample_name in v.get('samples', [])]
        
        # Second approach: try to extract from sample_effects if available
        if not sample_variants and sample_effects:
            variant_ids = set()
            for effect_type in ['homozygous', 'heterozygous', 'gene_compound_heterozygous', 'feature_compound_heterozygous']:
                for effect in sample_effects.get(effect_type, []):
                    for var in effect.get('variants', []):
                        if 'id' in var:
                            variant_ids.add(var['id'])
            
            # Filter variants by these IDs
            if variant_ids:
                sample_variants = [v for v in variants if v.get('id') in variant_ids]
        
        # Third approach: try to extract from feature effects
        if not sample_variants and feature_effects:
            variant_ids = set()
            for effect in feature_effects:
                for var in effect.get('variants', []):
                    if 'id' in var:
                        variant_ids.add(var['id'])
            
            # Filter variants by these IDs
            if variant_ids:
                sample_variants = [v for v in variants if v.get('id') in variant_ids]
    else:
        sample_variants = variants
    
    # Output variant summary
    if sample_variants:
        for variant in sample_variants:
            zygosity = ""
            # Determine zygosity if available
            if sample_effects:
                for effect_type in ['homozygous', 'heterozygous']:
                    for effect in sample_effects.get(effect_type, []):
                        if any(v.get('id') == variant.get('id') for v in effect.get('variants', [])):
                            zygosity = effect_type.upper()
                            break
                    if zygosity:
                        break
            
            # Display variant basic info
            var_type = variant.get('type', 'UNKNOWN')
            var_id = variant.get('id', 'unknown')
            out.write(f"Variant: {var_id} ({var_type}){' - ' + zygosity if zygosity else ''}\n")
            
            # Display position if available
            if 'pos' in variant and variant.get('pos', 0) > 0:
                out.write(f"Position: {variant.get('pos', 'unknown')}-{variant.get('end', variant.get('pos', 'unknown'))}\n")
            
            # Display length change with appropriate units and sign
            length_change = variant.get('length_change', 0)
            if length_change != 0:
                sign = "+" if length_change > 0 else ""
                out.write(f"Length Change: {sign}{length_change} bp\n")
            else:
                # For SNPs and other zero-length variants
                if var_type == "SNP" or var_type == "SNV":
                    out.write("Length Change: 0 bp (substitution)\n")
                else:
                    out.write("Length Change: 0 bp\n")
            
            # Add segment information
            out.write(f"Affected Segments: {', '.join(variant.get('segments', []))}\n")
            
            # Add repeat information if available
            if variant.get('type') in ['DUP', 'INS'] and variant.get('repeat_info'):
                repeat_info = variant['repeat_info']
                out.write(f"Repeat Units: {repeat_info.get('count', 'unknown')}\n")
                out.write(f"Repeat Sequence: {repeat_info.get('sequence', 'unknown')}\n")
                
            out.write("\n")
    else:
        out.write("No variants identified for this sample.\n\n")
    
    # Group features by type
    features_by_type = defaultdict(list)
    
    # Use sample_effects if provided, otherwise use feature_effects
    if sample_name and sample_effects and not sample_effects.get('incomplete'):
        effects_to_process = []
        effects_to_process.extend(sample_effects.get('homozygous', []))
        effects_to_process.extend(sample_effects.get('heterozygous', []))
    else:
        effects_to_process = feature_effects
    
    for effect in effects_to_process:
        features_by_type[effect['feature_type']].append(effect)
    
    # Write feature effect summary
    if sample_name:
        out.write(f"\n## Feature Effect Summary for Sample: {sample_name}\n")
    else:
        out.write("\n## Feature Effect Summary\n")
    out.write("-" * 80 + "\n")
    
    # Calculate statistics
    effect_counts = Counter()
    affected_feature_counts = Counter()
    zygosity_counts = Counter()
    
    for effect in effects_to_process:
        if effect['effects'] and effect['effects'] != ['no_change']:  # If feature is affected
            affected_feature_counts[effect['feature_type']] += 1
            for e in effect['effects']:
                effect_counts[e] += 1
            
            # Count zygosity if available
            if 'zygosity' in effect:
                zygosity_counts[effect['zygosity']] += 1
    
    # Write statistics
    out.write(f"Total Features Analyzed: {len(effects_to_process)}\n")
    out.write(f"Features Affected by Variants: {sum(affected_feature_counts.values())}\n\n")
    
    # Add zygosity statistics for sample-specific reports
    if sample_name and zygosity_counts:
        out.write("Zygosity Summary:\n")
        out.write(f"  Homozygous Effects: {zygosity_counts.get('homozygous', 0)}\n")
        out.write(f"  Heterozygous Effects: {zygosity_counts.get('heterozygous', 0)}\n\n")
    
    out.write("Affected Features by Type:\n")
    for feature_type, count in sorted(affected_feature_counts.items()):
        total = len(features_by_type[feature_type])
        percentage = (count / total * 100) if total > 0 else 0
        out.write(f"  {feature_type}: {count}/{total} ({percentage:.1f}%)\n")
    
    out.write("\nVariant Effects by Type:\n")
    for effect_type, count in sorted(effect_counts.items(), key=lambda x: x[1], reverse=True):
        out.write(f"  {effect_type}: {count}\n")
    
    # Write detailed feature effects
    if sample_name:
        out.write(f"\n\n## Detailed Feature Effects for Sample: {sample_name}\n")
    else:
        out.write("\n\n## Detailed Feature Effects\n")
    
    # Process each feature type
    for feature_type, effects in sorted(features_by_type.items()):
        affected_effects = [e for e in effects if e['effects'] and e['effects'] != ['no_change']]
        
        if not affected_effects:
            continue  # Skip feature types not affected by variants
        
        out.write(f"\n### {feature_type.upper()} Features\n")
        out.write("-" * 80 + "\n")
        
        for effect in affected_effects:
            feature = effect['feature']
            feature_id = feature['attributes'].get('ID', 'unknown')
            feature_name = feature['attributes'].get('Name', feature_id)
            
            out.write(f"\nFeature: {feature_name} ({feature_id})\n")
            out.write(f"Location: {feature['start']}-{feature['end']} ({feature['strand']})\n")
            
            # Add zygosity information if available
            if 'zygosity' in effect:
                out.write(f"Zygosity: {effect['zygosity'].upper()}\n")
            
            # List variants affecting this feature
            if effect.get('variants'):
                out.write("Affected by variants:\n")
                for variant in effect['variants']:
                    out.write(f"  - {variant.get('id', 'unknown')} ({variant.get('type', 'UNKNOWN')})\n")
            
            # List effects
            out.write("Effects:\n")
            
            # Use combined_effects for heterozygous variants if available
            if 'zygosity' in effect and effect['zygosity'] == 'heterozygous' and 'combined_effects' in effect:
                effect_list = effect['combined_effects']
            else:
                effect_list = effect['effects']
                
            for effect_type in sorted(effect_list):
                if effect_type == 'no_change':
                    out.write("  - No effect on this feature\n")
                    continue
                
                # Format effect with details if available
                detail_str = ""
                if effect_type in effect['details']:
                    detail = effect['details'][effect_type]
                    if isinstance(detail, dict):
                        detail_str = ": " + ", ".join(f"{k}={v}" for k, v in detail.items())
                    elif isinstance(detail, list):
                        detail_str = ": " + ", ".join(str(d) for d in detail)
                    else:
                        detail_str = f": {detail}"
                
                out.write(f"  - {effect_type}{detail_str}\n")
            
            # For heterozygous effects, show haplotype-specific differences
            if 'zygosity' in effect and effect['zygosity'] == 'heterozygous' and 'haplotype_effects' in effect:
                out.write("\nHaplotype-specific Effects:\n")
                for hap_name, hap_effect in effect['haplotype_effects'].items():
                    out.write(f"  {hap_name}:\n")
                    for e in sorted(hap_effect['effects']):
                        if e == 'no_change':
                            continue
                        out.write(f"    - {e}\n")
            
            # Add sequence details for CDS features
            if feature_type == 'CDS' and 'amino_acid_change' in effect['effects']:
                aa_changes = effect['details'].get('amino_acid_changes', {})
                if aa_changes:
                    out.write("\nAmino Acid Changes:\n")
                    for detail in aa_changes.get('details', []):
                        out.write(f"  - {detail}\n")
            
            # Write sequence changes
            if len(effect['ref_feature_seq']) <= 50 and len(effect['alt_feature_seq']) <= 50:
                out.write("\nSequence Changes:\n")
                out.write(f"  Reference: {effect['ref_feature_seq']}\n")
                out.write(f"  Alternate: {effect['alt_feature_seq']}\n")
            else:
                out.write("\nSequence Changes (truncated):\n")
                out.write(f"  Reference: {effect['ref_feature_seq'][:25]}...{effect['ref_feature_seq'][-25:]}\n")
                out.write(f"  Alternate: {effect['alt_feature_seq'][:25]}...{effect['alt_feature_seq'][-25:]}\n")
            
            # Add frame information for CDS features
            if feature_type == 'CDS':
                if 'frame_shift' in effect['effects']:
                    frame_shift = effect['details'].get('frame_shift', 'unknown')
                    out.write(f"\nFrame Shift: {frame_shift} bases\n")
                elif 'in_frame_change' in effect['effects']:
                    out.write("\nIn-frame change (multiple of 3 bases)\n")
                
                if 'premature_stop_codon' in effect['effects']:
                    stop_pos = effect['details'].get('premature_stop_position', 'unknown')
                    out.write(f"Premature stop codon introduced at position: {stop_pos}\n")
    
    # Write gene-level compound heterozygous effects if available
    if sample_name and sample_effects and 'gene_compound_heterozygous' in sample_effects:
        compound_het_effects = sample_effects['gene_compound_heterozygous']
        if compound_het_effects:
            out.write("\n\n## Gene-Level Compound Heterozygous Effects\n")
            out.write("-" * 80 + "\n")
            
            for effect in compound_het_effects:
                gene_id = effect.get('gene_id', 'unknown')
                gene_name = effect.get('gene_name', gene_id)
                
                out.write(f"\nGene: {gene_name} ({gene_id})\n")
                out.write("Compound heterozygous effect: Different variants affecting the same gene on different haplotypes\n")
                
                # List variants by haplotype
                out.write("\nVariants by Haplotype:\n")
                for hap_name, hap_variants in effect.get('haplotype_variants', {}).items():
                    out.write(f"  {hap_name}:\n")
                    for variant in hap_variants:
                        out.write(f"    - {variant.get('id', 'unknown')} ({variant.get('type', 'UNKNOWN')})\n")
                
                out.write("\n")
    
    # Write feature-level compound heterozygous effects if available
    if sample_name and sample_effects and 'feature_compound_heterozygous' in sample_effects:
        feature_compound_het_effects = sample_effects['feature_compound_heterozygous']
        if feature_compound_het_effects:
            out.write("\n\n## Feature-Level Compound Heterozygous Effects\n")
            out.write("-" * 80 + "\n")
            out.write("This section lists features affected by different variants on different haplotypes.\n\n")
            
            # Group by feature type
            effects_by_type = defaultdict(list)
            for effect in feature_compound_het_effects:
                effects_by_type[effect.get('feature_type', 'unknown')].append(effect)
            
            # Process each feature type
            for feature_type, effects in sorted(effects_by_type.items()):
                out.write(f"\n### {feature_type.upper()} Features\n")
                
                for effect in effects:
                    feature_id = effect.get('feature_id', 'unknown')
                    feature_name = effect.get('feature_name', feature_id)
                    
                    out.write(f"\nFeature: {feature_name} ({feature_id})\n")
                    
                    # Get feature info
                    feature = effect.get('feature', {})
                    if feature:
                        out.write(f"Location: {feature.get('start', 'unknown')}-{feature.get('end', 'unknown')} ({feature.get('strand', 'unknown')})\n")
                    
                    out.write("Compound heterozygous effect: Different variants affecting this feature on different haplotypes\n")
                    
                    # List variants by haplotype
                    out.write("\nVariants by Haplotype:\n")
                    for hap_name, hap_variants in effect.get('haplotype_variants', {}).items():
                        out.write(f"  {hap_name}:\n")
                        for variant in hap_variants:
                            out.write(f"    - {variant.get('id', 'unknown')} ({variant.get('type', 'UNKNOWN')})\n")
                    
                    out.write("\n")
    
    # Close output file if opened
    if outfile:
        out.close()
        logging.info(f"Report written to {outfile}")


def calculate_variant_length_changes(variants, segments):
    """
    Calculate length changes for each variant by comparing variant segments to their original segments.
    
    Args:
        variants: List of variant dictionaries
        segments: Dictionary of segment data
        
    Returns:
        Updated list of variants with accurate length change values
    """
    for variant in variants:
        # Skip if length change is already set and non-zero
        if variant.get('length_change', 0) != 0:
            continue
            
        var_segments = variant.get('segments', [])
        if not var_segments:
            continue
            
        total_length_change = 0
        
        for seg_id in var_segments:
            if seg_id not in segments:
                continue
                
            # Get variant segment length
            var_seg_length = segments[seg_id].get('length', 0)
            
            # Try to get original segment ID
            orig_seg_id = segments[seg_id].get('original_segment')
            
            # If not directly specified, try to extract from segment name
            if not orig_seg_id and '_VAR' in seg_id:
                parts = seg_id.split('_VAR')
                if len(parts) > 1:
                    orig_seg_id = parts[0]
            
            # Calculate length change if we found the original segment
            if orig_seg_id and orig_seg_id in segments:
                orig_seg_length = segments[orig_seg_id].get('length', 0)
                segment_length_change = var_seg_length - orig_seg_length
                total_length_change += segment_length_change
            # Get length change from the segment's tag if available
            elif 'LD' in segments[seg_id].get('tags', {}):
                tag_data = segments[seg_id]['tags']['LD']
                if tag_data[0] == 'i':  # Integer type
                    try:
                        segment_length_change = int(tag_data[1])
                        total_length_change += segment_length_change
                    except ValueError:
                        pass
        
        # Update the variant with the calculated length change
        variant['length_change'] = total_length_change
    
    return variants
        
def generate_sample_reports(segments, paths, variants, features, feature_by_id, children_by_parent, outfile_prefix=None):
    """Generate sample-specific variant effect reports including compound heterozygous effects."""
    # Identify samples with paths
    samples = defaultdict(list)
    haplotypes = defaultdict(dict)
    
    for path_name, path_data in paths.items():
        sample_name = path_data.get('sample')
        if sample_name:
            samples[sample_name].append(path_name)
            
            # Track haplotype paths
            haplotype = path_data.get('haplotype')
            if haplotype and haplotype.startswith('haplotype_'):
                haplotypes[sample_name][haplotype] = path_name
    
    if not samples:
        logging.warning("No sample-specific paths found in GFA")
        return {}
    
    logging.info(f"Found {len(samples)} samples with paths")
    sample_reports = {}
    
    # Get reference path
    ref_path_name = 'REF'
    if ref_path_name not in paths:
        logging.error("REF path not found in GFA")
        return {}
    
    ref_path_segments = paths[ref_path_name]['segments']
    ref_seq, _ = build_path_sequence(segments, ref_path_segments)
    
    # Collect variant segments for lookup
    variant_segments = {}
    for seg_id, seg_data in segments.items():
        if seg_data.get('variant_id'):
            variant_segments[seg_id] = seg_data
    
    # Process each sample
    for sample_name, path_list in samples.items():
        logging.info(f"Generating report for sample: {sample_name}")
        
        # Determine if we have phased haplotypes for this sample
        sample_haplotypes = haplotypes.get(sample_name, {})
        has_phased_haplotypes = len(sample_haplotypes) >= 2
        
        if has_phased_haplotypes:
            logging.info(f"Found {len(sample_haplotypes)} phased haplotypes for sample {sample_name}")
            
            # Process each haplotype
            feature_effects_by_haplotype = {}
            
            for hap_name, path_name in sample_haplotypes.items():
                logging.info(f"Processing haplotype: {hap_name} (path: {path_name})")
                
                # Build path sequence
                path_segments = paths[path_name]['segments']
                hap_seq, segment_offsets = build_path_sequence(segments, path_segments)
                
                # Analyze haplotype differences
                hap_effects = analyze_haplotype_differences(
                    ref_seq, 
                    hap_seq, 
                    features, 
                    segment_offsets, 
                    variant_segments,
                    segments
                )
                
                feature_effects_by_haplotype[hap_name] = hap_effects
            
            # Analyze homozygous vs heterozygous effects and compound heterozygous effects
            zygosity_effects = analyze_sample_haplotypes(
                feature_effects_by_haplotype, 
                features,
                feature_by_id,  # Add feature hierarchy info
                children_by_parent  # Add feature hierarchy info
            )
            
            # Generate report for this sample
            outfile = f"{outfile_prefix}_{sample_name}.txt" if outfile_prefix else None
            generate_variant_effect_report([], variants, outfile, sample_name, zygosity_effects)
            
            sample_reports[sample_name] = outfile
        else:
            # Process single haplotype
            logging.info(f"Only one haplotype found for sample {sample_name}, processing as single path")
            
            # Get the first path for this sample
            path_name = path_list[0]
            path_segments = paths[path_name]['segments']
            
            # Build path sequence
            alt_seq, segment_offsets = build_path_sequence(segments, path_segments)
            
            # Analyze differences
            feature_effects = analyze_haplotype_differences(
                ref_seq, 
                alt_seq, 
                features, 
                segment_offsets, 
                variant_segments,
                segments
            )
            
            # Generate report
            outfile = f"{outfile_prefix}_{sample_name}.txt" if outfile_prefix else None
            generate_variant_effect_report(feature_effects, variants, outfile, sample_name)
            
            sample_reports[sample_name] = outfile
    
    return sample_reports


def infer_variant_types(variants, segments):
    """
    Infer variant types based on segment information and length changes.
    
    Args:
        variants: List of variant dictionaries
        segments: Dictionary of segment data
        
    Returns:
        Updated list of variants with inferred types
    """
    for variant in variants:
        # Skip if type is already defined and not UNKNOWN
        if variant.get('type') and variant.get('type') != 'UNKNOWN':
            continue
            
        length_change = variant.get('length_change', 0)
        var_segments = variant.get('segments', [])
        
        # Look for variant type hints in segment names
        for seg_id in var_segments:
            if seg_id not in segments:
                continue
                
            # Check segment name for type hints
            if "SNP" in seg_id or "SNV" in seg_id:
                variant['type'] = "SNP"
                break
            elif "DEL" in seg_id:
                variant['type'] = "DEL"
                break
            elif "INS" in seg_id:
                variant['type'] = "INS"
                break
            elif "DUP" in seg_id:
                variant['type'] = "DUP"
                break
            elif "INV" in seg_id:
                variant['type'] = "INV"
                break
        
        # If no type found from names, infer from length change
        if (not variant.get('type') or variant.get('type') == 'UNKNOWN') and length_change != 0:
            if length_change > 0:
                # Check for duplication - if the variant segment is much longer
                for seg_id in var_segments:
                    if seg_id in segments:
                        seg_data = segments[seg_id]
                        orig_seg_id = seg_data.get('original_segment')
                        
                        if orig_seg_id and orig_seg_id in segments:
                            var_seq = seg_data.get('sequence', '')
                            orig_seq = segments[orig_seg_id].get('sequence', '')
                            
                            # If the variant contains the original sequence multiple times
                            if orig_seq and len(orig_seq) > 0 and var_seq.count(orig_seq) > 1:
                                variant['type'] = "DUP"
                                break
                
                # Default to insertion if not a duplication
                if not variant.get('type') or variant.get('type') == 'UNKNOWN':
                    variant['type'] = "INS"
            elif length_change < 0:
                variant['type'] = "DEL"
        
        # Default to SNP for small zero-length-change variants without a type
        if (not variant.get('type') or variant.get('type') == 'UNKNOWN') and length_change == 0:
            variant['type'] = "SNP"
    
    return variants

def get_shex_schema():
    """
    Return the ShEx schema for validating variant effect RDF data.
    
    The ShEx (Shape Expressions) schema defines the expected structure of the RDF data
    generated by the variant effect analyzer. It can be used with ShEx validation tools
    to verify that the RDF output conforms to the expected structure and semantics.
    
    The schema validates the following key components:
    
    1. Report structure with variants, features, and effects
    2. Variant information including type, length change, and position
    3. Feature information with proper genomic locations
    4. Effect information including type, affected feature, and zygosity
    5. Sample and haplotype information for complex genetic analyses
    6. Compound heterozygous effect representation
    
    Returns:
        A string containing the ShEx schema
    """
    return """PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
PREFIX faldo: <http://biohackathon.org/resource/faldo#>
PREFIX so: <http://purl.obolibrary.org/obo/SO_>
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX dc: <http://purl.org/dc/terms/>
PREFIX sio: <http://semanticscience.org/resource/>
PREFIX gff: <http://www.sequenceontology.org/gff3/1.0/>

# Define base URIs for our namespaces
PREFIX base: <http://example.org/genomics/>
PREFIX variant: <http://example.org/genomics/variant/>
PREFIX sample: <http://example.org/genomics/sample/>
PREFIX feature: <http://example.org/genomics/feature/>
PREFIX effect: <http://example.org/genomics/effect/>
PREFIX position: <http://example.org/genomics/position/>

# Define the report shape
<ReportShape> {
  rdf:type [sio:SequenceVariantAnalysisReport] ;
  rdfs:label xsd:string? ;
  dc:created xsd:dateTime ;
  variant:hasVariant @<VariantShape>+ ;
  feature:hasFeature @<FeatureShape>+ ;
  effect:hasEffect @<EffectShape>* ;
  sio:refersTo @<SampleShape>*
}

# Define the variant shape
<VariantShape> {
  rdf:type [so:sequence_variant] ;
  rdfs:label xsd:string ;
  variant:variantType IRI ;
  variant:lengthChange xsd:integer ;
  faldo:location @<RegionShape>? ;
  variant:affectsSegment @<SegmentShape>*
}

# Define the feature shape
<FeatureShape> {
  rdf:type IRI ;  # Any GFF3 feature type
  rdfs:label xsd:string ;
  faldo:location @<RegionShape> ;
  faldo:strand [faldo:ForwardStrand faldo:ReverseStrand faldo:BothStrand]?
}

# Define the effect shape
<EffectShape> {
  rdf:type [effect:VariantEffect] ;
  effect:effectType IRI ;
  effect:affectsFeature @<FeatureShape> ;
  variant:zygosity xsd:string? ;
  effect:inSample @<SampleShape>? ;
  effect:inHaplotype @<HaplotypeShape>* ;
  dc:description xsd:string? ;
  variant:referenceSequence xsd:string? ;
  variant:alternateSequence xsd:string? ;
  effect:causedBy @<VariantShape>* ;
  effect:variantInHaplotype @<HaplotypeVariantAssociationShape>*
}

# Define the sample shape
<SampleShape> {
  rdf:type [sio:Sample] ;
  rdfs:label xsd:string ;
  base:hasHaplotype @<HaplotypeShape>*
}

# Define the haplotype shape
<HaplotypeShape> {
  rdf:type [base:Haplotype] ;
  rdfs:label xsd:string
}

# Define the region shape (for genomic locations)
<RegionShape> {
  rdf:type [faldo:Region] ;
  faldo:begin @<PositionShape> ;
  faldo:end @<PositionShape>
}

# Define the position shape
<PositionShape> {
  rdf:type [faldo:ExactPosition] ;
  faldo:position xsd:integer
}

# Define the segment shape
<SegmentShape> {
  rdf:type [base:Segment] ;
  rdfs:label xsd:string
}

# Define the haplotype-variant association shape
<HaplotypeVariantAssociationShape> {
  rdf:type [effect:HaplotypeVariantAssociation] ;
  effect:haplotype @<HaplotypeShape> ;
  effect:variant @<VariantShape>
}
"""


def main():
    """Main function to run the variant effect report generation."""
    parser = argparse.ArgumentParser(description='Generate a variant effect report based on GFA and GFF3 files.')
    parser.add_argument('gfa_file', help='Input GFA file with variants')
    parser.add_argument('gff_file', help='Input GFF3 file with feature annotations')
    
    # Output options
    output_group = parser.add_argument_group('Output Options')
    output_group.add_argument('--output', '-o', help='Output file (default: stdout)')
    output_group.add_argument('--format', '-f', choices=['text', 'rdf'], default='text',
                         help='Output format: text (default) or RDF')
    output_group.add_argument('--rdf-format', choices=['turtle', 'n3', 'xml', 'json-ld', 'ntriples'], default='n3',
                         help='RDF serialization format (default: n3)')
    output_group.add_argument('--base-uri', default='http://example.org/genomics/',
                         help='Base URI for RDF output (default: http://example.org/genomics/)')
    output_group.add_argument('--consolidated', action='store_true',
                         help='Output a single consolidated RDF file for all samples')
    
    # Sample options
    sample_group = parser.add_argument_group('Sample Options')
    sample_group.add_argument('--sample-reports', action='store_true', 
                         help='Generate individual reports for each sample')
    sample_group.add_argument('--output-prefix', help='Prefix for sample report files')
    sample_group.add_argument('--samples', help='Comma-separated list of sample names to process (default: all)')
    
    # Analysis options
    analysis_group = parser.add_argument_group('Analysis Options')
    analysis_group.add_argument('--use-alignment', action='store_true',
                         help='Use sequence alignment for variant effect analysis')
    
    # Debug and logging options
    debug_group = parser.add_argument_group('Debug Options')
    debug_group.add_argument('--debug', action='store_true', help='Enable debug output')
    debug_group.add_argument('--verbose', action='store_true', help='Enable verbose output without full debug')
    debug_group.add_argument('--log-file', help='Write log to this file')
    debug_group.add_argument('--save-shex', help='Save ShEx validation schema to the specified file')
    
    args = parser.parse_args()
    
    # If --save-shex is specified, save the ShEx schema and exit
    if args.save_shex:
        with open(args.save_shex, 'w') as f:
            f.write(get_shex_schema())
        logging.info(f"ShEx schema saved to {args.save_shex}")
        return 0
    
    # Setup logging
    logger = setup_logging(debug=args.debug, log_file=args.log_file, verbose=args.verbose)
    
    try:
        # Parse input files
        segments, links, paths, variants, samples, haplotypes = parse_gfa(args.gfa_file)
        features, feature_by_id, children_by_parent = parse_gff3(args.gff_file)
        
        # Mark first CDS in each gene for start codon analysis
        features = identify_first_cds_in_genes(features, feature_by_id, children_by_parent)
        
        # Calculate accurate length changes for variants
        variants = calculate_variant_length_changes(variants, segments)
        
        # Infer variant types based on segment information and length changes
        variants = infer_variant_types(variants, segments)
        
        # Analyze repetitive sequences in variants and ensure variant info is complete
        for variant in variants:
            # Analyze repetitive sequences for INS and DUP variants
            if variant.get('type') in ['DUP', 'INS']:
                repeat_info = identify_repeat_sequences(variant, segments)
                if repeat_info:
                    variant['repeat_info'] = repeat_info
            
            # Ensure we have segment associations
            if 'segments' not in variant:
                variant['segments'] = []
                
            # Make sure all variants have position information
            if 'pos' not in variant or not variant['pos']:
                # Try to infer position from segment info
                for seg_id in variant.get('segments', []):
                    if seg_id in segments and segments[seg_id].get('original_segment'):
                        orig_seg = segments[seg_id]['original_segment']
                        if orig_seg and orig_seg.startswith('S'):
                            try:
                                pos = int(orig_seg[1:])
                                variant['pos'] = pos
                                variant['end'] = pos + len(segments[seg_id]['sequence']) - 1
                                break
                            except ValueError:
                                pass
        
        # Log variant information
        for var in variants:
            logging.info(f"Variant {var.get('id', 'unknown')}: {var.get('type', 'UNKNOWN')}, length change: {var.get('length_change', 0)} bp")
        
        # Filter samples if requested
        if args.samples and args.sample_reports:
            requested_samples = [s.strip() for s in args.samples.split(',')]
            filtered_samples = {name: paths for name, paths in samples.items() if name in requested_samples}
            if not filtered_samples:
                logging.warning(f"None of the requested samples found in GFA")
            samples = filtered_samples
        
        # Select the appropriate analysis function based on the --use-alignment flag
        if args.use_alignment:
            logging.info("Using sequence alignment for variant effect analysis")
            analyze_function = analyze_haplotype_differences_with_alignment
        else:
            logging.info("Using position-based approach for variant effect analysis")
            analyze_function = analyze_haplotype_differences
        
        # Handle consolidated output if requested (only for RDF format)
        if args.consolidated and args.format == 'rdf':
            logging.info("Generating consolidated RDF report for all samples")
            
            # Create consolidated graph
            consolidated_graph = create_consolidated_rdf_report(
                variants,
                features,
                samples,
                haplotypes,
                paths,
                segments,
                feature_by_id,
                children_by_parent,
                base_uri=args.base_uri
            )
            
            # Output consolidated graph
            output_rdf_report(consolidated_graph, args.output, args.rdf_format)
            
            # Skip other report generation
            return 0
            
        # Generate sample-specific reports if requested
        if args.sample_reports:
            if samples:
                sample_reports = {}
                
                for sample_name, sample_paths in samples.items():
                    logging.info(f"Processing sample: {sample_name}")
                    
                    # Get reference path
                    ref_path_name = 'REF'
                    if ref_path_name not in paths:
                        logging.error("REF path not found in GFA")
                        continue
                    
                    # Build reference sequence
                    ref_path_segments = paths[ref_path_name]['segments']
                    ref_seq, _ = build_path_sequence(segments, ref_path_segments)
                    
                    # Determine if we have phased haplotypes for this sample
                    sample_haplotypes = haplotypes.get(sample_name, {})
                    has_phased_haplotypes = len(sample_haplotypes) >= 2
                    
                    # Collect variant segments
                    variant_segments = {}
                    for seg_id, seg_data in segments.items():
                        if seg_data.get('variant_id'):
                            variant_segments[seg_id] = seg_data
                    
                    # Prepare output file path
                    if args.output_prefix:
                        if args.format == 'text':
                            outfile = f"{args.output_prefix}_{sample_name}.txt"
                        else:  # RDF format
                            outfile = f"{args.output_prefix}_{sample_name}.{args.rdf_format}"
                    else:
                        outfile = None  # Use stdout
                    
                    if has_phased_haplotypes:
                        # Process each haplotype
                        feature_effects_by_haplotype = {}
                        
                        for hap_name, path_name in sample_haplotypes.items():
                            # Build path sequence
                            path_segments = paths[path_name]['segments']
                            hap_seq, segment_offsets = build_path_sequence(segments, path_segments)
                            
                            # Analyze differences using the selected analysis method
                            hap_effects = analyze_function(
                                ref_seq, 
                                hap_seq, 
                                features, 
                                segment_offsets, 
                                variant_segments,
                                segments
                            )
                            
                            feature_effects_by_haplotype[hap_name] = hap_effects
                        
                        # Analyze homozygous vs heterozygous effects
                        zygosity_effects = analyze_sample_haplotypes(
                            feature_effects_by_haplotype, 
                            features,
                            feature_by_id,
                            children_by_parent
                        )
                        
                        # Generate output
                        if args.format == 'text':
                            # Generate text report
                            generate_variant_effect_report([], variants, outfile, sample_name, zygosity_effects)
                            sample_reports[sample_name] = outfile
                        else:  # RDF format
                            # Create RDF graph
                            graph = create_rdf_variant_report(
                                variants,
                                features,
                                feature_effects=None,
                                sample_name=sample_name,
                                sample_effects=zygosity_effects,
                                base_uri=args.base_uri
                            )
                            # Output RDF graph
                            output_rdf_report(graph, outfile, args.rdf_format)
                            sample_reports[sample_name] = outfile
                    else:
                        # Process single haplotype
                        path_name = sample_paths[0]
                        path_segments = paths[path_name]['segments']
                        
                        # Build path sequence
                        alt_seq, segment_offsets = build_path_sequence(segments, path_segments)
                        
                        # Analyze differences using the selected analysis method
                        feature_effects = analyze_function(
                            ref_seq, 
                            alt_seq, 
                            features, 
                            segment_offsets, 
                            variant_segments,
                            segments
                        )
                        
                        # Generate output
                        if args.format == 'text':
                            # Generate text report
                            generate_variant_effect_report(feature_effects, variants, outfile, sample_name)
                            sample_reports[sample_name] = outfile
                        else:  # RDF format
                            # Create RDF graph
                            graph = create_rdf_variant_report(
                                variants,
                                features,
                                feature_effects=feature_effects,
                                sample_name=sample_name,
                                sample_effects=None,
                                base_uri=args.base_uri
                            )
                            # Output RDF graph
                            output_rdf_report(graph, outfile, args.rdf_format)
                            sample_reports[sample_name] = outfile
                
                logging.info(f"Generated reports for {len(sample_reports)} samples")
            else:
                logging.warning("No samples found in GFA, cannot generate sample reports")
        
        # Generate main report if requested
        if not args.sample_reports and not args.consolidated:
            # Get reference path
            ref_path_name = 'REF'
            if ref_path_name not in paths:
                logging.error("REF path not found in GFA")
                return 1
            
            # Build reference sequence
            ref_path_segments = paths[ref_path_name]['segments']
            ref_seq, ref_offsets = build_path_sequence(segments, ref_path_segments)
            
            # Collect variant segments
            variant_segments = {}
            for seg_id, seg_data in segments.items():
                if seg_data.get('variant_id'):
                    variant_segments[seg_id] = seg_data
            
            # Process each sample or first available sample if none specified
            processed_sample = False
            
            for sample_name, sample_haplotypes in haplotypes.items():
                if not args.samples or sample_name in args.samples.split(','):
                    logging.info(f"Processing sample: {sample_name}")
                    
                    # Process all haplotypes for compound heterozygous detection
                    if len(sample_haplotypes) >= 2:
                        feature_effects_by_haplotype = {}
                        
                        for hap_name, path_name in sample_haplotypes.items():
                            # Build path sequence
                            path_segments = paths[path_name]['segments']
                            hap_seq, segment_offsets = build_path_sequence(segments, path_segments)
                            
                            # Analyze differences using the selected analysis method
                            hap_effects = analyze_function(
                                ref_seq, 
                                hap_seq, 
                                features, 
                                segment_offsets, 
                                variant_segments,
                                segments
                            )
                            
                            feature_effects_by_haplotype[hap_name] = hap_effects
                        
                        # Analyze homozygous vs heterozygous effects
                        zygosity_effects = analyze_sample_haplotypes(
                            feature_effects_by_haplotype, 
                            features,
                            feature_by_id,
                            children_by_parent
                        )
                        
                        # Generate output
                        if args.format == 'text':
                            # Generate text report
                            generate_variant_effect_report([], variants, args.output, sample_name, zygosity_effects)
                        else:  # RDF format
                            # Create RDF graph
                            graph = create_rdf_variant_report(
                                variants,
                                features,
                                feature_effects=None,
                                sample_name=sample_name,
                                sample_effects=zygosity_effects,
                                base_uri=args.base_uri
                            )
                            # Output RDF graph
                            output_rdf_report(graph, args.output, args.rdf_format)
                    else:
                        # Process single haplotype
                        hap_name, path_name = next(iter(sample_haplotypes.items()))
                        
                        # Build path sequence
                        path_segments = paths[path_name]['segments']
                        alt_seq, segment_offsets = build_path_sequence(segments, path_segments)
                        
                        # Analyze differences using the selected analysis method
                        feature_effects = analyze_function(
                            ref_seq, 
                            alt_seq, 
                            features, 
                            segment_offsets, 
                            variant_segments,
                            segments
                        )
                        
                        # Generate output
                        if args.format == 'text':
                            # Generate text report
                            generate_variant_effect_report(feature_effects, variants, args.output, sample_name)
                        else:  # RDF format
                            # Create RDF graph
                            graph = create_rdf_variant_report(
                                variants,
                                features,
                                feature_effects=feature_effects,
                                sample_name=sample_name,
                                sample_effects=None,
                                base_uri=args.base_uri
                            )
                            # Output RDF graph
                            output_rdf_report(graph, args.output, args.rdf_format)
                    
                    processed_sample = True
                    break
            
            # If no samples processed, use the first available path
            if not processed_sample:
                for path_name, path_data in paths.items():
                    if path_name != 'REF' and 'segments' in path_data:
                        # Build path sequence
                        path_segments = path_data['segments']
                        alt_seq, segment_offsets = build_path_sequence(segments, path_segments)
                        
                        # Analyze differences using the selected analysis method
                        feature_effects = analyze_function(
                            ref_seq, 
                            alt_seq, 
                            features, 
                            segment_offsets, 
                            variant_segments,
                            segments
                        )
                        
                        # Generate output
                        if args.format == 'text':
                            # Generate text report
                            generate_variant_effect_report(feature_effects, variants, args.output)
                        else:  # RDF format
                            # Create RDF graph
                            graph = create_rdf_variant_report(
                                variants,
                                features,
                                feature_effects=feature_effects,
                                sample_name=None,
                                sample_effects=None,
                                base_uri=args.base_uri
                            )
                            # Output RDF graph
                            output_rdf_report(graph, args.output, args.rdf_format)
                        break
        
        return 0
        
    except Exception as e:
        logging.error(f"Error during report generation: {e}")
        if args.debug:
            import traceback
            logging.error(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())

