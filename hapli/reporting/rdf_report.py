"""
RDF report generation for variant effects.
"""

import sys
import logging
import uuid
from datetime import datetime
import rdflib
from rdflib import Graph, Literal, BNode, Namespace, RDF, URIRef, XSD
from rdflib.namespace import FOAF, RDFS, SKOS

from hapli.models.genomic import GenomicNamespaces

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

# Import needed for consolidated report
from hapli.analysis.sequence_analysis import build_path_sequence, analyze_haplotype_differences
from hapli.analysis.feature_analysis import analyze_sample_haplotypes
