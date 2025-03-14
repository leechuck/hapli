"""
Convert variant effect data to RDF.
"""

import uuid
from datetime import datetime
from typing import Dict, List, Optional, Any

from rdflib import Graph, Literal, BNode, Namespace, RDF, URIRef, XSD

from variant_effect_report.rdf.namespaces import GenomicNamespaces
from variant_effect_report.core.models import Feature, Variant, FeatureEffect


class RDFConverter:
    """Convert variant effect data to RDF."""
    
    def __init__(self, base_uri="http://example.org/genomics/"):
        """
        Initialize the RDF converter.
        
        Args:
            base_uri: Base URI for the RDF graph
        """
        self.base_uri = base_uri
        self.ns = GenomicNamespaces(base_uri)
    
    def create_variant_report(self, variants: List[Dict], features: List[Feature], 
                             feature_effects: Optional[List[FeatureEffect]] = None, 
                             sample_name: Optional[str] = None, 
                             sample_effects: Optional[Dict] = None) -> Graph:
        """
        Create an RDF graph of variant effects.
        
        Args:
            variants: List of variant dictionaries
            features: List of feature objects
            feature_effects: List of feature effect objects (optional)
            sample_name: Name of the sample being analyzed (optional)
            sample_effects: Dictionary of sample-specific effects (optional)
            
        Returns:
            RDF graph object
        """
        # Initialize graph and namespaces
        g = Graph()
        self.ns.bind_to_graph(g)
        
        # Create report node
        report_id = str(uuid.uuid4())
        report_uri = URIRef(f"{self.base_uri}report/{report_id}")
        g.add((report_uri, RDF.type, self.ns.sio["SequenceVariantAnalysisReport"]))
        g.add((report_uri, self.ns.dc.created, Literal(datetime.now().isoformat(), datatype=XSD.dateTime)))
        
        # Add sample information if available
        if sample_name:
            sample_uri = URIRef(f"{self.ns.sample}{sample_name.replace(' ', '_')}")
            g.add((sample_uri, RDF.type, self.ns.sio.Sample))
            g.add((sample_uri, RDFS.label, Literal(sample_name)))
            g.add((report_uri, self.ns.sio.refersTo, sample_uri))
        
        # Process variants
        for variant in variants:
            var_id = variant.get('id', f"unknown_{str(uuid.uuid4())}")
            var_type = variant.get('type', 'UNKNOWN')
            
            # Create variant node
            var_uri = URIRef(f"{self.ns.variant}{var_id}")
            g.add((var_uri, RDF.type, self.ns.so.sequence_variant))
            g.add((var_uri, RDFS.label, Literal(var_id)))
            g.add((var_uri, self.ns.variant.variantType, self.ns.get_so_term(var_type)))
            g.add((var_uri, self.ns.variant.lengthChange, Literal(variant.get('length_change', 0), datatype=XSD.integer)))
            
            # Add variant to report
            g.add((report_uri, self.ns.variant.hasVariant, var_uri))
            
            # Add position information if available
            if 'pos' in variant and variant.get('pos', 0) > 0:
                pos_uri = URIRef(f"{self.ns.position}{var_id}")
                g.add((pos_uri, RDF.type, self.ns.faldo.Region))
                g.add((var_uri, self.ns.faldo.location, pos_uri))
                
                # Add start position
                start_uri = URIRef(f"{self.ns.position}{var_id}_start")
                g.add((start_uri, RDF.type, self.ns.faldo.ExactPosition))
                g.add((start_uri, self.ns.faldo.position, Literal(variant.get('pos', 0), datatype=XSD.integer)))
                g.add((pos_uri, self.ns.faldo.begin, start_uri))
                
                # Add end position
                end_uri = URIRef(f"{self.ns.position}{var_id}_end")
                g.add((end_uri, RDF.type, self.ns.faldo.ExactPosition))
                g.add((end_uri, self.ns.faldo.position, Literal(variant.get('end', variant.get('pos', 0)), datatype=XSD.integer)))
                g.add((pos_uri, self.ns.faldo.end, end_uri))
            
            # Add affected segments
            if 'segments' in variant and variant['segments']:
                for seg_id in variant['segments']:
                    seg_uri = URIRef(f"{self.ns.base}segment/{seg_id}")
                    g.add((seg_uri, RDF.type, self.ns.base.Segment))
                    g.add((seg_uri, RDFS.label, Literal(seg_id)))
                    g.add((var_uri, self.ns.variant.affectsSegment, seg_uri))
        
        # Process features
        for feature in features:
            feature_id = feature.attributes.get('ID', f"unknown_{str(uuid.uuid4())}")
            feature_name = feature.attributes.get('Name', feature_id)
            feature_type = feature.type
            
            # Create feature node
            feature_uri = URIRef(f"{self.ns.feature}{feature_id}")
            g.add((feature_uri, RDF.type, URIRef(f"{self.ns.gff}{feature_type}")))
            g.add((feature_uri, RDFS.label, Literal(feature_name)))
            
            # Add feature to report
            g.add((report_uri, self.ns.feature.hasFeature, feature_uri))
            
            # Add location information
            loc_uri = URIRef(f"{self.ns.position}feature_{feature_id}")
            g.add((loc_uri, RDF.type, self.ns.faldo.Region))
            g.add((feature_uri, self.ns.faldo.location, loc_uri))
            
            # Add start position
            start_uri = URIRef(f"{self.ns.position}feature_{feature_id}_start")
            g.add((start_uri, RDF.type, self.ns.faldo.ExactPosition))
            g.add((start_uri, self.ns.faldo.position, Literal(feature.start, datatype=XSD.integer)))
            g.add((loc_uri, self.ns.faldo.begin, start_uri))
            
            # Add end position
            end_uri = URIRef(f"{self.ns.position}feature_{feature_id}_end")
            g.add((end_uri, RDF.type, self.ns.faldo.ExactPosition))
            g.add((end_uri, self.ns.faldo.position, Literal(feature.end, datatype=XSD.integer)))
            g.add((loc_uri, self.ns.faldo.end, end_uri))
            
            # Add strand information
            if feature.strand == '+':
                g.add((feature_uri, self.ns.faldo.strand, self.ns.faldo.ForwardStrand))
            elif feature.strand == '-':
                g.add((feature_uri, self.ns.faldo.strand, self.ns.faldo.ReverseStrand))
            else:
                g.add((feature_uri, self.ns.faldo.strand, self.ns.faldo.BothStrand))
        
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
            if not effect.effects or effect.effects == ['no_change']:
                continue  # Skip effects with no change
                
            feature = effect.feature
            feature_id = feature.attributes.get('ID', 'unknown')
            feature_type = feature.type
            feature_uri = URIRef(f"{self.ns.feature}{feature_id}")
            
            # For each effect type
            for effect_type in effect.effects:
                if effect_type == 'no_change':
                    continue
                    
                effect_id = f"{feature_id}_{effect_type}_{str(uuid.uuid4())}"
                effect_uri = URIRef(f"{self.ns.effect}{effect_id}")
                
                # Create effect node
                g.add((effect_uri, RDF.type, self.ns.effect.VariantEffect))
                g.add((effect_uri, self.ns.effect.effectType, self.ns.get_effect_term(effect_type)))
                g.add((effect_uri, self.ns.effect.affectsFeature, feature_uri))
                
                # Add effect to report
                g.add((report_uri, self.ns.effect.hasEffect, effect_uri))
                
                # Add detailed description if available
                if effect_type in effect.details:
                    detail = effect.details[effect_type]
                    if isinstance(detail, dict):
                        desc = ", ".join(f"{k}={v}" for k, v in detail.items())
                    elif isinstance(detail, list):
                        desc = ", ".join(str(d) for d in detail)
                    else:
                        desc = str(detail)
                    g.add((effect_uri, self.ns.dc.description, Literal(desc)))
                
                # Add zygosity if available
                if effect.zygosity:
                    g.add((effect_uri, self.ns.variant.zygosity, Literal(effect.zygosity)))
                
                # Add sequence information
                if effect.ref_feature_seq and effect.alt_feature_seq:
                    g.add((effect_uri, self.ns.variant.referenceSequence, Literal(effect.ref_feature_seq)))
                    g.add((effect_uri, self.ns.variant.alternateSequence, Literal(effect.alt_feature_seq)))
                
                # Add variants causing this effect
                if effect.variants:
                    for var in effect.variants:
                        var_id = var.get('id', 'unknown')
                        var_uri = URIRef(f"{self.ns.variant}{var_id}")
                        g.add((effect_uri, self.ns.effect.causedBy, var_uri))
        
        return g
    
    def output_report(self, graph: Graph, output=None, format='turtle'):
        """
        Output an RDF graph in the specified format.
        
        Args:
            graph: RDF graph object
            output: Output file (default: stdout)
            format: RDF serialization format (default: turtle)
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
            import sys
            output_str = graph.serialize(format=rdf_format)
            sys.stdout.write(output_str.decode('utf-8') if isinstance(output_str, bytes) else output_str)
