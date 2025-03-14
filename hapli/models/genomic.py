"""
Data models for genomic entities like variants, features, and effects.
"""

import rdflib
from rdflib import Graph, Literal, BNode, Namespace, RDF, URIRef, XSD
from rdflib.namespace import FOAF, RDFS, SKOS
import uuid
from datetime import datetime

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
