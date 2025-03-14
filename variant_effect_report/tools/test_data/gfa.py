"""
GFA file generation utilities using gfapy.
"""

import os
import time
import logging
import random
from typing import Dict, List, Tuple, Optional

try:
    import gfapy
except ImportError:
    logging.error("gfapy module not found. Please install it with: pip install gfapy")
    raise

from .sequence import SequenceHandler


class GFAGenerator:
    """Generate GFA files from reference sequences"""
    
    def __init__(self, ref_seq: str, node_size: int = 100, sequence_id: str = "REF", seed: Optional[int] = None):
        """
        Initialize the generator
        
        Args:
            ref_seq: Reference sequence
            node_size: Average size of nodes
            sequence_id: ID for the reference sequence
            seed: Random seed for reproducibility
        """
        if not ref_seq:
            raise ValueError("Reference sequence cannot be empty")
            
        self.ref_seq = ref_seq
        self.seq_length = len(ref_seq)
        self.node_size = max(10, node_size)  # Minimum node size of 10
        self.sequence_id = sequence_id
        
        # Set random seed
        if seed is not None:
            random.seed(seed)
        
        # Create nodes and paths
        self.nodes = self._create_nodes()
        self.links = self._create_links()
        self.paths = {self.sequence_id: self._create_reference_path()}
        
        # Initialize sequence handler
        self.seq_handler = SequenceHandler(seed)
    
    def _create_nodes(self) -> Dict[str, str]:
        """
        Create nodes by segmenting the reference sequence
        
        Returns:
            Dictionary mapping node IDs to sequences
        """
        nodes = {}
        pos = 0
        node_id = 1
        
        while pos < self.seq_length:
            # Vary node size slightly for realism
            size_variation = min(self.node_size // 5, 10)  # Up to 20% variation or 10bp
            size = min(
                self.node_size + random.randint(-size_variation, size_variation), 
                self.seq_length - pos
            )
            
            if size <= 0:
                break
                
            node_seq = self.ref_seq[pos:pos+size]
            nodes[f"S{node_id}"] = node_seq
            pos += size
            node_id += 1
            
        return nodes
    
    def _create_links(self) -> List[Tuple[str, str, str, str, str]]:
        """
        Create links between consecutive nodes
        
        Returns:
            List of link tuples (from_node, from_orient, to_node, to_orient, cigar)
        """
        links = []
        node_ids = list(self.nodes.keys())
        
        for i in range(len(node_ids) - 1):
            # Link each node to the next one
            from_node = node_ids[i]
            to_node = node_ids[i+1]
            links.append((from_node, '+', to_node, '+', '0M'))
            
        return links
    
    def _create_reference_path(self) -> str:
        """
        Create reference path through the nodes
        
        Returns:
            Path string in GFA format
        """
        node_ids = list(self.nodes.keys())
        # Create a reference path that uses all nodes in order
        return ','.join([f"{node}+" for node in node_ids])
    
    def export_gfa(self, filename: str, add_metadata: bool = True) -> None:
        """
        Export the graph to a GFA file using gfapy
        
        Args:
            filename: Output filename
            add_metadata: Whether to add metadata to the header
            
        Raises:
            IOError: If there's an error writing the file
        """
        try:
            # Create GFA object
            gfa = gfapy.Gfa()
            
            # Add header with metadata
            if add_metadata:
                # Use current timestamp as integer (seconds since epoch)
                timestamp = int(time.time())
                gfa.add_line(f"H\tVN:Z:1.0\tTS:i:{timestamp}\tPG:Z:testdata_generator")
            else:
                gfa.add_line("H\tVN:Z:1.0")
            
            # Add segments
            for node_id, sequence in self.nodes.items():
                gfa.add_line(f"S\t{node_id}\t{sequence}")
            
            # Add links
            for from_node, from_orient, to_node, to_orient, overlap in self.links:
                gfa.add_line(f"L\t{from_node}\t{from_orient}\t{to_node}\t{to_orient}\t{overlap}")
            
            # Add paths
            for path_id, segment_list in self.paths.items():
                gfa.add_line(f"P\t{path_id}\t{segment_list}\t*")
            
            # Write to file
            gfa.to_file(filename)
                
        except Exception as e:
            logging.error(f"Error writing GFA file: {str(e)}")
            raise IOError(f"Failed to write GFA file: {str(e)}")
    
    def export_fasta(self, filename: str) -> None:
        """
        Export the reference sequence to a FASTA file
        
        Args:
            filename: Output filename
            
        Raises:
            IOError: If there's an error writing the file
        """
        try:
            self.seq_handler.save_fasta(self.ref_seq, self.sequence_id, filename)
        except Exception as e:
            logging.error(f"Error writing FASTA file: {str(e)}")
            raise IOError(f"Failed to write FASTA file: {str(e)}")
    
    def generate_vcf_variants(self, num_variants: int = 5, var_types: Optional[List[str]] = None, 
                            size_ranges: Optional[Dict[str, Tuple[int, int]]] = None, 
                            hgvs_prefix: str = "NC_000001.11") -> List[Dict]:
        """
        Generate VCF variants
        
        Args:
            num_variants: Number of variants to generate
            var_types: List of variant types to include (SNP, DEL, INS, DUP, INV)
            size_ranges: Dictionary with min/max size for each variant type
            hgvs_prefix: Prefix for HGVS notation
            
        Returns:
            List of variant dictionaries
        """
        from .variants import VariantGenerator
        
        # Default variant types if not specified
        if var_types is None:
            var_types = ['SNP', 'DEL', 'INS', 'DUP', 'INV']
        
        # Create variant generator
        variant_generator = VariantGenerator(self.ref_seq)
        
        # Generate variants
        return variant_generator.generate_variants(
            num_variants=num_variants,
            var_types=var_types,
            size_ranges=size_ranges,
            hgvs_prefix=hgvs_prefix
        )
    
    def export_vcf(self, variants: List[Dict], filename: str, contig_id: str = "1", 
                 samples: Optional[List[str]] = None, phased: bool = False) -> None:
        """
        Export variants to a VCF file
        
        Args:
            variants: List of variant dictionaries
            filename: Output VCF filename
            contig_id: Contig ID to use in VCF
            samples: List of sample names for multi-sample VCF (None = no samples)
            phased: Whether to generate phased genotypes (|) or unphased (/)
            
        Raises:
            IOError: If there's an error writing the file
        """
        from .vcf import VCFGenerator
        
        # Create VCF generator
        vcf_generator = VCFGenerator(self.seq_length, contig_id)
        
        # Export VCF
        vcf_generator.export_vcf(variants, filename, samples, phased)
