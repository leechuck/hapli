"""
Utility functions and classes for generating test data for variant effect analysis.
"""

import random
import os
import gzip
from typing import List, Tuple, Dict, Optional
from Bio import SeqIO
from Bio.Seq import Seq

class SimpleGFAGenerator:
    """Generates a simple GFA file with a reference path and node structure"""
    
    def __init__(self, ref_seq=None, seq_length=1000, node_size=100, seed=42, sequence_id="REF"):
        """
        Initialize the generator
        
        Args:
            ref_seq: Reference sequence (string) to use. If None, a random sequence is generated
            seq_length: Total length of the reference sequence if generating randomly
            node_size: Average size of nodes
            seed: Random seed for reproducibility
            sequence_id: ID for the reference sequence
        """
        random.seed(seed)
        self.sequence_id = sequence_id
        
        # Use provided reference or generate random sequence
        if ref_seq is not None:
            self.ref_seq = ref_seq
            self.seq_length = len(ref_seq)
        else:
            self.seq_length = seq_length
            self.ref_seq = self._generate_random_sequence(seq_length)
        
        # Create nodes and paths
        self.node_size = node_size
        self.nodes = self._create_nodes()
        self.links = self._create_links()
        self.paths = self._create_paths()
        
    def _generate_random_sequence(self, length: int) -> str:
        """Generate a random DNA sequence of specified length"""
        bases = ['A', 'C', 'G', 'T']
        return ''.join(random.choice(bases) for _ in range(length))
    
    def _create_nodes(self) -> Dict[str, str]:
        """Create nodes by segmenting the reference sequence"""
        nodes = {}
        pos = 0
        node_id = 1
        
        while pos < self.seq_length:
            # Vary node size slightly for realism
            size = min(self.node_size + random.randint(-10, 10), self.seq_length - pos)
            if size <= 0:
                break
                
            node_seq = self.ref_seq[pos:pos+size]
            nodes[f"S{node_id}"] = node_seq
            pos += size
            node_id += 1
            
        return nodes
    
    def _create_links(self) -> List[Tuple[str, str, str, str, str]]:
        """Create links between consecutive nodes"""
        links = []
        node_ids = list(self.nodes.keys())
        
        for i in range(len(node_ids) - 1):
            # Link each node to the next one
            from_node = node_ids[i]
            to_node = node_ids[i+1]
            links.append((from_node, '+', to_node, '+', '0M'))
            
        return links
    
    def _create_paths(self) -> Dict[str, str]:
        """Create reference path through the nodes"""
        node_ids = list(self.nodes.keys())
        # Create a reference path that uses all nodes in order
        ref_path = ','.join([f"{node}+" for node in node_ids])
        
        return {self.sequence_id: ref_path}
    
    def export_gfa(self, filename: str) -> None:
        """Export the graph to a GFA file"""
        with open(filename, 'w') as f:
            # Write header
            f.write("H\tVN:Z:1.0\n")
            
            # Write nodes
            for node_id, sequence in self.nodes.items():
                f.write(f"S\t{node_id}\t{sequence}\n")
            
            # Write links
            for from_node, from_orient, to_node, to_orient, overlap in self.links:
                f.write(f"L\t{from_node}\t{from_orient}\t{to_node}\t{to_orient}\t{overlap}\n")
            
            # Write paths
            for path_id, segment_list in self.paths.items():
                f.write(f"P\t{path_id}\t{segment_list}\t*\n")
    
    def export_fasta(self, filename: str) -> None:
        """Export the reference sequence to a FASTA file"""
        with open(filename, 'w') as f:
            f.write(f">{self.sequence_id}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(self.ref_seq), 80):
                f.write(f"{self.ref_seq[i:i+80]}\n")
    
    def generate_vcf_variants(self, num_variants=5, var_types=None, size_ranges=None, hgvs_prefix="NC_000001.11") -> List[Dict]:
        """
        Generate random structural variants based on the reference sequence
        
        Args:
            num_variants: Number of variants to generate
            var_types: List of variant types to include (SNP, DEL, INS, DUP, INV)
            size_ranges: Dictionary with min/max size for each variant type
                         e.g., {'DEL': (1, 50), 'INS': (1, 20)}
            hgvs_prefix: Prefix to use for HGVS nomenclature
        """
        variants = []
        ref_len = self.seq_length
        
        # Default variant types if none specified
        if var_types is None:
            var_types = ['SNP', 'DEL', 'INS', 'DUP', 'INV']
        
        # Default size ranges
        if size_ranges is None:
            size_ranges = {
                'DEL': (1, 50),
                'INS': (1, 20),
                'DUP': (1, 50),
                'INV': (3, 50)
            }
            
        for _ in range(num_variants):
            # Choose a random position
            pos = random.randint(10, ref_len - 20)
            
            # Choose variant type
            var_type = random.choice(var_types)
            
            if var_type == 'SNP':
                # Single nucleotide polymorphism
                ref_base = self.ref_seq[pos-1]
                alt_bases = [b for b in ['A', 'C', 'G', 'T'] if b != ref_base]
                alt_base = random.choice(alt_bases)
                
                variants.append({
                    'type': 'SNP',
                    'pos': pos,
                    'ref': ref_base,
                    'alt': alt_base,
                    'hgvs': f"{hgvs_prefix}:g.{pos}{ref_base}>{alt_base}"
                })
                
            elif var_type == 'DEL':
                # Deletion
                min_size, max_size = size_ranges.get('DEL', (1, 50))
                del_size = random.randint(min_size, max_size)
                if pos + del_size > ref_len:
                    del_size = ref_len - pos
                
                del_seq = self.ref_seq[pos-1:pos+del_size-1]
                
                variants.append({
                    'type': 'DEL',
                    'pos': pos,
                    'ref': del_seq,
                    'alt': del_seq[0],  # Keep first base
                    'hgvs': f"{hgvs_prefix}:g.{pos}_{pos+del_size-1}del"
                })
                
            elif var_type == 'INS':
                # Insertion
                min_size, max_size = size_ranges.get('INS', (1, 20))
                ins_size = random.randint(min_size, max_size)
                ins_seq = self._generate_random_sequence(ins_size)
                ref_base = self.ref_seq[pos-1]
                
                variants.append({
                    'type': 'INS',
                    'pos': pos,
                    'ref': ref_base,
                    'alt': ref_base + ins_seq,
                    'hgvs': f"{hgvs_prefix}:g.{pos}_{pos+1}ins{ins_seq}"
                })
                
            elif var_type == 'DUP':
                # Duplication
                min_size, max_size = size_ranges.get('DUP', (1, 50))
                dup_size = random.randint(min_size, max_size)
                if pos + dup_size > ref_len:
                    dup_size = ref_len - pos
                
                dup_seq = self.ref_seq[pos-1:pos+dup_size-1]
                
                variants.append({
                    'type': 'DUP',
                    'pos': pos,
                    'ref': dup_seq,
                    'alt': dup_seq + dup_seq,
                    'hgvs': f"{hgvs_prefix}:g.{pos}_{pos+dup_size-1}dup"
                })
                
            elif var_type == 'INV':
                # Inversion
                min_size, max_size = size_ranges.get('INV', (3, 50))
                inv_size = random.randint(min_size, max_size)
                if pos + inv_size > ref_len:
                    inv_size = ref_len - pos
                
                inv_seq = self.ref_seq[pos-1:pos+inv_size-1]
                rev_comp_seq = self._reverse_complement(inv_seq)
                
                variants.append({
                    'type': 'INV',
                    'pos': pos,
                    'ref': inv_seq,
                    'alt': rev_comp_seq,
                    'hgvs': f"{hgvs_prefix}:g.{pos}_{pos+inv_size-1}inv"
                })
        
        return variants
    
    def _reverse_complement(self, seq: str) -> str:
        """Return the reverse complement of a DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join(complement.get(base, 'N') for base in reversed(seq))
    
    def export_vcf(self, variants: List[Dict], filename: str, contig_id="1", 
                 samples=None, phased=False) -> None:
        """
        Export variants to a VCF file
        
        Args:
            variants: List of variant dictionaries
            filename: Output VCF filename
            contig_id: Contig ID to use in VCF
            samples: List of sample names for multi-sample VCF (None = no samples)
            phased: Whether to generate phased genotypes (|) or unphased (/)
        """
        with open(filename, 'w') as f:
            # Write VCF header
            f.write("##fileformat=VCFv4.2\n")
            f.write("##reference=ref.fa\n")
            f.write(f"##contig=<ID={contig_id},length={self.seq_length}>\n")
            
            # Add structural variant INFO fields
            f.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
            f.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
            f.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">\n")
            f.write("##INFO=<ID=HGVS,Number=1,Type=String,Description=\"HGVS notation for the variant\">\n")
            
            # Add ALT fields for SVs
            f.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
            f.write("##ALT=<ID=INS,Description=\"Insertion\">\n")
            f.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
            f.write("##ALT=<ID=INV,Description=\"Inversion\">\n")
            
            # Add FORMAT fields if we have samples
            if samples:
                f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
                f.write("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles\">\n")
                f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n")
            
            # Write column headers
            header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
            if samples:
                header += "\tFORMAT"
                for sample in samples:
                    header += f"\t{sample}"
            f.write(header + "\n")
            
            # Store the phase separator
            phase_sep = "|" if phased else "/"
            
            # Write variants
            for i, var in enumerate(variants):
                pos = var['pos']
                ref = var['ref']
                alt = var['alt']
                var_id = f"VAR{i+1}"
                var_type = var['type']
                hgvs = var['hgvs']
                
                # Calculate SV length and end position
                if var_type == 'DEL':
                    svlen = -(len(ref) - 1)
                    end = pos + len(ref) - 1
                elif var_type == 'INS':
                    svlen = len(alt) - 1
                    end = pos
                elif var_type == 'DUP':
                    svlen = len(alt) - len(ref)
                    end = pos + len(ref) - 1
                elif var_type == 'INV':
                    svlen = 0
                    end = pos + len(ref) - 1
                else:  # SNP
                    svlen = 0
                    end = pos
                
                # Build INFO field
                info = f"SVTYPE={var_type};SVLEN={svlen};END={end};HGVS={hgvs}"
                
                # Start the variant line
                line = f"{contig_id}\t{pos}\t{var_id}\t{ref}\t{alt}\t100\tPASS\t{info}"
                
                # Add genotype information if samples are specified
                if samples:
                    line += "\tGT:AD:DP"
                    
                    # Generate genotype for each sample
                    for _ in samples:
                        # Randomly assign genotype (0/0, 0/1, or 1/1)
                        # With 20% homozygous ref, 50% heterozygous, 30% homozygous alt
                        r = random.random()
                        if r < 0.2:  # Homozygous reference
                            gt = f"0{phase_sep}0"
                            ad = f"30,0"
                            dp = "30"
                        elif r < 0.7:  # Heterozygous
                            gt = f"0{phase_sep}1"
                            ad = f"15,15"
                            dp = "30"
                        else:  # Homozygous alternate
                            gt = f"1{phase_sep}1"
                            ad = f"0,30"
                            dp = "30"
                        
                        line += f"\t{gt}:{ad}:{dp}"
                
                f.write(line + "\n")


def load_reference_sequence(ref_file: str) -> Tuple[str, str]:
    """
    Load reference sequence from FASTA file (compressed or uncompressed)
    
    Args:
        ref_file: Path to reference FASTA file
        
    Returns:
        Tuple of (sequence_id, sequence)
    """
    # Check if file is gzipped
    is_gzipped = ref_file.endswith('.gz')
    
    # Open appropriate file handle
    if is_gzipped:
        handle = gzip.open(ref_file, 'rt')  # 'rt' for text mode
    else:
        handle = open(ref_file, 'r')
    
    try:
        # Read first sequence in the file
        for record in SeqIO.parse(handle, 'fasta'):
            sequence_id = record.id
            sequence = str(record.seq)
            return sequence_id, sequence
            
    finally:
        handle.close()
    
    raise ValueError(f"Could not load reference sequence from {ref_file}")


def parse_var_types(var_types_str: str) -> List[str]:
    """
    Parse comma-separated variant types
    
    Args:
        var_types_str: Comma-separated string of variant types
        
    Returns:
        List of valid variant types
    """
    valid_types = {'SNP', 'DEL', 'INS', 'DUP', 'INV'}
    requested_types = [t.strip().upper() for t in var_types_str.split(',')]
    
    # Filter out invalid types
    valid_requested = [t for t in requested_types if t in valid_types]
    
    if not valid_requested:
        raise ValueError(f"No valid variant types specified. Valid types are: {', '.join(valid_types)}")
    
    return valid_requested


def parse_size_range(size_range_str: str) -> Dict[str, Tuple[int, int]]:
    """
    Parse size range specifications for SV types
    Format: "DEL:1-50,INS:1-20"
    
    Args:
        size_range_str: String with size range specifications
        
    Returns:
        Dictionary mapping variant types to size ranges
    """
    size_ranges = {}
    
    if not size_range_str:
        return size_ranges
    
    for item in size_range_str.split(','):
        if ':' not in item:
            continue
        
        var_type, range_str = item.split(':', 1)
        var_type = var_type.strip().upper()
        
        if '-' in range_str:
            try:
                min_size, max_size = map(int, range_str.split('-', 1))
                size_ranges[var_type] = (min_size, max_size)
            except ValueError:
                continue
    
    return size_ranges
