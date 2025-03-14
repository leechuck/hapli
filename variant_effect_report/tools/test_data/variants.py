"""
Variant generation with biological consistency using hgvs library.
"""

import random
import logging
from typing import Dict, List, Tuple, Optional, Set
import hgvs.parser
import hgvs.validator
import hgvs.exceptions

from .sequence import SequenceHandler


class VariantGenerator:
    """Generate biologically consistent genomic variants"""
    
    def __init__(self, reference: str, seed: Optional[int] = None):
        """
        Initialize with reference sequence
        
        Args:
            reference: Reference DNA sequence
            seed: Random seed for reproducibility
        """
        if not reference:
            raise ValueError("Reference sequence cannot be empty")
            
        self.reference = reference
        self.seq_length = len(reference)
        self.seq_handler = SequenceHandler(seed)
        self.used_positions: Set[int] = set()  # Track positions to avoid overlaps
        
        # Initialize HGVS parser
        self.hp = hgvs.parser.Parser()
        
    def _is_position_available(self, start_pos: int, end_pos: int) -> bool:
        """
        Check if a position range is available (not overlapping with existing variants)
        
        Args:
            start_pos: Start position (1-based)
            end_pos: End position (inclusive)
            
        Returns:
            True if position range is available, False otherwise
        """
        # Check if any position in the range is already used
        for pos in range(start_pos, end_pos + 1):
            if pos in self.used_positions:
                return False
        return True
    
    def _reserve_positions(self, start_pos: int, end_pos: int) -> None:
        """
        Mark positions as used
        
        Args:
            start_pos: Start position (1-based)
            end_pos: End position (inclusive)
        """
        for pos in range(start_pos, end_pos + 1):
            self.used_positions.add(pos)
    
    def _create_hgvs(self, variant_type: str, position: int, ref: str, alt: str, 
                   contig_id: str = "NC_000001.11") -> str:
        """
        Create a valid HGVS expression
        
        Args:
            variant_type: Type of variant (SNP, DEL, INS, DUP, INV)
            position: Variant position (1-based)
            ref: Reference sequence
            alt: Alternate sequence
            contig_id: Contig ID for HGVS notation
            
        Returns:
            HGVS expression string
        """
        try:
            if variant_type == 'SNP':
                hgvs_str = f"{contig_id}:g.{position}{ref}>{alt}"
            elif variant_type == 'DEL':
                if len(ref) == 1:
                    hgvs_str = f"{contig_id}:g.{position}del"
                else:
                    hgvs_str = f"{contig_id}:g.{position}_{position+len(ref)-1}del"
            elif variant_type == 'INS':
                inserted_seq = alt[1:]  # Remove first base which is reference
                hgvs_str = f"{contig_id}:g.{position}_{position+1}ins{inserted_seq}"
            elif variant_type == 'DUP':
                if len(ref) == 1:
                    hgvs_str = f"{contig_id}:g.{position}dup"
                else:
                    hgvs_str = f"{contig_id}:g.{position}_{position+len(ref)-1}dup"
            elif variant_type == 'INV':
                hgvs_str = f"{contig_id}:g.{position}_{position+len(ref)-1}inv"
            else:
                raise ValueError(f"Unsupported variant type: {variant_type}")
            
            # Validate HGVS syntax
            var_c = self.hp.parse_hgvs_variant(hgvs_str)
            return hgvs_str
            
        except (hgvs.exceptions.HGVSError, ValueError) as e:
            logging.warning(f"Could not validate HGVS: {e}, using simplified format")
            return f"{contig_id}:g.{position}_{variant_type}"
    
    def generate_snp(self, hgvs_prefix: str = "NC_000001.11") -> Optional[Dict]:
        """
        Generate a single nucleotide polymorphism
        
        Args:
            hgvs_prefix: Prefix for HGVS notation
            
        Returns:
            Dictionary with variant information or None if no valid position found
        """
        # Try up to 10 times to find an available position
        for _ in range(10):
            pos = random.randint(1, self.seq_length)
            
            if not self._is_position_available(pos, pos):
                continue
                
            ref_base = self.reference[pos-1]
            alt_bases = [b for b in ['A', 'C', 'G', 'T'] if b != ref_base]
            
            if not alt_bases:  # Handle case of non-ACGT reference base
                continue
                
            alt_base = random.choice(alt_bases)
            
            # Reserve position
            self._reserve_positions(pos, pos)
            
            # Create HGVS notation
            hgvs = self._create_hgvs('SNP', pos, ref_base, alt_base, hgvs_prefix)
            
            return {
                'type': 'SNP',
                'pos': pos,
                'ref': ref_base,
                'alt': alt_base,
                'hgvs': hgvs
            }
            
        return None  # Could not find valid position
    
    def generate_deletion(self, min_size: int = 1, max_size: int = 50, 
                        hgvs_prefix: str = "NC_000001.11") -> Optional[Dict]:
        """
        Generate a deletion variant
        
        Args:
            min_size: Minimum deletion size
            max_size: Maximum deletion size
            hgvs_prefix: Prefix for HGVS notation
            
        Returns:
            Dictionary with variant information or None if no valid position found
        """
        if min_size < 1:
            raise ValueError("Minimum deletion size must be at least 1")
            
        # Try up to 20 times to find an available position
        for _ in range(20):
            del_size = random.randint(min_size, max_size)
            
            # Ensure we have enough sequence left
            if del_size >= self.seq_length - 1:
                del_size = self.seq_length - 2  # Leave at least one base
                
            if del_size < 1:
                return None  # Reference too short
                
            # Choose position
            max_start = self.seq_length - del_size
            if max_start < 1:
                return None  # Reference too short
                
            pos = random.randint(1, max_start)
            
            # Check if position range is available
            if not self._is_position_available(pos, pos + del_size - 1):
                continue
                
            # Get deletion sequence
            del_seq = self.reference[pos-1:pos+del_size-1]
            
            # Reserve positions
            self._reserve_positions(pos, pos + del_size - 1)
            
            # Create HGVS notation
            hgvs = self._create_hgvs('DEL', pos, del_seq, del_seq[0], hgvs_prefix)
            
            return {
                'type': 'DEL',
                'pos': pos,
                'ref': del_seq,
                'alt': del_seq[0],  # Keep first base
                'hgvs': hgvs
            }
            
        return None  # Could not find valid position
    
    def generate_insertion(self, min_size: int = 1, max_size: int = 20, 
                         hgvs_prefix: str = "NC_000001.11") -> Optional[Dict]:
        """
        Generate an insertion variant
        
        Args:
            min_size: Minimum insertion size
            max_size: Maximum insertion size
            hgvs_prefix: Prefix for HGVS notation
            
        Returns:
            Dictionary with variant information or None if no valid position found
        """
        if min_size < 1:
            raise ValueError("Minimum insertion size must be at least 1")
            
        # Try up to 10 times to find an available position
        for _ in range(10):
            ins_size = random.randint(min_size, max_size)
            
            # Choose position (must be at least 1 and less than sequence length)
            if self.seq_length < 2:
                return None  # Reference too short
                
            pos = random.randint(1, self.seq_length - 1)
            
            # Check if position is available
            if not self._is_position_available(pos, pos):
                continue
                
            # Generate insertion sequence
            ins_seq = self.seq_handler.generate_random_sequence(ins_size)
            ref_base = self.reference[pos-1]
            
            # Reserve position
            self._reserve_positions(pos, pos)
            
            # Create HGVS notation
            hgvs = self._create_hgvs('INS', pos, ref_base, ref_base + ins_seq, hgvs_prefix)
            
            return {
                'type': 'INS',
                'pos': pos,
                'ref': ref_base,
                'alt': ref_base + ins_seq,
                'hgvs': hgvs
            }
            
        return None  # Could not find valid position
    
    def generate_duplication(self, min_size: int = 1, max_size: int = 50, 
                           hgvs_prefix: str = "NC_000001.11") -> Optional[Dict]:
        """
        Generate a duplication variant
        
        Args:
            min_size: Minimum duplication size
            max_size: Maximum duplication size
            hgvs_prefix: Prefix for HGVS notation
            
        Returns:
            Dictionary with variant information or None if no valid position found
        """
        if min_size < 1:
            raise ValueError("Minimum duplication size must be at least 1")
            
        # Try up to 20 times to find an available position
        for _ in range(20):
            dup_size = random.randint(min_size, max_size)
            
            # Ensure we have enough sequence left
            if dup_size >= self.seq_length:
                dup_size = self.seq_length - 1  # Use almost all sequence
                
            if dup_size < 1:
                return None  # Reference too short
                
            # Choose position
            max_start = self.seq_length - dup_size
            if max_start < 1:
                return None  # Reference too short
                
            pos = random.randint(1, max_start)
            
            # Check if position range is available
            if not self._is_position_available(pos, pos + dup_size - 1):
                continue
                
            # Get duplication sequence
            dup_seq = self.reference[pos-1:pos+dup_size-1]
            
            # Reserve positions
            self._reserve_positions(pos, pos + dup_size - 1)
            
            # Create HGVS notation
            hgvs = self._create_hgvs('DUP', pos, dup_seq, dup_seq + dup_seq, hgvs_prefix)
            
            return {
                'type': 'DUP',
                'pos': pos,
                'ref': dup_seq,
                'alt': dup_seq + dup_seq,
                'hgvs': hgvs
            }
            
        return None  # Could not find valid position
    
    def generate_inversion(self, min_size: int = 3, max_size: int = 50, 
                         hgvs_prefix: str = "NC_000001.11") -> Optional[Dict]:
        """
        Generate an inversion variant
        
        Args:
            min_size: Minimum inversion size
            max_size: Maximum inversion size
            hgvs_prefix: Prefix for HGVS notation
            
        Returns:
            Dictionary with variant information or None if no valid position found
        """
        if min_size < 2:
            raise ValueError("Minimum inversion size must be at least 2")
            
        # Try up to 20 times to find an available position
        for _ in range(20):
            inv_size = random.randint(min_size, max_size)
            
            # Ensure we have enough sequence left
            if inv_size >= self.seq_length:
                inv_size = self.seq_length - 1  # Use almost all sequence
                
            if inv_size < 2:
                return None  # Reference too short
                
            # Choose position
            max_start = self.seq_length - inv_size
            if max_start < 1:
                return None  # Reference too short
                
            pos = random.randint(1, max_start)
            
            # Check if position range is available
            if not self._is_position_available(pos, pos + inv_size - 1):
                continue
                
            # Get inversion sequence
            inv_seq = self.reference[pos-1:pos+inv_size-1]
            rev_comp_seq = self.seq_handler.reverse_complement(inv_seq)
            
            # Reserve positions
            self._reserve_positions(pos, pos + inv_size - 1)
            
            # Create HGVS notation
            hgvs = self._create_hgvs('INV', pos, inv_seq, rev_comp_seq, hgvs_prefix)
            
            return {
                'type': 'INV',
                'pos': pos,
                'ref': inv_seq,
                'alt': rev_comp_seq,
                'hgvs': hgvs
            }
            
        return None  # Could not find valid position
    
    def generate_variants(self, num_variants: int, var_types: List[str], 
                        size_ranges: Optional[Dict[str, Tuple[int, int]]] = None,
                        hgvs_prefix: str = "NC_000001.11") -> List[Dict]:
        """
        Generate a list of non-overlapping variants
        
        Args:
            num_variants: Number of variants to generate
            var_types: List of variant types to include (SNP, DEL, INS, DUP, INV)
            size_ranges: Dictionary with min/max size for each variant type
            hgvs_prefix: Prefix for HGVS notation
            
        Returns:
            List of variant dictionaries
        """
        if not var_types:
            raise ValueError("At least one variant type must be specified")
            
        # Default size ranges
        if size_ranges is None:
            size_ranges = {
                'DEL': (1, 50),
                'INS': (1, 20),
                'DUP': (1, 50),
                'INV': (3, 50)
            }
            
        variants = []
        attempts = 0
        max_attempts = num_variants * 10  # Limit total attempts
        
        while len(variants) < num_variants and attempts < max_attempts:
            attempts += 1
            
            # Choose variant type
            var_type = random.choice(var_types)
            
            # Generate variant based on type
            variant = None
            
            if var_type == 'SNP':
                variant = self.generate_snp(hgvs_prefix)
                
            elif var_type == 'DEL':
                min_size, max_size = size_ranges.get('DEL', (1, 50))
                variant = self.generate_deletion(min_size, max_size, hgvs_prefix)
                
            elif var_type == 'INS':
                min_size, max_size = size_ranges.get('INS', (1, 20))
                variant = self.generate_insertion(min_size, max_size, hgvs_prefix)
                
            elif var_type == 'DUP':
                min_size, max_size = size_ranges.get('DUP', (1, 50))
                variant = self.generate_duplication(min_size, max_size, hgvs_prefix)
                
            elif var_type == 'INV':
                min_size, max_size = size_ranges.get('INV', (3, 50))
                variant = self.generate_inversion(min_size, max_size, hgvs_prefix)
            
            # Add variant if successfully generated
            if variant:
                variants.append(variant)
        
        return variants


def parse_var_types(var_types_str: str) -> List[str]:
    """
    Parse comma-separated variant types
    
    Args:
        var_types_str: Comma-separated string of variant types
        
    Returns:
        List of valid variant types
        
    Raises:
        ValueError: If no valid variant types are specified
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
        size_range_str: Size range specification string
        
    Returns:
        Dictionary mapping variant types to (min_size, max_size) tuples
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
                if min_size > max_size:
                    logging.warning(f"Invalid size range for {var_type}: min > max. Swapping values.")
                    min_size, max_size = max_size, min_size
                size_ranges[var_type] = (min_size, max_size)
            except ValueError:
                logging.warning(f"Could not parse size range for {var_type}: {range_str}")
                continue
    
    return size_ranges
