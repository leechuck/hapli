"""
Sequence generation and manipulation using BioPython.
"""

import random
import logging
import gzip
from typing import Dict, List, Tuple, Optional
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils import gc_fraction


class SequenceHandler:
    """Handle DNA sequence operations using BioPython"""
    
    def __init__(self, seed: Optional[int] = None):
        """Initialize with optional random seed"""
        if seed is not None:
            random.seed(seed)
    
    def generate_random_sequence(self, length: int) -> str:
        """
        Generate a random DNA sequence of specified length
        
        Args:
            length: Length of sequence to generate
            
        Returns:
            Random DNA sequence
        """
        if length <= 0:
            raise ValueError("Sequence length must be positive")
            
        bases = ['A', 'C', 'G', 'T']
        return ''.join(random.choice(bases) for _ in range(length))
    
    @staticmethod
    def reverse_complement(seq: str) -> str:
        """
        Return the reverse complement of a DNA sequence using BioPython
        
        Args:
            seq: DNA sequence
            
        Returns:
            Reverse complement of the sequence
        """
        try:
            return str(Seq(seq).reverse_complement())
        except Exception as e:
            logging.error(f"Error generating reverse complement: {e}")
            raise
    
    @staticmethod
    def load_fasta(filename: str) -> Tuple[str, str]:
        """
        Load a sequence from a FASTA file (compressed or uncompressed)
        
        Args:
            filename: Path to FASTA file (.fa or .fa.gz)
            
        Returns:
            Tuple of (sequence_id, sequence)
            
        Raises:
            FileNotFoundError: If the file doesn't exist
            ValueError: If the file is empty or not in FASTA format
        """
        try:
            # Check if file is gzipped
            is_gzipped = filename.endswith('.gz')
            
            # Open appropriate file handle
            if is_gzipped:
                handle = gzip.open(filename, 'rt')  # 'rt' for text mode
            else:
                handle = open(filename, 'r')
            
            # Read first sequence in the file
            for record in SeqIO.parse(handle, 'fasta'):
                sequence_id = record.id
                sequence = str(record.seq)
                handle.close()
                
                if not sequence:
                    raise ValueError(f"Empty sequence in {filename}")
                    
                return sequence_id, sequence
                
            # If we get here, no sequences were found
            handle.close()
            raise ValueError(f"No sequences found in {filename}")
            
        except FileNotFoundError:
            logging.error(f"File not found: {filename}")
            raise
        except Exception as e:
            logging.error(f"Error loading FASTA file: {e}")
            raise
    
    @staticmethod
    def save_fasta(sequence: str, sequence_id: str, filename: str) -> None:
        """
        Save a sequence to a FASTA file
        
        Args:
            sequence: DNA sequence
            sequence_id: Sequence identifier
            filename: Output filename
            
        Raises:
            IOError: If there's an error writing the file
        """
        try:
            record = SeqRecord(Seq(sequence), id=sequence_id, description="")
            SeqIO.write(record, filename, "fasta")
        except Exception as e:
            logging.error(f"Error saving FASTA file: {e}")
            raise IOError(f"Failed to write FASTA file: {str(e)}")
    
    @staticmethod
    def get_sequence_stats(sequence: str) -> Dict:
        """
        Get statistics for a DNA sequence
        
        Args:
            sequence: DNA sequence
            
        Returns:
            Dictionary with sequence statistics
        """
        return {
            'length': len(sequence),
            'gc_content': gc_fraction(sequence) * 100,  # Convert to percentage
            'molecular_weight': molecular_weight(sequence),
            'base_counts': {
                'A': sequence.upper().count('A'),
                'C': sequence.upper().count('C'),
                'G': sequence.upper().count('G'),
                'T': sequence.upper().count('T'),
                'N': sequence.upper().count('N')
            }
        }


def load_reference_sequence(ref_file: str) -> Tuple[str, str]:
    """
    Load reference sequence from FASTA file (wrapper for backward compatibility)
    
    Args:
        ref_file: Path to reference FASTA file
        
    Returns:
        Tuple of (sequence_id, sequence)
    """
    handler = SequenceHandler()
    return handler.load_fasta(ref_file)
