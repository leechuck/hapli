"""
VCF file generation utilities using cyvcf2.
"""

import os
import random
import logging
from typing import Dict, List, Tuple, Optional
import cyvcf2


class VCFGenerator:
    """Generate VCF files with structural variants"""
    
    def __init__(self, seq_length: int, contig_id: str = "1"):
        """
        Initialize the generator
        
        Args:
            seq_length: Length of the reference sequence
            contig_id: Contig ID to use in VCF
        """
        self.seq_length = seq_length
        self.contig_id = contig_id
    
    def export_vcf(self, variants: List[Dict], filename: str, 
                 samples: Optional[List[str]] = None, 
                 phased: bool = False) -> None:
        """
        Export variants to a VCF file using cyvcf2
        
        Args:
            variants: List of variant dictionaries
            filename: Output VCF filename
            samples: List of sample names for multi-sample VCF (None = no samples)
            phased: Whether to generate phased genotypes (|) or unphased (/)
            
        Raises:
            IOError: If there's an error writing the file
        """
        try:
            # Set up header lines
            header_lines = [
                "##fileformat=VCFv4.2",
                "##reference=ref.fa",
                f"##contig=<ID={self.contig_id},length={self.seq_length}>",
                "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
                "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
                "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">",
                "##INFO=<ID=HGVS,Number=1,Type=String,Description=\"HGVS notation for the variant\">",
                "##ALT=<ID=DEL,Description=\"Deletion\">",
                "##ALT=<ID=INS,Description=\"Insertion\">",
                "##ALT=<ID=DUP,Description=\"Duplication\">",
                "##ALT=<ID=INV,Description=\"Inversion\">"
            ]
            
            # Add FORMAT fields if we have samples
            if samples:
                header_lines.extend([
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles\">",
                    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">"
                ])
            
            # Create column headers
            columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
            if samples:
                columns.append("FORMAT")
                columns.extend(samples)
            
            header_lines.append("\t".join(columns))
            
            # Write header to file
            with open(filename, 'w') as f:
                for line in header_lines:
                    f.write(line + "\n")
            
            # Create a VCF writer
            # Note: cyvcf2 doesn't have a built-in writer, so we'll append to the file
            with open(filename, 'a') as f:
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
                    line = f"{self.contig_id}\t{pos}\t{var_id}\t{ref}\t{alt}\t100\tPASS\t{info}"
                    
                    # Add genotype information if samples are specified
                    if samples:
                        line += "\tGT:AD:DP"
                        
                        # Store the phase separator
                        phase_sep = "|" if phased else "/"
                        
                        # Generate genotype for each sample
                        for _ in samples:
                            # Randomly assign genotype (0/0, 0/1, or 1/1)
                            # With 20% homozygous ref, 50% heterozygous, 30% homozygous alt
                            r = random.random()
                            if r < 0.2:  # Homozygous reference
                                gt = f"0{phase_sep}0"
                                ad = "30,0"
                                dp = "30"
                            elif r < 0.7:  # Heterozygous
                                gt = f"0{phase_sep}1"
                                ad = "15,15"
                                dp = "30"
                            else:  # Homozygous alternate
                                gt = f"1{phase_sep}1"
                                ad = "0,30"
                                dp = "30"
                            
                            line += f"\t{gt}:{ad}:{dp}"
                    
                    f.write(line + "\n")
                    
        except Exception as e:
            logging.error(f"Error writing VCF file: {str(e)}")
            raise IOError(f"Failed to write VCF file: {str(e)}")
