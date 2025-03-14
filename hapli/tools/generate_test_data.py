#!/usr/bin/env python3
"""
Generate test data for hapli.

This script generates synthetic test data for testing hapli functionality:
- GFA files with reference and alternate paths
- GFF3 files with genomic features
- VCF files with variants

Usage:
  python -m hapli.tools.generate_test_data [options]
"""

import argparse
import logging
import os
import random
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
from collections import defaultdict
import cyvcf2
import tempfile

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

def generate_random_sequence(length, gc_content=0.5):
    """Generate a random DNA sequence with specified GC content."""
    bases = []
    for _ in range(length):
        if random.random() < gc_content:
            bases.append(random.choice(['G', 'C']))
        else:
            bases.append(random.choice(['A', 'T']))
    return ''.join(bases)

def generate_gfa(output_file, seq_length=10000, num_segments=5, num_variants=10, 
                num_paths=2, variant_types=None):
    """Generate a synthetic GFA file with segments and reference path only."""
    if variant_types is None:
        variant_types = ['SNP', 'INS', 'DEL']
    
    logging.info(f"Generating GFA file with {num_segments} segments")
    
    # Generate segments
    segments = {}
    segment_length = seq_length // num_segments
    for i in range(1, num_segments + 1):
        seg_id = f"s{i}"
        # Vary segment length slightly
        length_variation = random.randint(-int(segment_length * 0.1), int(segment_length * 0.1))
        actual_length = max(100, segment_length + length_variation)
        segments[seg_id] = generate_random_sequence(actual_length)
    
    # Generate reference path
    ref_path = [(seg_id, '+') for seg_id in segments.keys()]
    
    # Build reference sequence and segment offsets for accurate variant positioning
    reference_seq = ""
    segment_offsets = {}
    current_offset = 0
    
    for seg_id, orientation in ref_path:
        segment_seq = segments[seg_id]
        segment_offsets[seg_id] = (current_offset, current_offset + len(segment_seq) - 1, orientation)
        reference_seq += segment_seq
        current_offset += len(segment_seq)
    
    # Generate variants for VCF (not used in GFA, but returned for VCF generation)
    variants = []
    # Track positions and regions to avoid overlaps
    used_regions = []  # List of (start, end) tuples
    
    for i in range(1, num_variants + 1):
        # Choose a random segment
        seg_id = random.choice(list(segments.keys()))
        segment_seq = segments[seg_id]
        
        # Choose a random position in the segment that doesn't overlap with existing variants
        max_attempts = 100
        attempts = 0
        valid_pos = False
        
        while not valid_pos and attempts < max_attempts:
            attempts += 1
            pos = random.randint(1, len(segment_seq) - 10)
            
            # Determine potential variant size (for checking overlaps)
            var_type = random.choice(variant_types)
            if var_type == 'DEL':
                # For deletions, reserve more space
                del_length = random.randint(1, 10)
                if pos + del_length > len(segment_seq):
                    continue  # Skip if deletion would go beyond segment end
                region_start = pos
                region_end = pos + del_length
            elif var_type == 'INS':
                # For insertions, just reserve the position
                region_start = pos
                region_end = pos + 1
            else:  # SNP
                region_start = pos
                region_end = pos + 1
            
            # Add buffer zone to prevent nearby variants
            buffer = 20  # Buffer zone around variants
            check_start = region_start - buffer
            check_end = region_end + buffer
            
            # Check if this region overlaps with any existing variant
            overlaps = False
            for start, end in used_regions:
                if (check_start <= end and check_end >= start):
                    overlaps = True
                    break
            
            if not overlaps:
                valid_pos = True
                used_regions.append((region_start, region_end))
        
        if not valid_pos:
            logging.warning(f"Could not find non-overlapping position for variant {i} after {max_attempts} attempts")
            continue
        
        # Create variant based on type, ensuring reference sequence matches
        if var_type == 'SNP':
            ref_base = segment_seq[pos-1]
            alt_bases = [b for b in 'ACGT' if b != ref_base]
            alt_base = random.choice(alt_bases)
            variants.append({
                'id': f"var{i}",
                'type': 'SNP',
                'segment': seg_id,
                'position': pos,
                'ref': ref_base,
                'alt': alt_base
            })
        elif var_type == 'INS':
            ref_base = segment_seq[pos-1]
            ins_length = random.randint(1, 10)
            ins_seq = generate_random_sequence(ins_length)
            variants.append({
                'id': f"var{i}",
                'type': 'INS',
                'segment': seg_id,
                'position': pos,
                'ref': ref_base,
                'alt': ref_base + ins_seq
            })
        elif var_type == 'DEL':
            del_length = random.randint(1, min(10, len(segment_seq) - pos))
            ref_seq = segment_seq[pos-1:pos+del_length]
            alt_base = segment_seq[pos-1]
            variants.append({
                'id': f"var{i}",
                'type': 'DEL',
                'segment': seg_id,
                'position': pos,
                'ref': ref_seq,
                'alt': alt_base
            })
    
    # Write GFA file
    with open(output_file, 'w') as f:
        # Write header
        f.write("H\tVN:Z:1.0\n")
        
        # Write segments
        for seg_id, seq in segments.items():
            f.write(f"S\t{seg_id}\t{seq}\tLN:i:{len(seq)}\n")
        
        # Write links
        for i in range(len(ref_path) - 1):
            from_seg, from_orient = ref_path[i]
            to_seg, to_orient = ref_path[i + 1]
            f.write(f"L\t{from_seg}\t{from_orient}\t{to_seg}\t{to_orient}\t0M\n")
        
        # Write reference path only
        ref_path_str = ','.join(f"{seg_id}{orient}" for seg_id, orient in ref_path)
        f.write(f"P\tREF\t{ref_path_str}\t*\n")
    
    logging.info(f"Generated GFA file: {output_file}")
    logging.info(f"  Segments: {len(segments)}")
    logging.info(f"  Variants: {len(variants)}")
    logging.info(f"  Paths: 1")  # Only reference path
    
    return segments, variants, ref_path, []  # Empty list for alt_paths

def generate_gff3(output_file, sequence_length=10000, num_genes=5, num_exons_per_gene=3):
    """Generate a synthetic GFF3 file with gene features."""
    logging.info(f"Generating GFF3 file with {num_genes} genes")
    
    # Calculate gene distribution
    gene_length = sequence_length // (num_genes * 2)  # Leave space between genes
    
    with open(output_file, 'w') as f:
        # Write GFF3 header
        f.write("##gff-version 3\n")
        f.write(f"##sequence-region 1 1 {sequence_length}\n")
        
        # Generate genes
        for i in range(1, num_genes + 1):
            # Calculate gene position
            gene_start = i * (sequence_length // (num_genes + 1))
            gene_start = max(1, gene_start - gene_length // 2)
            gene_end = min(sequence_length, gene_start + gene_length)
            
            # Randomly choose strand
            strand = '+' if random.random() > 0.3 else '-'
            
            # Write gene feature
            gene_id = f"gene{i}"
            f.write(f"1\thapli\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\tID={gene_id};Name=gene_{i}\n")
            
            # Write mRNA feature
            mrna_id = f"mRNA{i}"
            f.write(f"1\thapli\tmRNA\t{gene_start}\t{gene_end}\t.\t{strand}\t.\tID={mrna_id};Parent={gene_id}\n")
            
            # Generate exons
            exon_length = (gene_end - gene_start) // num_exons_per_gene
            for j in range(1, num_exons_per_gene + 1):
                exon_start = gene_start + (j - 1) * exon_length
                exon_end = exon_start + exon_length - 10  # Leave small gaps between exons
                
                # Write exon feature
                exon_id = f"exon{i}.{j}"
                f.write(f"1\thapli\texon\t{exon_start}\t{exon_end}\t.\t{strand}\t.\tID={exon_id};Parent={mrna_id}\n")
                
                # Write CDS feature (same as exon for simplicity)
                cds_id = f"cds{i}.{j}"
                f.write(f"1\thapli\tCDS\t{exon_start}\t{exon_end}\t.\t{strand}\t0\tID={cds_id};Parent={mrna_id}\n")
    
    logging.info(f"Generated GFF3 file: {output_file}")
    logging.info(f"  Genes: {num_genes}")
    logging.info(f"  Exons per gene: {num_exons_per_gene}")
    
    return num_genes, num_exons_per_gene

def generate_hgvs_notation(variant):
    """Generate HGVS notation for a variant that's compatible with the hgvs library."""
    var_type = variant['type']
    pos = variant['pos']
    
    # Use NC_000001.11 as a dummy accession for chromosome 1
    accession = "NC_000001.11"
    
    if var_type == 'SNP':
        # Format: NC_000001.11:g.100A>G
        return f"{accession}:g.{pos}{variant['ref']}>{variant['alt']}"
    elif var_type == 'INS':
        # Format: NC_000001.11:g.100_101insACGT
        return f"{accession}:g.{pos}_{pos+1}ins{variant['alt'][1:]}"
    elif var_type == 'DEL':
        # Format: NC_000001.11:g.100_105del
        end_pos = pos + len(variant['ref']) - 1
        return f"{accession}:g.{pos}_{end_pos}del"
    return ""

def generate_vcf(output_file, sequence_length=10000, num_variants=20, variant_types=None,
                num_samples=2, phased=True, reference_seq=None):
    """Generate a synthetic VCF file with variants using cyvcf2.
    
    Args:
        output_file: Path to output VCF file
        sequence_length: Length of the reference sequence
        num_variants: Number of variants to generate
        variant_types: List of variant types to generate
        num_samples: Number of samples in the VCF
        phased: Whether to generate phased genotypes
        reference_seq: Optional reference sequence to ensure variants match
    """
    if variant_types is None:
        variant_types = ['SNP', 'INS', 'DEL']
    
    logging.info(f"Generating VCF file with {num_variants} variants for {num_samples} samples")
    
    # Generate reference sequence if not provided
    if reference_seq is None:
        reference_seq = generate_random_sequence(sequence_length)
        logging.info(f"Generated random reference sequence of length {len(reference_seq)}")
    else:
        sequence_length = len(reference_seq)
        logging.info(f"Using provided reference sequence of length {sequence_length}")
    
    # Generate variants
    variants = []
    # Track positions and regions to avoid overlaps
    used_regions = []  # List of (start, end) tuples
    
    for i in range(1, num_variants + 1):
        # Choose a random position that doesn't overlap with existing variants
        max_attempts = 100
        attempts = 0
        valid_pos = False
        
        while not valid_pos and attempts < max_attempts:
            attempts += 1
            pos = random.randint(10, sequence_length - 20)
            
            # Determine potential variant size (for checking overlaps)
            var_type = random.choice(variant_types)
            if var_type == 'DEL':
                # For deletions, reserve more space
                del_length = random.randint(1, min(10, sequence_length - pos - 5))
                region_start = pos
                region_end = pos + del_length
            elif var_type == 'INS':
                # For insertions, just reserve the position
                region_start = pos
                region_end = pos + 1
            else:  # SNP
                region_start = pos
                region_end = pos + 1
            
            # Add buffer zone to prevent nearby variants
            buffer = 20  # Buffer zone around variants
            check_start = region_start - buffer
            check_end = region_end + buffer
            
            # Check if this region overlaps with any existing variant
            overlaps = False
            for start, end in used_regions:
                if (check_start <= end and check_end >= start):
                    overlaps = True
                    break
            
            if not overlaps:
                valid_pos = True
                used_regions.append((region_start, region_end))
        
        if not valid_pos:
            logging.warning(f"Could not find non-overlapping position for variant {i} after {max_attempts} attempts")
            continue
        
        # Create variant based on type, ensuring it matches the reference sequence
        if var_type == 'SNP':
            ref_base = reference_seq[pos-1]  # Use actual reference base
            alt_bases = [b for b in 'ACGT' if b != ref_base]
            alt_base = random.choice(alt_bases)
            variants.append({
                'id': f"var{i}",
                'type': 'SNP',
                'pos': pos,
                'ref': ref_base,
                'alt': alt_base
            })
        elif var_type == 'INS':
            ref_base = reference_seq[pos-1]  # Use actual reference base
            ins_length = random.randint(1, 10)
            ins_seq = generate_random_sequence(ins_length)
            variants.append({
                'id': f"var{i}",
                'type': 'INS',
                'pos': pos,
                'ref': ref_base,
                'alt': ref_base + ins_seq
            })
        elif var_type == 'DEL':
            # Ensure deletion doesn't go beyond sequence end
            del_length = min(random.randint(1, 10), sequence_length - pos)
            ref_seq = reference_seq[pos-1:pos+del_length]  # Use actual reference sequence
            alt_base = reference_seq[pos-1]  # First base of the deletion
            variants.append({
                'id': f"var{i}",
                'type': 'DEL',
                'pos': pos,
                'ref': ref_seq,
                'alt': alt_base
            })
    
    # Sort variants by position
    variants.sort(key=lambda v: v['pos'])
    
    # Add HGVS notation to each variant
    for variant in variants:
        variant['hgvs'] = generate_hgvs_notation(variant)
    
    # Create a temporary VCF file with cyvcf2
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_file:
        temp_filename = temp_file.name
        
        # Write VCF header
        temp_file.write("##fileformat=VCFv4.2\n")
        temp_file.write("##source=hapli_test_data_generator\n")
        temp_file.write("##reference=synthetic\n")
        temp_file.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of variant">\n')
        temp_file.write('##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS notation">\n')
        temp_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        
        # Write column headers
        sample_names = [f"SAMPLE{i}" for i in range(1, num_samples + 1)]
        temp_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(sample_names) + "\n")
        
        # Write variants
        for variant in variants:
            # Generate genotypes for each sample
            genotypes = []
            for _ in range(num_samples):
                # Randomly assign genotype
                if random.random() < 0.7:  # 70% chance of having the variant
                    if random.random() < 0.5:  # 50% chance of heterozygous
                        gt = "0|1" if phased else "0/1"
                    else:  # Homozygous alt
                        gt = "1|1" if phased else "1/1"
                else:  # Homozygous ref
                    gt = "0|0" if phased else "0/0"
                genotypes.append(gt)
            
            # Write VCF line
            temp_file.write(f"1\t{variant['pos']}\t{variant['id']}\t{variant['ref']}\t{variant['alt']}\t.\tPASS\t")
            genotypes_str = "\t".join(genotypes)
            temp_file.write(f"TYPE={variant['type']};HGVS={variant['hgvs']}\tGT\t{genotypes_str}\n")
    
    try:
        # Use cyvcf2 to read and write the VCF file to ensure it's compliant
        vcf_reader = cyvcf2.VCF(temp_filename)
        vcf_writer = cyvcf2.Writer(output_file, vcf_reader)
        
        for variant in vcf_reader:
            vcf_writer.write_record(variant)
        
        vcf_writer.close()
        vcf_reader.close()
        
    finally:
        # Clean up the temporary file
        if os.path.exists(temp_filename):
            os.unlink(temp_filename)
    
    logging.info(f"Generated VCF file: {output_file}")
    logging.info(f"  Variants: {len(variants)}")
    logging.info(f"  Samples: {num_samples}")
    
    return variants, sample_names

def generate_fasta(output_file, sequence_length=10000, num_sequences=1):
    """Generate a synthetic FASTA file with random sequences."""
    logging.info(f"Generating FASTA file with {num_sequences} sequences")
    
    sequences = []
    for i in range(1, num_sequences + 1):
        seq_id = f"seq{i}"
        seq = generate_random_sequence(sequence_length)
        
        record = SeqRecord(
            Seq(seq),
            id=seq_id,
            name=seq_id,
            description=f"Synthetic sequence {i} of length {sequence_length}"
        )
        sequences.append(record)
    
    # Write to file
    with open(output_file, 'w') as f:
        SeqIO.write(sequences, f, "fasta")
    
    logging.info(f"Generated FASTA file: {output_file}")
    logging.info(f"  Sequences: {num_sequences}")
    logging.info(f"  Length: {sequence_length}")
    
    return sequences

def validate_variants_against_reference(variants, reference_seq):
    """Validate that variants match the reference sequence.
    
    Args:
        variants: List of variant dictionaries
        reference_seq: Reference sequence string
        
    Returns:
        tuple: (is_valid, list of error messages)
    """
    errors = []
    is_valid = True
    
    for variant in variants:
        pos = variant['pos'] - 1  # Convert to 0-based
        ref = variant['ref']
        
        # Check if position is within reference bounds
        if pos >= len(reference_seq):
            errors.append(f"Variant {variant['id']} position {variant['pos']} is beyond reference length {len(reference_seq)}")
            is_valid = False
            continue
            
        # For deletions, check if the deletion extends beyond reference
        if variant['type'] == 'DEL' and pos + len(ref) > len(reference_seq):
            errors.append(f"Deletion variant {variant['id']} extends beyond reference end")
            is_valid = False
            continue
        
        # Extract actual reference sequence at variant position
        actual_ref = reference_seq[pos:pos+len(ref)]
        
        # Compare with variant's reference allele
        if actual_ref != ref:
            errors.append(f"Reference mismatch for variant {variant['id']} at position {variant['pos']}: "
                         f"Expected '{ref}', Found '{actual_ref}'")
            is_valid = False
    
    return is_valid, errors

def generate_test_dataset(output_dir, prefix="test", sequence_length=10000):
    """Generate a complete test dataset with GFA, GFF3, VCF, and FASTA files."""
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate FASTA file
    fasta_file = os.path.join(output_dir, f"{prefix}.fa")
    generate_fasta(fasta_file, sequence_length)
    
    # Generate GFA file
    gfa_file = os.path.join(output_dir, f"{prefix}.gfa")
    generate_gfa(gfa_file, sequence_length)
    
    # Generate GFF3 file
    gff_file = os.path.join(output_dir, f"{prefix}.gff3")
    generate_gff3(gff_file, sequence_length)
    
    # Generate VCF file
    vcf_file = os.path.join(output_dir, f"{prefix}.vcf")
    generate_vcf(vcf_file, sequence_length)
    
    logging.info(f"Generated complete test dataset in {output_dir}")
    logging.info(f"  FASTA: {fasta_file}")
    logging.info(f"  GFA: {gfa_file}")
    logging.info(f"  GFF3: {gff_file}")
    logging.info(f"  VCF: {vcf_file}")
    
    return {
        'fasta': fasta_file,
        'gfa': gfa_file,
        'gff3': gff_file,
        'vcf': vcf_file
    }

def main():
    """Main function to run the test data generator."""
    # Check if we're being called from hapli.py
    if len(sys.argv) > 1 and sys.argv[1] == 'generate-test-data':
        # Remove the subcommand from sys.argv
        sys.argv.pop(1)
    
    # Check if cyvcf2 is available
    try:
        import cyvcf2
    except ImportError:
        logging.error("cyvcf2 is required for VCF generation. Please install it with: pip install cyvcf2")
        return 1
    
    parser = argparse.ArgumentParser(description='Generate test data for hapli.')
    parser.add_argument('--output-dir', default='testdata', help='Output directory for test files')
    parser.add_argument('--prefix', default='test', help='Prefix for output files')
    parser.add_argument('--length', type=int, default=10000, help='Length of reference sequence')
    parser.add_argument('--test', action='store_true', help='Run tests to verify variant generation')
    
    # File-specific options
    parser.add_argument('--gfa', action='store_true', help='Generate GFA file')
    parser.add_argument('--gff3', action='store_true', help='Generate GFF3 file')
    parser.add_argument('--vcf', action='store_true', help='Generate VCF file')
    parser.add_argument('--fasta', action='store_true', help='Generate FASTA file')
    
    # GFA options
    parser.add_argument('--segments', type=int, default=5, help='Number of segments in GFA')
    parser.add_argument('--variants', type=int, default=10, help='Number of variants')
    parser.add_argument('--paths', type=int, default=2, help='Number of paths in GFA')
    
    # GFF3 options
    parser.add_argument('--genes', type=int, default=5, help='Number of genes in GFF3')
    parser.add_argument('--exons', type=int, default=3, help='Number of exons per gene')
    
    # VCF options
    parser.add_argument('--samples', type=int, default=2, help='Number of samples in VCF')
    parser.add_argument('--unphased', action='store_true', help='Generate unphased genotypes')
    
    # Debug and logging options
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose output without full debug')
    parser.add_argument('--log-file', help='Write log to this file')
    parser.add_argument('--seed', type=int, help='Random seed for reproducible data generation')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(debug=args.debug, log_file=args.log_file, verbose=args.verbose)
    
    # Set random seed if provided
    if args.seed:
        random.seed(args.seed)
        logging.info(f"Using random seed: {args.seed}")
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        # Run tests if requested
        if args.test:
            logging.info("Running tests...")
            test_result = test_variant_reference_matching()
            if test_result:
                logging.info("All tests passed!")
                return 0
            else:
                logging.error("Tests failed!")
                return 1
        
        # Generate all files by default if no specific file type is requested
        generate_all = not (args.gfa or args.gff3 or args.vcf or args.fasta)
        
        if generate_all:
            generate_test_dataset(args.output_dir, args.prefix, args.length)
        else:
            # Generate requested file types
            reference_seq = None
            
            # Generate FASTA first if requested
            if args.fasta:
                fasta_file = os.path.join(args.output_dir, f"{args.prefix}.fa")
                sequences = generate_fasta(fasta_file, args.length)
                if sequences:
                    reference_seq = str(sequences[0].seq)
            
            # Generate GFA next
            if args.gfa:
                gfa_file = os.path.join(args.output_dir, f"{args.prefix}.gfa")
                segments, variants, ref_path, _ = generate_gfa(gfa_file, args.length, args.segments, args.variants, args.paths)
                
                # Build reference sequence from GFA if not already created
                if reference_seq is None and segments:
                    reference_seq = ""
                    for seg_id, orientation in ref_path:
                        segment_seq = segments[seg_id]
                        if orientation == '-':
                            segment_seq = reverse_complement(segment_seq)
                        reference_seq += segment_seq
            
            # Generate GFF3
            if args.gff3:
                gff_file = os.path.join(args.output_dir, f"{args.prefix}.gff3")
                generate_gff3(gff_file, args.length, args.genes, args.exons)
            
            # Generate VCF last, using the reference sequence
            if args.vcf:
                vcf_file = os.path.join(args.output_dir, f"{args.prefix}.vcf")
                generate_vcf(vcf_file, args.length, args.variants, None, args.samples, 
                            not args.unphased, reference_seq=reference_seq)
        
        logging.info("Test data generation completed successfully")
        return 0
        
    except Exception as e:
        logging.error(f"Error during test data generation: {e}")
        if args.debug:
            import traceback
            logging.error(traceback.format_exc())
        return 1

def test_variant_reference_matching():
    """Test function to verify variants match reference sequence."""
    logging.info("Running variant reference matching test...")
    
    # Generate a reference sequence
    seq_length = 5000
    reference_seq = generate_random_sequence(seq_length)
    
    # Generate variants using this reference
    num_variants = 20
    variants = []
    used_regions = []
    
    for i in range(1, num_variants + 1):
        pos = random.randint(10, seq_length - 20)
        
        # Ensure no overlaps
        valid_pos = True
        for start, end in used_regions:
            if pos >= start - 20 and pos <= end + 20:
                valid_pos = False
                break
        
        if not valid_pos:
            continue
            
        var_type = random.choice(['SNP', 'INS', 'DEL'])
        
        if var_type == 'SNP':
            ref_base = reference_seq[pos-1]
            alt_bases = [b for b in 'ACGT' if b != ref_base]
            alt_base = random.choice(alt_bases)
            variants.append({
                'id': f"var{i}",
                'type': 'SNP',
                'pos': pos,
                'ref': ref_base,
                'alt': alt_base
            })
        elif var_type == 'INS':
            ref_base = reference_seq[pos-1]
            ins_seq = generate_random_sequence(5)
            variants.append({
                'id': f"var{i}",
                'type': 'INS',
                'pos': pos,
                'ref': ref_base,
                'alt': ref_base + ins_seq
            })
        elif var_type == 'DEL':
            del_length = min(5, seq_length - pos)
            ref_seq = reference_seq[pos-1:pos+del_length]
            alt_base = reference_seq[pos-1]
            variants.append({
                'id': f"var{i}",
                'type': 'DEL',
                'pos': pos,
                'ref': ref_seq,
                'alt': alt_base
            })
            
        used_regions.append((pos, pos + (len(variants[-1]['ref']) - 1)))
    
    # Validate variants
    is_valid, errors = validate_variants_against_reference(variants, reference_seq)
    
    if is_valid:
        logging.info("✅ Test passed: All variants match reference sequence")
    else:
        logging.error("❌ Test failed: Variants do not match reference sequence")
        for error in errors:
            logging.error(f"  - {error}")
    
    # Test with deliberately mismatched variants
    bad_variants = variants.copy()
    if bad_variants:
        # Modify a variant to create a mismatch
        bad_variants[0]['ref'] = 'X' + bad_variants[0]['ref'][1:] if bad_variants[0]['ref'] else 'X'
        
        is_valid, errors = validate_variants_against_reference(bad_variants, reference_seq)
        
        if not is_valid:
            logging.info("✅ Validation correctly detected mismatched variants")
        else:
            logging.error("❌ Validation failed to detect mismatched variants")
    
    return is_valid

if __name__ == "__main__":
    # Run tests if requested with --test flag
    if len(sys.argv) > 1 and sys.argv[1] == '--test':
        setup_logging(debug=True)
        test_variant_reference_matching()
        sys.exit(0)
    
    sys.exit(main())
