"""
Text report generation for variant effects.
"""

import sys
import logging
from collections import Counter

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
