"""
Zygosity determination for variant effect reporting.
"""

from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict

from variant_effect_report.core.models import Feature, Variant, FeatureEffect


class ZygosityAnalyzer:
    """Analyze homozygous, heterozygous, and compound heterozygous effects."""
    
    def __init__(self):
        """Initialize the zygosity analyzer."""
        pass
    
    def analyze_sample_haplotypes(self, feature_effects_by_haplotype: Dict[str, List[FeatureEffect]], 
                                 features: List[Feature],
                                 feature_by_id: Optional[Dict[str, Feature]] = None, 
                                 children_by_parent: Optional[Dict[str, List[Feature]]] = None) -> Dict:
        """
        Analyze differences between sample haplotypes to identify homo/heterozygous effects
        and compound heterozygous variants.
        
        Args:
            feature_effects_by_haplotype: Dictionary mapping haplotype names to lists of feature effects
            features: List of all features
            feature_by_id: Dictionary mapping feature IDs to features
            children_by_parent: Dictionary mapping parent feature IDs to lists of child features
            
        Returns:
            Dictionary with homozygous, heterozygous, and compound heterozygous effects
        """
        if len(feature_effects_by_haplotype) < 2:
            # Cannot determine zygosity with less than 2 haplotypes
            return {'homozygous': [], 'heterozygous': [], 'incomplete': True}
        
        # Get haplotype names
        haplotype_names = list(feature_effects_by_haplotype.keys())
        
        # Map feature IDs to effects for each haplotype
        feature_map = {}
        for hap_name, effects in feature_effects_by_haplotype.items():
            for effect in effects:
                feature_id = effect.feature.attributes.get('ID', 'unknown')
                if feature_id not in feature_map:
                    feature_map[feature_id] = {}
                
                feature_map[feature_id][hap_name] = effect
        
        # Identify homozygous vs heterozygous effects
        homozygous = []
        heterozygous = []
        
        for feature_id, hap_effects in feature_map.items():
            # Check if effects are identical across haplotypes
            effect_signatures = {}
            for hap_name, effect in hap_effects.items():
                # Create a signature of the effect for comparison
                # Include sequence and effect types to determine if truly identical
                sig = (
                    tuple(sorted(effect.effects)),
                    effect.alt_feature_seq,
                    tuple(v.get('id', '') for v in effect.variants)
                )
                effect_signatures[hap_name] = sig
            
            # Skip features without effects on all haplotypes
            if len(hap_effects) < len(haplotype_names):
                continue
            
            # Check if all signatures are the same
            is_homozygous = len(set(effect_signatures.values())) == 1
            
            if is_homozygous:
                # Use the effect from the first haplotype
                first_hap = haplotype_names[0]
                effect = hap_effects[first_hap]
                effect.zygosity = 'homozygous'
                homozygous.append(effect)
            else:
                # Create a heterozygous effect entry
                # We'll use the effect from the first haplotype but note the differences
                first_hap = haplotype_names[0]
                effect = hap_effects[first_hap]
                effect.zygosity = 'heterozygous'
                effect.haplotype_effects = hap_effects
                
                # Analyze differences
                diff_effects = set()
                for hap_name, hap_effect in hap_effects.items():
                    for e in hap_effect.effects:
                        diff_effects.add(e)
                
                effect.combined_effects = list(diff_effects)
                heterozygous.append(effect)
        
        # Identify compound heterozygous effects if we have feature hierarchy info
        compound_heterozygous = []
        if feature_by_id and children_by_parent:
            compound_heterozygous = self._identify_compound_heterozygous_effects(
                feature_effects_by_haplotype, 
                features, 
                feature_by_id, 
                children_by_parent
            )
        
        return {
            'homozygous': homozygous,
            'heterozygous': heterozygous,
            'gene_compound_heterozygous': compound_heterozygous,
            'incomplete': False
        }
    
    def _identify_compound_heterozygous_effects(self, feature_effects_by_haplotype: Dict[str, List[FeatureEffect]], 
                                              features: List[Feature], 
                                              feature_by_id: Dict[str, Feature], 
                                              children_by_parent: Dict[str, List[Feature]]) -> List[Dict]:
        """
        Identify compound heterozygous effects where different variants affect the same gene in different haplotypes.
        
        Compound heterozygosity occurs when a gene has different variants on each chromosome copy.
        This is important for recessive conditions where two different variants can result in disease.
        
        Args:
            feature_effects_by_haplotype: Dictionary mapping haplotype names to lists of feature effects
            features: List of all features
            feature_by_id: Dictionary mapping feature IDs to features
            children_by_parent: Dictionary mapping parent feature IDs to lists of child features
            
        Returns:
            List of compound heterozygous effect dictionaries
        """
        compound_het_effects = []
        
        # Need at least 2 haplotypes
        if len(feature_effects_by_haplotype) < 2:
            return compound_het_effects
        
        # Get list of haplotypes
        haplotype_names = list(feature_effects_by_haplotype.keys())
        
        # Map genes to their CDS features
        gene_to_cds = {}
        for feature in features:
            if feature.type == 'gene':
                gene_id = feature.attributes.get('ID')
                if not gene_id:
                    continue
                    
                # Find all CDS features associated with this gene
                cds_features = []
                # First find mRNAs for this gene
                mrnas = children_by_parent.get(gene_id, [])
                for mrna in mrnas:
                    mrna_id = mrna.attributes.get('ID')
                    if mrna_id:
                        # Find CDS features for this mRNA
                        cds_list = [f for f in children_by_parent.get(mrna_id, []) if f.type == 'CDS']
                        cds_features.extend(cds_list)
                
                if cds_features:
                    gene_to_cds[gene_id] = cds_features
        
        # Create a map of gene_id -> haplotype -> variants
        gene_variants_by_haplotype = {}
        
        # For each haplotype's effects
        for hap_name, effects in feature_effects_by_haplotype.items():
            # For each effect
            for effect in effects:
                feature = effect.feature
                feature_id = feature.attributes.get('ID', 'unknown')
                feature_type = feature.type
                
                # Skip if no variants or feature is not affected
                if not effect.variants or effect.effects == ['no_change']:
                    continue
                    
                # Only concerned with CDS features
                if feature_type != 'CDS':
                    continue
                    
                # Find the gene this CDS belongs to
                gene_id = None
                for gene, cds_list in gene_to_cds.items():
                    if any(cds.attributes.get('ID') == feature_id for cds in cds_list):
                        gene_id = gene
                        break
                
                if not gene_id:
                    continue
                    
                # Add variants to the map
                if gene_id not in gene_variants_by_haplotype:
                    gene_variants_by_haplotype[gene_id] = {}
                
                if hap_name not in gene_variants_by_haplotype[gene_id]:
                    gene_variants_by_haplotype[gene_id][hap_name] = []
                
                # Add only variants we haven't seen for this haplotype and gene
                for variant in effect.variants:
                    # Skip if we've already added this variant
                    if any(v.get('id') == variant.get('id') for v in gene_variants_by_haplotype[gene_id][hap_name]):
                        continue
                    gene_variants_by_haplotype[gene_id][hap_name].append(variant)
        
        # Identify compound heterozygous effects
        for gene_id, haplotype_variants in gene_variants_by_haplotype.items():
            # Need variants in at least 2 haplotypes
            if len(haplotype_variants) < 2:
                continue
                
            # Check that all haplotypes have at least one variant
            if any(len(variants) == 0 for variants in haplotype_variants.values()):
                continue
            
            # Get variant IDs for each haplotype
            variant_ids_by_haplotype = {}
            for hap, variants in haplotype_variants.items():
                variant_ids_by_haplotype[hap] = [v.get('id') for v in variants if v.get('id')]
            
            # Check if different variants are affecting different haplotypes
            is_compound_het = False
            
            # Check if the variants are different between haplotypes
            variant_sets = [set(ids) for ids in variant_ids_by_haplotype.values()]
            
            # If any haplotype has a unique variant not shared with other haplotypes
            # This is a potential compound heterozygous situation
            if len(variant_sets) >= 2:
                for i, hap_variants in enumerate(variant_sets):
                    other_variants = set().union(*[s for j, s in enumerate(variant_sets) if j != i])
                    unique_variants = hap_variants - other_variants
                    if unique_variants:
                        is_compound_het = True
                        break
            
            if is_compound_het:
                # Create compound heterozygous effect entry
                gene_feature = next((f for f in features if f.type == 'gene' and f.attributes.get('ID') == gene_id), None)
                if not gene_feature:
                    continue
                    
                # Get gene name
                gene_name = gene_feature.attributes.get('Name', gene_id)
                
                # Get all variants affecting this gene across haplotypes
                all_variants = []
                for hap, variants in haplotype_variants.items():
                    for var in variants:
                        if not any(v.get('id') == var.get('id') for v in all_variants):
                            all_variants.append(var)
                
                # Create effect entry with detailed information
                compound_effect = {
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                    'haplotype_variants': haplotype_variants,
                    'variants': all_variants,
                    'compound_heterozygous': True,
                    'feature': gene_feature,
                    # Add more detailed analysis of the combined effect
                    'details': {
                        'variants_by_haplotype': {hap: [v.get('id') for v in vars] 
                                                 for hap, vars in haplotype_variants.items()},
                    }
                }
                
                compound_het_effects.append(compound_effect)
        
        return compound_het_effects
