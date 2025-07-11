import numpy as np
from scipy.stats import chi2_contingency

def simulate_gwas_ld():
    """
    Simulates a GWAS scenario to demonstrate how Tag SNPs in a haplotype block
    can lead to misleading associations due to Linkage Disequilibrium (LD).
    """
    # --- 1. Simulation Parameters ---
    n_individuals = 5000
    # Define a haplotype block with 5 SNPs. We'll have two common haplotypes.
    # Haplotype 1 (Risk):   T G A C G
    # Haplotype 2 (Normal): C A G T A
    haplotype_1 = np.array([1, 1, 1, 1, 1]) # Risk alleles coded as 1
    haplotype_2 = np.array([0, 0, 0, 0, 0]) # Normal alleles coded as 0
    
    # Let the 3rd SNP (index 2) be the true causal variant.
    causal_snp_idx_in_block = 2
    
    # We will genotype 3 "Tag SNPs" to represent this haplotype.
    tag_snp_indices_in_block = [0, 3, 4]

    # An unlinked SNP on another chromosome.
    unlinked_snp_maf = 0.3 # Minor allele frequency

    # --- 2. Generate Genotypes for the Population ---
    genotypes = np.zeros((n_individuals, 7), dtype=int) # 5 block SNPs + 1 unlinked SNP + 1 trait status
    
    # Assign haplotypes to each individual
    for i in range(n_individuals):
        # Each individual gets two haplotypes, one from each parent
        h1 = haplotype_1 if np.random.rand() < 0.5 else haplotype_2
        h2 = haplotype_1 if np.random.rand() < 0.5 else haplotype_2
        
        # Genotype is the sum of alleles from both haplotypes
        individual_block_genotype = h1 + h2
        genotypes[i, 0:5] = individual_block_genotype
        
        # Generate genotype for the unlinked SNP
        genotypes[i, 5] = np.random.binomial(2, unlinked_snp_maf)

    # --- 3. Generate Phenotype (Trait) based on the CAUSAL SNP ---
    trait_status = np.zeros(n_individuals, dtype=int)
    for i in range(n_individuals):
        causal_genotype = genotypes[i, causal_snp_idx_in_block]
        
        # Define penetrance: probability of trait given genotype at the causal SNP
        if causal_genotype == 2: # Two risk alleles
            prob_trait = 0.75
        elif causal_genotype == 1: # One risk allele
            prob_trait = 0.5
        else: # Zero risk alleles
            prob_trait = 0.25
            
        if np.random.rand() < prob_trait:
            trait_status[i] = 1 # Case (has the trait)
        else:
            trait_status[i] = 0 # Control (does not have the trait)

    # --- 4. Perform Association Tests (Chi-squared) ---
    def run_association_test(snp_genotypes, trait):
        """Calculates p-value from a chi-squared test."""
        # We test for an association between carrying at least one minor allele and the trait
        has_allele = snp_genotypes > 0
        
        # Create a 2x2 contingency table:
        #          [Has Trait, No Trait]
        # Has Allele  [  a,        b   ]
        # No Allele   [  c,        d   ]
        contingency_table = np.zeros((2, 2))
        contingency_table[0, 0] = np.sum((has_allele == 1) & (trait == 1))
        contingency_table[0, 1] = np.sum((has_allele == 1) & (trait == 0))
        contingency_table[1, 0] = np.sum((has_allele == 0) & (trait == 1))
        contingency_table[1, 1] = np.sum((has_allele == 0) & (trait == 0))
        
        chi2, p, dof, expected = chi2_contingency(contingency_table)
        return p

    # Test the true causal SNP
    causal_snp_p_value = run_association_test(genotypes[:, causal_snp_idx_in_block], trait_status)

    # Test the three Tag SNPs
    tag_snp_p_values = []
    for idx in tag_snp_indices_in_block:
        p_val = run_association_test(genotypes[:, idx], trait_status)
        tag_snp_p_values.append(p_val)
        
    # Test the unlinked SNP
    unlinked_snp_p_value = run_association_test(genotypes[:, 5], trait_status)
    
    # --- 5. Print Results ---
    print("GWAS Simulation Results:")
    print("-" * 40)
    print("This simulation shows how non-causal 'Tag SNPs' can show a strong association")
    print("with a trait because they are in Linkage Disequilibrium with the true 'Causal SNP'.")
    print("-" * 40)
    
    print(f"Association p-value for the TRUE CAUSAL SNP (SNP 3): {causal_snp_p_value:.4e}")
    print("\n--- Misleading Associations Due to LD ---")
    print("The following Tag SNPs are NOT causal, but are on the same inherited haplotype:")
    print(f"Association p-value for Tag SNP 1: {tag_snp_p_values[0]:.4e}")
    print(f"Association p-value for Tag SNP 4: {tag_snp_p_values[1]:.4e}")
    print(f"Association p-value for Tag SNP 5: {tag_snp_p_values[2]:.4e}")
    
    print("\n--- No Association for Unlinked SNP ---")
    print("This SNP is on a different 'chromosome' and shows no significant association:")
    print(f"Association p-value for Unlinked SNP: {unlinked_snp_p_value:.4f}")
    
if __name__ == '__main__':
    simulate_gwas_ld()
<<<E>>>