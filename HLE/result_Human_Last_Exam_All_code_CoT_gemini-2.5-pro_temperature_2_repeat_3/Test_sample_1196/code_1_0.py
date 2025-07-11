import numpy as np
from scipy.stats import pearsonr

def simulate_gwas_ld():
    """
    Simulates a scenario illustrating how Linkage Disequilibrium (LD)
    can cause misleading associations in a GWAS.
    """
    # Let's define two common haplotypes (blocks of linked SNPs) in a population.
    # Haplotype 1 has all 0s, Haplotype 2 has all 1s.
    # This represents a perfect LD block where all alleles are linked.
    # snp_positions =   [0, 1, 2, 3, 4]
    haplotype_A =       [0, 0, 0, 0, 0]
    haplotype_B =       [1, 1, 1, 1, 1]

    # Let's assume the SNP at index 2 is the TRUE CAUSAL variant for a trait.
    # The other SNPs have no biological effect. We'll choose SNP at index 0 as our "Tag SNP".
    true_causal_snp_idx = 2
    tag_snp_idx = 0

    print(f"Scenario Setup:")
    print(f" - A perfect haplotype block is simulated (like in choice E).")
    print(f" - The True Causal SNP is at index: {true_causal_snp_idx}")
    print(f" - We will genotype a 'Tag SNP' at index: {tag_snp_idx}\n")

    # Now, we create a population of 1000 individuals.
    # Each individual gets one of the two haplotypes randomly.
    num_individuals = 1000
    population_haplotypes = np.array([haplotype_A if np.random.rand() > 0.5 else haplotype_B for _ in range(num_individuals)])

    # Get the genotypes for our two SNPs of interest across the population
    true_causal_genotypes = population_haplotypes[:, true_causal_snp_idx]
    tag_genotypes = population_haplotypes[:, tag_snp_idx]

    # The trait is determined by the TRUE CAUSAL SNP, plus some random noise.
    # If the causal allele is 1, the trait value is higher.
    trait_values = 10 * true_causal_genotypes + np.random.normal(0, 1, num_individuals)

    # --- GWAS Analysis ---
    # Now, let's test the association between each SNP and the trait.
    # We use Pearson correlation as a simple association test.

    # 1. Association of the True Causal SNP with the trait
    corr_causal, p_val_causal = pearsonr(true_causal_genotypes, trait_values)
    print("Association with the True Causal SNP (the real cause):")
    print(f"  Correlation: {corr_causal:.4f}")
    print(f"  P-value: {p_val_causal:.2e}")
    print("  Result: As expected, a very strong and significant association!\n")

    # 2. Association of the Tag SNP with the trait
    corr_tag, p_val_tag = pearsonr(tag_genotypes, trait_values)
    print("Association with the Tag SNP (the proxy):")
    print(f"  Correlation: {corr_tag:.4f}")
    print(f"  P-value: {p_val_tag:.2e}")
    print("  Result: A very strong and significant association, even though this SNP is not causal.\n")
    
    print("Conclusion:")
    print("The Tag SNP shows a powerful association with the trait, almost identical to the True Causal SNP's association.")
    print("This is a 'misleading association' because the Tag SNP is just a proxy due to being on the same inherited haplotype (in perfect LD).")
    print("This perfectly illustrates option E, where Tag SNPs that predict an entire haplotype lead to such findings.")

simulate_gwas_ld()
<<<E>>>