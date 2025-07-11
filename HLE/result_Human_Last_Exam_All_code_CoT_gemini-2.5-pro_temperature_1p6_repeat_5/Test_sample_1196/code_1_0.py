import numpy as np
from scipy import stats

def demonstrate_ld_confounding():
    """
    Simulates a GWAS scenario to show how linkage disequilibrium (LD)
    can lead to misleading associations.
    """
    # Step 1: Define parameters for our simulation
    num_individuals = 2000
    num_snps_in_haplotype = 5
    # We will choose one SNP to be the causal one and another to be the "tag" SNP
    causal_snp_index = 3  # SNP #4 is causal
    tag_snp_index = 0     # SNP #1 is our tag SNP

    # Step 2: Create a population with two common haplotypes
    # Haplotype 0: [0, 0, 0, 0, 0] (all major alleles)
    # Haplotype 1: [1, 1, 1, 1, 1] (all minor alleles, perfectly linked in an LD block)
    # Each individual inherits two haplotypes, one from each parent.
    np.random.seed(42) # for reproducibility
    haplotype1_alleles = np.random.binomial(1, 0.3, num_individuals)
    haplotype2_alleles = np.random.binomial(1, 0.3, num_individuals)

    # Calculate genotypes. The genotype is the sum of alleles (0, 1, or 2).
    # Since all SNPs are perfectly linked, the genotype is the same for all SNPs.
    # genotypes will be a single vector of length num_individuals
    genotypes = haplotype1_alleles + haplotype2_alleles

    # We can represent the full genotype matrix, though all columns will be identical
    full_genotype_matrix = np.tile(genotypes, (num_snps_in_haplotype, 1)).T

    # Step 3: Simulate a quantitative trait
    # The trait is influenced ONLY by the causal SNP, plus some random noise.
    trait_baseline = 50.0
    causal_effect_size = 5.0
    noise_level = 4.0

    causal_snp_genotypes = full_genotype_matrix[:, causal_snp_index]
    random_noise = np.random.normal(0, noise_level, num_individuals)
    trait = trait_baseline + (causal_snp_genotypes * causal_effect_size) + random_noise

    # Step 4: Perform association tests (what a GWAS does)
    # We test the association between the trait and our measured "Tag SNP"
    tag_snp_genotypes = full_genotype_matrix[:, tag_snp_index]
    
    # Run linear regression: trait ~ tag_snp_genotype
    tag_snp_regr = stats.linregress(tag_snp_genotypes, trait)
    
    # For comparison, let's also test the true causal SNP
    causal_snp_regr = stats.linregress(causal_snp_genotypes, trait)

    # Step 5: Print and interpret the results
    print("--- Simulation of a Misleading GWAS Signal Due to a Tag SNP ---")
    print(f"We simulated a trait influenced ONLY by SNP #{causal_snp_index + 1}.")
    print(f"However, in our GWAS, we only measured the non-causal SNP #{tag_snp_index + 1}.")
    print("This Tag SNP is on the same inherited haplotype (in high LD) as the causal SNP.\n")

    print(f"Association results for the TRUE Causal SNP (SNP #{causal_snp_index + 1}):")
    print(f"  - R-squared: {causal_snp_regr.rvalue**2:.4f}")
    print(f"  - P-value:   {causal_snp_regr.pvalue:.2e}")
    print("  - Conclusion: A very strong, true association is detected, as expected.\n")
    
    print(f"Association results for the measured TAG SNP (SNP #{tag_snp_index + 1}):")
    print(f"  - R-squared: {tag_snp_regr.rvalue**2:.4f}")
    print(f"  - P-value:   {tag_snp_regr.pvalue:.2e}")
    print("  - Conclusion: A nearly identical strong association is found for the Tag SNP.\n")

    print("--- Final Interpretation ---")
    print("The GWAS result for the Tag SNP is highly significant, but this is a MISLEADING association.")
    print("It incorrectly suggests SNP #1 is functionally important.")
    print("The signal only exists because SNP #1 is a proxy for the true causal SNP #4, to which it is tightly linked.")
    print("This demonstrates the scenario in answer choice E, where tag SNPs for a haplotype can confound the identification of the true causal variant.")

demonstrate_ld_confounding()