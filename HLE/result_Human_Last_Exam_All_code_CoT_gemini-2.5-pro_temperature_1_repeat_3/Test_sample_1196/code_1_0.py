import numpy as np
from scipy import stats

def demonstrate_ld_misleading_association():
    """
    Simulates a GWAS scenario to demonstrate how linkage disequilibrium (LD)
    can cause misleading associations.
    """
    # 1. Simulation Parameters
    n_individuals = 2000
    n_snps = 10
    # The true causal SNP is at index 4 (the 5th SNP)
    causal_snp_index = 4
    # The effect size of carrying one copy of the causal allele
    genetic_effect_beta = 0.4
    # Frequency of the haplotype containing the causal 'risk' allele
    risk_haplotype_freq = 0.25

    # 2. Simulate Haplotypes to create a perfect LD block
    # Haplotype 0 (major): all alleles are 0
    # Haplotype 1 (minor/risk): all alleles are 1
    # This creates a block where all SNPs are perfectly correlated (r^2 = 1).
    haplotype_0 = np.zeros(n_snps, dtype=int)
    haplotype_1 = np.ones(n_snps, dtype=int)

    # 3. Simulate Genotypes for the population
    # Each individual's genotype is the sum of two randomly chosen haplotypes.
    genotypes = np.zeros((n_individuals, n_snps), dtype=int)
    for i in range(n_individuals):
        # Paternal haplotype
        paternal_hap = haplotype_1 if np.random.rand() < risk_haplotype_freq else haplotype_0
        # Maternal haplotype
        maternal_hap = haplotype_1 if np.random.rand() < risk_haplotype_freq else haplotype_0
        # Genotype (allele count: 0, 1, or 2) is the sum
        genotypes[i, :] = paternal_hap + maternal_hap

    # 4. Simulate a complex trait (phenotype)
    # The trait is determined by the causal SNP's genotype plus random noise.
    causal_genotype = genotypes[:, causal_snp_index]
    random_noise = np.random.normal(loc=0, scale=1.0, size=n_individuals)
    # Phenotype = baseline + genetic effect + noise
    phenotype = 10 + (genetic_effect_beta * causal_genotype) + random_noise

    # 5. Perform a mock GWAS: test each SNP for association with the trait
    print("--- Simulating GWAS Results for a Region with High LD ---")
    print(f"The true causal variant is SNP #{causal_snp_index + 1}.")
    print("We will test all 10 tightly linked SNPs in the LD block.")
    print("A strong association for any non-causal SNP would be 'misleading'.\n")
    print("SNP | P-value      | Is Causal? | Note")
    print("----|--------------|-----------|---------------------------------")
    
    p_values = []
    for j in range(n_snps):
        snp_genotype = genotypes[:, j]
        # Perform a linear regression to get the p-value for association
        slope, intercept, r_value, p_value, std_err = stats.linregress(snp_genotype, phenotype)
        p_values.append(p_value)

        is_causal = "Yes" if j == causal_snp_index else "No"
        note = ""
        # Highlight the group of SNPs from option A's scenario
        if j in [0, 1, 2]: # Example: SNP #1, #2, #3
            note = "<- Part of a 'tightly linked' group"
        
        print(f" {j+1:<3}| {p_value:<12.2e} | {is_causal:<9} | {note}")

    print("---------------------------------------------------------")
    print("\nExplanation of Results:")
    print("Because all 10 SNPs are in a perfect LD block, they are inherited together.")
    print("As a result, a GWAS finds a highly significant association (a very small p-value)")
    print("for EVERY SNP in the block, not just the true causal one (SNP #5).")
    print("An investigator looking at the strong signals for SNPs #1, #2, and #3 might")
    print("incorrectly conclude one of them influences the trait. This is a classic example")
    print("of a misleading association due to linkage disequilibrium.")

if __name__ == '__main__':
    demonstrate_ld_misleading_association()
