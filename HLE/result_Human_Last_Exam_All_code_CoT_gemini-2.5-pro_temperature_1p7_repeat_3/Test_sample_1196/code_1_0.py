import numpy as np
from scipy.stats import pearsonr

def simulate_gwas_ld():
    """
    Simulates a small dataset to demonstrate how Linkage Disequilibrium (LD)
    can create misleading associations in a Genome-Wide Association Study (GWAS).
    """
    print("### Simulation: LD Causing Misleading Association ###")
    print("We will simulate genotypes for 1000 individuals and 4 SNPs.")
    print("SNPs 1, 2, and 3 are tightly linked in an LD block.")
    print("SNP 4 is unlinked (on another 'chromosome').")
    print("The trait is directly caused by SNP 2 only.\n")

    # Parameters
    n_individuals = 1000
    
    # --- Generate Genotypes ---
    
    # For the LD block (SNPs 1, 2, 3), assume two haplotypes exist: '111' and '000'
    # This creates perfect LD (r^2 = 1.0)
    haplotype_choices = [0, 1] # 0 represents '000', 1 represents '111'
    # Each individual inherits two haplotypes
    h1 = np.random.choice(haplotype_choices, n_individuals)
    h2 = np.random.choice(haplotype_choices, n_individuals)
    
    # Genotypes are the sum of the two haplotype alleles (0, 1, or 2)
    snp1_genotype = h1 + h2
    snp2_genotype = h1 + h2 # Causal SNP
    snp3_genotype = h1 + h2

    # For the unlinked SNP 4, genotypes are generated independently
    # Allele frequency p=0.5
    snp4_allele1 = np.random.choice([0, 1], n_individuals)
    snp4_allele2 = np.random.choice([0, 1], n_individuals)
    snp4_genotype = snp4_allele1 + snp4_allele2

    # --- Simulate a Trait ---
    # Trait is determined by the genotype of the causal SNP (SNP2) plus some random noise
    noise = np.random.normal(0, 0.5, n_individuals)
    trait = snp2_genotype + noise

    print("--- Association Results ---")
    print("Calculating correlation between each SNP and the trait.")
    
    corr1, _ = pearsonr(snp1_genotype, trait)
    corr2, _ = pearsonr(snp2_genotype, trait)
    corr3, _ = pearsonr(snp3_genotype, trait)
    corr4, _ = pearsonr(snp4_genotype, trait)

    print(f"Correlation of SNP 1 (Linked, Non-Causal) with Trait: {corr1:.4f}")
    print(f"Correlation of SNP 2 (Linked, Causal)     with Trait: {corr2:.4f}")
    print(f"Correlation of SNP 3 (Linked, Non-Causal) with Trait: {corr3:.4f}")
    print(f"Correlation of SNP 4 (Unlinked, Non-Causal) with Trait: {corr4:.4f}\n")

    print("Observation: All three tightly linked SNPs show a strong association with the trait,")
    print("even though only SNP 2 is causal. The unlinked SNP shows no association.")
    print("This is 'misleading' because we cannot identify the true causal SNP from this result.\n")

    # --- LD Calculation Demonstration ---
    print("--- LD Calculation Example: r^2 between SNP1 and SNP2 ---")
    
    # In our simulation, the frequencies are known by design
    # Frequencies of allele '1' for SNP1 (pA) and SNP2 (pB)
    pA = 0.5
    pa = 1 - pA # freq of allele '0'
    pB = 0.5
    pb = 1 - pB # freq of allele '0'
    
    # Haplotype frequency of '11' (haplotype where SNP1=1 and SNP2=1)
    P_AB = 0.5
    
    # D = P(AB) - P(A)P(B)
    D = P_AB - (pA * pB)
    
    # r^2 = D^2 / (pA * pa * pB * pb)
    r_squared = (D**2) / (pA * pa * pB * pb)

    print("The formula for the LD statistic r^2 is: D^2 / (p(A) * p(a) * p(B) * p(b))")
    print(f"p(A) = {pA}, p(a) = {pa}")
    print(f"p(B) = {pB}, p(b) = {pb}")
    print(f"P(AB) = {P_AB}")
    print(f"D = P(AB) - p(A)p(B) = {P_AB} - ({pA} * {pB}) = {D}")
    print(f"r^2 = {D**2} / ({pA} * {pa} * {pB} * {pb}) = {r_squared:.2f}")

    print("\nThe r^2 value of 1.0 indicates perfect linkage, explaining the misleading association.")


if __name__ == '__main__':
    simulate_gwas_ld()
    print("\n<<<A>>>")
