import numpy as np
from scipy.stats import chi2_contingency

def simulate_gwas_ld():
    """
    Simulates a GWAS scenario to demonstrate how tightly linked SNPs
    can lead to misleading associations due to Linkage Disequilibrium (LD).
    """
    # --- 1. Simulation Parameters ---
    n_individuals = 2000  # Number of individuals in the study
    haplotype1_freq = 0.2  # Frequency of the risk haplotype
    
    # Define two haplotypes in an LD block.
    # Haplotype 1 (risk): Causal=1, Linked1=1, Linked2=1
    # Haplotype 2 (non-risk): Causal=0, Linked1=0, Linked2=0
    # This creates perfect LD (r^2 = 1) between the three SNPs.
    haplotypes = {
        1: np.array([1, 1, 1]), # Risk haplotype
        0: np.array([0, 0, 0])  # Non-risk haplotype
    }
    
    # --- 2. Generate Genotypes ---
    # Generate genotypes for the LD block
    hap1_draws = np.random.choice([1, 0], size=(n_individuals, 2), p=[haplotype1_freq, 1 - haplotype1_freq])
    
    # Each individual has two haplotypes (one from each parent)
    genotypes_ld_block = np.zeros((n_individuals, 3), dtype=int)
    for i in range(n_individuals):
        h1 = haplotypes[hap1_draws[i, 0]]
        h2 = haplotypes[hap1_draws[i, 1]]
        genotypes_ld_block[i, :] = h1 + h2 # Additive model (0, 1, or 2)

    # Generate genotypes for an unlinked SNP on another chromosome
    unlinked_snp_freq = 0.3
    genotype_unlinked = np.random.binomial(2, unlinked_snp_freq, n_individuals)
    
    # SNP names for clarity
    snp_names = ["Causal SNP", "Linked SNP 1", "Linked SNP 2", "Unlinked SNP"]
    all_genotypes = np.column_stack([genotypes_ld_block, genotype_unlinked])

    # --- 3. Generate Phenotypes (Trait Status) ---
    # The probability of having the trait depends ONLY on the causal SNP's genotype.
    phenotypes = np.zeros(n_individuals, dtype=int) # 0=control, 1=case
    base_risk = 0.10 # Base risk for genotype 0
    risk_per_allele = 0.15 # Additional risk per risk allele
    
    causal_genotypes = all_genotypes[:, 0]
    for i in range(n_individuals):
        risk = base_risk + causal_genotypes[i] * risk_per_allele
        if np.random.rand() < risk:
            phenotypes[i] = 1 # Assign as "case"

    # --- 4. Perform Association Test for each SNP ---
    print("Running Association Tests (Chi-Squared Test)...")
    print("-" * 50)
    
    for i in range(all_genotypes.shape[1]):
        snp_genotype = all_genotypes[:, i]
        
        # Create a 2x3 contingency table: (case/control) vs (genotype 0/1/2)
        contingency_table = np.zeros((2, 3), dtype=int)
        for genotype_val in [0, 1, 2]:
            is_genotype = (snp_genotype == genotype_val)
            contingency_table[0, genotype_val] = np.sum(is_genotype & (phenotypes == 0)) # Controls
            contingency_table[1, genotype_val] = np.sum(is_genotype & (phenotypes == 1)) # Cases

        # Perform Chi-Squared test
        chi2, p, dof, expected = chi2_contingency(contingency_table)

        print(f"SNP Name:       {snp_names[i]}")
        print(f"P-value:        {p:.4e}")
        if p < 0.001:
            print("Result:         Strongly associated with the trait.")
        else:
            print("Result:         Not associated with the trait.")
        print("-" * 50)
        
    print("\nConclusion:")
    print("As shown above, the 'Linked SNP 1' and 'Linked SNP 2' show a very strong (and misleading) association with the trait.")
    print("This is not because they are causal, but because they are tightly linked to the true 'Causal SNP'.")
    print("The 'Unlinked SNP' shows no association, as expected.")


if __name__ == '__main__':
    simulate_gwas_ld()