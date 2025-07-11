import numpy as np

def simulate_gwas_scenario():
    """
    Simulates genetic data to demonstrate how Linkage Disequilibrium (LD) can
    create misleading associations in a GWAS.
    """
    # Number of individuals in our study
    n_individuals = 1000

    # 1. Define the true CAUSAL variant and the Phenotype (trait)
    # The causal SNP is the "true" genetic cause of the trait.
    # Genotypes are 0 or 1.
    causal_snp = np.random.randint(0, 2, n_individuals)
    # The phenotype is determined by the causal SNP plus some random biological noise.
    noise = np.random.normal(0, 0.2, n_individuals)
    phenotype = causal_snp + noise

    # --- Scenario A: Three SNPs tightly linked within a single LD block ---
    # These SNPs are not the cause, but are highly correlated (in LD) with the causal SNP.
    # We simulate this by copying the causal SNP and flipping a few alleles (e.g., a 2% difference).
    def create_linked_snp(base_snp, divergence_rate=0.02):
        linked_snp = base_snp.copy()
        n_flips = int(len(linked_snp) * divergence_rate)
        flip_indices = np.random.choice(len(linked_snp), n_flips, replace=False)
        for idx in flip_indices:
            linked_snp[idx] = 1 - linked_snp[idx] # Flip the bit
        return linked_snp

    snp_A1 = create_linked_snp(causal_snp)
    snp_A2 = create_linked_snp(causal_snp)
    snp_A3 = create_linked_snp(causal_snp)


    # --- Scenario C: Two SNPs on different chromosomes ---
    # These SNPs are not in LD with the causal SNP. They assort independently.
    # We simulate this by generating them randomly.
    snp_C1 = np.random.randint(0, 2, n_individuals)
    snp_C2 = np.random.randint(0, 2, n_individuals)

    # 2. Perform a mock "Association Test"
    # In a real GWAS, this is a regression. Here, we'll use correlation as a proxy.
    # A high correlation means a strong (but potentially misleading) association.
    corr_causal = np.corrcoef(phenotype, causal_snp)[0, 1]
    corr_A1 = np.corrcoef(phenotype, snp_A1)[0, 1]
    corr_A2 = np.corrcoef(phenotype, snp_A2)[0, 1]
    corr_A3 = np.corrcoef(phenotype, snp_A3)[0, 1]
    corr_C1 = np.corrcoef(phenotype, snp_C1)[0, 1]
    corr_C2 = np.corrcoef(phenotype, snp_C2)[0, 1]

    # 3. Print the results
    print("--- Simulation of GWAS Association Signals ---")
    print(f"\nThis simulation shows how a SNP's correlation with a trait is measured.")
    print(f"The 'true' Causal SNP determines the trait.")
    print(f"Correlation of Phenotype with Causal SNP: {corr_causal:.4f}\n")

    print("--- Scenario A: Three SNPs in High LD with Causal SNP ---")
    print("These SNPs are NOT causal, but are physically close to the causal SNP on the chromosome.")
    print("Their signals are 'misleading' because they point strongly to the trait without being the cause.")
    print(f"Correlation with Phenotype for Linked SNP 1: {corr_A1:.4f}")
    print(f"Correlation with Phenotype for Linked SNP 2: {corr_A2:.4f}")
    print(f"Correlation with Phenotype for Linked SNP 3: {corr_A3:.4f}\n")

    print("--- Scenario C: Two SNPs on Different Chromosomes ---")
    print("These SNPs are unlinked to the Causal SNP.")
    print("As expected, they show little to no correlation with the trait.")
    print(f"Correlation with Phenotype for Unlinked SNP 1: {corr_C1:.4f}")
    print(f"Correlation with Phenotype for Unlinked SNP 2: {corr_C2:.4f}\n")

    print("Conclusion: The tightly linked SNPs in Scenario A show a strong, misleading association,")
    print("while the unlinked SNPs in Scenario C do not. This makes A the most likely answer.")


if __name__ == '__main__':
    simulate_gwas_scenario()
