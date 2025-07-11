import numpy as np
from scipy.stats import pearsonr

def simulate_gwas_ld():
    """
    Simulates a scenario to demonstrate how Linkage Disequilibrium (LD)
    can cause misleading associations in a GWAS.
    """
    # Number of individuals in the study
    n_individuals = 2000

    # --- Step 1: Simulate Genotypes for four SNPs ---

    # SNP1: The true causal SNP (e.g., allele '1' increases trait value)
    # Genotypes are 0, 1, or 2 (representing copies of the minor allele)
    causal_snp = np.random.choice([0, 1, 2], size=n_individuals, p=[0.49, 0.42, 0.09])

    # SNP2 & SNP3: Two SNPs in high Linkage Disequilibrium (LD) with the causal SNP
    # We simulate "tight linkage" by making them identical to the causal SNP most of the time.
    # 99% of the time, their genotype is the same as the causal SNP's.
    ld_level = 0.99
    linked_snp_1 = np.where(np.random.random(n_individuals) < ld_level, causal_snp, np.random.choice([0, 1, 2], n_individuals))
    linked_snp_2 = np.where(np.random.random(n_individuals) < ld_level, causal_snp, np.random.choice([0, 1, 2], n_individuals))

    # SNP4: A SNP not in LD with the causal SNP (on another chromosome)
    # Its genotypes are generated independently.
    unlinked_snp = np.random.choice([0, 1, 2], size=n_individuals, p=[0.25, 0.5, 0.25])

    # --- Step 2: Simulate a Phenotype (the Trait) ---
    # The trait value is directly influenced ONLY by the causal_snp, plus random noise.
    beta_causal = 0.5
    noise = np.random.normal(0, 1, n_individuals)
    trait = beta_causal * causal_snp + noise

    # --- Step 3: Perform Association Tests (using correlation as a proxy) ---
    # In a real GWAS, we would use linear regression and get a p-value.
    # For this illustration, Pearson correlation coefficient is sufficient.
    corr_causal, _ = pearsonr(causal_snp, trait)
    corr_linked_1, _ = pearsonr(linked_snp_1, trait)
    corr_linked_2, _ = pearsonr(linked_snp_2, trait)
    corr_unlinked, _ = pearsonr(unlinked_snp, trait)

    # --- Step 4: Print and Explain the Results ---
    print("GWAS Simulation: Linkage Disequilibrium (LD) Demonstration\n")
    print(f"The trait is DIRECTLY caused by 'Causal SNP' only.\n")
    print("Association between each SNP and the Trait (measured by correlation):")
    print("-" * 65)
    print(f"  Causal SNP:        {corr_causal:.4f}")
    print(f"  Linked SNP 1:       {corr_linked_1:.4f}  <-- Misleading association due to high LD")
    print(f"  Linked SNP 2:       {corr_linked_2:.4f}  <-- Misleading association due to high LD")
    print(f"  Unlinked SNP:       {corr_unlinked:.4f}  <-- No association, as expected")
    print("-" * 65)
    print("\nExplanation:")
    print("The simulation shows a strong correlation for both 'Linked SNP 1' and 'Linked SNP 2' with the trait,")
    print("even though they have no direct causal effect. Their association is 'misleading' because it is")
    print("entirely due to their tight linkage (high LD) with the 'Causal SNP'.")
    print("This perfectly illustrates the scenario described in Answer A.")


if __name__ == '__main__':
    simulate_gwas_ld()
<<<A>>>