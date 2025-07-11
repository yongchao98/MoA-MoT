import numpy as np
from scipy.stats import pearsonr

# This script simulates a scenario to demonstrate how linkage disequilibrium (LD)
# can lead to misleading associations in a genome-wide association study (GWAS).

# --- Plan ---
# 1. Define a "Causal SNP" (SNP_C) that directly influences a quantitative trait.
# 2. Define three "Tightly Linked SNPs" (SNP_A1, A2, A3) that are highly correlated with the Causal SNP, simulating scenario A.
# 3. Define an "Unlinked SNP" (SNP_U) that is not correlated with the Causal SNP or the trait.
# 4. Generate a population and their trait values, determined by the Causal SNP plus random noise.
# 5. Calculate the association score (correlation) for each SNP with the trait.
# 6. Print the results to show that the Tightly Linked SNPs have misleadingly high association scores,
#    making it appear as if they are all causing the trait.

def run_gwas_simulation():
    """
    Runs a simulation to show misleading GWAS associations due to tight linkage.
    """
    print("--- Simulation of GWAS Association ---")
    print("One SNP, the 'Causal SNP', truly influences a trait.")
    print("Three other SNPs are 'Tightly Linked' to it in an LD block (Scenario A).")
    print("An 'Unlinked SNP' has no connection.")
    print("We will calculate an association score (correlation) for each SNP with the trait.\n")

    # Number of individuals in the study
    n_individuals = 5000

    # Alleles are coded as 0 and 1.
    # The '1' allele of the Causal SNP will increase the trait value.
    causal_snp = np.random.randint(0, 2, n_individuals)

    # These three SNPs are tightly linked to the Causal SNP (99% correlation)
    ld_strength = 0.99
    linked_snp_1 = np.where(np.random.rand(n_individuals) < ld_strength, causal_snp, 1 - causal_snp)
    linked_snp_2 = np.where(np.random.rand(n_individuals) < ld_strength, causal_snp, 1 - causal_snp)
    linked_snp_3 = np.where(np.random.rand(n_individuals) < ld_strength, causal_snp, 1 - causal_snp)

    # This SNP is unlinked (genetically independent)
    unlinked_snp = np.random.randint(0, 2, n_individuals)

    # The trait is determined by the Causal SNP's effect plus some random biological noise
    trait_value = 100 + (causal_snp * 10) + np.random.normal(0, 5, n_individuals)

    # Calculate association scores (absolute Pearson correlation)
    # The 'r' value measures the strength of the linear association.
    score_causal, _ = pearsonr(causal_snp, trait_value)
    score_linked_1, _ = pearsonr(linked_snp_1, trait_value)
    score_linked_2, _ = pearsonr(linked_snp_2, trait_value)
    score_linked_3, _ = pearsonr(linked_snp_3, trait_value)
    score_unlinked, _ = pearsonr(unlinked_snp, trait_value)

    print("--- Results: Association Scores ---")
    print("A higher score means a stronger statistical association with the trait.")
    print("Notice how the linked, non-causal SNPs have scores almost as high as the true cause.\n")
    print(f"Association(Causal SNP, Trait)      = {abs(score_causal):.4f}")
    print("-------------------------------------------------")
    print("Association(Linked SNP 1, Trait)    = {abs(score_linked_1):.4f}  <-- Misleading Signal")
    print("Association(Linked SNP 2, Trait)    = {abs(score_linked_2):.4f}  <-- Misleading Signal")
    print("Association(Linked SNP 3, Trait)    = {abs(score_linked_3):.4f}  <-- Misleading Signal")
    print("-------------------------------------------------")
    print(f"Association(Unlinked SNP, Trait)    = {abs(score_unlinked):.4f}  <-- No Signal")

if __name__ == '__main__':
    run_gwas_simulation()