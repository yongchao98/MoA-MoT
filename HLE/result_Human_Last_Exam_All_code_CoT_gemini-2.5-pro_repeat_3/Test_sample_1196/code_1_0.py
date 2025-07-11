import numpy as np
from scipy import stats

def demonstrate_ld_effect():
    """
    Simulates a GWAS scenario to demonstrate how tightly linked SNPs
    (as in Option A) can lead to misleading associations.
    """
    # Step 1: Define haplotypes for a single LD block.
    # A haplotype is a set of alleles on a chromosome that are inherited together.
    # We have one true causal SNP and three tightly linked, non-causal SNPs.
    # risk_haplotype carries the allele that increases trait value.
    #
    # Haplotype format: [Causal_SNP, Linked_SNP_1, Linked_SNP_2, Linked_SNP_3]
    risk_haplotype =     ['T', 'A', 'G', 'C']
    non_risk_haplotype = ['C', 'G', 'T', 'A']
    snp_names = ["True Causal SNP", "Linked SNP 1", "Linked SNP 2", "Linked SNP 3"]

    # Step 2: Simulate a population of 1000 individuals.
    n_individuals = 1000
    population_genotypes = []
    for _ in range(n_individuals):
        # Each individual gets two haplotypes, one from each parent.
        # We assume a 30% frequency for the risk haplotype in the population.
        hap1 = risk_haplotype if np.random.rand() < 0.3 else non_risk_haplotype
        hap2 = risk_haplotype if np.random.rand() < 0.3 else non_risk_haplotype
        population_genotypes.append([hap1, hap2])

    # Step 3: Simulate a quantitative trait influenced ONLY by the causal SNP.
    trait_values = []
    for hap1, hap2 in population_genotypes:
        # Count the number of 'T' alleles at the causal position (index 0).
        risk_allele_count = (hap1[0] == 'T') + (hap2[0] == 'T')
        
        # Base trait value is 10. Each risk allele adds 5 to the value.
        # Add random noise to make the simulation more realistic.
        base_value = 10
        effect_size = 5
        noise = np.random.normal(loc=0, scale=2)
        trait = base_value + (risk_allele_count * effect_size) + noise
        trait_values.append(trait)

    # Step 4: Perform a mock "GWAS". Test each SNP for association with the trait.
    print("--- Simulating Association Test for Tightly Linked SNPs ---")
    print("The true cause of the trait is the 'True Causal SNP'. The others are non-causal but linked.\n")

    for i in range(len(snp_names)):
        risk_allele = risk_haplotype[i]
        
        # For each SNP, we create two groups: individuals with the risk allele vs. without.
        # A more advanced analysis would use a linear regression model, but this illustrates the point.
        group_with_risk_allele = []
        group_without_risk_allele = []

        for j in range(n_individuals):
            hap1, hap2 = population_genotypes[j]
            # Check if the individual carries the risk allele for the current SNP
            if (hap1[i] == risk_allele) or (hap2[i] == risk_allele):
                group_with_risk_allele.append(trait_values[j])
            else:
                group_without_risk_allele.append(trait_values[j])
        
        # Perform an independent t-test. A small p-value indicates a significant association.
        t_stat, p_value = stats.ttest_ind(group_with_risk_allele, group_without_risk_allele, equal_var=False)

        print(f"Testing Association for: '{snp_names[i]}'")
        # The key numbers from the "equation" of the t-test are the statistic and p-value.
        print(f"  T-statistic = {t_stat:.2f}")
        print(f"  P-value = {p_value:.2e}") # Using scientific notation for small p-values
        print("-" * 40)

    print("\nCONCLUSION:")
    print("All four SNPs show a highly significant association (a very small p-value) with the trait.")
    print("This happens because the three non-causal SNPs are in high LD with the true causal SNP.")
    print("A researcher seeing these results would not know which of the four is the true cause without more information.")
    print("This demonstrates how a combination of tightly linked SNPs (Option A) creates a misleading association.")

demonstrate_ld_effect()