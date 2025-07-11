import pandas as pd

def simulate_haplotype_association():
    """
    Simulates how Tag SNPs can lead to a misleading association due to LD.
    In this model, a trait is determined ONLY by the "Causal SNP" (SNP 4).
    However, we will test the "Tag SNPs" (SNP 1, 3, 7) for an association.
    """
    # Define two common haplotypes in a population.
    # Each list represents a haplotype with alleles for 8 SNPs (0 or 1).
    haplotype_A = [0, 1, 0, 0, 1, 0, 0, 1]  # Trait-negative haplotype
    haplotype_B = [1, 0, 1, 1, 0, 1, 1, 0]  # Trait-positive haplotype

    # The causal SNP is at index 3 (the 4th SNP).
    # The trait is present if the allele at SNP 4 is 1.
    causal_snp_index = 3

    # We choose Tag SNPs at indices 0, 2, and 6 (SNP 1, 3, and 7).
    tag_snp_indices = [0, 2, 6]

    # Create a simulated population of 100 individuals
    # 50 have Haplotype A, 50 have Haplotype B
    population_data = []
    for i in range(50):
        population_data.append({'haplotype': 'A', 'alleles': haplotype_A})
    for i in range(50):
        population_data.append({'haplotype': 'B', 'alleles': haplotype_B})

    df = pd.DataFrame(population_data)

    # Determine trait status based on the CAUSAL SNP
    df['trait_positive'] = df['alleles'].apply(lambda x: x[causal_snp_index] == 1)

    # Extract the alleles of the Tag SNPs
    for i, snp_idx in enumerate(tag_snp_indices):
        df[f'Tag_SNP_{snp_idx+1}_allele'] = df['alleles'].apply(lambda x: x[snp_idx])

    print("--- Simulation of Misleading Association due to Haplotype Tagging ---")
    print(f"The true CAUSAL SNP is SNP {causal_snp_index + 1}.")
    print(f"The trait is positive if and only if the allele at SNP {causal_snp_index + 1} is 1.\n")

    # Let's check the association for one of the Tag SNPs (SNP 1)
    tag_snp_to_test = f'Tag_SNP_{tag_snp_indices[0]+1}_allele'
    contingency_table = pd.crosstab(df[tag_snp_to_test], df['trait_positive'])

    print(f"Analyzing association between a Tag SNP (SNP {tag_snp_indices[0]+1}) and the trait:")
    print("Contingency Table:")
    print(contingency_table)
    
    allele_0_trait_neg = contingency_table.loc[0, False]
    allele_1_trait_pos = contingency_table.loc[1, True]
    
    print(f"\nResult of Association Test:")
    print(f"When Tag SNP {tag_snp_indices[0]+1} has allele 0, there are {allele_0_trait_neg} individuals without the trait.")
    print(f"When Tag SNP {tag_snp_indices[0]+1} has allele 1, there are {allele_1_trait_pos} individuals with the trait.")
    print("\nThis shows a perfect correlation!")
    print(f"The combination of Tag SNPs ({[i+1 for i in tag_snp_indices]}) perfectly predicts the allele of the Causal SNP ({causal_snp_index+1}).")
    print("This creates a strong, but MISLEADING, association because the Tag SNPs themselves do not cause the trait.")
    print("For a complex polygenic trait, this strong signal from one locus could obscure weaker, true signals from other parts of the genome.")


simulate_haplotype_association()