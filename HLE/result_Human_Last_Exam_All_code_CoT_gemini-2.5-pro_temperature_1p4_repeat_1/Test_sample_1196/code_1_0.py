import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

def demonstrate_ld_misleading_association():
    """
    Simulates a GWAS scenario to show how LD can create misleading associations.
    """
    print("--- Simulating Misleading GWAS Association due to Linkage Disequilibrium (LD) ---")

    # Step 1: Define two haplotypes (sets of linked alleles on a chromosome)
    # Haplotype 1 carries the 'risk' allele ('G') for the true causal SNP.
    # Haplotype 2 carries the 'protective' allele ('T').
    # The other SNPs (SNP1, SNP2, SNP3) are linked and have different alleles on each haplotype.
    risk_haplotype = {'SNP1': 'C', 'SNP2': 'T', 'SNP3': 'A', 'SNP_causal': 'G'}
    protective_haplotype = {'SNP1': 'A', 'SNP2': 'C', 'SNP3': 'G', 'SNP_causal': 'T'}

    # Step 2: Simulate a population of 1000 individuals (2000 haplotypes)
    # Assume the risk haplotype is present at a 20% frequency.
    n_individuals = 1000
    n_haplotypes = n_individuals * 2
    haplotypes = []
    for _ in range(n_haplotypes):
        if np.random.rand() < 0.2:
            haplotypes.append(risk_haplotype)
        else:
            haplotypes.append(protective_haplotype)

    np.random.shuffle(haplotypes)

    # Create diploid genotypes by pairing the haplotypes
    genotypes = []
    for i in range(0, n_haplotypes, 2):
        h1 = haplotypes[i]
        h2 = haplotypes[i+1]
        # Since LD is nearly perfect, recombination is ignored. Genotypes are combined.
        individual = {snp: "".join(sorted(h1[snp] + h2[snp])) for snp in risk_haplotype.keys()}
        genotypes.append(individual)

    df = pd.DataFrame(genotypes)

    # Step 3: Define the trait based ONLY on the causal SNP
    # The trait is present (1) if the causal genotype contains the 'G' allele.
    df['trait'] = df['SNP_causal'].apply(lambda genotype: 1 if 'G' in genotype else 0)

    print(f"\nSimulation setup is complete. The trait is caused ONLY by 'SNP_causal'.")
    print("Now, let's test the association for the other three tightly linked SNPs.\n")

    # Step 4 & 5: Perform association tests for the non-causal SNPs
    tested_snps = ['SNP1', 'SNP2', 'SNP3']
    for snp in tested_snps:
        print(f"--- Testing association for {snp} ---")
        
        # Create a contingency table (e.g., counts of AA, AC, CC vs. trait 0 or 1)
        contingency_table = pd.crosstab(df[snp], df['trait'])
        
        # Perform chi-squared test
        chi2, p_value, dof, expected = chi2_contingency(contingency_table)
        
        print(f"Result for {snp}: p-value = {p_value:.2e}")
        if p_value < 1e-6:
            print(f"Conclusion: A highly significant association is found for {snp}.")
            print("This is a 'misleading' association because this SNP is not causal; it is just linked to the true causal variant.\n")

    print("--- Final Summary ---")
    print("The simulation shows that all three tested SNPs have a strong statistical association with the trait.")
    print("This occurs because they are tightly linked (in high LD) with the real causal SNP.")
    print("This scenario, where multiple SNPs in a block are all associated, makes it hard to pinpoint the true cause and is the essence of choice A.")


if __name__ == '__main__':
    demonstrate_ld_misleading_association()