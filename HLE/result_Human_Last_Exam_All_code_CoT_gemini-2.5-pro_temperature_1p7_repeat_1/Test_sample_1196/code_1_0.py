import numpy as np
from scipy.stats import chi2_contingency

def simulate_gwas_ld():
    """
    Simulates how tightly linked SNPs can all show a misleading association
    with a trait, even when only one is causal.
    """
    # Define two haplotypes in perfect linkage disequilibrium (LD).
    # A haplotype is a set of SNP alleles on a single chromosome.
    # Let's say there are 4 SNPs in this LD block.
    # Haplotype 1: [Allele for SNP1, Allele for SNP2, ...]
    haplotype_1 = ['G', 'A', 'T', 'C']
    haplotype_2 = ['A', 'C', 'G', 'T']
    
    # Let SNP2 (index 1) be the causal SNP.
    # The 'A' allele of SNP2 causes the trait (dominant model).
    causal_snp_index = 1
    causal_allele = 'A'

    print("GWAS Simulation Demonstrating Linkage Disequilibrium\n")
    print(f"Scenario: 4 SNPs are in a tight LD block.")
    print(f"The two common haplotypes are: {haplotype_1} and {haplotype_2}")
    print(f"Only SNP {causal_snp_index + 1} is causal: allele '{causal_allele}' causes the trait.\n")

    # Generate a population of 1000 individuals
    num_individuals = 1000
    population = []
    trait_status = []

    # Each individual gets two haplotypes, one from each "parent".
    for _ in range(num_individuals):
        # In a real population, haplotype frequencies would vary.
        # For simplicity, we'll draw them with equal probability.
        h1_idx = np.random.randint(0, 2)
        h2_idx = np.random.randint(0, 2)
        h1 = haplotype_1 if h1_idx == 0 else haplotype_2
        h2 = haplotype_1 if h2_idx == 0 else haplotype_2
        
        individual_genotype = (h1, h2)
        population.append(individual_genotype)

        # Determine trait status based on the causal SNP
        allele1 = h1[causal_snp_index]
        allele2 = h2[causal_snp_index]
        if allele1 == causal_allele or allele2 == causal_allele:
            trait_status.append(1) # Case (has trait)
        else:
            trait_status.append(0) # Control (no trait)

    # Now, perform an association test for EACH SNP
    print("Performing association test (Chi-Squared) for each SNP:")
    print("-" * 55)
    print(f"{'SNP':<5} | {'Alleles':<12} | {'Chi-Squared':<14} | {'P-Value':<15} | {'Note'}")
    print("-" * 55)

    num_snps = len(haplotype_1)
    for i in range(num_snps):
        # Get the two possible alleles for the current SNP
        alleles = np.unique([h[i] for h in [haplotype_1, haplotype_2]])
        
        # Count allele occurrences in cases and controls
        # Contingency Table:
        #          Allele_1 | Allele_2
        #  Case    |          |
        #  Control |          |
        counts = {
            'case': {alleles[0]: 0, alleles[1]: 0},
            'control': {alleles[0]: 0, alleles[1]: 0}
        }

        for j, individual in enumerate(population):
            allele1 = individual[0][i]
            allele2 = individual[1][i]
            if trait_status[j] == 1: # Case
                counts['case'][allele1] += 1
                counts['case'][allele2] += 1
            else: # Control
                counts['control'][allele1] += 1
                counts['control'][allele2] += 1
        
        contingency_table = [
            [counts['case'][alleles[0]], counts['case'][alleles[1]]],
            [counts['control'][alleles[0]], counts['control'][alleles[1]]]
        ]

        # Perform chi-squared test
        chi2, p, _, _ = chi2_contingency(contingency_table)

        note = "Causal" if i == causal_snp_index else "Misleading Signal"
        allele_str = f"{alleles[0]}/{alleles[1]}"
        print(f"SNP {i+1:<4} | {allele_str:<12} | {chi2:<14.2f} | {p:<15.2e} | {note}")

    print("-" * 55)
    print("\nConclusion: Because all SNPs are in high LD, they are all strongly")
    print("associated with the trait, not just the true causal SNP.")
    print("This demonstrates how a signal from a non-causal SNP can be misleading.")

if __name__ == '__main__':
    simulate_gwas_ld()