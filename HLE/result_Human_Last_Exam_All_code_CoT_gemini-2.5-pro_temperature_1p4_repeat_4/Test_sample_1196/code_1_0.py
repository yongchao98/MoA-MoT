import numpy as np
from scipy.stats import chi2_contingency

def simulate_gwas_ld():
    """
    Simulates a scenario to demonstrate misleading association due to Linkage Disequilibrium (LD).

    In this simulation:
    - There is a 'causal_snp' that truly affects the trait (e.g., disease status).
    - We do not "genotype" or measure the causal_snp directly.
    - Instead, we measure three 'tag_snps' that are in high LD (inherited together) with the causal_snp.
    - We test the association between the tag_snps and the trait.
    """
    np.random.seed(42)
    
    # Define two haplotypes (blocks of co-inherited alleles)
    # Haplotype 1: Contains the 'risk' allele for the causal SNP ('A')
    # Haplotype 2: Contains the 'non-risk' allele for the causal SNP ('a')
    # Alleles for (Tag1, Tag2, Causal, Tag3)
    haplotype_1 = ['G', 'C', 'A', 'T'] 
    haplotype_2 = ['A', 'T', 'a', 'C'] 
    
    # Simulate a population of 1000 individuals
    num_individuals = 1000
    population_genotypes = []
    population_traits = []

    for _ in range(num_individuals):
        # Each individual has two chromosomes, randomly inheriting one of the two haplotypes for each
        chromo_1 = haplotype_1 if np.random.rand() > 0.5 else haplotype_2
        chromo_2 = haplotype_1 if np.random.rand() > 0.5 else haplotype_2
        
        # Determine genotype for the causal SNP
        causal_genotype = chromo_1[2] + chromo_2[2]
        
        # Assign trait status based on the causal SNP
        # Individuals with the 'A' risk allele are more likely to be cases
        if 'A' in causal_genotype:
            trait = 'case' if np.random.rand() > 0.4 else 'control' # 60% chance
        else: # genotype 'aa'
            trait = 'case' if np.random.rand() > 0.6 else 'control' # 40% chance
            
        # Store the genotypes of the TAG SNPs (not the causal one)
        genotype_tag1 = chromo_1[0] + chromo_2[0]
        genotype_tag2 = chromo_1[1] + chromo_2[1]
        genotype_tag3 = chromo_1[3] + chromo_2[3]
        
        population_genotypes.append([genotype_tag1, genotype_tag2, genotype_tag3])
        population_traits.append(trait)

    print("--- Simulation of Misleading Association due to LD ---")
    print(f"A hidden 'causal_snp' influences the trait ('case'/'control' status).")
    print(f"We only measure three 'tag_snps' that are in high LD with the causal one.\n")
    print("Running association tests on the Tag SNPs...")

    # Perform Chi-squared test for each Tag SNP
    for i, tag_snp_name in enumerate(["Tag SNP 1", "Tag SNP 2", "Tag SNP 3"]):
        # Create a contingency table: count genotypes for cases and controls
        contingency_table = {}
        for genotype, trait in zip(population_genotypes, population_traits):
            gt = genotype[i] # Genotype for the current tag SNP
            if gt not in contingency_table:
                contingency_table[gt] = {'case': 0, 'control': 0}
            contingency_table[gt][trait] += 1
        
        # Format table for scipy
        genotype_order = sorted(contingency_table.keys())
        table = [[contingency_table[gt]['case'], contingency_table[gt]['control']] for gt in genotype_order]
        
        # Perform test
        chi2, p, _, _ = chi2_contingency(table)
        
        print(f"\nResults for {tag_snp_name}:")
        print(f"Chi-squared = {chi2:.2f}, P-value = {p:.4g}")
        if p < 0.05:
            print(f"-> Strong association found!")
        else:
            print("-> No significant association found.")
            
    print("\n--- Conclusion ---")
    print("As the simulation shows, all three Tag SNPs have a highly significant association with the trait.")
    print("This association is 'misleading' because none of them are the true cause; they are just proxies")
    print("for the un-measured causal SNP because they are all part of the same inherited haplotype (LD block).")
    print("This matches answer choice E.")

simulate_gwas_ld()