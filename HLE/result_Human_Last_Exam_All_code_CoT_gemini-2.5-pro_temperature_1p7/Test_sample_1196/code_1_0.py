import numpy as np
from scipy.stats import pearsonr

def simulate_gwas_ld_scenario():
    """
    Simulates a scenario to demonstrate how Tag SNPs can lead to a misleading GWAS association
    due to Linkage Disequilibrium (LD).
    """
    population_size = 2000
    # Define a haplotype block from SNP index 10 to 29
    block_start, block_end = 10, 30
    
    # Let's say the SNP at index 20 is the true CAUSAL variant, but it's not genotyped.
    causal_snp_idx = 20
    
    # We only genotype three Tag SNPs within this block.
    tag_snp_indices = [12, 18, 25]

    print(f"Scenario Setup:")
    print(f"A hidden CAUSAL SNP exists at genomic position {causal_snp_idx}.")
    print(f"We only measure three Tag SNPs at positions: {tag_snp_indices}\n")
    
    # Create two primary haplotypes for the block. In a high LD block,
    # these patterns are rarely broken up by recombination.
    haplotype_A = np.zeros(block_end - block_start)  # e.g., all '0' alleles
    haplotype_B = np.ones(block_end - block_start)   # e.g., all '1' alleles
    
    # Simulate genotypes for the population
    # Each person gets two chromosomes, so two haplotypes.
    # We will assume a diploid organism and sum the alleles.
    genotypes = np.zeros((population_size, block_end - block_start))
    for i in range(population_size):
        # With high probability (98%), inherit a clean haplotype
        if np.random.rand() < 0.98:
            # Choose Haplotype A or B for each of the two parental chromosomes
            h1 = haplotype_A if np.random.rand() < 0.5 else haplotype_B
            h2 = haplotype_A if np.random.rand() < 0.5 else haplotype_B
        else:
            # With low probability, a recombination event creates a mixed haplotype
            h1 = np.random.randint(0, 2, size=block_end - block_start)
            h2 = np.random.randint(0, 2, size=block_end - block_start)
        genotypes[i, :] = h1 + h2 # Summing alleles for a diploid genotype
        
    # Get the allele counts for the specific SNPs we care about
    causal_snp_alleles = genotypes[:, causal_snp_idx - block_start]
    
    # Generate a trait that is directly influenced by the causal SNP, plus some noise
    noise = np.random.normal(0, 0.5, population_size)
    trait = 2 * causal_snp_alleles + noise
    
    print("Results:")
    print("We now test the association between the trait and our measured Tag SNPs.")
    print("-" * 50)
    
    # Check the association (correlation) for each Tag SNP
    for tag_idx in tag_snp_indices:
        tag_snp_alleles = genotypes[:, tag_idx - block_start]
        correlation, p_value = pearsonr(tag_snp_alleles, trait)
        
        print(f"Association for Tag SNP at position {tag_idx}:")
        print(f"  - Pearson Correlation: {correlation:.4f}")
        print(f"  - P-value: {p_value:0.3e}")
        # An extremely low p-value indicates a statistically significant association
        
        if correlation > 0.8:
            print("  - Verdict: Strong and significant association found!\n")
        else:
            print("  - Verdict: No strong association.\n")
            
    print("-" * 50)
    print("Conclusion:")
    print("Even though the trait is only directly caused by the SNP at position 20,")
    print("the Tag SNPs show a very strong statistical association. This is because")
    print("they are in high Linkage Disequilibrium with the causal SNP, making them")
    print("excellent predictors of the causal allele's presence. The association is real,")
    print("but it is MISLEADING because it points to the Tag SNPs, not the true cause.")

simulate_gwas_ld_scenario()
<<<E>>>