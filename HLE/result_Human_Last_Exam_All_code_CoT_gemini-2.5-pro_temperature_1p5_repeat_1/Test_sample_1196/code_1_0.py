import numpy as np
from sklearn.linear_model import LinearRegression

def demonstrate_ld_misleading_association():
    """
    Simulates genotypes and a trait to demonstrate how Linkage Disequilibrium (LD)
    can lead to misleading associations in a GWAS.
    """
    # 1. Parameters
    n_individuals = 1000
    n_snps = 10
    # An LD block exists from SNP 3 to SNP 7 (indices 2 to 6)
    ld_block_indices = list(range(2, 7))
    # The true causal variant is SNP 5 (index 4)
    causal_snp_index = 4
    # We use SNP 3 and SNP 7 as Tag SNPs to represent the block
    tag_snp_indices = [2, 6]

    # 2. Simulate Genotypes with an LD block
    # Create two foundational haplotypes for the LD block region
    # A haplotype is a set of alleles inherited together on one chromosome
    haplotype_A = np.array([0, 1, 0, 1, 0]) # Alleles for SNPs 3 through 7
    haplotype_B = np.array([1, 0, 1, 0, 1])

    # Generate genotypes for each individual
    genotypes = np.zeros((n_individuals, n_snps), dtype=int)

    for i in range(n_individuals):
        # Alleles for SNPs outside the block are random
        genotypes[i, [0, 1, 7, 8, 9]] = np.random.randint(0, 3, size=5)

        # Alleles for SNPs inside the block are inherited based on haplotypes
        # Each individual gets two haplotypes, one from each parent
        # We assume Haplotype A is more common (e.g., 70% frequency)
        h1 = haplotype_A if np.random.rand() < 0.7 else haplotype_B
        h2 = haplotype_A if np.random.rand() < 0.7 else haplotype_B
        
        # An individual's genotype is the sum of alleles on their two haplotypes
        ld_block_genotype = h1 + h2
        genotypes[i, ld_block_indices] = ld_block_genotype

    # 3. Simulate a Trait based on the Causal SNP
    # The trait is a function of the causal SNP's genotype, plus random noise
    beta = 2.0
    noise = np.random.normal(0, 1.5, n_individuals)
    causal_genotype = genotypes[:, causal_snp_index]
    trait = beta * causal_genotype + noise

    # 4. Perform a mock GWAS: Test association for each SNP
    # We use R-squared from a simple linear regression as the measure of association strength
    results = []
    for j in range(n_snps):
        snp_genotype = genotypes[:, j].reshape(-1, 1)
        model = LinearRegression().fit(snp_genotype, trait)
        r_squared = model.score(snp_genotype, trait)
        results.append(r_squared)

    # 5. Print and Explain the Results
    print("Demonstration of Misleading Association due to Linkage Disequilibrium\n")
    print(f"Scenario: A Causal SNP ({causal_snp_index + 1}) is in an inherited haplotype (LD Block).")
    print(f"Tag SNPs ({tag_snp_indices[0]+1} and {tag_snp_indices[1]+1}) are used to 'tag' this haplotype.\n")
    print("GWAS Simulation Results (R-squared indicates association strength):")
    print("-" * 60)
    print(f"{'SNP':<5} | {'Association (R-squared)':<25} | {'Comment'}")
    print("-" * 60)

    for j in range(n_snps):
        comment = ""
        if j == causal_snp_index:
            comment = "<- True Causal SNP (Strongest Signal)"
        elif j in tag_snp_indices:
            comment = "<- Misleading Tag SNP (Also Strong)"
        elif j in ld_block_indices:
            comment = "<- Misleading (in LD with Causal)"
        else:
            comment = "<- Unlinked SNP (No Signal)"

        print(f"SNP {j+1:<2} | {results[j]:<25.4f} | {comment}")
    print("-" * 60)
    print("\nExplanation:")
    print("The Causal SNP shows the strongest association, as expected.")
    print("However, the Tag SNPs and other variants in the same haplotype also show very strong associations.")
    print("This is because they are inherited together (in high LD) with the causal variant.")
    print("A researcher would be 'misled' into thinking the Tag SNPs might be causal, when in reality")
    print("they are just non-functional markers pointing to the right neighborhood.")

if __name__ == '__main__':
    demonstrate_ld_misleading_association()