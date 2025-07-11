import numpy as np

# Step 1: Define parameters and a population
num_individuals = 1000

# Step 2: Simulate Haplotypes based on Linkage Disequilibrium (LD)
# We define two fundamental haplotypes. In a real population, there would be more variation.
# Haplotype 1 (Causal): The causal variant (C) is 1, and the tag SNPs (T1, T2, T3) are also 1.
# This represents a block of alleles in perfect LD.
haplotype_causal =    {'C': 1, 'T1': 1, 'T2': 1, 'T3': 1}
# Haplotype 0 (Non-causal): All corresponding alleles are 0.
haplotype_non_causal = {'C': 0, 'T1': 0, 'T2': 0, 'T3': 0}

haplotypes = [haplotype_causal, haplotype_non_causal]

# Create individuals by giving them two randomly chosen haplotypes
population = []
for _ in range(num_individuals):
    # Each individual gets two haplotypes
    h1 = haplotypes[np.random.randint(0, 2)]
    h2 = haplotypes[np.random.randint(0, 2)]
    population.append((h1, h2))

# Step 3: Determine Genotypes and Phenotype for each individual
genotypes = {}
phenotypes = []

for snp in ['C', 'T1', 'T2', 'T3']:
    genotypes[snp] = []

for h1, h2 in population:
    # Genotype is the sum of alleles from both haplotypes (0, 1, or 2)
    for snp in ['C', 'T1', 'T2', 'T3']:
        genotypes[snp].append(h1[snp] + h2[snp])
    
    # CRUCIAL STEP: Phenotype is determined *only* by the causal SNP 'C'
    trait_value = h1['C'] + h2['C']
    phenotypes.append(trait_value)

phenotypes = np.array(phenotypes)

print("Simulating a GWAS scenario based on Option E...\n")
print(f"The trait is DIRECTLY caused by SNP 'C'. A higher genotype at 'C' leads to a higher trait value.")
print("Tag SNPs 'T1', 'T2', and 'T3' have NO direct effect on the trait.")
print("However, they are in perfect Linkage Disequilibrium (LD) with 'C'.\n")
print("Let's test the association for each SNP:\n")

# Step 4: Perform and print association tests for each SNP
# We calculate the average phenotype for each genotype group (0, 1, 2)
snp_list = {'C': 'Causal SNP', 'T1': 'Tag SNP 1', 'T2': 'Tag SNP 2', 'T3': 'Tag SNP 3'}

for snp, snp_name in snp_list.items():
    print(f"--- Association Test for {snp_name} ({snp}) ---")
    current_genotypes = np.array(genotypes[snp])
    
    # Calculate average phenotype for genotype 0 (e.g., c/c or t1/t1)
    mean_pheno_0 = np.mean(phenotypes[current_genotypes == 0])
    
    # Calculate average phenotype for genotype 1 (e.g., C/c or T1/t1)
    mean_pheno_1 = np.mean(phenotypes[current_genotypes == 1])
    
    # Calculate average phenotype for genotype 2 (e.g., C/C or T1/T1)
    mean_pheno_2 = np.mean(phenotypes[current_genotypes == 2])
    
    print(f"Average trait for genotype 0: {mean_pheno_0:.2f}")
    print(f"Average trait for genotype 1: {mean_pheno_1:.2f}")
    print(f"Average trait for genotype 2: {mean_pheno_2:.2f}")
    
    # Check if the association is statistically perfect in our simulation
    if (round(mean_pheno_0) == 0 and round(mean_pheno_1) == 1 and round(mean_pheno_2) == 2):
        print("Result: A perfect association is found.\n")
    else:
        # This case shouldn't happen in this perfect simulation
        print("Result: No clear association found.\n")

print("Conclusion: The association signal for the non-causal Tag SNPs (T1, T2, T3) is identical to the signal for the true Causal SNP (C).")
print("An investigator only seeing the test results would not know which SNP is causal.")
print("This is a 'misleading association' because the tag SNPs appear to be causal when they are merely proxies for the true causal variant due to LD.")
<<<E>>>