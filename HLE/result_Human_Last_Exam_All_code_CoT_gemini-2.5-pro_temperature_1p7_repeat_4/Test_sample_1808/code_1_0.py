# My plan is to demonstrate the Wahlund effect, which occurs when you mix
# genetically distinct populations, as in a hybrid zone. This effect leads to a
# high Fis (inbreeding coefficient), showing that option C is possible.

# Step 1: Define two parental populations at a single locus (Alleles 'A' and 'a').
# Population 1 is fixed for allele 'A'.
# Population 2 is fixed for allele 'a'.
pop1_genotypes = ['AA'] * 100  # 100 individuals, all homozygous AA
pop2_genotypes = ['aa'] * 100  # 100 individuals, all homozygous aa

# Step 2: Create a sample from the "hybrid zone" by mixing the two populations.
mixed_sample = pop1_genotypes + pop2_genotypes
total_individuals = len(mixed_sample)

# Step 3: Calculate the allele frequencies (p and q) in the mixed sample.
# Each individual has 2 alleles.
num_A_alleles = mixed_sample.count('AA') * 2
num_a_alleles = mixed_sample.count('aa') * 2
total_alleles = total_individuals * 2

p = num_A_alleles / total_alleles  # Frequency of 'A'
q = num_a_alleles / total_alleles  # Frequency of 'a'

# Step 4: Calculate the expected heterozygosity (He) under Hardy-Weinberg Equilibrium.
# He is the proportion of heterozygotes expected if the population were randomly mating.
expected_heterozygosity_He = 2 * p * q

# Step 5: Calculate the observed heterozygosity (Ho).
# Ho is the actual proportion of heterozygotes in our sample.
observed_heterozygotes_Ho = mixed_sample.count('Aa') / total_individuals

# Step 6: Calculate Fis (the inbreeding coefficient).
# Fis = (He - Ho) / He
# A high positive Fis indicates a deficit of heterozygotes.
Fis = (expected_heterozygosity_He - observed_heterozygotes_Ho) / expected_heterozygosity_He

# Step 7: Print the results to show how high Fis can be.
print("--- Demonstrating the Wahlund Effect in a Hybrid Zone ---")
print(f"Parental Population 1 is 100% 'AA'")
print(f"Parental Population 2 is 100% 'aa'")
print(f"A mixed sample of {total_individuals} individuals is created.\n")
print("CALCULATIONS FOR THE MIXED SAMPLE:")
print(f"Allele 'A' frequency (p): {p:.2f}")
print(f"Allele 'a' frequency (q): {q:.2f}")
print(f"Observed Heterozygosity (Ho): {observed_heterozygotes_Ho:.4f}")
print(f"Expected Heterozygosity (He = 2*p*q): {expected_heterozygosity_He:.4f}")
print("\nThe resulting inbreeding coefficient (Fis) is calculated as (He - Ho) / He:")
print(f"Fis = ({expected_heterozygosity_He:.2f} - {observed_heterozygotes_Ho:.2f}) / {expected_heterozygosity_He:.2f}")
print(f"Final Fis value = {Fis:.2f}\n")
print("This shows a very high Fis value of 1.00, demonstrating that high Fis (Option C) can definitely occur when populations are mixed, as they are in a hybrid zone.")
