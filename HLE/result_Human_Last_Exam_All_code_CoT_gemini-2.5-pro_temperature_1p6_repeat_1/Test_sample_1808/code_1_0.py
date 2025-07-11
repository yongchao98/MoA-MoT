import random

def calculate_pop_stats(genotypes, population_name):
    """
    Calculates and prints population genetics statistics for a list of genotypes.
    """
    n_total = len(genotypes)
    if n_total == 0:
        print(f"\n--- Stats for {population_name} ---")
        print("Population is empty.")
        return

    # Count genotypes
    n_AA = genotypes.count('AA')
    n_Aa = genotypes.count('Aa')
    n_aa = genotypes.count('aa')

    # Calculate allele frequencies
    # p is the frequency of 'A', q is the frequency of 'a'
    p = (2 * n_AA + n_Aa) / (2 * n_total)
    q = 1 - p

    # Calculate observed and expected heterozygosity
    h_obs = n_Aa / n_total
    h_exp = 2 * p * q

    # Calculate Fis
    if h_exp == 0:
        # If no heterozygotes are expected, Fis is typically considered 0 or undefined.
        # A non-zero h_obs would mean Fis is negative infinity,
        # but this case (h_exp=0, h_obs > 0) is impossible.
        fis = 0
    else:
        fis = (h_exp - h_obs) / h_exp

    print(f"\n--- Stats for {population_name} ---")
    print(f"Genotype Counts: AA={n_AA}, Aa={n_Aa}, aa={n_aa}")
    print(f"Total Individuals: {n_total}")
    print(f"Allele Frequencies: p(A)={p:.3f}, q(a)={q:.3f}")
    print(f"Observed Heterozygosity (H_obs): {h_obs:.3f}")
    print(f"Expected Heterozygosity (H_exp): {h_exp:.3f}")
    print(f"Inbreeding Coefficient (Fis): {fis:.3f}")
    if fis > 0.1:
        print("Result: High positive Fis, indicating a deficit of heterozygotes.")
    elif fis < -0.1:
        print("Result: Negative Fis, indicating an excess of heterozygotes.")
    else:
        print("Result: Fis is near zero, as expected in a randomly mating population.")


# --- Step 1: Define two distinct parental populations ---
P1 = ['AA'] * 50  # Population 1 is fixed for the 'A' allele
P2 = ['aa'] * 50  # Population 2 is fixed for the 'a' allele

# --- Step 2: Simulate a structured population by mechanically mixing P1 and P2 ---
# This is known as the Wahlund effect. There is population structure but no gene flow (no interbreeding).
mixed_population = P1 + P2
calculate_pop_stats(mixed_population, "Scenario 1: Mixed Population (Wahlund Effect)")
print("\nExplanation: When distinct populations are pooled without interbreeding,")
print("there is a strong deficit of heterozygotes, leading to a high positive Fis.")

# --- Step 3: Simulate the direct result of gene flow (interbreeding) ---
# For simplicity, we create an F1 generation where each individual is an offspring
# of one parent from P1 and one from P2.
f1_gene_flow_population = []
for _ in range(50):
    # One gamete 'A' from P1 and one gamete 'a' from P2
    f1_gene_flow_population.append('Aa')

calculate_pop_stats(f1_gene_flow_population, "Scenario 2: F1 Offspring from Gene Flow")
print("\nExplanation: The act of gene flow (interbreeding between distinct populations)")
print("creates an abundance of heterozygotes, leading to a negative Fis.")

print("\nFinal Conclusion: Since gene flow involves interbreeding that creates heterozygotes,")
print("it is fundamentally incompatible with a high positive Fis, which signals a deficit of heterozygotes.")
