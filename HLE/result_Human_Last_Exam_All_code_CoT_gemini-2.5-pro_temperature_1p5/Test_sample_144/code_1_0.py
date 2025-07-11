import collections

# Step 1: Define problem parameters
n_snps = 5
parental_a = "0" * n_snps
parental_b = "1" * n_snps
# A crossover can happen after SNP 1, 2, 3, or 4.
# These correspond to string indices 1, 2, 3, and 4.
crossover_indices = range(1, n_snps)

# Helper function to perform a single crossover
def perform_crossover(chrom1, chrom2, index):
    """Performs a single crossover, returning the two resulting chromosomes."""
    gamete1 = chrom1[:index] + chrom2[index:]
    gamete2 = chrom2[:index] + chrom1[index:]
    return gamete1, gamete2

# Helper function to count recombination junctions in a sequence
def count_junctions(sequence):
    """Counts the number of times adjacent bits are different."""
    junctions = 0
    for i in range(len(sequence) - 1):
        if sequence[i] != sequence[i+1]:
            junctions += 1
    return junctions

# Step 2: Generate the gametes from the F1 generation.
# These are the chromosomes that will be found in F2 individuals.
f1_gametes = set()
for idx in crossover_indices:
    g1, g2 = perform_crossover(parental_a, parental_b, idx)
    f1_gametes.add(g1)
    f1_gametes.add(g2)

# This list is the pool of chromosomes available to form F2 individuals.
f2_chromosome_pool = list(f1_gametes)

# Step 3: Generate all possible gametes produced by the F2 generation.
# These are the unique sequences that will be found in the F3 generation.
f3_sequences = set()

# An F2 individual's genotype is (chrom1, chrom2), where both are from the F2 pool.
# We iterate through all possible unique pairs of chromosomes.
for i in range(len(f2_chromosome_pool)):
    for j in range(i, len(f2_chromosome_pool)):
        chrom1 = f2_chromosome_pool[i]
        chrom2 = f2_chromosome_pool[j]
        
        # For each F2 parental pair, find all gametes they can produce.
        for idx in crossover_indices:
            g1, g2 = perform_crossover(chrom1, chrom2, idx)
            f3_sequences.add(g1)
            f3_sequences.add(g2)

# Step 4: Classify the resulting sequences by their number of junctions
junction_counts = collections.defaultdict(int)
for seq in f3_sequences:
    j = count_junctions(seq)
    junction_counts[j] += 1

# Step 5: Output the results and the final equation
print("The analysis of possible chromosome sequences in the F3 generation yields the following breakdown:")
print("-" * 70)
print(f"Number of parental-type sequences (0 junctions): {junction_counts[0]}")
print(f"Number of single-recombinant type sequences (1 junction): {junction_counts[1]}")
print(f"Number of double-recombinant type sequences (2 junctions): {junction_counts[2]}")
print(f"Number of triple-recombinant type sequences (3 junctions): {junction_counts[3]}")
print(f"Number of quadruple-recombinant type sequences (4 junctions): {junction_counts[4]}")
print("-" * 70)

# Build the final equation string as requested
equation_parts = []
# Iterate through sorted keys (0, 1, 2...) for a clean equation
for num_junctions in sorted(junction_counts.keys()):
    count = junction_counts[num_junctions]
    if count > 0:
        equation_parts.append(str(count))

total_count = len(f3_sequences)
equation_str = " + ".join(equation_parts) + f" = {total_count}"

print(f"The total number of unique sequences is derived from the sum of sequence types:")
print(f"Final Equation: {equation_str}")
<<<22>>>