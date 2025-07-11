import math

def combinations(n, k):
    """
    Calculates the number of combinations 'n choose k'.
    """
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve():
    """
    Calculates the number of unique autosomal sequences in the F3 generation.
    """
    num_snps = 5
    # The number of intervals between SNPs where a junction can occur.
    num_junction_locations = num_snps - 1

    print("Analyzing the number of unique sequences based on the number of 'junctions' between parental alleles.")
    print("-" * 60)

    # 0 Junctions: The two original parental haplotypes (e.g., 00000 and 11111).
    # These are present in all generations.
    count_0_junc = 2 * combinations(num_junction_locations, 0)

    # 1 Junction: Generated from the F1 cross. (e.g., 00111 or 11000).
    count_1_junc = 2 * combinations(num_junction_locations, 1)

    # 2 Junctions: Generated from F2 crosses. (e.g., 00100 or 11011).
    count_2_junc = 2 * combinations(num_junction_locations, 2)

    # 3 Junctions: Generated from F2 crosses. (e.g., 01010 - oops, length 5 means 01011 or 01101, etc.).
    count_3_junc = 2 * combinations(num_junction_locations, 3)
    
    # 4 Junctions: These haplotypes (e.g., 01010) cannot be created by a single crossover
    # between haplotypes that themselves only have 0 or 1 junctions (the F1 gametes).
    # So, they won't appear in the F2 gamete pool.

    total_sequences = count_0_junc + count_1_junc + count_2_junc + count_3_junc

    print("The total number of unique sequences is the sum of all constructible types:")
    print(f"Parental types (0 junctions)        : {count_0_junc}")
    print(f"F1 recombinant types (1 junction)   : {count_1_junc}")
    print(f"F2 recombinant types (2 junctions)  : {count_2_junc}")
    print(f"F2 recombinant types (3 junctions)  : {count_3_junc}")
    
    print("\nThe final equation is:")
    print(f"Total Unique Sequences = {count_0_junc} + {count_1_junc} + {count_2_junc} + {count_3_junc} = {total_sequences}")


solve()