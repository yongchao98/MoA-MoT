from math import comb

def solve_genetics_problem():
    """
    Calculates the number of unique autosome sequences in the F3 generation.
    """
    # The number of SNPs on the autosome.
    n_snps = 5
    
    # The number of intervals between SNPs where a crossover can occur.
    n_intervals = n_snps - 1

    # --- Step 1: Count sequences with 0 transitions ---
    # These are the two original parental haplotypes (e.g., AAAAA and BBBBB).
    num_0_trans = 2
    print(f"Number of 0-transition sequences (parentals): {num_0_trans}")

    # --- Step 2: Count sequences with 1 transition ---
    # A single crossover can occur in any of the 'n_intervals'.
    # For each location, the sequence can start with either parental type (A or B).
    # e.g., AABBB or BBAAA.
    num_1_trans = n_intervals * 2
    print(f"Number of 1-transition sequences: {n_intervals} locations * 2 starting types = {num_1_trans}")

    # --- Step 3: Count sequences with 2 transitions ---
    # These are formed by recombining 1-transition sequences in the F2 generation.
    # We need to choose 2 distinct intervals for the crossovers to occur.
    # The number of ways to choose 2 intervals from 'n_intervals' is C(n_intervals, 2).
    # For each pair of locations, the sequence can start with either parental type (A or B).
    # e.g., AABAA or BBABB.
    num_2_trans_locations = comb(n_intervals, 2)
    num_2_trans = num_2_trans_locations * 2
    print(f"Number of 2-transition sequences: C({n_intervals}, 2) * 2 starting types = {num_2_trans_locations} * 2 = {num_2_trans}")
    
    # --- Step 4: Calculate the total number of unique sequences ---
    # The set of sequences in the F3 generation is the union of all sequences with 0, 1, or 2 transitions.
    total_sequences = num_0_trans + num_1_trans + num_2_trans
    
    print("\nThe total number of possible unique sequences is the sum of these counts.")
    print(f"Total = {num_0_trans} (0-transition) + {num_1_trans} (1-transition) + {num_2_trans} (2-transition)")
    print(f"Final Answer: {total_sequences}")

solve_genetics_problem()