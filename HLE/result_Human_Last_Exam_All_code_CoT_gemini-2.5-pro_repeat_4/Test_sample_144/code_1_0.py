import itertools

def solve_recombination():
    """
    Calculates the number of unique autosome sequences in the F3 generation.
    """
    # --- Problem Setup ---
    # 5 SNPs means the haplotype string has length 5.
    # Crossovers can occur in the 4 intervals between SNPs.
    num_snps = 5
    crossover_locations = range(1, num_snps)

    # Represent the two inbred strains' haplotypes.
    # Strain A: 'AAAAA', Strain B: 'BBBBB'
    h_strain_A = 'A' * num_snps
    h_strain_B = 'B' * num_snps

    # --- F1 Gametes -> F2 Haplotypes ---
    # The F1 generation is AAAAA/BBBBB. Their gametes form the F2 haplotype pool.
    # This pool starts with the parental haplotypes.
    f2_haplotypes = {h_strain_A, h_strain_B}
    num_parental_haplotypes = len(f2_haplotypes)

    # A single crossover between AAAAA and BBBBB creates recombinant haplotypes.
    for c in crossover_locations:
        # A crossover event produces two reciprocal recombinant haplotypes.
        recombinant1 = h_strain_A[:c] + h_strain_B[c:]
        recombinant2 = h_strain_B[:c] + h_strain_A[c:]
        f2_haplotypes.add(recombinant1)
        f2_haplotypes.add(recombinant2)

    num_f2_haplotypes = len(f2_haplotypes)

    # --- F2 Gametes -> F3 Haplotypes ---
    # The F3 haplotype pool is formed by recombination between all possible pairs
    # of haplotypes present in the F2 generation.
    # We start with the F2 set and add any new ones.
    f3_haplotypes = set(f2_haplotypes)

    # Generate all unique pairs of haplotypes from the F2 pool for recombination.
    # We use combinations_with_replacement since a haplotype can recombine with itself.
    f2_haplotype_list = sorted(list(f2_haplotypes))
    haplotype_pairs = itertools.combinations_with_replacement(f2_haplotype_list, 2)

    for h1, h2 in haplotype_pairs:
        for c in crossover_locations:
            recombinant1 = h1[:c] + h2[c:]
            recombinant2 = h2[:c] + h1[c:]
            f3_haplotypes.add(recombinant1)
            f3_haplotypes.add(recombinant2)

    num_f3_haplotypes = len(f3_haplotypes)
    num_new_haplotypes = num_f3_haplotypes - num_f2_haplotypes

    # --- Final Output ---
    # The user requested the numbers in the final equation to be output.
    print(f"The process starts with {num_parental_haplotypes} parental haplotypes.")
    print(f"After one generation of recombination, the number of unique haplotypes in the F2 generation is {num_f2_haplotypes}.")
    print(f"These {num_f2_haplotypes} haplotypes can then recombine in the F2 generation.")
    print(f"This second round of recombination generates {num_new_haplotypes} new, unique haplotypes not seen before.")
    print("\nThe final number of unique sequences in the F3 generation is the sum of the haplotypes from the F2 generation and the newly created ones:")
    print(f"{num_f2_haplotypes} (from F2) + {num_new_haplotypes} (new) = {num_f3_haplotypes}")

solve_recombination()