import itertools

def solve():
    """
    Calculates the number of unique chromosome sequences possible in the F3 generation.
    """
    # Step 1: Define the set of unique gametes from the F1 generation (S1).
    # These are the parental types and all single-crossover recombinants between them.
    parental_a = "AAAAA"
    parental_b = "BBBBB"
    
    f1_gametes = {parental_a, parental_b}
    num_snps = 5
    # There are num_snps - 1 = 4 regions for crossover.
    for i in range(1, num_snps):
        recombinant_1 = parental_a[:i] + parental_b[i:]
        recombinant_2 = parental_b[:i] + parental_a[i:]
        f1_gametes.add(recombinant_1)
        f1_gametes.add(recombinant_2)

    # This is the set of chromosomes present in the F2 population.
    # Gametes produced by F2 individuals form the basis for the F3 generation.
    # We need to find the closure of this set under recombination.
    possible_sequences = f1_gametes.copy()
    
    print(f"Number of unique sequences in F1 gametes = {len(possible_sequences)}")

    iteration = 1
    while True:
        # We'll collect new sequences found in this pass.
        newly_generated = set()
        
        # Take all pairs of sequences from the current set.
        current_list = list(possible_sequences)
        for s1, s2 in itertools.product(current_list, repeat=2):
            # Recombine them at every possible locus.
            for i in range(1, num_snps):
                recombinant_1 = s1[:i] + s2[i:]
                recombinant_2 = s2[:i] + s1[i:]
                
                # If a new sequence is found, add it to our temporary set.
                if recombinant_1 not in possible_sequences:
                    newly_generated.add(recombinant_1)
                if recombinant_2 not in possible_sequences:
                    newly_generated.add(recombinant_2)
        
        # If new sequences were generated, add them to the main set and continue.
        if newly_generated:
            possible_sequences.update(newly_generated)
            print(f"After {iteration} round(s) of recombination, the number of unique sequences is {len(possible_sequences)}")
            iteration += 1
        else:
            # If no new sequences are found, the set is closed.
            print("The set of sequences is now closed.")
            break
            
    print("\n--- Final Calculation ---")
    print(f"The number of unique sequences possible in the F3 generation is the size of the closed set.")
    print(f"Final Count: {len(possible_sequences)}")

solve()
<<<32>>>