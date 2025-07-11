import itertools

def calculate_unique_sequences():
    """
    Calculates the number of unique autosome sequences possible in the F3 generation
    based on the specified genetic rules.
    """
    
    # Define the number of SNPs on the chromosome.
    num_snps = 5
    
    # Step 1: Represent the parental chromosomes.
    # '0' represents the allele from the first inbred strain (A).
    # '1' represents the allele from the second inbred strain (B).
    parent_A = '0' * num_snps
    parent_B = '1' * num_snps
    
    # Step 2: Generate the set of unique gametes from the F1 generation.
    # The F1 individuals are heterozygous ('00000' / '11111').
    # A single crossover occurs at one of the 4 possible locations between the 5 SNPs.
    f1_gametes = set()
    for i in range(1, num_snps):
        # A crossover between position i-1 and i creates two recombinant types.
        recombinant_1 = parent_A[:i] + parent_B[i:]
        recombinant_2 = parent_B[:i] + parent_A[i:]
        f1_gametes.add(recombinant_1)
        f1_gametes.add(recombinant_2)
        
    # Step 3: Generate the set of unique gametes from the F2 generation.
    # These gametes will form the F3 generation.
    # F2 individuals are formed by any pair of F1 gametes.
    f3_gamete_pool = set()
    f1_gamete_list = list(f1_gametes)
    
    # Iterate over all possible pairs of chromosomes in F2 individuals.
    # itertools.combinations_with_replacement ensures we consider both
    # heterozygous (c1, c2) and homozygous (c1, c1) F2 individuals.
    for c1, c2 in itertools.combinations_with_replacement(f1_gamete_list, 2):
        # For each F2 individual, generate their gametes via a single crossover.
        for i in range(1, num_snps):
            gamete_1 = c1[:i] + c2[i:]
            gamete_2 = c2[:i] + c1[i:]
            f3_gamete_pool.add(gamete_1)
            f3_gamete_pool.add(gamete_2)
            
    # Step 4: The result is the total number of unique sequences generated.
    num_unique_f3_sequences = len(f3_gamete_pool)
    
    # --- Output the results and reasoning ---
    print("Calculation Steps:")
    print(f"1. Number of SNPs: {num_snps}")
    print(f"2. Parental chromosomes: {parent_A} and {parent_B}")
    print(f"3. Number of unique gametes from F1 generation (single crossover products): {len(f1_gametes)}")
    
    # Creating pairs from the 8 F1 gametes:
    num_f1_gametes = len(f1_gametes)
    num_f2_genotypes = (num_f1_gametes * (num_f1_gametes + 1)) // 2

    print(f"4. Number of unique F2 genotypes to consider: {num_f2_genotypes}")
    print(f"5. Generating all possible gametes from these F2 individuals by applying another single crossover.")
    print("\nFinal Result:")
    print(f"The total number of possible unique sequences in the F3 generation is the size of the resulting gamete pool.")
    print(f"Final Count = {num_unique_f3_sequences}")

# Run the calculation and print the output
calculate_unique_sequences()
