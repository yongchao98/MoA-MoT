import itertools

def solve_chromosome_recombination():
    """
    Calculates the number of unique chromosome sequences in the F3 generation.
    """
    num_snps = 5

    # Step 1: Define the parental sequences (P generation)
    strain_A = "A" * num_snps
    strain_B = "B" * num_snps
    
    print(f"Parental (P) sequences: {strain_A}, {strain_B}")

    # Step 2: Generate the set of unique gametes from the F1 generation.
    # An F1 individual is AAAAA/BBBBB. Its gametes include parental types
    # (no crossover between SNPs) and recombinant types (one crossover).
    f1_gametes = {strain_A, strain_B}
    for i in range(num_snps - 1):
        # A single crossover occurs at one of the 4 locations between the 5 SNPs.
        # Location 'i' is between SNP i and SNP i+1.
        recombinant_1 = strain_A[:i+1] + strain_B[i+1:]
        recombinant_2 = strain_B[:i+1] + strain_A[i+1:]
        f1_gametes.add(recombinant_1)
        f1_gametes.add(recombinant_2)

    num_f1_gametes = len(f1_gametes)
    print(f"Number of unique gamete sequences from F1 generation: {num_f1_gametes}")

    # Step 3: Generate the set of unique gametes from the F2 generation.
    # These are the sequences found in the F3 generation.
    # F2 individuals are formed from pairs of F1 gametes. We must consider
    # recombination between all possible pairs of F1 gametes.
    f2_gametes = f1_gametes.copy()
    
    # We iterate through all pairs of F1 gametes, including pairing a sequence with itself.
    # The `itertools.product` function is suitable for this.
    for c1, c2 in itertools.product(f1_gametes, repeat=2):
        for i in range(num_snps - 1):
            recombinant_1 = c1[:i+1] + c2[i+1:]
            recombinant_2 = c2[:i+1] + c1[i+1:]
            f2_gametes.add(recombinant_1)
            f2_gametes.add(recombinant_2)

    # Step 4: The final answer is the total number of unique sequences generated.
    num_f2_gametes = len(f2_gametes)
    
    # The final equation is that the number of possible sequences is 2 to the power of the number of SNPs.
    final_equation_result = 2**num_snps

    print(f"Total possible unique sequences with {num_snps} SNPs is 2^{num_snps} = {final_equation_result}")
    print(f"Number of unique sequences found in the F3 generation: {num_f2_gametes}")

solve_chromosome_recombination()
<<<32>>>