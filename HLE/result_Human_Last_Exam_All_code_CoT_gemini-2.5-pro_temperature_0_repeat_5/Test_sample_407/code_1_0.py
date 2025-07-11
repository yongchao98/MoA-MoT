def calculate_gene_flow_cost():
    """
    This function simulates the calculation of the cost of gene flow in yeast
    by determining the selection coefficient (s) against F2 recombinants.

    The cost of gene flow (outbreeding depression) is often revealed after meiosis
    in the F2 generation. We measure this by comparing the fitness of the F2
    recombinants to the original parental lines.

    The selection coefficient (s) is calculated as: s = 1 - w
    where 'w' is the relative fitness of the genotype being selected against.
    """

    # Step 1: Define hypothetical fitness values (e.g., measured as competitive growth rate).
    # The parental line is the reference and is considered to have optimal fitness in its environment.
    fitness_parental = 1.0
    
    # The F2 generation results from mating F1 hybrids. Due to the breakdown of
    # co-adapted gene complexes during meiosis, their average fitness is often lower.
    fitness_f2_recombinant = 0.85

    print("--- Measuring Cost of Gene Flow in Yeast ---")
    print(f"Step 1: Measure fitness of the non-hybrid parental line.")
    print(f"Parental Fitness = {fitness_parental}\n")

    print(f"Step 2: Create F1 hybrids, then perform 'within mating' to produce F2 recombinants.")
    print(f"Measure the average fitness of the F2 recombinant population.")
    print(f"F2 Recombinant Fitness = {fitness_f2_recombinant}\n")

    # Step 2: Calculate the relative fitness (w) of the F2 recombinants.
    # This is the fitness of the F2s relative to the parental line.
    relative_fitness_w = fitness_f2_recombinant / fitness_parental
    
    print(f"Step 3: Calculate the relative fitness (w) of the F2 generation.")
    print(f"w = F2 Fitness / Parental Fitness")
    print(f"w = {fitness_f2_recombinant} / {fitness_parental} = {relative_fitness_w:.2f}\n")

    # Step 3: Calculate the selection coefficient (s), which represents the cost.
    selection_coefficient_s = 1 - relative_fitness_w

    print(f"Step 4: Calculate the selection coefficient (s) against the F2s.")
    print(f"This coefficient represents the 'cost of gene flow'.")
    print(f"s = 1 - w")
    # Final equation with numbers plugged in
    print(f"s = {1} - {relative_fitness_w:.2f} = {selection_coefficient_s:.2f}\n")

    print(f"The calculated cost due to gene flow (selection coefficient) is: {selection_coefficient_s:.2f}")
    print("This indicates a 15% fitness reduction in the F2 generation compared to the parental lines.")

# Run the calculation
calculate_gene_flow_cost()