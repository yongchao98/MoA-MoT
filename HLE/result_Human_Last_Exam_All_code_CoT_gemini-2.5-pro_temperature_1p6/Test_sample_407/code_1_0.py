def calculate_gene_flow_cost():
    """
    This function demonstrates how to calculate the cost of gene flow by measuring
    the selection coefficient (s) in F1 and F2 hybrid generations compared to a
    parental (no gene flow) line.

    The cost is calculated based on a fitness proxy, such as growth rate.
    - A positive selection coefficient (s > 0) indicates a fitness cost (outbreeding depression).
    - A negative selection coefficient (s < 0) indicates a fitness advantage (hybrid vigor).

    This example shows a scenario with slight outbreeding depression in the F1
    generation and more severe depression (hybrid breakdown) in the F2 generation,
    which can occur after meiosis and recombination.
    """

    # --- Step 1: Define hypothetical fitness data (e.g., growth rate) ---
    # Fitness of the original parental line (no gene flow)
    parental_fitness = 0.25  # e.g., growth rate in divisions per hour

    # Fitness of the F1 hybrid (direct result of gene flow)
    f1_hybrid_fitness = 0.24  # Slightly lower fitness than parent

    # Fitness of the F2 hybrid (after F1s mate, showing effects of meiosis)
    f2_hybrid_fitness = 0.20  # Significantly lower fitness due to hybrid breakdown

    print("--- Measuring Cost of Gene Flow in Yeast ---")
    print(f"Parental (No Gene Flow) Fitness: {parental_fitness}")
    print(f"F1 Hybrid Fitness: {f1_hybrid_fitness}")
    print(f"F2 Hybrid Fitness (Post-Meiosis): {f2_hybrid_fitness}\n")

    # --- Step 2: Calculate selection coefficient for F1 hybrids ---
    # Relative fitness (w) = Fitness of hybrid / Fitness of parent
    # Selection coefficient (s) = 1 - w
    relative_fitness_f1 = f1_hybrid_fitness / parental_fitness
    selection_coefficient_f1 = 1 - relative_fitness_f1

    print("--- F1 Hybrid Generation Analysis ---")
    print("This measures the immediate effect of hybridization.")
    print(f"Calculation: s = 1 - (F1_Fitness / Parental_Fitness)")
    print(f"Equation: s_F1 = 1 - ({f1_hybrid_fitness} / {parental_fitness})")
    print(f"Result: Selection coefficient (s) for F1 hybrids = {selection_coefficient_f1:.4f}\n")


    # --- Step 3: Calculate selection coefficient for F2 hybrids to check for effects of meiosis ---
    relative_fitness_f2 = f2_hybrid_fitness / parental_fitness
    selection_coefficient_f2 = 1 - relative_fitness_f2

    print("--- F2 Hybrid Generation Analysis (Post 'Within Mating') ---")
    print("This measures for hybrid breakdown, accounting for the effects of meiosis.")
    print(f"Calculation: s = 1 - (F2_Fitness / Parental_Fitness)")
    print(f"Equation: s_F2 = 1 - ({f2_hybrid_fitness} / {parental_fitness})")
    print(f"Result: Selection coefficient (s) for F2 hybrids = {selection_coefficient_f2:.4f}\n")
    
    if selection_coefficient_f2 > selection_coefficient_f1:
        print("Conclusion: The fitness cost increased in the F2 generation, demonstrating hybrid breakdown.")

# Run the calculation
calculate_gene_flow_cost()