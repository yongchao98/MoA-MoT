def calculate_gene_flow_cost():
    """
    This function simulates the calculation of the cost of gene flow in yeast
    as described in the most comprehensive answer.

    The cost is quantified by the selection coefficient (s), which measures the
    relative fitness disadvantage of hybrids compared to a parental line.
    
    s = 1 - (W_hybrid / W_parent)
    
    where W is fitness (e.g., growth rate). A positive 's' indicates a cost.
    """
    
    # --- Hypothetical Fitness Data (e.g., growth rate in doublings/hour) ---
    # Fitness of the non-hybridized parental line (the "no gene flow" control)
    fitness_parent = 0.50
    
    # Fitness of the first-generation (F1) hybrid
    fitness_hybrid_f1 = 0.48
    
    # Fitness of the second-generation (F2) hybrid, after meiosis and recombination.
    # The cost is often greater in the F2 generation due to outbreeding depression.
    fitness_hybrid_f2 = 0.41

    print("Measuring the cost of gene flow using the selection coefficient (s).\n")
    print(f"Parental (no gene flow) fitness: {fitness_parent}")
    print(f"F1 Hybrid (gene flow) fitness: {fitness_hybrid_f1}")
    print(f"F2 Hybrid (post-meiosis) fitness: {fitness_hybrid_f2}\n")

    # --- Calculation for the F1 Generation ---
    # This represents the initial cost of hybridization.
    s_f1 = 1 - (fitness_hybrid_f1 / fitness_parent)
    
    print("1. Calculate the selection coefficient (s) for the F1 hybrid generation.")
    print(f"Equation: s_F1 = 1 - (fitness_hybrid_f1 / fitness_parent)")
    print(f"Calculation: s_F1 = 1 - ({fitness_hybrid_f1} / {fitness_parent})")
    print(f"Result: The cost in the F1 generation is s_F1 = {s_f1:.4f}\n")
    
    # --- Calculation for the F2 Generation ---
    # This calculation after "within mating" accounts for the effects of meiosis.
    s_f2 = 1 - (fitness_hybrid_f2 / fitness_parent)

    print("2. Calculate the selection coefficient (s) for the F2 hybrid generation (post-meiosis).")
    print(f"Equation: s_F2 = 1 - (fitness_hybrid_f2 / fitness_parent)")
    print(f"Calculation: s_F2 = 1 - ({fitness_hybrid_f2} / {fitness_parent})")
    print(f"Result: The cost in the F2 generation is s_F2 = {s_f2:.4f}\n")
    
    print("Conclusion: The fitness cost is higher in the F2 generation,")
    print("demonstrating the importance of checking for effects after meiosis.")

# Run the calculation
calculate_gene_flow_cost()
