import math

def calculate_cost_of_gene_flow():
    """
    This function simulates the calculation of the cost of gene flow (selection coefficient)
    in yeast, following the principles in the best answer.

    We will use hypothetical growth rates (e.g., doublings per hour).
    - A 'no gene flow' line serves as the control.
    - An 'F1 hybrid' line from the initial cross.
    - An 'F2 hybrid' line from mating F1 hybrids ('within mating') to see effects of meiosis.
    """
    # Step 1: Define hypothetical fitness data (growth rate)
    no_gene_flow_fitness = 0.55  # Fitness of the control line
    f1_hybrid_fitness = 0.54     # Fitness of F1 hybrids (often similar to parents)
    f2_hybrid_fitness = 0.45     # Fitness of F2 hybrids (often lower due to meiotic breakdown)

    print("--- Measuring Cost of Gene Flow in Yeast ---")
    print(f"Control (No Gene Flow) Fitness: {no_gene_flow_fitness}")
    print(f"F1 Hybrid Fitness: {f1_hybrid_fitness}")
    print(f"F2 Hybrid Fitness (post-meiosis): {f2_hybrid_fitness}\n")

    # Step 2: Calculate relative fitness (w) compared to the control line.
    # The control line has a relative fitness of 1.
    print("--- Calculating F1 Hybrid Cost ---")
    w_f1 = f1_hybrid_fitness / no_gene_flow_fitness
    print(f"1. Calculate Relative Fitness (w) of F1 Hybrid:")
    print(f"   w_f1 = F1 Hybrid Fitness / Control Fitness")
    print(f"   w_f1 = {f1_hybrid_fitness} / {no_gene_flow_fitness} = {w_f1}")

    # Step 3: Calculate the selection coefficient (s), which represents the cost. s = 1 - w
    s_f1 = 1 - w_f1
    print("\n2. Calculate Selection Coefficient (s) for F1 Hybrid:")
    print(f"   s_f1 = 1 - w_f1")
    print(f"   s_f1 = {1} - {w_f1} = {s_f1}")
    print(f"The cost of gene flow in the F1 generation is {s_f1:.4f}\n")


    # Repeat for the F2 generation, where the cost is often more pronounced
    print("--- Calculating F2 Hybrid Cost (reveals effects of meiosis) ---")
    w_f2 = f2_hybrid_fitness / no_gene_flow_fitness
    print(f"1. Calculate Relative Fitness (w) of F2 Hybrid:")
    print(f"   w_f2 = F2 Hybrid Fitness / Control Fitness")
    print(f"   w_f2 = {f2_hybrid_fitness} / {no_gene_flow_fitness} = {w_f2}")

    s_f2 = 1 - w_f2
    print("\n2. Calculate Selection Coefficient (s) for F2 Hybrid:")
    print(f"   s_f2 = 1 - w_f2")
    print(f"   s_f2 = {1} - {w_f2} = {s_f2}")
    print(f"The cost of gene flow in the F2 generation is {s_f2:.4f}")

# Execute the function
calculate_cost_of_gene_flow()