import math

def calculate_gene_flow_cost():
    """
    This function demonstrates how to calculate the cost of gene flow in yeast,
    following the method described in the correct answer.

    The cost is quantified by the selection coefficient (s).
    """
    
    # Step 1: Measure fitness for both control and hybrid lines.
    # Fitness is often measured as the maximum growth rate in a competitive assay.
    # Let's assume we have measured these values from an experiment.

    # Fitness of the control line (parental strain, no gene flow).
    # We can normalize this to 1.0 for simplicity.
    fitness_no_gene_flow = 1.0 

    # Fitness of the hybrid line after it has undergone meiosis (e.g., F2 generation).
    # A lower value indicates a fitness cost due to gene flow.
    # Let's assume the hybrid line has a 5% lower fitness.
    fitness_hybrid_after_meiosis = 0.95
    
    print(f"Step 1: Define measured fitness values.")
    print(f"  - Fitness of No Gene Flow Line (W_control): {fitness_no_gene_flow}")
    print(f"  - Fitness of Hybrid Line after Meiosis (W_hybrid): {fitness_hybrid_after_meiosis}\n")
    
    # Step 2: Calculate the relative fitness of the hybrid line.
    # Relative fitness (W) is the fitness of the hybrid line in comparison to the control.
    relative_fitness = fitness_hybrid_after_meiosis / fitness_no_gene_flow
    
    print(f"Step 2: Calculate the relative fitness (W) of the hybrid.")
    print(f"  - Formula: W = W_hybrid / W_control")
    print(f"  - Calculation: W = {fitness_hybrid_after_meiosis} / {fitness_no_gene_flow} = {relative_fitness}\n")

    # Step 3: Calculate the selection coefficient (s), which represents the cost.
    # s = 1 - W
    selection_coefficient = 1 - relative_fitness

    print(f"Step 3: Calculate the selection coefficient (s).")
    print(f"  - This coefficient represents the 'cost' of gene flow.")
    print(f"  - Formula: s = 1 - W\n")

    # Final Output: Print the final equation with the numbers plugged in.
    print("--- Final Result ---")
    print("The final equation for the selection coefficient (s) is:")
    print(f"s = 1 - W_relative")
    # As requested, output the final equation showing each number.
    print(f"Final Equation: s = 1 - ({fitness_hybrid_after_meiosis} / {fitness_no_gene_flow})")
    print(f"Calculated cost of gene flow (s): {selection_coefficient:.2f}")
    if selection_coefficient > 0:
        print(f"This positive value indicates a fitness cost of {selection_coefficient:.0%}.")
    else:
        print("This value does not indicate a fitness cost.")

# Run the calculation
calculate_gene_flow_cost()