import numpy as np

def calculate_gene_flow_cost():
    """
    This function simulates the calculation of the cost of gene flow in yeast.

    It follows the principle of comparing the fitness of F2 hybrids (which have
    undergone gene flow and meiosis) to the fitness of the original parental lines
    (the 'no gene flow' control). The cost is quantified as a selection coefficient (s).
    """

    # Step 1: Define hypothetical fitness data (e.g., relative growth rates).
    # These represent measurements from a competition experiment.
    # A value of 1.0 represents the baseline fitness.
    parental_line_fitness = np.array([1.02, 0.99, 1.01, 1.03, 0.98, 1.00])
    
    # The F2 hybrids may show reduced fitness (outbreeding depression) after meiosis.
    hybrid_F2_fitness = np.array([0.95, 0.92, 0.99, 0.91, 0.96, 0.94])

    # Step 2: Calculate the mean fitness for each group.
    mean_parental_fitness = np.mean(parental_line_fitness)
    mean_hybrid_fitness = np.mean(hybrid_F2_fitness)

    print(f"Mean Parental Fitness (No Gene Flow Control): {mean_parental_fitness:.4f}")
    print(f"Mean F2 Hybrid Fitness (With Gene Flow + Meiosis): {mean_hybrid_fitness:.4f}\n")

    # Step 3: Calculate the relative fitness (w) of the hybrids.
    # This compares the hybrid fitness to the parental fitness.
    relative_fitness_w = mean_hybrid_fitness / mean_parental_fitness
    
    print("--- Calculating Relative Fitness (w) ---")
    print(f"w = Mean Hybrid Fitness / Mean Parental Fitness")
    print(f"w = {mean_hybrid_fitness:.4f} / {mean_parental_fitness:.4f}")
    print(f"w = {relative_fitness_w:.4f}\n")

    # Step 4: Calculate the selection coefficient (s).
    # s = w - 1. A negative value indicates a fitness cost.
    selection_coefficient_s = relative_fitness_w - 1

    print("--- Calculating Selection Coefficient (s) ---")
    print(f"s = Relative Fitness (w) - 1")
    print(f"s = {relative_fitness_w:.4f} - 1")
    print(f"s = {selection_coefficient_s:.4f}\n")

    print(f"Final Result:")
    print(f"The selection coefficient (s) against the hybrids is: {selection_coefficient_s:.4f}")
    print("A negative selection coefficient indicates a fitness cost due to gene flow.")

# Run the calculation
calculate_gene_flow_cost()