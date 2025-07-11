import math

def calculate_selection_coefficient():
    """
    This function demonstrates the calculation of the cost of gene flow
    by simulating fitness data and calculating the selection coefficient (s)
    as described in the correct answer choice.

    The cost of gene flow (outbreeding depression) is measured against a control
    population that did not experience gene flow. The selection coefficient, s,
    quantifies this cost.
    s = 1 - w, where w is the relative fitness of the hybrid line.
    """

    # --- Step 1: Define hypothetical fitness data ---
    # Fitness is measured (e.g., via growth rate) for different lines.
    # The parental line is the reference, so its fitness is normalized to 1.0.
    fitness_parental_line = 1.0  # Control (no gene flow)
    # Fitness of the F1 hybrids (first generation after gene flow)
    fitness_f1_hybrid = 0.95
    # Fitness of the F2 hybrids (after F1s mate, "within mating to account for effects of meiosis")
    # This value is often lower due to hybrid breakdown.
    fitness_f2_hybrid = 0.88

    print("Experimental Data (Hypothetical):")
    print(f"Fitness of No Gene Flow Line (Control): {fitness_parental_line}")
    print(f"Fitness of F1 Hybrid Line: {fitness_f1_hybrid}")
    print(f"Fitness of F2 Hybrid Line (post-meiosis): {fitness_f2_hybrid}\n")

    # --- Step 2: Calculate Relative Fitness (w) ---
    # Relative fitness w = fitness_of_line / fitness_of_reference
    w_f1 = fitness_f1_hybrid / fitness_parental_line
    w_f2 = fitness_f2_hybrid / fitness_parental_line

    # --- Step 3: Calculate Selection Coefficient (s) for each hybrid generation ---
    # s = 1 - w
    s_f1 = 1 - w_f1
    s_f2 = 1 - w_f2

    print("--- Calculation of Cost (Selection Coefficient, s) ---")
    print("\n1. Cost in F1 Generation:")
    print(f"   The relative fitness (w) of F1 hybrids is {fitness_f1_hybrid} / {fitness_parental_line} = {w_f1:.2f}")
    print(f"   The selection coefficient (s) = 1 - {w_f1:.2f} = {s_f1:.2f}")
    print(f"   This represents a cost of {s_f1:.2%} due to initial hybridization.")

    print("\n2. Cost in F2 Generation (after 'within mating'):")
    print(f"   The relative fitness (w) of F2 hybrids is {fitness_f2_hybrid} / {fitness_parental_line} = {w_f2:.2f}")
    print(f"   The selection coefficient (s) = 1 - {w_f2:.2f} = {s_f2:.2f}")
    print(f"   This represents a cost of {s_f2:.2%} due to hybrid breakdown from meiosis.")
    print("\nBy measuring both F1 and F2 fitness, we capture the full cost of gene flow.")

calculate_selection_coefficient()