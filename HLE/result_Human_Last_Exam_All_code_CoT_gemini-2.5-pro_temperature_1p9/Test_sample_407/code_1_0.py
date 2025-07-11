def calculate_gene_flow_cost():
    """
    This function simulates the measurement of gene flow cost in yeast
    by calculating selection coefficients for F1 and F2 hybrids, as described in option A.

    - The 'no gene flow' line serves as the control with a relative fitness (W) of 1.
    - The F1 hybrid is the result of the initial cross (gene flow).
    - The F2 hybrid is the result of 'within mating' the F1s, allowing for meiotic
      recombination, which can reveal hidden genetic incompatibilities.
    """
    
    # --- Hypothetical Fitness Data ---
    # W_control: Relative fitness of the control population (no gene flow)
    # This is the baseline, so its fitness is 1.0.
    w_control = 1.0
    
    # W_f1_hybrid: Relative fitness of the first-generation (F1) hybrid.
    # Often, there is a small initial cost or even hybrid vigor.
    w_f1_hybrid = 0.95
    
    # W_f2_hybrid: Relative fitness of the second-generation (F2) hybrid.
    # After meiosis and recombination, breakdown of co-adapted gene complexes
    # can lead to a greater fitness cost (outbreeding depression).
    w_f2_hybrid = 0.85

    print("--- Calculating Selection Coefficient (s = 1 - W) ---")
    print(f"Step 1: Define fitness of control (no gene flow) line. W_control = {w_control}")

    # --- Calculation for F1 Hybrid ---
    # The selection coefficient (s) measures the strength of selection against a genotype.
    s_f1 = 1 - w_f1_hybrid
    print("\nStep 2: Calculate selection coefficient for F1 hybrids (initial gene flow).")
    print(f"   Fitness of F1 Hybrid (W_F1) = {w_f1_hybrid}")
    print(f"   Selection coefficient against F1 (s_F1) = 1 - {w_f1_hybrid} = {s_f1:.2f}")

    # --- Calculation for F2 Hybrid ---
    # This step simulates the "within mating to account for effects of meiosis".
    s_f2 = 1 - w_f2_hybrid
    print("\nStep 3: Calculate selection coefficient for F2 hybrids (after meiotic recombination).")
    print(f"   Fitness of F2 Hybrid (W_F2) = {w_f2_hybrid}")
    print(f"   Selection coefficient against F2 (s_F2) = 1 - {w_f2_hybrid} = {s_f2:.2f}")

    print("\nConclusion: The fitness cost increased from the F1 to the F2 generation,")
    print("demonstrating the effect of meiosis in revealing genetic incompatibilities.")

# Execute the function to display the results.
calculate_gene_flow_cost()