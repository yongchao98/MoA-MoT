def analyze_simpsons_diversity():
    """
    Demonstrates the calculation of Simpson's Diversity Index for different scenarios
    and explains the mathematical and ecological validity of a D=0 result.
    """

    print("--- Analysis of Simpson's Diversity Index (D) ---")
    print("Formula: D = 1 - (Σn(n-1) / N(N-1))\n")

    # --- Case 1: The student's actual finding (hypothetical numbers) ---
    # This scenario leads to a diversity index of 0.
    # Let's assume the student found 85 bats, all of the same species.
    one_species_counts = [85]
    n_case1 = one_species_counts[0]
    N_case1 = sum(one_species_counts)

    print("--- Scenario 1: Mathematically Valid Case (D=0) ---")
    print("This represents a sample where only ONE species was found.")
    print(f"Species counts (n): {one_species_counts}")
    print(f"Total individuals (N): {N_case1}")

    # Calculate terms for the equation
    sum_n_n_minus_1_case1 = n_case1 * (n_case1 - 1)
    N_N_minus_1_case1 = N_case1 * (N_case1 - 1)
    
    # Calculate D
    if N_N_minus_1_case1 == 0:
        D_case1 = 0.0
    else:
        D_case1 = 1 - (sum_n_n_minus_1_case1 / N_N_minus_1_case1)

    print("\nCalculating the equation parts:")
    print(f"Σn(n-1) = {n_case1}*({n_case1}-1) = {sum_n_n_minus_1_case1}")
    print(f"N(N-1)   = {N_case1}*({N_case1}-1) = {N_N_minus_1_case1}")
    
    print("\nFinal Equation:")
    print(f"D = 1 - ({sum_n_n_minus_1_case1} / {N_N_minus_1_case1})")
    print(f"D = 1 - {sum_n_n_minus_1_case1 / N_N_minus_1_case1 if N_N_minus_1_case1 != 0 else 'undefined'}")
    print(f"Resulting Simpson's Index (D): {D_case1:.2f}\n")
    print("Conclusion: A value of D=0 is mathematically valid and occurs when the sample has no diversity.")

    # --- Case 2: A more ecologically realistic finding for this island ---
    # Let's assume a sample with three species.
    multiple_species_counts = [40, 30, 15]
    N_case2 = sum(multiple_species_counts)

    print("\n\n--- Scenario 2: Ecologically Realistic Case (D>0) ---")
    print("This represents a sample where MULTIPLE species were found.")
    print(f"Species counts (n): {multiple_species_counts}")
    print(f"Total individuals (N): {N_case2}")
    
    # Calculate terms for the equation
    sum_n_n_minus_1_case2 = sum([n * (n - 1) for n in multiple_species_counts])
    N_N_minus_1_case2 = N_case2 * (N_case2 - 1)

    # Calculate D
    D_case2 = 1 - (sum_n_n_minus_1_case2 / N_N_minus_1_case2)
    
    print("\nCalculating the equation parts:")
    print(f"Σn(n-1) = (40*39) + (30*29) + (15*14) = 1560 + 870 + 210 = {sum_n_n_minus_1_case2}")
    print(f"N(N-1)   = {N_case2}*({N_case2}-1) = {N_N_minus_1_case2}")

    print("\nFinal Equation:")
    print(f"D = 1 - ({sum_n_n_minus_1_case2} / {N_N_minus_1_case2})")
    print(f"D = 1 - {(sum_n_n_minus_1_case2 / N_N_minus_1_case2):.4f}")
    print(f"Resulting Simpson's Index (D): {D_case2:.4f}\n")
    print("Conclusion: A value of D>0 indicates the presence of diversity in the sample.")

# Execute the analysis
analyze_simpsons_diversity()