def solve_formula():
    """
    This function derives and prints the formula for the second voltage plateau
    of a graphite anode based on a physical approximation.
    """

    # Define the numbers that appear in the formula as variables
    stage_index_1 = 1
    stage_index_3 = 3
    denominator_coefficient = 2

    # Print the explanation for the formula
    print("The simple formula that best approximates the second plateau voltage (V_2) is based on the physical approximation that the chemical potential of stage 2 (μ_2) is the average of the chemical potentials of the adjacent stages, stage 1 (μ_1) and stage 3 (μ_3).")
    print("\nThe relationship is: μ_2 ≈ (μ_1 + μ_3) / 2")
    print("\nSince the voltage for a given stage k is V_k = -μ_k / e (where μ_k is the chemical potential relative to the Li/Li+ reference), we can substitute the approximation to find the formula for V_2.")

    # Construct and print the final formula using an f-string
    final_formula = f"-(μ_{stage_index_1} + μ_{stage_index_3}) / ({denominator_coefficient}*e)"
    print("\nThe final formula is:")
    print(f"V_2 ≈ {final_formula}")

    # Explicitly output each number in the final equation as requested
    print("\nIn this formula:")
    print(f"- The number '{stage_index_1}' is the index for the chemical potential of stage 1 (μ_1).")
    print(f"- The number '{stage_index_3}' is the index for the chemical potential of stage 3 (μ_3).")
    print(f"- The number '{denominator_coefficient}' is the coefficient in the denominator.")

solve_formula()