def solve_electron_energy_problem():
    """
    This program derives the minimum energy for an electron to cause an
    inter-band transition based on energy and momentum conservation.
    """
    print("Derivation for the Minimum Energy of Electron 1\n")

    # --- Define variables for clarity in print statements ---
    E_g = "E_g"
    C = "C" # Represents ħ²/(2m*)

    # --- Step 1: State Conservation Laws ---
    print("Step 1: Conservation Laws")
    print("-------------------------")
    print("Energy: E_1 + E_2 = E'_1 + E'_2")
    print("Momentum: k_1 + k_2 = k'_1 + k'_2\n")

    # --- Step 2: Substitute Energy Expressions ---
    print("Step 2: Substitute Energy Expressions into Energy Conservation")
    print("---------------------------------------------------------")
    print(f"Using E_1=E_g+C*k_1², E_2=-C*k_2², E'_1=E_g+C*(k'_1)², E'_2=E_g+C*(k'_2)²:")
    print(f"E_g + C*k_1² - C*k_2² = (E_g + C*(k'_1)²) + (E_g + C*(k'_2)²)")
    print(f"After simplifying: C*k_1² = {E_g} + C*k_2² + C*((k'_1)² + (k'_2)²)\n")

    # --- Step 3: Apply Threshold Conditions ---
    print("Step 3: Apply Threshold (Minimum Energy) Conditions")
    print("----------------------------------------------------")
    print("To minimize E_1 (and thus k_1), we must minimize the other terms.")
    print("We choose the most favorable initial state for electron 2, which is k_2 = 0.")
    print("The conservation laws become:")
    print(f"  Energy:   C*k_1² = {E_g} + C*((k'_1)² + (k'_2)²)  (Equation A)")
    print(f"  Momentum: k_1 = k'_1 + k'_2                     (Equation B)\n")

    # --- Step 4: Solve for the Final State ---
    print("Step 4: Combine the Equations")
    print("-----------------------------")
    print("Substitute momentum (B) into energy (A):")
    print(f"C*(k'_1 + k'_2)² = {E_g} + C*((k'_1)² + (k'_2)²) ")
    print("Expanding and simplifying yields a condition on the final state:")
    print(f"  2*C*k'_1.k'_2 = {E_g}\n")

    # --- Step 5: Express E_1 and Minimize ---
    print("Step 5: Express E_1 in terms of the final state and find the minimum")
    print("----------------------------------------------------------------------")
    print(f"From E_1 = E_g + C*k_1² and Equation A:")
    print(f"E_1 = E_g + ({E_g} + C*((k'_1)² + (k'_2)²)) = 2*{E_g} + C*((k'_1)² + (k'_2)²)")
    print("To minimize E_1, we must minimize the final kinetic energy C*((k'_1)² + (k'_2)²),")
    print(f"subject to the constraint 2*C*k'_1.k'_2 = {E_g}.")
    print("This minimum is achieved when k'_1 = k'_2.")
    print(f"In this case, the minimum final kinetic energy is exactly {E_g}.\n")

    # --- Step 6: Final Result ---
    print("Step 6: Calculate the Final Result")
    print("----------------------------------")
    print(f"Substituting the minimum final kinetic energy back into the expression for E_1:")
    print(f"E_1_min = 2*{E_g} + {E_g}")
    
    coefficient = 3
    print("\nThe final equation for the minimum energy of electron 1 is:")
    print(f"E_1_min = {coefficient} * {E_g}")

# Execute the derivation
solve_electron_energy_problem()