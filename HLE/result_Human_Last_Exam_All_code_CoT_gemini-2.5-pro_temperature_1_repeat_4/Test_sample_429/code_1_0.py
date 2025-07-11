import math

def solve_chemistry_riddle():
    """
    This script solves the chemistry problem by systematically searching for the
    unknown metals based on the data provided.
    """
    # 1. Define constants from the problem statement.
    m_solution_initial = 10.0  # g
    w_salt_initial = 0.10      # 10%
    delta_m_plate = 0.172      # g
    w_salt_final = 0.1152      # 11.52%
    AR_CL = 35.5               # Atomic mass of Chlorine

    # 2. Calculate key mass values based on the provided data.
    # The mass of the solution increases by the amount the plate mass decreases.
    m_solution_final = m_solution_initial + delta_m_plate

    # Mass of the initial salt (MCl_n) and final salt (ACl_2)
    m_salt_initial = m_solution_initial * w_salt_initial
    m_salt_final = m_solution_final * w_salt_final

    # 3. Check for data consistency and select the most reliable values.
    # The difference in the mass of the salts is due to the difference in the mass of the metals.
    # delta_m_salts = m_metal_A - m_metal_M
    # This value should be equal to the change in the plate's mass.
    delta_m_salts = m_salt_final - m_salt_initial
    
    # print(f"Note: Plate mass change given is {delta_m_plate} g.")
    # print(f"Note: Mass change calculated from salts is {delta_m_salts:.5f} g.")
    # There is a small inconsistency. We will use the value derived from the salt masses
    # as it makes the underlying stoichiometric equations consistent.
    
    # Let m_A be the mass of metal A that reacted, and m_M be the mass of metal M deposited.
    # From conservation of chlorine:
    # m_salt_final = m_A + m_chlorine
    # m_salt_initial = m_M + m_chlorine
    # Subtracting gives: m_salt_final - m_salt_initial = m_A - m_M.
    # So, m_A - m_M = delta_m_salts
    
    m_A_minus_m_M = delta_m_salts

    # 4. Define a database of common metals with their properties.
    # 'val' is a list of common valences (oxidation states).
    metals = {
        "K": {"Ar": 39.1, "val": [1]}, "Na": {"Ar": 23.0, "val": [1]},
        "Ca": {"Ar": 40.1, "val": [2]}, "Mg": {"Ar": 24.3, "val": [2]},
        "Al": {"Ar": 27.0, "val": [3]}, "Zn": {"Ar": 65.4, "val": [2]},
        "Fe": {"Ar": 55.8, "val": [2, 3]}, "Pb": {"Ar": 207.2, "val": [2]},
        "Cu": {"Ar": 63.5, "val": [1, 2]}, "Ag": {"Ar": 107.9, "val": [1]},
        "Ni": {"Ar": 58.7, "val": [2]}, "Mn": {"Ar": 54.9, "val": [2, 4, 7]},
        "Sr": {"Ar": 87.6, "val": [2]}, "Ba": {"Ar": 137.3, "val": [2]},
        "Cd": {"Ar": 112.4, "val": [2]}, "Ga": {"Ar": 69.7, "val": [3]},
    }
    # Simplified reactivity series to check if the reaction is possible.
    reactivity_series = ["K", "Ca", "Na", "Mg", "Al", "Mn", "Zn", "Fe", "Ni", "Pb", "Cu", "Ag"]

    # 5. Iterate through possibilities to find the solution.
    # Loop through all potential metals for A (must be divalent).
    for A_name, A_props in metals.items():
        if 2 not in A_props['val']:
            continue
        Ar_A = A_props['Ar']

        # Loop through possible valences (n) for metal M.
        for n in [1, 2, 3]:
            # From stoichiometric relationships, we can derive a formula for the atomic mass of M (Ar_M)
            # based on the atomic mass of A (Ar_A) and the valence of M (n).
            # Derived from: 2*(m_A/Ar_A) = n*(m_M/Ar_M) and relations for m_A, m_M.
            # Ar_M = 35.5 * n * (m_A - delta_m_salts) / (1 - (m_A - delta_m_salts))
            # Ar_A = 71 * m_A / (m_salt_final - m_A)
            # After substitution and simplification:
            numerator = AR_CL * n * ((m_salt_final - m_A_minus_m_M) * Ar_A - 2 * AR_CL * m_A_minus_m_M)
            denominator = 2 * AR_CL * m_salt_final
            Ar_M_target = numerator / denominator
            
            # Loop through all potential metals for M.
            for M_name, M_props in metals.items():
                if n not in M_props['val']:
                    continue

                # Check if the calculated atomic mass matches a known metal within a tolerance.
                if abs(M_props['Ar'] - Ar_M_target) < 1.5:
                    # Check if A is more reactive than M.
                    try:
                        if reactivity_series.index(A_name) < reactivity_series.index(M_name):
                            # Solution found. Print the results.
                            print(f"Analysis complete. The metals have been identified.\n")
                            print(f"Metal A (divalent) is: {A_name} (Atomic mass: {A_props['Ar']})")
                            print(f"The unknown chloride contains metal M: {M_name} (Atomic mass: {M_props['Ar']}, Valence: {n})\n")
                            
                            # Determine coefficients for the balanced equation: n A + 2 MCl_n -> n ACl_2 + 2 M
                            common_divisor = math.gcd(n, 2)
                            coeff_A = n // common_divisor
                            coeff_M = 2 // common_divisor
                            
                            print("The balanced chemical equation for the reaction is:")
                            print(f"{coeff_A} {A_name} + {coeff_M} {M_name}Cl_{n} -> {coeff_A} {A_name}Cl_2 + {coeff_M} {M_name}")
                            
                            global final_answer
                            final_answer = A_name
                            return

solve_chemistry_riddle()
<<<Calcium>>>