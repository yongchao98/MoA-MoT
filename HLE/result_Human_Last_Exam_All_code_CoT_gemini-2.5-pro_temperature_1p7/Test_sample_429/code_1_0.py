import math

def solve_chemistry_problem():
    """
    Solves the chemistry problem by searching for the best-fit metals.
    """
    # Step 1: Define a list of common divalent metals with their properties
    # Reactivity is simplified: higher number is more reactive.
    metals = {
        'Mg': {'mass': 24.31, 'reactivity': 8},
        'Al': {'mass': 26.98, 'reactivity': 7}, # Often trivalent, but for completeness
        'Zn': {'mass': 65.38, 'reactivity': 6},
        'Fe': {'mass': 55.85, 'reactivity': 5},
        'Cd': {'mass': 112.41, 'reactivity': 4},
        'Ni': {'mass': 58.69, 'reactivity': 3},
        'Pb': {'mass': 207.2, 'reactivity': 2},
        'Cu': {'mass': 63.55, 'reactivity': 1},
        'Sr': {'mass': 87.62, 'reactivity': 9} # A less common but possible candidate
    }

    # Problem data
    m_initial_MCl2 = 1.0  # g
    target_delta_plate = -0.172  # g
    target_w_final = 0.1152 # 11.52%

    best_fit = None
    min_error = float('inf')

    # Step 2: Iterate through all possible pairs for metal A and metal M
    for name_A, props_A in metals.items():
        for name_M, props_M in metals.items():
            if name_A == name_M:
                continue

            # Step 3: Apply physical and chemical constraints
            # Constraint 1: Metal A is more reactive than Metal M
            is_more_reactive = props_A['reactivity'] > props_M['reactivity']
            # Constraint 2: Mass of plate decreased, so Ar(A) > Ar(M)
            is_heavier = props_A['mass'] > props_M['mass']

            if is_more_reactive and is_heavier:
                # Step 4: For this valid pair, calculate the expected outcomes
                ar_A = props_A['mass']
                ar_M = props_M['mass']

                # Molar mass of MCl2
                mr_MCl2 = ar_M + 2 * 35.45
                
                # Moles of MCl2 reacted (from 1g initial)
                n = m_initial_MCl2 / mr_MCl2

                # Calculated change in plate mass
                calc_delta_plate = n * (ar_M - ar_A)

                # Calculated mass of ACl2 produced
                mr_ACl2 = ar_A + 2 * 35.45
                calc_m_ACl2 = n * mr_ACl2

                # Calculated final solution weight
                calc_sol_final = 10.0 - calc_delta_plate
                
                # Calculated final mass fraction
                calc_w_final = calc_m_ACl2 / calc_sol_final

                # Step 5: Calculate error compared to problem data
                error_delta = abs(calc_delta_plate - target_delta_plate) / abs(target_delta_plate)
                error_w = abs(calc_w_final - target_w_final) / target_w_final
                total_error = error_delta + error_w # Simple error sum

                if total_error < min_error:
                    min_error = total_error
                    best_fit = {
                        'A': name_A,
                        'M': name_M,
                        'Ar(A)': ar_A,
                        'Ar(M)': ar_M
                    }
    
    # Step 6: Print the determined metal and the reaction equation
    if best_fit:
        metal_A_name = best_fit['A']
        metal_M_name = best_fit['M']
        
        print(f"The analysis indicates that Metal A is {metal_A_name} and the unknown chloride was {metal_M_name}(II) Chloride.")
        print("The reaction equation is:")
        print(f"{metal_A_name} + {metal_M_name}Cl2 -> {metal_A_name}Cl2 + {metal_M_name}")
    else:
        print("No suitable pair of metals found.")

solve_chemistry_problem()
<<<Sr + CuCl2 -> SrCl2 + Cu>>>