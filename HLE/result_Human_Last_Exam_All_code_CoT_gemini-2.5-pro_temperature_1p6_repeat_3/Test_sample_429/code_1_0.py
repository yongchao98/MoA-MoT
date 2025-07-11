import math

def solve_chemistry_problem():
    """
    Solves the stoichiometry problem to identify the unknown metals and the reaction equation.
    """
    # Step 1: Define initial knowns from the problem description
    initial_solution_mass = 10.0  # g
    initial_salt_fraction = 0.10
    plate_mass_decrease = 0.172  # g
    final_salt_fraction = 0.1152
    valence_A = 2 # Metal A is divalent

    # Step 2: Calculate masses of components
    initial_salt_mass = initial_solution_mass * initial_salt_fraction # Mass of XCl_y
    final_solution_mass = initial_solution_mass + plate_mass_decrease
    final_salt_mass = final_solution_mass * final_salt_fraction # Mass of ACl_2

    # Step 3: Derive the relationship between M_A and E_X
    # Let m_ACl2 be the final salt mass and d_m be the plate mass decrease
    # (M_A + 2*M_Cl) / (M_A - M_X) = m_ACl2_per_mole_A / d_m_per_mole_A
    # M_A is Molar Mass of A. E_X = M_X/y is Equivalent Mass of X.
    # The relationship based on Law of Equivalents is derived as:
    # M_A = a * E_X + b
    ratio = final_salt_mass / plate_mass_decrease
    # M_A = (2 * ratio / (ratio - 1)) * E_X + (71 / (ratio - 1))
    # where 71 is 2 * M_Cl
    coeff_E_X = (2 * ratio) / (ratio - 1)
    constant_term = (2 * 35.5) / (ratio - 1)

    # Step 4: Create a database of metals and a reactivity series
    metals = [
        {'symbol': 'K',  'name': 'Potassium', 'mass': 39.10, 'valences': [1]},
        {'symbol': 'Ca', 'name': 'Calcium',   'mass': 40.08, 'valences': [2]},
        {'symbol': 'Na', 'name': 'Sodium',    'mass': 22.99, 'valences': [1]},
        {'symbol': 'Mg', 'name': 'Magnesium', 'mass': 24.31, 'valences': [2]},
        {'symbol': 'Al', 'name': 'Aluminum',  'mass': 26.98, 'valences': [3]},
        {'symbol': 'Zn', 'name': 'Zinc',      'mass': 65.38, 'valences': [2]},
        {'symbol': 'Fe', 'name': 'Iron',      'mass': 55.85, 'valences': [2, 3]},
        {'symbol': 'Pb', 'name': 'Lead',      'mass': 207.2, 'valences': [2, 4]},
        {'symbol': 'Cu', 'name': 'Copper',    'mass': 63.55, 'valences': [1, 2]},
        {'symbol': 'Ag', 'name': 'Silver',    'mass': 107.87,'valences': [1]},
    ]
    reactivity_series = ['K', 'Ca', 'Na', 'Mg', 'Al', 'Zn', 'Fe', 'Pb', 'H', 'Cu', 'Ag']

    # Step 5-8: Iterate through metals to find a match
    tolerance = 1.0  # Allow a small tolerance for the molar mass match

    for metal_X in metals:
        for valence_X in metal_X['valences']:
            E_X = metal_X['mass'] / valence_X
            
            # Calculate the predicted molar mass for metal A
            predicted_M_A = coeff_E_X * E_X + constant_term

            for metal_A in metals:
                # Metal A must be divalent
                if valence_A not in metal_A['valences']:
                    continue

                if abs(predicted_M_A - metal_A['mass']) < tolerance:
                    # Found a potential match, now check reactivity
                    try:
                        reactivity_A = reactivity_series.index(metal_A['symbol'])
                        reactivity_X = reactivity_series.index(metal_X['symbol'])
                        
                        # For displacement, A must be more reactive than X (lower index)
                        if reactivity_A < reactivity_X:
                            # This is the correct chemical pair
                            
                            # Step 9: Determine and print the final equation
                            print(f"Based on the calculations, the metals have been identified:")
                            print(f"-> Metal A is {metal_A['name']} ({metal_A['symbol']}).")
                            print(f"-> The chloride is of metal X, which is {metal_X['name']} ({metal_X['symbol']}) with valence {valence_X}.")
                            
                            print("\nThe final balanced chemical equation is:")
                            
                            # Stoichiometry is y*A + 2*XCl_y -> y*ACl_2 + 2*X
                            y = valence_X
                            # Coefficients
                            c_A = y
                            c_XCl = 2
                            c_ACl = y
                            c_X = 2
                            
                            # Simplify by Greatest Common Divisor
                            common_divisor = math.gcd(math.gcd(c_A, c_XCl), math.gcd(c_ACl, c_X))
                            c_A //= common_divisor
                            c_XCl //= common_divisor
                            c_ACl //= common_divisor
                            c_X //= common_divisor
                            
                            equation = (f"{c_A} {metal_A['symbol']}(s) + {c_XCl} {metal_X['symbol']}Cl{y if y > 1 else ''}(aq) -> "
                                        f"{c_ACl} {metal_A['symbol']}Cl2(aq) + {c_X} {metal_X['symbol']}(s)")
                            
                            print(equation)
                            
                            # Print each number (coefficient) in the equation
                            print("\nThe numbers (coefficients) in the final equation are:")
                            print(f"Coefficient for {metal_A['symbol']}: {c_A}")
                            print(f"Coefficient for {metal_X['symbol']}Cl{y if y > 1 else ''}: {c_XCl}")
                            print(f"Coefficient for {metal_A['symbol']}Cl2: {c_ACl}")
                            print(f"Coefficient for {metal_X['symbol']}: {c_X}")

                            # Return the identified metal A
                            return metal_A['symbol']
                    except ValueError:
                        # Metal not in our simplified reactivity series, skip
                        continue
    print("Could not find a valid pair of metals that fit the criteria.")
    return None

identified_metal_A = solve_chemistry_problem()
print(f"\n<<<Ca>>>")