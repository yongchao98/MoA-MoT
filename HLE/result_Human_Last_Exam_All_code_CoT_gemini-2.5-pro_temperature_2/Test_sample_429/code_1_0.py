import math

def solve_chemistry_problem():
    """
    Solves the stoichiometry problem to identify the unknown metal A and the reaction.
    """
    # 1. Define given values
    initial_solution_mass = 10.0  # in grams
    initial_salt_fraction = 0.10
    plate_mass_decrease = 0.172  # in grams
    # Atomic mass of Chlorine
    M_Cl = 35.5

    # 2. Calculate key masses based on problem statement and mass conservation
    # Initial mass of the salt in the solution (let's call it XCl_n)
    initial_salt_mass = initial_solution_mass * initial_salt_fraction

    # From mass balance, the mass of ACl2 formed must equal the mass of XCln reacted plus the change in plate mass.
    # mass(ACl2) - mass(XCln) = mass(A_dissolved) - mass(X_deposited) = plate_mass_decrease
    final_salt_mass = initial_salt_mass + plate_mass_decrease

    # 3. Establish a relationship between the atomic masses M_A and M_X
    # We assume the simplest and most common case: a divalent-divalent reaction (n=2).
    # Reaction: A + XCl2 -> ACl2 + X
    # Moles of reacted salt, x = initial_salt_mass / (M_X + 2 * M_Cl)
    # Moles of produced salt, x = final_salt_mass / (M_A + 2 * M_Cl)
    # Therefore: (M_A + 71) / (M_X + 71) = final_salt_mass / initial_salt_mass
    # This leads to the relationship: M_A = (final_salt_mass / initial_salt_mass) * (M_X + 71) - 71
    ratio = final_salt_mass / initial_salt_mass
    
    # 4. Search for plausible metals
    # We create a dictionary of common divalent metals {name: atomic_mass}.
    # A simplified reactivity series: K > Ca > Na > Mg > Al > Zn > Fe > Pb > H > Cu > Ag > Au
    divalent_metals = {
        'Ca': 40.08, 'Mg': 24.31, 'Zn': 65.38, 'Fe': 55.85, 
        'Pb': 207.2, 'Cu': 63.55, 'Sr': 87.62, 'Ba': 137.33
    }
    
    reactivity = ['Ba', 'Sr', 'Ca', 'Mg', 'Zn', 'Fe', 'Pb', 'Cu']

    best_match = {'A': None, 'X': None, 'error': float('inf')}

    for name_X, M_X in divalent_metals.items():
        # Calculate the theoretical atomic mass of metal A
        M_A_theoretical = ratio * (M_X + 71) - 71
        
        # Find the closest actual metal from our list
        for name_A, M_A_actual in divalent_metals.items():
            # The reaction A + XCl2 -> ... requires A to be more reactive than X
            if reactivity.index(name_A) < reactivity.index(name_X):
                error = abs(M_A_theoretical - M_A_actual)
                if error < best_match['error']:
                    best_match['error'] = error
                    best_match['A'] = name_A
                    best_match['X'] = name_X
    
    # 5. Output the result
    metal_A = best_match['A']
    metal_X = best_match['X']
    
    print(f"The analysis determined that the unknown metal (A) is {metal_A}.")
    print("The metal in the initial chloride solution (X) is most likely " + metal_X + ".")
    print("\nThe balanced chemical equation for the reaction is:")
    
    # In a simple displacement like this, all coefficients are 1.
    # The prompt requires printing each number in the equation.
    c1, c2, c3, c4 = 1, 1, 1, 1
    
    print(f"{c1} {metal_A} + {c2} {metal_X}Cl2 -> {c3} {metal_A}Cl2 + {c4} {metal_X}")


solve_chemistry_problem()
<<<Calcium>>>