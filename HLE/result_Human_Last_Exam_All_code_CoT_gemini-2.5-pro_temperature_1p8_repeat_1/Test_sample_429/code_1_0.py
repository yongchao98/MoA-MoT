import math

def solve_chemistry_problem():
    """
    This function solves the chemistry problem by finding the identity of metal A
    and writing the reaction equation.
    """

    # Dictionary of common divalent metals with their molar masses and a reactivity score (lower is more reactive)
    # A must be more reactive than B, so reactivity_A < reactivity_B.
    metals = {
        'Sr': {'mass': 87.62, 'reactivity': 2},
        'Ca': {'mass': 40.08, 'reactivity': 3},
        'Mg': {'mass': 24.31, 'reactivity': 4},
        'Zn': {'mass': 65.38, 'reactivity': 5},
        'Fe': {'mass': 55.85, 'reactivity': 6},
        'Ni': {'mass': 58.69, 'reactivity': 7},
        'Pb': {'mass': 207.2, 'reactivity': 8},
        'Cu': {'mass': 63.55, 'reactivity': 9}
    }

    best_match = {
        'metal_A': None,
        'metal_B': None,
        'min_diff': float('inf')
    }

    # Derived relationship: M_A = 1.172 * M_B + 12.212
    # We will iterate through possible metals for B and find the best matching A.
    for b_name, b_props in metals.items():
        m_b = b_props['mass']
        # Calculate the theoretical molar mass for metal A
        m_a_theoretical = 1.172 * m_b + 12.212

        # Find the best real metal A that matches the theoretical mass and reactivity
        for a_name, a_props in metals.items():
            if a_name == b_name:
                continue

            # Check reactivity condition: A must be more reactive than B
            if a_props['reactivity'] < b_props['reactivity']:
                m_a_real = a_props['mass']
                diff = abs(m_a_real - m_a_theoretical)

                if diff < best_match['min_diff']:
                    best_match['min_diff'] = diff
                    best_match['metal_A'] = a_name
                    best_match['metal_B'] = b_name

    metal_A_name = best_match['metal_A']
    metal_B_name = best_match['metal_B']

    if metal_A_name and metal_B_name:
        # Recalculate and print final values to verify the choice
        m_a = metals[metal_A_name]['mass']
        m_b = metals[metal_B_name]['mass']
        
        # Moles of BCl2 consumed
        moles_reacted = 1.0 / (m_b + 71)
        
        # Verify plate mass decrease
        plate_mass_decrease = moles_reacted * (m_a - m_b)
        
        # Verify final concentration
        mass_ACl2_formed = moles_reacted * (m_a + 71)
        final_solution_mass = 10 + plate_mass_decrease
        final_concentration = (mass_ACl2_formed / final_solution_mass) * 100

        print(f"The program analyzed common divalent metals based on the problem's data.")
        print(f"The best fit is found with Metal A = {metal_A_name} and Metal B = {metal_B_name}.")
        print("-" * 30)
        print("Verification with this pair:")
        print(f"  Calculated plate mass decrease: {plate_mass_decrease:.3f} g (Problem value: 0.172 g)")
        print(f"  Calculated final concentration: {final_concentration:.2f}% (Problem value: 11.52%)")
        print("-" * 30)
        print(f"The determined metal (A) is: {metal_A_name}")
        print("\nThe equation for the reaction is:")
        
        # Print the final equation with each part separated
        print(f"{metal_A_name} + {metal_B_name}Cl₂ → {metal_A_name}Cl₂ + {metal_B_name}")

    else:
        print("Could not determine the metals based on the provided data and list of common metals.")

solve_chemistry_problem()
<<<Sr>>>