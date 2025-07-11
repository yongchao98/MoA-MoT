import sys

def solve_chemistry_problem():
    """
    This function solves the stoichiometry problem to identify two metals and write the reaction equation.
    """
    # Step 1 & 2: Define variables and initial calculations
    m_solution_initial = 10.0  # g
    w_salt_initial = 0.10      # 10%
    plate_mass_decrease = 0.172 # g
    w_salt_final = 0.1152      # 11.52%
    Ar_Cl = 35.5               # Atomic mass of Chlorine

    # Mass of the initial unknown chloride salt (MCl2)
    m_MCl2 = m_solution_initial * w_salt_initial

    # The mass of the solution increases by the amount the plate mass decreases
    m_solution_final = m_solution_initial + plate_mass_decrease

    # Mass of the final salt (ACl2) formed
    m_ACl2 = m_solution_final * w_salt_final
    
    # Step 3 & 4: Establish the relationship between atomic masses
    # The relationship is derived from n = mass / molar_mass, where n is constant for both salts.
    # m_ACl2 / (Ar_A + 71) = m_MCl2 / (Ar_M + 71)
    # Ar_A = (m_ACl2 / m_MCl2) * (Ar_M + 71) - 71
    # Let's calculate the ratio m_ACl2 / m_MCl2
    mass_ratio = m_ACl2 / m_MCl2

    def calculate_Ar_A(Ar_M):
        """Calculates the theoretical atomic mass of A given the atomic mass of M."""
        return mass_ratio * (Ar_M + 2 * Ar_Cl) - (2 * Ar_Cl)

    # Step 5: Create a list of candidate metals with their atomic mass and a reactivity score.
    # Higher score means more reactive. A must have a higher score than M.
    # Reactivity series: K > Ca > Mg > Al > Zn > Fe > Pb > Cu
    # We will use common divalent metals.
    candidate_metals = {
        'Magnesium (Mg)': {'mass': 24.3, 'reactivity': 6},
        'Calcium (Ca)':   {'mass': 40.1, 'reactivity': 7},
        'Strontium (Sr)': {'mass': 87.6, 'reactivity': 8},
        'Barium (Ba)':    {'mass': 137.3,'reactivity': 9},
        'Iron (Fe)':      {'mass': 55.8, 'reactivity': 4},
        'Zinc (Zn)':      {'mass': 65.4, 'reactivity': 5},
        'Lead (Pb)':      {'mass': 207.2,'reactivity': 3},
        'Copper (Cu)':    {'mass': 63.5, 'reactivity': 2},
        'Nickel (Ni)':    {'mass': 58.7, 'reactivity': 4},
        'Manganese (Mn)': {'mass': 54.9, 'reactivity': 5}
    }

    best_fit = {
        'metal_A': None,
        'metal_M': None,
        'min_diff': sys.float_info.max
    }

    # Iterate through all possibilities for metal M
    for m_name, m_props in candidate_metals.items():
        Ar_M_actual = m_props['mass']
        reactivity_M = m_props['reactivity']
        
        # Calculate the theoretical atomic mass for A
        Ar_A_theoretical = calculate_Ar_A(Ar_M_actual)

        # Find the best match for A from the candidate list
        for a_name, a_props in candidate_metals.items():
            Ar_A_actual = a_props['mass']
            reactivity_A = a_props['reactivity']
            
            # Check for two conditions:
            # 1. Metal A must be more reactive than Metal M.
            # 2. The pair provides the closest match found so far.
            if reactivity_A > reactivity_M:
                diff = abs(Ar_A_theoretical - Ar_A_actual)
                if diff < best_fit['min_diff']:
                    best_fit['min_diff'] = diff
                    best_fit['metal_A'] = a_name
                    best_fit['metal_M'] = m_name

    # Step 6: Output the results
    metal_A_name = best_fit['metal_A']
    metal_M_name = best_fit['metal_M']
    
    # Extract symbols for the equation
    metal_A_symbol = metal_A_name.split('(')[1].replace(')', '')
    metal_M_symbol = metal_M_name.split('(')[1].replace(')', '')
    
    print(f"The analysis indicates that metal A is {metal_A_name}.")
    print(f"The unknown chloride was a salt of {metal_M_name}.")
    print("\nThe equation for the reaction is:")
    
    # The stoichiometry is 1:1 for this reaction type
    coeff_A = 1
    coeff_MCl2 = 1
    coeff_ACl2 = 1
    coeff_M = 1
    
    print(f"{coeff_A} {metal_A_symbol} + {coeff_MCl2} {metal_M_symbol}Cl2 -> {coeff_ACl2} {metal_A_symbol}Cl2 + {coeff_M} {metal_M_symbol}")

solve_chemistry_problem()
<<<Metal A: Strontium (Sr), Unknown metal: Copper (Cu), Equation: 1 Sr + 1 CuCl2 -> 1 SrCl2 + 1 Cu>>>