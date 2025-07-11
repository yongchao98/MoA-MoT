import numpy as np

def solve_chemistry_problem():
    """
    Solves the chemistry problem by testing common divalent metals for substance A.
    """
    # Given constants
    initial_solution_mass = 10.0  # g
    initial_salt_fraction = 0.10  # 10%
    final_salt_fraction = 0.1152 # 11.52%
    plate_mass_decrease = 0.172 # g
    E_Cl = 35.5

    # Calculated values
    initial_salt_mass = initial_solution_mass * initial_salt_fraction
    # The final solution mass must equal the initial mass plus the mass lost by the plate
    final_solution_mass = initial_solution_mass + plate_mass_decrease
    final_salt_mass = final_solution_mass * final_salt_fraction
    
    # Let E_A be the equivalent mass of metal A, and E_X be the equivalent mass of metal X.
    # We have two independent relationships to find E_X if E_A is known.
    
    # Equation 1 from law of equivalents on the salts:
    # Eq = initial_salt_mass / (E_X + E_Cl) = final_salt_mass / (E_A + E_Cl)
    # E_A + E_Cl = (final_salt_mass / initial_salt_mass) * (E_X + E_Cl)
    def get_Ex_from_salts(E_A):
        # E_A + E_Cl = (final_salt_mass / initial_salt_mass) * (E_X + E_Cl)
        # (E_A + E_Cl) * (initial_salt_mass / final_salt_mass) = E_X + E_Cl
        Ex = (E_A + E_Cl) * (initial_salt_mass / final_salt_mass) - E_Cl
        return Ex

    # Equation 2 from the plate mass change:
    # plate_mass_decrease = Eq * (E_A - E_X)
    # Substitute Eq = final_salt_mass / (E_A + E_Cl)
    # plate_mass_decrease = (final_salt_mass / (E_A + E_Cl)) * (E_A - E_X)
    def get_Ex_from_plate(E_A):
        # (plate_mass_decrease / final_salt_mass) * (E_A + E_Cl) = E_A - E_X
        # E_X = E_A - (plate_mass_decrease / final_salt_mass) * (E_A + E_Cl)
        Ex = E_A - (plate_mass_decrease / final_salt_mass) * (E_A + E_Cl)
        return Ex

    # Dictionary of common divalent metals {Symbol: Molar Mass}
    divalent_metals = {
        'Fe': 55.8,
        'Ni': 58.7,
        'Zn': 65.4,
        'Sn': 118.7,
        'Pb': 207.2,
        'Cd': 112.4
    }

    best_fit = {'metal': None, 'diff': float('inf')}

    print("--- Testing common divalent metals for A ---")
    for symbol, M_A in divalent_metals.items():
        E_A = M_A / 2.0
        Ex1 = get_Ex_from_salts(E_A)
        Ex2 = get_Ex_from_plate(E_A)
        diff = abs(Ex1 - Ex2)
        print(f"Testing A = {symbol} (M_A = {M_A:.1f}, E_A = {E_A:.2f}):")
        print(f"  E_X from salt data: {Ex1:.2f}")
        print(f"  E_X from plate data: {Ex2:.2f}")
        print(f"  Difference: {diff:.4f}")
        
        if diff < best_fit['diff']:
            best_fit['metal'] = symbol
            best_fit['M_A'] = M_A
            best_fit['E_A'] = E_A
            best_fit['E_X'] = np.mean([Ex1, Ex2])
            best_fit['diff'] = diff

    print("\n--- Results ---")
    metal_A_symbol = best_fit['metal']
    M_A = best_fit['M_A']
    E_X = best_fit['E_X']

    print(f"The metal 'A' that provides the most consistent result is {metal_A_symbol}.")
    print(f"Molar mass of A (M_A) = {M_A:.1f} g/mol")
    print(f"The calculated equivalent mass of unknown metal X (E_X) is approximately {E_X:.2f} g/mol.")
    
    # Find a plausible metal X. E_X = M_X / v. So M_X = E_X * v
    # Try common valencies v = 1, 2, 3
    print("\nFinding potential identity for metal X (M_X = E_X * v):")
    v1_mass = E_X * 1
    v2_mass = E_X * 2
    v3_mass = E_X * 3
    print(f"If v=1, M_X = {v1_mass:.1f} g/mol (No common stable metal)")
    print(f"If v=2, M_X = {v2_mass:.1f} g/mol (Close to Rubidium, but Rb is monovalent)")
    print(f"If v=3, M_X = {v3_mass:.1f} g/mol (Close to Tellurium, M=127.6 g/mol)")

    metal_X_symbol = 'Te'
    metal_X_valency = 3
    M_X = 127.6

    print(f"\nConclusion: Metal X is likely Tellurium (Te) with valency 3.")
    
    # The reaction involves Metal A (divalent) and Metal X (trivalent).
    # To balance the electrons (charge):
    # A -> A^2+ + 2e-
    # X^3+ + 3e- -> X
    # We need 3 moles of A for every 2 moles of X.
    # 3 A + 2 XCl3 -> 3 ACl2 + 2 X
    
    print("\nFinal Determined Reaction Equation:")
    print(f"Metal A: {metal_A_symbol}")
    print(f"Metal X: {metal_X_symbol}")
    
    a_coeff = 3
    b_coeff = 2
    c_coeff = 3
    d_coeff = 2
    
    print(f"{a_coeff} {metal_A_symbol} + {b_coeff} {metal_X_symbol}Cl{metal_X_valency} -> {c_coeff} {metal_A_symbol}Cl2 + {d_coeff} {metal_X_symbol}")

solve_chemistry_problem()
<<<3 Cd + 2 TeCl3 -> 3 CdCl2 + 2 Te>>>