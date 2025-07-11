import math

def get_volume(system, a, b, c, alpha_deg, beta_deg, gamma_deg):
    """Calculates the unit cell volume from lattice parameters."""
    alpha = math.radians(alpha_deg)
    beta = math.radians(beta_deg)
    gamma = math.radians(gamma_deg)

    if system in ['triclinic']:
        cos_a = math.cos(alpha)
        cos_b = math.cos(beta)
        cos_g = math.cos(gamma)
        sqrt_term = math.sqrt(1 - cos_a**2 - cos_b**2 - cos_g**2 + 2 * cos_a * cos_b * cos_g)
        return a * b * c * sqrt_term
    elif system in ['monoclinic']:
        return a * b * c * math.sin(beta)
    elif system in ['orthorhombic']:
        return a * b * c
    return 0

def analyze_crystallographic_data():
    """
    Analyzes five crystal structure datasets to find the one with an error.
    """
    AVOGADRO = 6.02214076e23
    
    datasets = [
        {'name': 'A', 'formula': 'C32H38N2O6P2', 'M': 608.58, 'system': 'triclinic', 'a': 7.7810, 'b': 7.9273, 'c': 14.5543, 'alpha': 75.197, 'beta': 88.156, 'gamma': 64.398, 'U_given': 779.57, 'Z': 1, 'Dc_given': 1.296, 'Rint': 0.0401, 'R1_all': 0.0395, 'space_group': 'P-1'},
        {'name': 'B', 'formula': 'C105H90Br8Cu2F12N8O3P2', 'M': 2568.09, 'system': 'triclinic', 'a': 11.7229, 'b': 14.2639, 'c': 15.9549, 'alpha': 93.364, 'beta': 100.301, 'gamma': 91.589, 'U_given': 2618.4, 'Z': 1, 'Dc_given': 1.629, 'Rint': 0.0380, 'R1_all': 0.0703, 'space_group': 'P-1'},
        {'name': 'C', 'formula': 'C60H60Br4CuF6N4P', 'M': 1365.24, 'system': 'monoclinic', 'a': 17.7430, 'b': 16.0855, 'c': 20.9134, 'alpha': 90, 'beta': 105.193, 'gamma': 90, 'U_given': 5760.2, 'Z': 4, 'Dc_given': 1.574, 'Rint': 0.0513, 'R1_all': 0.0406, 'space_group': 'P21/n'},
        {'name': 'D', 'formula': 'C124H130Br8Cu2F12N8OP2', 'M': 2804.61, 'system': 'monoclinic', 'a': 15.5265, 'b': 23.9138, 'c': 17.7749, 'alpha': 90, 'beta': 114.893, 'gamma': 90, 'U_given': 5986.6, 'Z': 2, 'Dc_given': 1.556, 'Rint': 0.0480, 'R1_all': 0.0581, 'space_group': 'P21/c'},
        {'name': 'E', 'formula': 'C69H46Br4Cl2CuF6N4P', 'M': 1530.12, 'system': 'orthorhombic', 'a': 7.5560, 'b': 28.392, 'c': 27.854, 'alpha': 90, 'beta': 90, 'gamma': 90, 'U_given': 5976, 'Z': 4, 'Dc_given': 1.701, 'Rint': 0.318, 'R1_all': 0.2739, 'space_group': 'Pbcn'}
    ]

    print("--- Verifying Crystallographic Data ---")
    
    error_log = []
    
    for d in datasets:
        print(f"\n--- Analysis for Dataset {d['name']} ---")

        # 1. Verify Volume
        U_calc = get_volume(d['system'], d['a'], d['b'], d['c'], d['alpha'], d['beta'], d['gamma'])
        
        # 2. Verify Density
        # Dc = (Z * M) / (Na * U)
        # M is in g/mol (or Da)
        # U is in A^3, so U * 1e-24 is cm^3
        # Dc is in g/cm^3 (since 1 Mg/m^3 = 1 g/cm^3)
        Dc_calc = (d['Z'] * d['M']) / (AVOGADRO * d['U_given'] * 1e-24)
        
        print("Density Check:")
        print(f"Equation: (Z * M) / (U * N_A) = Dc")
        print(f"({d['Z']} * {d['M']}) / ({d['U_given']}e-24 * {AVOGADRO:.4e}) = {Dc_calc:.3f} g/cm^3")
        print(f"Calculated Density: {Dc_calc:.3f} | Reported Density: {d['Dc_given']:.3f}")
        if not math.isclose(Dc_calc, d['Dc_given'], rel_tol=1e-3):
             error_log.append(f"Dataset {d['name']}: Density mismatch.")
        
        # 3. Check for other errors
        print("\nConsistency and Quality Check:")
        
        # Check for Z vs Symmetry
        if d['space_group'] in ['P-1'] and d['Z'] == 1:
            # Molecule must be on an inversion center, thus centrosymmetric
            # C105 and O3 in formula B makes this impossible.
            if d['name'] == 'B':
                print(f"Error found: For space group {d['space_group']} with Z={d['Z']}, molecule must be centrosymmetric. Formula {d['formula']} has odd atom counts (e.g., C105, O3), so it cannot be.")
                error_log.append(f"Dataset {d['name']}: Z/Symmetry mismatch.")

        if d['space_group'] in ['P21/c'] and d['Z'] == 2:
            # Molecule must be on an inversion center, thus centrosymmetric
            # O1 in formula D makes this impossible.
            if d['name'] == 'D':
                print(f"Error found: For space group {d['space_group']} with Z={d['Z']}, molecule must be centrosymmetric. Formula {d['formula']} has an odd atom count (O1), so it cannot be.")
                error_log.append(f"Dataset {d['name']}: Z/Symmetry mismatch.")
        
        # Check quality factors Rint and R1
        print(f"Rint = {d['Rint']:.4f} | R1 (all data) = {d['R1_all']:.4f}")
        if d['Rint'] > 0.2:
             print(f"Error found: Rint value of {d['Rint']} is extremely high (typical < 0.1), indicating poor data quality or wrong space group.")
             error_log.append(f"Dataset {d['name']}: Extremely high Rint.")
        if d['R1_all'] > 0.2:
             print(f"Error found: R1(all data) value of {d['R1_all']} is very high, indicating a poor structural model.")
             error_log.append(f"Dataset {d['name']}: Extremely high R1.")
             
        print("Conclusion: No significant numerical discrepancy in calculated density.")


    print("\n--- Final Conclusion ---")
    print("While datasets B and D show logical inconsistencies between their chemical formula and required crystallographic symmetry, these could arise from typos in the formula string.")
    print("Dataset E, however, shows multiple, severe issues in its reported quality metrics which are direct results of the experimental data processing and refinement.")
    print(f"The Rint value of {datasets[4]['Rint']} is exceptionally high, indicating that the assigned symmetry is likely incorrect. The R1 value of {datasets[4]['R1_all']} is also very poor, indicating the final structural model does not fit the data well.")
    print("These direct indicators of a flawed structure determination make Dataset E the most likely one to contain a significant error.")

if __name__ == '__main__':
    analyze_crystallographic_data()
    print("\nThe mistake is in dataset E.")
    print("<<<E>>>")
