import math

def solve_crystallography_error():
    """
    Calculates crystallographic densities from lattice parameters and compares
    them to reported values to find the dataset with an error.
    """
    # Avogadro's number
    NA = 6.02214076e23

    # Store all data in a list of dictionaries
    datasets = [
        {'name': 'A', 'system': 'triclinic', 'M': 608.58, 'Z': 1, 'Dc_rep': 1.296,
         'a': 7.7810, 'b': 7.9273, 'c': 14.5543, 'alpha': 75.197, 'beta': 88.156, 'gamma': 64.398},
        {'name': 'B', 'system': 'triclinic', 'M': 2568.09, 'Z': 1, 'Dc_rep': 1.629,
         'a': 11.7229, 'b': 14.2639, 'c': 15.9549, 'alpha': 93.364, 'beta': 100.301, 'gamma': 91.589},
        {'name': 'C', 'system': 'monoclinic', 'M': 1365.24, 'Z': 4, 'Dc_rep': 1.574,
         'a': 17.7430, 'b': 16.0855, 'c': 20.9134, 'beta': 105.193},
        {'name': 'D', 'system': 'monoclinic', 'M': 2804.61, 'Z': 2, 'Dc_rep': 1.556,
         'a': 15.5265, 'b': 23.9138, 'c': 17.7749, 'beta': 114.893},
        {'name': 'E', 'system': 'orthorhombic', 'M': 1530.12, 'Z': 4, 'Dc_rep': 1.701,
         'a': 7.5560, 'b': 28.392, 'c': 27.854}
    ]

    max_relative_error = -1.0
    error_dataset_name = None
    
    print("Checking for inconsistencies in crystal structure data...")
    print("Plan: Recalculate density (Dc) from lattice parameters and molar mass, then compare to the reported Dc.")
    print("-" * 70)

    for data in datasets:
        name = data['name']
        system = data['system']
        a, b, c = data['a'], data['b'], data['c']
        
        # Step 1: Calculate unit cell volume (U) from lattice parameters
        U_calc = 0.0
        if system == 'triclinic':
            alpha_r, beta_r, gamma_r = math.radians(data['alpha']), math.radians(data['beta']), math.radians(data['gamma'])
            cos_a, cos_b, cos_g = math.cos(alpha_r), math.cos(beta_r), math.cos(gamma_r)
            U_calc = a * b * c * math.sqrt(1 - cos_a**2 - cos_b**2 - cos_g**2 + 2 * cos_a * cos_b * cos_g)
        elif system == 'monoclinic':
            beta_r = math.radians(data['beta'])
            U_calc = a * b * c * math.sin(beta_r)
        elif system == 'orthorhombic':
            U_calc = a * b * c

        # Step 2: Calculate density (Dc) using the calculated volume
        # Formula: Dc = (Z * M) / (U * NA)
        # Units: M(g/mol), U(Å³ -> cm³ by 1e-24), NA(mol⁻¹), Dc(g/cm³ == Mg/m⁻³)
        M = data['M']
        Z = data['Z']
        Dc_calc = (Z * M) / (U_calc * 1e-24 * NA)

        # Step 3: Compare calculated Dc with reported Dc
        Dc_rep = data['Dc_rep']
        relative_error = abs(Dc_calc - Dc_rep) / Dc_rep
        
        print(f"Dataset {name}:")
        print(f"  Calculated Volume = {U_calc:.2f} Å³")
        print(f"  M = {M}, Z = {Z}")
        print(f"  Resulting Calculated Density Dc = ({Z} * {M}) / ({U_calc:.2f} * {NA:.4e} * 1e-24) = {Dc_calc:.4f} Mg m⁻³")
        print(f"  Reported Density Dc = {Dc_rep:.4f} Mg m⁻³")
        print(f"  Relative Error = {relative_error:.4%}")
        print("-" * 70)

        if relative_error > max_relative_error:
            max_relative_error = relative_error
            error_dataset_name = name

    print(f"\nCONCLUSION:")
    print(f"The dataset with the largest relative error between its reported density and the density calculated from its fundamental parameters is '{error_dataset_name}'.")
    print(f"This indicates a significant inconsistency in the data provided for dataset {error_dataset_name}, suggesting it contains the mistake.")

solve_crystallography_error()
<<<A>>>