import numpy as np

def solve_puzzle():
    """
    Solves the chemistry puzzle by deciphering the molecules, calculating their
    energies, finding the molecule with the lowest energy, and then computing
    the product of this energy with its min and max Moran's I values.
    """
    
    # Step 1 & 2: Decipher the mapping and define element properties
    # This mapping is the solution to the cipher part of the puzzle.
    X_map = {
        'X1': 'P', 'X2': 'I', 'X3': 'S', 'X4': 'C', 'X5': 'H',
        'X6': 'N', 'X7': 'O', 'X8': 'B', 'X9': 'Se', 'X10': 'Br',
        'X11': 'Cl', 'X12': 'As', 'X13': 'Ge', 'X14': 'F'
    }
    
    # Atomic masses (M) and atomic numbers (Z) for the elements
    properties = {
        'H': {'M': 1.008, 'Z': 1}, 'C': {'M': 12.011, 'Z': 6},
        'N': {'M': 14.007, 'Z': 7}, 'O': {'M': 15.999, 'Z': 8},
        'F': {'M': 18.998, 'Z': 9}, 'P': {'M': 30.974, 'Z': 15},
        'S': {'M': 32.06, 'Z': 16}, 'Cl': {'M': 35.45, 'Z': 17},
        'B': {'M': 10.81, 'Z': 5}, 'Ge': {'M': 72.630, 'Z': 32},
        'As': {'M': 74.922, 'Z': 33}, 'Se': {'M': 78.971, 'Z': 34},
        'Br': {'M': 79.904, 'Z': 35}, 'I': {'M': 126.90, 'Z': 53},
    }

    # Step 3: Define the molecules as sequences of X variables
    Y_formulas = {
        'Y1': ['X1', 'X2', 'X3', 'X4', 'X5', 'X6'],
        'Y2': ['X8', 'X9', 'X10', 'X11', 'X5', 'X6', 'X12', 'X8', 'X9', 'X8', 'X12', 'X13'],
        'Y3': ['X11', 'X3', 'X4', 'X14', 'X3', 'X4', 'X5', 'X6'],
        'Y4': ['X12', 'X4', 'X13', 'X5', 'X6', 'X3'],
        'Y5': ['X8', 'X9', 'X10', 'X11', 'X9', 'X14', 'X5', 'X6', 'X3'],
        'Y6': ['X1', 'X10', 'X5', 'X1', 'X9', 'X4', 'X3'],
        'Y7': ['X8', 'X9', 'X10', 'X11', 'X12', 'X4', 'X5', 'X6'],
        'Y8': ['X10', 'X2', 'X5', 'X13', 'X9', 'X4', 'X12', 'X4', 'X3'],
        'Y9': ['X9', 'X14', 'X5', 'X11', 'X3', 'X4', 'X14', 'X3', 'X4', 'X3'],
        'Y10': ['X1', 'X12', 'X1', 'X3', 'X10', 'X12', 'X13', 'X12', 'X4', 'X3'],
    }

    energies = {}

    # Step 4: Calculate Mass-Weighted Barysz Graph Energy for each molecule
    for name, formula in Y_formulas.items():
        atom_list = [X_map[x] for x in formula]
        n_atoms = len(atom_list)
        masses = np.array([properties[atom]['M'] for atom in atom_list])
        
        # Build the Barysz matrix B
        B = np.zeros((n_atoms, n_atoms))
        mass_sqrt_matrix = np.sqrt(np.outer(masses, masses))
        
        for i in range(n_atoms):
            B[i, i] = masses[i]
            if i < n_atoms - 1:
                B[i, i + 1] = mass_sqrt_matrix[i, i + 1]
                B[i + 1, i] = mass_sqrt_matrix[i + 1, i]

        eigenvalues, _ = np.linalg.eig(B)
        energy = np.sum(np.abs(eigenvalues))
        energies[name] = energy

    # Step 5: Identify the target molecule with the lowest energy
    min_energy_name = min(energies, key=energies.get)
    min_energy_value = energies[min_energy_name]
    
    # Step 6: Calculate Mass-Weighted Moran's I for the target molecule
    formula_min = Y_formulas[min_energy_name]
    atom_list_min = [X_map[x] for x in formula_min]
    masses_min = np.array([properties[atom]['M'] for atom in atom_list_min])
    n_atoms_min = len(atom_list_min)
    
    mean_mass = np.mean(masses_min)
    mass_deviations = masses_min - mean_mass
    sum_sq_dev = np.sum(mass_deviations ** 2)

    # Calculate topological distance matrix for the chain
    dist_matrix = np.abs(np.arange(n_atoms_min)[:, None] - np.arange(n_atoms_min)[None, :])
    
    moran_I_values = []
    max_dist = n_atoms_min - 1

    for d in range(1, max_dist + 1):
        numerator = 0
        weight_sum = 0
        
        # Find pairs of atoms at distance d
        pairs = np.argwhere(dist_matrix == d)
        if len(pairs) == 0:
            continue
            
        weight_sum = len(pairs)
        for i, j in pairs:
            numerator += mass_deviations[i] * mass_deviations[j]

        if sum_sq_dev > 1e-9 and weight_sum > 0:
            moran_I = (n_atoms_min / weight_sum) * (numerator / sum_sq_dev)
            moran_I_values.append(moran_I)

    # Step 7: Find the min and max Moran's I and the final product
    min_moran_I = min(moran_I_values)
    max_moran_I = max(moran_I_values)
    
    final_product = min_energy_value * min_moran_I * max_moran_I

    # Final output
    print(f"Identified molecule with the lowest energy: {min_energy_name}")
    print(f"The lowest Mass-Weighted Barysz Graph Energy is: {min_energy_value}")
    print(f"Minimum Mass-Weighted Moran's I is: {min_moran_I}")
    print(f"Maximum Mass-Weighted Moran's I is: {max_moran_I}")
    print("The final result (Energy * I_min * I_max) is:")
    print(f"{min_energy_value} * {min_moran_I} * {max_moran_I} = {final_product}")
    
    # Return final answer in specified format
    print(f"\n<<<{final_product}>>>")

solve_puzzle()