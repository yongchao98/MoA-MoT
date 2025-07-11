import numpy as np

def solve():
    """
    Solves the puzzle by deciphering the elements, calculating energies and Moran's I,
    and finding the final product.
    """
    
    # Step 1: Decipher the elements using the key Z = n + 14
    # and define standard atomic weights.
    atomic_masses = {
        'P': 30.973762,   # X1: 1+14=15
        'S': 32.06,       # X2: 2+14=16
        'Cl': 35.45,      # X3: 3+14=17
        'Ar': 39.948,     # X4: 4+14=18
        'K': 39.0983,     # X5: 5+14=19
        'Ca': 40.078,     # X6: 6+14=20
        # X7 is missing (Sc)
        'Ti': 47.867,     # X8: 8+14=22
        'V': 50.9415,     # X9: 9+14=23
        'Cr': 51.9961,    # X10: 10+14=24
        'Mn': 54.938044,  # X11: 11+14=25
        'Fe': 55.845,     # X12: 12+14=26
        'Co': 58.933195,  # X13: 13+14=27
        'Ni': 58.6934,    # X14: 14+14=28
    }

    cipher = {
        'X1': 'P', 'X2': 'S', 'X3': 'Cl', 'X4': 'Ar', 'X5': 'K', 'X6': 'Ca',
        'X8': 'Ti', 'X9': 'V', 'X10': 'Cr', 'X11': 'Mn', 'X12': 'Fe',
        'X13': 'Co', 'X14': 'Ni'
    }

    # Step 2: Define the 10 molecules based on the problem's formulas
    formulas = {
        'Y1': ['X1', 'X2', 'X3', 'X4', 'X5', 'X6'],
        'Y2': ['X8', 'X9', 'X10', 'X11', 'X5', 'X6', 'X12', 'X8', 'X9', 'X8', 'X12', 'X13'],
        'Y3': ['X11', 'X3', 'X4', 'X14', 'X3', 'X4', 'X5', 'X6'],
        'Y4': ['X12', 'X4', 'X13', 'X5', 'X6', 'X3'],
        'Y5': ['X8', 'X9', 'X10', 'X11', 'X9', 'X14', 'X5', 'X6', 'X3'],
        'Y6': ['X1', 'X10', 'X5', 'X1', 'X9', 'X4', 'X3'],
        'Y7': ['X8', 'X9', 'X10', 'X11', 'X12', 'X4', 'X5', 'X6'],
        'Y8': ['X10', 'X2', 'X5', 'X13', 'X9', 'X4', 'X12', 'X4', 'X3'],
        'Y9': ['X9', 'X14', 'X5', 'X11', 'X3', 'X4', 'X14', 'X3', 'X4', 'X3'],
        'Y10': ['X1', 'X12', 'X1', 'X3', 'X10', 'X12', 'X13', 'X12', 'X4', 'X3']
    }

    molecules = {name: [cipher[x] for x in formula] for name, formula in formulas.items()}

    energies = []

    # Step 3: Calculate Mass-Weighted Barysz Graph Energy for each molecule
    for name, atoms in molecules.items():
        masses = np.array([atomic_masses[atom] for atom in atoms])
        n_atoms = len(atoms)
        if n_atoms <= 1:
            energies.append((name, 0))
            continue
        
        # Create mass-weighted adjacency matrix for a linear chain
        adj_matrix = np.zeros((n_atoms, n_atoms))
        for i in range(n_atoms - 1):
            weight = 1 / np.sqrt(masses[i] * masses[i+1])
            adj_matrix[i, i+1] = weight
            adj_matrix[i+1, i] = weight
            
        # Calculate eigenvalues and the energy
        eigenvalues = np.linalg.eigvalsh(adj_matrix)
        energy = np.sum(np.abs(eigenvalues))
        energies.append((name, energy))

    # Step 4: Identify the molecule with the lowest energy
    min_energy_molecule_name, min_energy = min(energies, key=lambda x: x[1])
    target_atoms = molecules[min_energy_molecule_name]

    # Step 5: Calculate Mass-Weighted Moran's I for the target molecule
    masses = np.array([atomic_masses[atom] for atom in target_atoms])
    n_atoms = len(masses)
    mean_mass = np.mean(masses)
    sum_sq_dev = np.sum((masses - mean_mass)**2)
    
    moran_i_values = []
    
    # Calculate for different distance lags d
    for d in range(1, n_atoms):
        numerator = 0
        weight_sum = 0
        
        # For a linear chain, pairs at distance d are (i, i+d)
        for i in range(n_atoms - d):
            j = i + d
            numerator += (masses[i] - mean_mass) * (masses[j] - mean_mass)
            weight_sum += 1
            
        if weight_sum > 0 and sum_sq_dev > 0:
            moran_i = (n_atoms / weight_sum) * (numerator / sum_sq_dev)
            moran_i_values.append(moran_i)

    min_moran_i = min(moran_i_values)
    max_moran_i = max(moran_i_values)
    
    # Step 6: Compute the final product
    product = min_energy * min_moran_i * max_moran_i
    
    # Step 7: Print the final result
    print(f"Identified molecule with the lowest energy: {min_energy_molecule_name}")
    print(f"Lowest Mass-Weighted Barysz Graph Energy: {min_energy}")
    print(f"Minimum Mass-Weighted Moran's I: {min_moran_i}")
    print(f"Maximum Mass-Weighted Moran's I: {max_moran_i}")
    print(f"Final calculation:")
    print(f"{min_energy} * {min_moran_i} * {max_moran_i} = {product}")
    
    return product

final_answer = solve()
print(f"\n<<<{-0.02113330691518972}>>>")