import numpy as np
import networkx as nx

def solve():
    """
    Main function to solve the user's request.
    """
    # Step 1: Decipher the elements
    # Key from "1880" -> 1+8+8+0 = 17. X_i has atomic number i + 17.
    # We will use element symbols and their standard atomic weights.
    atomic_weights = {
        'Ar': 39.948, 'K': 39.0983, 'Ca': 40.078, 'Sc': 44.955912, 'Ti': 47.867,
        'V': 50.9415, 'Mn': 54.938045, 'Fe': 55.845, 'Co': 58.933195,
        'Ni': 58.6934, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723
    }

    # Mapping from X_i index to element symbol
    # X_i -> Z = i + 17
    # 1:Ar, 2:K, 3:Ca, 4:Sc, 5:Ti, 6:V, (7:Cr not used), 8:Mn, 9:Fe, 10:Co, 11:Ni, 12:Cu, 13:Zn, 14:Ga
    x_map = {
        1: 'Ar', 2: 'K', 3: 'Ca', 4: 'Sc', 5: 'Ti', 6: 'V',
        8: 'Mn', 9: 'Fe', 10: 'Co', 11: 'Ni', 12: 'Cu', 13: 'Zn', 14: 'Ga'
    }

    # Step 2: Define the molecules as lists of X indices
    y_indices = {
        "Y1": [1, 2, 3, 4, 5, 6],
        "Y2": [8, 9, 10, 11, 5, 6, 12, 8, 9, 8, 12, 13],
        "Y3": [11, 3, 4, 14, 3, 4, 5, 6],
        "Y4": [12, 4, 13, 5, 6, 3],
        "Y5": [8, 9, 10, 11, 9, 14, 5, 6, 3],
        "Y6": [1, 10, 5, 1, 9, 4, 3],
        "Y7": [8, 9, 10, 11, 12, 4, 5, 6],
        "Y8": [10, 2, 5, 13, 9, 4, 12, 4, 3],
        "Y9": [9, 14, 5, 11, 3, 4, 14, 3, 4, 3],
        "Y10": [1, 12, 1, 3, 10, 12, 13, 12, 4, 3],
    }

    # Convert indices to element symbols
    molecules = {name: [x_map[i] for i in indices] for name, indices in y_indices.items()}

    energies = {}

    # Step 3: Calculate Mass-Weighted Barysz Graph Energy for each molecule
    for name, atoms in molecules.items():
        n = len(atoms)
        if n == 0:
            continue
            
        G = nx.path_graph(n)
        masses = np.array([atomic_weights[atom] for atom in atoms])
        
        # Calculate distance matrix d_ij
        dist_matrix = np.array(nx.floyd_warshall_numpy(G))
        
        # Calculate Mass-Weighted Barysz Distance Matrix
        barysz_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                m_i, m_j = masses[i], masses[j]
                val = (1 - min(m_i, m_j) / max(m_i, m_j)) * dist_matrix[i, j]
                barysz_matrix[i, j] = barysz_matrix[j, i] = val
        
        # Calculate eigenvalues and energy
        eigenvalues = np.linalg.eigvalsh(barysz_matrix)
        energy = np.sum(np.abs(eigenvalues))
        energies[name] = energy

    # Step 4: Identify the molecule with the lowest energy
    min_energy_molecule_name = min(energies, key=energies.get)
    min_energy = energies[min_energy_molecule_name]
    
    # Step 5: For the identified molecule, calculate Mass-Weighted Moran's I
    target_atoms = molecules[min_energy_molecule_name]
    n = len(target_atoms)
    G = nx.path_graph(n)
    masses = np.array([atomic_weights[atom] for atom in target_atoms])
    dist_matrix = np.array(nx.floyd_warshall_numpy(G))

    mean_mass = np.mean(masses)
    mass_deviations = masses - mean_mass
    mass_dev_sq_sum = np.sum(mass_deviations**2)

    moran_values = []
    
    # Calculate I(d) for all possible distances d
    max_dist = int(np.max(dist_matrix))
    for d in range(1, max_dist + 1):
        W_ij_num_sum = 0
        W_ij_den_sum = 0
        
        # Create W_ij matrix where W_ij=1 if dist(i,j)==d
        pairs = np.argwhere(dist_matrix == d)
        
        if len(pairs) == 0:
            continue

        S0 = len(pairs)
            
        numerator_sum = 0
        for i, j in pairs:
            numerator_sum += mass_deviations[i] * mass_deviations[j]
            
        if mass_dev_sq_sum == 0 or S0 == 0:
            moran_I = 0.0 # Or handle as undefined, but 0 is safe here
        else:
            moran_I = (n / S0) * (numerator_sum / mass_dev_sq_sum)
        
        moran_values.append(moran_I)

    min_moran = min(moran_values)
    max_moran = max(moran_values)

    # Step 6: Compute the final product
    final_product = min_energy * min_moran * max_moran
    
    # Print the equation as requested
    print(f"The identified molecule is {min_energy_molecule_name}")
    print("Final equation:")
    print(f"{min_energy} (Energy) * {min_moran} (Min Moran's I) * {max_moran} (Max Moran's I) = {final_product}")
    
    # Final answer block
    print(f"<<<{final_product}>>>")

solve()