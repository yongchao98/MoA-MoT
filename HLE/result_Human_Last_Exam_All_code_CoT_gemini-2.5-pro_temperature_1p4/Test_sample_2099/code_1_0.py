import numpy as np

def solve_chemistry_puzzle():
    """
    Solves the puzzle by deciphering the elements, calculating graph energies,
    finding the minimum, and then computing Moran's I values for the final product.
    """
    # Step 1: Decipher the system using the standard key and define atomic properties
    ELEMENT_DATA = {
        'H':  {'Z': 1,  'A': 1.008}, 'C':  {'Z': 6,  'A': 12.011}, 'N':  {'Z': 7,  'A': 14.007},
        'O':  {'Z': 8,  'A': 15.999}, 'F':  {'Z': 9,  'A': 18.998}, 'P':  {'Z': 15, 'A': 30.974},
        'V':  {'Z': 23, 'A': 50.942}, 'Cr': {'Z': 24, 'A': 51.996}, 'Mn': {'Z': 25, 'A': 54.938},
        'Fe': {'Z': 26, 'A': 55.845}, 'Co': {'Z': 27, 'A': 58.933}, 'Ni': {'Z': 28, 'A': 58.693},
        'Cu': {'Z': 29, 'A': 63.546}
    }
    X_MAP = {
        1: 'H', 2: 'C', 3: 'N', 4: 'O', 5: 'F', 6: 'P', 8: 'V', 9: 'Cr',
        10: 'Mn', 11: 'Fe', 12: 'Co', 13: 'Ni', 14: 'Cu'
    }

    Y_FORMULAS = {
        'Y1': [1, 2, 3, 4, 5, 6],
        'Y2': [8, 9, 10, 11, 5, 6, 12, 8, 9, 8, 12, 13],
        'Y3': [11, 3, 4, 14, 3, 4, 5, 6],
        'Y4': [12, 4, 13, 5, 6, 3],
        'Y5': [8, 9, 10, 11, 9, 14, 5, 6, 3],
        'Y6': [1, 10, 5, 1, 9, 4, 3],
        'Y7': [8, 9, 10, 11, 12, 4, 5, 6],
        'Y8': [10, 2, 5, 13, 9, 4, 12, 4, 3],
        'Y9': [9, 14, 5, 11, 3, 4, 14, 3, 4, 3],
        'Y10': [1, 12, 1, 3, 10, 12, 13, 12, 4, 3]
    }

    def parse_formula(x_indices):
        """Builds graph structure from a formula sequence."""
        atom_sequence = [X_MAP[i] for i in x_indices]
        unique_atoms = sorted(list(set(atom_sequence)))
        atom_map = {name: i for i, name in enumerate(unique_atoms)}
        n_atoms = len(unique_atoms)
        adj_matrix = np.zeros((n_atoms, n_atoms))
        
        for i in range(len(atom_sequence) - 1):
            u_name, v_name = atom_sequence[i], atom_sequence[i+1]
            if u_name != v_name:
                u_idx, v_idx = atom_map[u_name], atom_map[v_name]
                adj_matrix[u_idx, v_idx] = 1
                adj_matrix[v_idx, u_idx] = 1
        return unique_atoms, adj_matrix

    def calculate_barysz_energy(atoms, adj_matrix):
        """Calculates the Mass-Weighted Barysz Graph Energy."""
        n = len(atoms)
        masses = np.array([ELEMENT_DATA[atom]['A'] for atom in atoms])
        barysz_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i + 1, n):
                if adj_matrix[i, j] == 1:
                    val = (masses[i] * masses[j])**(-0.5)
                    barysz_matrix[i, j] = val
                    barysz_matrix[j, i] = val
                    
        eigenvalues = np.linalg.eigvalsh(barysz_matrix)
        energy = np.sum(np.abs(eigenvalues))
        return energy

    def calculate_moran_i(atoms, adj_matrix, prop):
        """Calculates Mass-Weighted Moran's I for a given atomic property."""
        n = len(atoms)
        prop_values = np.array([ELEMENT_DATA[atom][prop] for atom in atoms])
        prop_mean = np.mean(prop_values)
        z = prop_values - prop_mean
        
        s_w = np.sum(adj_matrix)
        if s_w == 0: return 0

        numerator = 0
        for i in range(n):
            for j in range(n):
                numerator += adj_matrix[i, j] * z[i] * z[j]
        
        denominator = np.sum(z**2)
        if denominator == 0: return 0

        moran_i = (n / s_w) * (numerator / denominator)
        return moran_i

    # Step 2 & 3: Calculate E_BG for all molecules and find the minimum
    energies = {}
    print("Calculating Mass-Weighted Barysz Graph Energy for each molecule:")
    for name, formula in Y_FORMULAS.items():
        atoms, adj_matrix = parse_formula(formula)
        energy = calculate_barysz_energy(atoms, adj_matrix)
        energies[name] = energy
        print(f"{name}: {energy:.4f}")

    min_energy_y = min(energies, key=energies.get)
    min_energy_val = energies[min_energy_y]
    print(f"\nMolecule with the lowest energy is {min_energy_y} with E = {min_energy_val:.4f}\n")

    # Step 4: Calculate Moran's I for the identified molecule
    min_y_atoms, min_y_adj_matrix = parse_formula(Y_FORMULAS[min_energy_y])
    
    i_mass = calculate_moran_i(min_y_atoms, min_y_adj_matrix, 'A')
    i_atomic_number = calculate_moran_i(min_y_atoms, min_y_adj_matrix, 'Z')
    
    i_min = min(i_mass, i_atomic_number)
    i_max = max(i_mass, i_atomic_number)

    print(f"For {min_energy_y}:")
    print(f"Moran's I (using Atomic Mass): {i_mass:.4f}")
    print(f"Moran's I (using Atomic Number): {i_atomic_number:.4f}")
    print(f"I_min = {i_min:.4f}")
    print(f"I_max = {i_max:.4f}\n")

    # Step 5: Compute the final product
    final_product = min_energy_val * i_min * i_max
    
    print("Final calculation:")
    print(f"{min_energy_val:.4f} * {i_min:.4f} * {i_max:.4f} = {final_product:.4f}")
    
    return final_product

final_answer = solve_chemistry_puzzle()
print(f"<<<{final_answer:.4f}>>>")