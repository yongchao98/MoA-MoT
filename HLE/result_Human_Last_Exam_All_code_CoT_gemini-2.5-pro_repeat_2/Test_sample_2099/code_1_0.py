import numpy as np

def solve_chemistry_puzzle():
    """
    Solves the puzzle by deciphering the elements, calculating energies
    and Moran's I values, and finding the final product.
    """

    # Step 1: Define element data (Z: (Symbol, Mass, Electronegativity))
    # Pauling electronegativity for noble gases is considered 0.
    element_data = {
        1: ('H', 1.008, 2.20), 2: ('He', 4.0026, 0), 3: ('Li', 6.94, 0.98),
        4: ('Be', 9.012, 1.57), 5: ('B', 10.81, 2.04), 6: ('C', 12.011, 2.55),
        7: ('N', 14.007, 3.04), 8: ('O', 15.999, 3.44), 9: ('F', 18.998, 3.98),
        10: ('Ne', 20.180, 0), 11: ('Na', 22.990, 0.93), 12: ('Mg', 24.305, 1.31),
        13: ('Al', 26.982, 1.61), 14: ('Si', 28.085, 1.90)
    }

    # Step 2: Define the molecules as lists of atomic numbers (Z) based on Xn -> Z=n
    # The space in Y2 from the prompt is assumed to be a typo and ignored.
    molecules_z = {
        "Y1": [1, 2, 3, 4, 5, 6],
        "Y2": [8, 9, 10, 11, 5, 6, 12, 8, 9, 8, 12, 13],
        "Y3": [11, 3, 4, 14, 3, 4, 5, 6],
        "Y4": [12, 4, 13, 5, 6, 3],
        "Y5": [8, 9, 10, 11, 9, 14, 5, 6, 3],
        "Y6": [1, 10, 5, 1, 9, 4, 3],
        "Y7": [8, 9, 10, 11, 12, 4, 5, 6],
        "Y8": [10, 2, 5, 13, 9, 4, 12, 4, 3],
        "Y9": [9, 14, 5, 11, 3, 4, 14, 3, 4, 3],
        "Y10": [1, 12, 1, 3, 10, 12, 13, 12, 4, 3]
    }

    def calculate_barysz_energy(atom_z_list):
        n = len(atom_z_list)
        if n <= 1: return 0
        
        masses = [element_data[z][1] for z in atom_z_list]
        
        # Assume linear chain connectivity
        adj_matrix = np.zeros((n, n))
        for i in range(n - 1):
            adj_matrix[i, i + 1] = adj_matrix[i + 1, i] = 1
            
        # Build Mass-Weighted Adjacency Matrix
        mw_adj_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if adj_matrix[i, j] == 1:
                    mw_adj_matrix[i, j] = 1 / np.sqrt(masses[i] * masses[j])
        
        eigenvalues = np.linalg.eigvals(mw_adj_matrix)
        energy = np.sum(np.abs(eigenvalues))
        return energy

    # Step 3: Find the molecule with the lowest energy
    min_energy = float('inf')
    min_energy_molecule_name = None
    min_energy_molecule_z_list = None

    for name, z_list in molecules_z.items():
        energy = calculate_barysz_energy(z_list)
        if energy < min_energy:
            min_energy = energy
            min_energy_molecule_name = name
            min_energy_molecule_z_list = z_list

    # Step 4: Calculate Moran's I for the identified molecule
    def calculate_morans_i(atom_z_list, property_name):
        n = len(atom_z_list)
        if n <= 1: return 0

        if property_name == 'z':
            x = np.array(atom_z_list, dtype=float)
        elif property_name == 'mass':
            x = np.array([element_data[z][1] for z in atom_z_list])
        elif property_name == 'en':
            x = np.array([element_data[z][2] for z in atom_z_list])
        
        x_bar = np.mean(x)
        if np.var(x) == 0: return 0 # No variance in property, Moran's I is undefined/0

        masses = [element_data[z][1] for z in atom_z_list]
        adj_matrix = np.zeros((n, n))
        for i in range(n - 1):
            adj_matrix[i, i + 1] = adj_matrix[i + 1, i] = 1
        
        W = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if adj_matrix[i, j] == 1:
                    W[i, j] = 1 / np.sqrt(masses[i] * masses[j])
        
        sum_W = np.sum(W)
        if sum_W == 0: return 0

        numerator = np.sum(W * np.outer(x - x_bar, x - x_bar))
        denominator = np.sum((x - x_bar)**2)
        
        return (n / sum_W) * (numerator / denominator)

    properties = ['z', 'mass', 'en']
    moran_i_values = [calculate_morans_i(min_energy_molecule_z_list, p) for p in properties]
    
    min_moran_i = min(moran_i_values)
    max_moran_i = max(moran_i_values)

    # Step 5: Final calculation and output
    final_product = min_energy * min_moran_i * max_moran_i
    
    print(f"The molecule with the lowest Mass-Weighted Barysz Graph Energy is {min_energy_molecule_name}.")
    print(f"Identified Energy (E): {min_energy}")
    print(f"Moran's I values for properties (Z, Mass, EN): {moran_i_values}")
    print(f"Minimum Moran's I (I_min): {min_moran_i}")
    print(f"Maximum Moran's I (I_max): {max_moran_i}")
    print("\nFinal Calculation:")
    print(f"E * I_min * I_max = {min_energy} * {min_moran_i} * {max_moran_i} = {final_product}")
    print(f"<<<{final_product}>>>")

solve_chemistry_puzzle()