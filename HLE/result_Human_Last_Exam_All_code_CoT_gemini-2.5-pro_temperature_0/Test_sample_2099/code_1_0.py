import numpy as np

def solve_saber_puzzle():
    """
    Solves the entire puzzle from deciphering to final calculation.
    """
    # Step 1: Deciphering
    atomic_mass = {'Fe': 55.845, 'C': 12.011}
    x_map = {i: 'Fe' if i % 2 == 1 else 'C' for i in range(1, 15) if i != 7}

    y_definitions = {
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

    y_sequences = {name: [x_map[i] for i in indices] for name, indices in y_definitions.items()}

    # Step 2: Calculate Mass-Weighted Barysz Graph Energy
    def calculate_barysz_energy(sequence):
        n = len(sequence)
        if n == 0: return 0
        masses = np.array([atomic_mass[atom] for atom in sequence])
        b_matrix = np.zeros((n, n))
        for i in range(n):
            b_matrix[i, i] = masses[i]
            for j in range(i + 1, n):
                dist = float(j - i)
                val = -1.0 / (np.sqrt(masses[i] * masses[j]) * dist)
                b_matrix[i, j] = b_matrix[j, i] = val
        eigenvalues = np.linalg.eigvalsh(b_matrix)
        return np.sum(np.abs(eigenvalues))

    energies = {name: calculate_barysz_energy(seq) for name, seq in y_sequences.items()}

    # Step 3: Identify the element Y with the lowest energy
    min_energy = float('inf')
    identified_y_name = None
    for name, energy in energies.items():
        if energy < min_energy:
            min_energy = energy
            identified_y_name = name
    
    identified_energy = min_energy
    identified_y_sequence = y_sequences[identified_y_name]

    # Step 4: Calculate min and max values of Mass-Weighted Moran's I
    def calculate_moran_i_bounds(sequence):
        n = len(sequence)
        if n <= 1: return 0, 0
        w_matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i != j:
                    w_matrix[i, j] = 1.0 / (float(abs(i - j))**2)
        
        identity_matrix = np.identity(n)
        ones_matrix = np.ones((n, n))
        centering_matrix = identity_matrix - (ones_matrix / n)
        m_matrix = centering_matrix @ w_matrix @ centering_matrix
        eigenvalues = np.linalg.eigvalsh(m_matrix)
        return np.min(eigenvalues), np.max(eigenvalues)

    min_moran_i, max_moran_i = calculate_moran_i_bounds(identified_y_sequence)

    # Step 5: Final Calculation and Output
    final_product = identified_energy * min_moran_i * max_moran_i

    print(f"The identified element is {identified_y_name}, which has the lowest energy.")
    print("\nFinal equation components:")
    print(f"Lowest Mass-Weighted Barysz Graph Energy = {identified_energy}")
    print(f"Minimum value of Mass-Weighted Moran's I = {min_moran_i}")
    print(f"Maximum value of Mass-Weighted Moran's I = {max_moran_i}")
    print("\nFinal calculation:")
    print(f"{identified_energy} * {min_moran_i} * {max_moran_i} = {final_product}")
    
    # Output the final answer in the required format
    print(f"\n<<<{final_product}>>>")

solve_saber_puzzle()