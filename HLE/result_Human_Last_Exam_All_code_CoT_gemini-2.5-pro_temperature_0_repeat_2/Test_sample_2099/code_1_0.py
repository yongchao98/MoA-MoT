import numpy as np

def calculate_barysz_energy(masses):
    """
    Calculates the Mass-Weighted Barysz Graph Energy for a molecule.
    Assumes a clique structure and a Barysz matrix B_ij = m_i * m_j.
    """
    n = len(masses)
    if n == 0:
        return 0
    m = np.array(masses, dtype=float)
    
    # Create the Barysz matrix B_ij = m_i * m_j
    barysz_matrix = np.outer(m, m)
    # Set diagonal to zero
    np.fill_diagonal(barysz_matrix, 0)
    
    # Energy is the sum of the absolute values of the eigenvalues
    # Use eigvalsh for symmetric matrices, it's faster and returns real eigenvalues.
    eigenvalues = np.linalg.eigvalsh(barysz_matrix)
    energy = np.sum(np.abs(eigenvalues))
    return energy

def calculate_moran_min_max(masses):
    """
    Calculates the min and max values of the term (m_i - m_bar)(m_j - m_bar)
    for the Moran's I calculation.
    """
    n = len(masses)
    if n < 2:
        return 0, 0
    m = np.array(masses, dtype=float)
    m_bar = np.mean(m)
    
    deviations = m - m_bar
    
    # The product term matrix contains (m_i - m_bar)(m_j - m_bar) for all pairs i, j
    term_matrix = np.outer(deviations, deviations)
    
    min_val = np.min(term_matrix)
    max_val = np.max(term_matrix)
    
    return min_val, max_val

def solve_puzzle():
    """
    Main function to solve the puzzle.
    """
    # Step 1: Decipher X_i as having mass 'i'. Define Y_k molecules.
    # The space in Y2 is interpreted as a typo and ignored.
    formulas = {
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

    min_energy = float('inf')
    min_energy_y_name = None
    
    # Step 2 & 3: Calculate energy for each Y and find the minimum.
    for name, masses in formulas.items():
        energy = calculate_barysz_energy(masses)
        if energy < min_energy:
            min_energy = energy
            min_energy_y_name = name

    # The identified element Y is the one with the lowest energy
    identified_y_masses = formulas[min_energy_y_name]
    
    # Step 4: Calculate Moran's I min/max for the identified Y
    min_moran_val, max_moran_val = calculate_moran_min_max(identified_y_masses)
    
    # Step 5: Calculate the final product
    final_product = min_energy * min_moran_val * max_moran_val
    
    # Print the final equation as requested
    print(f"Identified element: {min_energy_y_name}")
    print(f"Lowest Mass-Weighted Barysz Graph Energy: {min_energy}")
    print(f"Minimum Moran's I term value: {min_moran_val}")
    print(f"Maximum Moran's I term value: {max_moran_val}")
    print("Final calculation:")
    print(f"{min_energy} * {min_moran_val} * {max_moran_val} = {final_product}")
    
    # Return the final answer in the specified format
    print(f"<<<{final_product}>>>")

solve_puzzle()