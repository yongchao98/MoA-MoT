import numpy as np
from mendeleev import element

def solve_chemistry_puzzle():
    """
    Solves the puzzle by deciphering the elements, calculating the specified
    physicochemical properties for each molecule, finding the molecule with the
    lowest energy, and computing the final product.
    """
    
    # Step 1: Decipher the elements. Assume Xn is the element with atomic number n.
    elements = {i: element(i) for i in range(1, 15)}

    # The molecular formulas as sequences of atomic numbers
    formulas_str = [
        "1,2,3,4,5,6",
        "8,9,10,11,5,6,12,8,9,8,12,13",
        "11,3,4,14,3,4,5,6",
        "12,4,13,5,6,3",
        "8,9,10,11,9,14,5,6,3",
        "1,10,5,1,9,4,3",
        "8,9,10,11,12,4,5,6",
        "10,2,5,13,9,4,12,4,3",
        "9,14,5,11,3,4,14,3,4,3",
        "1,12,1,3,10,12,13,12,4,3"
    ]
    
    formulas = [[int(x) for x in s.split(',')] for s in formulas_str]

    results = []

    # Step 2 & 3: Iterate through each molecule and calculate its properties
    for y_idx, atom_indices in enumerate(formulas):
        masses = np.array([elements[i].mass for i in atom_indices])
        n_atoms = len(masses)

        # Calculate Mass-Weighted Graph Energy
        adj_matrix_weighted = np.zeros((n_atoms, n_atoms))
        for i in range(n_atoms - 1):
            # Assuming a linear chain structure
            weight = np.sqrt(masses[i] * masses[i+1])
            adj_matrix_weighted[i, i+1] = weight
            adj_matrix_weighted[i+1, i] = weight
        
        eigenvalues = np.linalg.eigvalsh(adj_matrix_weighted)
        energy = np.sum(np.abs(eigenvalues))

        # Calculate Min/Max of Mass-Weighted Moran's I
        mean_mass = np.mean(masses)
        deviations = masses - mean_mass
        
        local_morans = []
        for i in range(n_atoms):
            neighbor_dev_sum = 0
            # Sum deviations of neighbors in the chain
            if i > 0:
                neighbor_dev_sum += deviations[i-1]
            if i < n_atoms - 1:
                neighbor_dev_sum += deviations[i+1]
            
            # Local Moran's I is proportional to this product
            local_i = deviations[i] * neighbor_dev_sum
            local_morans.append(local_i)
            
        min_moran = min(local_morans)
        max_moran = max(local_morans)
        
        results.append({
            'y_index': y_idx + 1,
            'energy': energy,
            'min_moran': min_moran,
            'max_moran': max_moran
        })

    # Step 4: Find the molecule with the lowest energy
    min_energy_result = min(results, key=lambda x: x['energy'])
    
    # Step 5: Calculate the final product for the identified molecule
    final_product = min_energy_result['energy'] * min_energy_result['min_moran'] * min_energy_result['max_moran']

    # As requested, output the components of the final calculation
    print(f"The molecule with the lowest energy is Y{min_energy_result['y_index']}.")
    print(f"Identified Energy: {min_energy_result['energy']}")
    print(f"Minimum Moran's I: {min_energy_result['min_moran']}")
    print(f"Maximum Moran's I: {min_energy_result['max_moran']}")
    
    # Final Answer
    print(f"The final product is: {final_product}")
    
    print(f"<<<{final_product}>>>")


solve_chemistry_puzzle()