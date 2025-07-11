import numpy as np

def solve_chemistry_puzzle():
    """
    Solves the puzzle by identifying a target molecule from a list of candidates
    based on minimum graph energy, then calculating a product of its properties.
    """
    # Step 1: Define atomic data (Symbol: (Atomic Number, Standard Atomic Weight))
    ATOMIC_DATA = {
        'H': (1, 1.008), 'B': (5, 10.81), 'C': (6, 12.011), 'N': (7, 14.007),
        'O': (8, 15.999), 'F': (9, 18.998), 'Li': (3, 6.94), 'Al': (13, 26.982),
        'Si': (14, 28.085), 'P': (15, 30.974), 'S': (16, 32.06), 'Ar': (18, 39.948),
        'Ni': (28, 58.693), 'Se': (34, 78.971), 'Ti': (22, 47.867), 'Mn': (25, 54.938),
        'Zn': (30, 65.38), 'Y': (39, 88.906), 'Ru': (44, 101.07), 'I': (53, 126.90),
        'La': (57, 138.91), 'Sm': (62, 150.36), 'Ho': (67, 164.93), 'W': (74, 183.84),
        'Th': (90, 232.04), 'U': (92, 238.03), 'Am': (95, 243.0), 'Ra': (88, 226.0),
        'Bi': (83, 208.98)
    }

    # Step 2: Define the 10 candidate molecules, whose names are spelled by element symbols
    MOLECULES = {
        "TUNGSTEN": ['W', 'O', 'La', 'F', 'Ra', 'Am'],
        "BISMUTH": ['Bi', 'Sm', 'U', 'Th'],
        "CARBON": ['C', 'Ar', 'B', 'O', 'N'],
        "YTTRIUM": ['Y', 'Th', 'Ru', 'I', 'U', 'Mn'],
        "HAFNIUM": ['H', 'Al', 'F', 'Ni', 'U', 'Mn'],
        "SILICON": ['Si', 'Li', 'C', 'O', 'N'],
        "PHOSPHORUS": ['P', 'H', 'O', 'S', 'P', 'Ho', 'Ru', 'S'],
        "ARSENIC": ['Ar', 'Se', 'Ni', 'C'],
        "TIN": ['Ti', 'N'],
        "ZINC": ['Zn', 'C']
    }

    results = []

    # Step 3: Calculate energy for each molecule
    for name, atoms in MOLECULES.items():
        n = len(atoms)
        if n == 0:
            continue
        
        masses = np.array([ATOMIC_DATA[atom][1] for atom in atoms])
        
        # Construct Mass-Weighted Barysz Matrix assuming a linear chain
        b_matrix = np.zeros((n, n))
        for i in range(n):
            b_matrix[i, i] = masses[i]**2
            if i < n - 1:
                b_matrix[i, i+1] = masses[i] * masses[i+1]
                b_matrix[i+1, i] = masses[i+1] * masses[i]
        
        eigvals_b = np.linalg.eigvalsh(b_matrix)
        energy = np.sum(np.abs(eigvals_b))
        
        results.append({'name': name, 'energy': energy, 'atoms': atoms})

    # Step 4: Find the molecule with the lowest energy
    min_energy_molecule = min(results, key=lambda x: x['energy'])

    # Step 5: Calculate Moran's I stats for the identified molecule
    name = min_energy_molecule['name']
    energy = min_energy_molecule['energy']
    atoms = min_energy_molecule['atoms']
    n = len(atoms)
    masses = np.array([ATOMIC_DATA[atom][1] for atom in atoms])

    # Construct the mass-weighting matrix W
    w_matrix = np.zeros((n, n))
    for i in range(n-1):
        w_matrix[i, i+1] = masses[i] * masses[i+1]
        w_matrix[i+1, i] = masses[i+1] * masses[i]
    
    s0 = w_matrix.sum()
    
    if s0 == 0:
        i_min, i_max = 0.0, 0.0
    else:
        # Construct the centered weighting matrix Wc
        identity = np.identity(n)
        one_vector = np.ones((n, 1))
        centering_matrix = identity - (one_vector @ one_vector.T) / n
        wc_matrix = centering_matrix @ w_matrix @ centering_matrix
        
        eigvals_wc = np.linalg.eigvalsh(wc_matrix)
        
        i_min = (n / s0) * np.min(eigvals_wc)
        i_max = (n / s0) * np.max(eigvals_wc)

    # Step 6: Print the components of the final calculation
    print(f"The identified element is {name}, spelled by the molecule {'-'.join(atoms)}.")
    print("The components for the final product are:")
    print(f"Mass-Weighted Barysz Graph Energy: {energy}")
    print(f"Minimum Mass-Weighted Moran's I: {i_min}")
    print(f"Maximum Mass-Weighted Moran's I: {i_max}")

    final_product = energy * i_min * i_max
    # This print is for the final answer requested by the prompt format.
    # print(f"Final Answer: {final_product}")

solve_chemistry_puzzle()