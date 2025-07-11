import numpy as np

def solve_chemistry_puzzle():
    """
    Solves the puzzle by deciphering the elements, calculating graph energies and indices,
    and finding the final product.
    """

    # Step 1: Decipher the Mapping
    # List of 14 common steel alloying elements
    elements_list = ['C', 'Al', 'Si', 'P', 'S', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Mo', 'W']
    
    # Atomic masses for these elements
    atomic_masses = {
        'C': 12.011, 'Al': 26.982, 'Si': 28.085, 'P': 30.974, 'S': 32.06,
        'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938, 'Fe': 55.845,
        'Co': 58.933, 'Ni': 58.693, 'Mo': 95.96, 'W': 183.84
    }

    # Sort elements by atomic number to create the mapping
    # Z: C(6), Al(13), Si(14), P(15), S(16), Ti(22), V(23), Cr(24), Mn(25), Fe(26), Co(27), Ni(28), Mo(42), W(74)
    elements_by_z = sorted(elements_list, key=lambda el: list(atomic_masses.keys()).index(el))
    
    # Mapping from X_i to element symbol
    x_map = {f"X{i+1}": el for i, el in enumerate(elements_by_z)}

    # Step 2: Define Molecules
    y_formulas_str = {
        "Y1": "X1X2X3X4X5X6",
        "Y2": "X8X9X10X11X5X6X12X8X9X8X12X13",
        "Y3": "X11X3X4X14X3X4X5X6",
        "Y4": "X12X4X13X5X6X3",
        "Y5": "X8X9X10X11X9X14X5X6X3",
        "Y6": "X1X10X5X1X9X4X3",
        "Y7": "X8X9X10X11X12X4X5X6",
        "Y8": "X10X2X5X13X9X4X12X4X3",
        "Y9": "X9X14X5X11X3X4X14X3X4X3",
        "Y10": "X1X12X1X3X10X12X13X12X4X3",
    }

    molecules = {}
    for name, formula_str in y_formulas_str.items():
        # Split string like "X1X2X10" into ["X1", "X2", "X10"]
        import re
        x_vars = re.findall(r'X\d+', formula_str)
        molecules[name] = [x_map[x] for x in x_vars]

    energies = {}
    moran_indices = {}

    # Step 3: Calculate Chemical Properties
    for name, atoms in molecules.items():
        n = len(atoms)
        masses = np.array([atomic_masses[atom] for atom in atoms])

        # Calculate Mass-Weighted Moran's I for a clique
        # I_m = -1 / (n - 1)
        if n > 1:
            moran_indices[name] = -1 / (n - 1)
        else:
            moran_indices[name] = 0 # Undefined, but 0 is a neutral choice

        # Calculate Mass-Weighted Barysz Graph Energy
        # Construct the mass-weighted adjacency matrix A*
        sqrt_masses = np.sqrt(masses)
        A_star = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if i != j:
                    A_star[i, j] = 1 / (sqrt_masses[i] * sqrt_masses[j])
        
        # Calculate eigenvalues and sum their absolute values
        eigenvalues = np.linalg.eigvalsh(A_star)
        energies[name] = np.sum(np.abs(eigenvalues))

    # Step 4: Identify Key Values
    e_min = min(energies.values())
    i_m_min = min(moran_indices.values())
    i_m_max = max(moran_indices.values())
    
    # Step 5: Final Calculation
    result = e_min * i_m_min * i_m_max

    # Output the results for clarity
    print(f"Deciphering Assumption: X_i maps to steel elements sorted by atomic number.")
    print(f"Identified lowest Mass-Weighted Barysz Graph Energy (E_min): {e_min}")
    print(f"Identified minimum Mass-Weighted Moran's I (I_m_min): {i_m_min}")
    print(f"Identified maximum Mass-Weighted Moran's I (I_m_max): {i_m_max}")
    print("\nFinal Calculation:")
    print(f"Product = E_min * I_m_min * I_m_max")
    print(f"Product = {e_min} * {i_m_min} * {i_m_max}")
    print(f"Final Answer: {result}")
    print(f"<<<{result}>>>")

solve_chemistry_puzzle()