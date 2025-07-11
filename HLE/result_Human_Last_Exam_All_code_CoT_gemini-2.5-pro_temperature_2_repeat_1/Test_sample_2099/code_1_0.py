import numpy as np
import re

def solve_chemistry_puzzle():
    """
    Solves the user's chemistry puzzle by performing graph-based calculations
    on an interpreted list of carbide molecules.
    """

    # Atomic masses of the elements involved in common hard carbides
    atomic_masses = {
        'W': 183.84, 'Mo': 95.96, 'Cr': 51.996, 'V': 50.942, 'Fe': 55.845,
        'Si': 28.085, 'B': 10.81, 'Ti': 47.867, 'Zr': 91.224, 'Hf': 178.49,
        'C': 12.011
    }

    # Interpreted list of molecules based on the puzzle's theme ("weapon preparation").
    molecule_formulas = {
        "Y1": "WC", "Y2": "Mo2C", "Y3": "Cr3C2", "Y4": "VC",
        "Y5": "Fe3C", "Y6": "SiC", "Y7": "B4C", "Y8": "TiC",
        "Y9": "ZrC", "Y10": "HfC"
    }

    def parse_formula(formula_str):
        """Parses a formula like 'Cr3C2' into a list of atom symbols, e.g., ['Cr', 'Cr', 'Cr', 'C', 'C']."""
        parts = re.findall(r'([A-Z][a-z]?)(\d*)', formula_str)
        atom_list = []
        for symbol, count_str in parts:
            count = int(count_str) if count_str else 1
            atom_list.extend([symbol] * count)
        return atom_list

    def calculate_moran_i(d, N, w, w_bar, ss):
        """Calculates Mass-Weighted Moran's I for a given distance `d` in a path graph."""
        if d >= N: return 0.0
        
        # S0 is the sum of weights, which for an unweighted adjacency matrix is the number of edges.
        # For a path graph at distance d, there are 2*(N-d) connections (i->j and j->i).
        S0 = 2.0 * (N - d)
        if S0 == 0: return 0.0

        # Numerator term of Moran's I, simplified for a path graph structure.
        num_cross_product = np.sum([(w[i] - w_bar) * (w[i+d] - w_bar) for i in range(N - d)])
        num_term = 2 * num_cross_product
        
        if ss == 0: return 0.0  # Avoid division by zero if all weights are the same.
        
        return (N / S0) * (num_term / ss)

    results = []

    # Iterate through each molecule to perform calculations
    for name, formula in molecule_formulas.items():
        atom_list = parse_formula(formula)
        N = len(atom_list)
        if N <= 1: continue

        # Vector of atomic masses for the current molecule
        w = np.array([atomic_masses[atom] for atom in atom_list])

        # --- Mass-Weighted Barysz Graph Energy Calculation ---
        # We assume a linear path graph. The distance matrix D is built from this assumption.
        D = np.abs(np.arange(N)[:, None] - np.arange(N))
        
        # The Barysz Matrix BM has atomic masses on the diagonal and 1/distance off-diagonal.
        BM = np.diag(w)
        # Set off-diagonal elements where D_ij is not zero to avoid division by zero.
        off_diag_indices = (D > 0)
        BM[off_diag_indices] = 1.0 / D[off_diag_indices]

        # The energy is the sum of the absolute values of the eigenvalues.
        eigenvalues = np.linalg.eigvalsh(BM) # Use eigvalsh for symmetric matrices
        energy = np.sum(np.abs(eigenvalues))

        # --- Mass-Weighted Moran's I Calculation ---
        # Calculate Moran's I for all possible topological distances d in the path graph.
        w_bar = np.mean(w)
        ss = np.sum((w - w_bar)**2) # Sum of squared deviations from the mean
        d_max = N - 1
        
        moran_values = [calculate_moran_i(d, N, w, w_bar, ss) for d in range(1, d_max + 1)]
            
        min_moran = min(moran_values) if moran_values else 0
        max_moran = max(moran_values) if moran_values else 0

        results.append({
            "name": name, "formula": formula, "energy": energy,
            "min_moran": min_moran, "max_moran": max_moran
        })

    # Identify the molecule Y with the lowest Barysz energy
    lowest_energy_molecule = min(results, key=lambda x: x['energy'])

    # Extract the final values for this molecule and compute the product
    E_final = lowest_energy_molecule['energy']
    min_I_final = lowest_energy_molecule['min_moran']
    max_I_final = lowest_energy_molecule['max_moran']
    final_product = E_final * min_I_final * max_I_final
    
    # Print the required output, including each number in the final equation
    print(f"The identified element Y with the lowest energy is {lowest_energy_molecule['name']} ({lowest_energy_molecule['formula']}).\n")
    print(f"Identified Energy = {E_final}")
    print(f"Minimum Mass-Weighted Moran's I = {min_I_final}")
    print(f"Maximum Mass-Weighted Moran's I = {max_I_final}\n")
    print("The final result is the product of these three values:")
    print(f"{E_final} * {min_I_final} * {max_I_final} = {final_product}")

solve_chemistry_puzzle()
>>>-17.585773173511874