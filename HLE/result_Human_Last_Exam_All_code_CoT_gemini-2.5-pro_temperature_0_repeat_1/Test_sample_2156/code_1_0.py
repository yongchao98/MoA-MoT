import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

def solve_chemistry_problem():
    """
    Solves the entire multi-step computational chemistry problem.
    """

    # --- Step 1: Define Molecules and Constants ---
    # Sanderson Electronegativity values for relevant atoms
    SANDERSON_EN = {1: 2.592, 6: 2.746, 8: 3.654}  # H, C, O

    # Homologous series of aldehydes from C1 to C10
    ALDEHYDE_NAMES = [
        "Formaldehyde", "Acetaldehyde", "Propionaldehyde", "Butyraldehyde",
        "Pentanal", "Hexanal", "Heptanal", "Octanal", "Nonanal", "Decanal"
    ]
    ALDEHYDE_SMILES = [
        "C=O", "CC=O", "CCC=O", "CCCC=O", "CCCCC=O",
        "CCCCCC=O", "CCCCCCC=O", "CCCCCCCC=O", "CCCCCCCCC=O", "CCCCCCCCCC=O"
    ]

    # --- Step 2: Function for Geary Autocorrelation ---
    def calculate_geary_autocorrelations(mol):
        """
        Calculates Geary autocorrelation weighted by Sanderson EN for all lags.
        The formula used is based on the implementation in popular chemoinformatics libraries.
        """
        N = mol.GetNumAtoms()
        if N == 0:
            return {}

        try:
            weights = np.array([SANDERSON_EN[atom.GetAtomicNum()] for atom in mol.GetAtoms()])
        except KeyError as e:
            print(f"Error: Atom with atomic number {e.args[0]} not found in EN dictionary.")
            return {}

        avg_w = np.mean(weights)
        denom = np.sum((weights - avg_w)**2)

        if denom == 0:
            return {}

        dist_matrix = Chem.GetDistanceMatrix(mol)
        max_lag = int(np.max(dist_matrix))
        
        geary_values = {}

        for d in range(1, max_lag + 1):
            numerator_sum = 0
            num_pairs = 0
            for i in range(N):
                for j in range(i + 1, N):  # Sum over unique pairs (i < j)
                    if dist_matrix[i, j] == d:
                        numerator_sum += (weights[i] - weights[j])**2
                        num_pairs += 1
            
            if num_pairs > 0:
                # Formula: G_d = (N - 1) * sum_{i<j,d} (w_i-w_j)^2 / (num_pairs_d * sum_i (w_i-w_bar)^2)
                geary_val = (N - 1) * numerator_sum / (num_pairs * denom)
                geary_values[d] = geary_val
                
        return geary_values

    # --- Step 4: Function for Chi Index Difference ---
    def calculate_avg_chi_difference(mol):
        """
        Calculates the difference between the average valence and simple path chi indices.
        Averages are taken over orders 0 to 10.
        """
        orders = range(11)  # Orders 0 to 10
        valence_chis = []
        simple_chis = []

        for k in orders:
            val_chi_func = getattr(Descriptors, f'Chi{k}v')
            sim_chi_func = getattr(Descriptors, f'Chi{k}n')
            valence_chis.append(val_chi_func(mol))
            simple_chis.append(sim_chi_func(mol))
        
        avg_valence_chi = np.mean(valence_chis)
        avg_simple_chi = np.mean(simple_chis)
        
        return avg_valence_chi - avg_simple_chi

    # --- Main Logic: Steps 3, 5, 6 ---
    min_product = float('inf')
    final_result_info = {}

    for name, smiles in zip(ALDEHYDE_NAMES, ALDEHYDE_SMILES):
        mol = Chem.MolFromSmiles(smiles)
        # Add hydrogens for autocorrelation calculation, as they contribute to distances and properties
        mol_with_hs = Chem.AddHs(mol)

        # Calculate Geary values
        geary_values = calculate_geary_autocorrelations(mol_with_hs)
        if not geary_values:
            continue

        max_geary_val = max(geary_values.values())
        imax = max(geary_values, key=geary_values.get)

        # Step 3: Filter homologs based on the max Geary value
        if 2.0 <= max_geary_val <= 3.0:
            # Step 4: Calculate Chi difference for the found homolog
            # Note: RDKit Descriptors work on the original mol object (hydrogens are implicit)
            avg_chi_diff = calculate_avg_chi_difference(mol)

            # Step 5: Calculate the final product
            product = imax * avg_chi_diff

            # Step 6: Check if this is the new minimum product
            if product < min_product:
                min_product = product
                final_result_info = {
                    "name": name,
                    "product": product,
                    "imax": imax,
                    "chi_diff": avg_chi_diff,
                    "max_geary": max_geary_val
                }

    # --- Output the results ---
    if not final_result_info:
        print("No homologs found matching the Geary autocorrelation criteria (max value between 2 and 3).")
    else:
        imax_val = final_result_info['imax']
        chi_diff_val = final_result_info['chi_diff']
        
        print(f"Found homolog with minimum product: {final_result_info['name']}")
        print(f"Its maximum Geary autocorrelation value is {final_result_info['max_geary']:.4f} (which is between 2 and 3).")
        print("\n--- Calculation of the Minimum Product ---")
        print(f"Lag of maximum Geary value (i_max): {imax_val}")
        print(f"Difference between average chi indices (valence - simple): {chi_diff_val:.4f}")
        print(f"Final Equation: {imax_val} * {chi_diff_val:.4f} = {min_product:.4f}")
        print(f"\nThe minimum product is {min_product:.4f}")
        print(f"<<<{min_product:.4f}>>>")

# Execute the function to solve the problem
solve_chemistry_problem()