import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import numpy as np

def solve_chemoinformatics_task():
    """
    This function performs the entire task as requested.
    1. It identifies homologs of formaldehyde (aliphatic aldehydes).
    2. For each homolog, it calculates the Geary autocorrelation weighted by Sanderson electronegativities for all possible lags.
    3. It finds the maximum Geary value (G_max) and the corresponding lag (i_max).
    4. It filters for homologs where 2 <= G_max <= 3.
    5. For these found homologs, it calculates the average simple path chi index and the average valence path chi index.
    6. It computes the difference between these two averages (delta_chi).
    7. It then calculates the product of i_max and delta_chi.
    8. Finally, it determines the minimum product among all found homologs and prints the result.
    """

    def calculate_geary_autocorrelation(mol):
        """Calculates Geary autocorrelation for all lags for a given molecule."""
        en_map = {'C': 2.55, 'O': 3.44, 'H': 2.20}
        
        mol_h = Chem.AddHs(mol)
        num_atoms = mol_h.GetNumAtoms()

        if num_atoms <= 1:
            return {}

        weights = np.array([en_map.get(atom.GetSymbol()) for atom in mol_h.GetAtoms()])
        avg_weight = np.mean(weights)
        
        denom_term_sum_sq_dev = np.sum((weights - avg_weight)**2)

        if denom_term_sum_sq_dev == 0:
            return {}
            
        dist_matrix = Chem.GetDistanceMatrix(mol_h)
        max_lag = int(np.max(dist_matrix))
        geary_correlations = {}

        for lag in range(1, max_lag + 1):
            indices = np.where(np.triu(dist_matrix) == lag)
            num_pairs = len(indices[0])
            
            if num_pairs > 0:
                numerator_sum_sq_diff = 0
                for i, j in zip(indices[0], indices[1]):
                    numerator_sum_sq_diff += (weights[i] - weights[j])**2
                
                denominator = 2 * num_pairs * denom_term_sum_sq_dev
                if denominator > 0:
                    geary_val = ((num_atoms - 1) * numerator_sum_sq_diff) / denominator
                    geary_correlations[lag] = geary_val
        
        return geary_correlations

    def calculate_chi_indices_and_diff(mol):
        """Calculates average path chi indices and their difference."""
        num_heavy_atoms = mol.GetNumAtoms()
        if num_heavy_atoms <= 1:
            return 0, 0, 0
            
        max_order = num_heavy_atoms - 1
        
        p_vals = []
        pv_vals = []
        
        # We average over all path orders from 1 to the maximum possible path length
        for i in range(1, max_order + 1):
            # Type 0: simple path chi index
            p_val = rdMolDescriptors.CalcChi(mol, i, 0)
            # Type 1: valence path chi index
            pv_val = rdMolDescriptors.CalcChi(mol, i, 1)
            p_vals.append(p_val)
            pv_vals.append(pv_val)

        avg_chi_p = np.mean(p_vals) if p_vals else 0
        avg_chi_pv = np.mean(pv_vals) if pv_vals else 0
        # The difference of averages is the same as the average of differences
        avg_diff = avg_chi_pv - avg_chi_p
        
        return avg_chi_p, avg_chi_pv, avg_diff

    # --- Main Execution ---
    aldehyde_smiles_list = [
        'C=O', 'CC=O', 'CCC=O', 'CCCC=O', 'CCCCC=O',
        'CCCCCC=O', 'CCCCCCC=O', 'CCCCCCCC=O', 'CCCCCCCCC=O', 'CCCCCCCCCC=O'
    ]

    found_homologs_data = []

    for smiles in aldehyde_smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        
        geary_correlations = calculate_geary_autocorrelation(mol)
        
        if not geary_correlations:
            continue
        
        i_max = max(geary_correlations, key=geary_correlations.get)
        g_max = geary_correlations[i_max]

        if 2 <= g_max <= 3:
            avg_chi_p, avg_chi_pv, delta_chi = calculate_chi_indices_and_diff(mol)
            product = i_max * delta_chi
            
            found_homologs_data.append({
                'smiles': smiles,
                'product': product,
                'i_max': i_max,
                'avg_chi_p': avg_chi_p,
                'avg_chi_pv': avg_chi_pv,
                'delta_chi': delta_chi
            })

    if not found_homologs_data:
        print("No homologs were found that match the specified criteria.")
        return

    min_result = min(found_homologs_data, key=lambda x: x['product'])
    
    final_product = min_result['product']
    final_i_max = min_result['i_max']
    final_avg_chi_p = min_result['avg_chi_p']
    final_avg_chi_pv = min_result['avg_chi_pv']
    final_delta_chi = min_result['delta_chi']

    print(f"The minimum product is found for the homolog with SMILES representation: {min_result['smiles']}")
    print(f"The lag with maximum Geary autocorrelation (i_max) is: {final_i_max}")
    print(f"The average valence path chi index for this homolog is: {final_avg_chi_pv}")
    print(f"The average simple path chi index for this homolog is: {final_avg_chi_p}")
    print(f"The difference between the average chi indices (average valence - average simple) is: {final_delta_chi}")
    print("\nThe final equation is:")
    print(f"{final_i_max} * ({final_avg_chi_pv} - {final_avg_chi_p}) = {final_product}")
    
    # Final answer in the required format
    print(f"<<<{final_product}>>>")

# Execute the main function to solve the problem
solve_chemoinformatics_task()