import numpy as np
from rdkit import Chem
from rdkit.Chem import GraphDescriptors

def solve():
    """
    This script finds the minimum product of i_max and a chi index difference
    for formaldehyde homologs that meet a specific Geary autocorrelation criterion.
    """

    # Step 1: Define homologs and constants
    aldehydes = {
        'Formaldehyde': 'C=O', 'Acetaldehyde': 'CC=O', 'Propanal': 'CCC=O',
        'Butanal': 'CCCC=O', 'Pentanal': 'CCCCC=O', 'Hexanal': 'CCCCCC=O',
        'Heptanal': 'CCCCCCC=O', 'Octanal': 'CCCCCCCC=O', 'Nonanal': 'CCCCCCCCC=O',
        'Decanal': 'CCCCCCCCCC=O', 'Undecanal': 'CCCCCCCCCCC=O', 'Dodecanal': 'CCCCCCCCCCCC=O'
    }
    EN = {'C': 2.55, 'O': 3.44, 'H': 2.20}

    # Step 2: Geary autocorrelation calculation function
    def calculate_geary(mol, weights):
        n_atoms = mol.GetNumAtoms()
        if n_atoms <= 1:
            return {}
        
        weights = np.array(weights)
        dist_matrix = Chem.GetDistanceMatrix(mol)
        
        w_mean = np.mean(weights)
        var_sum = np.sum((weights - w_mean)**2)
        if var_sum == 0:
            return {}
        denominator = var_sum / (n_atoms - 1)
        
        geary_results = {}
        max_lag = int(np.max(dist_matrix))
        
        for k in range(1, max_lag + 1):
            indices_i, indices_j = np.where(dist_matrix == k)
            unique_pairs = {(i, j) for i, j in zip(indices_i, indices_j) if i < j}
            pair_count = len(unique_pairs)
            
            if pair_count == 0:
                continue
            
            sum_sq_diff = sum((weights[i] - weights[j])**2 for i, j in unique_pairs)
            numerator = sum_sq_diff / pair_count
            geary_k = numerator / denominator
            geary_results[k] = geary_k
            
        return geary_results

    # Steps 3 & 4: Loop through homologs, filter, and perform calculations
    found_homologs_data = []
    for name, smiles in aldehydes.items():
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        weights = [EN[atom.GetSymbol()] for atom in mol.GetAtoms()]
        
        geary_values = calculate_geary(mol, weights)
        if not geary_values:
            continue
            
        i_max = max(geary_values, key=geary_values.get)
        g_max = geary_values[i_max]
        
        if 2.0 <= g_max <= 3.0:
            num_orders = 11  # Orders 0 to 10
            simple_chis = [GraphDescriptors.ChiNp(mol, n) for n in range(num_orders)]
            valence_chis = [GraphDescriptors.ChiNv(mol, n) for n in range(num_orders)]
            
            avg_simple_chi = np.mean(simple_chis)
            avg_valence_chi = np.mean(valence_chis)
            
            diff_chi = avg_valence_chi - avg_simple_chi
            product = i_max * diff_chi
            
            found_homologs_data.append({
                'name': name,
                'product': product,
                'i_max': i_max,
                'avg_valence_chi': avg_valence_chi,
                'avg_simple_chi': avg_simple_chi
            })

    # Steps 5 & 6: Find minimum product and print results
    if not found_homologs_data:
        print("No homologs found with a maximum Geary autocorrelation between 2.0 and 3.0.")
        print("<<<NO_CANDIDATES_FOUND>>>")
    else:
        min_homolog = min(found_homologs_data, key=lambda x: x['product'])
        
        print(f"Found {len(found_homologs_data)} homolog(s) that satisfy the Geary autocorrelation condition.")
        print(f"The homolog with the minimum product is '{min_homolog['name']}'.")
        print("\nThe product is calculated as: i_max * (average_valence_chi - average_simple_chi)")
        print("\nFinal Equation:")
        print(f"{min_homolog['i_max']} * ({min_homolog['avg_valence_chi']:.4f} - {min_homolog['avg_simple_chi']:.4f}) = {min_homolog['product']:.4f}")
        
        print(f"<<<{min_homolog['product']}>>>")

solve()