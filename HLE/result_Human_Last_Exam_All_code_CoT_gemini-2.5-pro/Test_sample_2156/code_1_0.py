import numpy as np
from rdkit import Chem
from rdkit.Chem import GraphDescriptors

def solve():
    """
    Solves the entire problem as described in the user prompt.
    """

    # Step 1: Define constants and generate homologs
    sanderson_en = {
        'H': 2.592,
        'C': 2.746,
        'O': 3.654,
    }
    # Generate SMILES for aldehydes from C1 (formaldehyde) to C16 (hexadecanal)
    smiles_list = ['C' * k + '=O' for k in range(1, 17)]

    def calculate_geary_autocorrelation(mol, atom_weights):
        """
        Calculates Geary autocorrelation for all lags, finds the max value and its lag.
        """
        num_atoms = mol.GetNumAtoms()
        if num_atoms <= 1:
            return 0, 0

        dist_matrix = Chem.GetDistanceMatrix(mol)
        max_lag = int(np.max(dist_matrix))
        
        weights = np.array([atom_weights.get(atom.GetSymbol()) for atom in mol.GetAtoms()])
        
        mean_weight = np.mean(weights)
        var_sum_sq = np.sum((weights - mean_weight)**2)
        
        if var_sum_sq == 0 or num_atoms <= 1:
            return 0, 0
        
        denominator = var_sum_sq / (num_atoms - 1)
        if denominator == 0:
            return 0, 0
            
        geary_values = {}
        
        for lag in range(1, max_lag + 1):
            sum_sq_diff = 0
            pair_count = 0
            for i in range(num_atoms):
                for j in range(i + 1, num_atoms):
                    if dist_matrix[i, j] == lag:
                        sum_sq_diff += (weights[i] - weights[j])**2
                        pair_count += 1
            
            if pair_count > 0:
                numerator = sum_sq_diff / pair_count
                geary = numerator / denominator
                geary_values[lag] = geary

        if not geary_values:
            return 0, 0
        
        i_max = max(geary_values, key=geary_values.get)
        max_geary = geary_values[i_max]
                
        return max_geary, i_max

    def calculate_avg_chi_diff(mol):
        """
        Calculates the difference between average valence and simple path chi indices.
        """
        dist_matrix = Chem.GetDistanceMatrix(mol)
        max_order = int(np.max(dist_matrix))
        
        # RDKit's GraphDescriptors have functions from Chi0p/v to Chi10p/v
        calc_limit = min(max_order, 10)

        simple_chi_vals = []
        valence_chi_vals = []

        for i in range(calc_limit + 1):
            func_name_p = f"Chi{i}p"
            chi_func_p = getattr(GraphDescriptors, func_name_p)
            simple_chi_vals.append(chi_func_p(mol))
            
            func_name_pv = f"Chi{i}pv"
            chi_func_pv = getattr(GraphDescriptors, func_name_pv)
            valence_chi_vals.append(chi_func_pv(mol))
            
        avg_spc = np.mean(simple_chi_vals) if simple_chi_vals else 0
        avg_vpc = np.mean(valence_chi_vals) if valence_chi_vals else 0
        
        return abs(avg_vpc - avg_spc)

    # Main analysis loop
    valid_homologs_results = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        
        # Step 2: Calculate Geary autocorrelation
        max_geary, i_max = calculate_geary_autocorrelation(mol, sanderson_en)
        
        # Step 3: Filter homologs based on max_geary
        if 2.0 <= max_geary <= 3.0:
            # Step 4: Calculate chi indices difference
            avg_chi_diff = calculate_avg_chi_diff(mol)
            
            # Step 5: Calculate the product
            product = i_max * avg_chi_diff
            
            valid_homologs_results.append({
                'smiles': smiles,
                'i_max': i_max,
                'avg_chi_diff': avg_chi_diff,
                'product': product
            })

    # Step 6: Find the minimum product and print the results
    if not valid_homologs_results:
        print("No homologs found within the specified Geary autocorrelation range [2, 3].")
    else:
        min_result = min(valid_homologs_results, key=lambda x: x['product'])
        
        i_max_val = min_result['i_max']
        chi_diff_val = min_result['avg_chi_diff']
        product_val = min_result['product']

        print("The homolog that gives the minimum product is found.")
        print(f"SMILES: {min_result['smiles']}")
        print(f"The equation for the product is: i_max * |avg_vpc - avg_spc|")
        print(f"The values are:")
        print(f"i_max = {i_max_val}")
        print(f"|avg_vpc - avg_spc| = {chi_diff_val}")
        print(f"Minimum product = {product_val}")
        
        print(f"<<<{product_val}>>>")

solve()