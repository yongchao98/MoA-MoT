import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def find_minimum_product():
    """
    Solves the main problem by finding the minimum product of i_max and the chi index difference
    for formaldehyde homologs satisfying the given conditions.
    """
    # Sanderson electronegativities for relevant atoms
    sanderson_en = {1: 2.20, 6: 2.55, 8: 3.44}

    def calculate_gats(mol):
        """Calculates the maximum Geary autocorrelation and the corresponding lag."""
        mol_h = Chem.AddHs(mol)
        n_atoms = mol_h.GetNumAtoms()
        if n_atoms < 2:
            return None, None

        weights = [sanderson_en.get(atom.GetAtomicNum(), 0) for atom in mol_h.GetAtoms()]
        w_bar = sum(weights) / n_atoms
        
        denominator_term = sum([(w - w_bar) ** 2 for w in weights])
        if denominator_term == 0: return None, None
        denominator = denominator_term / (n_atoms - 1)
        if denominator == 0: return None, None

        dist_matrix = Chem.GetDistanceMatrix(mol_h)
        max_lag = int(np.max(dist_matrix))
        gats_values = {}

        for lag in range(1, max_lag + 1):
            sum_sq_diff = 0.0
            delta_count = 0
            for i in range(n_atoms):
                for j in range(i + 1, n_atoms):
                    if dist_matrix[i, j] == lag:
                        sum_sq_diff += (weights[i] - weights[j]) ** 2
                        delta_count += 1
            
            if delta_count > 0:
                numerator = sum_sq_diff / (2 * delta_count)
                gats = numerator / denominator
                gats_values[lag] = gats

        if not gats_values:
            return None, None

        i_max = max(gats_values, key=gats_values.get)
        max_gats = gats_values[i_max]
        return max_gats, i_max

    def calculate_chi_components(mol):
        """Calculates average valence and simple path chi indices."""
        mol_no_hs = Chem.RemoveHs(mol)
        num_heavy_atoms = mol_no_hs.GetNumAtoms()
        
        vpc_values = []
        spc_values = []
        
        # Calculate chi indices for path orders k=1, 2, ...
        # The average is taken over all orders for which paths exist.
        for k in range(1, num_heavy_atoms):
            try:
                vpc = rdMolDescriptors.CalcChiNpathv(mol_no_hs, k)
                spc = rdMolDescriptors.CalcChiNpath(mol_no_hs, k)
                # Ensure a path of this order exists by checking for non-zero sum,
                # as some RDKit versions might return 0 for non-existent paths.
                if vpc > 0 or spc > 0:
                    vpc_values.append(vpc)
                    spc_values.append(spc)
                else:
                    break
            except RuntimeError:
                break # Stop if path of length k doesn't exist
        
        if not vpc_values:
            return 0, 0, 0

        avg_vpc = sum(vpc_values) / len(vpc_values)
        avg_spc = sum(spc_values) / len(spc_values)
        diff = avg_vpc - avg_spc
        return diff, avg_vpc, avg_spc

    # A list of aldehydes (formaldehyde homologs and their isomers) to test
    aldehydes = {
        "Formaldehyde": "C=O",
        "Acetaldehyde": "CC=O",
        "Propanal": "CCC=O",
        "Butanal": "CCCC=O",
        "Isobutyraldehyde": "CC(C)C=O",
        "Pentanal": "CCCCC=O",
        "2,2-Dimethylpropanal": "C(C)(C)C=O",
        "Hexanal": "CCCCCC=O",
        "2,2-Dimethylbutanal": "CCC(C)(C)C=O",
        "2,3-Dimethylbutanal": "CC(C)C(C)C=O",
    }

    min_product = float('inf')
    best_result = {}

    for name, smiles in aldehydes.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        max_gats, i_max = calculate_gats(mol)

        if max_gats is not None and 2 <= max_gats <= 3:
            chi_diff, avg_vpc, avg_spc = calculate_chi_components(mol)
            product = i_max * chi_diff

            if product < min_product:
                min_product = product
                best_result = {
                    "name": name,
                    "i_max": i_max,
                    "avg_vpc": avg_vpc,
                    "avg_spc": avg_spc,
                    "product": product,
                }

    if not best_result:
        print("No homolog found that satisfies the Geary autocorrelation condition.")
    else:
        # Print the final equation with each number, as requested
        i_max_val = best_result['i_max']
        avg_vpc_val = best_result['avg_vpc']
        avg_spc_val = best_result['avg_spc']
        product_val = best_result['product']
        
        print(f"The minimum product is found for {best_result['name']}.")
        print("The final calculation is:")
        print(f"{i_max_val} * ({avg_vpc_val:.4f} - {avg_spc_val:.4f}) = {product_val:.4f}")

find_minimum_product()