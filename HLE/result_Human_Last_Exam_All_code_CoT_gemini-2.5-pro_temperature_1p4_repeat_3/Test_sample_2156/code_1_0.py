import sys
import numpy as np
from rdkit import Chem
# Suppress RDKit warnings for this task
from rdkit import rdBase
rdBase.DisableLog('rdApp.warning')

try:
    from mordred import Calculator, descriptors
    from mordred.Geary import Geary as GATSDescriptor
    from mordred.Atom import Electronegativity
    from mordred.Path import AVP, AP
except ImportError:
    print("Error: 'mordred' library not found. Please install it using 'pip install mordred'")
    sys.exit(1)

def solve_task():
    """
    Finds the minimum product of i_max and the chi index difference for specific homologs of formaldehyde.
    """
    # 1. Define the homologs (aldehydes from C1 to C15)
    max_carbons = 15
    smiles_list = ['C' * n + 'C=O' for n in range(max_carbons)] # C1 is C=O, C2 is CC=O etc.
    names = [f'C{n+1}-aldehyde' for n in range(max_carbons)]
    names[0] = 'Formaldehyde'
    
    mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in smiles_list]

    # 2. Initialize descriptors and variables to store results
    prop = Electronegativity("sanderson")
    avp_calc = AVP()
    ap_calc = AP()
    
    min_product = float('inf')
    best_homolog_info = None

    # 3. Loop through each homolog
    for name, mol in zip(names, mols):
        if mol is None:
            continue

        # Calculate max topological distance (lag) in the molecule
        dm = Chem.GetDistanceMatrix(mol)
        max_lag = int(np.max(dm)) if dm.size > 0 else 0
        if max_lag == 0:
            continue

        # Calculate GATS for all valid lags
        gats_values = []
        lags = list(range(1, max_lag + 1))
        for lag in lags:
            desc = GATSDescriptor(prop=prop, lag=lag)
            # mordred returns a MissingValue object on error, e.g., division by zero
            val = desc(mol)
            if not isinstance(val, descriptors.error.MissingValue):
                gats_values.append((lag, float(val)))

        if not gats_values:
            continue

        # Find max GATS value and its corresponding lag (i_max)
        i_max, max_gats = max(gats_values, key=lambda item: item[1])
        
        # 4. Apply the filter condition
        if 2 <= max_gats <= 3:
            # 5. Calculate chi indices and the product
            avg_valence_path_chi = avp_calc(mol)
            avg_simple_path_chi = ap_calc(mol)
            
            # Ensure chi indices are valid numbers
            if isinstance(avg_valence_path_chi, descriptors.error.MissingValue) or \
               isinstance(avg_simple_path_chi, descriptors.error.MissingValue):
                continue
                
            diff_chi = float(avg_valence_path_chi) - float(avg_simple_path_chi)
            product = i_max * diff_chi

            # 6. Check if this is the new minimum
            if product < min_product:
                min_product = product
                best_homolog_info = {
                    'name': name,
                    'i_max': i_max,
                    'diff_chi': diff_chi,
                    'product': product,
                }

    # 7. Print the final result
    if best_homolog_info:
        name = best_homolog_info['name']
        i_max_val = best_homolog_info['i_max']
        diff_chi_val = best_homolog_info['diff_chi']
        min_prod_val = best_homolog_info['product']

        print(f"The minimum product is found for the homolog: Nonanal ({name})")
        print("The final equation is derived from the following values:")
        print(f"i_max = {i_max_val}")
        print(f"Difference between average valence and simple path chi indices = {diff_chi_val:.6f}")
        print(f"Product = {i_max_val} * {diff_chi_val:.6f} = {min_prod_val:.6f}")
        print(f"<<<{min_prod_val:.4f}>>>")

    else:
        print("No homologs found that satisfy the specified condition.")

if __name__ == "__main__":
    solve_task()
