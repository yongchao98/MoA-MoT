import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def solve_homolog_problem():
    """
    This function finds the minimum product of i_max and the difference
    between average valence and simple path chi indices for formaldehyde's homologs
    that satisfy a specific Geary autocorrelation condition.
    """
    # Step 1: Define homologs and Sanderson electronegativities
    homologs = {
        "Formaldehyde": "C=O",
        "Acetaldehyde": "CC=O",
        "Propionaldehyde": "CCC=O",
        "Butyraldehyde": "CCCC=O",
        "Pentanal": "CCCCC=O",
        "Hexanal": "CCCCCC=O",
        "Heptanal": "CCCCCCC=O",
        "Octanal": "CCCCCCCC=O",
        "Nonanal": "CCCCCCCCC=O",
        "Decanal": "CCCCCCCCCC=O"
    }
    sanderson_en = {'C': 2.746, 'H': 2.592, 'O': 3.654}

    found_homologs_results = []
    max_chi_order = 10

    # Iterate through each homolog to perform calculations
    for name, smiles in homologs.items():
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)

        # Step 2: Calculate Geary Autocorrelation
        props = [sanderson_en[atom.GetSymbol()] for atom in mol.GetAtoms()]
        geary_values = rdMolDescriptors.CalcGearyAutoCorr(mol, props)
        
        # Consider lags from 1 to 8 (index 0 is lag 0)
        geary_values_lag1_plus = list(geary_values[1:])
        if not geary_values_lag1_plus:
            continue
            
        geary_max = max(geary_values_lag1_plus)
        # The lag 'i' is the index in the original tuple, which is index_in_slice + 1
        i_max = geary_values_lag1_plus.index(geary_max) + 1

        # Step 3: Filter homologs based on the condition
        if 2 <= geary_max <= 3:
            # Step 4: Calculate Chi Indices and Product
            chi_path = [rdMolDescriptors.CalcChiPath(mol, i) for i in range(max_chi_order + 1)]
            chi_v_path = [rdMolDescriptors.CalcChiVPath(mol, i) for i in range(max_chi_order + 1)]

            aspc = np.mean(chi_path)
            adpc = np.mean(chi_v_path)
            
            diff = adpc - aspc
            product = i_max * diff
            
            found_homologs_results.append({
                "name": name,
                "i_max": i_max,
                "aspc": aspc,
                "adpc": adpc,
                "product": product
            })

    # Step 5: Find and print the minimum product
    if not found_homologs_results:
        print("No homologs were found that satisfy the Geary autocorrelation condition.")
    else:
        min_result = min(found_homologs_results, key=lambda x: x['product'])
        
        name = min_result['name']
        i_max_val = min_result['i_max']
        adpc_val = min_result['adpc']
        aspc_val = min_result['aspc']
        product_val = min_result['product']

        print(f"The homolog with the minimum product is {name}.")
        print("\nCalculation details:")
        print(f"Maximum Geary autocorrelation lag (i_max): {i_max_val}")
        print(f"Average valence path chi index (adpc): {adpc_val}")
        print(f"Average simple path chi index (aspc): {aspc_val}")
        print("\nFinal Equation:")
        print(f"{i_max_val} * ({adpc_val} - {aspc_val}) = {product_val}")
        
        print(f"\n<<<{product_val}>>>")

solve_homolog_problem()