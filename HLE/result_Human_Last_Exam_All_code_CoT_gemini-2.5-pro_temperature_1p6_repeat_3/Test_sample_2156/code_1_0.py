import numpy as np
from rdkit import Chem
from mordred import Calculator, descriptors

def solve_homolog_problem():
    """
    Finds formaldehyde's homologs that meet specific descriptor criteria and 
    determines the minimum product of i_max and a difference of chi indices.
    """
    # 1. Define the list of formaldehyde homologs (aldehydes)
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
        "Decanal": "CCCCCCCCCC=O",
        "Undecanal": "CCCCCCCCCCC=O",
        "Dodecanal": "CCCCCCCCCCCC=O",
    }

    # 2. Create a Mordred calculator for the required descriptors
    calc = Calculator(descriptors, ignore_3D=True)
    # Geary autocorrelation with Sanderson electronegativity (lags 1 to 8)
    calc.register(descriptors.geary.Geary(prop='se', lag=range(1, 9)))
    # Valence and Simple path chi indices (orders 1 to 10)
    calc.register(descriptors.path_count.PathCount(type='valence', order=range(1, 11)))
    calc.register(descriptors.path_count.PathCount(type='path', order=range(1, 11)))

    found_homologs_results = []

    # Iterate through each homolog to perform calculations
    for name, smiles in homologs.items():
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue
        mol = Chem.AddHs(mol)

        # Calculate all registered descriptors for the molecule
        mol_descriptors = calc(mol)

        # 3. Extract GATS values and find the maximum
        gats_values = {}
        for i in range(1, 9):
            desc_name = f'GATS{i}se'
            value = mol_descriptors[desc_name]
            # Mordred may return non-float error values
            if isinstance(value, float) and np.isfinite(value):
                gats_values[i] = value
        
        if not gats_values:
            continue

        i_max = max(gats_values, key=gats_values.get)
        max_gats = gats_values[i_max]

        # 4. Filter homologs where 2 <= max_gats <= 3
        if 2 <= max_gats <= 3:
            # 5. Calculate average WPC and SPC
            wpc_values = [mol_descriptors[f'WPC{i:02d}'] for i in range(1, 11) if isinstance(mol_descriptors[f'WPC{i:02d}'], float)]
            spc_values = [mol_descriptors[f'SPC{i:02d}'] for i in range(1, 11) if isinstance(mol_descriptors[f'SPC{i:02d}'], float)]

            if not wpc_values or not spc_values:
                continue

            avg_wpc = np.mean(wpc_values)
            avg_spc = np.mean(spc_values)
            
            # 6. Compute the product
            diff = avg_wpc - avg_spc
            product = i_max * diff
            
            found_homologs_results.append({
                "name": name,
                "i_max": i_max,
                "diff": diff,
                "product": product,
            })

    # 7. Determine the minimum product
    if not found_homologs_results:
        print("No homologs were found that met the specified criteria.")
        return

    min_result = min(found_homologs_results, key=lambda x: x['product'])

    print(f"A total of {len(found_homologs_results)} homologs were found matching the criteria.")
    print(f"The homolog with the minimum product is: {min_result['name']}")
    print("\nFinal Equation:")
    print(f"Product = i_max * (Average WPC - Average SPC)")
    print(f"Product = {min_result['i_max']} * ({min_result['diff']:.6f})")
    print(f"Minimum Product = {min_result['product']:.6f}")
    
    # Print the final answer in the specified format
    print(f"\n<<<{min_result['product']}>>>")


if __name__ == "__main__":
    solve_homolog_problem()