import sys
from rdkit import rdBase, Chem
from rdkit.Chem import AllChem
from mordred import Calculator, descriptors

def find_minimum_product():
    """
    Finds formaldehyde's homologs that satisfy a specific Geary autocorrelation criterion
    and then determines the minimum product of i_max and a difference of chi indices.
    """
    # Suppress RDKit warnings for cleaner output
    rdBase.DisableLog('rdApp.warning')

    # 1. Generate homologs: Aldehyde series from C1 to C12
    homologs_data = []
    # Start with formaldehyde (n=1)
    homologs_data.append({"name": "Formaldehyde", "smiles": "C=O"})
    # Generate C2 to C12 aldehydes
    names = [
        "Acetaldehyde", "Propanal", "Butanal", "Pentanal", "Hexanal",
        "Heptanal", "Octanal", "Nonanal", "Decanal", "Undecanal", "Dodecanal"
    ]
    for i, name in enumerate(names, start=2):
        smiles = 'C' * (i-1) + 'C=O'
        homologs_data.append({"name": name, "smiles": smiles})

    # 2. Create a mordred calculator
    # Using the default 2D descriptors is efficient and includes what we need.
    calc = Calculator(descriptors, ignore_3D=True)

    found_homologs = []

    for homolog in homologs_data:
        mol = Chem.MolFromSmiles(homolog['smiles'])
        if not mol:
            continue
        # Descriptors are calculated on molecules with explicit hydrogens
        mol_h = Chem.AddHs(mol)

        # Calculate all descriptors
        all_desc = calc(mol_h)

        # 3. Extract Geary values and find the maximum and corresponding lag (i_max)
        geary_values = []
        for i in range(1, 9):  # Lags from 1 to 8
            desc_name = f'GATS{i}e'
            val = all_desc[desc_name]
            # mordred returns an error object on failure
            if isinstance(val, (float, int)):
                geary_values.append(val)
            else:
                # Stop if a lag cannot be calculated (molecule is too small)
                break
        
        if not geary_values:
            continue

        max_geary_auto = max(geary_values)
        # The lag is the index in the list + 1
        i_max = geary_values.index(max_geary_auto) + 1

        # 4. Check if the max Geary value is within the specified range [2, 3]
        if 2.0 <= max_geary_auto <= 3.0:
            
            # 5. Retrieve the average chi path indices
            avg_simple_chi = all_desc['APC']
            avg_valence_chi = all_desc['AVPC']

            # Ensure chi indices were calculated successfully
            if not isinstance(avg_simple_chi, (float, int)) or not isinstance(avg_valence_chi, (float, int)):
                continue

            # 6. Calculate the difference and the final product
            diff = avg_valence_chi - avg_simple_chi
            product = i_max * diff
            
            found_homologs.append({
                "name": homolog['name'],
                "i_max": i_max,
                "avg_valence_chi": avg_valence_chi,
                "avg_simple_chi": avg_simple_chi,
                "diff": diff,
                "product": product
            })

    if not found_homologs:
        print("No homologs satisfying the criteria were found.")
        return

    # 7. Find the homolog with the minimum product
    min_result_homolog = min(found_homologs, key=lambda x: x['product'])

    # Print the detailed breakdown of the final answer
    print(f"The homolog with the minimum product is {min_result_homolog['name']}.")
    print("\nCalculation details for this homolog:")
    print(f"  - Lag for max Geary Autocorrelation (i_max): {min_result_homolog['i_max']}")
    print(f"  - Average Valence Path Chi Index: {min_result_homolog['avg_valence_chi']:.4f}")
    print(f"  - Average Simple Path Chi Index: {min_result_homolog['avg_simple_chi']:.4f}")
    print(f"  - Difference in Chi Indices: {min_result_homolog['diff']:.4f}")
    
    print("\nThe final product is the result of the following equation:")
    print(f"{min_result_homolog['i_max']} * ({min_result_homolog['avg_valence_chi']:.4f} - {min_result_homolog['avg_simple_chi']:.4f}) = {min_result_homolog['product']:.4f}")
    
    # Final answer in the required format
    print(f"\n<<<{min_result_homolog['product']}>>>")

if __name__ == '__main__':
    # Ensure necessary packages are installed
    try:
        import rdkit
        import mordred
    except ImportError:
        print("This script requires the 'rdkit-pypi' and 'mordred' packages.")
        print("Please install them using: pip install rdkit-pypi mordred")
        sys.exit(1)
        
    find_minimum_product()
