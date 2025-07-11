import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def solve_chemoinformatics_problem():
    """
    Finds formaldehyde homologs based on Geary autocorrelation and calculates a minimum product
    of graph theory descriptors for them.
    """
    # Step 1: Define formaldehyde's homologs
    aldehyde_smiles = {
        'Formaldehyde': 'C=O',
        'Acetaldehyde': 'CC=O',
        'Propanal': 'CCC=O',
        'Butanal': 'CCCC=O',
        'Pentanal': 'CCCCC=O',
        'Hexanal': 'CCCCCC=O',
        'Heptanal': 'CCCCCCC=O',
        'Octanal': 'CCCCCCCC=O',
        'Nonanal': 'CCCCCCCCC=O',
        'Decanal': 'CCCCCCCCCC=O'
    }

    # Step 2: Define atomic property (Sanderson electronegativity)
    sanderson_en = {'C': 2.55, 'H': 2.20, 'O': 3.44}

    # List to store results for homologs that meet the criteria
    found_homologs_results = []

    # Step 3: Iterate through each homolog and perform calculations
    for name, smi in aldehyde_smiles.items():
        mol = Chem.MolFromSmiles(smi)
        mol_with_hs = Chem.AddHs(mol)

        # 3a: Create a list of atom weights based on Sanderson electronegativity
        weights = [sanderson_en.get(atom.GetSymbol(), 0) for atom in mol_with_hs.GetAtoms()]

        # 3b: Calculate Geary autocorrelations for lags 0-8
        geary_vals = rdMolDescriptors.CalcAUTOCORR2D(mol_with_hs, atomprop=weights)
        
        # 3c: Find the maximum Geary value and its lag (i_max) for lags > 0
        # We analyze lags from 1 to 8
        if len(geary_vals) > 1:
            geary_vals_from_lag1 = geary_vals[1:]
            max_geary_val = max(geary_vals_from_lag1)
            # The lag is the index in the sublist + 1
            i_max = geary_vals_from_lag1.index(max_geary_val) + 1
        else:
            continue

        # Step 4: Filter homologs where max_geary_val is between 2 and 3
        if 2 < max_geary_val < 3:
            
            # Step 5: Calculate average Chi indices
            orders = range(1, 8)  # Calculate for orders 1 through 7
            
            simple_chi_paths = []
            valence_chi_paths = []
            
            for n in orders:
                try:
                    # RDKit may raise an error if a path of length n does not exist.
                    # In such cases, the chi index is conceptually zero.
                    simple_chi_paths.append(rdMolDescriptors.CalcChiNpath(mol_with_hs, n))
                    valence_chi_paths.append(rdMolDescriptors.CalcChiNvalence(mol_with_hs, n))
                except RuntimeError:
                    simple_chi_paths.append(0.0)
                    valence_chi_paths.append(0.0)

            avg_simple_chi = sum(simple_chi_paths) / len(simple_chi_paths)
            avg_valence_chi = sum(valence_chi_paths) / len(valence_chi_paths)
            diff_chi = avg_valence_chi - avg_simple_chi

            # Step 6: Calculate the final product
            product = i_max * diff_chi

            # Store the detailed results for this homolog
            found_homologs_results.append({
                'name': name,
                'i_max': i_max,
                'max_geary': max_geary_val,
                'avg_simple_chi': avg_simple_chi,
                'avg_valence_chi': avg_valence_chi,
                'diff_chi': diff_chi,
                'product': product
            })

    # Step 7 & 8: Determine the minimum product and output the results
    if not found_homologs_results:
        print("No homologs were found with a maximum Geary autocorrelation value between 2 and 3.")
        return

    print("Found homologs satisfying the Geary autocorrelation condition:")
    for res in found_homologs_results:
        print(f"\n- For {res['name']}:")
        print(f"  - Maximum Geary Autocorrelation = {res['max_geary']:.4f} (at lag i_max = {res['i_max']})")
        print(f"  - Average Simple Path Chi Index = {res['avg_simple_chi']:.4f}")
        print(f"  - Average Valence Path Chi Index = {res['avg_valence_chi']:.4f}")
        print(f"  - Difference in Chi Indices = {res['avg_valence_chi']:.4f} - {res['avg_simple_chi']:.4f} = {res['diff_chi']:.4f}")
        print(f"  - Calculated Product = {res['i_max']} * {res['diff_chi']:.4f} = {res['product']:.4f}")

    # Find the result with the minimum product
    min_result = min(found_homologs_results, key=lambda x: x['product'])
    min_product_val = min_result['product']

    print("\n" + "="*50)
    print(f"The minimum product is for the homolog: {min_result['name']}")
    print(f"Final Minimum Product Value: {min_product_val:.4f}")
    print("="*50)

    # Print the final answer in the required format
    print(f"<<<{min_product_val:.4f}>>>")

# Execute the function
solve_chemoinformatics_problem()