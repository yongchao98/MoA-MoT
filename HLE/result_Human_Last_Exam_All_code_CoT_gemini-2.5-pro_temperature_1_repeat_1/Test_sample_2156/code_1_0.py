import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, GraphDescriptors

def solve_chemtask():
    """
    Finds the minimum product of i_max and a chi index difference
    for specific homologs of formaldehyde.
    """
    min_product = float('inf')
    final_equation_components = None

    # We will test aldehydes with 1 to 20 carbon atoms.
    for n_carbons in range(1, 21):
        if n_carbons == 1:
            smiles = 'C=O' # Formaldehyde
        else:
            # SMILES for other straight-chain aldehydes (e.g., CCC=O for propanal)
            smiles = 'C' * (n_carbons - 1) + 'C=O'

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue
        # Add hydrogens, as they are required for these descriptor calculations.
        mol = Chem.AddHs(mol)

        # Step 2: Calculate 2D autocorrelation descriptors.
        all_autocorr = rdMolDescriptors.CalcAUTOCORR2D(mol)
        
        # GATS weighted by Sanderson EN for lags 1-8 are at indices 80-87.
        gats_en_values = all_autocorr[80:88]

        # Step 3: Find the maximum GATS value and its corresponding lag (i_max).
        max_gats = -float('inf')
        i_max = -1
        for i, gats_val in enumerate(gats_en_values):
            lag = i + 1  # Lags are 1-based, indices are 0-based.
            if gats_val > max_gats:
                max_gats = gats_val
                i_max = lag
        
        # Filter homologs where 2 <= max_gats <= 3.
        if 2 <= max_gats <= 3:
            # Step 4: Calculate the difference in average chi indices.
            # We interpret "average ... chi index" as the first-order average chi index.
            avg_chi_v = GraphDescriptors.GetAverageValenceIpc(mol)
            avg_chi_s = GraphDescriptors.GetAverageConnectivityIpc(mol)
            diff = avg_chi_v - avg_chi_s

            # Step 5: Calculate the product.
            product = i_max * diff
            
            # Step 6: Keep track of the minimum product found.
            if product < min_product:
                min_product = product
                final_equation_components = {
                    'i_max': i_max,
                    'diff': diff,
                    'product': product,
                }

    # Step 7 & 8: Output the final result.
    if final_equation_components:
        i_max_final = final_equation_components['i_max']
        diff_final = final_equation_components['diff']
        min_prod_final = final_equation_components['product']
        
        print("The final calculation uses the following numbers:")
        print(f"i_max: {i_max_final}")
        print(f"Difference between average chi indices: {diff_final}")
        print("\nFinal Equation:")
        print(f"{i_max_final} * {diff_final} = {min_prod_final}")
    else:
        print("No homologs were found that met the specified criteria.")

if __name__ == "__main__":
    # Ensure RDKit is installed
    try:
        from rdkit import Chem
    except ImportError:
        print("RDKit is not installed. Please install it using 'pip install rdkit-pypi'", file=sys.stderr)
        sys.exit(1)
        
    solve_chemtask()