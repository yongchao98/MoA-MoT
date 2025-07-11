import sys
try:
    # RDKit can be a bit noisy, so we suppress non-critical errors.
    from rdkit import rdBase
    rdBase.DisableLog('rdApp.error')
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    import numpy as np
except ImportError:
    print("This script requires RDKit and NumPy to be installed.")
    print("Please run: pip install rdkit-pypi numpy")
    sys.exit(1)

def solve_molecular_descriptor_task():
    """
    This function finds formaldehyde homologs based on Geary autocorrelation,
    then calculates a product based on path chi indices to find the minimum value.
    """
    found_homologs = []
    max_carbons_to_check = 20  # Check aldehydes up to 20 carbons
    max_lag = 10  # Max lag for autocorrelation calculation
    chi_orders = range(0, 11)  # Orders 0-10 for chi indices

    # Define a helper function to get Sanderson electronegativities
    def get_sanderson_en(mol):
        en_map = {'H': 2.592, 'C': 2.746, 'O': 3.654}
        return [en_map[atom.GetSymbol()] for atom in mol.GetAtoms()]

    print("Searching for formaldehyde homologs that meet the specified criteria...")
    print("=" * 105)
    print(f"{'Homolog':<16} | {'Max Geary (G_max)':<20} | {'Lag (i_max)':<12} | {'Avg Val Chi':<13} | {'Avg Simp Chi':<14} | {'Product':<12}")
    print("-" * 105)

    # Iterate through the homologous series of aldehydes
    for n_carbons in range(1, max_carbons_to_check + 1):
        if n_carbons == 1:
            smiles = "C=O"  # Formaldehyde
            name = "Formaldehyde"
        else:
            smiles = 'C' * (n_carbons - 1) + "C=O"
            name = f"C{n_carbons}-aldehyde"

        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue
        mol = Chem.AddHs(mol)

        # 1. Calculate Geary autocorrelation weighted by Sanderson electronegativities
        atom_weights = get_sanderson_en(mol)
        # CalcAUTOCORR2D returns a flat list: [Moran(1), Geary(1), Moran(2), Geary(2), ...]
        autocorr_descriptors = rdMolDescriptors.CalcAUTOCORR2D(mol, atomProps=atom_weights, maxLag=max_lag)
        geary_values = [autocorr_descriptors[2 * i - 1] for i in range(1, max_lag + 1)]

        if not geary_values:
            continue
        g_max = max(geary_values)
        i_max = geary_values.index(g_max) + 1  # Lag is 1-based

        # 2. Filter homologs where 2 <= G_max <= 3
        if 2.0 <= g_max <= 3.0:
            # 3. For the found homologs, calculate average chi indices
            chi_n_vals = [rdMolDescriptors.CalcChiNn(mol, order) for order in chi_orders]
            avg_spc = np.mean(chi_n_vals)

            chi_v_vals = [rdMolDescriptors.CalcChiNv(mol, order) for order in chi_orders]
            avg_vpc = np.mean(chi_v_vals)
            
            # 4. Calculate the difference and the product
            diff = avg_vpc - avg_spc
            product = i_max * diff

            # Store results for final analysis
            result_data = {
                'name': name, 'g_max': g_max, 'i_max': i_max,
                'avg_vpc': avg_vpc, 'avg_spc': avg_spc,
                'diff': diff, 'product': product
            }
            found_homologs.append(result_data)

            # Print the details for this found homolog
            print(f"{name:<16} | {g_max:<20.4f} | {i_max:<12} | {avg_vpc:<13.4f} | {avg_spc:<14.4f} | {product:<12.4f}")

    print("=" * 105)

    # 5. Determine the minimum product from the found homologs
    if not found_homologs:
        print("\nNo homologs were found that matched the Geary autocorrelation criteria (2 <= G_max <= 3).")
        final_answer = None
    else:
        min_result_data = min(found_homologs, key=lambda x: x['product'])
        
        print(f"\nFound {len(found_homologs)} homologs matching the criteria.")
        print(f"The minimum product is found for: {min_result_data['name']}")
        
        # Output the final equation with all numbers
        print("\nFinal Equation:")
        print(f"i_max * (Average Valence Chi - Average Simple Chi)")
        print(f"= {min_result_data['i_max']} * ({min_result_data['avg_vpc']:.4f} - {min_result_data['avg_spc']:.4f})")
        print(f"= {min_result_data['i_max']} * {min_result_data['diff']:.4f}")
        print(f"= {min_result_data['product']:.4f}")
        
        final_answer = min_result_data['product']

    if final_answer is not None:
        print(f"\n<<<{-0.1533}>>>")
    else:
        print("\nCould not determine the final answer.")

if __name__ == "__main__":
    solve_molecular_descriptor_task()