import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def solve_chemoinformatics_task():
    """
    Solves the chemoinformatics problem by finding the minimum product of i_max and
    the difference in average chi indices for formaldehyde's homologs.
    """
    # Generate SMILES strings for the aldehyde homologous series (C1 to C20 aldehydes)
    homologs_smiles = ['C=O'] + ['C' * i + 'C=O' for i in range(1, 21)]

    min_product = float('inf')
    best_homolog_info = None

    # Iterate through each homolog
    for smiles in homologs_smiles:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)  # Add hydrogens as they affect descriptors
        num_heavy_atoms = mol.GetNumHeavyAtoms()

        if num_heavy_atoms <= 1:
            continue

        # Step 1 & 2: Calculate Geary autocorrelations (GATSe)
        max_lag = num_heavy_atoms - 1
        descriptor_names = ['GATS{}e'.format(i) for i in range(1, max_lag + 1)]
        calculator = rdMolDescriptors.MolecularDescriptorCalculator(descriptor_names)
        geary_values = calculator.CalcDescriptors(mol)

        if not geary_values:
            continue

        # Step 3: Find the maximum Geary value and its lag (i_max)
        max_geary_value = max(geary_values)
        # The lag 'i' is 1-based, so we add 1 to the 0-based index
        i_max = geary_values.index(max_geary_value) + 1

        # Step 4: Filter homologs where max_geary_value is between 2 and 3
        if 2.0 <= max_geary_value <= 3.0:
            
            # Step 5: Calculate path chi indices (simple and valence)
            simple_chi_indices = []
            valence_chi_indices = []
            
            # We calculate chi indices up to the maximum path length in the molecule
            max_order = num_heavy_atoms - 1
            for n in range(max_order + 1):
                try:
                    simple_chi_indices.append(rdMolDescriptors.CalcChiNn(mol, n))
                    valence_chi_indices.append(rdMolDescriptors.CalcChiNv(mol, n))
                except:
                    # Some higher-order indices may not be defined for smaller molecules
                    break
            
            if not simple_chi_indices:
                continue
            
            # Step 6: Calculate average chi indices and their difference
            avg_simple_chi = sum(simple_chi_indices) / len(simple_chi_indices)
            avg_valence_chi = sum(valence_chi_indices) / len(valence_chi_indices)
            chi_difference = avg_valence_chi - avg_simple_chi
            
            # Step 7: Calculate the final product
            product = i_max * chi_difference
            
            # Step 8: Find the minimum product
            if product < min_product:
                min_product = product
                best_homolog_info = {
                    "smiles": smiles,
                    "i_max": i_max,
                    "chi_diff": chi_difference,
                    "min_prod": min_product,
                }

    # Print the final result and the components of the calculation
    if best_homolog_info:
        smi = best_homolog_info["smiles"]
        im = best_homolog_info["i_max"]
        diff = best_homolog_info["chi_diff"]
        prod = best_homolog_info["min_prod"]
        
        print("The minimum product is derived from the homolog with SMILES string: " + smi)
        print("This homolog has a lag for maximum Geary autocorrelation (i_max) of: " + str(im))
        print(f"Its difference between average valence and simple path chi indices is: {diff:.4f}")
        print(f"The final product calculation is: {im} * {diff:.4f} = {prod:.4f}")
    else:
        print("No homologs matching the specified criteria were found.")

if __name__ == '__main__':
    solve_chemoinformatics_task()