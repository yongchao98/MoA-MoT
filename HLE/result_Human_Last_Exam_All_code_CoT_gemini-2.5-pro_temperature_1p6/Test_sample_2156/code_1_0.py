import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def solve_molecular_properties():
    """
    Finds formaldehyde homologs based on Geary autocorrelation and calculates a specific property product.
    
    The function identifies n-aldehydes where the maximum Geary autocorrelation (weighted by
    Sanderson electronegativity) is between 2 and 3. For these homologs, it computes the
    product of the lag at which the maximum occurs (i_max) and the average difference
    between valence and simple path chi indices. The minimum of these products is then reported.
    """

    # Sanderson electronegativity values (J. Am. Chem. Soc. 1989, 111, 5, 1731–1735)
    sanderson_en = {
        1: 2.592,  # H
        6: 2.746,  # C
        8: 3.654,  # O
    }

    # Helper function to generate SMILES strings for n-aldehydes
    def get_aldehyde_smiles(n_carbons):
        if n_carbons == 1:
            return "C=O"  # Formaldehyde
        return "C" * (n_carbons - 1) + "C=O"

    # Helper function to get the list of EN values for atoms in a molecule
    def get_en_list(mol):
        return [sanderson_en.get(atom.GetAtomicNum(), 0) for atom in mol.GetAtoms()]

    found_homologs = []

    # Iterate through n-aldehydes from C1 to C15 to find candidates
    for n in range(1, 16):
        smiles = get_aldehyde_smiles(n)
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue
        mol_h = Chem.AddHs(mol)

        # Get atom properties (Sanderson EN) for the molecule
        en_values = get_en_list(mol_h)
        
        # Determine the maximum possible lag (graph diameter)
        try:
            max_lag = int(np.max(Chem.GetDistanceMatrix(mol_h)))
        except:
            continue # Skip if distance matrix fails

        geary_values = []
        for i in range(1, max_lag + 1):
            try:
                # Calculate Geary autocorrelation for lag i
                geary = rdMolDescriptors.CalcGearyAutocorrelation(mol_h, lag=i, prop=en_values)
                geary_values.append(geary)
            except:
                geary_values.append(np.nan)

        # Find the maximum Geary value and its corresponding lag
        clean_geary_values = [g for g in geary_values if not np.isnan(g)]
        if not clean_geary_values:
            continue
        
        max_geary = max(clean_geary_values)
        i_max = geary_values.index(max_geary) + 1  # Lag is 1-indexed

        # Check if the homolog meets the specified Geary autocorrelation criteria
        if 2 <= max_geary <= 3:
            # If criteria met, calculate average chi index difference
            chi_diffs = []
            # Use orders from 1 to 8 for the average calculation
            for order in range(1, 9):
                try:
                    chi_v = rdMolDescriptors.CalcChiNv(mol_h, order)
                    chi_n = rdMolDescriptors.CalcChiNn(mol_h, order)
                    chi_diffs.append(chi_v - chi_n)
                except:
                    # This order might not have corresponding paths in the molecule
                    pass
            
            avg_chi_diff = np.mean(chi_diffs) if chi_diffs else 0
            
            # Calculate the final product for this homolog
            product = i_max * avg_chi_diff
            
            found_homologs.append({
                "name": f"{'n' * (n > 3)}{['formaldehyde', 'acetaldehyde', 'propionaldehyde', 'butyraldehyde'][n-1] if n <= 4 else f'C{n}-aldehyde'} ({smiles})",
                "i_max": i_max,
                "max_geary": max_geary,
                "avg_chi_diff": avg_chi_diff,
                "product": product
            })

    # Identify the homolog with the minimum product
    if not found_homologs:
        print("No homologs were found that satisfy the Geary autocorrelation condition.")
    else:
        min_homolog = min(found_homologs, key=lambda x: x["product"])

        print(f"The minimum product was found for the molecule: {min_homolog['name']}.")
        print("This molecule meets the criteria with the following properties:")
        print(f"  - Maximum Geary Autocorrelation = {min_homolog['max_geary']:.4f}")
        print(f"  - Lag at Maximum (i_max) = {min_homolog['i_max']}")
        print(f"  - Average Chi Index Difference = {min_homolog['avg_chi_diff']:.4f}")
        
        print("\nThe final equation to determine the minimum product is:")
        print(f"{min_homolog['i_max']} * {min_homolog['avg_chi_diff']:.4f} = {min_homolog['product']:.4f}")
        
        print(f"\n<<<{min_homolog['product']}>>>")

if __name__ == "__main__":
    solve_molecular_properties()