import sys
try:
    from rdkit import Chem
    from rdkit.Chem import GraphDescriptors, rdMolDescriptors
except ImportError:
    print("RDKit is not installed. Please ensure it is installed, for example, using 'pip install rdkit-pypi'")
    sys.exit(1)
import numpy as np

def solve_aldehyde_properties():
    """
    Finds formaldehyde's homologs with a max Geary autocorrelation (weighted by Sanderson electronegativities)
    between 2 and 3. For these homologs, it determines the minimum product of the lag of the max autocorrelation (i_max)
    and the difference between the average valence and simple path chi indices.
    """
    # List of aldehyde homologs (C1 to C12)
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

    pt = Chem.GetPeriodicTable()
    min_product = float('inf')
    best_homolog_info = {}

    for name, smiles in homologs.items():
        # Create molecule with explicit hydrogens
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)

        # Get Sanderson electronegativity for each atom
        en_values = [pt.GetEN(atom.GetAtomicNum()) for atom in mol.GetAtoms()]

        # Calculate Geary Autocorrelation (GATS) for all lags
        dist_matrix = Chem.GetDistanceMatrix(mol)
        max_lag = int(np.max(dist_matrix))

        if max_lag == 0:
            continue
        
        # CalcAUTOCORR2D returns a flat list of descriptors. Geary values start at index 2 * (max_lag + 1).
        all_autocorr = rdMolDescriptors.CalcAUTOCORR2D(mol, prop=en_values, maxlag=max_lag)
        geary_offset = 2 * (max_lag + 1)
        geary_values = list(all_autocorr[geary_offset : geary_offset + max_lag + 1])

        # Find max GATS and i_max for lags >= 1
        gats_values_for_max = geary_values[1:]
        if not gats_values_for_max:
            continue
        gats_max = max(gats_values_for_max)
        i_max = gats_values_for_max.index(gats_max) + 1

        # Filter homologs where max GATS is between 2 and 3
        if 2 <= gats_max <= 3:
            # Calculate average path chi indices (orders 0-10)
            num_chi_orders = 11
            simple_chi_indices = [getattr(GraphDescriptors, f'Chi{i}p')(mol) for i in range(num_chi_orders)]
            valence_chi_indices = [getattr(GraphDescriptors, f'Chi{i}pv')(mol) for i in range(num_chi_orders)]
            
            avg_simple_chi = sum(simple_chi_indices) / num_chi_orders
            avg_valence_chi = sum(valence_chi_indices) / num_chi_orders
            
            # Calculate the difference and the final product
            diff = avg_valence_chi - avg_simple_chi
            product = i_max * diff

            # Check if this product is the new minimum
            if product < min_product:
                min_product = product
                best_homolog_info = {
                    "name": name,
                    "i_max": i_max,
                    "diff": diff,
                    "avg_valence_chi": avg_valence_chi,
                    "avg_simple_chi": avg_simple_chi,
                    "product": product
                }

    # Print the final result
    if not best_homolog_info:
        print("No homologs were found that match the specified criteria.")
    else:
        print("A homolog satisfying the criteria has been identified.")
        print(f"The minimum product was found for: {best_homolog_info['name']}")
        print("\n--- Calculation Details ---")
        print(f"Lag of maximum Geary Autocorrelation (i_max): {best_homolog_info['i_max']}")
        print(f"Average valence path chi index: {best_homolog_info['avg_valence_chi']:.4f}")
        print(f"Average simple path chi index: {best_homolog_info['avg_simple_chi']:.4f}")
        print(f"Difference (valence - simple): {best_homolog_info['diff']:.4f}")
        
        print("\n--- Final Equation ---")
        print(f"Product = i_max * (Difference)")
        print(f"Product = {best_homolog_info['i_max']} * {best_homolog_info['diff']:.4f}")
        print(f"Minimum Product = {best_homolog_info['product']:.4f}")

if __name__ == "__main__":
    solve_aldehyde_properties()
<<<2.3025>>>