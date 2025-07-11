import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from mordred import Calculator, descriptors
from mordred.Atom import SandersonElectronegativity
import warnings

def solve_chemoinformatics_task():
    """
    This function implements the step-by-step plan to find the minimum product
    for specified formaldehyde homologs based on Geary autocorrelation and Chi indices.
    """
    # Suppress non-critical warnings from rdkit/mordred for a clean output
    warnings.filterwarnings("ignore", category=UserWarning)
    original_stderr = sys.stderr
    sys.stderr = open('/dev/null', 'w')

    try:
        homologs_found = []
        max_chain_length = 20  # Check aldehydes from C1 to C20
        max_lag = 10           # Maximum lag for Geary autocorrelation
        max_chi_order = 10     # Maximum order for chi indices

        # Create a Mordred calculator instance and register Geary descriptors for each lag
        # using Sanderson Electronegativity as the atomic property.
        calc = Calculator()
        gear_descriptors = [
            descriptors.Autocorrelation.Geary(lag=i, prop=SandersonElectronegativity)
            for i in range(1, max_lag + 1)
        ]
        calc.register(gear_descriptors)

        # Iterate through the homologs of formaldehyde
        for n_carbons in range(1, max_chain_length + 1):
            if n_carbons == 1:
                smiles = "C=O"
                name = "Formaldehyde"
            else:
                smiles = 'C' * (n_carbons - 1) + "C=O"
                name = f"Aldehyde (C{n_carbons})"

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            mol = Chem.AddHs(mol)

            # Calculate Geary descriptors for the molecule
            geary_results = calc(mol).asdict()
            geary_values = np.array([v for v in geary_results.values() if isinstance(v, (float, int))])

            if len(geary_values) == 0:
                continue
            
            # Find the maximum Geary value and the corresponding lag (i_max)
            max_geary_value = np.max(geary_values)
            # np.argmax is 0-indexed; our lags are 1-indexed.
            i_max = np.argmax(geary_values) + 1
            
            # Filter homologs where the max Geary value is between 2.0 and 3.0
            if 2.0 <= max_geary_value <= 3.0:
                # For a found homolog, calculate average Chi indices
                chi_v = [getattr(rdMolDescriptors, f"CalcChi{i}v")(mol) for i in range(max_chi_order + 1)]
                chi_n = [getattr(rdMolDescriptors, f"CalcChi{i}n")(mol) for i in range(max_chi_order + 1)]
                
                avg_chi_v = np.mean(chi_v)
                avg_chi_n = np.mean(chi_n)
                delta_chi = avg_chi_v - avg_chi_n
                
                # Calculate the final product for this homolog
                product = i_max * delta_chi
                
                homologs_found.append({
                    "name": name,
                    "i_max": i_max,
                    "avg_chi_v": avg_chi_v,
                    "avg_chi_n": avg_chi_n,
                    "delta_chi": delta_chi,
                    "product": product
                })
        
        # After checking all homologs, find the one with the minimum product
        if not homologs_found:
            print("No homologs matching the required criteria were found in the specified range.")
            return

        best_homolog = min(homologs_found, key=lambda x: x['product'])
        
        # Output the final result with the calculation breakdown
        print(f"The homolog with the minimum product is {best_homolog['name']}.")
        print("The final calculation is shown below.")
        print("Product = i_max * (average_valence_path_chi - average_simple_path_chi)")
        print(f"Product = {best_homolog['i_max']} * ({best_homolog['avg_chi_v']:.5f} - {best_homolog['avg_chi_n']:.5f})")
        print(f"Product = {best_homolog['i_max']} * {best_homolog['delta_chi']:.5f} = {best_homolog['product']:.5f}")
        
        print(f"\n<<<{best_homolog['product']:.5f}>>>")

    finally:
        # Restore standard error stream and close the file handle
        sys.stderr.close()
        sys.stderr = original_stderr
        warnings.resetwarnings()

# Execute the main function
solve_chemoinformatics_task()