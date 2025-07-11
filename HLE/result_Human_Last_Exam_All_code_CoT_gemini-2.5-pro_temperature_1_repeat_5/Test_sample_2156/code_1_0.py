import sys
import numpy as np
from rdkit import Chem
from mordred import Calculator, descriptors

def find_minimum_product_of_homologs():
    """
    This script solves the user's request by following these steps:

    1.  It defines a list of formaldehyde homologs (aldehydes from C1 to C12).
    2.  It initializes the `mordred` calculator with the required descriptors:
        - Geary autocorrelation weighted by Sanderson electronegativities.
        - Path Chi indices (both average valence and average simple).
    3.  It then iterates through each homolog:
        a. Calculates all required descriptors for the molecule.
        b. Finds the maximum Geary autocorrelation value and its corresponding lag (`i_max`).
        c. Checks if this maximum value is between 2.0 and 3.0.
        d. If the condition is met, it calculates the difference between the sum of
           average valence path chi indices and the sum of average simple path chi indices (`delta_chi`).
        e. It computes the product of `i_max` and `delta_chi`.
        f. It keeps track of the homolog that produces the minimum product.
    4.  Finally, it prints the details of the calculation for the homolog with the minimum
        product and the final answer in the specified format.
    """
    try:
        # 1. Define homologs
        homologs = {
            'formaldehyde': 'C=O',
            'acetaldehyde': 'CC=O',
            'propanal':     'CCC=O',
            'butanal':      'CCCC=O',
            'pentanal':     'CCCCC=O',
            'hexanal':      'CCCCCC=O',
            'heptanal':     'CCCCCCC=O',
            'octanal':      'CCCCCCCC=O',
            'nonanal':      'CCCCCCCCC=O',
            'decanal':      'CCCCCCCCCC=O',
            'undecanal':    'CCCCCCCCCCC=O',
            'dodecanal':    'CCCCCCCCCCCC=O'
        }

        # 2. Initialize mordred calculator
        calc = Calculator(descriptors, ignore_3D=True)
        calc.register([descriptors.Geary, descriptors.PathChi])

        min_product = float('inf')
        result_details = {}

        # 3. Iterate through homologs
        for name, smiles in homologs.items():
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            mol_h = Chem.AddHs(mol)

            d_vals = calc(mol_h)

            # Extract Geary values (weighted by Sanderson electronegativity 'pe')
            gats_values = []
            for i in range(1, 9): # Lags from 1 to 8
                val = d_vals[f'GATS{i}pe']
                gats_values.append(float(val) if val is not None and not np.isnan(val) else -1e9)

            if not gats_values or max(gats_values) <= -1e9:
                continue

            # Find max GATS and i_max
            max_gats = max(gats_values)
            i_max = gats_values.index(max_gats) + 1 # Lag is index + 1

            # Filter by GATS value
            if 2.0 <= max_gats <= 3.0:
                # Calculate Path Chi difference
                avp_chi_sum = sum(float(d_vals[f'AVPChi{i}']) for i in range(11) if d_vals[f'AVPChi{i}'] is not None and not np.isnan(d_vals[f'AVPChi{i}']))
                asp_chi_sum = sum(float(d_vals[f'ASPChi{i}']) for i in range(11) if d_vals[f'ASPChi{i}'] is not None and not np.isnan(d_vals[f'ASPChi{i}']))
                
                delta_chi = avp_chi_sum - asp_chi_sum
                product = i_max * delta_chi

                if product < min_product:
                    min_product = product
                    result_details = {
                        'name': name,
                        'i_max': i_max,
                        'delta_chi': delta_chi,
                        'product': product,
                        'max_gats': max_gats
                    }

        # 4. Print the final result
        if not result_details:
            print("No homologs found that match the Geary autocorrelation criteria (max value between 2 and 3).")
        else:
            print(f"The homolog with the minimum product is '{result_details['name']}'.")
            print("Its maximum Geary autocorrelation value is {:.4f}.".format(result_details['max_gats']))
            print("\nThe final equation is derived from the following numbers:")
            print(f"  Lag (i_max) for the maximum Geary autocorrelation: {result_details['i_max']}")
            print(f"  Difference between average path chi indices (delta_chi): {result_details['delta_chi']:.4f}")
            print(f"\nMinimum Product = i_max * delta_chi")
            print(f"Minimum Product = {result_details['i_max']} * {result_details['delta_chi']:.4f} = {result_details['product']:.4f}")
            print(f"<<<{result_details['product']:.4f}>>>")

    except ImportError:
        print("Error: RDKit or Mordred library not found.")
        print("Please install them using: pip install rdkit-pypi mordred")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    find_minimum_product_of_homologs()