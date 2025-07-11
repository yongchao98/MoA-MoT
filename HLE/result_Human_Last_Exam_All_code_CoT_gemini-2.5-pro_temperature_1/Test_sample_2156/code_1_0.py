import sys
import warnings
from rdkit import Chem
from mordred import Calculator, descriptors
import numpy as np

# Suppress warnings from RDKit and Mordred for cleaner output
warnings.filterwarnings("ignore")

def solve_chemistry_problem():
    """
    Solves the described chemical informatics problem.
    """
    try:
        # 1. Define homologs of formaldehyde (aldehydes C1-C12)
        homologs = {
            'Formaldehyde': 'C=O',
            'Acetaldehyde': 'CC=O',
            'Propionaldehyde': 'CCC=O',
            'Butyraldehyde': 'CCCC=O',
            'Pentanal': 'CCCCC=O',
            'Hexanal': 'CCCCCC=O',
            'Heptanal': 'CCCCCCC=O',
            'Octanal': 'CCCCCCCC=O',
            'Nonanal': 'CCCCCCCCC=O',
            'Decanal': 'CCCCCCCCCC=O',
            'Undecanal': 'CCCCCCCCCCC=O',
            'Dodecanal': 'CCCCCCCCCCCC=O',
        }

        # 2. Initialize Mordred calculator
        calc = Calculator(descriptors, ignore_3D=True)
        
        found_homologs_results = []

        # 3. Iterate through each homolog
        for name, smiles in homologs.items():
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            mol = Chem.AddHs(mol)

            # Calculate all descriptors
            all_descs = calc(mol)
            # Replace errors/missing values with NaN for safe processing
            all_descs.fill_missing(value=np.nan)

            # Extract Geary autocorrelation descriptors (GATS)
            geary_descs = {k: v for k, v in all_descs.items() if str(k).startswith('GATS') and str(k).endswith('se')}
            
            if not geary_descs or all(np.isnan(v) for v in geary_descs.values()):
                continue

            # Find max Geary value and its corresponding lag (i_max)
            max_geary_val = -np.inf
            i_max = 0
            for desc_name, value in geary_descs.items():
                if not np.isnan(value) and value > max_geary_val:
                    max_geary_val = value
                    i_max = int(str(desc_name).replace('GATS', '').replace('se', ''))

            # 4. Filter homologs where 2 <= max_geary_val <= 3
            if 2 <= max_geary_val <= 3:
                
                # 5. Calculate average path chi indices
                # Simple path chi indices (X0, X1, ..., X9)
                spc_keys = {f'X{i}' for i in range(10)}
                spc_vals = [v for k, v in all_descs.items() if str(k) in spc_keys and not np.isnan(v)]
                avg_spc = np.mean(spc_vals) if spc_vals else 0

                # Valence path chi indices (X0v, X1v, ..., X9v)
                vpc_keys = {f'X{i}v' for i in range(10)}
                vpc_vals = [v for k, v in all_descs.items() if str(k) in vpc_keys and not np.isnan(v)]
                avg_vpc = np.mean(vpc_vals) if vpc_vals else 0
                
                # Calculate the difference and the final product
                diff = avg_vpc - avg_spc
                product = i_max * diff
                
                found_homologs_results.append({
                    'name': name,
                    'i_max': i_max,
                    'avg_vpc': avg_vpc,
                    'avg_spc': avg_spc,
                    'product': product
                })

        # 6. Find the homolog with the minimum product
        if not found_homologs_results:
            print("No homologs were found that satisfy the Geary autocorrelation condition.")
            return

        min_result = min(found_homologs_results, key=lambda x: x['product'])

        # 7. Print the final result in the required format
        i = min_result['i_max']
        vpc = min_result['avg_vpc']
        spc = min_result['avg_spc']
        p = min_result['product']
        
        print(f"{i} * ({vpc:.4f} - {spc:.4f}) = {p:.4f}")
        print(f"<<<{p}>>>")

    except ImportError:
        print("Required libraries (rdkit, mordred, numpy) are not installed.", file=sys.stderr)
        print("Please install them using: pip install rdkit mordred numpy", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

if __name__ == "__main__":
    solve_chemistry_problem()