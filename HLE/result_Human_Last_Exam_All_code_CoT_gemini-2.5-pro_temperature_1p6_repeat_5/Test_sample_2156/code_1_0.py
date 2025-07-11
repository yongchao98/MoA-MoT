import numpy as np
import math
from rdkit import Chem
from rdkit.Chem import Descriptors

def solve():
    """
    Solves the problem by finding formaldehyde's homologs that satisfy certain
    molecular descriptor conditions and calculating a target value from them.
    """

    # Sanderson Electronegativities for relevant atoms
    SANDERSON_EN = {'H': 2.59, 'C': 2.75, 'O': 3.65}

    # List of alkanal homologs (n=1 to 14) by SMILES and common name
    HOMOLOGS_DATA = {
        'Formaldehyde': 'C=O',
        'Acetaldehyde': 'CC=O',
        'Propanal': 'CCC=O',
        'Butanal': 'CCCC=O',
        'Pentanal': 'CCCCC=O',
        'Hexanal': 'CCCCCC=O',
        'Heptanal': 'CCCCCCC=O',
        'Octanal': 'CCCCCCCC=O',
        'Nonanal': 'CCCCCCCCC=O',
        'Decanal': 'CCCCCCCCCC=O',
        'Undecanal': 'CCCCCCCCCCC=O',
        'Dodecanal': 'CCCCCCCCCCCC=O',
        'Tridecanal': 'CCCCCCCCCCCCC=O',
        'Tetradecanal':'CCCCCCCCCCCCCC=O',
    }

    def calculate_geary_autocorrelation(mol):
        """Calculates Geary autocorrelation for all lags for a molecule."""
        num_atoms = mol.GetNumAtoms()
        if num_atoms < 2:
            return {}, -1, -1

        weights = [SANDERSON_EN[atom.GetSymbol()] for atom in mol.GetAtoms()]
        dist_matrix = Chem.GetDistanceMatrix(mol)
        avg_w = np.mean(weights)
        variance_sum = sum([(w - avg_w)**2 for w in weights])
        
        if variance_sum == 0:
            return {}, -1, -1

        # Denominator of the Geary formula
        denom = variance_sum / (num_atoms - 1)
        if denom == 0:
             return {}, -1, -1

        max_lag = int(np.max(dist_matrix))
        geary_by_lag = {}

        for lag in range(1, max_lag + 1):
            sum_sq_diff = 0
            num_pairs = 0
            for i in range(num_atoms):
                for j in range(i + 1, num_atoms):
                    if dist_matrix[i, j] == lag:
                        num_pairs += 1
                        sum_sq_diff += (weights[i] - weights[j])**2
            
            if num_pairs > 0:
                # Numerator is the average squared difference at the given lag
                numerator = sum_sq_diff / num_pairs
                geary_by_lag[lag] = numerator / denom
            else:
                geary_by_lag[lag] = 0

        if not geary_by_lag:
            return {}, -1, -1

        max_geary_val = max(geary_by_lag.values())
        # In case of ties for max value, this gets the first max lag
        max_lag_val = max(geary_by_lag, key=geary_by_lag.get)
        
        return geary_by_lag, max_geary_val, max_lag_val

    def get_simple_avg_ipc(mol):
        """
        Calculates the Average Simple Path Chi Index.
        This is analogous to RDKit's AvgIpc but uses simple atom connectivities.
        The calculation is based on the information content of path length partitions.
        """
        num_atoms = mol.GetNumAtoms()
        if num_atoms < 2:
            return 0.0

        path_counts = []
        # Max path length is num_atoms - 1
        for i in range(1, num_atoms):
            # FindAllPathsOfLengthN uses atom indices
            paths = Chem.FindAllPathsOfLengthN(mol, i, useBonds=False, useHs=True)
            if paths:
                path_counts.append(len(paths))
            else:
                break
        
        if not path_counts:
            return 0.0
            
        total_paths = sum(path_counts)
        if total_paths == 0:
            return 0.0

        # Information content calculation
        entropy_sum = 0.0
        for count in path_counts:
            if count > 0:
                entropy_sum += count * math.log(count, 2)
        
        if total_paths == 1 and entropy_sum == 0:
            ipc = 0.0 # handle log2(1) = 0 case
        else:
            ipc = total_paths * math.log(total_paths, 2) - entropy_sum
        
        # Return average IPC
        return ipc / num_atoms

    def calculate_ipc_difference(mol):
        """
        Calculates the difference between average valence and simple path chi indices.
        """
        # RDKit's built-in average valence path chi index
        avg_valence_ipc = Descriptors.AvgIpc(mol)
        
        # Custom-calculated average simple path chi index
        avg_simple_ipc = get_simple_avg_ipc(mol)
        
        return avg_valence_ipc - avg_simple_ipc

    found_homologs = []

    for name, smiles in HOMOLOGS_DATA.items():
        mol = Chem.MolFromSmiles(smiles)
        mol_h = Chem.AddHs(mol)

        _, max_geary, i_max = calculate_geary_autocorrelation(mol_h)
        
        if 2 <= max_geary <= 3:
            diff_ipc = calculate_ipc_difference(mol_h)
            product = i_max * diff_ipc
            found_homologs.append({
                'name': name,
                'i_max': i_max,
                'diff_ipc': diff_ipc,
                'product': product,
            })
    
    if not found_homologs:
        print("No homologs were found that satisfy the condition 2 <= max(Geary) <= 3.")
        return

    # Find the homolog with the minimum product
    min_entry = min(found_homologs, key=lambda x: x['product'])

    # Output the result
    print("The minimum product is found among the homologs whose maximum Geary autocorrelation (weighted by Sanderson electronegativity) is between 2 and 3.")
    print("The found homologs are:")
    for entry in found_homologs:
        print(f"- {entry['name']}")

    print("\nThe minimum product belongs to:")
    print(f"{min_entry['name']}")
    
    print("\nThe final equation is based on the following values:")
    print(f"Lag at maximum Geary autocorrelation (i_max): {min_entry['i_max']}")
    print(f"Difference between average chi indices (AvgValenceIpc - AvgSimpleIpc): {min_entry['diff_ipc']:.8f}")
    
    print(f"\nResulting equation:")
    print(f"{min_entry['i_max']} * {min_entry['diff_ipc']:.8f} = {min_entry['product']:.8f}")

    print(f"<<<{min_entry['product']}>>>")

solve()