import numpy as np
from rdkit import Chem
from rdkit.Chem import rdmolops, rdMolDescriptors

def solve():
    """
    Finds the minimum product of i_max and a chi index difference for formaldehyde homologs.
    """
    
    # Plan Step 1 & 2: Define molecules and atomic properties
    aldehyde_names = [
        "Formaldehyde", "Acetaldehyde", "Propanal", "Butanal", 
        "Pentanal", "Hexanal", "Heptanal", "Octanal", "Nonanal", "Decanal"
    ]
    aldehyde_smiles = [
        'C=O', 'CC=O', 'CCC=O', 'CCCC=O', 'CCCCC=O', 'CCCCCC=O',
        'CCCCCCC=O', 'CCCCCCCC=O', 'CCCCCCCCC=O', 'CCCCCCCCCC=O'
    ]
    sanderson_en = {'H': 2.592, 'C': 2.746, 'O': 3.654}

    # Helper function to calculate Geary Autocorrelation
    def calculate_geary_autocorrelation(mol, lag, atomic_props):
        """Calculates the Geary autocorrelation for a given molecule, lag, and atomic property list."""
        num_atoms = mol.GetNumAtoms()
        if num_atoms < 2 or lag == 0:
            return 0.0

        dist_matrix = rdmolops.GetDistanceMatrix(mol)
        prop_avg = np.mean(atomic_props)
        
        denominator_sum_sq_dev = np.sum([(prop - prop_avg)**2 for prop in atomic_props])
        
        if denominator_sum_sq_dev == 0:
            return 0.0

        numerator_sum_sq_diff = 0.0
        num_pairs_at_lag = 0
        
        for i in range(num_atoms):
            for j in range(i + 1, num_atoms):
                if dist_matrix[i, j] == lag:
                    numerator_sum_sq_diff += (atomic_props[i] - atomic_props[j])**2
                    num_pairs_at_lag += 1

        if num_pairs_at_lag == 0:
            return 0.0
            
        numerator = (num_atoms - 1) * numerator_sum_sq_diff
        denominator = (2 * num_pairs_at_lag) * denominator_sum_sq_dev
        
        return numerator / denominator if denominator != 0 else 0.0

    found_homologs_data = []

    print("Searching for homologs with max Geary Autocorrelation between 2 and 3...")
    
    # Plan Step 3: Iterate and Calculate Geary Autocorrelation
    for name, smiles in zip(aldehyde_names, aldehyde_smiles):
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        
        atomic_props = [sanderson_en[atom.GetSymbol()] for atom in mol.GetAtoms()]
        
        dist_matrix = rdmolops.GetDistanceMatrix(mol)
        max_lag = int(np.max(dist_matrix))
        
        geary_values = [calculate_geary_autocorrelation(mol, lag, atomic_props) for lag in range(1, max_lag + 1)]

        if not geary_values:
            continue

        max_geary = max(geary_values)
        i_max = geary_values.index(max_geary) + 1
        
        # Plan Step 4: Filter Homologs
        if 2 <= max_geary <= 3:
            
            # Plan Step 5: Calculate Chi Indices and Product
            print(f"\nFound homolog: {name}")
            print(f"  - Max Geary Autocorrelation = {max_geary:.4f} (at lag i_max = {i_max})")

            chi_funcs_n = [
                rdMolDescriptors.CalcChi0n, rdMolDescriptors.CalcChi1n,
                rdMolDescriptors.CalcChi2n, rdMolDescriptors.CalcChi3n,
                rdMolDescriptors.CalcChi4n, rdMolDescriptors.CalcChi5n,
                rdMolDescriptors.CalcChi6n, rdMolDescriptors.CalcChi7n,
                rdMolDescriptors.CalcChi8n, rdMolDescriptors.CalcChi9n,
                rdMolDescriptors.CalcChi10n,
            ]
            chi_funcs_v = [
                rdMolDescriptors.CalcChi0v, rdMolDescriptors.CalcChi1v,
                rdMolDescriptors.CalcChi2v, rdMolDescriptors.CalcChi3v,
                rdMolDescriptors.CalcChi4v, rdMolDescriptors.CalcChi5v,
                rdMolDescriptors.CalcChi6v, rdMolDescriptors.CalcChi7v,
                rdMolDescriptors.CalcChi8v, rdMolDescriptors.CalcChi9v,
                rdMolDescriptors.CalcChi10v,
            ]
        
            simple_chi_vals = [f(mol) for f in chi_funcs_n]
            valence_chi_vals = [f(mol) for f in chi_funcs_v]

            avg_chi_simple = np.mean(simple_chi_vals)
            avg_chi_valence = np.mean(valence_chi_vals)
            
            delta_chi = avg_chi_valence - avg_chi_simple
            product = i_max * delta_chi
            
            found_homologs_data.append({
                'name': name,
                'i_max': i_max,
                'avg_chi_simple': avg_chi_simple,
                'avg_chi_valence': avg_chi_valence,
                'product': product
            })

    print("\n--------------------------------------------------")

    # Plan Step 6: Determine Minimum Product
    if not found_homologs_data:
        print("No homologs were found that satisfy the Geary autocorrelation condition.")
    else:
        min_result = min(found_homologs_data, key=lambda x: x['product'])
        
        # Plan Step 7: Output the Result
        print("Calculation for the homolog with the minimum product:")
        print(f"Homolog: {min_result['name']}")
        
        print("\nFinal equation: i_max * (average_valence_chi - average_simple_chi) = minimum_product")
        print("Values:")
        print(f"i_max = {min_result['i_max']}")
        print(f"average_valence_chi = {min_result['avg_chi_valence']:.4f}")
        print(f"average_simple_chi = {min_result['avg_chi_simple']:.4f}")
        
        # Explicitly show the equation with numbers as requested
        print(f"\n{min_result['i_max']} * ({min_result['avg_chi_valence']:.4f} - {min_result['avg_chi_simple']:.4f}) = {min_result['product']:.4f}")

        print(f"\n<<<__{min_result['product']:.4f}__>>>")

solve()