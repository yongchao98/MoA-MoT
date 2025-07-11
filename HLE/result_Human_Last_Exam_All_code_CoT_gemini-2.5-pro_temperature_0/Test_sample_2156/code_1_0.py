import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def solve_cheminformatics_problem():
    """
    This script finds formaldehyde's homologs with a specific Geary autocorrelation range,
    and then determines the minimum product of the max-value lag (i_max) and the
    difference between average valence and simple path chi indices.
    """

    # --- Helper Functions ---

    # Sanderson electronegativities for relevant atoms
    s_en = {'H': 2.592, 'C': 2.746, 'O': 3.654}

    def get_s_en(atom):
        """Returns the Sanderson electronegativity for a given RDKit atom."""
        return s_en.get(atom.GetSymbol(), 0)

    def calculate_geary_autocorrelation(mol):
        """
        Calculates Geary autocorrelation weighted by Sanderson electronegativities.
        Returns the max value and the lag at which it occurs.
        """
        num_atoms = mol.GetNumAtoms()
        if num_atoms < 2:
            return 0, 0

        # Get atomic properties (electronegativity) and calculate the average
        props = np.array([get_s_en(atom) for atom in mol.GetAtoms()])
        avg_prop = np.mean(props)
        
        # Denominator of the Geary formula
        denominator = np.sum((props - avg_prop)**2)
        if denominator == 0:
            return 0, 0

        dist_matrix = Chem.GetDistanceMatrix(mol)
        max_lag = int(np.max(dist_matrix))
        geary_values = {}

        for lag in range(1, max_lag + 1):
            sum_sq_diff = 0
            count_pairs = 0
            # Sum squared differences for all atom pairs at the current lag
            for i in range(num_atoms):
                for j in range(i + 1, num_atoms):
                    if dist_matrix[i, j] == lag:
                        sum_sq_diff += (props[i] - props[j])**2
                        count_pairs += 1
            
            if count_pairs > 0:
                # Standard formula for Geary Autocorrelation of Topological Structure (GATS)
                geary_d = ((num_atoms - 1) / (2.0 * count_pairs)) * (sum_sq_diff / denominator)
                geary_values[lag] = geary_d

        if not geary_values:
            return 0, 0

        # Find the maximum Geary value and the corresponding lag (i_max)
        max_geary = max(geary_values.values())
        i_max = max(geary_values, key=geary_values.get)

        return max_geary, i_max

    def calculate_avg_chi_indices(mol):
        """
        Calculates the average simple and valence path chi indices.
        """
        num_atoms = mol.GetNumAtoms()
        if num_atoms == 0:
            return 0, 0

        # The maximum order for chi indices is the maximum path length (graph diameter)
        max_order = int(np.max(Chem.GetDistanceMatrix(mol)))

        # Sum chi indices from order 0 to max_order
        sum_path = sum(rdMolDescriptors.CalcChiNpath(mol, i) for i in range(max_order + 1))
        sum_wpath = sum(rdMolDescriptors.CalcChiNpathv(mol, i) for i in range(max_order + 1))

        # Calculate the average by dividing by the number of atoms
        avg_path = sum_path / num_atoms
        avg_wpath = sum_wpath / num_atoms

        return avg_path, avg_wpath

    # --- Main Execution ---

    # Generate SMILES strings for aldehyde homologs (C1 to C20)
    homologs_smiles = ['C=O'] + [('C' * n) + 'C=O' for n in range(1, 20)]
    
    found_homologs_data = []

    for smiles in homologs_smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        mol = Chem.AddHs(mol)

        # Step 1 & 2: Calculate Geary autocorrelation and find max value/lag
        max_geary, i_max = calculate_geary_autocorrelation(mol)

        # Step 3: Filter homologs based on the max Geary value
        if 2.0 <= max_geary <= 3.0:
            # Step 4: Calculate average chi indices
            avg_path, avg_wpath = calculate_avg_chi_indices(mol)
            
            # Step 5: Calculate the final product
            diff = avg_wpath - avg_path
            product = i_max * diff
            
            found_homologs_data.append({
                'i_max': i_max,
                'avg_wpath': avg_wpath,
                'avg_path': avg_path,
                'product': product
            })

    if not found_homologs_data:
        print("Error: No homologs were found that match the specified criteria.")
        return

    # Step 6: Find the homolog with the minimum product
    min_result = min(found_homologs_data, key=lambda x: x['product'])

    # Output the final equation with the numbers that produce the minimum product
    print(f"{min_result['i_max']} * ({min_result['avg_wpath']} - {min_result['avg_path']}) = {min_result['product']}")

solve_cheminformatics_problem()