import rdkit
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem import Descriptors
import numpy as np
import sys

def solve_chemoinformatics_problem():
    """
    This script finds formaldehyde's homologs that satisfy a specific Geary autocorrelation
    condition and then determines the minimum product of i_max and a chi index difference
    among them.
    """

    # --- Helper function to calculate Geary Autocorrelation (Full Formula) ---
    def calculate_geary_autocorrelation(mol):
        # Sanderson electronegativities for relevant atoms
        sander_en = {'H': 2.592, 'C': 2.746, 'O': 3.654}
        num_atoms = mol.GetNumAtoms()
        if num_atoms <= 1:
            return 0, 0
        
        try:
            props = np.array([sander_en[atom.GetSymbol()] for atom in mol.GetAtoms()])
        except KeyError as e:
            print(f"Error: Atom type {e} not found in Sanderson EN dictionary.", file=sys.stderr)
            return 0, 0
        
        # Calculate denominator term: sum of squared deviations from the mean
        prop_mean = np.mean(props)
        prop_dev_sq_sum = np.sum((props - prop_mean)**2)
        if prop_dev_sq_sum == 0:
            return 0, 0

        dist_matrix = Chem.GetDistanceMatrix(mol)
        max_lag = int(np.max(dist_matrix))
        
        gats_values = {}
        
        # Constant scaling factor part from the Geary formula
        scaling_factor1 = (num_atoms - 1) / prop_dev_sq_sum

        for lag in range(1, max_lag + 1):
            # Find pairs of atoms at the specified lag
            rows, cols = np.where(dist_matrix == lag)
            valid_pairs_mask = rows < cols
            rows, cols = rows[valid_pairs_mask], cols[valid_pairs_mask]

            pair_count = len(rows)
            if pair_count > 0:
                # Numerator term: sum of squared differences for pairs at distance `lag`
                prop_diffs_sq = (props[rows] - props[cols])**2
                numerator_term = np.sum(prop_diffs_sq)

                # Final scaling factor for this lag
                scaling_factor2 = 1 / (2 * pair_count)
                
                gats_values[lag] = scaling_factor1 * scaling_factor2 * numerator_term

        if not gats_values:
            return 0, 0

        # Find the maximum GATS value and the corresponding lag
        i_max = max(gats_values, key=gats_values.get)
        max_gats = gats_values[i_max]
        
        return max_gats, i_max

    # --- Helper function to calculate Chi Index information ---
    def calculate_chi_index_info(mol_with_h):
        mol_no_h = Chem.RemoveHs(mol_with_h)
        num_heavy_atoms = mol_no_h.GetNumAtoms()
        
        if num_heavy_atoms == 0:
            return 0, 0, 0

        # Simple vertex degrees (v) on the H-suppressed graph
        v = np.array([atom.GetDegree() for atom in mol_no_h.GetAtoms()], dtype=float)

        # Valence vertex degrees (v^v) on the H-suppressed graph
        Z_v = {'C': 4, 'O': 6}
        v_v = np.array([Z_v[atom.GetSymbol()] - atom.GetTotalNumHs() for atom in mol_no_h.GetAtoms()], dtype=float)

        # Determine the maximum path length (k)
        max_k = int(np.max(Chem.GetDistanceMatrix(mol_no_h))) if num_heavy_atoms > 1 else 0

        chi_p_values = []
        chi_pv_values = []
        
        # Order k=0
        with np.errstate(divide='ignore', invalid='ignore'):
            inv_sqrt_v = np.power(v, -0.5)
            inv_sqrt_v_v = np.power(v_v, -0.5)
        inv_sqrt_v[np.isinf(inv_sqrt_v)] = 0
        inv_sqrt_v_v[np.isinf(inv_sqrt_v_v)] = 0
        chi_p_values.append(np.sum(inv_sqrt_v))
        chi_pv_values.append(np.sum(inv_sqrt_v_v))

        # Orders k > 0
        for k in range(1, max_k + 1):
            paths = rdmolops.FindAllPathsOfLengthN(mol_no_h, k + 1, useBonds=False, useHs=False)
            chi_p_k = 0.0
            chi_pv_k = 0.0
            for path in paths:
                path_indices = list(path)
                
                prod_p = np.prod(v[path_indices])
                if prod_p > 0:
                    chi_p_k += prod_p**-0.5

                prod_pv = np.prod(v_v[path_indices])
                if prod_pv > 0:
                    chi_pv_k += prod_pv**-0.5
            
            chi_p_values.append(chi_p_k)
            chi_pv_values.append(chi_pv_k)

        num_orders = len(chi_p_values)
        if num_orders == 0:
            return 0, 0, 0
        
        avg_chi_p = sum(chi_p_values) / num_orders
        avg_chi_pv = sum(chi_pv_values) / num_orders

        return avg_chi_pv, avg_chi_p, (avg_chi_pv - avg_chi_p)

    # --- Main Loop ---
    min_product = float('inf')
    best_homolog_info = {}

    # Check n-aldehydes from propanal (n=3) up to n=30
    for n_carbons in range(3, 31): 
        smiles = 'C' * n_carbons + '=O'
        mol = Chem.MolFromSmiles(smiles)
        mol_h = Chem.AddHs(mol)

        max_gats, i_max = calculate_geary_autocorrelation(mol_h)

        if 2.0 <= max_gats <= 3.0:
            avg_chi_pv, avg_chi_p, diff = calculate_chi_index_info(mol_h)
            product = i_max * diff
            
            if product < min_product:
                min_product = product
                best_homolog_info = {
                    'formula': Descriptors.CalcMolFormula(mol_h),
                    'n_carbons': n_carbons,
                    'i_max': i_max,
                    'avg_chi_pv': avg_chi_pv,
                    'avg_chi_p': avg_chi_p,
                    'diff': diff,
                    'product': product
                }

    if not best_homolog_info:
        print("No homologs found in the specified Geary autocorrelation range [2.0, 3.0].")
        return

    # Print the detailed breakdown for the homolog with the minimum product
    info = best_homolog_info
    name = f"n-C{info['n_carbons']} aldehyde (heptacosanal)" if info['n_carbons'] == 27 else f"n-C{info['n_carbons']} aldehyde"
    
    print(f"Homolog with minimum product: {name} ({info['formula']})")
    print(f"\nThe final equation is: i_max * (average_valence_chi - average_simple_chi)")
    print(f"The numbers for the final equation are:")
    print(f"i_max = {info['i_max']}")
    print(f"average_valence_chi = {info['avg_chi_pv']:.5f}")
    print(f"average_simple_chi = {info['avg_chi_p']:.5f}")
    print(f"\nResulting in:")
    print(f"Difference = {info['diff']:.5f}")
    print(f"Minimum Product = {info['product']:.5f}")

    # Final answer in specified format
    print(f"\n<<<{info['product']:.5f}>>>")

solve_chemoinformatics_problem()