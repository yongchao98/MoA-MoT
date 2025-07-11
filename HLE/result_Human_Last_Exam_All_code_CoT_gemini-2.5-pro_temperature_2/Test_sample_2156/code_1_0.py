import sys
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    from mordred import Calculator, descriptors
except ImportError:
    print("This script requires rdkit and mordred.", file=sys.stderr)
    print("Please install them using: pip install rdkit-pypi mordred", file=sys.stderr)
    sys.exit(1)

def solve_molecular_property_task():
    """
    Solves the user's request by finding the minimum product of i_max and
    a delta chi value for a specific set of formaldehyde homologs.
    """
    # 1. Define formaldehyde's homologs (alkanals) from C1 to C20
    alkanals = []
    for i in range(20):
        if i == 0:
            name = "Formaldehyde (C1)"
            smiles = "C=O"
        else:
            name = f"C{i + 1}-anal"
            smiles = 'C' * (i + 1) + '=O'
        alkanals.append({"name": name, "smiles": smiles})

    # 2. Setup Mordred calculator for Geary Autocorrelation (GATS)
    # We will use descriptors weighted by Sanderson electronegativity ('s') for lags 1-10.
    max_lag = 10
    gats_instances = [descriptors.GearyAuto(lag=i, prop='En_s') for i in range(1, max_lag + 1)]
    calc_gats = Calculator(gats_instances, ignore_3D=True)

    found_homologs = []

    # 3. Iterate through homologs, calculate GATS, and filter them
    for alkanal in alkanals:
        mol = Chem.MolFromSmiles(alkanal["smiles"])
        if mol is None:
            continue
            
        mol_with_hs = Chem.AddHs(mol)
        
        gats_values = {}
        try:
            mordred_results = calc_gats(mol_with_hs)
            for i in range(1, max_lag + 1):
                key = f'GATS{i}s'
                val = mordred_results[key]
                # Ensure the result is a valid number
                if val is not None and not np.isnan(val) and not isinstance(val, Exception):
                    gats_values[i] = val
        except Exception:
            continue

        if not gats_values:
            continue
            
        # Find the maximum GATS value and the corresponding lag (i_max)
        gats_max = max(gats_values.values())
        i_max = max(gats_values, key=gats_values.get)

        # Apply the filtering condition: 2 < max(GATS) < 3
        if 2 < gats_max < 3:
            found_homologs.append({
                'name': alkanal["name"],
                'mol_with_hs': mol_with_hs,
                'i_max': i_max,
                'gats_max': gats_max
            })
    
    if not found_homologs:
        print("No homologs found matching the criteria 2 < max(GATS) < 3.")
        return

    # 4. Calculate Chi indices and the final product for each found homolog
    results = []
    chi_n_funcs = [getattr(rdMolDescriptors, f'CalcChi{i}n') for i in range(11)]
    chi_v_funcs = [getattr(rdMolDescriptors, f'CalcChi{i}v') for i in range(11)]

    for homolog in found_homologs:
        mol = homolog['mol_with_hs']
        
        sum_n = sum(f(mol) for f in chi_n_funcs)
        sum_v = sum(f(mol) for f in chi_v_funcs)
        
        num_orders = float(len(chi_n_funcs))
        avg_simple_chi = sum_n / num_orders
        avg_valence_chi = sum_v / num_orders
        
        delta_chi = avg_valence_chi - avg_simple_chi
        product = homolog['i_max'] * delta_chi
        
        homolog['avg_simple_chi'] = avg_simple_chi
        homolog['avg_valence_chi'] = avg_valence_chi
        homolog['product'] = product
        results.append(homolog)

    if not results:
        print("Calculation failed for all found homologs.")
        return

    # 5. Find the minimum product among all qualifying homologs
    min_result = min(results, key=lambda x: x['product'])

    # 6. Output the final result as requested
    print("--- Calculation for the minimum product ---")
    print(f"Molecule with minimum product: {min_result['name']}")
    print(f"Maximum Geary autocorrelation lag (i_max): {min_result['i_max']}")
    print(f"Average valence path chi index: {min_result['avg_valence_chi']:.4f}")
    print(f"Average simple path chi index: {min_result['avg_simple_chi']:.4f}")
    
    print("\nFinal Equation:")
    print(f"{min_result['i_max']} * ({min_result['avg_valence_chi']:.4f} - {min_result['avg_simple_chi']:.4f}) = {min_result['product']:.4f}")
    
    print(f"\n<<<{min_result['product']:.4f}>>>")

if __name__ == "__main__":
    solve_molecular_property_task()