import sys
import io
from contextlib import redirect_stdout

# This suppresses noisy RDKit logs to keep the output clean
old_stdout = sys.stdout
sys.stdout = io.StringIO()
try:
    from rdkit import Chem
    from rdkit import RDLogger
    # Disable RDKit's console output
    RDLogger.DisableLog('rdApp.*')
finally:
    # Restore standard output
    sys.stdout = old_stdout

from mordred import Calculator, descriptors

def solve_chemistry_problem():
    """
    This function executes the plan to find the minimum product of i_max and the
    difference between average valence and simple path chi indices for specific
    formaldehyde homologs.
    """
    # 1. Define the homologous series of aldehydes
    homologs = {
        'Formaldehyde': 'C=O',
        'Ethanal': 'CC=O',
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
        'Tetradecanal': 'CCCCCCCCCCCCCC=O',
        'Pentadecanal': 'CCCCCCCCCCCCCCC=O'
    }

    # 2. Set up the mordred calculator with the required descriptors
    descriptor_list = []
    # Geary autocorrelations (lag 1-8) weighted by Sanderson EN
    for i in range(1, 9):
        descriptor_list.append(descriptors.Autocorrelation.Geary(lag=i, prop='electronegativity_sanderson'))
    
    # Simple and Valence Path chi indices (order 0-10)
    for i in range(11):
        descriptor_list.append(descriptors.PathChi.PathChi(order=i, valence=False))
        descriptor_list.append(descriptors.PathChi.PathChi(order=i, valence=True))

    calc = Calculator(descriptor_list, ignore_3D=True)

    # 3. & 4. Iterate through homologs, calculate descriptors, and filter
    found_homologs_results = []

    for name, smiles in homologs.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        # Descriptors are sensitive to hydrogens
        mol = Chem.AddHs(mol)
        
        # Calculate all descriptors for the molecule
        desc_values = calc(mol)
        
        # Extract Geary values and find the maximum and corresponding lag (i_max)
        geary_vals = {}
        for i in range(1, 9):
            val = desc_values[f'GATSe{i}']
            # Handle potential calculation errors from mordred
            if val and not isinstance(val, (float, int)) and val.is_missing:
                val = 0.0
            geary_vals[i] = float(val)

        if not geary_vals:
            continue
        
        i_max = max(geary_vals, key=geary_vals.get)
        max_geary_val = geary_vals[i_max]

        # Check if the homolog meets the Geary autocorrelation condition
        if 2.0 <= max_geary_val <= 3.0:
            # If it does, calculate the average chi indices
            chi_p_vals = []
            chi_pv_vals = []
            for i in range(11):
                p_val = desc_values[f'chi{i}p']
                pv_val = desc_values[f'chi{i}pv']
                
                # Sanitize values in case of calculation errors
                p_val = 0.0 if (p_val and not isinstance(p_val, (float, int)) and p_val.is_missing) else float(p_val)
                pv_val = 0.0 if (pv_val and not isinstance(pv_val, (float, int)) and pv_val.is_missing) else float(pv_val)

                chi_p_vals.append(p_val)
                chi_pv_vals.append(pv_val)

            avg_chi_p = sum(chi_p_vals) / len(chi_p_vals)
            avg_chi_pv = sum(chi_pv_vals) / len(chi_pv_vals)
            
            chi_diff = avg_chi_pv - avg_chi_p
            product = i_max * chi_diff

            # Store all relevant data for the qualifying homolog
            found_homologs_results.append({
                "name": name,
                "i_max": i_max,
                "max_geary": max_geary_val,
                "avg_chi_p": avg_chi_p,
                "avg_chi_pv": avg_chi_pv,
                "chi_diff": chi_diff,
                "product": product,
                "chi_p_vals_raw": chi_p_vals,
                "chi_pv_vals_raw": chi_pv_vals
            })

    # 5. Find the homolog with the minimum product among those found
    if not found_homologs_results:
        print("No formaldehyde homologs were found with a maximum Geary autocorrelation (weighted by Sanderson EN) between 2.0 and 3.0.")
        return

    min_result = min(found_homologs_results, key=lambda x: x['product'])

    # 6. Print the detailed, step-by-step calculation for the minimum case
    print(f"The homolog with the minimum product is {min_result['name']}.")
    print(f"This was determined from the set of aldehydes from C1 to C15.")
    print("\n--- Calculation Breakdown ---")
    
    print(f"\n1. Maximum Geary Autocorrelation:")
    print(f"   - The maximum Geary value for {min_result['name']} is {min_result['max_geary']:.4f}.")
    print(f"   - This maximum occurs at lag i_max = {min_result['i_max']}.")
    
    print("\n2. Average Path Chi Indices (calculated over orders 0-10):")
    
    chi_p_str = " + ".join([f"{v:.4f}" for v in min_result['chi_p_vals_raw']])
    print(f"   - Average simple path chi index = ({chi_p_str}) / 11 = {min_result['avg_chi_p']:.4f}")
    
    chi_pv_str = " + ".join([f"{v:.4f}" for v in min_result['chi_pv_vals_raw']])
    print(f"   - Average valence path chi index = ({chi_pv_str}) / 11 = {min_result['avg_chi_pv']:.4f}")
    
    print(f"\n3. Difference between Average Chi Indices:")
    print(f"   - Difference = {min_result['avg_chi_pv']:.4f} - {min_result['avg_chi_p']:.4f} = {min_result['chi_diff']:.4f}")
    
    print("\n4. Final Product Calculation:")
    print(f"   - Final Product = i_max * (Difference)")
    print(f"   - Final Product = {min_result['i_max']} * {min_result['chi_diff']:.4f} = {min_result['product']:.4f}")

    # 7. Provide the final answer in the required format
    print(f"\n<<<{min_result['product']:.4f}>>>")

# Execute the main function to solve the problem
solve_chemistry_problem()