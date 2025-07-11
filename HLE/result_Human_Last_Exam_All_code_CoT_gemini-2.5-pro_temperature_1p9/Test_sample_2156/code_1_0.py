import numpy as np
from rdkit import Chem
from mordred import Calculator, descriptors

def solve_homologs_problem():
    """
    Solves the described cheminformatics problem by:
    1. Generating aldehydes as homologs of formaldehyde.
    2. Filtering them based on the maximum Geary autocorrelation value.
    3. Calculating the minimum product of i_max and a chi index difference.
    """
    # 1. Generate a series of aldehydes (homologs of formaldehyde)
    # SMILES strings for C1 to C20 aldehydes.
    # Formaldehyde is 'C=O'. Homologs add CH2 groups, e.g., Acetaldehyde is 'CC=O'.
    smiles_list = ['C=O'] + ['C' * n + 'C=O' for n in range(1, 20)]

    # 2. Set up mordred calculators for the required descriptors
    max_lag = 20
    
    # Geary autocorrelations with Sanderson electronegativity for lags 1 to 20
    geary_descs = [descriptors.GearyAuto(prop='se', lag=i) for i in range(1, max_lag + 1)]
    # Average simple and valence path chi indices
    chi_descs = [descriptors.Chi.Apc(), descriptors.Chi.Avpc()]

    calc_geary = Calculator(geary_descs, ignore_3D=True)
    calc_chi = Calculator(chi_descs, ignore_3D=True)
    
    filtered_homologs = []

    # Iterate through each aldehyde to find candidates
    for smi in smiles_list:
        # Create an RDKit molecule object with explicit hydrogens
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        if mol is None:
            continue
        
        # Calculate Geary autocorrelation values for all lags
        # mordred may return error objects for invalid calculations (e.g., lag > diameter)
        geary_results = calc_geary(mol)
        geary_values = [
            float(v) if not isinstance(v, descriptors.error.Error) else np.nan
            for v in geary_results
        ]
        
        # Find the maximum valid Geary value and its corresponding lag
        valid_geary_values = [v for v in geary_values if not np.isnan(v)]
        if not valid_geary_values:
            continue
            
        max_geary = max(valid_geary_values)
        # The lag is 1-based, corresponding to the index in the original list + 1
        i_max = geary_values.index(max_geary) + 1

        # Filter homologs where 2.0 <= max_geary <= 3.0
        if 2.0 <= max_geary <= 3.0:
            # For qualifying homologs, calculate chi indices
            # Result order matches descriptor list: [Apc, Avpc]
            chi_results = calc_chi(mol)
            apc = float(chi_results[0])
            avpc = float(chi_results[1])
            
            # Calculate the difference and the final product
            diff = avpc - apc
            product = i_max * diff
            
            # Store the results for this homolog
            filtered_homologs.append({
                'i_max': i_max,
                'apc': apc,
                'avpc': avpc,
                'product': product,
            })

    # Find the homolog with the minimum product
    if not filtered_homologs:
        print("No homologs were found that satisfy the Geary autocorrelation condition.")
        return

    min_homolog = min(filtered_homologs, key=lambda x: x['product'])

    # Output the final result in the requested format
    final_i_max = min_homolog['i_max']
    final_avpc = min_homolog['avpc']
    final_apc = min_homolog['apc']
    final_product = min_homolog['product']
    
    # Print the equation with all its numeric components
    print(f"{final_i_max} * ({final_avpc} - {final_apc}) = {final_product}")
    
    # Print the final answer in the specified format
    print(f"<<<{final_product}>>>")

# Execute the function
solve_homologs_problem()
