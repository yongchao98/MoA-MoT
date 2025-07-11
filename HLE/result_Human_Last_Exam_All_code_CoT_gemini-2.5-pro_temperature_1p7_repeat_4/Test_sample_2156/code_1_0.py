import sys
import numpy as np
from rdkit import Chem
from mordred import Calculator, descriptors, error
from mordred.geary import Geary
from mordred.path_chi import AverageValencePathChi, AverageSimplePathChi

def find_minimum_product():
    """
    Finds formaldehyde homologs satisfying a Geary autocorrelation condition,
    then determines the minimum product of i_max and a path chi index difference.
    """
    # Define the series of aldehyde homologs (from formaldehyde C1 to decanal C10)
    homolog_names = [
        "Formaldehyde", "Acetaldehyde", "Propanal", "Butanal",
        "Pentanal", "Hexanal", "Heptanal", "Octanal", "Nonanal", "Decanal"
    ]
    homolog_smiles = [
        "C=O", "CC=O", "CCC=O", "CCCC=O", "CCCCC=O",
        "CCCCCC=O", "CCCCCCC=O", "CCCCCCCC=O", "CCCCCCCCC=O", "CCCCCCCCCC=O"
    ]

    # Convert SMILES to RDKit molecule objects, adding hydrogens is good practice
    mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in homolog_smiles]

    # Create a calculator for the path chi indices we will need for the filtered molecules
    calc_chi = Calculator([AverageValencePathChi(), AverageSimplePathChi()], ignore_3D=True)

    found_homologs_data = []

    # Maximum lag to check for Geary autocorrelation
    max_lag = 10

    print("Analyzing formaldehyde homologs to find those meeting the criteria...")
    print("-" * 60)

    # 1. Iterate through molecules to find homologs that meet the criteria
    for name, mol in zip(homolog_names, mols):
        if mol is None:
            continue

        geary_values = []
        
        # Calculate Geary autocorrelation for all relevant lags
        for i in range(1, max_lag + 1):
            geary_desc = Geary(lag=i, prop='en_sanderson')
            calc_geary = Calculator([geary_desc], ignore_3D=True)
            
            try:
                result = calc_geary(mol)
                geary_val = result.values()[0]
                # Mordred returns a Missing object on failure for a given lag
                if isinstance(geary_val, error.Missing):
                    geary_val = np.nan
            except Exception:
                geary_val = np.nan
            
            geary_values.append(geary_val)
        
        # Ensure we have valid calculations before finding the max
        if np.all(np.isnan(geary_values)):
            continue

        # Find the maximum Geary value and its corresponding lag (ignoring NaNs)
        g_max = np.nanmax(geary_values)
        i_max = np.nanargmax(geary_values) + 1  # Add 1 as lags are 1-indexed

        # 2. Check if the homolog meets the condition 2.0 <= G_max <= 3.0
        print(f"Checking {name}: G_max = {g_max:.4f} at lag i_max = {i_max}")
        if 2.0 <= g_max <= 3.0:
            print(f"-> Found homolog: {name}")
            # 3. If it meets the condition, calculate path chi indices
            try:
                chi_indices = calc_chi(mol).asdict()
                avg_val_chi = chi_indices['Xp-Av']
                avg_sim_chi = chi_indices['Xp-As']
                
                if isinstance(avg_val_chi, error.Missing) or isinstance(avg_sim_chi, error.Missing):
                    continue
                    
                difference = avg_val_chi - avg_sim_chi
                product = i_max * difference

                # Store the results for this valid homolog
                found_homologs_data.append({
                    "name": name,
                    "i_max": i_max,
                    "difference": difference,
                    "product": product,
                })
            except Exception as e:
                print(f"-> Error calculating Path Chi indices for {name}: {e}")

    print("-" * 60)

    # 4. Find the minimum product among the found homologs
    if not found_homologs_data:
        print("No homologs were found that satisfy the Geary autocorrelation condition.")
        return

    print("Calculating final product for all found homologs:")
    for data in found_homologs_data:
        print(f"  {data['name']:<10}: i_max * difference = {data['i_max']} * {data['difference']:.4f} = {data['product']:.4f}")

    min_data = min(found_homologs_data, key=lambda x: x['product'])

    print("\n" + "=" * 60)
    print("                Final Answer Calculation")
    print("=" * 60)
    print(f"The minimum product is from the homolog '{min_data['name']}'.")
    print("\nThe final equation with the corresponding values is:")
    print(f"Minimum Product = {min_data['i_max']} * {min_data['difference']:.4f} = {min_data['product']:.4f}")

    global final_answer
    final_answer = min_data['product']


if __name__ == '__main__':
    # A global variable to hold the final numeric answer
    final_answer = None
    find_minimum_product()

