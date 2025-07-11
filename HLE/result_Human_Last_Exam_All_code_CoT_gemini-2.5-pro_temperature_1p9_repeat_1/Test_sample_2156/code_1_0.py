import sys
from rdkit import Chem
from mordred import Calculator, descriptors
from mordred.error import MissingValueBase
import warnings

# Suppress warnings from mordred for missing values, which are handled in the code
warnings.filterwarnings("ignore", category=UserWarning)

def solve_homolog_problem():
    """
    Analyzes aldehyde homologs to find the minimum product of i_max and a chi index difference.
    """
    # 1. Define the list of aldehyde homologs to be tested
    homolog_data = []
    base_names = [
        "Formaldehyde", "Acetaldehyde", "Propanal", "Butanal", "Pentanal",
        "Hexanal", "Heptanal", "Octanal", "Nonanal", "Decanal",
        "Undecanal", "Dodecanal", "Tridecanal", "Tetradecanal", "Pentadecanal"
    ]
    for i in range(len(base_names)):
        # SMILES string for C(n)H(2n+1)CHO is C...C=O with n 'C's
        smiles = "C" * i + "C=O"
        name = base_names[i]
        homolog_data.append({"name": name, "smiles": smiles})

    found_homologs = []
    
    # Create a calculator for Avp and Aap, which can be reused
    calc_chi = Calculator(descriptors.Avp, descriptors.Aap, ignore_3D=True)

    # 2. Loop through each homolog and perform calculations
    for homolog in homolog_data:
        mol = Chem.MolFromSmiles(homolog['smiles'])
        if mol is None:
            continue
        mol = Chem.AddHs(mol)

        # 3. Calculate Geary autocorrelation (weighted by S-ene) for lags 1-10
        geary_values = []
        max_lag = 10
        for lag in range(1, max_lag + 1):
            try:
                # A new calculator is created for each lag-dependent descriptor
                calc_geary = Calculator(descriptors.Geary(prop="S-ene", lag=lag), ignore_3D=True)
                val = next(iter(calc_geary(mol).values()))
                if isinstance(val, MissingValueBase):
                    # Use a very small number for non-calculable values to ignore them in max()
                    geary_values.append(float('-inf')) 
                else:
                    geary_values.append(float(val))
            except Exception:
                # Handle any other calculation errors
                geary_values.append(float('-inf'))
        
        # Proceed only if at least one Geary value was calculated
        if not any(v > float('-inf') for v in geary_values):
            continue

        g_max = max(geary_values)
        # The lag is the 1-based index of the max value
        i_max = geary_values.index(g_max) + 1  

        # 4. Filter homologs where 2 < G_max < 3
        if 2 < g_max < 3:
            try:
                # 5. Calculate Avp, Aap, and their difference
                chi_indices = calc_chi(mol)
                avp_val = chi_indices['Avp']
                aap_val = chi_indices['Aap']

                if isinstance(avp_val, MissingValueBase) or isinstance(aap_val, MissingValueBase):
                    continue

                diff = float(avp_val) - float(aap_val)
                product = i_max * diff

                # Store results for qualifying homologs
                found_homologs.append({
                    "name": homolog['name'],
                    "i_max": i_max,
                    "g_max": g_max,
                    "diff": diff,
                    "product": product
                })

            except Exception:
                # Skip if chi indices fail (e.g., for formaldehyde with 0 paths)
                continue

    # 6. Report findings and determine the minimum product
    if not found_homologs:
        print("No homologs were found that match the criteria (2 < G_max < 3).")
        return

    print("Found homologs with maximum Geary autocorrelation (G_max) between 2 and 3:")
    for h in found_homologs:
        print(f"- {h['name']}: i_max = {h['i_max']}, G_max = {h['g_max']:.4f}, Product = {h['product']:.4f}")
    
    print("\n" + "="*50 + "\n")
    
    min_homolog = min(found_homologs, key=lambda x: x['product'])
    min_product_val = min_homolog['product']
    
    print(f"The minimum product is found for {min_homolog['name']}.")
    print("The final calculation is the product of its i_max and the difference between its average valence and simple path chi indices:")
    # Print the final equation with all its numbers
    print(f"{min_homolog['i_max']} * {min_homolog['diff']:.4f} = {min_product_val:.4f}")
    
    # Final answer in the required format
    print(f"<<<{min_product_val:.4f}>>>")

# Execute the main function
solve_homolog_problem()