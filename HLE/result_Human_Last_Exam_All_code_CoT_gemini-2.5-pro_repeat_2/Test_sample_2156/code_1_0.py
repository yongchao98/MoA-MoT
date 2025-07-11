import sys
import io
from rdkit import Chem
from mordred import Calculator, descriptors

# Suppress potential RDKit warnings for cleaner output
from rdkit import rdBase
rdBase.DisableLog('rdApp.warning')

def solve_homolog_problem():
    """
    Finds the minimum product of i_max and the difference between average 
    valence and simple path chi indices for formaldehyde homologs satisfying a 
    specific Geary autocorrelation condition.
    """
    
    # Step 1: Define the homologous series of aldehydes
    aldehyde_names = [
        "Formaldehyde", "Acetaldehyde", "Propanal", "Butanal",
        "Pentanal", "Hexanal", "Heptanal", "Octanal"
    ]
    aldehyde_smiles = [
        "C=O", "CC=O", "CCC=O", "CCCC=O", "CCCCC=O",
        "CCCCCC=O", "CCCCCCC=O", "CCCCCCCC=O"
    ]

    # Step 2: Set up the mordred calculator for required descriptors
    # GATS with Sanderson Electronegativity (se) for lags 1 to 8
    gats_descriptor_names = [f"GATS{i}se" for i in range(1, 9)]
    gats_descriptors = [getattr(descriptors.Geary, name) for name in gats_descriptor_names]
    
    chi_descriptors = [
        descriptors.Path.ChiPath.AvgVPath,
        descriptors.Path.ChiPath.AvgPath
    ]
    
    # Create the mordred calculator instance
    calc = Calculator(gats_descriptors + chi_descriptors, ignore_3D=True)

    found_homologs_data = []

    # Iterate through each aldehyde in the series
    for name, smiles in zip(aldehyde_names, aldehyde_smiles):
        mol = Chem.MolFromSmiles(smiles)
        # Add hydrogens, as they are crucial for these descriptors
        mol_h = Chem.AddHs(mol)

        # Calculate all descriptors for the molecule
        results = calc(mol_h)

        # Step 3: Find the maximum Geary autocorrelation value and its lag (i_max)
        gats_values = {}
        for i in range(1, 9):
            try:
                # Descriptor values can be errors, so we handle them by converting to float
                val = float(results[f"GATS{i}se"])
                gats_values[i] = val
            except (KeyError, TypeError, ValueError):
                # This lag might not be applicable for the molecule (e.g., max path is shorter)
                continue
        
        if not gats_values:
            continue

        i_max = max(gats_values, key=gats_values.get)
        gats_max = gats_values[i_max]

        # Filter homologs based on the condition 2 < GATS_max < 3
        if 2 < gats_max < 3:
            # Step 4: For found homologs, calculate the target product
            try:
                avg_vpath = float(results["AvgVPath"])
                avg_path = float(results["AvgPath"])
            except (KeyError, TypeError, ValueError):
                # Skip if chi indices cannot be calculated
                continue

            product = i_max * (avg_vpath - avg_path)
            
            found_homologs_data.append({
                "name": name,
                "i_max": i_max,
                "gats_max": gats_max,
                "avg_vpath": avg_vpath,
                "avg_path": avg_path,
                "product": product
            })

    # Step 5: Determine the minimum product among the found homologs
    if not found_homologs_data:
        print("No homologs were found that satisfy the condition 2 < max GATS < 3.")
        return

    min_data = min(found_homologs_data, key=lambda x: x['product'])
    
    final_name = min_data['name']
    final_i_max = min_data['i_max']
    final_avg_vpath = min_data['avg_vpath']
    final_avg_path = min_data['avg_path']
    final_product = min_data['product']

    print("Homologs satisfying the condition (2 < max GATS < 3):")
    for data in found_homologs_data:
        print(f"- {data['name']} (max GATS = {data['gats_max']:.4f} at lag {data['i_max']})")
    print("-" * 50)

    print(f"The homolog with the minimum product is '{final_name}'.")
    print("The final calculation is the product of its i_max and the difference between its average valence path chi index and average simple path chi index.")
    
    # Final output showing the numbers in the equation
    print(f"\nFinal Equation: {final_i_max} * ({final_avg_vpath:.4f} - {final_avg_path:.4f}) = {final_product:.4f}")

    # The required answer format
    print(f"\n<<<{final_product}>>>")

# Execute the function
solve_homolog_problem()