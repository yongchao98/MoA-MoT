import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors

def find_minimum_product():
    """
    This script solves the user's request by:
    1. Identifying formaldehyde homologs (aldehydes).
    2. Filtering them based on the maximum Geary autocorrelation value (weighted by Sanderson EN).
    3. For the filtered homologs, calculating the product of i_max and the difference between
       average valence and simple path chi indices.
    4. Determining the minimum of these products.
    """
    # Step 1: Define the homologs of formaldehyde (aldehydes)
    aldehyde_names = [
        "Formaldehyde", "Acetaldehyde", "Propanal", "Butanal",
        "Pentanal", "Hexanal", "Heptanal", "Octanal", "Nonanal",
        "Decanal", "Undecanal", "Dodecanal"
    ]
    aldehyde_smiles = [
        "C=O", "CC=O", "CCC=O", "CCCC=O", "CCCCC=O", "CCCCCC=O",
        "CCCCCCC=O", "CCCCCCCC=O", "CCCCCCCCC=O", "CCCCCCCCCC=O",
        "CCCCCCCCCCC=O", "CCCCCCCCCCCC=O"
    ]

    # Step 2: Find homologs that meet the Geary autocorrelation criteria
    found_homologs = []
    print("--- Step 1: Finding homologs with max Geary Autocorrelation between 2 and 3 ---")
    
    for name, smiles in zip(aldehyde_names, aldehyde_smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            
            # Add hydrogens for the Geary calculation
            mol_h = Chem.AddHs(mol)
            num_atoms = mol_h.GetNumAtoms()

            if num_atoms <= 1:
                continue

            dist_matrix = Chem.GetDistanceMatrix(mol_h)
            max_lag = int(np.max(dist_matrix))

            geary_values = {i: Descriptors.GearyAutocorrelation(mol_h, lag=i, prop='sanderson') 
                            for i in range(1, max_lag + 1)}

            if not geary_values:
                continue

            i_max = max(geary_values, key=geary_values.get)
            G_max = geary_values[i_max]

            print(f"Molecule: {name:<12} | G_max = {G_max:.4f} @ lag i_max = {i_max}")

            # Check if the homolog meets the condition 2 <= G_max <= 3
            if 2 <= G_max <= 3:
                print(f"  -> MATCH FOUND: {name}")
                found_homologs.append({"name": name, "smiles": smiles, "i_max": i_max})
        except Exception as e:
            print(f"An error occurred while processing {name}: {e}")
    
    if not found_homologs:
        print("\nNo homologs were found matching the specified criteria.")
        return

    # Step 3: Calculate the product for each found homolog
    results = []
    print("\n--- Step 2: Calculating indices and product for found homologs ---")

    for homolog in found_homologs:
        name, smiles, i_max = homolog["name"], homolog["smiles"], homolog["i_max"]
        
        # Chi indices are calculated on the heavy-atom graph
        mol = Chem.MolFromSmiles(smiles)
        num_heavy_atoms = mol.GetNumAtoms()
        max_order = min(10, num_heavy_atoms - 1)

        orders = list(range(max_order + 1))
        simple_chi = [Descriptors.ChiNn(mol, n) for n in orders]
        valence_chi = [Descriptors.ChiNv(mol, n) for n in orders]
        
        avg_X_s = np.mean(simple_chi)
        avg_X_v = np.mean(valence_chi)
        
        difference = avg_X_v - avg_X_s
        product = i_max * difference
        
        results.append({
            "name": name, "product": product, "i_max": i_max,
            "avg_X_v": avg_X_v, "avg_X_s": avg_X_s, "diff": difference
        })
        
        print(f"\nCalculating for: {name}")
        print(f"  i_max = {i_max}")
        print(f"  Avg Simple Path Chi (orders 0-{max_order})   = {avg_X_s:.4f}")
        print(f"  Avg Valence Path Chi (orders 0-{max_order}) = {avg_X_v:.4f}")
        print(f"  Product = {i_max} * ({avg_X_v:.4f} - {avg_X_s:.4f}) = {product:.4f}")

    # Step 4: Find and report the minimum product
    if not results:
        print("\nCould not calculate the final product for any found homolog.")
        return

    min_result = min(results, key=lambda x: x["product"])
    
    print("\n--- Step 3: Determining the minimum product ---")
    print(f"The minimum product is for the molecule: {min_result['name']}")
    
    # Final formatted output
    print("\nFinal Equation for the Minimum Product:")
    print(f"min_product = i_max * (avg_valence_chi - avg_simple_chi)")
    print(f"min_product = {min_result['i_max']} * ({min_result['avg_X_v']:.4f} - {min_result['avg_X_s']:.4f})")
    print(f"min_product = {min_result['i_max']} * {min_result['diff']:.4f}")
    
    final_answer = min_result['product']
    print(f"min_product = {final_answer:.4f}")
    
    print(f"\n<<<{final_answer:.4f}>>>")

if __name__ == '__main__':
    find_minimum_product()