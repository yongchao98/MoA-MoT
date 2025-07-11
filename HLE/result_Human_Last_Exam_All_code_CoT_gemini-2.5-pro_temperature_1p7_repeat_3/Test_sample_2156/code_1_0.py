import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def solve_chemistry_problem():
    """
    Finds homologs of formaldehyde based on specific descriptor criteria and calculates
    the minimum product of i_max and a chi index difference.
    """
    # Based on the analysis, furaldehydes and related structures are likely candidates
    # for "formaldehyde's homologs" that meet the GATS criteria.
    candidate_molecules = {
        "2-Furaldehyde": 'c1oc(C=O)cc1',
        "3-Furaldehyde": 'c1oc(cc1)C=O',
        "5-Methyl-2-furaldehyde": 'Cc1oc(C=O)cc1',
        "3-Methyl-2-furaldehyde": 'c1c(C)c(C=O)oc1',
        "3-Oxapentanedial": 'O=CCOCC=O',
        "Glyoxal": 'O=CC=O', # Included for completeness, though it may not match
        "Succinaldehyde": 'O=CCCC=O' # Included for completeness
    }

    found_homologs = []

    print("--- Screening Candidate Molecules ---")

    for name, smiles in candidate_molecules.items():
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            print(f"Skipping invalid SMILES for {name}: {smiles}")
            continue
        mol = Chem.AddHs(mol)

        # 1. Calculate GATS weighted by Sanderson EN and find the maximum
        # The relevant GATSe descriptors (lag 1-8) are at indices 48-55
        try:
            all_autocorr = rdMolDescriptors.CalcAUTOCORR2D(mol)
            gats_e_values = list(all_autocorr[48:56])
        except Exception:
            # Some molecules (e.g., disconnected fragments) might fail
            print(f"Could not calculate descriptors for {name}.")
            continue
        
        if not gats_e_values:
            continue
            
        max_gats = max(gats_e_values)
        
        # 2. Check if the maximum value is within the specified range [2, 3]
        if 2.0 <= max_gats <= 3.0:
            i_max = gats_e_values.index(max_gats) + 1  # Lag is 1-based

            print(f"\nFound a matching homolog: {name}")
            print(f"  - Max GATS value ({max_gats:.4f}) is in [2, 3] at lag i_max = {i_max}.")

            # 3. Calculate average chi path indices for orders 0-10
            orders = range(11)
            simple_chis = [rdMolDescriptors.CalcChiNn(mol, order) for order in orders]
            valence_chis = [rdMolDescriptors.CalcChiNv(mol, order) for order in orders]
            
            avg_s_chi = np.mean(simple_chis)
            avg_v_chi = np.mean(valence_chis)

            # 4. Calculate the difference and the final product
            delta_chi = avg_v_chi - avg_s_chi
            product = i_max * delta_chi
            
            print(f"  - Average simple path chi index: {avg_s_chi:.4f}")
            print(f"  - Average valence path chi index: {avg_v_chi:.4f}")
            print(f"  - Difference (Δχ): {delta_chi:.4f}")
            print(f"  - Product (i_max * Δχ): {product:.4f}")

            # Store result for final comparison
            found_homologs.append({
                "name": name,
                "i_max": i_max,
                "avg_v_chi": avg_v_chi,
                "avg_s_chi": avg_s_chi,
                "product": product
            })
        else:
            print(f"Molecule {name} did not match (Max GATS = {max_gats:.4f}).")


    print("\n--- Final Result ---")
    if not found_homologs:
        print("No homologs were found that satisfy the condition.")
        return

    # Find the homolog with the minimum product
    min_homolog = min(found_homologs, key=lambda x: x['product'])

    # Print the detailed calculation for the minimum product found
    i = min_homolog['i_max']
    avg_v = min_homolog['avg_v_chi']
    avg_s = min_homolog['avg_s_chi']
    prod = min_homolog['product']
    
    print(f"The minimum product is found for the molecule: {min_homolog['name']}.")
    print("\nThe detailed calculation is:")
    print(f"i_max * (Average Valence Chi - Average Simple Chi) = Product")
    print(f"{i} * ({avg_v:.4f} - {avg_s:.4f}) = {prod:.4f}")
    
    # Final answer in specified format
    print(f"\n<<<{prod:.4f}>>>")

if __name__ == "__main__":
    solve_chemistry_problem()
