import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import numpy as np

def solve_aldehyde_properties():
    """
    Finds the minimum product of i_max and the difference between average chi indices
    for formaldehyde homologs that meet a specific Geary autocorrelation criterion.
    """
    # Step 1: Generate Homologs (aldehydes C1 to C20)
    max_carbons = 20
    # C1 is Formaldehyde, C2 is Acetaldehyde, etc.
    aldehydes = {f'C{n} Aldehyde': smiles for n, smiles in enumerate(
        ['C=O'] + ['C' * i + 'C=O' for i in range(1, max_carbons)], 1)}

    # Step 2: Prepare for descriptor calculation
    # Get the list of all 2D autocorrelation descriptor names from RDKit
    all_descriptor_names = list(rdMolDescriptors.AUTOCORR2D_DESCRIPTORS)
    # Find the indices corresponding to Geary autocorrelation with Sanderson electronegativity
    # weighting (GATSde) for lags d=1 to 8.
    geary_en_indices = [all_descriptor_names.index(f"GATS{i}e") for i in range(1, 9)]

    found_homologs_data = []

    # Step 3 & 4: Iterate through homologs, calculate descriptors, and filter
    for name, smiles in aldehydes.items():
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue
        # Add hydrogens, as they are required for some descriptor calculations
        mol = Chem.AddHs(mol)

        # Calculate all 2D autocorrelation descriptors
        all_descriptors = rdMolDescriptors.CalcAUTOCORR2D(mol)
        
        # Extract the 8 Geary EN values
        geary_en_values = np.array([all_descriptors[i] for i in geary_en_indices])
        
        # Find the maximum Geary value and the corresponding lag (i_max)
        max_geary_value = np.max(geary_en_values)
        # The lag is the 1-based index of the maximum value
        i_max = np.argmax(geary_en_values) + 1

        # Check if the homolog meets the criterion
        if 2 <= max_geary_value <= 3:
            # Step 5: Calculate Chi Indices for the found homolog
            # We average over orders 0 to 10 (11 values total)
            num_orders = 11
            
            # Average simple path chi index
            chi_p_values = [rdMolDescriptors.CalcChiNn(mol, order) for order in range(num_orders)]
            avg_chi_p = sum(chi_p_values) / num_orders
            
            # Average valence path chi index
            chi_vp_values = [rdMolDescriptors.CalcChiNv(mol, order) for order in range(num_orders)]
            avg_chi_vp = sum(chi_vp_values) / num_orders

            # Step 6: Calculate the final product
            diff_chi = avg_chi_vp - avg_chi_p
            product = i_max * diff_chi

            found_homologs_data.append({
                "name": name,
                "smiles": smiles,
                "max_geary_value": max_geary_value,
                "i_max": i_max,
                "avg_chi_p": avg_chi_p,
                "avg_chi_vp": avg_chi_vp,
                "diff_chi": diff_chi,
                "product": product
            })

    # Step 7: Find the minimum product and report the results
    if not found_homologs_data:
        print("No homologs were found that satisfy the criterion.")
        return

    print("Found homologs satisfying the condition (Max Geary Autocorrelation between 2 and 3):")
    for data in found_homologs_data:
        print(f"- {data['name']} ({data['smiles']}): Max Geary = {data['max_geary_value']:.4f}, Product = {data['product']:.4f}")
    
    print("\n" + "="*80)

    # Find the homolog with the minimum product
    min_data = min(found_homologs_data, key=lambda x: x['product'])

    print(f"The minimum product is found for {min_data['name']}.")
    print("Calculation details for this homolog:")
    print(f"  - Lag with maximum Geary autocorrelation (i_max): {min_data['i_max']}")
    print(f"  - Average valence path chi index (χ_vp): {min_data['avg_chi_vp']:.4f}")
    print(f"  - Average simple path chi index (χ_p): {min_data['avg_chi_p']:.4f}")
    
    # Final equation printing as requested
    print("\nFinal Equation:")
    print(f"{min_data['i_max']} * ({min_data['avg_chi_vp']:.4f} - {min_data['avg_chi_p']:.4f}) = {min_data['product']:.4f}")
    print(f"{min_data['i_max']} * {min_data['diff_chi']:.4f} = {min_data['product']:.4f}")
    
    final_answer = min_data['product']
    print(f"\n<<<ANSER>>>")
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    solve_aldehyde_properties()