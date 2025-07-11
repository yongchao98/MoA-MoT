import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def find_minimum_descriptor_product():
    """
    This function solves the described cheminformatics problem.
    It identifies homologs of formaldehyde that satisfy a specific Geary autocorrelation criterion,
    then calculates a product based on their chi indices and i_max values,
    and finally determines the minimum possible value of this product.
    """
    # Define Sanderson electronegativity values for the relevant atoms.
    sanderson_en = {'C': 2.746, 'O': 3.654, 'H': 2.592}
    # A common maximum order for calculating average Chi indices.
    MAX_CHI_ORDER = 10
    
    # Initialize variables to track the minimum product and the details of the molecule that produces it.
    min_product_so_far = float('inf')
    min_details = {}

    # Iterate through the aldehyde homologous series from 1 to 20 carbons.
    for n_carbons in range(1, 21):
        # Generate the SMILES string for the linear aldehyde.
        smiles = 'C' * (n_carbons - 1) + 'C=O'
        
        # Create an RDKit molecule object from SMILES for Chi index calculations (heavy atoms only).
        mol_heavy = Chem.MolFromSmiles(smiles)
        if mol_heavy is None:
            continue
        
        # Create a second molecule object with explicit hydrogens for GATS calculation.
        mol_with_H = Chem.AddHs(mol_heavy)

        # Assign Sanderson electronegativity to each atom as its weight.
        try:
            weights = [sanderson_en[atom.GetSymbol()] for atom in mol_with_H.GetAtoms()]
        except KeyError:
            # Skip this molecule if it contains an atom for which we don't have an EN value.
            continue

        # Determine the maximum possible lag (topological distance) in the molecule.
        dist_matrix = Chem.GetDistanceMatrix(mol_with_H)
        max_lag = int(np.max(dist_matrix)) if dist_matrix.size > 0 else 0
        if max_lag == 0:
            continue

        # Calculate Geary autocorrelation for all lags from 1 to the maximum.
        gats_values = [rdMolDescriptors.GATS(mol_with_H, weights, lag=i) for i in range(1, max_lag + 1)]

        # Find the maximum GATS value and the corresponding lag (i_max).
        if not gats_values:
            continue
        max_gats = np.max(gats_values)
        i_max = np.argmax(gats_values) + 1  # Lags are 1-based, index is 0-based.

        # Filter for homologs where the maximum GATS value is between 2 and 3.
        if 2 < max_gats < 3:
            # For filtered homologs, calculate average chi path indices.
            chi_p_vals = [rdMolDescriptors.CalcChiNpath(mol_heavy, i) for i in range(MAX_CHI_ORDER + 1)]
            avg_chi_p = np.mean(chi_p_vals)

            chi_vp_vals = [rdMolDescriptors.CalcChiNpathv(mol_heavy, i) for i in range(MAX_CHI_ORDER + 1)]
            avg_chi_vp = np.mean(chi_vp_vals)
            
            # Calculate the difference between the average chi indices.
            diff = avg_chi_vp - avg_chi_p

            # Calculate the final product.
            product = i_max * diff

            # If this product is the smallest found so far, store its value and details.
            if product < min_product_so_far:
                min_product_so_far = product
                min_details = {
                    'smiles': smiles,
                    'i_max': i_max,
                    'avg_chi_vp': avg_chi_vp,
                    'avg_chi_p': avg_chi_p,
                    'product': product
                }

    # After checking all homologs, print the results.
    if not min_details:
        print("No formaldehyde homologs were found with a maximum Geary autocorrelation between 2 and 3.")
    else:
        print(f"The homolog with the minimum product is {min_details['smiles']}.")
        print("The calculation is: product = i_max * (average_valence_path_chi - average_simple_path_chi)")
        print("\nThe numbers for the final equation that result in the minimum product are:")
        print(f"i_max = {min_details['i_max']}")
        print(f"average_valence_path_chi = {min_details['avg_chi_vp']}")
        print(f"average_simple_path_chi = {min_details['avg_chi_p']}")

        # The final answer in the requested format.
        print(f"\n<<<{min_details['product']}>>>")

if __name__ == '__main__':
    find_minimum_descriptor_product()