import numpy as np
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdPeriodicTable

def solve_chemistry_problem():
    """
    Solves the computational chemistry problem as described.
    1.  Generates aldehyde homologs.
    2.  For each, finds the max Geary autocorrelation (weighted by Sanderson EN) and its lag (i_max).
    3.  Filters homologs where 2 <= max_geary <= 3.
    4.  For filtered homologs, calculates avg simple and valence chi path indices.
    5.  Calculates the product of i_max and the difference in chi indices.
    6.  Finds the homolog with the minimum product.
    """
    # Generate aldehyde homologous series (formaldehyde to dodecanal)
    smiles_list = ['C=O'] + ['C' * n + 'C=O' for n in range(1, 12)]
    homolog_names = [
        'Formaldehyde', 'Acetaldehyde', 'Propanal', 'Butanal', 'Pentanal',
        'Hexanal', 'Heptanal', 'Octanal', 'Nonanal', 'Decanal',
        'Undecanal', 'Dodecanal'
    ]

    candidate_homologs = []

    # Maximum lag and chi order to consider
    max_lag = 10
    max_chi_order = 10

    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue

        num_heavy_atoms = mol.GetNumHeavyAtoms()

        # 1. Get Sanderson Electronegativity for each heavy atom
        sanderson_en = [rdPeriodicTable.GetEN(atom.GetAtomicNum()) for atom in mol.GetAtoms()]

        # 2. Find the maximum Geary autocorrelation and its corresponding lag (i_max)
        max_geary_val = -float('inf')
        i_max = -1

        for lag in range(1, max_lag + 1):
            if lag >= num_heavy_atoms:
                break
            # Geary autocorrelation weighted by Sanderson electronegativities
            geary_val = rdMolDescriptors.GearyAuto(mol, lag, sanderson_en)
            if geary_val > max_geary_val:
                max_geary_val = geary_val
                i_max = lag

        # 3. Filter homologs where max_geary is between 2 and 3
        if 2 <= max_geary_val <= 3:
            # 4. Calculate average chi indices
            chi_s_values = []
            chi_v_values = []
            for order in range(max_chi_order + 1):
                chi_s_values.append(rdMolDescriptors.CalcChiNn(mol, order))
                chi_v_values.append(rdMolDescriptors.CalcChiNv(mol, order))

            avg_chi_s = np.mean(chi_s_values)
            avg_chi_v = np.mean(chi_v_values)
            
            # 5. Calculate the product
            product = i_max * (avg_chi_v - avg_chi_s)

            candidate_homologs.append({
                'name': homolog_names[i],
                'smiles': smi,
                'i_max': i_max,
                'avg_chi_s': avg_chi_s,
                'avg_chi_v': avg_chi_v,
                'product': product,
                'max_geary': max_geary_val
            })

    # 6. Find the homolog with the minimum product
    if not candidate_homologs:
        print("No homologs found matching the criteria.")
        return

    min_product_homolog = min(candidate_homologs, key=lambda x: x['product'])

    # Print the detailed result for the best homolog
    print(f"Homolog with minimum product: {min_product_homolog['name']} ({min_product_homolog['smiles']})")
    print(f"  - Max Geary Autocorrelation: {min_product_homolog['max_geary']:.4f} (within [2, 3])")
    print("\nFinal Equation Components:")
    
    i_val = min_product_homolog['i_max']
    avg_v_val = min_product_homolog['avg_chi_v']
    avg_s_val = min_product_homolog['avg_chi_s']
    prod_val = min_product_homolog['product']
    
    print(f"  - i_max = {i_val}")
    print(f"  - Average Valence Path Chi Index = {avg_v_val:.4f}")
    print(f"  - Average Simple Path Chi Index  = {avg_s_val:.4f}")
    print(f"\nCalculation: {i_val} * ({avg_v_val:.4f} - {avg_s_val:.4f}) = {prod_val:.4f}")

    # Print the final answer in the required format
    print(f"\n<<<{prod_val:.4f}>>>")

# Execute the function
solve_chemistry_problem()