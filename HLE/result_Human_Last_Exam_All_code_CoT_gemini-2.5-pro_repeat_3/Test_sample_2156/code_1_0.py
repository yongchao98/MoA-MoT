import sys
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors, GraphDescriptors
except ImportError:
    print("RDKit is not installed. Please install it using 'pip install rdkit-pypi'")
    sys.exit(1)
import numpy as np

def solve_chemoinformatics_problem():
    """
    Solves the chemoinformatics problem as described.
    1. Finds aldehyde homologs with max GATS value between 2 and 3.
    2. For these homologs, calculates a product involving i_max and chi indices.
    3. Determines the minimum product among them.
    """
    # Generate formaldehyde's homologous series (n-aldehydes)
    homologs_smiles = ['C' * n + '=O' for n in range(1, 11)]
    homolog_names = [
        "Formaldehyde", "Acetaldehyde", "Propanal", "Butanal", "Pentanal",
        "Hexanal", "Heptanal", "Octanal", "Nonanal", "Decanal"
    ]

    # Lists of chi index functions from RDKit for orders 0-10
    chi_s_funcs = [
        GraphDescriptors.Chi0p, GraphDescriptors.Chi1p, GraphDescriptors.Chi2p,
        GraphDescriptors.Chi3p, GraphDescriptors.Chi4p, GraphDescriptors.Chi5p,
        GraphDescriptors.Chi6p, GraphDescriptors.Chi7p, GraphDescriptors.Chi8p,
        GraphDescriptors.Chi9p, GraphDescriptors.Chi10p
    ]
    chi_v_funcs = [
        GraphDescriptors.Chi0pv, GraphDescriptors.Chi1pv, GraphDescriptors.Chi2pv,
        GraphDescriptors.Chi3pv, GraphDescriptors.Chi4pv, GraphDescriptors.Chi5pv,
        GraphDescriptors.Chi6pv, GraphDescriptors.Chi7pv, GraphDescriptors.Chi8pv,
        GraphDescriptors.Chi9pv, GraphDescriptors.Chi10pv
    ]

    found_homologs_results = []

    # Process each homolog
    for name, smiles in zip(homolog_names, homologs_smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        
        # Add explicit hydrogens for correct descriptor calculation
        mol_h = Chem.AddHs(mol)

        # Calculate 2D Autocorrelation descriptors
        autocorr_descriptors = rdMolDescriptors.CalcAUTOCORR2D(mol_h)
        # GATS descriptors (lags 1-8) are at indices 96-103
        gats_values = list(autocorr_descriptors[96:104])
        
        max_gats = max(gats_values)
        
        # Filter based on max GATS value
        if 2 < max_gats < 3:
            # lag = index + 1
            i_max = gats_values.index(max_gats) + 1

            # Calculate simple path chi indices and their average
            chi_s = [f(mol_h) for f in chi_s_funcs]
            avg_X_s = sum(chi_s) / len(chi_s)

            # Calculate valence path chi indices and their average
            chi_v = [f(mol_h) for f in chi_v_funcs]
            avg_X_v = sum(chi_v) / len(chi_v)
            
            # Calculate the difference and the final product
            diff_chi = avg_X_v - avg_X_s
            product = i_max * diff_chi

            found_homologs_results.append({
                "name": name,
                "smiles": smiles,
                "i_max": i_max,
                "avg_X_v": avg_X_v,
                "avg_X_s": avg_X_s,
                "product": product
            })

    if not found_homologs_results:
        print("No homologs found with max Geary autocorrelation between 2 and 3.")
        return

    print("Found homologs satisfying the condition. Calculating the product for each:")
    print("-" * 70)

    min_product = float('inf')
    min_result = None

    for result in found_homologs_results:
        print(f"Homolog: {result['name']} ({result['smiles']})")
        # Output each number in the final equation
        print(f"  i_max = {result['i_max']}")
        print(f"  Average Valence Path Chi Index (avg_X_v) = {result['avg_X_v']}")
        print(f"  Average Simple Path Chi Index (avg_X_s)  = {result['avg_X_s']}")
        print(f"  Product = {result['i_max']} * ({result['avg_X_v']} - {result['avg_X_s']}) = {result['product']}")
        print("-" * 70)

        if result['product'] < min_product:
            min_product = result['product']
            min_result = result
            
    print(f"\nThe minimum product is {min_product}, which corresponds to {min_result['name']}.")
    
    # Final answer in the required format
    print(f"\n<<<{min_product}>>>")


if __name__ == "__main__":
    solve_chemoinformatics_problem()