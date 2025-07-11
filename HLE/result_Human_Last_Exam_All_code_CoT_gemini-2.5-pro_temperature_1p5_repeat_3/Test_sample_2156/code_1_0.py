import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, GraphDescriptors

def solve():
    """
    Solves the problem by finding the minimum product of i_max and a difference of chi indices
    for formaldehyde homologs that meet a specific Geary autocorrelation criterion.
    """

    # 1. Define Homologs
    homologs = {}
    for i in range(1, 21):
        name = f'C{i}-Aldehyde'
        if i == 1:
            name = "Formaldehyde"
            smiles = "O=C"
        else:
            smiles = 'C' * i + '=O'
            if i == 2: name = "Acetaldehyde"
            if i == 3: name = "Propanal"
            if i == 4: name = "Butanal"
        homologs[smiles] = name
    
    # Sanderson electronegativity values
    en_map = {'H': 2.20, 'C': 2.55, 'O': 3.44}
    
    # Store results for filtered homologs
    found_homologs = []
    
    # 2. Loop through homologs to perform calculations
    max_lag = 30
    for smiles, name in homologs.items():
        # Create molecule object with explicit hydrogens for GATS calculation
        mol_with_hs = Chem.MolFromSmiles(smiles)
        mol_with_hs = Chem.AddHs(mol_with_hs)

        # Create list of Sanderson EN values
        try:
            en_values = [en_map[atom.GetSymbol()] for atom in mol_with_hs.GetAtoms()]
        except KeyError as e:
            print(f"Skipping {name}: Atom {e} not in electronegativity map.")
            continue
            
        # 2c. Calculate GATS vector
        gats_vector = Descriptors.GATS(mol_with_hs, weights=en_values, lag=max_lag)
        
        # We are interested in lags > 0
        gats_lags_1_plus = list(gats_vector[1:])
        
        if not gats_lags_1_plus:
            continue
            
        # 3. Filter based on max GATS value
        max_gats = max(gats_lags_1_plus)
        
        if 2 <= max_gats <= 3:
            # Find i_max (lag is 1-based index)
            i_max = gats_lags_1_plus.index(max_gats) + 1
            
            # 4. Calculate Chi indices on H-suppressed graph
            mol_no_hs = Chem.MolFromSmiles(smiles)
            avg_val_chi = GraphDescriptors.GetAverageValenceConnectivityI(mol_no_hs)
            avg_sim_chi = GraphDescriptors.GetAverageConnectivityI(mol_no_hs)
            diff_chi = avg_val_chi - avg_sim_chi

            # 5. Calculate product
            product = i_max * diff_chi
            
            # Store results
            found_homologs.append({
                "name": name,
                "i_max": i_max,
                "max_gats": max_gats,
                "diff_chi": diff_chi,
                "product": product
            })

    # 6. Find and report the minimum product
    if not found_homologs:
        print("No homologs were found that match the specified criteria.")
        return

    min_product_homolog = min(found_homologs, key=lambda x: x['product'])
    
    name = min_product_homolog['name']
    i_max_val = min_product_homolog['i_max']
    diff_chi_val = min_product_homolog['diff_chi']
    min_prod = min_product_homolog['product']
    
    print(f"The minimum product was found for the homolog: {name}")
    print("\nThe calculation is: i_max * (AverageValenceChi - AverageSimpleChi)")
    print("The values for this molecule are:")
    print(f"i_max = {i_max_val}")
    print(f"Difference of chi indices = {diff_chi_val:.4f}")
    print(f"Product = {i_max_val} * {diff_chi_val:.4f} = {min_prod:.4f}")

    print(f"\n<<<{min_prod}>>>")

# Execute the solver function
solve()