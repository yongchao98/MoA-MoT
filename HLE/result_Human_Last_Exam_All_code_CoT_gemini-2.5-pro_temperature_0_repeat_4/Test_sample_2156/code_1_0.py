import sys
try:
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    import numpy as np
except ImportError:
    print("Error: RDKit or NumPy is not installed. Please install them to run this script.")
    print("You can install them using pip: pip install rdkit-pypi numpy")
    sys.exit(1)

def calculate_geary_sane(mol):
    """
    Calculates Geary autocorrelation for all possible lags, weighted by
    Sanderson electronegativity.
    The Geary Autocorrelation G(d) for a lag d is the average squared difference
    of the property for all atom pairs at that topological distance.
    """
    # Sanderson electronegativity values for C, O, H.
    sane = {'C': 2.746, 'O': 3.654, 'H': 2.592}
    
    atoms = mol.GetAtoms()
    N = len(atoms)
    if N <= 1:
        return {}

    # Create a list of electronegativity values for each atom.
    try:
        w_list = [sane[atom.GetSymbol()] for atom in atoms]
    except KeyError as e:
        print(f"Error: Atom type {e} not found in Sanderson electronegativity dictionary.")
        return {}

    # Get the topological distance matrix.
    dist_matrix = Chem.GetDistanceMatrix(mol)
    
    geary_values = {}
    max_lag = int(np.max(dist_matrix))

    # Iterate through all possible lags (from 1 to the max distance).
    for lag in range(1, max_lag + 1):
        numerator_sum = 0
        pair_count = 0
        # Iterate over all unique pairs of atoms.
        for i in range(N):
            for j in range(i + 1, N):
                if dist_matrix[i, j] == lag:
                    numerator_sum += (w_list[i] - w_list[j])**2
                    pair_count += 1
        
        if pair_count > 0:
            geary_d = numerator_sum / pair_count
            geary_values[lag] = geary_d
            
    return geary_values

def calculate_avg_chi_diff(mol):
    """
    Calculates the difference between the average valence path chi index and
    the average simple path chi index for orders 0 to 10.
    """
    vpc_list = []
    spc_list = []
    num_orders = 11  # Corresponds to orders 0, 1, ..., 10

    for n in range(num_orders):
        # Calculate valence path chi index (sensitive to heteroatoms).
        vpc_list.append(rdMolDescriptors.CalcChiNv(mol, n))
        # Calculate simple path chi index (based on connectivity only).
        spc_list.append(rdMolDescriptors.CalcChiNn(mol, n))
        
    avg_vpc = sum(vpc_list) / len(vpc_list)
    avg_spc = sum(spc_list) / len(spc_list)
    
    return avg_vpc - avg_spc

def solve_task():
    """
    Main function to find the minimum product based on the problem's criteria.
    """
    # Define the homologous series of aldehydes using SMILES strings.
    aldehydes = {
        "Formaldehyde": "C=O", "Acetaldehyde": "CC=O", "Propanal": "CCC=O",
        "Butanal": "CCCC=O", "Pentanal": "CCCCC=O", "Hexanal": "CCCCCC=O",
        "Heptanal": "CCCCCCC=O", "Octanal": "CCCCCCCC=O", "Nonanal": "CCCCCCCCC=O",
        "Decanal": "CCCCCCCCCC=O", "Undecanal": "CCCCCCCCCCC=O", "Dodecanal": "CCCCCCCCCCCC=O"
    }

    found_homologs = []

    # Process each aldehyde in the series.
    for name, smi in aldehydes.items():
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)

        geary_values = calculate_geary_sane(mol)

        if not geary_values:
            continue

        # Find the maximum Geary value and its corresponding lag (i_max).
        i_max = max(geary_values, key=geary_values.get)
        max_geary = geary_values[i_max]

        # Filter homologs where max_geary is between 2 and 3.
        if 2 <= max_geary <= 3:
            chi_diff = calculate_avg_chi_diff(mol)
            product = i_max * chi_diff
            
            found_homologs.append({
                "name": name,
                "i_max": i_max,
                "chi_diff": chi_diff,
                "product": product
            })

    if not found_homologs:
        print("No homologs were found with a maximum Geary autocorrelation between 2 and 3.")
        return

    # Find the homolog with the minimum product.
    best_homolog = min(found_homologs, key=lambda x: x['product'])
    
    final_imax = best_homolog['i_max']
    final_chi_diff = best_homolog['chi_diff']
    final_product = best_homolog['product']

    # Output the numbers involved in the final calculation.
    print(f"The homolog that yields the minimum product is {best_homolog['name']}.")
    print(f"The lag for its maximum Geary autocorrelation is i_max = {final_imax}.")
    print(f"The difference between its average chi indices is {final_chi_diff:.6f}.")
    print(f"The final equation is: {final_imax} * {final_chi_diff:.6f} = {final_product:.6f}")
    
    # Output the final answer in the specified format.
    print(f"<<<{final_product:.6f}>>>")

if __name__ == "__main__":
    solve_task()