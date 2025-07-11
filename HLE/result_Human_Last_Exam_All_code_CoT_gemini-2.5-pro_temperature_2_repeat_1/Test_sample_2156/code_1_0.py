import sys

try:
    from rdkit import Chem
    from mordred import Calculator, descriptors, error
except ImportError:
    print("Error: RDKit or mordred library is not installed.", file=sys.stderr)
    print("Please install them using: pip install rdkit-pypi mordred", file=sys.stderr)
    sys.exit(1)

def find_minimum_product():
    """
    Finds the minimum product of i_max and a path chi index difference for formaldehyde homologs
    that satisfy a specific Geary autocorrelation criterion.
    """
    # Initialize the Mordred descriptor calculator. We only need the GearyAuto and Path sets.
    # ignore_3D=True is used as we are working with 2D representations (SMILES).
    calc = Calculator(descriptors.GearyAuto, descriptors.Path, ignore_3D=True)

    # Define the formaldehyde homologs to be analyzed using their common names and SMILES strings.
    homologs = {
        'Formaldehyde': 'C=O',
        'Acetaldehyde': 'CC=O',
        'Propionaldehyde': 'CCC=O',
        'Butyraldehyde': 'CCCC=O',
        'Pentanal': 'CCCCC=O',
        'Hexanal': 'CCCCCC=O',
        'Heptanal': 'CCCCCCC=O',
        'Octanal': 'CCCCCCCC=O',
        'Nonanal': 'CCCCCCCCC=O',
        'Decanal': 'CCCCCCCCCC=O'
    }

    # This list will store the calculated results for homologs that meet the criteria.
    found_homologs = []

    # These are the specific descriptor names we're interested in.
    gats_descriptor_names = [f'GATS{i}e' for i in range(1, 9)]

    # Iterate over the defined homologs.
    for name, smiles in homologs.items():
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            # This is a safeguard against invalid SMILES strings.
            continue
        
        # Add explicit hydrogen atoms, as they are required for accurate descriptor calculation.
        mol = Chem.AddHs(mol)
        
        # Calculate descriptors. The 'quiet=True' argument suppresses per-descriptor error logging.
        all_descriptors = calc(mol, quiet=True)

        # Extract the GATS values.
        gats_values = [all_descriptors[d_name] for d_name in gats_descriptor_names]

        # If any descriptor failed to calculate, skip this molecule.
        if any(isinstance(val, error.Error) for val in gats_values):
            continue

        # Find the maximum GATS value and its corresponding lag 'i'.
        max_gats_value = -float('inf')
        i_max = 0
        for i, value in enumerate(gats_values, 1):
            if value > max_gats_value:
                max_gats_value = value
                i_max = i

        # Apply the filter: maximum GATS value must be between 2 and 3.
        if 2.0 < max_gats_value < 3.0:
            # Retrieve the necessary path chi indices.
            avg_valence_path = all_descriptors['Xpaval']
            avg_simple_path = all_descriptors['Xp_avg']

            # Ensure path chi indices were calculated successfully.
            if isinstance(avg_valence_path, error.Error) or isinstance(avg_simple_path, error.Error):
                continue
            
            # Calculate the difference and the final product.
            delta_chi = avg_valence_path - avg_simple_path
            product = i_max * delta_chi
            
            # Store the results for this valid homolog.
            found_homologs.append({
                "name": name,
                "i_max": i_max,
                "delta_chi": delta_chi,
                "product": product
            })

    if not found_homologs:
        print("No homologs were found that met the specified criteria.")
        return

    # Determine which of the found homologs has the minimum product.
    min_result_homolog = min(found_homologs, key=lambda x: x['product'])

    # Present the final result as requested.
    print(f"The minimum product is found for the homolog: {min_result_homolog['name']}")
    print("The final equation is:")
    
    i = min_result_homolog['i_max']
    d_chi = min_result_homolog['delta_chi']
    prod = min_result_homolog['product']
    
    # Output each number in the final equation.
    print(f"{i} * {d_chi} = {prod}")
    
    # Finally, output the answer in the specified format.
    print(f"\n<<<{prod}>>>")

if __name__ == "__main__":
    find_minimum_product()