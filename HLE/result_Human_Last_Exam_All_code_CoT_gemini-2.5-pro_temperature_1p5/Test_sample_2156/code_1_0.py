import sys
import subprocess

# Install required packages if they are not already installed
try:
    import pandas as pd
    from rdkit import Chem
    from mordred import Calculator, descriptors
except ImportError:
    print("One or more required libraries are not installed. Attempting to install them...")
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas", "rdkit", "mordred"])
        print("Libraries installed successfully. Please run the script again.")
    except Exception as e:
        print(f"Failed to install libraries. Please install them manually using 'pip install pandas rdkit mordred'. Error: {e}")
    sys.exit()

def calculate_properties():
    """
    This script finds formaldehyde homologs that meet specific chemical criteria
    and calculates a target value based on their molecular descriptors.
    """
    print("Starting the analysis...")
    print("Step 1: Defining the set of formaldehyde homologs (aldehydes C1-C12).")
    aldehydes = {
        'formaldehyde': 'C=O',
        'acetaldehyde': 'CC=O',
        'propanal': 'CCC=O',
        'butanal': 'CCCC=O',
        'pentanal': 'CCCCC=O',
        'hexanal': 'CCCCCC=O',
        'heptanal': 'CCCCCCC=O',
        'octanal': 'CCCCCCCC=O',
        'nonanal': 'CCCCCCCCC=O',
        'decanal': 'CCCCCCCCCC=O',
        'undecanal': 'CCCCCCCCCCC=O',
        'dodecanal': 'CCCCCCCCCCCC=O',
    }

    print("Step 2: Initializing descriptor calculator from the 'mordred' library.")
    # Define all necessary descriptors
    gats_descs = [getattr(descriptors.Autocorrelation, f'GATS{i}e') for i in range(1, 9)]
    chi_p_descs = [getattr(descriptors.Chi, f'Chi{i}p') for i in range(11)]
    chi_pv_descs = [getattr(descriptors.Chi, f'Chi{i}pv') for i in range(11)]
    
    # Create a single calculator for efficiency
    calc = Calculator(gats_descs + chi_p_descs + chi_pv_descs, ignore_3D=True)

    found_homologs = []
    print("\nStep 3: Processing each aldehyde to find those with 2 < GATS_max < 3.")
    # Iterate through molecules, calculate properties, and filter
    for name, smiles in aldehydes.items():
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)  # Add hydrogens for correct calculations

        # Mordred returns a Series; fill any calculation errors/missing values
        results = calc(mol).fill_missing(0.0)

        # Extract GATS values and find the maximum
        gats_values = {i: float(results[f'GATS{i}e']) for i in range(1, 9)}
        gats_max = max(gats_values.values())
        i_max = max(gats_values, key=gats_values.get)

        # Filter by GATS_max value
        if 2 < gats_max < 3:
            print(f"- Found a matching homolog: {name.capitalize()} (GATS_max = {gats_max:.4f})")
            
            # Step 4: Calculate average chi indices
            simple_chi_values = [results[f'Chi{i}p'] for i in range(11)]
            valence_chi_values = [results[f'Chi{i}pv'] for i in range(11)]

            avg_simple_chi = sum(simple_chi_values) / len(simple_chi_values)
            avg_valence_chi = sum(valence_chi_values) / len(valence_chi_values)

            # Step 5: Compute the final product
            chi_diff = avg_valence_chi - avg_simple_chi
            product = i_max * chi_diff

            found_homologs.append({
                'name': name,
                'gats_max': gats_max,
                'i_max': i_max,
                'avg_simple_chi': avg_simple_chi,
                'avg_valence_chi': avg_valence_chi,
                'chi_diff': chi_diff,
                'product': product,
            })

    # Step 6: Determine the minimum product
    if not found_homologs:
        print("\nNo homologs were found that satisfy the GATS criteria (2 < GATS_max < 3).")
        print("Cannot determine the minimum product.")
        print("\n<<<No matching homologs found>>>")
        return

    print("\nStep 6: Identifying the minimum product among the found homologs.")
    min_homolog = min(found_homologs, key=lambda x: x['product'])

    print("\n--- Calculation Details ---")
    print(f"The following homologs satisfied the GATS criteria:")
    for homolog in found_homologs:
        print(f"  - {homolog['name'].capitalize()}: Product = {homolog['product']:.4f}")

    print(f"\nThe minimum product is for {min_homolog['name'].capitalize()}.")
    print("\nThe final calculation is based on this molecule:")
    print(f"  - Lag of maximum GATS (i_max): {min_homolog['i_max']}")
    print(f"  - Average Valence Path Chi Index: {min_homolog['avg_valence_chi']:.4f}")
    print(f"  - Average Simple Path Chi Index: {min_homolog['avg_simple_chi']:.4f}")
    
    print("\nMinimum product = i_max * (Average Valence Path Chi Index - Average Simple Path Chi Index)")
    print(f"Minimum product = {min_homolog['i_max']} * ({min_homolog['avg_valence_chi']:.4f} - {min_homolog['avg_simple_chi']:.4f})")
    print(f"Minimum product = {min_homolog['i_max']} * {min_homolog['chi_diff']:.4f}")
    print(f"Minimum product = {min_homolog['product']:.4f}")
    
    # Final answer in the specified format
    print(f"\n<<<{min_homolog['product']:.4f}>>>")

calculate_properties()
