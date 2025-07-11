import sys
import subprocess

# Install necessary libraries if they are not already installed
try:
    from rdkit import Chem
    from mordred import Calculator, descriptors
    from mordred.AtomProperty import ene_s
except ImportError:
    print("Installing required libraries: rdkit-pypi and mordred...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "rdkit-pypi", "mordred"])
    from rdkit import Chem
    from mordred import Calculator, descriptors
    from mordred.AtomProperty import ene_s

def solve_homologs_problem():
    """
    Finds formaldehyde homologs based on a Geary autocorrelation criterion and
    determines the minimum product of i_max and a path chi index difference.
    """
    # Define the homologs of formaldehyde to be analyzed
    aldehydes = {
        "Formaldehyde": "C=O",
        "Acetaldehyde": "CC=O",
        "Propionaldehyde": "CCC=O",
        "Butyraldehyde": "CCCC=O",
        "Pentanal": "CCCCC=O",
        "Hexanal": "CCCCCC=O",
        "Heptanal": "CCCCCCC=O",
        "Octanal": "CCCCCCCC=O",
        "Nonanal": "CCCCCCCCC=O",
        "Decanal": "CCCCCCCCCC=O"
    }

    # Initialize the mordred calculator
    calc = Calculator(descriptors, ignore_3D=True)

    # We need Geary autocorrelation (GATS) with Sanderson electronegativity ('e').
    # Mordred calculates this for lags 1-8 by default.
    calc.register(descriptors.GearyAuto(prop=ene_s))

    # We need the average path chi indices, simple (mASP) and valence (mAVP).
    calc.register(descriptors.Path.AvPath(valence=True))
    calc.register(descriptors.Path.AvPath(valence=False))

    found_homologs_products = []
    print("Analyzing formaldehyde homologs...")
    print("-" * 60)

    # Iterate through each aldehyde
    for name, smiles in aldehydes.items():
        # Create and prepare the molecule object
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)

        # Calculate all registered descriptors
        results = calc(mol)

        # Extract Geary autocorrelation values and find the maximum
        gats_values = [results[f"GATS{i}e"] for i in range(1, 9)]
        
        # Handle cases where calculation fails for a molecule
        if any(isinstance(v, Exception) for v in gats_values):
            print(f"Skipping {name} due to a calculation error.\n" + "-"*60)
            continue
            
        max_gats = max(gats_values)
        i_max = gats_values.index(max_gats) + 1  # Lag is 1-based index

        print(f"Molecule: {name} ({smiles})")
        print(f"  Maximum Geary Autocorrelation is GATS{i_max}e = {max_gats:.4f}")

        # Check if the homolog meets the specified condition
        if 2.0 <= max_gats <= 3.0:
            print(f"  CONDITION MET: This homolog is selected.")

            # Get the path chi indices
            mAVP = results['mAVP']
            mASP = results['mASP']
            
            # Calculate the difference and the final product
            delta_chi = mAVP - mASP
            product = i_max * delta_chi
            found_homologs_products.append(product)

            # Print the components of the final equation as requested
            print(f"\n  Calculating the product for {name}:")
            print(f"  i_max = {i_max}")
            print(f"  Average Valence Path Chi Index (mAVP) = {mAVP:.4f}")
            print(f"  Average Simple Path Chi Index (mASP) = {mASP:.4f}")
            print(f"  Difference (delta_chi) = {mAVP:.4f} - {mASP:.4f} = {delta_chi:.4f}")
            print(f"  Product = i_max * delta_chi = {i_max} * {delta_chi:.4f} = {product:.4f}")
        else:
            print("  CONDITION NOT MET: This homolog is not selected.")
        
        print("-" * 60)

    # Determine the minimum product from all found homologs
    if not found_homologs_products:
        print("\nNo homologs were found that satisfy the specified condition.")
        final_answer = "N/A"
    else:
        min_product = min(found_homologs_products)
        print(f"\nSummary: Found {len(found_homologs_products)} homologs meeting the criteria.")
        print(f"Calculated products: {[f'{p:.4f}' for p in found_homologs_products]}")
        print(f"\nThe minimum product among the found homologs is: {min_product:.4f}")
        final_answer = f"{min_product:.4f}"
    
    return final_answer

if __name__ == '__main__':
    result = solve_homologs_problem()
    print(f"<<<{result}>>>")
