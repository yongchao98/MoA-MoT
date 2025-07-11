import sys
from rdkit import Chem
from mordred import Calculator, descriptors

def solve_molecular_properties():
    """
    Finds formaldehyde homologs based on specific molecular properties and calculates
    the minimum product of i_max and a chi index difference.
    """
    # Create a calculator with the specific descriptors we need.
    # This is much faster than calculating all descriptors.
    calc = Calculator([
        descriptors.Geary(weighting='sanderson'),
        descriptors.Chi.AVPC(),
        descriptors.Chi.ASPC()
    ])

    min_product = float('inf')
    best_homolog = None

    # Name mapping for aldehydes
    aldehyde_names = {
        1: "Methanal (Formaldehyde)", 2: "Ethanal", 3: "Propanal", 4: "Butanal",
        5: "Pentanal", 6: "Hexanal", 7: "Heptanal", 8: "Octanal", 9: "Nonanal",
        10: "Decanal", 11: "Undecanal", 12: "Dodecanal", 13: "Tridecanal",
        14: "Tetradecanal", 15: "Pentadecanal", 16: "Hexadecanal", 17: "Heptadecanal",
        18: "Octadecanal", 19: "Nonadecanal", 20: "Icosanal"
    }

    print("Analyzing formaldehyde homologs (aldehydes C1-C30)...\n")

    # Iterate through aldehydes with 1 to 30 carbon atoms
    for n in range(1, 31):
        # Generate SMILES string for the aldehyde homolog
        # 'C=O' for Methanal (n=1), 'CC=O' for Ethanal (n=2), etc.
        smiles = 'C' * n + '=O'
        
        # Create an RDKit molecule object
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            continue
        
        # Add explicit hydrogens, important for descriptor calculation
        mol = Chem.AddHs(mol)

        # Calculate the descriptors
        try:
            results = calc(mol)
        except Exception:
            continue
        
        # Extract Geary Autocorrelation values (GATS1e to GATS8e)
        geary_values = {}
        for i in range(1, 9):
            key = f'GATS{i}e'
            if key in results:
                # Mordred can return missing values as special objects, filter them out
                if hasattr(results[key], 'real'):
                     geary_values[i] = results[key].real
        
        if not geary_values:
            continue

        # Find the maximum Geary value and the lag (i_max) at which it occurs
        i_max = max(geary_values, key=geary_values.get)
        max_geary = geary_values[i_max]

        # Check if the homolog meets the condition
        if 2 <= max_geary <= 3:
            homolog_name = aldehyde_names.get(n, f"{n}-Carbon Aldehyde")
            print(f"Found a qualifying homolog: {homolog_name} ({n} carbons)")
            print(f" - Max Geary Autocorrelation (GATS): {max_geary:.4f} (at lag i_max = {i_max})")

            # Get the average chi indices
            avpc = results['AVPC'].real
            aspc = results['ASPC'].real
            delta_chi = avpc - aspc
            
            print(f" - Average Valence Path Chi (AVPC): {avpc:.4f}")
            print(f" - Average Simple Path Chi (ASPC): {aspc:.4f}")
            print(f" - Difference (AVPC - ASPC): {avpc:.4f} - {aspc:.4f} = {delta_chi:.4f}")

            # Calculate the product of i_max and the chi difference
            product = i_max * delta_chi
            print(f" - Product (i_max * Difference): {i_max} * {delta_chi:.4f} = {product:.4f}")
            print("-" * 20)

            # Check if this product is the new minimum
            if product < min_product:
                min_product = product
                best_homolog = (homolog_name, product)

    if best_homolog:
        name, min_val = best_homolog
        print(f"\nCalculation complete.")
        print(f"The minimum product of {min_val:.4f} was found for {name}.")
        sys.stdout.write(f"<<<{min_val:.4f}>>>\n")
    else:
        print("\nNo homologs found that satisfy the specified Geary autocorrelation criteria.")
        sys.stdout.write("<<<No answer found>>>\n")

if __name__ == "__main__":
    solve_molecular_properties()