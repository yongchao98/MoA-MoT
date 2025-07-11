import sys
try:
    from rdkit import Chem
    from mordred import Calculator, descriptors
    from mordred.Geary import Geary
    from mordred.Chi import AVPChi, APChi
except ImportError:
    print("Please install rdkit-pypi and mordred first, e.g., 'pip install rdkit-pypi mordred'", file=sys.stderr)
    sys.exit(1)

def solve_chemoinformatics_task():
    """
    Solves the user's request by finding the minimum product of i_max and a chi index difference
    for specific formaldehyde homologs.
    """
    print("Step 1: Analyzing homologs of formaldehyde to find candidates.")
    print("The condition is that the maximum Geary Autocorrelation (weighted by Sanderson electronegativity) must be between 2.0 and 3.0.")
    print("-" * 70)

    # Generate homologs of formaldehyde (aldehydes)
    homolog_names = [
        "Formaldehyde", "Acetaldehyde", "Propanal", "Butanal", "Pentanal",
        "Hexanal", "Heptanal", "Octanal", "Nonanal", "Decanal",
        "Undecanal", "Dodecanal", "Tridecanal", "Tetradecanal", "Pentadecanal"
    ]
    # SMILES for C1 to C15 aldehydes
    smiles_list = ['C' * n + '=O' for n in range(1, len(homolog_names) + 1)]
    mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in smiles_list]

    found_homologs = []
    
    # Iterate through homologs to find i_max and check the Geary value
    for name, mol in zip(homolog_names, mols):
        if not mol:
            continue

        max_geary_val = -float('inf')
        i_max = -1
        # Test lags from 1 up to 10 (a common upper limit for path-based descriptors)
        # The maximum possible lag is mol.GetNumAtoms() - 1
        max_lag_to_test = min(mol.GetNumAtoms() - 1, 10)

        for lag in range(1, max_lag_to_test + 1):
            try:
                # mordred's Geary requires instantiation per lag
                calc = Calculator([Geary(lag=lag, prop='e')])
                # Key format is 'GATS' + lag + property, e.g., 'GATS4e'
                geary_val = calc(mol).as_dict()[f'GATS{lag}e']

                if geary_val > max_geary_val:
                    max_geary_val = geary_val
                    i_max = lag
            except Exception:
                # Calculation may fail for some lags on small molecules
                continue
        
        print(f"Molecule: {name:<12} | Max Geary Value = {max_geary_val:.4f} (at lag i_max = {i_max})")

        # Filter for homologs where the max Geary value is between 2 and 3
        if 2.0 <= max_geary_val <= 3.0:
            found_homologs.append({
                'name': name,
                'mol': mol,
                'i_max': i_max
            })
            print(f"  -> Found a candidate: {name}")

    print("-" * 70)
    print("Step 2: Calculating indices and the final product for the candidate homologs.")

    if not found_homologs:
        print("No homologs were found that satisfy the Geary autocorrelation condition.")
        return

    min_product = float('inf')
    final_details = {}

    chi_calc = Calculator([AVPChi(), APChi()])
    
    for homolog in found_homologs:
        name = homolog['name']
        mol = homolog['mol']
        i_max = homolog['i_max']
        
        print(f"\nProcessing candidate: {name}")
        
        # Calculate chi indices
        chi_results = chi_calc(mol).as_dict()
        avp_chi = chi_results['AVPChi']
        ap_chi = chi_results['APChi']
        
        # Calculate the difference and the final product
        chi_diff = avp_chi - ap_chi
        product = i_max * chi_diff

        print(f"  i_max = {i_max}")
        print(f"  Average Valence Path Chi (AVPChi) = {avp_chi:.4f}")
        print(f"  Average Simple Path Chi (APChi) = {ap_chi:.4f}")
        print(f"  Difference (AVPChi - APChi) = {chi_diff:.4f}")
        print(f"  Product (i_max * Difference) = {product:.4f}")
        
        if product < min_product:
            min_product = product
            final_details = {
                'name': name,
                'product': product,
                'i_max': i_max,
                'avp_chi': avp_chi,
                'ap_chi': ap_chi
            }
            
    print("-" * 70)
    print("Step 3: Determining the minimum product.")

    if 'name' in final_details:
        result = final_details
        print(f"The minimum product is for the homolog '{result['name']}'.")
        print(f"The calculation is:")
        # Output each number in the final equation as requested
        print(f"{result['i_max']} * ({result['avp_chi']:.4f} - {result['ap_chi']:.4f}) = {result['product']:.4f}")
        
        # Final answer in the specified format
        print(f"\nFinal Answer Rounded to 4 Decimal Places: {result['product']:.4f}")
        print(f"\n<<<{result['product']:.4f}>>>")
    else:
        # This case should not be reached if candidates were found
        print("Could not determine the minimum product.")

if __name__ == '__main__':
    solve_chemoinformatics_task()