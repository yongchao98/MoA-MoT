def solve_molecular_formula():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    """
    # Step 1: Define the molecular formula of the starting material.
    # Base structure (2-azabicyclo[2.2.1]hept-5-en-3-one): C6H7NO
    # Add CF3 group (replaces one H): +C, -H, +3F -> C7H6F3NO
    # Add PMB group (C8H9O) (replaces one H on N): +C8, +H9, +O, -H -> C15H14F3NO2
    atoms = {'C': 15, 'H': 14, 'F': 3, 'N': 1, 'O': 2}
    print("Starting Material Formula: C15H14F3NO2")
    
    # Step 2: Reaction 1 -> Intermediate 1 (PMB deprotection)
    # Remove PMB group (C8H9O) and add one H.
    # Net change: -C8, -H8, -O1
    atoms['C'] -= 8
    atoms['H'] -= 8
    atoms['O'] -= 1
    # Intermediate 1: C7H6F3NO
    
    # Step 3: Reaction 2 -> Intermediate 2 (Hydrogenation)
    # Add two H atoms to reduce the C=C bond.
    atoms['H'] += 2
    # Intermediate 2: C7H8F3NO

    # Step 4: Reaction 3 -> Product (Lactam Hydrolysis)
    # Add one molecule of water (H2O).
    atoms['H'] += 2
    atoms['O'] += 1
    # Product: C7H10F3NO2
    
    print("\nThe final product's molecular formula is calculated as follows:")
    # Using Hill system order (C, H, then alphabetical)
    print(f"Number of Carbon atoms: {atoms['C']}")
    print(f"Number of Hydrogen atoms: {atoms['H']}")
    print(f"Number of Fluorine atoms: {atoms['F']}")
    print(f"Number of Nitrogen atoms: {atoms['N']}")
    print(f"Number of Oxygen atoms: {atoms['O']}")

    final_formula = f"C{atoms['C']}H{atoms['H']}F{atoms['F']}N{atoms['N']}O{atoms['O']}"
    print(f"\nFinal Molecular Formula: {final_formula}")

solve_molecular_formula()