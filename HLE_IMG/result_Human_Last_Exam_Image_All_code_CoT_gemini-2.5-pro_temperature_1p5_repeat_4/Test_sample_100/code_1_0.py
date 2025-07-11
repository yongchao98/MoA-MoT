def calculate_product_formula():
    """
    This function calculates the molecular formula of the product
    by tracking the changes in atomic composition through the given reaction scheme.
    """
    # Step 1: Define the composition of the starting material.
    # Starting Material Formula: C16H16F3NO2
    atoms = {'C': 16, 'H': 16, 'F': 3, 'N': 1, 'O': 2}

    # Step 2: Reaction 1 - PMB (C8H9O) deprotection with CAN.
    # A hydrogen atom replaces the PMB group on the nitrogen.
    # Net change: Lose C8H8O.
    atoms['C'] -= 8
    atoms['H'] -= 9
    atoms['O'] -= 1
    atoms['H'] += 1 # Add H to Nitrogen

    # Step 3: Reaction 2 - Hydrogenation of the C=C double bond.
    # Adds two hydrogen atoms.
    atoms['H'] += 2

    # Step 4: Reaction 3 - Hydrolysis of the lactam (cyclic amide).
    # Adds one molecule of water (H2O).
    atoms['H'] += 2
    atoms['O'] += 1

    # Print the final counts for each element
    print("The final product is composed of:")
    print(f"Carbon (C) atoms: {atoms['C']}")
    print(f"Hydrogen (H) atoms: {atoms['H']}")
    print(f"Fluorine (F) atoms: {atoms['F']}")
    print(f"Nitrogen (N) atoms: {atoms['N']}")
    print(f"Oxygen (O) atoms: {atoms['O']}")

    # Format the final molecular formula string
    final_formula = f"C{atoms['C']}H{atoms['H']}F{atoms['F']}N{atoms['N']}O{atoms['O']}"
    print(f"\nThe molecular formula of the product is: {final_formula}")

calculate_product_formula()