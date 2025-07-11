def calculate_product_formula():
    """
    This function calculates the molecular formula of the final product in the given reaction scheme.
    """
    # Step 1: Define the molecular formula of the starting material.
    # The starting material is 2-(p-methoxybenzyl)-5-(trifluoromethyl)-2-azabicyclo[2.2.1]hept-5-en-3-one.
    # Its formula is C15H14F3NO2.
    atoms = {
        'C': 15,
        'H': 14,
        'F': 3,
        'N': 1,
        'O': 2
    }

    # Step 2: Account for Reaction 1 (CAN deprotection).
    # This reaction removes the PMB group (C8H9O) and adds an H atom to the nitrogen.
    # Net change: remove C8H8O.
    atoms['C'] -= 8
    atoms['H'] -= 8
    atoms['O'] -= 1
    # Intermediate 1 formula: C7H6F3NO

    # Step 3: Account for Reaction 2 (Hydrogenation).
    # Pd/C, H2 reduces the C=C double bond, adding two hydrogen atoms (H2).
    atoms['H'] += 2
    # Intermediate 2 formula: C7H8F3NO

    # Step 4: Account for Reaction 3 (Lactam Hydrolysis).
    # Acid hydrolysis with 4 N HCl cleaves the lactam by adding one molecule of water (H2O).
    atoms['H'] += 2
    atoms['O'] += 1
    # Final product formula: C7H10F3NO2

    # Step 5: Print the final atom counts and the molecular formula.
    # Following the prompt to "output each number in the final equation".
    c = atoms['C']
    h = atoms['H']
    f = atoms['F']
    n = atoms['N']
    o = atoms['O']

    print("Final atom counts for the product:")
    print(f"Carbon (C): {c}")
    print(f"Hydrogen (H): {h}")
    print(f"Fluorine (F): {f}")
    print(f"Nitrogen (N): {n}")
    print(f"Oxygen (O): {o}")

    # Construct the molecular formula string according to chemical conventions (Hill system).
    # C, then H, then other elements alphabetically.
    # For a single atom of an element (like Nitrogen here), the number '1' is typically omitted.
    formula = f"C{c}H{h}F{f}NO{o}"
    print(f"\nThe molecular formula of the product is: {formula}")

calculate_product_formula()