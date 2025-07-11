def calculate_product_formula():
    """
    Calculates the molecular formula of the final product through a three-step reaction.
    """
    # Step 1: Define the molecular formula of the starting material.
    # The starting material is 2-(4-methoxybenzyl)-5-(trifluoromethyl)-2-azabicyclo[2.2.1]hept-5-en-3-one.
    # Core (C6H5NO) + PMB group (C8H9O) + CF3 group (CF3)
    # C: 6 + 8 + 1 = 15
    # H: 5 + 9 = 14
    # F: 3
    # N: 1
    # O: 1 + 1 = 2
    # Formula: C15H14F3NO2
    atoms = {
        'C': 15,
        'H': 14,
        'F': 3,
        'N': 1,
        'O': 2
    }
    print(f"Starting material formula: C{atoms['C']}H{atoms['H']}F{atoms['F']}NO{atoms['O']}")

    # Step 2: PMB group deprotection with CAN.
    # The PMB group (C8H9O) is replaced by a hydrogen (H) atom.
    # The net change is the removal of a C8H8O fragment.
    atoms['C'] -= 8
    atoms['H'] -= 8
    atoms['O'] -= 1
    print(f"Intermediate 1 formula (after PMB deprotection): C{atoms['C']}H{atoms['H']}F{atoms['F']}NO{atoms['O']}")

    # Step 3: Hydrogenation of the C=C double bond with Pd/C and H2.
    # This adds two hydrogen atoms.
    atoms['H'] += 2
    print(f"Intermediate 2 formula (after hydrogenation): C{atoms['C']}H{atoms['H']}F{atoms['F']}NO{atoms['O']}")

    # Step 4: Hydrolysis of the lactam with acid.
    # This reaction adds one molecule of water (H2O).
    atoms['H'] += 2
    atoms['O'] += 1
    print(f"Final product formula (after hydrolysis): C{atoms['C']}H{atoms['H']}F{atoms['F']}NO{atoms['O']}")

    # Format the final result into a standard molecular formula string.
    # The convention is to omit the number '1'.
    final_formula = ""
    # Standard order of elements: C, H, then alphabetical.
    elements_in_order = ['C', 'H', 'F', 'N', 'O']
    for element in elements_in_order:
        count = atoms.get(element, 0)
        if count > 0:
            final_formula += element
            if count > 1:
                final_formula += str(count)

    print("\nThe molecular formula of the final product is:")
    # The requested output format is the full equation with each number.
    print(f"C{atoms['C']}H{atoms['H']}F{atoms['F']}N{atoms['N']}O{atoms['O']}")


calculate_product_formula()