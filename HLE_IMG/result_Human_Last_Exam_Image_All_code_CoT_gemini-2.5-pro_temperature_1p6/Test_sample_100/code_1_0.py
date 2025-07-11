def calculate_product_formula():
    """
    Calculates the molecular formula of the final product by tracking the atomic changes
    throughout the provided three-step chemical reaction.
    """

    # Step 1: Define the molecular formula of the starting material.
    # The structure is a complex bicyclic lactam. Let's break it down to count the atoms.
    # - Bicyclic lactam core with a C=C bond and a C=O group: C6H5NO
    # - Trifluoromethyl group (CF3) replacing one H on the C=C bond: C1F3
    # - p-Methoxybenzyl (PMB) group on the Nitrogen: C8H9O
    # So, the initial counts are:
    atoms = {
        'C': 6 + 1 + 8,  # Sum of carbons from core, CF3, and PMB
        'H': 5 + 9,      # Sum of hydrogens from core and PMB
        'F': 3,          # From the CF3 group
        'N': 1,          # From the lactam core
        'O': 1 + 1,      # From the lactam C=O and PMB's methoxy group
    }

    # Step 2: The first reaction is PMB group deprotection using CAN.
    # This reaction removes the PMB group (C8H9O) and adds one Hydrogen atom to the Nitrogen.
    # Net change: -C8, -H9, -O1, +H1  =>  -C8, -H8, -O1
    atoms['C'] -= 8
    atoms['H'] -= 8
    atoms['O'] -= 1
    # Formula of Intermediate 1: C7H6F3NO

    # Step 3: The second reaction is the hydrogenation of the C=C double bond.
    # This adds two Hydrogen atoms across the double bond.
    # Net change: +H2
    atoms['H'] += 2
    # Formula of Intermediate 2: C7H8F3NO

    # Step 4: The third reaction is the acid-catalyzed hydrolysis of the lactam.
    # This reaction opens the ring by adding one molecule of water (H2O).
    # Net change: +H2, +O1
    atoms['H'] += 2
    atoms['O'] += 1
    # Formula of the final product: C7H10F3NO2

    # Step 5: Print the composition and the final molecular formula.
    print("The atomic composition of the final product is:")
    element_order = ['C', 'H', 'F', 'N', 'O']
    formula_string = ""
    for element in element_order:
        count = atoms.get(element, 0)
        if count > 0:
            print(f"{element}: {count}")
            formula_string += element
            if count > 1:
                formula_string += str(count)

    print(f"\nThe determined molecular formula of the product is: {formula_string}")


calculate_product_formula()
<<<C7H10F3NO2>>>