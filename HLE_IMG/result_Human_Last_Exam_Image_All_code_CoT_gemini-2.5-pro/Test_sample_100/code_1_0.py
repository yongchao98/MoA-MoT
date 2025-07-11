def calculate_product_formula():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    """
    # Step 1: Determine the molecular formula of the starting material.
    # The starting material is 2-(p-methoxybenzyl)-5-(trifluoromethyl)-2-azabicyclo[2.2.1]hept-5-en-3-one.
    # Core (C6H5NO) + CF3 (CF3) + PMB (C8H9O)
    # C: 6 + 1 + 8 = 15
    # H: 5 (from core) + 9 (from PMB) = 14
    # F: 3
    # N: 1
    # O: 1 (from core) + 1 (from PMB) = 2
    # Formula: C15H14F3NO2
    formula = {'C': 15, 'H': 14, 'F': 3, 'N': 1, 'O': 2}
    
    # Step 2: First reaction - Deprotection of PMB group with CAN.
    # The PMB group (C8H9O) is removed, and a hydrogen atom is added to the nitrogen.
    # Change: -C8, -H9, -O1, +H1
    formula['C'] -= 8
    formula['H'] -= 9
    formula['O'] -= 1
    formula['H'] += 1
    # Formula of Intermediate 1: C7H6F3NO

    # Step 3: Second reaction - Hydrogenation with Pd/C, H2.
    # The C=C double bond is reduced, which adds two hydrogen atoms (H2).
    formula['H'] += 2
    # Formula of Intermediate 2: C7H8F3NO

    # Step 4: Third reaction - Lactam hydrolysis with 4 N HCl.
    # The lactam ring is opened by adding one molecule of water (H2O).
    formula['H'] += 2
    formula['O'] += 1
    # Final Product Formula: C7H10F3NO2

    # Step 5: Format and print the final molecular formula.
    # The convention is C, H, then other elements in alphabetical order.
    elements_order = ['C', 'H', 'F', 'N', 'O']
    final_formula_str = ""
    print("The final product has the following elemental composition:")
    
    for elem in elements_order:
        count = formula.get(elem, 0)
        if count > 0:
            print(f"{elem}: {count}")
            final_formula_str += elem
            # As requested, explicitly output each number in the formula.
            final_formula_str += str(count)

    print("\nThe molecular formula of the product is:")
    print(final_formula_str)

calculate_product_formula()