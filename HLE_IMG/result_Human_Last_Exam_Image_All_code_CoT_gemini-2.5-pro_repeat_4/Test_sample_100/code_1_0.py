def calculate_product_formula():
    """
    Calculates the molecular formula of the product from the given reaction sequence.
    """
    # Plan Explanation
    print("This script determines the molecular formula of the final product by tracking atomic changes through the three reaction steps.")
    print("-----------------------------------------------------------------------------------------------------------------")

    # Step 0: Determine the formula of the starting material.
    # Structure: 5-(trifluoromethyl)-2-(4-methoxybenzyl)-2-azabicyclo[2.2.1]hept-5-en-3-one
    # Core (C6H5): 6 C, 5 H
    # Groups: C=O (1 O), N, CF3 (1 C, 3 F), PMB (p-methoxybenzyl, C8H9O)
    # Total C = 6 (core) + 1 (CF3) + 8 (PMB) = 15
    # Total H = 5 (core) + 9 (PMB) = 14
    # Total F = 3
    # Total N = 1
    # Total O = 1 (carbonyl) + 1 (PMB) = 2
    atoms = {'C': 15, 'H': 14, 'F': 3, 'N': 1, 'O': 2}
    print("Step 0: The starting material is C15H14F3NO2.")
    print(f"Initial atom counts: {atoms}")
    print("-----------------------------------------------------------------------------------------------------------------")


    # Step 1: Deprotection of the PMB group with CAN.
    # The PMB group (p-methoxybenzyl, C8H9O) is removed from the Nitrogen.
    # A Hydrogen atom is added to the Nitrogen.
    # Net change: -C8, -H9, -O1, +H1 => -C8, -H8, -O1
    atoms['C'] -= 8
    atoms['H'] -= 9
    atoms['O'] -= 1
    atoms['H'] += 1
    print("Step 1: PMB group is removed (-C8H9O) and an H atom is added.")
    print(f"Atom counts after Step 1: {atoms}")
    print("-----------------------------------------------------------------------------------------------------------------")


    # Step 2: Catalytic hydrogenation with Pd/C and H2.
    # The C=C double bond in the ring is reduced to a single bond.
    # Two Hydrogen atoms are added.
    # Net change: +H2
    atoms['H'] += 2
    print("Step 2: C=C double bond is hydrogenated (+H2).")
    print(f"Atom counts after Step 2: {atoms}")
    print("-----------------------------------------------------------------------------------------------------------------")


    # Step 3: Acidic hydrolysis of the lactam.
    # The lactam (amide) ring is opened by adding one molecule of water (H2O).
    # Net change: +H2O
    atoms['H'] += 2
    atoms['O'] += 1
    print("Step 3: Lactam ring is hydrolyzed by adding water (+H2O).")
    print(f"Atom counts after Step 3: {atoms}")
    print("-----------------------------------------------------------------------------------------------------------------")

    # Final product formula calculation
    c = atoms['C']
    h = atoms['H']
    f = atoms['F']
    n = atoms['N']
    o = atoms['O']

    print("The final count for each element in the product is:")
    print(f"Carbon (C): {c}")
    print(f"Hydrogen (H): {h}")
    print(f"Fluorine (F): {f}")
    print(f"Nitrogen (N): {n}")
    print(f"Oxygen (O): {o}")

    # Constructing the molecular formula string in the standard order (C, H, then alphabetical)
    # Numbers equal to 1 are omitted.
    formula_parts = []
    formula_parts.append(f"C{c}")
    formula_parts.append(f"H{h}")
    formula_parts.append(f"F{f}")
    formula_parts.append("N" if n == 1 else f"N{n}")
    formula_parts.append(f"O{o}")
    final_formula = "".join(formula_parts)

    print("\nThe molecular formula of the product is:")
    print(final_formula)


calculate_product_formula()