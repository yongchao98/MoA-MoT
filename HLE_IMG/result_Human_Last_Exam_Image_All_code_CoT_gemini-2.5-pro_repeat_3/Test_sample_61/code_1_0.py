def calculate_product_formula():
    """
    Calculates the molecular formula of product A based on the reaction described.
    """
    # Step 1: Define the molecular formulas of the reactants
    # Compound 1: methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate
    # Indanone (C9H8O) + COOMe (C2H3O2) - H = C11H10O3
    compound_1 = {'C': 11, 'H': 10, 'O': 3}
    
    # Compound 2: benzyl bromide (C6H5CH2Br)
    compound_2 = {'C': 7, 'H': 7, 'Br': 1}

    print("--- Reaction Analysis ---")
    print(f"Formula of Compound 1: C{compound_1['C']}H{compound_1['H']}O{compound_1['O']}")
    print(f"Formula of Compound 2: C{compound_2['C']}H{compound_2['H']}Br")
    
    # Step 2: Calculate the formula of the alkylated intermediate
    # This involves replacing an H in compound 1 with a benzyl group (C7H7) from compound 2.
    # Intermediate = (Compound 1 - H) + (Compound 2 - Br)
    intermediate = {}
    intermediate['C'] = compound_1['C'] + compound_2['C']
    intermediate['H'] = compound_1['H'] - 1 + compound_2['H']
    intermediate['O'] = compound_1['O']
    print(f"\nAlkylation step adds a benzyl group (C7H7) for an acidic H.")
    print(f"Formula of Intermediate: C{intermediate['C']}H{intermediate['H']}O{intermediate['O']}")

    # Step 3: Calculate the formula of the final product A
    # Saponification and decarboxylation replaces the -COOMe group with an H atom.
    # Atoms in -COOMe group:
    coome_group = {'C': 2, 'H': 3, 'O': 2}
    
    product_A = {}
    product_A['C'] = intermediate['C'] - coome_group['C']
    product_A['H'] = intermediate['H'] - coome_group['H'] + 1
    product_A['O'] = intermediate['O'] - coome_group['O']

    print(f"\nSaponification and decarboxylation replaces -COOMe (C{coome_group['C']}H{coome_group['H']}O{coome_group['O']}) with H.")
    print("\n--- Final Formula Calculation for Product A ---")
    print(f"Carbon atoms: {intermediate['C']} - {coome_group['C']} = {product_A['C']}")
    print(f"Hydrogen atoms: {intermediate['H']} - {coome_group['H']} + 1 = {product_A['H']}")
    print(f"Oxygen atoms: {intermediate['O']} - {coome_group['O']} = {product_A['O']}")
    
    print(f"\nThe molecular formula of compound A is C{product_A['C']}H{product_A['H']}O{product_A['O']}.")

calculate_product_formula()