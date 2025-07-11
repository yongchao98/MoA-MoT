def calculate_molecular_formula():
    """
    Calculates the molecular formula of product A based on the reaction stoichiometry.

    The reaction is the alkylation of methyl 1-oxoindane-2-carboxylate with benzyl bromide,
    followed by saponification and decarboxylation.
    The net transformation is the replacement of the -COOMe group with a benzyl group.
    """
    # Molecular formula of Compound 1 (Methyl 1-oxoindane-2-carboxylate)
    # C: 9 in indanone skeleton + 1 in ester C=O + 1 in O-Me = 11
    # H: 8 in indanone skeleton + 2 on C3 and C2? Let's count properly:
    # C6H4 (benzo) + C3H3O (cyclopentanone part) + COOCH3
    # C: 6+3+1+1 = 11
    # H: 4 (benzo) + 1 (C2) + 2 (C3) + 3 (Me) = 10
    # O: 1 (ketone) + 2 (ester) = 3
    compound_1 = {'C': 11, 'H': 10, 'O': 3}

    # Atoms in the group being removed: -COOMe (methyl carboxylate)
    # C: 1 (carbonyl) + 1 (methyl) = 2
    # H: 3 (methyl)
    # O: 2
    group_removed = {'C': 2, 'H': 3, 'O': 2}

    # Atoms in the group being added: -C7H7 (benzyl group)
    group_added = {'C': 7, 'H': 7, 'O': 0}

    # Calculate the molecular formula of product A
    product_A = {}
    product_A['C'] = compound_1['C'] - group_removed['C'] + group_added['C']
    product_A['H'] = compound_1['H'] - group_removed['H'] + group_added['H']
    product_A['O'] = compound_1['O'] - group_removed['O'] + group_added['O']

    # Print the explanation of the calculation
    print("Calculation Steps:")
    print(f"1. Start with Compound 1 formula: C{compound_1['C']}H{compound_1['H']}O{compound_1['O']}")
    print(f"2. The net reaction replaces the -COOMe group (C{group_removed['C']}H{group_removed['H']}O{group_removed['O']}) with a benzyl group (C{group_added['C']}H{group_added['H']}).")
    print("\nFinal Atom Count Calculation:")
    print(f"Carbon (C): {compound_1['C']} - {group_removed['C']} + {group_added['C']} = {product_A['C']}")
    print(f"Hydrogen (H): {compound_1['H']} - {group_removed['H']} + {group_added['H']} = {product_A['H']}")
    print(f"Oxygen (O): {compound_1['O']} - {group_removed['O']} + {group_added['O']} = {product_A['O']}")

    # Print the final molecular formula
    print("\nThe molecular formula of compound A is:")
    print(f"C{product_A['C']}H{product_A['H']}O{product_A['O']}")

calculate_molecular_formula()