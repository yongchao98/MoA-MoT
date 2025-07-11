def calculate_molecular_formula():
    """
    Calculates the molecular formula of product A based on the described reaction.
    """
    # Molecular formula of Compound 1 (Methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate)
    compound_1 = {'C': 11, 'H': 10, 'O': 3}

    # Step 1: Alkylation
    # A hydrogen atom from compound 1 is replaced by a benzyl group (C7H7).
    # We add the atoms of the benzyl group and subtract one hydrogen atom.
    intermediate = {
        'C': compound_1['C'] + 7,
        'H': compound_1['H'] - 1 + 7,
        'O': compound_1['O']
    }
    # The formula of the alkylated intermediate is C18H16O3.

    # Step 2: Saponification and Decarboxylation
    # The ester group (-COOCH3) is replaced by a hydrogen atom.
    # Atoms in -COOCH3 group: C=2, H=3, O=2.
    # Net change: Subtract C2H3O2 and add H, which is equivalent to subtracting C2H2O2.
    product_A = {
        'C': intermediate['C'] - 2,
        'H': intermediate['H'] - 2,
        'O': intermediate['O'] - 2
    }

    # Extract the final atom counts
    C_count = product_A['C']
    H_count = product_A['H']
    O_count = product_A['O']
    
    # Print the final molecular formula, including the numbers for each element.
    # The final product is 2-benzyl-2,3-dihydro-1H-inden-1-one.
    print(f"C{C_count}H{H_count}O{'' if O_count == 1 else O_count}")

calculate_molecular_formula()