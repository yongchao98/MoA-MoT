def calculate_molecular_formula():
    """
    Calculates the molecular formula of compound A based on the reaction steps.
    """
    # Step 1: Define the atomic composition of the starting material, Reactant 1.
    # Reactant 1 is methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate: C11H10O3
    # C: 6(benzo) + 1(C=O) + 1(CH) + 1(CH2) + 1(ester C=O) + 1(methyl) = 11
    # H: 4(benzo) + 1(CH) + 2(CH2) + 3(methyl) = 10
    # O: 1(ketone) + 2(ester) = 3
    reactant1 = {'C': 11, 'H': 10, 'O': 3}

    # The reaction involves alkylation with a benzyl group followed by decarboxylation.

    # Transformation 1: Alkylation
    # The acidic alpha-hydrogen is replaced by a benzyl group (C7H7).
    # Net change: -H + C7H7
    benzyl_group = {'C': 7, 'H': 7}
    
    # Calculate the formula of the alkylated intermediate
    intermediate = reactant1.copy()
    intermediate['C'] += benzyl_group['C']
    intermediate['H'] += benzyl_group['H'] - 1  # remove one H, add the benzyl group

    # Transformation 2: Saponification and Decarboxylation
    # The methoxycarbonyl group (-COOMe) is replaced by a hydrogen atom.
    # Atoms in -COOMe: C=2, H=3, O=2
    coome_group = {'C': 2, 'H': 3, 'O': 2}
    
    # Calculate the formula of the final product A
    product_A = intermediate.copy()
    product_A['C'] -= coome_group['C']
    product_A['H'] -= coome_group['H']
    product_A['O'] -= coome_group['O']
    product_A['H'] += 1 # Add one H

    c_count = product_A['C']
    h_count = product_A['H']
    o_count = product_A['O']

    print("The molecular formula of compound A is CxHyOz.")
    print(f"The number of Carbon atoms (x) is: {c_count}")
    print(f"The number of Hydrogen atoms (y) is: {h_count}")
    print(f"The number of Oxygen atoms (z) is: {o_count}")
    print(f"\nFinal Molecular Formula: C{c_count}H{h_count}O{o_count if o_count > 1 else ''}")

calculate_molecular_formula()