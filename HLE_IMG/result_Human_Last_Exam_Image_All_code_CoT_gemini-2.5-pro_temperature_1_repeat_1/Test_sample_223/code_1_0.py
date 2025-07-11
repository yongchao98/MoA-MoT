def calculate_molecular_formula():
    """
    Calculates the molecular formula of compound B based on the structure of compound A
    and the reagents used in the reaction.
    """
    # Molecular formula of the cation of compound A (N-propyl-tetramethoxyacridinium)
    # C: 13 (acridinium core) + 4 (methoxy) + 3 (propyl) = 20
    # H: 5 (acridinium core) + 12 (methoxy) + 7 (propyl) = 24
    # N: 1 (acridinium core)
    # O: 4 (methoxy)
    formula_A = {'C': 20, 'H': 24, 'N': 1, 'O': 4}

    # Atomic composition of the group removed (n-propyl group, -C3H7)
    group_removed = {'C': 3, 'H': 7, 'N': 0, 'O': 0}

    # Atomic composition of the group added (from methyl-3-aminopropionate, -CH2CH2COOCH3)
    # C: 1+1+1+1 = 4
    # H: 2+2+3 = 7
    # O: 2
    group_added = {'C': 4, 'H': 7, 'N': 0, 'O': 2}

    # Calculate the molecular formula of the cation of compound B
    formula_B = {}
    
    # Calculate Carbon atoms
    c_b = formula_A['C'] - group_removed['C'] + group_added['C']
    print(f"Calculation for Carbon (C): {formula_A['C']} - {group_removed['C']} + {group_added['C']} = {c_b}")
    formula_B['C'] = c_b

    # Calculate Hydrogen atoms
    h_b = formula_A['H'] - group_removed['H'] + group_added['H']
    print(f"Calculation for Hydrogen (H): {formula_A['H']} - {group_removed['H']} + {group_added['H']} = {h_b}")
    formula_B['H'] = h_b

    # Calculate Nitrogen atoms
    n_b = formula_A['N'] - group_removed['N'] + group_added['N']
    print(f"Calculation for Nitrogen (N): {formula_A['N']} - {group_removed['N']} + {group_added['N']} = {n_b}")
    formula_B['N'] = n_b
    
    # Calculate Oxygen atoms
    o_b = formula_A['O'] - group_removed['O'] + group_added['O']
    print(f"Calculation for Oxygen (O): {formula_A['O']} - {group_removed['O']} + {group_added['O']} = {o_b}")
    formula_B['O'] = o_b

    # Construct the final formula string
    final_formula = f"C{formula_B['C']}H{formula_B['H']}N{formula_B['N']}O{formula_B['O']}"
    
    print("\nThe molecular formula of compound B is:")
    print(final_formula)

calculate_molecular_formula()