def solve_molecular_formula():
    """
    Calculates the molecular formula of compound B based on the provided reaction scheme.
    """
    # Step 1: Define the atomic composition of the starting reactants.
    
    # The starting cation is 9-(2,6-dimethoxyphenyl)-1,8-dimethoxyxanthenylium.
    # Xanthene core (C13) + 2 methoxy groups (-OCH3) = C15H12O3 (as 1,8-dimethoxyxanthene)
    # The cation is formed at C9, and a 2,6-dimethoxyphenyl group is attached.
    # 9-Aryl-1,8-dimethoxyxanthenylium cation:
    # C: 13 (xanthene core) + 6 (phenyl ring) + 4 (four -CH3) = 23
    # H: 6 (xanthene core rings) + 3 (phenyl ring) + 12 (four -CH3) = 21
    # O: 1 (xanthene ether) + 4 (four -OCH3) = 5
    starting_cation = {'C': 23, 'H': 21, 'N': 0, 'O': 5}

    # The nucleophile is methyl-3-aminopropionate: NH2-CH2-CH2-COOCH3
    # C: 1+1+1+1 = 4
    # H: 2(NH2) + 2(CH2) + 2(CH2) + 3(CH3) = 9
    # N: 1
    # O: 2
    nucleophile = {'C': 4, 'H': 9, 'N': 1, 'O': 2}

    # Step 2: Define the molecule eliminated during the condensation reaction.
    # The reaction forms an acridinium salt by replacing the xanthene ether oxygen
    # with the nitrogen from the amine, eliminating one molecule of water.
    water = {'C': 0, 'H': 2, 'N': 0, 'O': 1}

    # Step 3: Calculate the atomic composition of the cation of compound B.
    product_cation_B = {}
    product_cation_B['C'] = starting_cation['C'] + nucleophile['C'] - water['C']
    product_cation_B['H'] = starting_cation['H'] + nucleophile['H'] - water['H']
    product_cation_B['N'] = starting_cation['N'] + nucleophile['N'] - water['N']
    product_cation_B['O'] = starting_cation['O'] + nucleophile['O'] - water['O']

    # Step 4: Print the calculation and the final result.
    print("Determining the molecular formula for the cation of compound B.")
    print("-" * 50)
    print(f"Formula of starting cation: C{starting_cation['C']}H{starting_cation['H']}O{starting_cation['O']}")
    print(f"Formula of methyl-3-aminopropionate: C{nucleophile['C']}H{nucleophile['H']}N{nucleophile['N']}O{nucleophile['O']}")
    print(f"Formula of eliminated water: H{water['H']}O{water['O']}")
    print("\nCalculation: (Starting Cation) + (Nucleophile) - (Water)")
    print(f"Carbon (C): {starting_cation['C']} + {nucleophile['C']} - {water['C']} = {product_cation_B['C']}")
    print(f"Hydrogen (H): {starting_cation['H']} + {nucleophile['H']} - {water['H']} = {product_cation_B['H']}")
    print(f"Nitrogen (N): {starting_cation['N']} + {nucleophile['N']} - {water['N']} = {product_cation_B['N']}")
    print(f"Oxygen (O): {starting_cation['O']} + {nucleophile['O']} - {water['O']} = {product_cation_B['O']}")
    print("-" * 50)
    
    # Format the final formula string, omitting elements with a count of 1.
    formula_B_str = f"C{product_cation_B['C']}H{product_cation_B['H']}"
    if product_cation_B['N'] > 0:
        formula_B_str += f"N" if product_cation_B['N'] == 1 else f"N{product_cation_B['N']}"
    if product_cation_B['O'] > 0:
        formula_B_str += f"O{product_cation_B['O']}"
    
    print(f"The molecular formula of the cation of compound B is: {formula_B_str}")

solve_molecular_formula()