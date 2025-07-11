def calculate_product_formula():
    """
    Calculates the molecular formula of the final product based on the reaction scheme.
    """
    # Step 1: Define the formula of the intermediate
    # Intermediate formula is C9H10N2O2S
    intermediate_atoms = {'C': 9, 'H': 10, 'N': 2, 'O': 2, 'S': 1}

    # Step 2: Define the formulas of the reagent added and byproduct removed in the second step
    # Reagent added: Benzylamine (C7H9N)
    benzylamine_atoms = {'C': 7, 'H': 9, 'N': 1, 'O': 0, 'S': 0}
    # Byproduct removed: Ethanol (C2H6O)
    ethanol_atoms = {'C': 2, 'H': 6, 'O': 1, 'N': 0, 'S': 0}

    # Step 3: Calculate the formula of the final product
    # Product = Intermediate + Benzylamine - Ethanol
    product_atoms = {
        'C': intermediate_atoms['C'] + benzylamine_atoms['C'] - ethanol_atoms['C'],
        'H': intermediate_atoms['H'] + benzylamine_atoms['H'] - ethanol_atoms['H'],
        'N': intermediate_atoms['N'] + benzylamine_atoms['N'] - ethanol_atoms['N'],
        'O': intermediate_atoms['O'] + benzylamine_atoms['O'] - ethanol_atoms['O'],
        'S': intermediate_atoms['S'] + benzylamine_atoms['S'] - ethanol_atoms['S']
    }

    # Step 4: Print the calculation and the final formula
    print("To find the molecular formula of the product, we perform atom balance on the second reaction:")
    print("Product = Intermediate + Benzylamine - Ethanol")
    print(f"C: {intermediate_atoms['C']} + {benzylamine_atoms['C']} - {ethanol_atoms['C']} = {product_atoms['C']}")
    print(f"H: {intermediate_atoms['H']} + {benzylamine_atoms['H']} - {ethanol_atoms['H']} = {product_atoms['H']}")
    print(f"N: {intermediate_atoms['N']} + {benzylamine_atoms['N']} - {ethanol_atoms['N']} = {product_atoms['N']}")
    print(f"O: {intermediate_atoms['O']} + {benzylamine_atoms['O']} - {ethanol_atoms['O']} = {product_atoms['O']}")
    print(f"S: {intermediate_atoms['S']} + {benzylamine_atoms['S']} - {ethanol_atoms['S']} = {product_atoms['S']}")
    
    # Format the final formula string
    formula_str = (f"C{product_atoms['C']}H{product_atoms['H']}"
                   f"N{product_atoms['N']}O{product_atoms['O']}S{product_atoms['S']}")
    # Adjust for atoms with count 1
    formula_str = formula_str.replace('O1', 'O').replace('S1', 'S')
    
    print(f"\nThe molecular formula of the product is: {formula_str}")

calculate_product_formula()