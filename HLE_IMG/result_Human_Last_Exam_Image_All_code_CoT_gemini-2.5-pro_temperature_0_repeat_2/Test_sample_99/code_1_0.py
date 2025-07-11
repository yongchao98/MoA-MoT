def calculate_molecular_formula():
    """
    Calculates the molecular formula of the final product based on the reaction scheme.
    """
    # Molecular formula of the intermediate: C9H10N2O2S
    intermediate = {'C': 9, 'H': 10, 'N': 2, 'O': 2, 'S': 1}

    # Molecular formula of benzylamine (C6H5CH2NH2): C7H9N
    benzylamine = {'C': 7, 'H': 9, 'N': 1, 'O': 0, 'S': 0}

    # Molecular formula of the byproduct, ethanol (CH3CH2OH): C2H6O
    ethanol = {'C': 2, 'H': 6, 'O': 1, 'S': 0}

    # Calculate the molecular formula of the final product
    product = {}
    print("Calculating the molecular formula of the final product:")
    print("Product = Intermediate + Benzylamine - Ethanol")
    print("-" * 30)

    # Carbon
    product['C'] = intermediate['C'] + benzylamine['C'] - ethanol['C']
    print(f"Carbon (C): {intermediate['C']} + {benzylamine['C']} - {ethanol['C']} = {product['C']}")

    # Hydrogen
    product['H'] = intermediate['H'] + benzylamine['H'] - ethanol['H']
    print(f"Hydrogen (H): {intermediate['H']} + {benzylamine['H']} - {ethanol['H']} = {product['H']}")

    # Nitrogen
    product['N'] = intermediate['N'] + benzylamine['N'] - ethanol['N']
    print(f"Nitrogen (N): {intermediate['N']} + {benzylamine['N']} - {ethanol['N']} = {product['N']}")

    # Oxygen
    product['O'] = intermediate['O'] + benzylamine['O'] - ethanol['O']
    print(f"Oxygen (O): {intermediate['O']} + {benzylamine['O']} - {ethanol['O']} = {product['O']}")

    # Sulfur
    product['S'] = intermediate['S'] + benzylamine['S'] - ethanol['S']
    print(f"Sulfur (S): {intermediate['S']} + {benzylamine['S']} - {ethanol['S']} = {product['S']}")

    # Construct the final formula string
    formula = ""
    for atom in ['C', 'H', 'N', 'O', 'S']:
        count = product[atom]
        if count > 0:
            formula += atom
            if count > 1:
                formula += str(count)
    
    print("-" * 30)
    print(f"The final molecular formula is: {formula}")

calculate_molecular_formula()