def calculate_molecular_formula():
    """
    Calculates the molecular formula of the product based on the reaction scheme.
    """

    # --- Step 1: Calculate the formula of the Intermediate ---
    print("Step 1: Calculating the molecular formula of the Intermediate.")

    # Molecular formulas of the starting materials
    aminothiazole = {'C': 3, 'H': 4, 'N': 2, 'O': 0, 'S': 1}
    chloro_keto_ester = {'C': 6, 'H': 9, 'Cl': 1, 'O': 3, 'S': 0, 'N': 0}

    # Small molecules eliminated during the reaction
    hcl = {'H': 1, 'Cl': 1}
    h2o = {'H': 2, 'O': 1}

    print("Formula of 2-aminothiazole: C3H4N2S")
    print("Formula of ethyl 2-chloro-3-oxobutanoate: C6H9ClO3")
    print("The reaction is a cyclocondensation, which eliminates one molecule of HCl and one molecule of H2O.")
    
    # Calculate intermediate formula
    intermediate_c = aminothiazole['C'] + chloro_keto_ester['C']
    intermediate_h = aminothiazole['H'] + chloro_keto_ester['H'] - hcl['H'] - h2o['H']
    intermediate_n = aminothiazole['N'] + chloro_keto_ester['N']
    intermediate_o = aminothiazole['O'] + chloro_keto_ester['O'] - h2o['O']
    intermediate_s = aminothiazole['S'] + chloro_keto_ester['S']

    print(f"Intermediate C = {aminothiazole['C']} + {chloro_keto_ester['C']} = {intermediate_c}")
    print(f"Intermediate H = {aminothiazole['H']} + {chloro_keto_ester['H']} - {hcl['H']} (from HCl) - {h2o['H']} (from H2O) = {intermediate_h}")
    print(f"Intermediate N = {aminothiazole['N']} = {intermediate_n}")
    print(f"Intermediate O = {chloro_keto_ester['O']} - {h2o['O']} (from H2O) = {intermediate_o}")
    print(f"Intermediate S = {aminothiazole['S']} = {intermediate_s}")
    
    intermediate_formula = f"C{intermediate_c}H{intermediate_h}N{intermediate_n}O{intermediate_o}S{intermediate_s}"
    print(f"Formula of the Intermediate product is: {intermediate_formula}\n")

    # --- Step 2: Calculate the formula of the Final Product ---
    print("Step 2: Calculating the molecular formula of the final product.")
    
    # Atoms in groups being swapped
    oet_group = {'C': 2, 'H': 5, 'N': 0, 'O': 1, 'S': 0}
    # Benzylamino group is -NH-CH2-Ph
    # H count = 1(NH) + 2(CH2) + 5(Ph) = 8
    # C count = 1(CH2) + 6(Ph) = 7
    nhbn_fragment = {'C': 7, 'H': 8, 'N': 1, 'O': 0, 'S': 0}
    
    print("The second reaction is an amidation, replacing an ethoxy group (-OEt) with an N-benzylamino group (-NHBn).")
    print(f"Lost group (-OEt): C{oet_group['C']}H{oet_group['H']}O{oet_group['O']}")
    print(f"Gained group (-NHBn): C{nhbn_fragment['C']}H{nhbn_fragment['H']}N{nhbn_fragment['N']}")

    # Calculate final product formula
    product_c = intermediate_c - oet_group['C'] + nhbn_fragment['C']
    product_h = intermediate_h - oet_group['H'] + nhbn_fragment['H']
    product_n = intermediate_n - oet_group['N'] + nhbn_fragment['N']
    product_o = intermediate_o - oet_group['O'] + nhbn_fragment['O']
    product_s = intermediate_s - oet_group['S'] + nhbn_fragment['S']

    print(f"Final Product C = {intermediate_c} (from Intermediate) - {oet_group['C']} + {nhbn_fragment['C']} = {product_c}")
    print(f"Final Product H = {intermediate_h} (from Intermediate) - {oet_group['H']} + {nhbn_fragment['H']} = {product_h}")
    print(f"Final Product N = {intermediate_n} (from Intermediate) - {oet_group['N']} + {nhbn_fragment['N']} = {product_n}")
    print(f"Final Product O = {intermediate_o} (from Intermediate) - {oet_group['O']} + {nhbn_fragment['O']} = {product_o}")
    print(f"Final Product S = {intermediate_s} (from Intermediate) - {oet_group['S']} + {nhbn_fragment['S']} = {product_s}")

    final_formula = f"C{product_c}H{product_h}N{product_n}O{product_o}S{product_s if product_s > 1 else ''}"
    print(f"\nThe molecular formula of the product is: {final_formula}")

calculate_molecular_formula()